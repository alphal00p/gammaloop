//! Collect metric/vector components before materializing their products.

mod factorized;
pub(crate) use factorized::{ContractionStatus, FactorizedContraction};

use super::{Endpoint, SlotContraction};
use crate::{
    shorthands::schoonschip::SimplificationCandidates,
    tensor::{SymbolicTensor, inference::InterfaceInference},
};
use ahash::AHashMap;
use spenso::shadowing::{TermLeaf, TermTape};
use spenso::structure::{
    abstract_index::AbstractIndex,
    partial::{PartialStructure, PartialStructureExt},
    representation::{LibraryRep, RepName},
    slot::{SlotMatch, SlotMatcher},
};
use symbolica::{
    atom::{Atom, AtomCore, AtomView, FunctionBuilder, Symbol, representation::FunView},
    domains::{
        Ring,
        rational::{Q, Rational},
    },
    poly::{PolyVariable, polynomial::MultivariatePolynomial},
};

#[derive(Clone, Copy, PartialEq, Eq, Hash)]
enum Argument<'a> {
    Original(AtomView<'a>),
    Vector(usize, usize),
}

#[derive(Clone, Copy, PartialEq, Eq, Hash)]
enum TensorSource<'a> {
    Function(Symbol),
    Literal(AtomView<'a>),
    Factor(usize),
}

#[derive(Clone, PartialEq, Eq, Hash)]
enum Variable<'a> {
    Scalar(AtomView<'a>),
    Subtree(usize),
    Dot(usize, [usize; 2]),
    Vector(usize, AtomView<'a>),
    Metric([AtomView<'a>; 2]),
    Tensor(TensorSource<'a>, Vec<Argument<'a>>),
}

#[derive(Clone, Copy)]
enum Terminal<'a> {
    Vector(usize),
    Free(AtomView<'a>),
    Tensor(usize, usize),
}

struct Node<'a> {
    parent: usize,
    degree: u8,
    space: usize,
    slot: AtomView<'a>,
    terminals: Vec<Terminal<'a>>,
    metrics: usize,
}

// Only a key for the existing tensor plan: exact non-internal arguments and
// (argument position, space) for internal ports. No new tensor semantics.
type TensorKey<'a> = (
    TensorSource<'a>,
    Vec<(usize, Argument<'a>)>,
    Vec<(usize, usize)>,
);

#[derive(PartialEq, Eq, Hash)]
struct MonomialKey {
    factors: Vec<(usize, u16)>,
    tensors: Vec<usize>,
    connections: Vec<(usize, [(usize, usize); 2])>,
}

// Leaves identify borrowed payloads or an opaque subtree of the existing graph.
// Only explicit expansion admits Atom arithmetic directly.
#[derive(Clone, Copy, PartialEq, Eq, Hash)]
enum InputLeaf<'a> {
    Atom(AtomView<'a>),
    // A completed graph scope with an established empty boundary. Its literal
    // value, rather than its occurrence ID, identifies an opaque coefficient.
    Scalar(AtomView<'a>),
    Subtree(usize),
}

#[derive(Clone, Copy, PartialEq, Eq)]
pub(super) enum Intake {
    Contraction,
    ScalarExpansion,
    #[cfg(test)]
    ExpandedContraction,
    #[cfg(test)]
    FactoredContraction,
}

// Scratch storage for this one coefficient-list emission, not a tensor interface.
struct ComponentSum<'a, 'b> {
    contractor: &'b SlotContraction,
    intake: Intake,
    rank_one: bool,
    slots: &'b mut SlotMatcher,
    heads: AHashMap<Symbol, bool>,
    vectors: Vec<FunView<'a>>,
    vector_sources: AHashMap<AtomView<'a>, usize>,
    vector_positions: AHashMap<(Symbol, Vec<AtomView<'a>>), usize>,
    spaces: Vec<(LibraryRep, AtomView<'a>)>,
    // Resolved syntax survives term resets; incidence counts and node IDs do not.
    endpoints: AHashMap<AtomView<'a>, (usize, AtomView<'a>)>,
    compact_spaces: AHashMap<AtomView<'a>, usize>,
    variables: Vec<Variable<'a>>,
    positions: AHashMap<Variable<'a>, usize>,
    nodes: Vec<Node<'a>>,
    occurrences: AHashMap<(usize, AtomView<'a>), usize>,
    factors: Vec<(usize, u16)>,
    coefficient: Rational,
    contracted: bool,
    inference: InterfaceInference,
    tensor_ports: AHashMap<AtomView<'a>, Vec<(usize, AtomView<'a>)>>,
    tensors: Vec<(TensorSource<'a>, Vec<Argument<'a>>)>,
    literal_relabellings: std::cell::RefCell<Vec<(Atom, Atom)>>,
    opaque_factors: Vec<(usize, Vec<AtomView<'a>>, PartialStructure)>,
    factor_roots: Vec<Option<usize>>,
    factor_emission: Option<&'b dyn Fn(usize) -> Option<Atom>>,
    overrides: Vec<(AtomView<'a>, Argument<'a>)>,
    metrics: usize,
    alpha_tensors: AHashMap<TensorKey<'a>, usize>,
    alpha_terms: AHashMap<MonomialKey, Vec<(usize, u16)>>,
    alpha_bytes: usize,
    coefficients: Vec<Rational>,
    monomials: Vec<Vec<(usize, u16)>>,
    input: TermTape<InputLeaf<'a>>,
    scalar_spectators: AHashMap<AtomView<'a>, bool>,
}

impl SlotContraction {
    pub(super) fn materialize_scalar_sum(
        &self,
        value: AtomView<'_>,
        slots: &mut SlotMatcher,
    ) -> Option<Atom> {
        if !matches!(value, AtomView::Add(_))
            || value.needs_normalization()
            || self.metric.get_evaluation_info().is_some()
        {
            return None;
        }
        let mut state = ComponentSum::new(self, Intake::ScalarExpansion, slots);
        let (node, _) = state.compile_input(value, 0)?;
        state
            .input
            .distribute(&mut vec![(node, 1)], &Rational::one(), 0)?;
        let terms = state.input.take_terms();
        state.coefficients.reserve(terms.len());
        state.monomials.reserve(terms.len());
        for (factors, coefficient) in terms {
            state.factors.clear();
            for (factor, exponent) in factors {
                // Compile-time admission proved every literal scalar. This
                // explicit expansion neither contracts nor rebuilds functions.
                let InputLeaf::Atom(value) = state.input.leaves()[factor] else {
                    return None;
                };
                state.variable(Variable::Scalar(value), exponent);
            }
            state.coefficients.push(coefficient);
            state.monomials.push(state.factors.clone());
        }
        state.emit_polynomial()
    }

    /// Explicit test oracle for the flat reducer. No production contraction
    /// dispatch reaches polynomial emission; the factorized engine owns that.
    #[cfg(test)]
    fn materialize_test_sum(
        &self,
        value: AtomView<'_>,
        slots: &mut SlotMatcher,
        intake: Intake,
    ) -> Option<Atom> {
        if intake == Intake::ScalarExpansion {
            return self.materialize_scalar_sum(value, slots);
        }
        let expand_sums = intake != Intake::ExpandedContraction;
        if !matches!(value, AtomView::Add(_)) && intake != Intake::FactoredContraction {
            return None;
        }
        if value.needs_normalization() || self.metric.get_evaluation_info().is_some() {
            return None;
        }
        let mut state = ComponentSum::new(self, intake, slots);
        let mut spectators = Vec::new();
        if expand_sums {
            let mut pending = Vec::new();
            if let AtomView::Mul(product) = value {
                for factor in product.iter() {
                    if state.scalar_spectator(factor) {
                        spectators.push(factor);
                    } else {
                        pending.push((factor, 1));
                    }
                }
            } else {
                pending.push((value, 1));
            }
            if pending.is_empty() {
                return None;
            }
            let mut count = (1usize, 0usize);
            let mut compiled = Vec::with_capacity(pending.len());
            for (factor, exponent) in pending {
                let (node, next) = state.compile_input(factor, 0)?;
                count = TermTape::<InputLeaf<'_>>::product_size(count, next)?;
                compiled.push((node, exponent));
            }
            state.input.distribute(&mut compiled, &Rational::one(), 0)?;
            // Match expansion's exact collection before checking graph degree
            // or opaque powers. Cancelled terms cannot create incidences.
            let terms = state.input.take_terms();
            state.coefficients.reserve(terms.len());
            state.monomials.reserve(terms.len());
            for (factors, coefficient) in terms {
                state.nodes.clear();
                state.occurrences.clear();
                state.factors.clear();
                state.tensors.clear();
                state.metrics = 0;
                state.coefficient = coefficient;
                for (factor, exponent) in factors {
                    let InputLeaf::Atom(value) = state.input.leaves()[factor] else {
                        return None;
                    };
                    state.factor(value, exponent)?;
                }
                state.reduce_components()?;
                state.coefficients.push(state.coefficient.clone());
                state.monomials.push(state.factors.clone());
            }
        } else {
            let AtomView::Add(sum) = value else {
                unreachable!()
            };
            state.coefficients.reserve(sum.iter().len());
            state.monomials.reserve(sum.iter().len());
            for term in sum.iter() {
                state.nodes.clear();
                state.occurrences.clear();
                state.factors.clear();
                state.tensors.clear();
                state.metrics = 0;
                state.coefficient = Rational::one();
                if let AtomView::Mul(product) = term {
                    for factor in product.iter() {
                        state.factor(factor, 1)?;
                    }
                } else {
                    state.factor(term, 1)?;
                }
                state.reduce_components()?;
                state.coefficients.push(state.coefficient.clone());
                state.monomials.push(state.factors.clone());
            }
        }
        let result = state.emit_polynomial()?;
        if spectators.is_empty() {
            Some(result)
        } else {
            Some(Atom::mul_many(
                spectators
                    .into_iter()
                    .map(|value| value.to_owned())
                    .chain([result]),
            ))
        }
    }
}

impl<'a, 'b> ComponentSum<'a, 'b> {
    fn emit_polynomial(self) -> Option<Atom> {
        if !self.contracted && !self.input.distributed() {
            return None;
        }
        if self.coefficients.is_empty() {
            return Some(Atom::Zero);
        }
        // Canonical variable aliases are resolved by the existing emitter.
        // A conservative degree bound proves even a complete alias collision
        // cannot overflow after materialization starts.
        for factors in &self.monomials {
            factors
                .iter()
                .try_fold(0u16, |degree, (_, exponent)| degree.checked_add(*exponent))?;
        }
        // Dense polynomial exponents suit the few scalar variables in a ladder,
        // but not high-entropy sums with a new variable in every monomial.
        const MAX_EXPONENT_BYTES: usize = 64 * 1024 * 1024;
        if self
            .coefficients
            .len()
            .checked_mul(self.variables.len())?
            .checked_mul(std::mem::size_of::<u16>())?
            > MAX_EXPONENT_BYTES
        {
            return None;
        }
        // Admission has finished: no callback-bearing or unsupported factor can
        // reach a builder. Canonical Atom equality also merges any variable alias.
        let mut atoms = AHashMap::new();
        let mut variables: Vec<PolyVariable> = Vec::new();
        let mut remap = Vec::with_capacity(self.variables.len());
        for variable in &self.variables {
            let atom = self.emit_variable(variable)?;
            let next = variables.len();
            let position = *atoms.entry(atom.clone()).or_insert_with(|| {
                variables.push(atom.try_into().unwrap());
                next
            });
            remap.push(position);
        }
        let count = variables.len();
        let size = self.coefficients.len().checked_mul(count)?;
        let mut exponents = vec![0u16; size];
        for (term, factors) in self.monomials.iter().enumerate() {
            for &(factor, exponent) in factors {
                let entry = &mut exponents[term * count + remap[factor]];
                *entry = entry.checked_add(exponent)?;
            }
        }
        Some(
            MultivariatePolynomial::<_, u16>::from_coefficient_list(
                self.coefficients,
                exponents,
                variables.into(),
                &Q,
            )
            .to_expression(),
        )
    }

    fn new(contractor: &'b SlotContraction, intake: Intake, slots: &'b mut SlotMatcher) -> Self {
        Self {
            contractor,
            intake,
            rank_one: true,
            slots,
            heads: AHashMap::new(),
            vectors: Vec::new(),
            vector_sources: AHashMap::new(),
            vector_positions: AHashMap::new(),
            spaces: Vec::new(),
            endpoints: AHashMap::new(),
            compact_spaces: AHashMap::new(),
            variables: Vec::new(),
            positions: AHashMap::new(),
            nodes: Vec::new(),
            occurrences: AHashMap::new(),
            factors: Vec::new(),
            coefficient: Rational::one(),
            contracted: false,
            inference: InterfaceInference::default(),
            tensor_ports: AHashMap::new(),
            tensors: Vec::new(),
            literal_relabellings: Default::default(),
            opaque_factors: Vec::new(),
            factor_roots: Vec::new(),
            factor_emission: None,
            overrides: Vec::new(),
            metrics: 0,
            alpha_tensors: AHashMap::new(),
            alpha_terms: AHashMap::new(),
            alpha_bytes: 0,
            coefficients: Vec::new(),
            monomials: Vec::new(),
            input: TermTape::default(),
            scalar_spectators: AHashMap::new(),
        }
    }

    fn emit_variable(&self, variable: &Variable<'a>) -> Option<Atom> {
        Some(match variable {
            Variable::Scalar(value) => (*value).to_owned(),
            Variable::Subtree(position) => {
                (self.factor_emission?)(self.opaque_factors[*position].0)?
            }
            Variable::Vector(vector, slot) => self.emit_vector(*vector, *slot),
            Variable::Metric([first, second]) => FunctionBuilder::new(self.contractor.metric)
                .add_arg(*first)
                .add_arg(*second)
                .finish(),
            Variable::Dot(space, [first, second]) => {
                let (representation, dimension) = self.spaces[*space];
                let compact = representation.to_symbolic([dimension]);
                let first = self.emit_vector(*first, compact.as_view());
                let second = self.emit_vector(*second, compact.as_view());
                FunctionBuilder::new(self.contractor.metric)
                    .add_arg(first)
                    .add_arg(second)
                    .finish()
            }
            Variable::Tensor(
                source @ (TensorSource::Function(_) | TensorSource::Literal(_)),
                arguments,
            ) => {
                let head = match source {
                    TensorSource::Function(head) => *head,
                    TensorSource::Literal(AtomView::Fun(function)) => function.get_symbol(),
                    _ => unreachable!("literal tensor source is a function"),
                };
                let mut function = FunctionBuilder::new(head);
                for argument in arguments {
                    function = match *argument {
                        Argument::Original(value) => function.add_arg(value),
                        Argument::Vector(space, head) => {
                            let (representation, dimension) = self.spaces[space];
                            let compact = representation.to_symbolic([dimension]);
                            function.add_arg(self.emit_vector(head, compact.as_view()))
                        }
                    };
                }
                let result = function.finish();
                if let TensorSource::Literal(source) = source
                    && *source != result.as_view()
                {
                    self.literal_relabellings
                        .borrow_mut()
                        .push((source.to_owned(), result.clone()));
                }
                result
            }
            Variable::Tensor(TensorSource::Factor(position), arguments) => {
                return self.emit_factor(*position, arguments);
            }
        })
    }

    // Count borrowed factor occurrences, not expanded Atom bytes. Both the
    // Cartesian work and the existing coefficient/exponent storage are bounded.
    const MAX_GENERATED_TERMS: usize = 1_000_000;
    const MAX_EXPANSION_DEPTH: usize = 256;

    fn compile_group(
        &mut self,
        children: impl ExactSizeIterator<Item = AtomView<'a>>,
        sum: bool,
        depth: usize,
    ) -> Option<(usize, (usize, usize))> {
        let children = children
            .map(|child| self.compile_input(child, depth + 1))
            .collect::<Option<Vec<_>>>()?;
        self.input.group(children.into_iter(), sum)
    }

    fn compile_input(
        &mut self,
        value: AtomView<'a>,
        depth: usize,
    ) -> Option<(usize, (usize, usize))> {
        if depth > Self::MAX_EXPANSION_DEPTH || value.needs_normalization() {
            return None;
        }
        match value {
            AtomView::Add(sum) => self.compile_group(sum.iter(), true, depth),
            AtomView::Mul(product) => self.compile_group(product.iter(), false, depth),
            AtomView::Pow(power) => {
                let (base, exponent) = power.get_base_exp();
                let exponent = u16::try_from(i64::try_from(exponent).ok()?).ok()?;
                if exponent == 0 {
                    return None;
                }
                let base = self.compile_input(base, depth + 1)?;
                self.input.power(base, exponent)
            }
            _ => {
                if let Some(node) = self.input.known_leaf(&InputLeaf::Atom(value)) {
                    return Some((node, (1, 1)));
                }
                if self.intake == Intake::ScalarExpansion {
                    match value {
                        AtomView::Num(_) => {}
                        AtomView::Var(variable) => {
                            let head = variable.get_symbol();
                            let tags = self.contractor.tags;
                            if head.get_wildcard_level() != 0
                                || head.has_tag(&tags.tensor)
                                || head.has_tag(&tags.rank1)
                                || head.has_tag(&tags.representation)
                                || head == tags.chain_in
                                || head == tags.chain_out
                            {
                                return None;
                            }
                        }
                        // The metric evaluation guard and compact-dot grammar
                        // prove every head intrinsic, including both operands.
                        AtomView::Fun(function) if self.compact_dot(function).is_some() => {}
                        _ => return None,
                    }
                } else if matches!(value, AtomView::Fun(_))
                    && !InterfaceInference::normalization_is_intrinsic(value)
                {
                    return None;
                }
                // Cache only certified original leaves. Compilation finishes
                // before distribution, so cancellation cannot hide a callback.
                if let AtomView::Num(_) = value {
                    self.input
                        .number(InputLeaf::Atom(value), Rational::try_from(value).ok()?)
                } else {
                    self.input.leaf(InputLeaf::Atom(value))
                }
            }
        }
    }

    fn scalar_spectator(&mut self, value: AtomView<'a>) -> bool {
        if let Some(&scalar) = self.scalar_spectators.get(&value) {
            return scalar;
        }
        // Do not mistake an internally contracted tensor product for an inert
        // scalar, or hide pending work inside scalar metadata. Use the shared
        // interface owner rather than a second tensor/metadata grammar.
        let candidates = SimplificationCandidates::scan(value, [], || true);
        let compact_power = if let AtomView::Pow(power) = value
            && let (AtomView::Fun(base), exponent) = power.get_base_exp()
            && i64::try_from(exponent).is_ok_and(|n| n > 0)
        {
            self.compact_dot(base).is_some()
        } else {
            false
        };
        let mut exact = true;
        value.visitor(&mut |node| {
            exact &= match node {
                AtomView::Num(_) => Rational::try_from(node).is_ok(),
                AtomView::Var(variable) => variable.get_symbol().get_wildcard_level() == 0,
                _ => true,
            };
            exact
        });
        // Omitted root factors do not enter compile_input. This complete
        // observer certificate covers all of their function metadata instead.
        let scalar = candidates.complete
            && candidates.intrinsic
            && (candidates.normalized() || compact_power)
            && exact
            && SymbolicTensor::validate_interface(
                &value.to_owned(),
                &PartialStructure::from_logical_slots([]),
            )
            .is_ok();
        if self.scalar_spectators.len() < 256 {
            self.scalar_spectators.insert(value, scalar);
        }
        scalar
    }

    // Bound both stored keys and transient signature work independently of the
    // existing dense-polynomial cap. Exhaustion declines the entire sum before
    // materialization; partial collection could expose more work on a rerun.
    // Exercise saturation through the actual public pipeline with a small
    // bounded fixture; no process-global mutable test configuration is needed.
    const MAX_ALPHA_BYTES: usize = if cfg!(test) {
        16 * 1024
    } else {
        16 * 1024 * 1024
    };

    fn reserve_alpha(&mut self, bytes: usize) -> Result<(), ()> {
        let total = self.alpha_bytes.checked_add(bytes).ok_or(())?;
        if total > Self::MAX_ALPHA_BYTES {
            return Err(());
        }
        self.alpha_bytes = total;
        Ok(())
    }

    fn monomial_key(&mut self) -> Result<Option<MonomialKey>, ()> {
        if self.tensors.is_empty()
            || self
                .tensors
                .iter()
                .any(|(source, _)| matches!(source, TensorSource::Factor(_)))
        {
            return Ok(None);
        }
        // These estimates include Vec headers and generous map-capacity slack;
        // borrowed Atoms retain the input's storage rather than copying bytes.
        let argument_count = self
            .tensors
            .iter()
            .try_fold(0usize, |total, (_, args)| total.checked_add(args.len()))
            .ok_or(())?;
        let transient = (|| {
            argument_count
                .checked_mul(4 * std::mem::size_of::<(usize, Argument<'a>)>())?
                .checked_add(self.factors.len().checked_mul(32)?)?
                .checked_add(self.tensors.len().checked_mul(256)?)
        })()
        .ok_or(())?;
        if transient > Self::MAX_ALPHA_BYTES {
            return Err(());
        }
        let mut ports = vec![Vec::new(); self.tensors.len()];
        let mut connections = Vec::new();
        for (index, node) in self.nodes.iter().enumerate() {
            if node.parent != index {
                continue;
            }
            if let &[Terminal::Tensor(first, a), Terminal::Tensor(second, b)] =
                node.terminals.as_slice()
            {
                ports[first].push((a, node.space));
                ports[second].push((b, node.space));
                connections.push((node.space, [(first, a), (second, b)]));
            }
        }
        if connections.is_empty() {
            return Ok(None);
        }
        let mut tensors = Vec::with_capacity(self.tensors.len());
        for (tensor, mut internal) in ports.into_iter().enumerate() {
            internal.sort_unstable();
            let mut arguments = Vec::new();
            let mut next_internal = 0;
            for position in 0..self.tensors[tensor].1.len() {
                if internal
                    .get(next_internal)
                    .is_some_and(|&(p, _)| p == position)
                {
                    next_internal += 1;
                    continue;
                }
                let mut argument = self.tensors[tensor].1[position];
                // Existing compact bindings and newly contracted vectors emit
                // the same Atom and must have the same signature.
                if let Argument::Original(AtomView::Fun(function)) = argument
                    && let Some(head) = self.vector(function)
                    && let Some(space) = self.compact_space(function.iter().last().unwrap())
                {
                    argument = Argument::Vector(space, head);
                }
                arguments.push((position, argument));
            }
            let key = (self.tensors[tensor].0, arguments, internal);
            let id = if let Some(&id) = self.alpha_tensors.get(&key) {
                id
            } else {
                let bytes = (|| {
                    key.1
                        .capacity()
                        .checked_mul(std::mem::size_of::<(usize, Argument<'a>)>())?
                        .checked_add(
                            key.2
                                .capacity()
                                .checked_mul(std::mem::size_of::<(usize, usize)>())?,
                        )?
                        .checked_add(4 * std::mem::size_of::<(TensorKey<'a>, usize)>())
                })()
                .ok_or(())?;
                self.reserve_alpha(bytes)?;
                let id = self.alpha_tensors.len();
                self.alpha_tensors.insert(key, id);
                id
            };
            tensors.push(id);
        }
        for (_, endpoints) in &mut connections {
            for (tensor, _) in endpoints.iter_mut() {
                *tensor = tensors[*tensor];
            }
            endpoints.sort_unstable();
        }
        tensors.sort_unstable();
        // Repeated indistinguishable factors require general graph labeling.
        // Keep their existing literal contraction result instead.
        if tensors.windows(2).any(|pair| pair[0] == pair[1]) {
            return Ok(None);
        }
        connections.sort_unstable();
        let mut factors = self.factors.clone();
        factors.sort_unstable();
        let mut combined: Vec<(usize, u16)> = Vec::with_capacity(factors.len());
        for (factor, exponent) in factors {
            if let Some((last, power)) = combined.last_mut()
                && *last == factor
            {
                *power = power.checked_add(exponent).ok_or(())?;
            } else {
                combined.push((factor, exponent));
            }
        }
        Ok(Some(MonomialKey {
            factors: combined,
            tensors,
            connections,
        }))
    }

    fn plain(&mut self, head: Symbol) -> bool {
        *self.heads.entry(head).or_insert_with(|| {
            head.get_wildcard_level() == 0
                && !head.is_symmetric()
                && !head.is_antisymmetric()
                && !head.is_cyclesymmetric()
                && !head.is_linear()
                && !head.is_flat()
                && head.get_normalization_function().is_none()
                && head.get_evaluation_info().is_none()
        })
    }

    fn vector(&mut self, function: FunView<'a>) -> Option<usize> {
        if let Some(&vector) = self.vector_sources.get(&function.as_view()) {
            return Some(vector);
        }
        // The shared matcher keeps representation-shaped metadata opaque. The
        // full function template, excluding its final port, identifies a vector;
        // momentum labels such as p(1, slot) must survive state merging.
        self.slots.vector_argument(function)?;
        let head = function.get_symbol();
        if head.is_scalar()
            || !head.has_tag(&self.contractor.tags.rank1)
            || !self.plain(head)
            || (self.intake != Intake::Contraction && function.get_nargs() != 1)
        {
            return None;
        }
        let metadata = function
            .iter()
            .take(function.get_nargs() - 1)
            .collect::<Vec<_>>();
        for &argument in &metadata {
            let observed = SimplificationCandidates::scan(argument, [], || true);
            if !observed.complete || !observed.intrinsic || !observed.normalized() {
                return None;
            }
        }
        let key = (head, metadata);
        let vector = if let Some(&vector) = self.vector_positions.get(&key) {
            vector
        } else {
            self.input.reserve_bytes(
                key.1
                    .capacity()
                    .checked_mul(std::mem::size_of::<AtomView<'a>>())?
                    .checked_add(4 * std::mem::size_of::<(Symbol, Vec<AtomView<'a>>, usize)>())?,
            )?;
            let vector = self.vectors.len();
            self.vectors.push(function);
            self.vector_positions.insert(key, vector);
            vector
        };
        self.input
            .reserve_bytes(4 * std::mem::size_of::<(AtomView<'a>, usize)>())?;
        self.vector_sources.insert(function.as_view(), vector);
        Some(vector)
    }

    fn emit_vector(&self, vector: usize, port: AtomView<'_>) -> Atom {
        let source = self.vectors[vector];
        source
            .iter()
            .take(source.get_nargs() - 1)
            .fold(
                FunctionBuilder::new(source.get_symbol()),
                |builder, argument| builder.add_arg(argument),
            )
            .add_arg(port)
            .finish()
    }

    fn atomic(value: AtomView<'_>) -> bool {
        match value {
            AtomView::Var(v) => v.get_symbol().get_wildcard_level() == 0,
            AtomView::Num(_) => i64::try_from(value).is_ok_and(|v| v >= 0),
            _ => false,
        }
    }

    fn space(
        &mut self,
        value: AtomView<'a>,
        representation: LibraryRep,
        dimension: AtomView<'a>,
    ) -> Option<usize> {
        let AtomView::Fun(function) = value else {
            return None;
        };
        if !representation.is_self_dual()
            || !representation.is_base()
            || representation == LibraryRep::Dummy
            || !Self::atomic(dimension)
            || function.get_symbol() != representation.symbol()
            || !self.plain(function.get_symbol())
        {
            return None;
        }
        let key = (representation, dimension);
        if let Some(position) = self.spaces.iter().position(|&space| space == key) {
            return Some(position);
        }
        self.spaces.push(key);
        Some(self.spaces.len() - 1)
    }

    fn compact_space(&mut self, value: AtomView<'a>) -> Option<usize> {
        if let Some(&space) = self.compact_spaces.get(&value) {
            return Some(space);
        }
        let view = self.slots.compact_representation(value)?;
        let representation = self
            .slots
            .parse_representation::<LibraryRep>(value)
            .ok()?
            .rep;
        let space = self.space(value, representation, view.dimension())?;
        self.compact_spaces.insert(value, space);
        Some(space)
    }

    fn compact_dot(&mut self, function: FunView<'a>) -> Option<(usize, usize, usize)> {
        if function.get_symbol() != self.contractor.metric || function.get_nargs() != 2 {
            return None;
        }
        let mut arguments = function.iter();
        let (AtomView::Fun(first), AtomView::Fun(second)) = (arguments.next()?, arguments.next()?)
        else {
            return None;
        };
        let first_head = self.vector(first)?;
        let second_head = self.vector(second)?;
        let space = self.compact_space(first.iter().last()?)?;
        (space == self.compact_space(second.iter().last()?)?).then_some((
            space,
            first_head,
            second_head,
        ))
    }

    fn resolve_endpoint(&mut self, value: AtomView<'a>) -> Option<(usize, AtomView<'a>)> {
        if let Some(&endpoint) = self.endpoints.get(&value) {
            return Some(endpoint);
        }
        let endpoint = Endpoint::parse(value, self.slots)?;
        if !Self::atomic(endpoint.index) {
            // Named and scoped indices are explicit identities admitted by the
            // shared slot grammar, not tensor payloads to distribute or rewrite.
            let index = self
                .slots
                .parse::<LibraryRep, AbstractIndex>(value)
                .ok()?
                .aind;
            let mut base = index;
            while let AbstractIndex::Scoped(scope) = base {
                base = scope.index();
            }
            if matches!(base, AbstractIndex::Open { .. })
                || !matches!(index, AbstractIndex::Named(..) | AbstractIndex::Scoped(..))
                || endpoint
                    .index
                    .get_all_symbols(true)
                    .iter()
                    .any(|symbol| symbol.get_wildcard_level() != 0)
            {
                return None;
            }
        }
        let space = self.space(value, endpoint.representation, endpoint.dimension)?;
        let resolved = (space, endpoint.index);
        self.endpoints.insert(value, resolved);
        Some(resolved)
    }

    fn endpoint(
        &mut self,
        value: AtomView<'a>,
        space: usize,
        index: AtomView<'a>,
    ) -> Option<usize> {
        let next = self.nodes.len();
        let node = *self.occurrences.entry((space, index)).or_insert_with(|| {
            self.nodes.push(Node {
                parent: next,
                degree: 0,
                space,
                slot: value,
                terminals: Vec::new(),
                metrics: 0,
            });
            next
        });
        self.nodes[node].degree += 1;
        if self.nodes[node].degree > 2 {
            return None;
        }
        Some(node)
    }

    fn root(&mut self, mut node: usize) -> usize {
        while self.nodes[node].parent != node {
            let parent = self.nodes[node].parent;
            self.nodes[node].parent = self.nodes[parent].parent;
            node = parent;
        }
        node
    }

    fn variable(&mut self, variable: Variable<'a>, exponent: u16) {
        let position = if let Some(&position) = self.positions.get(&variable) {
            position
        } else {
            let position = self.variables.len();
            self.variables.push(variable.clone());
            self.positions.insert(variable, position);
            position
        };
        self.factors.push((position, exponent));
    }

    fn dot(&mut self, space: usize, first: usize, second: usize, exponent: u16) {
        let heads = if first <= second {
            [first, second]
        } else {
            [second, first]
        };
        self.variable(Variable::Dot(space, heads), exponent);
    }

    fn tensor(&mut self, function: FunView<'a>) -> Option<()> {
        let head = function.get_symbol();
        if !head.has_tag(&self.contractor.tags.tensor)
            || (self.rank_one && head.has_tag(&self.contractor.tags.rank1))
            || head.is_scalar()
            || !self.plain(head)
            || !matches!(self.slots.classify(function.as_view()), SlotMatch::Other)
        {
            return None;
        }
        let value = function.as_view();
        let ports = if let Some(ports) = self.tensor_ports.get(&value) {
            ports.clone()
        } else {
            if !self.inference.algebra_preserves_leaf_interfaces(value) {
                return None;
            }
            let mut ports = Vec::new();
            for (position, argument) in function.iter().enumerate() {
                // Metric-only vectors use the same explicit tensor terminal as
                // other tensors. Their preceding arguments are opaque metadata,
                // as established by vector admission; only the final port participates.
                if head.has_tag(&self.contractor.tags.rank1) && position + 1 != function.get_nargs()
                {
                    continue;
                }
                if matches!(self.slots.classify(argument), SlotMatch::Explicit(_)) {
                    self.resolve_endpoint(argument)?;
                    ports.push((position, argument));
                } else if let AtomView::Fun(vector) = argument
                    && self.vector(vector).is_some()
                {
                    self.compact_space(vector.iter().last().unwrap())?;
                } else {
                    // Scalar metadata stays opaque to the port graph, but the
                    // ordinary contractor would still visit hidden products.
                    let metadata = SimplificationCandidates::scan(argument, [], || true);
                    if !metadata.complete || !metadata.intrinsic || !metadata.normalized() {
                        return None;
                    }
                }
            }
            self.tensor_ports.insert(value, ports.clone());
            ports
        };
        let tensor = self.tensors.len();
        self.tensors.push((
            if head == spenso::tensor_symbol!("idenso::tensor_alias") {
                // Different literal definitions may share an owner/head. Keep
                // the exact registered use in the key and its rewrite record.
                TensorSource::Literal(value)
            } else {
                TensorSource::Function(head)
            },
            function.iter().map(|arg| self.argument(arg)).collect(),
        ));
        for (position, slot) in ports {
            if let Argument::Original(slot) = self.argument(slot) {
                let (space, index) = self.resolve_endpoint(slot)?;
                let node = self.endpoint(slot, space, index)?;
                self.nodes[node]
                    .terminals
                    .push(Terminal::Tensor(tensor, position));
            }
        }
        Some(())
    }

    fn factor(&mut self, value: AtomView<'a>, exponent: u16) -> Option<()> {
        match value {
            AtomView::Num(_) => {
                self.coefficient *= Q.pow(&Rational::try_from(value).ok()?, u64::from(exponent));
            }
            AtomView::Var(v) if v.get_symbol().get_wildcard_level() == 0 => {
                self.variable(Variable::Scalar(value), exponent);
            }
            AtomView::Pow(power) if exponent == 1 => {
                let (base, power) = power.get_base_exp();
                let power = u16::try_from(i64::try_from(power).ok()?).ok()?;
                if power == 0 {
                    return None;
                }
                // Match DotNormalizer's local positive-power pairing before
                // counting residual explicit ports in the surrounding product.
                self.factor(base, power)?;
            }
            AtomView::Fun(function) if function.get_symbol() == self.contractor.metric => {
                let mut args = function.iter();
                if args.len() != 2 {
                    return None;
                }
                let first = args.next().unwrap();
                let second = args.next().unwrap();
                if let Some((space, first, second)) = self.compact_dot(function) {
                    self.dot(space, first, second, exponent);
                } else {
                    self.metric(self.argument(first), self.argument(second), exponent)?;
                }
            }
            AtomView::Fun(function) => {
                let Some(head) = self.vector(function).filter(|_| self.rank_one) else {
                    return (exponent == 1)
                        .then_some(())
                        .and_then(|()| self.tensor(function));
                };
                let argument = self.argument(function.iter().last().unwrap());
                let Argument::Original(argument) = argument else {
                    let Argument::Vector(space, other) = argument else {
                        unreachable!()
                    };
                    self.dot(space, head, other, exponent);
                    return Some(());
                };
                let (space, index) = self.resolve_endpoint(argument)?;
                if exponent > 1 {
                    self.dot(space, head, head, exponent / 2);
                    self.contracted = true;
                }
                if exponent % 2 == 1 {
                    let node = self.endpoint(argument, space, index)?;
                    self.nodes[node].terminals.push(Terminal::Vector(head));
                }
            }
            _ => return None,
        }
        Some(())
    }

    fn reduce_components(&mut self) -> Option<()> {
        for node in &mut self.nodes {
            // Degree-one indices survive as the exact original external slots.
            // A single vector or metric is already canonical and does no work.
            if node.degree == 1 {
                node.terminals.push(Terminal::Free(node.slot));
            }
        }
        for node in 0..self.nodes.len() {
            let root = self.root(node);
            if root != node {
                let terminals = std::mem::take(&mut self.nodes[node].terminals);
                self.nodes[root].terminals.extend(terminals);
                self.nodes[root].metrics += self.nodes[node].metrics;
                if AtomView::cmp(&self.nodes[root].slot, &self.nodes[node].slot).is_lt() {
                    self.nodes[root].slot = self.nodes[node].slot;
                }
            }
        }
        for node in 0..self.nodes.len() {
            if self.nodes[node].parent != node {
                continue;
            }
            let space = self.nodes[node].space;
            let edges = self.nodes[node].metrics;
            self.contracted |= edges > 1;
            match self.nodes[node].terminals.as_slice() {
                [] => {
                    self.contracted = true;
                    self.factor(self.spaces[space].1, 1)?;
                }
                &[Terminal::Vector(first), Terminal::Vector(second)] => {
                    self.contracted = true;
                    self.dot(space, first, second, 1);
                }
                &[Terminal::Vector(head), Terminal::Free(slot)]
                | &[Terminal::Free(slot), Terminal::Vector(head)] => {
                    self.contracted |= edges != 0;
                    self.variable(Variable::Vector(head, slot), 1);
                }
                &[Terminal::Free(first), Terminal::Free(second)] => {
                    let slots = if AtomView::cmp(&first, &second).is_le() {
                        [first, second]
                    } else {
                        [second, first]
                    };
                    self.variable(Variable::Metric(slots), 1);
                }
                &[Terminal::Tensor(tensor, position), Terminal::Vector(head)]
                | &[Terminal::Vector(head), Terminal::Tensor(tensor, position)] => {
                    self.contracted = true;
                    self.tensors[tensor].1[position] = Argument::Vector(space, head);
                }
                &[Terminal::Tensor(tensor, position), Terminal::Free(slot)]
                | &[Terminal::Free(slot), Terminal::Tensor(tensor, position)] => {
                    self.contracted |= edges != 0;
                    self.tensors[tensor].1[position] = Argument::Original(slot);
                }
                &[Terminal::Tensor(first, a), Terminal::Tensor(second, b)] => {
                    self.contracted |= edges != 0;
                    let slot = if self.metrics >= 3 && edges > 1 {
                        // The existing metric-path shortcut removes internal
                        // labels before its residual metric is substituted.
                        let Argument::Original(a) = self.tensors[first].1[a] else {
                            unreachable!()
                        };
                        let Argument::Original(b) = self.tensors[second].1[b] else {
                            unreachable!()
                        };
                        if AtomView::cmp(&a, &b).is_lt() { b } else { a }
                    } else {
                        // Ordered symmetric metrics eliminate their smaller
                        // endpoint while both ends still have a tensor partner.
                        self.nodes[node].slot
                    };
                    self.tensors[first].1[a] = Argument::Original(slot);
                    self.tensors[second].1[b] = Argument::Original(slot);
                }
                _ => return None,
            }
        }
        // Equal keys prove only a renaming of internal labels. Reusing the
        // first complete literal plan avoids fresh names, capture, and changes
        // to free labels or opaque metadata. All input admission still runs.
        let key = self.monomial_key().ok()?;
        if let Some(representative) = key.as_ref().and_then(|key| self.alpha_terms.get(key)) {
            self.factors.clone_from(representative);
            self.tensors.clear();
            self.contracted = true;
            return Some(());
        }
        let mut tensors = std::mem::take(&mut self.tensors);
        for (head, arguments) in tensors.drain(..) {
            self.variable(Variable::Tensor(head, arguments), 1);
        }
        self.tensors = tensors;
        if let Some(key) = key {
            let bytes = key
                .factors
                .capacity()
                .checked_mul(std::mem::size_of::<(usize, u16)>())
                .and_then(|bytes| {
                    bytes.checked_add(
                        key.tensors
                            .capacity()
                            .checked_mul(std::mem::size_of::<usize>())?,
                    )
                })
                .and_then(|bytes| {
                    bytes.checked_add(
                        key.connections
                            .capacity()
                            .checked_mul(std::mem::size_of::<(usize, [(usize, usize); 2])>())?,
                    )
                })
                .and_then(|bytes| {
                    bytes.checked_add(
                        self.factors
                            .len()
                            .checked_mul(std::mem::size_of::<(usize, u16)>())?,
                    )
                })
                .and_then(|bytes| {
                    bytes.checked_add(4 * std::mem::size_of::<(MonomialKey, Vec<(usize, u16)>)>())
                });
            self.reserve_alpha(bytes?).ok()?;
            self.alpha_terms.insert(key, self.factors.clone());
        }
        Some(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::shorthands::schoonschip::Schoonschip;
    use spenso::network::tags::SPENSO_TAG;
    use std::sync::{Arc, Mutex};
    use symbolica::parser::ParseSettings;

    fn parse(value: &str) -> Atom {
        Atom::parse(value, "closed_component_test", ParseSettings::symbolica()).unwrap()
    }
    #[test]
    fn metric_only_disconnected_components_keep_separate_dummy_representatives() {
        setup();
        let source = parse(
            "p(spenso::mink(D,mu))*q(spenso::mink(D,mu))*p(spenso::mink(D,nu))*q(spenso::mink(D,nu))",
        );
        let metric_source = parse(
            "spenso::g(spenso::mink(D,mu),spenso::mink(D,a))*p(spenso::mink(D,a))*q(spenso::mink(D,mu))*spenso::g(spenso::mink(D,nu),spenso::mink(D,b))*p(spenso::mink(D,b))*q(spenso::mink(D,nu))",
        );
        let spectator = parse("(x+y)^30");
        let dot_squared = parse("spenso::dot(p(spenso::mink(D)),q(spenso::mink(D)))^2");
        for (position, input) in [source, metric_source].into_iter().enumerate() {
            for power in [1, 2, -1] {
                // A scalar sum keeps its internal pair scope under Atom powers.
                // A raw product power would first distribute into vector powers.
                let base = if power == 1 {
                    input.clone()
                } else {
                    &input + Atom::one()
                };
                let expected = if power == 1 {
                    dot_squared.clone()
                } else {
                    &dot_squared + Atom::one()
                };
                let expression = &spectator * base.pow(power);
                let restricted = SymbolicTensor::infer(expression.clone())
                    .unwrap()
                    .contract(
                        crate::tensor::ContractionSettings::default().without_rank_one_tensors(),
                    )
                    .unwrap()
                    .resolved()
                    .unwrap();
                if position == 0 {
                    assert_eq!(restricted.expression(), &expression);
                }
                assert!(restricted.expression().contains(&spectator));
                let complete = restricted
                    .contract(Default::default())
                    .unwrap()
                    .resolved()
                    .unwrap()
                    .to_dots()
                    .unwrap();
                assert_eq!(complete.expression(), &(&spectator * expected.pow(power)));
            }
        }
    }

    #[test]
    fn factorized_contraction_preserves_named_and_scoped_index_identity() {
        let contractor = setup();
        let head = symbolica::symbol!(
            "closed_component_test::port",
            tags = [SPENSO_TAG.index.clone()]
        );
        let first = AbstractIndex::Named(head.into(), 90, 0);
        let second = AbstractIndex::Named(head.into(), 91, 0);
        let scope = symbolica::symbol!("closed_component_test::bra");
        let source = input(
            "(p(mink(4,a))*q(mink(4,b))+q(mink(4,a))*p(mink(4,b)))
             *(r(mink(4,a))*s(mink(4,b))+s(mink(4,a))*r(mink(4,b)))",
        );
        let expected = input(
            "2*g(p(mink(4)),r(mink(4)))*g(q(mink(4)),s(mink(4)))
             +2*g(p(mink(4)),s(mink(4)))*g(q(mink(4)),r(mink(4)))",
        );
        for (a, b) in [
            (first, second),
            (first.scoped(scope), second.scoped(scope)),
            (first, first.scoped(scope)),
        ] {
            use spenso::structure::slot::ParseableAind;
            let source = source
                .replace(input("a"))
                .with(a.to_atom())
                .replace(input("b"))
                .with(b.to_atom());
            let result = contractor
                .contract_factorized(source.as_view(), None, true)
                .expect("registered explicit indices admit the same factorized frontier");
            assert!(result.status == ContractionStatus::Complete);
            let mut resolved = symbolica::atom::AliasedAtom::from(result.root);
            for (handle, body) in result.aliases {
                resolved.register_alias(handle, body);
            }
            assert_eq!(resolved.into_inner(), expected);
            let result = SymbolicTensor::infer(source)
                .unwrap()
                .contract(Default::default())
                .unwrap();
            assert!(result.contraction_complete());
            assert_eq!(result.resolved().unwrap().expression(), &expected);
        }
    }

    #[test]
    fn factorized_named_index_admission_retains_wildcard_refusal() {
        let contractor = setup();
        let head = symbolica::symbol!(
            "closed_component_test::wild_port_",
            tags = [SPENSO_TAG.index.clone()]
        );
        let index = symbolica::function!(head, Atom::num(90), Atom::num(0));
        let a = input("a");
        let source =
            input("(p(mink(4,a))+q(mink(4,a)))*r(mink(4,a))").replace_map(|value, _, output| {
                if value == a.as_view() {
                    **output = index.clone();
                }
            });
        assert!(
            contractor
                .contract_factorized(source.as_view(), None, true)
                .is_none()
        );
    }

    #[test]
    fn factorized_named_index_admission_retains_open_and_dimension_refusal() {
        use spenso::structure::slot::ParseableAind;
        let contractor = setup();
        let scope = symbolica::symbol!("closed_component_test::guard_scope");
        let open = AbstractIndex::Open {
            owner: 98300,
            axis: 0,
        };
        for index in [open, open.scoped(scope)] {
            let a = input("a");
            let index = index.to_atom();
            let source = input("(p(mink(4,a))+q(mink(4,a)))*r(mink(4,a))").replace_map(
                |value, _, output| {
                    if value == a.as_view() {
                        **output = index.clone();
                    }
                },
            );
            assert!(
                contractor
                    .contract_factorized(source.as_view(), None, true)
                    .is_none()
            );
        }
        let head = symbolica::symbol!(
            "closed_component_test::dimension_port",
            tags = [SPENSO_TAG.index.clone()]
        );
        let a = input("a");
        let index = AbstractIndex::Named(head.into(), 98301, 0).to_atom();
        let source = input("(p(mink(D+1,a))+q(mink(D+1,a)))*r(mink(D+1,a))").replace_map(
            |value, _, output| {
                if value == a.as_view() {
                    **output = index.clone();
                }
            },
        );
        assert!(
            contractor
                .contract_factorized(source.as_view(), None, true)
                .is_none()
        );
    }

    pub(super) fn setup() -> SlotContraction {
        crate::representations::initialize();
        for head in ["p", "q", "r", "s"] {
            SPENSO_TAG.rank_one_tensor_symbol(&format!("closed_component_test::{head}"));
        }
        for head in ["t", "u"] {
            SPENSO_TAG.tensor_symbol(&format!("closed_component_test::{head}"));
        }
        let _ = symbolica::symbol!("closed_component_test::routing"; Scalar);
        SlotContraction::new()
    }
    pub(super) fn input(value: &str) -> Atom {
        let mut qualified = value.to_owned();
        for head in ["g(", "mink(", "cof("] {
            let mut output = String::with_capacity(qualified.len());
            let mut start = 0;
            for (position, _) in qualified.match_indices(head) {
                output.push_str(&qualified[start..position]);
                if !qualified[..position]
                    .ends_with(|c: char| c.is_alphanumeric() || c == '_' || c == ':')
                {
                    output.push_str("spenso::");
                }
                output.push_str(head);
                start = position + head.len();
            }
            output.push_str(&qualified[start..]);
            qualified = output;
        }
        parse(&qualified)
    }

    #[test]
    fn input_leaf_intrinsic_certificates_preserve_admission_and_results() {
        let contractor = setup();
        for (intake, source, admitted) in [
            (
                Intake::ScalarExpansion,
                "1+(x+y)^3*g(p(mink(4)),q(mink(4)))",
                true,
            ),
            (
                Intake::ScalarExpansion,
                "(x+y)*g(p(mink(4)),q(mink(4)))-x*g(p(mink(4)),q(mink(4)))-y*g(p(mink(4)),q(mink(4)))",
                true,
            ),
            (
                Intake::ScalarExpansion,
                "1+(x+y)*routing(p(mink(4)))",
                false,
            ),
            (
                Intake::ScalarExpansion,
                "1+(x+y)*g(p(mink(4)),q(mink(6)))",
                false,
            ),
            (Intake::ScalarExpansion, "1+(x+y)^(-2)", false),
            (
                Intake::FactoredContraction,
                "(x+y)*(p(mink(4,a))+q(mink(4,a)))*r(mink(4,a))",
                true,
            ),
            (
                Intake::FactoredContraction,
                "routing(g(p(mink(4)),q(mink(4))))*(p(mink(4,a))+q(mink(4,a)))*r(mink(4,a))",
                true,
            ),
            (
                Intake::FactoredContraction,
                "g(mink(4,a),mink(4,b))*(t(mink(4,a))+u(mink(4,a)))",
                true,
            ),
            (
                Intake::FactoredContraction,
                "p(mink(4,a))*(t(mink(4,a),mink(4,b))+u(mink(4,a),mink(4,b)))-p(mink(4,a))*t(mink(4,a),mink(4,b))-p(mink(4,a))*u(mink(4,a),mink(4,b))",
                true,
            ),
            (
                Intake::FactoredContraction,
                "(p(mink(4,a))+q(mink(4,a)))*r(mink(4,a))*t(mink(4,a))",
                false,
            ),
        ] {
            let source = input(source);
            let current = contractor.materialize_test_sum(
                source.as_view(),
                &mut SlotMatcher::default(),
                intake,
            );
            assert_eq!(current.is_some(), admitted, "{source}");
            if let Some(result) = current {
                assert!(InterfaceInference::normalization_is_intrinsic(
                    source.as_view()
                ));
                match intake {
                    Intake::ScalarExpansion => assert_eq!(result, source.expand()),
                    Intake::FactoredContraction => {
                        // Root scalar spectators remain factored in the new
                        // intake; compare ordinary fully expanded results.
                        let expected = source.expand().schoonschip();
                        assert_eq!(result.expand(), expected.expand(), "{source}");
                    }
                    Intake::ExpandedContraction | Intake::Contraction => unreachable!(),
                }
            }
        }
    }

    #[test]
    fn input_leaf_certification_precedes_cancellation_and_spectator_omission() {
        let contractor = setup();
        let calls = Arc::new(Mutex::new(Vec::new()));
        let observed = Arc::clone(&calls);
        let hook = spenso::vector_symbol!(
            "closed_component_test::input_leaf_callback",
            norm = move |value, _| {
                observed.lock().unwrap().push(value.to_owned());
            }
        );
        let hooked = FunctionBuilder::new(hook)
            .add_arg(input("mink(4)"))
            .finish();
        let metadata = FunctionBuilder::new(symbolica::symbol!(
            "closed_component_test::input_leaf_metadata"; Scalar
        ))
        .add_arg(&hooked)
        .finish();
        let core = input("(p(mink(4,a))+q(mink(4,a)))*r(mink(4,a))");
        let x = input("x");
        let y = input("y");
        let cancelled = (&x + &y) * &hooked - &x * &hooked - &y * &hooked;
        let sources = [
            &metadata * &core,
            &cancelled + &core,
            &cancelled + Atom::one(),
            (&x + &y) * &metadata - &x * &metadata - &y * &metadata,
        ];
        for source in sources {
            assert!(!InterfaceInference::normalization_is_intrinsic(
                source.as_view()
            ));
            for intake in [Intake::FactoredContraction, Intake::ScalarExpansion] {
                calls.lock().unwrap().clear();
                assert!(
                    contractor
                        .materialize_test_sum(
                            source.as_view(),
                            &mut SlotMatcher::default(),
                            intake,
                        )
                        .is_none(),
                    "{source}"
                );
                assert!(calls.lock().unwrap().is_empty(), "{source}");
            }
        }
        assert_eq!(cancelled.expand(), Atom::Zero);
    }

    #[test]
    fn input_leaf_callback_decline_keeps_metric_rank_loss_validation() {
        let contractor = setup();
        let a = input("mink(4,76401)");
        let b = input("mink(4,76403)");
        let calls = Arc::new(Mutex::new(Vec::new()));
        let observed = Arc::clone(&calls);
        let target = b.clone();
        let head = spenso::tensor_symbol!(
            "closed_component_test::input_leaf_rank_loss",
            norm = move |value, output| {
                observed.lock().unwrap().push(value.to_owned());
                if let AtomView::Fun(function) = value
                    && function.iter().any(|argument| argument == target.as_view())
                {
                    **output = Atom::one();
                }
            }
        );
        let source = FunctionBuilder::new(contractor.metric)
            .add_arg(&a)
            .add_arg(&b)
            .finish()
            * FunctionBuilder::new(head).add_arg(&a).finish();
        let typed = SymbolicTensor::<PartialStructure>::infer(source.clone()).unwrap();
        calls.lock().unwrap().clear();
        assert!(
            contractor
                .materialize_test_sum(
                    source.as_view(),
                    &mut SlotMatcher::default(),
                    Intake::FactoredContraction,
                )
                .is_none()
        );
        assert!(calls.lock().unwrap().is_empty());
        let rewritten = source.schoonschip();
        assert_eq!(rewritten, Atom::one());
        assert!(!calls.lock().unwrap().is_empty());
        assert!(typed.with_rewritten_expression(rewritten).is_err());
        let zero = typed.with_rewritten_expression(Atom::Zero).unwrap();
        assert_eq!(zero.structure, typed.structure);
        assert_eq!(typed.expression, source);
    }

    #[test]
    fn scalar_expansion_collects_without_contracting_or_rebuilding_functions() {
        let contractor = setup();
        for source in [
            "1/3*(g(p(mink(4)),q(mink(4)))+x)^3+2*(x+y)*(x-y)",
            "1+(x+y)^2*(x-y)^2",
            "g(p(mink(4)),q(mink(4)))*(x+y)-x*g(p(mink(4)),q(mink(4)))-y*g(p(mink(4)),q(mink(4)))",
            "2*(x+y)-2*x-2*y",
        ] {
            let source = input(source);
            assert!(matches!(source.as_view(), AtomView::Add(_)), "{source}");
            assert!(!source.is_expanded::<Atom>(None), "{source}");
            let result = contractor
                .materialize_test_sum(
                    source.as_view(),
                    &mut SlotMatcher::default(),
                    Intake::ScalarExpansion,
                )
                .expect("supported scalar sum must use the new intake");
            assert_eq!(result, source.expand());
            assert_eq!(result.expand(), result);
            assert_eq!(
                source.schoonschip(),
                source,
                "scalar Schoonschip stays a no-op"
            );
        }
    }

    #[test]
    fn scalar_expansion_checks_unsupported_leaves_before_cancellation() {
        let contractor = setup();
        for source in [
            "1+(x+y)*g(mink(4,a),mink(4,b))*t(mink(4,a))*u(mink(4,b))",
            "(x+y)*t(mink(4,a))^3-x*t(mink(4,a))^3-y*t(mink(4,a))^3",
            "(x+y)*t(p(mink(4)))-x*t(p(mink(4)))-y*t(p(mink(4)))",
            "1+(x+y)*p(mink(4))",
            "1+(x+y)*routing(p(mink(4)))",
            "1+(x+y)*g(p(mink(4)),q(mink(6)))",
        ] {
            let source = input(source);
            assert!(
                contractor
                    .materialize_test_sum(
                        source.as_view(),
                        &mut SlotMatcher::default(),
                        Intake::ScalarExpansion,
                    )
                    .is_none(),
                "{source}"
            );
        }
        let calls = Arc::new(Mutex::new(Vec::new()));
        let observed = Arc::clone(&calls);
        let hook = symbolica::symbol!(
            "scalar_intake_hook",
            norm = move |value, _| {
                observed.lock().unwrap().push(value.to_owned());
            }
        );
        let hidden = FunctionBuilder::new(symbolica::symbol!("scalar_intake_meta"; Scalar))
            .add_arg(FunctionBuilder::new(hook).add_arg(1).finish())
            .finish();
        let source = input("x+y") * &hidden - input("x") * &hidden - input("y") * &hidden;
        assert!(matches!(source.as_view(), AtomView::Add(_)));
        calls.lock().unwrap().clear();
        assert!(
            contractor
                .materialize_test_sum(
                    source.as_view(),
                    &mut SlotMatcher::default(),
                    Intake::ScalarExpansion,
                )
                .is_none()
        );
        assert!(calls.lock().unwrap().is_empty());
        assert_eq!(source.expand(), Atom::Zero);
    }

    #[test]
    fn scalar_expansion_no_work_powers_and_resource_refusal_keep_ordinary_fallback() {
        let contractor = setup();
        for source in [
            "1+x*g(p(mink(4)),q(mink(4)))",
            "1+(x+y)^-1",
            "1+(x+y)^(1/2)",
            "1+(x+y)^z",
            "1+x^65536*(x+y)",
            // Virtual Cartesian expansion would exceed its bound, although
            // ordinary binomial expansion has only 21 terms.
            "1+(x+y)^20",
        ] {
            let source = input(source);
            assert!(
                contractor
                    .materialize_test_sum(
                        source.as_view(),
                        &mut SlotMatcher::default(),
                        Intake::ScalarExpansion,
                    )
                    .is_none(),
                "{source}"
            );
            let typed = SymbolicTensor::<PartialStructure>::checked_parts(
                source.clone(),
                PartialStructure::from_logical_slots([]),
            )
            .unwrap();
            let expanded = typed.expanded(None, false).unwrap();
            assert_eq!(expanded.expression, source.expand());
            assert_eq!(expanded.structure, typed.structure);
        }
        let rounded = Atom::num(0.25f64) * input("x+y") + Atom::num(1);
        assert!(
            contractor
                .materialize_test_sum(
                    rounded.as_view(),
                    &mut SlotMatcher::default(),
                    Intake::ScalarExpansion,
                )
                .is_none()
        );
        let typed = SymbolicTensor::<PartialStructure>::checked_parts(
            rounded.clone(),
            PartialStructure::from_logical_slots([]),
        )
        .unwrap();
        assert_eq!(
            typed.expanded(None, false).unwrap().expression,
            rounded.expand()
        );
        assert!(SlotContraction::expand_scalar_sum(input("(x+y)*(x-y)").as_view()).is_none());
    }

    #[test]
    fn component_sum_collects_paths_loops_and_existing_dot_powers() {
        let contractor = setup();
        for (source, expected) in [
            (
                "p(mink(4,0))*q(mink(4,0))+p(mink(6,0))*q(mink(6,0))",
                "g(p(mink(4)),q(mink(4)))+g(p(mink(6)),q(mink(6)))",
            ),
            (
                "p(mink(4,a))*q(mink(4,a))+2*p(mink(4,b))*q(mink(4,b))",
                "3*g(p(mink(4)),q(mink(4)))",
            ),
            (
                "p(mink(4,a))*q(mink(4,a))*g(q(mink(4)),p(mink(4)))-g(p(mink(4)),q(mink(4)))^2",
                "0",
            ),
            (
                "1/3*g(mink(4,a),mink(4,b))*p(mink(4,a))*q(mink(4,b))+1/2*p(mink(4,c))*q(mink(4,c))",
                "5/6*g(p(mink(4)),q(mink(4)))",
            ),
            (
                "g(mink(D,a),mink(D,b))*g(mink(D,b),mink(D,c))*g(mink(D,c),mink(D,a))+7",
                "D+7",
            ),
            (
                "g(mink(D,a),mink(D,b))*g(mink(D,b),mink(D,c))*g(mink(D,c),mink(D,a))*g(mink(D,d),mink(D,e))*g(mink(D,e),mink(D,f))*g(mink(D,f),mink(D,d))+7",
                "D^2+7",
            ),
            (
                "g(mink(0,a),mink(0,b))*g(mink(0,b),mink(0,c))*g(mink(0,c),mink(0,a))+7",
                "7",
            ),
            (
                "p(mink(4,a))*q(mink(4,a))*r(mink(4,b))*s(mink(4,b))+g(p(mink(4)),q(mink(4)))*g(r(mink(4)),s(mink(4)))",
                "2*g(p(mink(4)),q(mink(4)))*g(r(mink(4)),s(mink(4)))",
            ),
        ] {
            let source = input(source);
            let result = contractor
                .materialize_test_sum(
                    source.as_view(),
                    &mut SlotMatcher::default(),
                    Intake::ExpandedContraction,
                )
                .expect("closed sum admitted");
            assert_eq!(result, input(expected), "{source}");
            assert_eq!(source.schoonschip(), result);
            assert_eq!(result.schoonschip(), result);
        }
    }

    #[test]
    fn component_sum_compact_spaces_keep_complete_representation_keys() {
        let contractor = setup();
        for (source, expected) in [
            (
                "p(mink(4,a))*q(mink(4,a))*g(p(mink(D)),q(mink(D)))+p(mink(6,a))*q(mink(6,a))*g(p(mink(D)),q(mink(D)))",
                "g(p(mink(4)),q(mink(4)))*g(p(mink(D)),q(mink(D)))+g(p(mink(6)),q(mink(6)))*g(p(mink(D)),q(mink(D)))",
            ),
            (
                "p(mink(4,a))*t(q(mink(6)),mink(4,a))+p(mink(6,a))*t(q(mink(4)),mink(6,a))",
                "t(q(mink(6)),p(mink(4)))+t(q(mink(4)),p(mink(6)))",
            ),
        ] {
            let source = input(source);
            let result = contractor
                .materialize_test_sum(
                    source.as_view(),
                    &mut SlotMatcher::default(),
                    Intake::ExpandedContraction,
                )
                .expect("repeated compact spaces admitted");
            assert_eq!(result, input(expected));
            assert_eq!(source.schoonschip(), result);
            assert_eq!(result.schoonschip(), result);
        }
        // A successful compact lookup cannot admit another dimension, indexed
        // or malformed syntax, or unsupported variance sharing the same head.
        for compact in [
            "mink(6)",
            "mink(4,a)",
            "mink(4,a,b)",
            "mink(D_)",
            "mink(D+1)",
            "cof(3)",
        ] {
            let source = input(&format!(
                "p(mink(4,a))*q(mink(4,a))*g(p(mink(4)),q({compact}))+1"
            ));
            assert!(
                contractor
                    .materialize_test_sum(
                        source.as_view(),
                        &mut SlotMatcher::default(),
                        Intake::ExpandedContraction
                    )
                    .is_none(),
                "{source}"
            );
        }
    }

    #[test]
    fn component_sum_rejects_uncontracted_ambiguous_opaque_wildcard_and_inexact_inputs() {
        let contractor = setup();
        for source in [
            "1+p(mink(4,a))*q(mink(4,b))",
            "1+p(mink(4,a))*q(mink(6,a))",
            "1+p(mink(4,a))*q(mink(4,a))*r(mink(4,a))",
            "1+p(cof(3,a))*q(cof(3,a))",
            "1+p(17,mink(4,a))*q(mink(4,a))",
            "1+f(p(mink(4,a))*q(mink(4,a)))",
            "1+p(mink(D_,a))*q(mink(D_,a))",
            "1+p(mink(4,a_))*q(mink(4,a_))",
            "1+p(mink(4,spenso::cind(0)))*q(mink(4,spenso::cind(0)))",
            "1+p(mink(4,a))*q(mink(4,a))*g(p(mink(4)),q(mink(4)))^-1",
            "1+p(mink(4,a))*q(mink(4,a))*g(p(mink(4)),q(mink(4)))^(1/2)",
            "1+p(mink(4,a))*q(mink(4,a))*g(p(mink(4)),q(mink(4)))^65535",
        ] {
            let source = input(source);
            assert!(
                contractor
                    .materialize_test_sum(
                        source.as_view(),
                        &mut SlotMatcher::default(),
                        Intake::ExpandedContraction
                    )
                    .is_none(),
                "{source}"
            );
        }
        let conflicting =
            symbolica::atom::SymbolBuilder::new(symbolica::wrap_symbol!("closed_sum_rep_vector"))
                .with_tags([
                    &SPENSO_TAG.rank1,
                    &SPENSO_TAG.tensor,
                    &SPENSO_TAG.representation,
                ])
                .build()
                .unwrap();
        let conflicting = FunctionBuilder::new(conflicting)
            .add_arg(input("mink(4,a)"))
            .finish()
            * input("q(mink(4,a))")
            + Atom::num(1);
        assert!(
            contractor
                .materialize_test_sum(
                    conflicting.as_view(),
                    &mut SlotMatcher::default(),
                    Intake::ExpandedContraction
                )
                .is_none()
        );
        let rounded = Atom::num(0.25f64) * input("p(mink(4,a))*q(mink(4,a))") + Atom::num(1);
        assert!(
            contractor
                .materialize_test_sum(
                    rounded.as_view(),
                    &mut SlotMatcher::default(),
                    Intake::ExpandedContraction
                )
                .is_none()
        );
    }

    #[test]
    fn component_sum_positive_powers_preserve_local_pairing() {
        let contractor = setup();
        for (source, expected) in [
            ("1+p(mink(4,a))^2", "1+g(p(mink(4)),p(mink(4)))"),
            (
                "1+p(mink(4,a))^2*q(mink(4,a))",
                "1+g(p(mink(4)),p(mink(4)))*q(mink(4,a))",
            ),
            (
                "1+g(mink(4,a),mink(4,b))^2*p(mink(4,a))*q(mink(4,b))",
                "1+4*p(mink(4,a))*q(mink(4,b))",
            ),
            ("1+p(mink(4,a))^4", "1+g(p(mink(4)),p(mink(4)))^2"),
            (
                "1+p(mink(4,a))^3*q(mink(4,a))",
                "1+g(p(mink(4)),p(mink(4)))*g(p(mink(4)),q(mink(4)))",
            ),
            (
                "1+p(mink(4,a))^2*q(mink(4,a))*r(mink(4,a))",
                "1+g(p(mink(4)),p(mink(4)))*g(q(mink(4)),r(mink(4)))",
            ),
            ("1+g(mink(D,a),mink(D,b))^2", "1+D"),
            ("1+g(mink(D,a),mink(D,b))^4", "1+D^2"),
            (
                "1+g(mink(D,a),mink(D,b))^3*p(mink(D,a))*q(mink(D,b))",
                "1+D*g(p(mink(D)),q(mink(D)))",
            ),
            ("1+g(mink(4,a),mink(4,b))^6", "65"),
            (
                "1+g(mink(4,a),mink(4,b))^2*p(mink(4,a))*q(mink(4,a))",
                "1+4*g(p(mink(4)),q(mink(4)))",
            ),
        ] {
            let source = input(source);
            let result = contractor
                .materialize_test_sum(
                    source.as_view(),
                    &mut SlotMatcher::default(),
                    Intake::ExpandedContraction,
                )
                .expect("locally paired closed power");
            assert_eq!(result, input(expected), "{source}");
            assert_eq!(source.schoonschip(), result);
            assert_eq!(result.schoonschip(), result);
        }
        for source in [
            "1+p(mink(4,a))^-2",
            "1+p(mink(4,a))^-3",
            "1+p(mink(4,a))^(1/2)",
            "1+g(mink(4,a),mink(4,b))^-2",
        ] {
            let source = input(source);
            assert!(
                contractor
                    .materialize_test_sum(
                        source.as_view(),
                        &mut SlotMatcher::default(),
                        Intake::ExpandedContraction
                    )
                    .is_none(),
                "{source}"
            );
        }
        let settings =
            crate::shorthands::schoonschip::SchoonschipSettings::default().without_rank1_tensors();
        assert_eq!(
            input("1+p(mink(4,a))^3*q(mink(4,a))").schoonschip_with_settings(&settings),
            input("1+g(p(mink(4)),p(mink(4)))*p(mink(4,a))*q(mink(4,a))")
        );
        assert_eq!(
            input("1+p(mink(4,a))^2").schoonschip_with_settings(&settings),
            input("1+g(p(mink(4)),p(mink(4)))")
        );
    }

    #[test]
    fn component_sum_preserves_free_terminals_and_collects_equal_paths() {
        let contractor = setup();
        for (source, expected) in [
            (
                "g(mink(4,a),mink(4,b))*p(mink(4,a))+p(mink(4,b))",
                "2*p(mink(4,b))",
            ),
            (
                "g(mink(4,a),mink(4,b))*g(mink(4,b),mink(4,c))+g(mink(4,a),mink(4,c))",
                "2*g(mink(4,a),mink(4,c))",
            ),
            (
                "g(mink(4,a),mink(4,b))*g(mink(4,b),mink(4,c))*p(mink(4,a))+p(mink(4,c))",
                "2*p(mink(4,c))",
            ),
            ("g(mink(4,a),mink(4,b))*p(mink(4,a))-p(mink(4,b))", "0"),
            (
                "g(mink(4,a),mink(4,b))*p(mink(4,a))*q(mink(4,c))+p(mink(4,b))*q(mink(4,c))",
                "2*p(mink(4,b))*q(mink(4,c))",
            ),
            (
                "g(mink(4,a),mink(4,b))*p(mink(4,a))*q(mink(6,a))+p(mink(4,b))*q(mink(6,a))",
                "2*p(mink(4,b))*q(mink(6,a))",
            ),
            (
                "g(mink(4,a),mink(4,b))*p(mink(4,a))*g(mink(6,a),mink(6,c))*q(mink(6,a))+p(mink(4,b))*q(mink(6,c))",
                "2*p(mink(4,b))*q(mink(6,c))",
            ),
            (
                "g(mink(4,0),mink(4,1))*p(mink(4,0))+p(mink(4,1))",
                "2*p(mink(4,1))",
            ),
            (
                "g(mink(D,a),mink(D,b))*g(mink(D,b),mink(D,c))*g(mink(D,c),mink(D,a))*p(mink(D,d))+p(mink(D,d))",
                "D*p(mink(D,d))+p(mink(D,d))",
            ),
            (
                "g(mink(4,a),mink(4,b))*p(mink(4,a))+q(mink(4,c))",
                "p(mink(4,b))+q(mink(4,c))",
            ),
            ("g(mink(4,a),mink(4,b))*p(mink(4,a))+1", "p(mink(4,b))+1"),
            (
                "g(mink(4,a),mink(4,b))*p(mink(4,a))+g(mink(4,c),mink(4,d))",
                "p(mink(4,b))+g(mink(4,c),mink(4,d))",
            ),
        ] {
            let source = input(source);
            let expected = input(expected);
            let result = contractor
                .materialize_test_sum(
                    source.as_view(),
                    &mut SlotMatcher::default(),
                    Intake::ExpandedContraction,
                )
                .expect("metric/vector components admitted");
            assert_eq!(result, expected, "{source}");
            assert_eq!(source.schoonschip(), expected);
            assert_eq!(result.schoonschip(), result);
        }
    }

    #[test]
    fn component_sum_preserves_open_local_power_residuals() {
        let contractor = setup();
        for (source, expected) in [
            (
                "1+p(mink(4,a))^3",
                "1+g(p(mink(4)),p(mink(4)))*p(mink(4,a))",
            ),
            (
                "g(mink(4,a),mink(4,b))*p(mink(4,a))^3+p(mink(4,b))",
                "g(p(mink(4)),p(mink(4)))*p(mink(4,b))+p(mink(4,b))",
            ),
            (
                "g(mink(D,a),mink(D,b))^3*p(mink(D,a))+p(mink(D,b))",
                "D*p(mink(D,b))+p(mink(D,b))",
            ),
            (
                "g(mink(D,a),mink(D,b))^2*p(mink(D,a))*q(mink(D,b))+p(mink(D,a))*q(mink(D,b))",
                "D*p(mink(D,a))*q(mink(D,b))+p(mink(D,a))*q(mink(D,b))",
            ),
            (
                "g(mink(4,a),mink(4,b))^2*p(mink(4,a))^3+4*p(mink(4,a))",
                "4*g(p(mink(4)),p(mink(4)))*p(mink(4,a))+4*p(mink(4,a))",
            ),
            (
                "g(mink(4,a),mink(4,b))*p(mink(4,a))*p(mink(4,b))+g(p(mink(4)),p(mink(4)))",
                "2*g(p(mink(4)),p(mink(4)))",
            ),
        ] {
            let source = input(source);
            let result = contractor
                .materialize_test_sum(
                    source.as_view(),
                    &mut SlotMatcher::default(),
                    Intake::ExpandedContraction,
                )
                .expect("locally paired powers with free residuals");
            assert_eq!(result, input(expected), "{source}");
            assert_eq!(source.schoonschip(), result);
            assert_eq!(result.schoonschip(), result);
        }
        // Missing work, occurrence-local ports, or multiplicity stay on the old route.
        for source in [
            "p(mink(4,a))+q(mink(4,a))",
            "g(mink(4,a),mink(4,b))+g(mink(4,a),mink(4,c))",
            "g(mink(4,a),mink(4,b))*p(mink(4,a))*q(mink(4,a))+p(mink(4,b))",
            "g(mink(4,a),mink(4,b))*p(mink(4,a))+q(mink(4))",
            "g(mink(4,a),mink(4,b))*p(mink(4,a))+q(mink(4,spenso::open(2,0)))",
            "g(mink(4,a),mink(4,b))*p(mink(4,a))+q(mink(4,spenso::cind(0)))",
        ] {
            let source = input(source);
            assert!(
                contractor
                    .materialize_test_sum(
                        source.as_view(),
                        &mut SlotMatcher::default(),
                        Intake::ExpandedContraction
                    )
                    .is_none(),
                "{source}"
            );
        }
    }

    #[test]
    fn component_sum_preserves_typed_order_and_zero() {
        use crate::tensor::SymbolicTensor;
        use spenso::structure::partial::{PartialStructure, PartialStructureExt};
        let contractor = setup();
        for source in [
            "g(mink(4,a),mink(4,b))*p(mink(4,a))+p(mink(4,b))",
            "g(mink(4,a),mink(4,b))*p(mink(4,a))*q(mink(6,c))+p(mink(4,b))*q(mink(6,c))",
            "g(mink(4,a),mink(4,b))*p(mink(4,a))-p(mink(4,b))",
        ] {
            let source = input(source);
            let mut typed = SymbolicTensor::<PartialStructure>::infer(source.clone()).unwrap();
            typed.structure = PartialStructure::from_logical_slots(
                typed.structure.logical_slots().into_iter().rev(),
            );
            let result = contractor
                .materialize_test_sum(
                    source.as_view(),
                    &mut SlotMatcher::default(),
                    Intake::ExpandedContraction,
                )
                .unwrap();
            SymbolicTensor::<PartialStructure>::validate_interface(&result, &typed.structure)
                .unwrap();
            let result = typed.with_rewritten_expression(result).unwrap();
            assert_eq!(result.structure, typed.structure);
            assert_eq!(result.expression, source.schoonschip());
        }
    }

    #[test]
    fn component_sum_callback_fallback_does_not_speculate_with_builders() {
        let contractor = setup();
        let calls = Arc::new(Mutex::new(Vec::new()));
        let observed = Arc::clone(&calls);
        let head = spenso::vector_symbol!(
            "closed_sum_custom",
            norm = move |value, _| {
                observed.lock().unwrap().push(value.to_owned());
            }
        );
        let custom = FunctionBuilder::new(head)
            .add_arg(input("mink(4,a)"))
            .finish();
        let source = custom * input("q(mink(4,a))") + input("p(mink(4,b))*q(mink(4,b))");
        calls.lock().unwrap().clear();
        assert!(
            contractor
                .materialize_test_sum(
                    source.as_view(),
                    &mut SlotMatcher::default(),
                    Intake::ExpandedContraction
                )
                .is_none()
        );
        assert!(calls.lock().unwrap().is_empty());
        // The unsuccessful admission must leave the exact existing callback route.
        let result = source.schoonschip();
        let transcript = calls.lock().unwrap().clone();
        assert!(!transcript.is_empty());
        calls.lock().unwrap().clear();
        let result_again = source.schoonschip();
        assert_eq!(result, result_again);
        assert_eq!(*calls.lock().unwrap(), transcript);
    }

    #[test]
    fn component_sum_binds_opaque_ports_and_retains_metadata() {
        let contractor = setup();
        let metadata = input("routing(1)");
        let AtomView::Fun(metadata) = metadata.as_view() else {
            panic!("metadata fixture must be a function");
        };
        assert!(metadata.get_symbol().is_scalar());
        for (source, expected) in [
            ("1+p(mink(4,a))*t(mink(4,a))", "1+t(p(mink(4)))"),
            (
                "p(mink(4,a))*q(mink(4,b))*t(r(mink(4)),mink(4,a),routing(p(mink(4))-q(mink(4))),mink(4,b))+t(r(mink(4)),p(mink(4)),routing(p(mink(4))-q(mink(4))),q(mink(4)))",
                "2*t(r(mink(4)),p(mink(4)),routing(p(mink(4))-q(mink(4))),q(mink(4)))",
            ),
            (
                "1+p(mink(4,a))*t(routing(mink(4,a)),mink(4,a))",
                "1+t(routing(mink(4,a)),p(mink(4)))",
            ),
            (
                "1+g(mink(4,a),mink(4,b))*p(mink(4,a))*t(mink(4,b),mink(6,c))",
                "1+t(p(mink(4)),mink(6,c))",
            ),
            (
                "1+p(mink(4,a))^3*t(mink(4,a))",
                "1+g(p(mink(4)),p(mink(4)))*t(p(mink(4)))",
            ),
            (
                "1+g(mink(4,a),mink(4,b))^3*t(mink(4,a))*u(mink(4,b))",
                "1+4*t(mink(4,b))*u(mink(4,b))",
            ),
            // Raw Schoonschip preserves differing additive interfaces; the
            // typed constructor continues to reject this source independently.
            (
                "p(mink(4,a))*t(mink(4,a))+u(mink(4,b))",
                "t(p(mink(4)))+u(mink(4,b))",
            ),
            (
                "1+p(mink(4,c))*q(mink(4,c))*t(mink(4,a))*u(mink(4,a))",
                "1+g(p(mink(4)),q(mink(4)))*t(mink(4,a))*u(mink(4,a))",
            ),
        ] {
            let source = input(source);
            let result = contractor
                .materialize_test_sum(
                    source.as_view(),
                    &mut SlotMatcher::default(),
                    Intake::ExpandedContraction,
                )
                .unwrap_or_else(|| panic!("plain opaque ports not admitted: {source}"));
            assert_eq!(result, input(expected), "{source}");
            assert_eq!(result.schoonschip(), result);
        }
        for source in ["1+t(mink(4,a))*u(mink(4,a))", "1+t(mink(4,a),mink(4,a))"] {
            assert!(
                contractor
                    .materialize_test_sum(
                        input(source).as_view(),
                        &mut SlotMatcher::default(),
                        Intake::ExpandedContraction
                    )
                    .is_none(),
                "opaque dummy pairs alone perform no contraction"
            );
        }
    }

    #[test]
    fn component_sum_opaque_metric_paths_are_alpha_equivalent_to_ordered_substitution() {
        use itertools::Itertools;

        let contractor = setup();
        for edges in 1..=4 {
            for indices in (0..=edges).permutations(edges + 1) {
                for same_tensor in [false, true] {
                    for extra in 0..=2 {
                        let mut factors = indices
                            .windows(2)
                            .map(|pair| {
                                input(&format!("g(mink(4,i{}),mink(4,i{}))", pair[0], pair[1]))
                            })
                            .collect::<Vec<_>>();
                        let first = indices[0];
                        let last = indices[edges];
                        factors.push(input(&if same_tensor {
                            format!("t(mink(4,i{first}),mink(4,i{last}))")
                        } else {
                            format!("t(mink(4,i{first}))*u(mink(4,i{last}))")
                        }));
                        for i in 0..extra {
                            factors.push(input(&format!("g(mink(4,j{i}),mink(4,k{i}))")));
                        }
                        let product = Atom::mul_many(factors);
                        // The flat and ordered reducers may retain different
                        // original dummy representatives. Compare their exact
                        // alpha key, including all free ports and metadata.
                        let expected =
                            SlotContraction::run(product.as_view(), false, true, &mut Vec::new())
                                .normalize_dots()
                                + Atom::num(1);
                        let source = product + Atom::num(1);
                        let result = contractor
                            .materialize_test_sum(
                                source.as_view(),
                                &mut SlotMatcher::default(),
                                Intake::ExpandedContraction,
                            )
                            .expect("plain metric path admitted");
                        let difference = (&result - &expected).expand();
                        assert!(
                            difference.is_zero()
                                || contractor
                                    .materialize_test_sum(
                                        difference.as_view(),
                                        &mut SlotMatcher::default(),
                                        Intake::ExpandedContraction,
                                    )
                                    .is_some_and(|value| value.is_zero()),
                            "{source}: {result} !=alpha {expected}"
                        );
                    }
                }
            }
        }
    }

    #[test]
    fn component_sum_opaque_fallback_keeps_scopes_and_callback_rank_changes() {
        let contractor = setup();
        for source in [
            "1+p(mink(4,a))*t(mink(4,a))^2",
            "1+p(mink(4,a))*(t(mink(4,a))+u(mink(4,b)))",
            "1+p(mink(4,a))*t(opaque(g(mink(4,b),mink(4,c))*u(mink(4,b))),mink(4,a))",
            "1+p(mink(4,a))*t(opaque(q(mink(4,b))*r(mink(4,b))),mink(4,a))",
            // Raw Schoonschip can normalize products inside Scalar functions;
            // the collector must defer to that route when metadata needs work.
            "1+p(mink(4,a))*t(routing(g(mink(4,b),mink(4,c))*u(mink(4,b))),mink(4,a))",
            "1+p(mink(4,a))*t(routing(q(mink(4,b))*r(mink(4,b))),mink(4,a))",
            "1+p(mink(4,a))*t(q(r(mink(4))),mink(4,a))",
            "1+p(mink(4,a))*t(q(mink(4,b)),mink(4,a))",
            "1+p(mink(4,a))*t(mink(4,a),mink(4,a))",
            "1+p(mink(4,a))*t(mink(4,a),mink(4))",
            "1+p(mink(4,a))*t(mink(4,a),mink(4,spenso::open(1,0)))",
        ] {
            let source = input(source);
            assert!(
                contractor
                    .materialize_test_sum(
                        source.as_view(),
                        &mut SlotMatcher::default(),
                        Intake::ExpandedContraction
                    )
                    .is_none(),
                "{source}"
            );
        }
        let calls = Arc::new(Mutex::new(Vec::new()));
        let observed = Arc::clone(&calls);
        let compact = input("p(mink(4))");
        let callback = spenso::tensor_symbol!(
            "component_opaque_rank_loss",
            norm = move |value, out| {
                observed.lock().unwrap().push(value.to_owned());
                if let AtomView::Fun(function) = value
                    && function.iter().next() == Some(compact.as_view())
                {
                    **out = Atom::num(1);
                }
            }
        );
        let leaf = FunctionBuilder::new(callback)
            .add_arg(input("mink(4,a)"))
            .finish();
        let source = leaf * input("p(mink(4,a))") + Atom::num(1);
        calls.lock().unwrap().clear();
        assert!(
            contractor
                .materialize_test_sum(
                    source.as_view(),
                    &mut SlotMatcher::default(),
                    Intake::ExpandedContraction
                )
                .is_none()
        );
        assert!(calls.lock().unwrap().is_empty());
        assert_eq!(source.schoonschip(), Atom::num(2));
        assert_eq!(
            calls.lock().unwrap().len(),
            2,
            "discovery and actual replacement retained"
        );
    }

    #[test]
    fn component_sum_opaque_ports_preserve_typed_order_and_zero() {
        use crate::tensor::SymbolicTensor;
        use spenso::structure::partial::{PartialStructure, PartialStructureExt};

        let contractor = setup();
        for source in [
            "p(mink(4,a))*t(mink(4,a),mink(4,b),mink(6,c))+t(p(mink(4)),mink(4,b),mink(6,c))",
            "p(mink(4,a))*t(mink(4,a),mink(4,b),mink(6,c))-t(p(mink(4)),mink(4,b),mink(6,c))",
        ] {
            let source = input(source);
            let mut tensor = SymbolicTensor::<PartialStructure>::infer(source.clone()).unwrap();
            tensor.structure = PartialStructure::from_logical_slots(
                tensor.structure.logical_slots().into_iter().rev(),
            );
            let expression = contractor
                .materialize_test_sum(
                    source.as_view(),
                    &mut SlotMatcher::default(),
                    Intake::ExpandedContraction,
                )
                .unwrap();
            SymbolicTensor::<PartialStructure>::validate_interface(&expression, &tensor.structure)
                .unwrap();
            let result = tensor.with_rewritten_expression(expression).unwrap();
            assert_eq!(result.structure, tensor.structure);
            assert_eq!(result.expression, source.schoonschip());
        }
        assert!(
            SymbolicTensor::<PartialStructure>::infer(input(
                "p(mink(4,a))*t(mink(4,a))+u(mink(4,b))"
            ))
            .is_err()
        );
    }

    #[test]
    fn component_sum_alpha_collects_connections_without_new_labels() {
        let contractor = setup();
        let _ = symbolica::symbol!("closed_component_test::alpha_metadata"; Scalar);
        for source in [
            "t(mink(4,a))*u(mink(4,a))-t(mink(4,b))*u(mink(4,b))",
            "t(mink(4,a),mink(4,a))-t(mink(4,b),mink(4,b))",
            "t(mink(4,a),mink(4,b))*u(mink(4,a),mink(4,b))-t(mink(4,c),mink(4,d))*u(mink(4,c),mink(4,d))",
            "t(mink(4,a),mink(6,b))*u(mink(4,a),mink(6,b))-t(mink(4,c),mink(6,d))*u(mink(4,c),mink(6,d))",
            "t(mink(4,0))*u(mink(4,0))-t(mink(4,1))*u(mink(4,1))",
            "t(1,mink(4,a))*t(2,mink(4,a))*t(3,mink(4,b),mink(4,b))-t(1,mink(4,c))*t(2,mink(4,c))*t(3,mink(4,d),mink(4,d))",
            "1/3*t(mink(4,a))*u(mink(4,a))-1/3*t(mink(4,b))*u(mink(4,b))",
            // The literal metadata slot is never renamed with the dummy.
            "t(alpha_metadata(mink(4,a)),mink(4,a))*u(mink(4,a))-t(alpha_metadata(mink(4,a)),mink(4,b))*u(mink(4,b))",
        ] {
            let source = input(source);
            assert_ne!(source, Atom::num(0));
            assert_eq!(
                contractor.materialize_test_sum(
                    source.as_view(),
                    &mut SlotMatcher::default(),
                    Intake::ExpandedContraction
                ),
                Some(Atom::num(0)),
                "{source}"
            );
            assert_eq!(
                crate::test_support::contracted_atom(source.as_view()).unwrap(),
                Atom::num(0)
            );
        }
        let source = input("t(mink(4,a))*u(mink(4,a))+t(mink(4,b))*u(mink(4,b))");
        let result = crate::test_support::contracted_atom(source.as_view()).unwrap();
        assert!(
            result == input("2*t(mink(4,a))*u(mink(4,a))")
                || result == input("2*t(mink(4,b))*u(mink(4,b))"),
            "reuse one whole original representative, never a fresh label"
        );
        assert_eq!(
            crate::test_support::contracted_atom(result.as_view()).unwrap(),
            result
        );
    }

    #[test]
    fn component_sum_alpha_preserves_free_metadata_and_connection_identity() {
        let contractor = setup();
        let _ = symbolica::symbol!("closed_component_test::alpha_metadata"; Scalar);
        for source in [
            "t(mink(4,a),mink(4,f))*u(mink(4,a))-t(mink(4,b),mink(4,h))*u(mink(4,b))",
            "t(alpha_metadata(mink(4,a)),mink(4,a))*u(mink(4,a))-t(alpha_metadata(mink(4,b)),mink(4,b))*u(mink(4,b))",
            "t(alpha_metadata(p(mink(4))),mink(4,a))*u(mink(4,a))-t(alpha_metadata(q(mink(4))),mink(4,b))*u(mink(4,b))",
            "t(alpha_metadata(x),mink(4,a))*u(mink(4,a))-t(mink(4,b),alpha_metadata(x))*u(mink(4,b))",
            "t(p(mink(4)),mink(4,a))*u(mink(4,a))-t(q(mink(4)),mink(4,b))*u(mink(4,b))",
            "t(mink(4,a))*u(mink(4,a))-t(mink(6,b))*u(mink(6,b))",
            "t(mink(4,a),mink(4,b))*u(mink(4,a),mink(4,b))-t(mink(4,c),mink(4,d))*u(mink(4,d),mink(4,c))",
            "t(mink(4,a),mink(4,a))*u(mink(4,b),mink(4,b))-t(mink(4,c),mink(4,d))*u(mink(4,c),mink(4,d))",
            // Same factor colors require a general graph-labeling algorithm.
            "t(mink(4,a),mink(4,b))*t(mink(4,c),mink(4,d))*u(mink(4,a),mink(4,c),mink(4,b),mink(4,d))-t(mink(4,e),mink(4,f))*t(mink(4,h),mink(4,i))*u(mink(4,e),mink(4,h),mink(4,f),mink(4,i))",
            // Raw missing-index branches remain separate; no sum shape check.
            "t(mink(4,a))*u(mink(4,a))+t(mink(4,b),mink(4,f))*u(mink(4,b))",
        ] {
            let source = input(source);
            assert!(
                contractor
                    .materialize_test_sum(
                        source.as_view(),
                        &mut SlotMatcher::default(),
                        Intake::ExpandedContraction
                    )
                    .is_none(),
                "no alpha-equivalent pair: {source}"
            );
            assert_eq!(source.schoonschip(), source);
        }
    }

    #[test]
    fn component_sum_alpha_preserves_local_powers_and_compact_binding_aliases() {
        let contractor = setup();
        for source in [
            "g(mink(4,a),mink(4,b))^3*t(mink(4,a))*u(mink(4,b))-4*t(mink(4,c))*u(mink(4,c))",
            "p(mink(4,a))^3*t(mink(4,a),mink(4,b))*u(mink(4,b))-g(p(mink(4)),p(mink(4)))*t(p(mink(4)),mink(4,c))*u(mink(4,c))",
            "g(p(mink(4)),q(mink(4)))^3*t(mink(4,a))*u(mink(4,a))-g(q(mink(4)),p(mink(4)))^3*t(mink(4,b))*u(mink(4,b))",
        ] {
            let source = input(source);
            assert_eq!(
                contractor.materialize_test_sum(
                    source.as_view(),
                    &mut SlotMatcher::default(),
                    Intake::ExpandedContraction
                ),
                Some(Atom::num(0)),
                "{source}"
            );
            assert_eq!(
                crate::test_support::contracted_atom(source.normalize_dots().as_view()).unwrap(),
                Atom::num(0)
            );
        }
        let scalar_power = input(
            "(x+g(p(mink(4)),q(mink(4))))^3*t(mink(4,a))*u(mink(4,a))-(x+g(q(mink(4)),p(mink(4))))^3*t(mink(4,b))*u(mink(4,b))",
        );
        assert_eq!(
            crate::test_support::contracted_atom(scalar_power.as_view()).unwrap(),
            Atom::Zero
        );
        // Derived equal factors still form their old power. Repeated identical
        // tensor identities deliberately do not get an alpha representative.
        let product = input("g(mink(4,a),mink(4,b))*t(mink(4,a))*t(mink(4,b))");
        let expected = SlotContraction::run(product.as_view(), false, true, &mut Vec::new())
            .normalize_dots()
            + Atom::num(1);
        let source = product + Atom::num(1);
        assert_eq!(
            crate::test_support::contracted_atom(source.normalize_dots().as_view()).unwrap(),
            expected
        );
        for source in [
            "t(mink(4,a))^2+t(mink(4,b))^2",
            "g(mink(4,a),mink(4,b))^-2*t(mink(4,a))*u(mink(4,b))+1",
        ] {
            assert!(
                contractor
                    .materialize_test_sum(
                        input(source).as_view(),
                        &mut SlotMatcher::default(),
                        Intake::ExpandedContraction
                    )
                    .is_none()
            );
        }
    }

    #[test]
    fn component_sum_alpha_stage_three_cancellation_retains_typed_free_ports() {
        use crate::tensor::SymbolicTensor;
        use spenso::structure::partial::{PartialStructure, PartialStructureExt};

        let contractor = setup();
        let _ = symbolica::symbol!("closed_component_test::alpha_metadata"; Scalar);
        // Portable form of the measured reverse-stage-3 mu6/mu11 mismatch:
        // same two vertices, ordered routing metadata, and free mu2/mu7 ports.
        for name in ["mu5", "mu6", "mu11"] {
            let _ = parse(name);
        }
        let vertex = |middle: &str| {
            format!(
                "t(3,alpha_metadata(-p(mink(4))),alpha_metadata(q(mink(4))),alpha_metadata(p(mink(4))-q(mink(4))),mink(4,mu2),mink(4,{middle}),mink(4,mu10))*t(7,alpha_metadata(r(mink(4))-q(mink(4))),alpha_metadata(-p(mink(4))+q(mink(4))),alpha_metadata(-r(mink(4))+p(mink(4))),mink(4,mu6),mink(4,mu10),mink(4,mu7))"
            )
        };
        let first = input(&format!("g(mink(4,mu5),mink(4,mu6))*{}", vertex("mu5")));
        let second = input(&format!("g(mink(4,mu6),mink(4,mu11))*{}", vertex("mu11")));
        let old_first =
            SlotContraction::run(first.as_view(), false, true, &mut Vec::new()).normalize_dots();
        let old_second =
            SlotContraction::run(second.as_view(), false, true, &mut Vec::new()).normalize_dots();
        assert_ne!(
            old_first, old_second,
            "legacy results intentionally retain different dummy labels"
        );
        let source = first - second;
        let mut tensor = SymbolicTensor::<PartialStructure>::infer(source.clone()).unwrap();
        assert_eq!(tensor.structure.logical_slots().len(), 2);
        tensor.structure = PartialStructure::from_logical_slots(
            tensor.structure.logical_slots().into_iter().rev(),
        );
        let expression = contractor
            .materialize_test_sum(
                source.as_view(),
                &mut SlotMatcher::default(),
                Intake::ExpandedContraction,
            )
            .unwrap();
        assert_eq!(expression, Atom::num(0));
        let result = tensor.with_rewritten_expression(expression).unwrap();
        assert_eq!(result.structure, tensor.structure);
        assert_eq!(
            crate::test_support::contracted_atom(source.as_view()).unwrap(),
            Atom::num(0)
        );
    }

    #[test]
    fn component_sum_alpha_budget_retains_the_bounded_frontier_and_exact_reference() {
        use crate::shorthands::schoonschip::SchoonschipSettings;

        let contractor = setup();
        assert_eq!(ComponentSum::MAX_ALPHA_BYTES, 16 * 1024);
        // These groups cancel only after the ordinary metric contraction. Their
        // distinct tensor metadata exceeds the small compile-time test budget,
        // even though the final fallback expression has only two small terms.
        let a = input("mink(4,budget_a)");
        let b = input("mink(4,budget_b)");
        assert!(AtomView::cmp(&a.as_view(), &b.as_view()).is_lt());
        let survivor = input("t(999,mink(4,c))*u(mink(4,c))-t(999,mink(4,d))*u(mink(4,d))");
        let mut terms = vec![survivor.clone()];
        for label in 0..64 {
            terms.push(input(&format!(
                "g(mink(4,budget_a),mink(4,budget_b))*t({label},mink(4,budget_a))*u(mink(4,budget_b))"
            )));
            terms.push(-input(&format!(
                "t({label},mink(4,budget_b))*u(mink(4,budget_b))"
            )));
        }
        let source = Atom::add_many(&terms);
        assert!(
            contractor
                .materialize_test_sum(
                    source.as_view(),
                    &mut SlotMatcher::default(),
                    Intake::ExpandedContraction
                )
                .is_none(),
            "resource exhaustion must discard all partial alpha collection"
        );
        let fallback = Atom::add_many(
            terms
                .iter()
                .map(|term| {
                    SlotContraction::run(term.as_view(), false, false, &mut Vec::new())
                        .normalize_dots()
                })
                .collect::<Vec<_>>(),
        );
        assert_eq!(
            fallback, survivor,
            "ordinary contractions actually shrink the declined sum"
        );
        assert_ne!(fallback, source);
        assert_ne!(fallback, Atom::num(0));
        // The public bounded frontier retains the exact source rather than
        // restarting the declined work through an unbounded fallback.
        let value = SymbolicTensor::infer(source.clone()).unwrap();
        let pending = Arc::new(value.contract(Default::default()).unwrap());
        assert!(!pending.contraction_complete());
        assert_eq!(pending.resolved().unwrap().expression, source);
        let again = pending.contract(Default::default()).unwrap();
        assert!(!again.contraction_complete());
        assert_eq!(again.resolved().unwrap(), pending.resolved().unwrap());
        // Keep the exact zero oracle on the existing explicit per-term
        // reference computed above; its mathematical expectation is unchanged.
        assert_eq!(
            crate::test_support::contracted_atom(fallback.as_view()).unwrap(),
            Atom::num(0)
        );
        assert_eq!(
            crate::test_support::contracted_atom(
                crate::test_support::contracted_atom(fallback.as_view())
                    .unwrap()
                    .as_view()
            )
            .unwrap(),
            Atom::num(0)
        );
        assert_eq!(
            source
                .schoonschip_with_settings(&SchoonschipSettings::default().without_rank1_tensors()),
            fallback,
            "disabled rank-one mode keeps its existing route"
        );
    }

    #[test]
    fn component_sum_alpha_finishes_metadata_admitted_after_dot_cleanup() {
        let contractor = setup();
        let source = input(
            "t(p(mink(4,x))*q(mink(4,x)),mink(4,a))*u(mink(4,a))-t(g(p(mink(4)),q(mink(4))),mink(4,b))*u(mink(4,b))",
        );
        assert!(
            contractor
                .materialize_test_sum(
                    source.as_view(),
                    &mut SlotMatcher::default(),
                    Intake::ExpandedContraction
                )
                .is_none(),
            "pending work in a tensor argument must retain the initial fallback"
        );
        let AtomView::Add(terms) = source.as_view() else {
            panic!("two metadata alternatives")
        };
        let after_contraction = Atom::add_many(
            terms
                .iter()
                .map(|term| SlotContraction::run(term, false, true, &mut Vec::new()))
                .collect::<Vec<_>>(),
        );
        let after_dots = after_contraction.normalize_dots();
        assert_ne!(
            after_contraction, after_dots,
            "final dot cleanup changes admission"
        );
        assert_ne!(
            after_dots,
            Atom::num(0),
            "literal dummy labels still differ"
        );
        assert_eq!(
            contractor.materialize_test_sum(
                after_dots.as_view(),
                &mut SlotMatcher::default(),
                Intake::ExpandedContraction
            ),
            Some(Atom::num(0))
        );
        assert_eq!(
            crate::test_support::contracted_atom(source.as_view()).unwrap(),
            Atom::num(0)
        );
        assert_eq!(
            crate::test_support::contracted_atom(
                crate::test_support::contracted_atom(source.as_view())
                    .unwrap()
                    .as_view()
            )
            .unwrap(),
            Atom::num(0)
        );
    }

    #[test]
    fn fused_expansion_contracts_opaque_ports_and_coalesced_powers() {
        let contractor = setup();
        for source in [
            "(p(mink(4,a))+q(mink(4,a)))*t(routing(p(mink(4))-q(mink(4))),mink(4,a),mink(4,b))*u(mink(4,b))",
            "(g(mink(4,a),mink(4,c))*p(mink(4,c))+g(mink(4,a),mink(4,d))*q(mink(4,d)))*t(mink(4,a),mink(4,b))*u(mink(4,b))",
            "(p(mink(4,a))+q(mink(4,a)))*p(mink(4,a))^2*r(mink(4,a))",
            "(p(mink(4,a))+q(mink(4,a)))^2",
            "(g(mink(4,a),mink(4,b))+g(mink(4,a),mink(4,c)))*g(mink(4,a),mink(4,b))^2*p(mink(4,a))*q(mink(4,b))",
            "(p(mink(4,a))+1)*t(mink(4,a))",
            "(g(p(mink(4)),q(mink(4)))-g(q(mink(4)),q(mink(4))))*g(mink(4,a),mink(4,b))*t(mink(4,a))*u(mink(4,b))+(p(mink(4,a))-q(mink(4,a)))*t(mink(4,a))",
        ] {
            let source = input(source);
            let expected = source.expand().schoonschip();
            let actual = contractor
                .materialize_test_sum(
                    source.as_view(),
                    &mut SlotMatcher::default(),
                    Intake::FactoredContraction,
                )
                .expect("the fused intake must admit the opaque-vertex grammar");
            assert_eq!(actual.expand(), expected.expand(), "{source}");
            assert_eq!(actual.schoonschip(), actual, "{source}");
        }
    }

    #[test]
    fn fused_expansion_keeps_common_scalar_spectators_factored() {
        let contractor = setup();
        let core = input("(p(mink(4,a))+q(mink(4,a)))*(r(mink(4,a))+s(mink(4,a)))");
        let contracted = core.expand().schoonschip();
        for spectator in [
            "(x+y)^8",
            "g(p(mink(4)),q(mink(4)))^3",
            "(g(p(mink(4)),q(mink(4)))+x)^3",
            "routing(mink(4,a))",
        ] {
            let spectator = input(spectator);
            let source = &spectator * &core;
            let actual = contractor
                .materialize_test_sum(
                    source.as_view(),
                    &mut SlotMatcher::default(),
                    Intake::FactoredContraction,
                )
                .unwrap();
            assert_eq!(actual, &spectator * &contracted, "{source}");
            assert_eq!(actual.schoonschip(), actual, "{source}");
        }
    }

    #[test]
    fn fused_expansion_preserves_zero_and_exact_branch_coefficients() {
        let contractor = setup();
        for (source, expected) in [
            (
                "(p(mink(4,a))+q(mink(4,a)))*t(mink(4,a))-t(p(mink(4)))-t(q(mink(4)))",
                "0",
            ),
            (
                "(p(mink(4,a))+2*q(mink(4,a)))*(p(mink(4,a))-2*q(mink(4,a)))",
                "g(p(mink(4)),p(mink(4)))-4*g(q(mink(4)),q(mink(4)))",
            ),
            (
                "(p(mink(4,a))/3+q(mink(4,a))/2)*t(mink(4,a))",
                "t(p(mink(4)))/3+t(q(mink(4)))/2",
            ),
        ] {
            let source = input(source);
            assert_eq!(
                contractor
                    .materialize_test_sum(
                        source.as_view(),
                        &mut SlotMatcher::default(),
                        Intake::FactoredContraction
                    )
                    .unwrap(),
                input(expected),
                "{source}"
            );
        }
    }

    #[test]
    fn fused_expansion_declines_whole_attempt_before_callbacks() {
        use std::sync::{
            Arc,
            atomic::{AtomicUsize, Ordering},
        };
        let contractor = setup();
        let calls = Arc::new(AtomicUsize::new(0));
        let observed = Arc::clone(&calls);
        let callback = spenso::tensor_symbol!(
            "closed_component_test::fused_callback",
            norm = move |_, _| {
                observed.fetch_add(1, Ordering::Relaxed);
            }
        );
        let callback = FunctionBuilder::new(callback)
            .add_arg(spenso::mink!(4, 75321))
            .finish();
        let source = input("(p(mink(4,a))+q(mink(4,a)))*t(mink(4,a))") + callback;
        calls.store(0, Ordering::Relaxed);
        assert!(
            contractor
                .materialize_test_sum(
                    source.as_view(),
                    &mut SlotMatcher::default(),
                    Intake::FactoredContraction
                )
                .is_none()
        );
        assert_eq!(calls.load(Ordering::Relaxed), 0);

        for source in [
            "(p(mink(4,a))+q(mink(4,a)))^20",
            "(p(mink(4,a))^65535+q(mink(4,a)))*p(mink(4,a))",
            "(t(mink(4,a))+u(mink(4,a)))^2",
            "(p(mink(4,a))+q(mink(4,a)))*t(routing(p(mink(4,b))*q(mink(4,b))),mink(4,a))",
            "(p(mink(4,a))+q(mink(4,a)))*t(mink(6,a))^2",
            "(p(mink(4,a))+q(mink(4,a)))*r(mink(4,a))*t(mink(4,a))",
            "(p(mink(4,a))+q(mink(4,a)))^(-1)",
        ] {
            let source = input(source);
            assert!(
                contractor
                    .materialize_test_sum(
                        source.as_view(),
                        &mut SlotMatcher::default(),
                        Intake::FactoredContraction
                    )
                    .is_none(),
                "{source}"
            );
        }
    }

    #[test]
    fn fused_input_collection_cancels_before_incidence_and_opaque_power_admission() {
        let contractor = setup();
        for source in [
            "(p(mink(4,a))+q(mink(4,a)))*(r(mink(4,a))+s(mink(4,a)))*t(mink(4,a))-p(mink(4,a))*r(mink(4,a))*t(mink(4,a))-p(mink(4,a))*s(mink(4,a))*t(mink(4,a))-q(mink(4,a))*r(mink(4,a))*t(mink(4,a))-q(mink(4,a))*s(mink(4,a))*t(mink(4,a))+p(mink(4,b))*q(mink(4,b))",
            "(t(mink(4,a))+u(mink(4,a)))^2-t(mink(4,a))^2-2*t(mink(4,a))*u(mink(4,a))-u(mink(4,a))^2+p(mink(4,b))*q(mink(4,b))",
        ] {
            let source = input(source);
            let expected = input("g(p(mink(4)),q(mink(4)))");
            assert_eq!(source.expand().schoonschip(), expected);
            assert_eq!(
                contractor.materialize_test_sum(
                    source.as_view(),
                    &mut SlotMatcher::default(),
                    Intake::FactoredContraction
                ),
                Some(expected),
            );
        }
        // A surviving malformed graph or opaque power must still decline.
        for source in [
            "(p(mink(4,a))+q(mink(4,a)))*r(mink(4,a))*t(mink(4,a))-p(mink(4,a))*r(mink(4,a))*t(mink(4,a))",
            "(t(mink(4,a))+u(mink(4,a)))^2-2*t(mink(4,a))*u(mink(4,a))-u(mink(4,a))^2",
        ] {
            assert!(
                contractor
                    .materialize_test_sum(
                        input(source).as_view(),
                        &mut SlotMatcher::default(),
                        Intake::FactoredContraction
                    )
                    .is_none()
            );
        }
    }

    #[test]
    fn fused_input_collection_numeric_distribution_can_cancel_every_row() {
        let contractor = setup();
        let source = input(
            "2*(t(mink(4,a),mink(4,a))+u(mink(4,a),mink(4,a)))-2*t(mink(4,a),mink(4,a))-2*u(mink(4,a),mink(4,a))",
        );
        assert_ne!(source, Atom::num(0), "the source still needs distribution");
        assert_eq!(source.expand(), Atom::num(0));
        assert_eq!(
            contractor.materialize_test_sum(
                source.as_view(),
                &mut SlotMatcher::default(),
                Intake::FactoredContraction
            ),
            Some(Atom::num(0)),
        );
        assert_eq!(
            crate::test_support::contracted_atom(source.as_view())
                .unwrap()
                .expand(),
            Atom::num(0),
        );

        // This cancellation happens before graph planning, while the typed
        // receiver still owns the surviving logical port b.
        let source = input(
            "p(mink(4,a))*(t(mink(4,a),mink(4,b))+u(mink(4,a),mink(4,b)))-p(mink(4,a))*t(mink(4,a),mink(4,b))-p(mink(4,a))*u(mink(4,a),mink(4,b))",
        );
        assert_ne!(source, Atom::Zero);
        let inferred = SymbolicTensor::<PartialStructure>::infer(source.clone()).unwrap();
        let typed = SymbolicTensor::checked_parts(inferred.expression, inferred.structure).unwrap();
        assert_eq!(typed.expression, source);
        let expected = SymbolicTensor::<PartialStructure>::infer(input("t(mink(4,b))")).unwrap();
        assert_eq!(typed.structure.logical_slots().len(), 1);
        assert_eq!(
            typed.structure.logical_slots(),
            expected.structure.logical_slots()
        );
        let result = typed
            .contract(Default::default())
            .unwrap()
            .expanded()
            .unwrap();
        assert_eq!(result.expression, Atom::Zero);
        assert_eq!(result.structure, typed.structure);
        assert_eq!(typed.expression, source);
    }

    #[test]
    fn fused_input_collection_rolls_back_powers_and_exact_coefficients() {
        let contractor = setup();
        for source in [
            "(p(mink(4,a))/3+q(mink(4,a))/2)^3*(r(mink(4,a))+s(mink(4,a)))",
            "(p(mink(4,a))+q(mink(4,a)))^2-(p(mink(4,a))-q(mink(4,a)))^2",
            "(p(mink(4,a))+q(mink(4,a)))*t(routing(mink(4,a)),mink(4,a))-(q(mink(4,a))+p(mink(4,a))/2)*t(routing(mink(4,a)),mink(4,a))",
        ] {
            let source = input(source);
            let actual = contractor
                .materialize_test_sum(
                    source.as_view(),
                    &mut SlotMatcher::default(),
                    Intake::FactoredContraction,
                )
                .unwrap();
            assert_eq!(
                actual.expand(),
                source.expand().schoonschip().expand(),
                "{source}"
            );
            assert_eq!(actual.schoonschip(), actual, "{source}");
        }
    }
}
