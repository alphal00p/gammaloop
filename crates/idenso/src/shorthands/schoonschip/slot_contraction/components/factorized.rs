//! State merging on the collector's existing port graph.
//!
//! Remaining factors are opaque terminals of the same component reducer used
//! for a completed monomial. Only the selected factor contributes alternatives;
//! substitutions remain port bindings until a distinct output variable is emitted.

use super::super::SlotContraction;
use super::{
    Argument, ComponentSum, InputLeaf, Intake, ScalarPolynomial, TensorSource, Terminal,
    TerminalPorts, Variable,
};
#[cfg(feature = "reference-cases")]
use crate::reference_cases::timing;
#[cfg(feature = "reference-cases")]
use crate::reference_cases::timing::{Frontier, Phase, ScopeOpenReason};
use crate::shorthands::bracket::BracketNormalizer;
use crate::shorthands::schoonschip::{
    DotNormalizer, SchoonschipSettings, SchoonschipWithSettings, SimplificationCandidates,
};
use crate::tensor::composition;
use crate::tensor::inference::{InterfaceInference, TensorInferenceError};
use crate::tensor::planning::{self, PlanningLeaf, PlanningNetwork, PlanningSource};
use crate::tensor::simplification::observation::DomainObservations;
use crate::tensor::{ContractSettings, SymbolicTensor, TensorProofs};
use ahash::{AHashMap, AHashSet};
use linnet::half_edge::NodeIndex;
use linnet::half_edge::subgraph::{Inclusion, SuBitGraph};
use linnet::half_edge::tree::SimpleTraversalTree;
use linnet::tree::child_vec::ChildVecStore;
use spenso::network::graph::{NetworkEdge, NetworkGraph, NetworkLeaf, NetworkNode, NetworkOp};
use spenso::network::library::DummyKey;
use spenso::network::parsing::{
    AtomStructureExt, LeafInterfaceCache, ParseSettings, ShorthandParsing,
};
use spenso::network::store::{NetworkStore, TensorScalarStore};
use spenso::shadowing::{TermLeaf, TermTape};
use spenso::structure::abstract_index::AbstractIndex;
use spenso::structure::partial::{PartialIndex, PartialStructure, PartialStructureExt};
use spenso::structure::representation::{LibraryRep, RepName};
use spenso::structure::slot::{IsAbstractSlot, SlotMatch, SlotMatcher};
use spenso::structure::{OrderedStructure, TensorStructure};
use std::borrow::Cow;
use std::cell::OnceCell;
use std::cmp::Reverse;
use std::collections::{BTreeSet, HashMap};
use std::mem;
use std::sync::Arc;
use symbolica::atom::{Atom, AtomCore, AtomOrView, AtomView, Symbol};
use symbolica::coefficient::CoefficientView;
use symbolica::domains::rational::{Q, Rational};
use symbolica::domains::{Ring, RingOps};
use symbolica::poly::PolyVariable;

// A state contains occurrence-local port bindings, not rebuilt expressions.
// Variables and the incidence scratch still belong to ComponentSum.
#[derive(Clone, PartialEq, Eq, Hash)]
struct Remaining<'a> {
    factors: Vec<(usize, Vec<Argument<'a>>)>,
    residual: Vec<(usize, u32)>,
}

impl Remaining<'_> {
    #[cfg(feature = "reference-cases")]
    fn storage_bytes(&self) -> usize {
        mem::size_of::<Self>()
            + self.factors.capacity() * mem::size_of::<(usize, Vec<Argument<'_>>)>()
            + self
                .factors
                .iter()
                .map(|(_, arguments)| arguments.capacity() * mem::size_of::<Argument<'_>>())
                .sum::<usize>()
            + self.residual.capacity() * mem::size_of::<(usize, u32)>()
    }
}

// Merge similarly sized coefficients instead of copying a growing polynomial
// for every contribution to one frontier state. There is at most one partial
// per logarithmic term-count bin, hence at most usize::BITS partials. The common
// one-coefficient case needs no additional allocation. Symbolica still owns
// variable-map unification, exact cancellation, and monomial ordering.
enum CoefficientSum {
    One(ScalarPolynomial),
    Binned(Vec<Option<ScalarPolynomial>>),
}

impl CoefficientSum {
    /// Restore the min-heap order of stream `heap[index]` by its current monomial.
    fn sift_down(heap: &mut [usize], mut index: usize, current: &[u32], nvars: usize) {
        let key = |stream: usize| &current[stream * nvars..(stream + 1) * nvars];
        loop {
            let left = 2 * index + 1;
            if left >= heap.len() {
                return;
            }
            let right = left + 1;
            let smaller = if right < heap.len() && key(heap[right]) < key(heap[left]) {
                right
            } else {
                left
            };
            if key(heap[smaller]) < key(heap[index]) {
                heap.swap(index, smaller);
                index = smaller;
            } else {
                return;
            }
        }
    }

    /// Sum the contributions of one frontier target in a single pass.
    ///
    /// Multiplying a sorted polynomial by a monomial keeps it sorted, and every
    /// source map is a prefix of the shared map, so the target is a k-way merge of
    /// lazily shifted and scaled sources: each output monomial is written once,
    /// and no product of a source with a monomial is materialized. A source with
    /// no use left is released, and one used unchanged by a single target moves.
    fn frontier_sum(
        sources: &mut [Option<ScalarPolynomial>],
        uses: &mut [usize],
        contributions: &[Contribution],
        map: &Arc<Vec<PolyVariable>>,
    ) -> Option<ScalarPolynomial> {
        let nvars = map.len();
        for contribution in contributions {
            let source = sources.get(contribution.source)?.as_ref()?;
            if source.nvars() > nvars
                || source.get_vars_ref() != &map[..source.nvars()]
                || contribution.shift.len() > nvars
            {
                return None;
            }
        }
        for contribution in contributions {
            uses[contribution.source] -= 1;
        }
        if let [contribution] = contributions {
            let source = sources[contribution.source].as_ref()?;
            if uses[contribution.source] == 0
                && source.nvars() == nvars
                && Q.is_one(&contribution.coefficient)
                && contribution.shift.iter().all(|&exponent| exponent == 0)
            {
                return sources[contribution.source].take();
            }
        }
        let sum = {
            let streams = contributions
                .iter()
                .map(|contribution| Some((sources[contribution.source].as_ref()?, contribution)))
                .collect::<Option<Vec<_>>>()?;
            if streams
                .iter()
                .any(|(source, _)| source.get_vars_ref() != &map[..source.nvars()])
            {
                return None;
            }
            // The input term count bounds the output: reserving it avoids the
            // repeated copies of a growing buffer. Pages never written are not
            // touched, so the unused tail costs address space only.
            let capacity = streams
                .iter()
                .map(|(source, _)| source.nterms())
                .sum::<usize>();
            let mut sum = ScalarPolynomial::new(&Q, Some(capacity), Arc::clone(map));
            let mut cursors = vec![0usize; streams.len()];
            let mut current = vec![0u32; streams.len() * nvars];
            let load = |stream: usize, cursor: usize, current: &mut [u32]| -> Option<()> {
                let (source, contribution) = streams[stream];
                let width = source.nvars();
                let exponents = &source.exponents[cursor * width..(cursor + 1) * width];
                let target = &mut current[stream * nvars..(stream + 1) * nvars];
                for (variable, entry) in target.iter_mut().enumerate() {
                    *entry = exponents
                        .get(variable)
                        .copied()
                        .unwrap_or(0)
                        .checked_add(contribution.shift.get(variable).copied().unwrap_or(0))?;
                }
                Some(())
            };
            let mut heap = Vec::with_capacity(streams.len());
            for (stream, (source, _)) in streams.iter().enumerate() {
                if !source.is_zero() {
                    load(stream, 0, &mut current)?;
                    heap.push(stream);
                }
            }
            for index in (0..heap.len() / 2).rev() {
                CoefficientSum::sift_down(&mut heap, index, &current, nvars);
            }
            while let Some(&top) = heap.first() {
                let (source, contribution) = streams[top];
                let cursor = cursors[top];
                let coefficient = if Q.is_one(&contribution.coefficient) {
                    source.coefficients[cursor].clone()
                } else {
                    Q.mul(&source.coefficients[cursor], &contribution.coefficient)
                };
                sum.append_monomial_back(coefficient, &current[top * nvars..(top + 1) * nvars]);
                cursors[top] += 1;
                if cursors[top] < source.nterms() {
                    load(top, cursors[top], &mut current)?;
                } else {
                    heap.swap_remove(0);
                }
                CoefficientSum::sift_down(&mut heap, 0, &current, nvars);
            }
            sum
        };
        for contribution in contributions {
            if uses[contribution.source] == 0 {
                sources[contribution.source] = None;
            }
        }
        Some(sum)
    }

    /// Largest exponent of every variable of `weight`.
    fn max_degrees(weight: &ScalarPolynomial) -> Vec<u32> {
        let nvars = weight.nvars();
        let mut degrees = vec![0u32; nvars];
        if nvars > 0 {
            for exponents in weight.exponents.chunks_exact(nvars) {
                for (degree, &exponent) in degrees.iter_mut().zip(exponents) {
                    *degree = (*degree).max(exponent);
                }
            }
        }
        degrees
    }

    fn check_product_degrees(degrees: &[u32], factor: &ScalarPolynomial) -> Option<()> {
        for term in factor.into_iter() {
            for position in 0..degrees.len().max(term.exponents.len()) {
                let degree = degrees.get(position).copied().unwrap_or(0);
                let exponent = term.exponents.get(position).copied().unwrap_or(0);
                // Symbolica's expression writer represents powers as i32.
                // Decline before multiplication can produce an unwriteable
                // exponent, including variables newly added to the shared map.
                i32::try_from(degree.checked_add(exponent)?).ok()?;
            }
        }
        Some(())
    }

    fn push(
        &mut self,
        mut coefficient: ScalarPolynomial,
        #[cfg(feature = "reference-cases")] surrounding_bytes: usize,
    ) {
        if coefficient.is_zero() {
            return;
        }
        if let Self::One(old) = self {
            if old.is_zero() {
                *old = coefficient;
                return;
            }
            let first = old.nterms().ilog2() as usize;
            let second = coefficient.nterms().ilog2() as usize;
            let old = mem::replace(old, ScalarPolynomial::new_zero(&Q));
            if first == second {
                *self = Self::One(Self::merge(
                    old,
                    coefficient,
                    #[cfg(feature = "reference-cases")]
                    surrounding_bytes.saturating_add(mem::size_of::<Self>()),
                ));
            } else {
                let mut bins = Vec::new();
                bins.resize_with(first.max(second) + 1, || None);
                bins[first] = Some(old);
                bins[second] = Some(coefficient);
                *self = Self::Binned(bins);
            }
            return;
        }
        let Self::Binned(bins) = self else {
            unreachable!("the single-polynomial case returned above")
        };
        while !coefficient.is_zero() {
            let bin = coefficient.nterms().ilog2() as usize;
            if bins.len() <= bin {
                bins.resize_with(bin + 1, || None);
            }
            let Some(old) = bins[bin].take() else {
                bins[bin] = Some(coefficient);
                return;
            };
            coefficient = Self::merge(
                old,
                coefficient,
                #[cfg(feature = "reference-cases")]
                surrounding_bytes.saturating_add(Self::binned_bytes(bins)),
            );
            // Cancellation can move the result to a smaller bin. Each merge
            // consumes an occupied bin, so the loop also terminates then.
        }
    }

    fn finish(
        self,
        #[cfg(feature = "reference-cases")] surrounding_bytes: usize,
    ) -> ScalarPolynomial {
        let bins = match self {
            Self::One(coefficient) => return coefficient,
            Self::Binned(bins) => bins,
        };
        #[cfg(feature = "reference-cases")]
        let mut remaining_bytes = Self::binned_bytes(&bins);
        let mut sum = ScalarPolynomial::new_zero(&Q);
        for coefficient in bins.into_iter().flatten() {
            #[cfg(feature = "reference-cases")]
            {
                // The iterator retains the allocation while each moved
                // coefficient leaves the unconsumed payload charge. Keep a
                // saturated bound conservative after the exact total is lost.
                if remaining_bytes != usize::MAX {
                    remaining_bytes -= ComponentSum::coefficient_bytes(&coefficient)
                        - mem::size_of::<ScalarPolynomial>();
                }
            }
            // Smallest first: the final fold never repeatedly copies the
            // largest partial while adding the remaining small ones.
            sum = Self::merge(
                sum,
                coefficient,
                #[cfg(feature = "reference-cases")]
                surrounding_bytes.saturating_add(remaining_bytes),
            );
        }
        sum
    }

    fn merge(
        first: ScalarPolynomial,
        second: ScalarPolynomial,
        #[cfg(feature = "reference-cases")] surrounding_bytes: usize,
    ) -> ScalarPolynomial {
        #[cfg(feature = "reference-cases")]
        let _merge_timing = timing::scope(Phase::CoefficientMerge);
        #[cfg(feature = "reference-cases")]
        {
            // Addition may widen both variable maps and allocate a merged
            // output. Include all other retained partials through the caller's
            // surrounding storage, not just the two inputs being merged.
            let variables = first.nvars().saturating_add(second.nvars());
            let inputs = first
                .nterms()
                .saturating_add(second.nterms())
                .saturating_mul(
                    mem::size_of::<Rational>()
                        .saturating_add(variables.saturating_mul(mem::size_of::<u32>())),
                );
            let maps = variables.saturating_mul(mem::size_of::<PolyVariable>());
            let retained = ComponentSum::coefficient_bytes(&first)
                .saturating_add(ComponentSum::coefficient_bytes(&second));
            timing::frontier_workspace(
                surrounding_bytes
                    .saturating_add(retained)
                    .saturating_add(inputs.saturating_add(maps).saturating_mul(2)),
            );
        }
        if first.is_zero() {
            return second;
        }
        if second.is_zero() {
            return first;
        }
        // Symbolica inherits the first operand's variable order. Preserve the
        // larger map so a different ordering only reconstructs the smaller
        // polynomial; reconstructing the larger one can require quadratic
        // monomial insertion.
        if first.nterms() < second.nterms() {
            second + first
        } else {
            first + second
        }
    }

    #[cfg(feature = "reference-cases")]
    fn binned_bytes(bins: &Vec<Option<ScalarPolynomial>>) -> usize {
        mem::size_of::<Self>()
            .saturating_add(
                bins.capacity()
                    .saturating_mul(mem::size_of::<Option<ScalarPolynomial>>()),
            )
            .saturating_add(
                bins.iter()
                    .flatten()
                    .map(|coefficient| {
                        ComponentSum::coefficient_bytes(coefficient)
                            - mem::size_of::<ScalarPolynomial>()
                    })
                    .fold(0usize, usize::saturating_add),
            )
    }

    #[cfg(all(test, feature = "reference-cases"))]
    fn storage_bytes(&self) -> usize {
        match self {
            Self::One(coefficient) => {
                mem::size_of::<Self>().saturating_add(ComponentSum::coefficient_bytes(coefficient))
                    - mem::size_of::<ScalarPolynomial>()
            }
            Self::Binned(bins) => Self::binned_bytes(bins),
        }
    }
}

/// Sources with fewer terms are multiplied out and summed at once; larger
/// ones are merged lazily (see [`CoefficientSum::frontier_sum`]).
const LAZY_SOURCE_TERMS: usize = 32;

/// One contribution to a frontier target: a source state times one monomial,
/// `coefficient` times the variables of the shared map raised to `shift`.
struct Contribution {
    source: usize,
    coefficient: Rational,
    shift: Vec<u32>,
}

/// The scalar variables of one frontier, in one append-only polynomial
/// variable map: products and merges of frontier coefficients then never
/// reorder and rebuild a polynomial. When enabled, small composite scalar
/// spectators (such as `2-dim` or a short closed trace) are expanded over their
/// base variables, so equal monomials merge instead of carrying an opaque copy
/// of the sum.
struct SharedScalars<'a> {
    contractor: &'a SlotContraction,
    atoms: Vec<Atom>,
    positions: AHashMap<Atom, usize>,
    dot_positions: AHashMap<Atom, usize>,
    /// Dot-notation spelling of every position.
    dots: Vec<Atom>,
    map: Arc<Vec<PolyVariable>>,
    expansions: AHashMap<usize, Expansion>,
}

/// Composite scalars with at most this many polynomial terms are expanded
/// over their base variables; larger ones stay opaque variables.
const EXPANSION_TERM_LIMIT: usize = 64;

/// Composites whose structural term bound (products of sum lengths, powers
/// expanded) exceeds this are not converted at all, so that a large spectator
/// such as `(x+y)^1000` is never expanded only to be discarded.
const EXPANSION_WORK_LIMIT: usize = 1024;

enum Expansion {
    Base(usize),
    Composite(ScalarPolynomial),
}

struct ExpansionBound {
    terms: usize,
    degree: u32,
}

impl<'a> SharedScalars<'a> {
    // Bound conversion work and every intermediate monomial degree before
    // asking Symbolica to perform unchecked polynomial arithmetic.
    fn bound(view: AtomView<'_>) -> Option<ExpansionBound> {
        match view {
            AtomView::Add(sum) => sum.iter().try_fold(
                ExpansionBound {
                    terms: 0,
                    degree: 0,
                },
                |total, term| {
                    let next = Self::bound(term)?;
                    Some(ExpansionBound {
                        terms: total.terms.checked_add(next.terms)?,
                        degree: total.degree.max(next.degree),
                    })
                },
            ),
            AtomView::Mul(product) => product.iter().try_fold(
                ExpansionBound {
                    terms: 1,
                    degree: 0,
                },
                |total, factor| {
                    let next = Self::bound(factor)?;
                    Some(ExpansionBound {
                        terms: total.terms.checked_mul(next.terms)?,
                        degree: total.degree.checked_add(next.degree)?,
                    })
                },
            ),
            AtomView::Pow(power) => {
                let (base, exponent) = power.get_base_exp();
                let exponent = u32::try_from(i64::try_from(exponent).ok()?).ok()?;
                let base = Self::bound(base)?;
                Some(ExpansionBound {
                    terms: base.terms.checked_pow(exponent)?,
                    degree: base.degree.checked_mul(exponent)?,
                })
            }
            AtomView::Num(_) => Some(ExpansionBound {
                terms: 1,
                degree: 0,
            }),
            _ => Some(ExpansionBound {
                terms: 1,
                degree: 1,
            }),
        }
    }

    fn composite(view: AtomView<'_>) -> bool {
        match view {
            AtomView::Add(_) => true,
            AtomView::Mul(product) => product.iter().any(Self::composite),
            AtomView::Pow(power) => {
                let (base, exponent) = power.get_base_exp();
                matches!(exponent, AtomView::Num(_)) && Self::composite(base)
            }
            _ => false,
        }
    }

    /// Whether `value` has no explicit index slot: a composite with internal dummy
    /// indices must not be distributed, since that repeats its index labels.
    fn index_free(slots: &mut SlotMatcher, value: AtomView<'_>) -> bool {
        let mut free = true;
        value.visitor(&mut |node| {
            if matches!(slots.classify(node), SlotMatch::Explicit(_)) {
                free = false;
            }
            free
        });
        free
    }

    /// Whether `terms^exponent` stays within the expansion limit.
    fn small_power(terms: usize, exponent: u32) -> bool {
        terms
            .checked_pow(exponent)
            .is_some_and(|terms| terms <= EXPANSION_TERM_LIMIT)
    }

    fn new(contractor: &'a SlotContraction) -> Self {
        Self {
            contractor,
            atoms: Vec::new(),
            positions: AHashMap::new(),
            dot_positions: AHashMap::new(),
            dots: Vec::new(),
            map: Arc::new(Vec::new()),
            expansions: AHashMap::new(),
        }
    }

    /// One position per scalar up to dot spelling: `g(p,q)` from this frontier
    /// and `dot(p,q)` from a nested contraction are one variable, so their
    /// exponents add instead of forming distinct monomials.
    fn position(&mut self, atom: Atom) -> usize {
        if let Some(&position) = self.positions.get(&atom) {
            return position;
        }
        let dots = DotNormalizer::notation(atom.as_view(), self.contractor.normalization_matcher());
        if let Some(&position) = self.dot_positions.get(&dots) {
            self.positions.insert(atom, position);
            return position;
        }
        let position = self.atoms.len();
        self.atoms.push(atom.clone());
        self.dots.push(dots.clone());
        self.dot_positions.insert(dots, position);
        self.positions.insert(atom, position);
        Arc::make_mut(&mut self.map).push(PolyVariable::Temporary(position));
        position
    }

    /// The expansion of the scalar variable `index`, whose atom is `atom`.
    /// `index_free` is asked once, and only when expansion is enabled.
    fn expansion(
        &mut self,
        index: usize,
        atom: &Atom,
        index_free: impl FnOnce(&Atom) -> bool,
    ) -> &Expansion {
        if !self.expansions.contains_key(&index) {
            let expansion = self.expand(atom, index_free);
            self.expansions.insert(index, expansion);
        }
        &self.expansions[&index]
    }

    fn expand(&mut self, atom: &Atom, index_free: impl FnOnce(&Atom) -> bool) -> Expansion {
        if self.contractor.expand_scalars
            && Self::composite(atom.as_view())
            && Self::bound(atom.as_view()).is_some_and(|bound| {
                bound.terms <= EXPANSION_WORK_LIMIT && i32::try_from(bound.degree).is_ok()
            })
            && index_free(atom)
            && let Ok(polynomial) = atom.try_to_polynomial::<_, u32>(&Q, None)
            && polynomial.nterms() <= EXPANSION_TERM_LIMIT
        {
            let positions = polynomial
                .variables()
                .iter()
                .map(|variable| self.position(variable.to_atom()))
                .collect::<Vec<_>>();
            let mut expanded =
                ScalarPolynomial::new(&Q, Some(polynomial.nterms()), Arc::clone(&self.map));
            let mut exponents = vec![0u32; self.map.len()];
            let mut overflow = false;
            for term in polynomial.into_iter() {
                exponents.fill(0);
                // Two variables of the composite (`g(p,q)` and `dot(p,q)`)
                // can share one position; their exponents add.
                for (&position, &exponent) in positions.iter().zip(term.exponents) {
                    match exponents[position].checked_add(exponent) {
                        Some(sum) => exponents[position] = sum,
                        None => overflow = true,
                    }
                }
                expanded.append_monomial(term.coefficient.clone(), &exponents);
            }
            if !overflow {
                return Expansion::Composite(expanded);
            }
        }
        Expansion::Base(self.position(atom.clone()))
    }

    /// Emit each variable in the spelling it first had in this frontier, or
    /// in dot notation when the request ends there. The polynomial is
    /// relabelled with those atoms and written into one expression buffer,
    /// which is normalized once.
    fn emit(&self, mut polynomial: ScalarPolynomial) -> Option<Atom> {
        let spellings = if self.contractor.dot_notation {
            &self.dots
        } else {
            &self.atoms
        };
        let variables = polynomial
            .get_vars_ref()
            .iter()
            .map(|variable| {
                let PolyVariable::Temporary(position) = variable else {
                    return None;
                };
                let atom = spellings.get(*position)?.clone();
                // Only the expression of a variable is read when writing.
                Some(match atom.as_view() {
                    AtomView::Var(variable) => PolyVariable::Symbol(variable.get_symbol()),
                    AtomView::Fun(function) => PolyVariable::Function(function.get_symbol(), atom),
                    _ => PolyVariable::Power(atom),
                })
            })
            .collect::<Option<Vec<_>>>()?;
        let mut labelled = ScalarPolynomial::new(&Q, None, Arc::new(variables));
        labelled.coefficients = mem::take(&mut polynomial.coefficients);
        labelled.exponents = mem::take(&mut polynomial.exponents);
        Some(
            labelled.to_expression_with_coeff_map(|_, coefficient, out| {
                out.to_num(coefficient.clone());
            }),
        )
    }
}

type FactorTerm = (Vec<(usize, u32)>, Rational);
// Selected expression/interface and retained factors.
type PrerequisiteFactors = (Atom, PartialStructure, Vec<Atom>);

pub(crate) use crate::tensor::simplification::ReductionStatus;

/// Exact output may retain unsupported contractions for subsequent operations.
/// Algebraic validity and completion are separate facts.
#[derive(Clone)]
pub(crate) struct FactorizedContraction<Root = Atom> {
    pub(crate) root: Root,
    pub(crate) status: ReductionStatus,
    pub(crate) observations: Option<DomainObservations>,
    /// A published Atom result owns only facts reachable from its current root.
    /// Typed roots carry their facts directly and leave this scratch cell empty.
    pub(crate) proofs: TensorProofs,
}

impl FactorizedContraction {
    /// Attach the producer's retained facts before accounting for generated indices.
    pub(crate) fn into_tensor(
        self,
        source: &SymbolicTensor<PartialStructure>,
    ) -> Result<FactorizedContraction<SymbolicTensor<PartialStructure>>, TensorInferenceError> {
        let root = source.with_published_materialization(
            "factorized slot contraction",
            self.root,
            self.proofs,
        )?;
        let root = if let Some(observations) = self.observations {
            root.with_identity_result(root.expression.clone(), Some(observations))?
        } else {
            root
        };
        Ok(FactorizedContraction {
            root,
            status: self.status,
            observations: None,
            proofs: TensorProofs::default(),
        })
    }

    /// An open terminal certificate describes its root, never scalar syntax.
    fn scalar_facts(observed: &DomainObservations) -> bool {
        observed.is_closed_scalar()
    }
}

// A planned domain already owns its topology, structures and leaf proofs.
// Scratch parses keep their borrowed source views for selectively opening leaves.
struct ScopeNetwork<'src> {
    graph: Cow<'src, NetworkGraph<DummyKey, Symbol>>,
    store: ScopeStore<'src>,
    admission: Option<Arc<TensorProofs>>,
}

type ScopeLeafStore<E> = NetworkStore<PlanningLeaf<E>, E>;

enum ScopeStore<'src> {
    Planned(&'src ScopeLeafStore<Atom>),
    Parsed(ScopeLeafStore<AtomOrView<'src>>),
}

impl<'src> ScopeNetwork<'src> {
    fn planned(network: &'src PlanningNetwork<Atom>) -> Self {
        Self {
            graph: Cow::Borrowed(&network.graph),
            store: ScopeStore::Planned(&network.store),
            admission: network.admission.clone(),
        }
    }

    fn parsed(network: PlanningNetwork<AtomOrView<'src>>) -> Self {
        Self {
            graph: Cow::Owned(network.network.graph),
            store: ScopeStore::Parsed(network.network.store),
            admission: network.admission,
        }
    }

    fn structure(&self, index: usize) -> &OrderedStructure {
        match &self.store {
            ScopeStore::Planned(store) => &store.tensors[index].structure,
            ScopeStore::Parsed(store) => &store.tensors[index].structure,
        }
    }

    fn leaf<'a, E>(
        &self,
        store: &'a ScopeLeafStore<E>,
        node: NodeIndex,
    ) -> Option<(&'a E, Option<&'a PlanningLeaf<E>>)> {
        match self.graph.graph[node] {
            NetworkNode::Leaf(NetworkLeaf::LocalTensor(index)) => {
                let tensor = &store.tensors[index];
                Some((&tensor.expression, Some(tensor)))
            }
            NetworkNode::Leaf(NetworkLeaf::Scalar(index)) => {
                Some((store.get_scalar_ref(index), None))
            }
            _ => None,
        }
    }

    fn literal(&self, node: NodeIndex) -> Option<AtomView<'_>> {
        match &self.store {
            ScopeStore::Planned(store) => Some(self.leaf(store, node)?.0.as_view()),
            ScopeStore::Parsed(store) => Some(self.leaf(store, node)?.0.as_view()),
        }
    }

    fn openable_leaf(
        &self,
        node: NodeIndex,
    ) -> Option<(AtomView<'src>, Option<PlanningSource<'src>>)> {
        match &self.store {
            ScopeStore::Planned(store) => {
                let (value, admission) = self.leaf(*store, node)?;
                Some((
                    value.as_view(),
                    admission.map(|leaf| leaf.source(value.as_view(), self.admission.as_ref())),
                ))
            }
            ScopeStore::Parsed(store) => {
                // Owned coefficients of a scratch parse have no tensor interior.
                let (AtomOrView::View(value), admission) = self.leaf(store, node)? else {
                    return None;
                };
                Some((
                    *value,
                    admission.map(|leaf| leaf.source(*value, self.admission.as_ref())),
                ))
            }
        }
    }
}

// A parse belongs to one arithmetic scope. Opening an occurrence is cached on
// that scope, so borrowed payloads and its port atoms outlive every DP state.
// This is parser scratch, not another symbolic tensor representation.
struct ContractionScope<'src> {
    source: AtomView<'src>,
    network: ScopeNetwork<'src>,
    tree: SimpleTraversalTree<ChildVecStore<()>>,
    root: NodeIndex,
    ports: AHashMap<NodeIndex, (Vec<Atom>, PartialStructure)>,
    opened: AHashMap<NodeIndex, OnceCell<Option<Box<Self>>>>,
    reductions: AHashMap<NodeIndex, OnceCell<Option<FactorizedContraction>>>,
    // Each scope serves one fixed contraction or admission policy.
    work: AHashMap<NodeIndex, OnceCell<bool>>,
    interfaces: Option<LeafInterfaceCache>,
}

impl<'src> ContractionScope<'src> {
    #[cfg(test)]
    fn parse(source: AtomView<'src>) -> Option<Self> {
        let network = SlotContraction::shallow_network(source, None, None).ok()?;
        Self::from_network(source, ScopeNetwork::parsed(network), None)
    }

    fn from_network(
        source: AtomView<'src>,
        network: ScopeNetwork<'src>,
        interfaces: Option<LeafInterfaceCache>,
    ) -> Option<Self> {
        let graph = &network.graph;
        let tree = graph.expr_tree().cast::<ChildVecStore<()>>();
        let root = graph.graph.node_id(graph.head());
        let nodes = tree
            .iter_preorder_tree_nodes(&graph.graph, root)
            .collect::<Vec<_>>();
        let external: SuBitGraph = graph.graph.external_filter();
        // Preorder makes each arithmetic subtree a contiguous interval. Reuse
        // those intervals for every boundary instead of allocating a descendant
        // set and traversing the same tree for each node.
        let positions = nodes
            .iter()
            .copied()
            .enumerate()
            .map(|(position, node)| (node, position))
            .collect::<AHashMap<_, _>>();
        let mut ends = (1..=nodes.len()).collect::<Vec<_>>();
        for (position, &node) in nodes.iter().enumerate().rev() {
            for child in tree.iter_children(node, &graph.graph) {
                ends[position] = ends[position].max(ends[positions[&child]]);
            }
        }
        // The parser already observed each leaf's written ports. Retain their
        // literal index payloads when describing connected boundaries instead
        // of rebuilding names through user-defined normalizers.
        let mut written_slots = AHashMap::new();
        if let Some(interfaces) = &interfaces {
            let interfaces = interfaces.lock().unwrap();
            let mut matcher = SlotMatcher::default();
            for &node in &nodes {
                if matches!(
                    &graph.graph[node],
                    NetworkNode::Leaf(NetworkLeaf::LocalTensor(_))
                ) && let Some((_, ports)) = interfaces.get(network.literal(node)?.get_data())
                {
                    for atom in ports {
                        if let Ok(slot) = matcher.parse::<LibraryRep, AbstractIndex>(atom.as_view())
                        {
                            written_slots.entry(slot).or_insert_with(|| atom.clone());
                        }
                    }
                }
            }
        }
        let mut ports = AHashMap::with_capacity(nodes.len());
        for (begin, &node) in nodes.iter().enumerate() {
            let end = ends[begin];
            let mut boundary = Vec::new();
            for &member in &nodes[begin..end] {
                for hedge in graph.graph.iter_crown(member) {
                    let NetworkEdge::Slot(slot) = graph.graph[[&hedge]] else {
                        continue;
                    };
                    if !external.includes(&hedge)
                        && positions
                            .get(&graph.graph.node_id(graph.graph.inv(hedge)))
                            .is_some_and(|&other| begin <= other && other < end)
                    {
                        continue;
                    }
                    // Sewing shares one edge label; the occurrence's storage
                    // axis still owns this endpoint's variance.
                    let slot = match &graph.graph[member] {
                        NetworkNode::Leaf(NetworkLeaf::LocalTensor(index)) => network
                            .structure(*index)
                            .get_slot(usize::from(graph.slot_order[hedge.0]))?,
                        _ => slot,
                    };
                    boundary.push(slot);
                }
            }
            boundary.sort();
            ports.insert(
                node,
                (
                    boundary
                        .iter()
                        .map(|slot| {
                            written_slots
                                .get(slot)
                                .cloned()
                                .unwrap_or_else(|| slot.to_atom())
                        })
                        .collect(),
                    PartialStructure::from_logical_slots(
                        boundary
                            .into_iter()
                            .map(|slot| slot.rep().slot(PartialIndex::Explicit(slot.aind()))),
                    ),
                ),
            );
        }
        Some(Self {
            source,
            network,
            tree,
            root,
            ports,
            opened: nodes
                .iter()
                .map(|&node| (node, Default::default()))
                .collect(),
            reductions: nodes
                .iter()
                .map(|&node| (node, Default::default()))
                .collect(),
            work: nodes
                .iter()
                .map(|&node| (node, Default::default()))
                .collect(),
            interfaces,
        })
    }

    fn literal(&self, node: NodeIndex) -> Option<AtomView<'_>> {
        self.network
            .literal(node)
            .or_else(|| (node == self.root).then_some(self.source))
    }

    fn open(&self, node: NodeIndex) -> Option<&Self> {
        self.opened[&node]
            .get_or_init(|| {
                let (value, admission) = self.network.openable_leaf(node)?;
                #[cfg(test)]
                crate::tensor::SELECTED_OPENS.with(|count| count.set(count.get() + 1));
                let network =
                    SlotContraction::shallow_network(value, self.interfaces.as_ref(), admission)
                        .ok()?;
                let opened = Self::from_network(
                    value,
                    ScopeNetwork::parsed(network),
                    self.interfaces.clone(),
                )?;
                if matches!(
                    opened.network.graph.graph[opened.root],
                    NetworkNode::Leaf(_)
                ) && opened.literal(opened.root)? == value
                {
                    // Unsupported powers and opaque scalar syntax do not expose
                    // arithmetic. Decline instead of reopening the same leaf.
                    return None;
                }
                Some(Box::new(opened))
            })
            .as_deref()
    }

    fn contains_work(
        &self,
        node: NodeIndex,
        contractor: &SlotContraction,
        settings: ContractSettings<'_>,
    ) -> bool {
        *self.work[&node].get_or_init(|| {
            let Some(value) = self.literal(node) else {
                return false;
            };
            if let Some(observed) = &contractor.observations {
                observed.region(value).is_none_or(|region| {
                    // Internal colour/gamma connections are not metric/vector
                    // work. Keep these regions opaque, including their open
                    // interfaces: incident sources can bind those ports without
                    // opening the payload. Minimal notation admission still
                    // needs its independent chain/trace candidates. Expanding
                    // collection also merges alpha-equivalent tensor terms.
                    (settings.expand
                        || settings.collect_chains
                        || settings.collect_traces
                        || !region.excludes_contraction_sources(settings))
                        && (region.has_indices() || region.counts[7] > 0)
                })
            } else {
                let observed = SimplificationCandidates::scan(value, [contractor.metric], || true);
                observed.symbols[0] || observed.dots || observed.repeated_indices
            }
        })
    }

    fn factors(&self, root: NodeIndex) -> Vec<NodeIndex> {
        if matches!(
            self.network.graph.graph[root],
            NetworkNode::Op(NetworkOp::Product)
        ) {
            let mut nodes = self
                .tree
                .iter_children(root, &self.network.graph.graph)
                .collect::<Vec<_>>();
            nodes.reverse();
            nodes
        } else {
            vec![root]
        }
    }
}

impl SlotContraction {
    fn contraction_settings(&self, rank_one: bool) -> ContractSettings<'_> {
        ContractSettings {
            metrics: self.metrics,
            representations: self.representations.as_deref(),
            rank_one,
            collect_chains: false,
            collect_traces: false,
            expand: self.expand,
            ..Default::default()
        }
    }

    /// Budget admission uses the same shallow scopes and boundary edges as the
    /// contraction. It does not normalize a tensor or invoke user callbacks.
    pub(crate) fn minimal_work_permitted(
        &self,
        source: AtomView<'_>,
        settings: ContractSettings<'_>,
    ) -> Option<bool> {
        let scope = self.scope(source)?;
        self.minimal_scope_work(&scope, scope.root, settings)
    }

    fn minimal_scope_work(
        &self,
        scope: &ContractionScope<'_>,
        root: NodeIndex,
        settings: ContractSettings<'_>,
    ) -> Option<bool> {
        match scope.network.graph.graph[root] {
            NetworkNode::Op(NetworkOp::Sum) => {
                for child in scope.tree.iter_children(root, &scope.network.graph.graph) {
                    if self.minimal_scope_work(scope, child, settings)? {
                        return Some(true);
                    }
                }
                return Some(false);
            }
            NetworkNode::Op(NetworkOp::Power(_)) => {
                let child = scope
                    .tree
                    .iter_children(root, &scope.network.graph.graph)
                    .next()?;
                if !scope.ports[&child].0.is_empty()
                    && !self.tensor_alternatives(scope.literal(child)?)
                {
                    // An atomic tensor power can connect copies directly.
                    let sources = self
                        .observations
                        .as_ref()
                        .and_then(|observed| observed.region(scope.literal(child)?))
                        .map(|region| region.contraction_sources(settings))?;
                    return Some(
                        (sources || (settings.collect_traces && scope.ports[&child].0.len() == 2))
                            && scope.ports[&child]
                                .1
                                .logical_slots()
                                .iter()
                                .any(|slot| self.permits(slot.rep().rep)),
                    );
                }
                return self.minimal_scope_work(scope, child, settings);
            }
            _ => {}
        }
        let nodes = scope.factors(root);
        let external: SuBitGraph = scope.network.graph.graph.external_filter();
        let components = scope.network.graph.slot_components(&scope.tree, &nodes);
        for &node in &nodes {
            let literal = scope.literal(node)?;
            let component = components
                .iter()
                .find(|component| component.contains(&node))?;
            let alternatives = component
                .iter()
                .copied()
                .filter(|&member| {
                    !scope.ports[&member].0.is_empty()
                        && scope
                            .literal(member)
                            .is_some_and(|source| self.tensor_alternatives(source))
                })
                .collect::<Vec<_>>();
            let open_component = component.iter().any(|&member| {
                scope.network.graph.graph.iter_crown(member).any(|edge| {
                    matches!(scope.network.graph.graph[[&edge]], NetworkEdge::Slot(_))
                        && external.includes(&edge)
                })
            });
            let matrix = scope.ports[&node].0.len() == 2
                && (settings.collect_traces || (settings.collect_chains && open_component));
            let connected = scope.network.graph.graph.iter_crown(node).any(|edge| {
                matches!(scope.network.graph.graph[[&edge]], NetworkEdge::Slot(slot)
                    if self.permits(slot.rep().rep))
                    && !external.includes(&edge)
                    && component.contains(
                        &scope
                            .network
                            .graph
                            .graph
                            .node_id(scope.network.graph.graph.inv(edge)),
                    )
            });
            if alternatives.len() == 1
                && alternatives[0] == node
                && connected
                && (matrix
                    || self
                        .observations
                        .as_ref()
                        .and_then(|observed| observed.region(literal))?
                        .contraction_sources(settings))
            {
                // One selected existing sum may supply a metric/vector to an
                // atomic partner. It never combines independent alternatives.
                return Some(true);
            }
            if let AtomView::Fun(function) = literal {
                let head = function.get_symbol();
                let source = (settings.metrics && head == self.metric)
                    || (settings.rank_one && head.has_tag(&self.tags.rank1));
                if (source || matrix) && connected {
                    return Some(true);
                }
                if !source
                    && scope.ports[&node].0.is_empty()
                    && settings.collect_traces
                    && self
                        .observations
                        .as_ref()
                        .and_then(|observed| observed.region(literal))?
                        .has_indices()
                {
                    // An internally closed matrix leaf has no graph boundary;
                    // its notation owner establishes trace eligibility.
                    return None;
                }
                if head == self.tags.chain || head == self.tags.trace || head == self.tags.bracket {
                    // These notation owners carry ports inside their words.
                    // Keep their existing admission and checked rewrite rules.
                    return None;
                }
            } else if scope.contains_work(node, self, settings) {
                if !matches!(scope.network.graph.graph[node], NetworkNode::Leaf(_)) {
                    if self.minimal_scope_work(scope, node, settings)? {
                        return Some(true);
                    }
                } else {
                    let opened = scope.open(node)?;
                    if self.minimal_scope_work(opened, opened.root, settings)? {
                        return Some(true);
                    }
                }
            }
        }
        Some(false)
    }

    fn tensor_alternatives(&self, source: AtomView<'_>) -> bool {
        if self
            .observations
            .as_ref()
            .and_then(|observed| observed.region(source))
            .is_some_and(|region| !region.has_indices())
        {
            return false;
        }
        match source {
            AtomView::Add(_) => true,
            AtomView::Mul(product) => product
                .iter()
                .any(|factor| self.tensor_alternatives(factor)),
            AtomView::Pow(power) => self.tensor_alternatives(power.get_base_exp().0),
            _ => false,
        }
    }

    fn scope<'a>(&'a self, source: AtomView<'a>) -> Option<ContractionScope<'a>> {
        let network = if let Some(planned) = &self.graph
            && planned.matches(source)
        {
            ScopeNetwork::planned(&planned.graph)
        } else {
            ScopeNetwork::parsed(
                Self::shallow_network(
                    source,
                    self.interfaces.as_ref(),
                    self.proofs.as_ref().map(PlanningSource::root),
                )
                .ok()?,
            )
        };
        ContractionScope::from_network(source, network, self.interfaces.clone())
    }

    /// The shared depth-one structural view keeps arithmetic leaves borrowed.
    /// Family collection opens selected leaves through its existing tape reader.
    pub(crate) fn shallow_network<'src>(
        source: AtomView<'src>,
        interfaces: Option<&LeafInterfaceCache>,
        admission: Option<PlanningSource<'src>>,
    ) -> Result<PlanningNetwork<AtomOrView<'src>>, TensorInferenceError> {
        #[cfg(feature = "reference-cases")]
        let _phase = timing::scope(Phase::ShallowGraph);
        let settings = ParseSettings {
            precontract_scalars: false,
            depth_limit: Some(1),
            shorthand_parsing: ShorthandParsing::Opaque,
            parse_composite_scalars_as_tensors: true,
            ..ParseSettings::default()
        };
        planning::parse(source, interfaces, admission, &settings)
    }

    /// Collect notation on connected scopes of the same shallow graph. Scalar
    /// spectators and components without a selected representation are reused.
    pub(crate) fn collect_matrix_notation(
        &self,
        source: AtomView<'_>,
        representations: &[LibraryRep],
        collect: bool,
        traces: bool,
    ) -> Result<Option<Atom>, TensorInferenceError> {
        let Some(scope) = self.scope(source) else {
            return Ok(None);
        };
        let nodes = scope.factors(scope.root);
        let components = scope.network.graph.slot_components(&scope.tree, &nodes);
        let external: SuBitGraph = scope.network.graph.graph.external_filter();
        let mut output = Vec::new();
        let mut changed = false;
        for component in components {
            let Some(factors) = component
                .iter()
                .map(|&node| scope.literal(node))
                .collect::<Option<Vec<_>>>()
            else {
                return Ok(None);
            };
            let mut expression = if factors.len() == 1 {
                factors[0].to_owned()
            } else {
                Atom::mul_many(factors.iter())
            };
            let mut component_value = None;
            for &representation in representations {
                let connection = component.len() > 1
                    // An opaque arithmetic leaf can own connections that the
                    // outer graph intentionally does not expose. Regional
                    // representation facts below select its local notation.
                    || factors.iter().any(|factor| {
                        (matches!(factor, AtomView::Add(_) | AtomView::Mul(_) | AtomView::Pow(_))
                            && factor.has_repeated_explicit_indices())
                            || matches!(factor, AtomView::Fun(function)
                                if function.get_symbol() == self.tags.chain)
                    })
                    || (traces
                        && component.iter().any(|&node| {
                            scope
                                .tree
                                .iter_preorder_tree_nodes(&scope.network.graph.graph, node)
                                .any(|member| {
                                    scope.network.graph.graph.iter_crown(member).any(|edge| {
                                        let NetworkEdge::Slot(slot) =
                                            scope.network.graph.graph[[&edge]]
                                        else {
                                            return false;
                                        };
                                        slot.rep().rep.base() == representation.base()
                                            && !external.includes(&edge)
                                            && scope
                                                .network
                                                .graph
                                                .graph
                                                .node_id(scope.network.graph.graph.inv(edge))
                                                == member
                                    })
                                })
                        }));
                let relevant = connection
                    && factors.iter().any(|&factor| {
                        self.observations
                            .as_ref()
                            .and_then(|observed| observed.region(factor))
                            .is_none_or(|region| region.has_chain_channel(representation))
                    });
                if relevant {
                    let value = match component_value.take() {
                        Some(value) => value,
                        None => {
                            let interfaces = component
                                .iter()
                                .map(|node| scope.ports[node].1.clone())
                                .collect::<Vec<_>>();
                            let interface =
                                InterfaceInference::merge_explicit_interface_sequence(&interfaces)?;
                            self.proofs.as_ref().and_then(|proofs| proofs.region_with_interface(expression.as_view(), interface)).ok_or_else(|| TensorInferenceError::invalid("matrix collection requires its admitted component boundary"))?
                        }
                    };
                    let next = value.collect_chain_words(representation, collect, traces)?;
                    changed |= next.expression != expression;
                    expression = next.expression.clone();
                    component_value = Some(next);
                }
            }
            if let Some(value) = component_value {
                self.generated.borrow_mut().push(value.proofs);
            }
            output.push(expression);
        }
        Ok(Some(if changed {
            Atom::mul_many(output)
        } else {
            source.to_owned()
        }))
    }

    /// Select only a family's factors and metric/vector paths incident to them.
    /// Connectivity authorizes these prerequisites, never a foreign tensor body.
    pub(crate) fn prerequisite_factors(
        &self,
        source: AtomView<'_>,
        mut selected: impl FnMut(AtomView<'_>) -> bool,
        allowed_vector: &dyn Fn(LibraryRep) -> bool,
    ) -> Option<PrerequisiteFactors> {
        let scope = self.scope(source)?;
        let nodes = scope.factors(scope.root);
        let graph = &scope.network.graph.graph;
        let mut owners = AHashMap::new();
        for (position, &node) in nodes.iter().enumerate() {
            for member in scope.tree.iter_preorder_tree_nodes(graph, node) {
                owners.insert(member, position);
            }
        }
        let mut active = nodes
            .iter()
            .map(|&node| scope.literal(node).is_some_and(&mut selected))
            .collect::<Vec<_>>();
        if !active.iter().any(|active| *active) {
            return None;
        }
        let eligible = nodes
            .iter()
            .map(|&node| {
                matches!(scope.literal(node), Some(AtomView::Fun(function))
                if (self.metrics && function.get_symbol() == self.metric)
                    || (function.get_symbol().has_tag(&self.tags.rank1)
                        && scope.ports[&node].1.logical_slots().iter()
                            .all(|slot| allowed_vector(slot.rep().rep))))
            })
            .collect::<Vec<_>>();
        let mut joined = false;
        let mut touched = vec![false; nodes.len()];
        loop {
            let mut changed = false;
            for (position, &node) in nodes.iter().enumerate() {
                if active[position] || !eligible[position] {
                    continue;
                }
                let mut attached = false;
                for member in scope.tree.iter_preorder_tree_nodes(graph, node) {
                    for edge in graph.iter_crown(member) {
                        let NetworkEdge::Slot(slot) = graph[[&edge]] else {
                            continue;
                        };
                        if self.permits(slot.rep().rep)
                            && let Some(&other) = owners.get(&graph.node_id(graph.inv(edge)))
                            && other != position
                            && active[other]
                        {
                            touched[other] = true;
                            attached = true;
                        }
                    }
                }
                if attached {
                    active[position] = true;
                    touched[position] = true;
                    changed = true;
                    joined = true;
                }
            }
            if !changed {
                break;
            }
        }
        if !joined {
            return None;
        }
        let mut factors = Vec::new();
        let mut interfaces = Vec::new();
        let mut spectators = Vec::new();
        for (position, &node) in nodes.iter().enumerate() {
            let value = scope.literal(node)?.to_owned();
            if active[position] && touched[position] {
                // A selected but disconnected scalar child is still opaque:
                // another family's occurrence is no reason to open it here.
                factors.push(value);
                interfaces.push(scope.ports[&node].1.clone());
            } else {
                spectators.push(value);
            }
        }
        let interface = InterfaceInference::merge_explicit_interface_sequence(&interfaces).ok()?;
        Some((Atom::mul_many(factors), interface, spectators))
    }

    /// Plan an admitted interface. Raw Atom adapters enter admission here once;
    /// typed callers supply their retained metadata with `admitted`.
    pub(crate) fn contract_factorized(
        &self,
        source: AtomView<'_>,
        order: Option<&[usize]>,
        rank_one: bool,
    ) -> Option<FactorizedContraction> {
        self.generated.borrow_mut().clear();
        if self.proofs.is_none() {
            let source = SymbolicTensor::observe_raw_expression(source.to_owned()).ok()?;
            return self.clone().admitted(&source).contract_factorized(
                source.expression().as_view(),
                order,
                rank_one,
            );
        }
        if rank_one && self.observations.is_none() && source.contains_symbol(self.tags.chain) {
            self.dummy_state().reserve_indices(source);
        }
        let intrinsic = self
            .graph
            .as_ref()
            .filter(|planned| planned.matches(source))
            .and(self.observations.as_ref())
            .map_or_else(
                || InterfaceInference::normalization_is_intrinsic(source),
                |observed| {
                    observed.candidates.intrinsic
                        || (!observed.candidates.complete
                            && InterfaceInference::normalization_is_intrinsic(source))
                },
            );
        if source.needs_normalization() || self.metric.get_evaluation_info().is_some() || !intrinsic
        {
            return None;
        }
        let scope = self.scope(source)?;
        let result = self.contract_graph(&scope, scope.root, order, rank_one);
        let Some(mut result) = result else {
            self.generated.borrow_mut().clear();
            return None;
        };
        result.proofs =
            self.take_generated(result.root.as_view(), scope.ports[&scope.root].1.clone());
        // This producer promises a globally terminal polynomial. Validate that
        // promise here: generic publication also carries temporary certificates
        // scoped to an identity family's policy (notably colour invariants).
        #[cfg(debug_assertions)]
        if let Some(observed) = &result.observations {
            observed.assert_terminal_result(
                &result.root,
                &scope.ports[&scope.root].1,
                SlotMatcher::with_leaf_interfaces([result.proofs.parser_leaf_interfaces()]),
            );
        }
        Some(result)
    }

    fn power_scope(
        &self,
        scope: &ContractionScope<'_>,
        node: NodeIndex,
    ) -> Option<(NodeIndex, i8, bool)> {
        let NetworkNode::Op(NetworkOp::Power(power)) = scope.network.graph.graph[node] else {
            return None;
        };
        let mut children = scope.tree.iter_children(node, &scope.network.graph.graph);
        let child = children.next()?;
        if children.next().is_some() {
            return None;
        }
        let source = scope.literal(child)?;
        if !scope.ports[&child].0.is_empty() && (self.expand || !self.tensor_alternatives(source)) {
            return None;
        }
        let internal = self
            .observations
            .as_ref()
            .and_then(|observed| observed.region(source))
            .map_or_else(
                || source.has_repeated_explicit_indices(),
                |region| region.has_indices(),
            );
        Some((child, power, internal))
    }

    fn scoped_power<'a>(
        &self,
        scope: &'a ContractionScope<'_>,
        node: NodeIndex,
        rank_one: bool,
    ) -> Option<&'a FactorizedContraction> {
        let (child, power, internal) = self.power_scope(scope, node)?;
        scope.reductions[&node]
            .get_or_init(|| {
                if !internal {
                    // A transparent bracket can remain after local dot pairing
                    // has removed its last explicit index. Lower just that
                    // selected wrapper, retaining ordinary scalar powers verbatim.
                    let value = scope.literal(child)?;
                    if matches!(value, AtomView::Fun(function)
                        if function.get_symbol() == self.tags.bracket)
                        && scope.contains_work(child, self, self.contraction_settings(rank_one))
                        && self
                            .observations
                            .as_ref()
                            .and_then(|observed| observed.region(value))
                            .and_then(|region| {
                                region.any_leaf(|leaf| {
                                    // The bracket itself marks this region as structural.
                                    // Only its tensor operands can supply contraction work.
                                    leaf.get_symbol() != Some(self.tags.bracket)
                                        && InterfaceInference::has_structured_syntax(leaf)
                                })
                            })
                            .unwrap_or(true)
                    {
                        let opened = scope.open(child)?;
                        let mut result =
                            self.contract_graph(opened, opened.root, None, rank_one)?;
                        result.root = result.root.pow(power);
                        result.observations = result
                            .observations
                            .filter(FactorizedContraction::scalar_facts);
                        return Some(result);
                    }
                    return Some(FactorizedContraction {
                        proofs: Default::default(),
                        root: scope.literal(node)?.to_owned(),
                        status: ReductionStatus::Complete,
                        observations: self.observations.as_ref().and_then(|observed| {
                            observed.certify_scalar_region(scope.literal(node)?)
                        }),
                    });
                }
                // Contract the base in its own scope before applying the
                // exponent. Closed bases keep independent dummy copies; minimal
                // open sums retain their alternatives and external connections.
                let mut result = self.contract_graph(scope, child, None, rank_one)?;
                result.root = result.root.pow(power);
                // A power of an open base connects its ports.
                result.observations = result
                    .observations
                    .filter(FactorizedContraction::scalar_facts);
                Some(result)
            })
            .as_ref()
    }

    fn factor_order(
        &self,
        scope: &ContractionScope<'_>,
        nodes: &[NodeIndex],
        mut selected: Vec<usize>,
        order: Option<&[usize]>,
    ) -> Option<Vec<usize>> {
        if let Some(order) = order {
            let factors = match scope.source {
                AtomView::Mul(product) => product.iter().collect::<Vec<_>>(),
                value => vec![value],
            };
            let mut positions = order.to_vec();
            positions.sort_unstable();
            if positions != (0..factors.len()).collect::<Vec<_>>() {
                return None;
            }
            let ranks = selected
                .iter()
                .map(|&position| {
                    let literal = scope.literal(nodes[position])?;
                    let ordinal = factors.iter().position(|&factor| factor == literal)?;
                    Some((position, order.iter().position(|&entry| entry == ordinal)?))
                })
                .collect::<Option<AHashMap<_, _>>>()?;
            selected.sort_by_key(|position| ranks[position]);
            return Some(selected);
        }
        // Minimize the live boundary of the accumulated factor set. Incidence
        // comes entirely from the depth-one graph; scoring never opens a leaf.
        let graph = &scope.network.graph.graph;
        let edges = nodes
            .iter()
            .map(|&node| {
                scope
                    .tree
                    .iter_preorder_tree_nodes(graph, node)
                    .flat_map(|member| graph.iter_crown(member))
                    .filter(|hedge| scope.network.graph.graph[[hedge]].is_slot())
                    .filter(|&hedge| {
                        graph.inv(hedge) == hedge || graph.node_id(graph.inv(hedge)) != node
                    })
                    .map(|hedge| hedge.0.min(graph.inv(hedge).0))
                    .collect::<BTreeSet<_>>()
            })
            .collect::<Vec<_>>();
        let mut boundary = BTreeSet::new();
        let mut planned = Vec::with_capacity(selected.len());
        while !selected.is_empty() {
            let best = (0..selected.len()).min_by_key(|&i| {
                let position = selected[i];
                (
                    boundary.symmetric_difference(&edges[position]).count(),
                    Reverse(boundary.intersection(&edges[position]).count()),
                    position,
                )
            })?;
            let position = selected.remove(best);
            boundary = boundary
                .symmetric_difference(&edges[position])
                .copied()
                .collect();
            planned.push(position);
        }
        Some(planned)
    }

    /// Compile selected, admitted arithmetic on the existing contraction tape.
    /// The outer shallow graph already supplied topology and factor order.
    /// Selected sums and products use the same tape as their parsed graph;
    /// indexed powers and transparent brackets retain the scoped graph path.
    fn compile_arithmetic<'a>(
        &self,
        value: AtomView<'a>,
        tape: &mut TermTape<InputLeaf<'a>>,
        sum: &mut ComponentSum<'a, '_>,
    ) -> Option<(usize, (usize, usize))> {
        if let AtomView::Num(number) = value {
            return if let Ok(coefficient) = Rational::try_from(value) {
                tape.number(InputLeaf::Atom(value), coefficient)
            } else if matches!(
                number.get_coeff_view(),
                CoefficientView::Natural(..) | CoefficientView::Large(..)
            ) {
                tape.leaf(InputLeaf::Scalar(value))
            } else {
                None
            };
        }
        if self
            .observations
            .as_ref()
            .and_then(|observed| observed.region(value))
            .is_some_and(|region| region.certifies_scalar_interface_syntax())
        {
            sum.scalar_interfaces.insert(value, true);
            return tape.leaf(InputLeaf::Scalar(value));
        }
        match value {
            AtomView::Add(add) => {
                let children = add
                    .iter()
                    .map(|child| self.compile_arithmetic(child, tape, sum))
                    .collect::<Option<Vec<_>>>()?;
                tape.group(children.into_iter(), true)
            }
            AtomView::Mul(product) => {
                let children = product
                    .iter()
                    .map(|child| self.compile_arithmetic(child, tape, sum))
                    .collect::<Option<Vec<_>>>()?;
                tape.group(children.into_iter(), false)
            }
            // Identity kernels may publish a terminal polynomial without
            // enumerating its scalar leaves. A new coefficient power therefore
            // need not have a regional observation yet. Reuse the existing
            // intrinsic scalar certificate and keep the whole power opaque;
            // indexed or callback-sensitive powers retain the scoped path.
            AtomView::Pow(_) if sum.scalar_spectator(value) => tape.leaf(InputLeaf::Scalar(value)),
            AtomView::Pow(_) => None,
            _ => tape.leaf(InputLeaf::Atom(value)),
        }
    }

    fn compile_scope<'a>(
        &self,
        scope: &'a ContractionScope<'_>,
        node: NodeIndex,
        tape: &mut TermTape<InputLeaf<'a>>,
        sum: &mut ComponentSum<'a, '_>,
        positions: &mut AHashMap<(usize, NodeIndex), usize>,
        scoped: &mut AHashMap<(usize, NodeIndex), &'a FactorizedContraction>,
    ) -> Option<(usize, (usize, usize))> {
        tape.compile_graph(
            &scope.network.graph,
            &scope.tree,
            node,
            &mut |tape, node, value| {
                if let Some(literal) = scope.literal(node)
                    && matches!(literal, AtomView::Add(_) | AtomView::Mul(_) | AtomView::Pow(_))
                    && !scope.contains_work(node, self, self.contraction_settings(sum.rank_one))
                {
                    // The graph already established this boundary. A completed
                    // structural region can still contain identity candidates;
                    // retain its payload without inventing terminal algebra facts.
                    sum.scalar_interfaces.insert(literal, scope.ports[&node].0.is_empty());
                    return Some(if scope.ports[&node].0.is_empty() {
                        TermLeaf::Value(InputLeaf::Scalar(literal))
                    } else {
                        let position = sum.register_factor(scope, node, positions)?;
                        TermLeaf::Value(InputLeaf::Subtree(position))
                    });
                }
                if let Some((child, power, _)) = self.power_scope(scope, node) {
                    let result = self.scoped_power(scope, node, sum.rank_one)?;
                    let unchanged = scope.literal(node).map_or_else(
                        || {
                            scope
                                .literal(child)
                                .is_some_and(|base| result.root == base.pow(power))
                        },
                        |source| result.root.as_view() == source,
                    );
                    // Visiting an unchanged scalar power is bookkeeping, not
                    // a contraction. It must not disable exact source reuse.
                    if !unchanged
                        || result.status != ReductionStatus::Complete
                    {
                        scoped.insert((scope as *const _ as usize, node), result);
                    }
                    if let Some(observed) = &result.observations {
                        sum.scalar_forms.insert(
                            result.root.as_view(),
                            observed.needs_dot_normalization(),
                        );
                    }
                    return Some(if scope.ports[&node].0.is_empty() {
                        TermLeaf::Value(InputLeaf::Scalar(result.root.as_view()))
                    } else {
                        let position = sum.register_factor(scope, node, positions)?;
                        sum.retain_factor_result(position, result);
                        TermLeaf::Value(InputLeaf::Subtree(position))
                    });
                }
                if !matches!(value, NetworkNode::Leaf(_)) {
                    return None;
                }
                let literal = scope.literal(node)?;
                sum.scalar_interfaces
                    .insert(literal, scope.ports[&node].0.is_empty());
                if matches!(
                    literal,
                    AtomView::Add(_) | AtomView::Mul(_) | AtomView::Pow(_)
                ) || matches!(literal, AtomView::Fun(function) if function.get_symbol() == self.tags.bracket) {
                    if scope.contains_work(node, self, self.contraction_settings(sum.rank_one)) {
                        let region = self.observations.as_ref()
                            .and_then(|observed| observed.region(literal));
                        if self.expand
                            && matches!(literal, AtomView::Add(_) | AtomView::Mul(_))
                            && region.as_ref()
                                .is_some_and(|region| !region.has_indexed_powers() && region.counts[7] == 0)
                        {
                            let (node, size) = self.compile_arithmetic(literal, tape, sum)?;
                            return Some(TermLeaf::Reference(node, size));
                        }
                        #[cfg(feature = "reference-cases")]
                        {
                            let reason = match region.as_ref() {
                                None => ScopeOpenReason::Unobserved,
                                Some(region) if region.contraction_sources(self.contraction_settings(sum.rank_one)) => ScopeOpenReason::MetricOrVector,
                                Some(_) => ScopeOpenReason::OpaqueOrBracket,
                            };
                            timing::scope_open(reason, region.as_ref().is_some_and(|region| region.has_indexed_powers()));
                        }
                        let opened = scope.open(node)?;
                        if !self.expand && sum.exposed_sum != Some((scope as *const _ as usize, node)) && self.tensor_alternatives(literal) {
                            // Reduce this existing scope independently; its
                            // alternatives remain inside one opaque factor of
                            // the surrounding product. Boundary bindings are
                            // applied by the same positional emitter below.
                            let result = scope.reductions[&node]
                                .get_or_init(|| self.contract_graph(opened, opened.root, None, sum.rank_one))
                                .as_ref()?;
                            if result.root.as_view() != literal || result.status != ReductionStatus::Complete {
                                scoped.insert((scope as *const _ as usize, node), result);
                            }
                            if let Some(observed) = result.observations.as_ref().filter(|observed| FactorizedContraction::scalar_facts(observed)) {
                                sum.scalar_forms.insert(result.root.as_view(), observed.needs_dot_normalization());
                            }
                            return Some(if scope.ports[&node].0.is_empty() {
                                TermLeaf::Value(InputLeaf::Scalar(result.root.as_view()))
                            } else {
                                let position = sum.register_factor(scope, node, positions)?;
                                sum.retain_factor_result(position, result);
                                TermLeaf::Value(InputLeaf::Subtree(position))
                            });
                        }
                        let (node, size) =
                            self.compile_scope(opened, opened.root, tape, sum, positions, scoped)?;
                        Some(TermLeaf::Reference(node, size))
                    } else if scope.ports[&node].0.is_empty() {
                        // Equal finished coefficients share expression identity
                        // across branches, independently of graph occurrence.
                        Some(TermLeaf::Value(InputLeaf::Scalar(literal)))
                    } else {
                        let position = sum.register_factor(scope, node, positions)?;
                        Some(TermLeaf::Value(InputLeaf::Subtree(position)))
                    }
                } else if let AtomView::Num(number) = literal {
                    if let Ok(coefficient) = Rational::try_from(literal) {
                        Some(TermLeaf::Number(InputLeaf::Atom(literal), coefficient))
                    } else if matches!(
                        number.get_coeff_view(),
                        CoefficientView::Natural(..) | CoefficientView::Large(..)
                    ) {
                        Some(TermLeaf::Value(InputLeaf::Scalar(literal)))
                    } else {
                        None
                    }
                } else {
                    Some(TermLeaf::Value(InputLeaf::Atom(literal)))
                }
            },
        )
    }

    fn contract_graph(
        &self,
        scope: &ContractionScope<'_>,
        root: NodeIndex,
        order: Option<&[usize]>,
        rank_one: bool,
    ) -> Option<FactorizedContraction> {
        if !self.expand
            && let Some(source) = scope.literal(root)
            && self
                .observations
                .as_ref()
                .and_then(|observed| observed.region(source))
                .is_some_and(|region| {
                    region.excludes_contraction_sources(self.contraction_settings(rank_one))
                })
        {
            return Some(FactorizedContraction {
                proofs: Default::default(),
                root: source.to_owned(),
                status: ReductionStatus::Complete,
                observations: self.observations.as_ref().and_then(|observed| {
                    observed.certify_terminal_region(source, &scope.ports[&root].1)
                }),
            });
        }
        if self.power_scope(scope, root).is_some() {
            return self.scoped_power(scope, root, rank_one).cloned();
        }
        let scalar_boundary = scope.ports[&root].0.is_empty();
        let root_node = root;
        if !self.expand
            && matches!(
                scope.network.graph.graph[root],
                NetworkNode::Op(NetworkOp::Sum)
            )
        {
            // Existing additive branches already own their dummy scopes. Visit
            // them independently instead of exposing them to a product tape.
            let mut terms = Vec::new();
            let mut status = ReductionStatus::Complete;
            let mut scalar_form = Some(false);
            for node in scope.tree.iter_children(root, &scope.network.graph.graph) {
                let result = self.contract_graph(scope, node, None, rank_one)?;
                status = status.max(result.status);
                scalar_form = scalar_form
                    .zip(result.observations.as_ref())
                    .filter(|(_, observed)| observed.is_terminal_polynomial())
                    .map(|(needs_dots, observed)| needs_dots || observed.needs_dot_normalization());
                terms.push(result.root);
            }
            return Some(FactorizedContraction {
                proofs: Default::default(),
                root: Atom::add_many(terms),
                status,
                observations: scalar_form
                    .filter(|_| status == ReductionStatus::Complete)
                    .map(|needs_dots| {
                        DomainObservations::terminal_polynomial(&scope.ports[&root].1, needs_dots)
                    }),
            });
        }
        let nodes = scope.factors(root);
        let selected = nodes
            .iter()
            .enumerate()
            .filter_map(|(i, &node)| {
                scope
                    .contains_work(node, self, self.contraction_settings(rank_one))
                    .then_some(i)
            })
            .collect();
        let mut selected = self.factor_order(scope, &nodes, selected, order)?;
        if !self.expand && order.is_none() {
            // Apply outer atomic sources before visiting opaque arithmetic
            // scopes. A scalar or tensor sum is never a list of alternatives
            // in the surrounding contraction state.
            selected.sort_by_key(|&position| {
                matches!(
                    scope.literal(nodes[position]),
                    Some(AtomView::Add(_) | AtomView::Mul(_) | AtomView::Pow(_))
                )
            });
        }
        let components = scope.network.graph.slot_components(&scope.tree, &nodes);
        let mut exposed_sums = AHashSet::new();
        if !self.expand {
            for component in &components {
                let alternatives = component
                    .iter()
                    .copied()
                    .filter(|&node| {
                        !scope.ports[&node].0.is_empty()
                            && scope
                                .literal(node)
                                .is_some_and(|source| self.tensor_alternatives(source))
                    })
                    .collect::<Vec<_>>();
                if let [node] = alternatives.as_slice() {
                    exposed_sums.insert(*node);
                }
            }
        }
        let ranks = selected
            .into_iter()
            .enumerate()
            .map(|(rank, position)| (nodes[position], rank))
            .collect::<AHashMap<_, _>>();
        let mut roots = Vec::with_capacity(components.len());
        let mut status = ReductionStatus::Complete;
        let mut scalar_form = Some(false);
        // Open slots of components that are scalar polynomials times vectors
        // and metrics on explicit slots, with their dot-normalization need.
        let mut terminal = Some(TerminalPorts {
            slots: Vec::new(),
            needs_dots: false,
        });
        let mut changed = false;
        let mut refused = Vec::new();
        for component in components {
            let mut selected = component
                .iter()
                .enumerate()
                .filter_map(|(position, node)| ranks.get(node).map(|&rank| (position, rank)))
                .collect::<Vec<_>>();
            selected.sort_unstable_by_key(|&(_, rank)| rank);
            let selected = selected
                .into_iter()
                .map(|(position, _)| position)
                .collect::<Vec<_>>();
            let result =
                self.contract_component(scope, &component, &selected, &exposed_sums, rank_one);
            let (result, component_changed, open_slots) = match result {
                Some(result) => result,
                None => {
                    // A refusal can follow tentative port substitutions or
                    // coefficient emission. Retain the exact source, discarding
                    // this component's scratch and tentative observations.
                    let root = Atom::mul_many(
                        component
                            .iter()
                            .map(|&node| scope.literal(node).map(|value| value.to_owned()))
                            .collect::<Option<Vec<_>>>()?,
                    );
                    refused.push(roots.len());
                    (
                        FactorizedContraction {
                            proofs: Default::default(),
                            root,
                            status: ReductionStatus::Deferred,
                            observations: None,
                        },
                        false,
                        None,
                    )
                }
            };
            changed |= component_changed;
            status = status.max(result.status);
            terminal = match (terminal, result.observations.as_ref(), open_slots) {
                (Some(mut terminal), Some(observed), None) if observed.is_closed_scalar() => {
                    terminal.needs_dots |= observed.needs_dot_normalization();
                    Some(terminal)
                }
                (Some(mut terminal), None, Some(open)) => {
                    terminal.slots.extend(open.slots);
                    terminal.needs_dots |= open.needs_dots;
                    Some(terminal)
                }
                _ => None,
            };
            scalar_form = scalar_form
                .zip(result.observations.as_ref())
                .map(|(needs_dots, observed)| needs_dots || observed.needs_dot_normalization());
            roots.push(result.root);
        }
        if !refused.is_empty() {
            if !changed {
                // Preserve the checked whole-source fallback when factorized
                // planning made no independent progress. In particular, an
                // inexact coefficient must not become an admitted polynomial.
                return None;
            }
            for position in refused {
                // Independent work is already exact. Rewrite only the refused
                // components, with the same dummy owner as the whole domain.
                roots[position] = self.contract_ordered(
                    roots[position].as_view(),
                    rank_one,
                    self.observations
                        .as_ref()
                        .map(|observed| &observed.candidates),
                );
            }
        }
        let retained_source = status == ReductionStatus::Complete && !changed;
        let root = if retained_source {
            let source = scope.literal(root)?;
            // The bracket normalizer owns whether a product can be exposed:
            // closed scalar wrappers disappear, while AUTO operand order stays.
            if matches!(source, AtomView::Fun(function)
                if function.get_symbol() == self.tags.bracket)
            {
                BracketNormalizer::with_matcher(source, self.normalization_matcher())
            } else {
                source.to_owned()
            }
        } else {
            Atom::mul_many(roots)
        };
        let observations = if retained_source {
            // Keeping the source's factorization also keeps its syntax facts.
            // A tentative coefficient sum may have cancelled only some indexed
            // terms; its terminal certificate cannot describe the old payload.
            self.observations.as_ref().and_then(|observed| {
                observed.certify_terminal_region(root.as_view(), &scope.ports[&root_node].1)
            })
        } else if scalar_boundary {
            scalar_form
                .filter(|_| status == ReductionStatus::Complete)
                .map(DomainObservations::closed_scalar)
        } else {
            // Every component is a certified scalar polynomial, and the open
            // vectors and metrics sit exactly on the root's ports: no internal
            // connection and no identity candidate remains. Certify the open
            // output as terminal, so that its replacement need not be rescanned.
            terminal
                .filter(|_| status == ReductionStatus::Complete)
                .and_then(|mut terminal| {
                    terminal.slots.sort();
                    let mut ports = scope.ports[&root_node].1.slots().ok()?;
                    ports.sort();
                    (terminal.slots == ports).then(|| {
                        DomainObservations::terminal_polynomial(
                            &scope.ports[&root_node].1,
                            terminal.needs_dots,
                        )
                    })
                })
        };
        Some(FactorizedContraction {
            proofs: Default::default(),
            root,
            status,
            observations,
        })
    }

    /// Execute one connected component on the admitted graph. Scratch bindings,
    /// compiled factors and tentative proofs cannot escape a refused component.
    fn contract_component<'a>(
        &self,
        scope: &'a ContractionScope<'_>,
        nodes: &[NodeIndex],
        selected: &[usize],
        exposed_sums: &AHashSet<NodeIndex>,
        rank_one: bool,
    ) -> Option<(FactorizedContraction, bool, Option<TerminalPorts>)> {
        let mut slots = SlotMatcher::default();
        let mut sum = ComponentSum::new(self, Intake::Contraction, &mut slots);
        sum.rank_one = rank_one;
        let mut positions = AHashMap::new();
        let mut remaining = Vec::with_capacity(nodes.len());
        for &node in nodes {
            let position = sum.register_factor(scope, node, &mut positions)?;
            remaining.push((
                position,
                sum.opaque_factors[position]
                    .1
                    .iter()
                    .copied()
                    .map(Argument::Original)
                    .collect::<Vec<_>>(),
            ));
        }
        if selected.is_empty() {
            // Untouched components need neither a polynomial frontier nor a
            // unit coefficient. Keep the admitted scalar facts and exact leaves.
            let mut scalar_form = Some(false);
            let mut terminal = Some(TerminalPorts {
                slots: Vec::new(),
                needs_dots: false,
            });
            let mut factors = Vec::with_capacity(remaining.len());
            for (&node, (position, arguments)) in nodes.iter().zip(remaining) {
                terminal = terminal.and_then(|mut terminal| {
                    let source = sum.opaque_factors[position].0;
                    let interface = &scope.ports[&node].1;
                    let observed = self
                        .observations
                        .as_ref()?
                        .certify_terminal_region(source, interface)?;
                    terminal.slots.extend(interface.slots().ok()?);
                    terminal.needs_dots |= observed.needs_dot_normalization();
                    Some(terminal)
                });
                scalar_form = scalar_form
                    .zip(sum.scalar_form(sum.opaque_factors[position].0))
                    .filter(|_| sum.opaque_factors[position].1.is_empty())
                    .map(|(first, second)| first || second);
                factors.push(sum.emit_factor(position, &arguments)?);
            }
            return Some((
                FactorizedContraction {
                    proofs: Default::default(),
                    root: Atom::mul_many(factors),
                    status: ReductionStatus::Complete,
                    observations: scalar_form.map(DomainObservations::closed_scalar),
                },
                false,
                terminal.filter(|terminal| !terminal.slots.is_empty()),
            ));
        }
        let mut scoped = AHashMap::new();
        let mut compile = |sum: &mut ComponentSum<'a, '_>, position: usize| {
            let mut tape = mem::take(&mut sum.input);
            sum.exposed_sum = exposed_sums
                .contains(&nodes[position])
                .then_some((scope as *const _ as usize, nodes[position]));
            let result = self.compile_scope(
                scope,
                nodes[position],
                &mut tape,
                sum,
                &mut positions,
                &mut scoped,
            );
            sum.exposed_sum = None;
            sum.input = tape;
            result.map(|(node, _)| node)
        };
        let mut result = sum.contract_states(
            Remaining {
                factors: remaining,
                residual: Vec::new(),
            },
            selected,
            &mut compile,
        )?;
        for replacement in scoped.values() {
            result.status = result.status.max(replacement.status);
        }
        // Coefficient cancellation can empty the frontier without contracting
        // a port. Publish that zero instead of pairing its facts with the old
        // source syntax, which can still contain internal indices.
        let changed = sum.contracted || !scoped.is_empty() || result.root.is_zero();
        Some((result, changed, sum.terminal_slots.take()))
    }

    /// The ordered rewrite is shared by refused intrinsic components and the
    /// checked callback path. Keep existing additive branches independently scoped.
    pub(crate) fn contract_ordered<const N: usize>(
        &self,
        source: AtomView<'_>,
        rank_one: bool,
        observed: Option<&SimplificationCandidates<N>>,
    ) -> Atom {
        let mut settings = SchoonschipSettings::default().with_chain_like_functions();
        settings.schoonschip_rank1_tensors = rank_one;
        settings.metrics = self.metrics;
        settings.representations = self.representations.clone();
        let rewrite = SchoonschipWithSettings {
            settings: &settings,
        };
        let rewrite_term = |term| {
            let prepared = self.contract_leaf_metadata(term, &rewrite);
            rewrite.run_observed(prepared.as_view(), observed, self)
        };
        match source {
            AtomView::Add(sum) => Atom::add_many(sum.iter().map(rewrite_term)),
            _ => rewrite_term(source),
        }
    }
}

impl<'a> ComponentSum<'a, '_> {
    /// Storage of one [`CoefficientSum::frontier_sum`]: its contribution records and the
    /// output it reserves at the input term count over the shared map.
    #[cfg(feature = "reference-cases")]
    fn merge_bytes(
        sources: &[Option<ScalarPolynomial>],
        contributions: &[Contribution],
        nvars: usize,
    ) -> usize {
        let terms = contributions
            .iter()
            .filter_map(|contribution| sources[contribution.source].as_ref())
            .map(ScalarPolynomial::nterms)
            .fold(0usize, usize::saturating_add);
        let records = contributions
            .iter()
            .map(|contribution| {
                mem::size_of::<Contribution>().saturating_add(
                    contribution
                        .shift
                        .capacity()
                        .saturating_mul(mem::size_of::<u32>()),
                )
            })
            .fold(0usize, usize::saturating_add);
        terms
            .saturating_mul(
                mem::size_of::<Rational>()
                    .saturating_add(nvars.saturating_mul(mem::size_of::<u32>())),
            )
            .saturating_add(records)
    }

    /// The explicit slots of a vector or metric variable, whose only open
    /// indices are these slots; `None` for any other variable.
    fn terminal_variable(&mut self, position: usize) -> Option<&TerminalPorts> {
        if !self.terminal_variables.contains_key(&position) {
            let certificate = self.certify_terminal_variable(position);
            self.terminal_variables.insert(position, certificate);
        }
        self.terminal_variables[&position].as_ref()
    }

    fn certify_terminal_variable(&mut self, position: usize) -> Option<TerminalPorts> {
        #[cfg(test)]
        tests::TERMINAL_VARIABLE_OBSERVATIONS.with(|count| count.set(count.get() + 1));
        match self.variables[position] {
            Variable::Vector(vector, slot) => {
                let source = self.vectors[vector];
                let needs_dots = source.iter().take(source.get_nargs() - 1).try_fold(
                    false,
                    |needs_dots, argument| {
                        self.scalar_form(argument).map(|next| needs_dots || next)
                    },
                )?;
                Some(TerminalPorts {
                    slots: vec![self.slots.parse(slot).ok()?],
                    needs_dots,
                })
            }
            Variable::Metric([first, second]) => Some(TerminalPorts {
                slots: vec![
                    self.slots.parse(first).ok()?,
                    self.slots.parse(second).ok()?,
                ],
                needs_dots: false,
            }),
            _ => None,
        }
    }

    /// Reuse the admitted inventory; never inspect a completed coefficient to
    /// rediscover whether its source contains an identity or an index.
    fn scalar_form(&self, value: AtomView<'a>) -> Option<bool> {
        if let Some(&needs_dots) = self.scalar_forms.get(&value) {
            return Some(needs_dots);
        }
        if matches!(value, AtomView::Num(_)) {
            return Some(false);
        }
        if let AtomView::Var(variable) = value {
            let head = variable.get_symbol();
            let tags = self.contractor.tags;
            if head.get_wildcard_level() == 0
                && !head.has_tag(&tags.tensor)
                && !head.has_tag(&tags.rank1)
                && !head.has_tag(&tags.representation)
                && head != tags.chain_in
                && head != tags.chain_out
            {
                return Some(false);
            }
        }
        let observed = self.contractor.observations.as_ref();
        let known = observed.is_some_and(|observed| observed.region(value).is_some());
        if known {
            observed
                .and_then(|observed| observed.certify_scalar_region(value))
                .map(|observed| observed.needs_dot_normalization())
        } else {
            // The shallow parser can compose scalar coefficients that were
            // separate factors at admission. Combine only these new arithmetic
            // connectors, stopping at every already-observed subtree.
            match value {
                AtomView::Add(sum) => sum.iter().try_fold(false, |needs_dots, term| {
                    self.scalar_form(term).map(|next| needs_dots || next)
                }),
                AtomView::Mul(product) => product.iter().try_fold(false, |needs_dots, factor| {
                    self.scalar_form(factor).map(|next| needs_dots || next)
                }),
                AtomView::Pow(power) => {
                    let (base, exponent) = power.get_base_exp();
                    self.scalar_form(base)
                        .zip(self.scalar_form(exponent))
                        .map(|(base, exponent)| base || exponent)
                }
                _ => None,
            }
        }
    }

    fn variable_scalar_form(&mut self, variable: &Variable<'a>) -> Option<bool> {
        match variable {
            Variable::Scalar(value @ AtomView::Fun(function))
                if self.scalar_form(*value).is_none() =>
            {
                // A compact dot can come from the shallow parser's scalar
                // store, after notation normalization changed its source key.
                // Reuse the same admitted vector/space parser as dot emission.
                let (space, first, second) = self.compact_dot(*function)?;
                let needs_dots =
                    self.variable_scalar_form(&Variable::Dot(space, [first, second]))?;
                self.scalar_forms.insert(*value, needs_dots);
                Some(needs_dots)
            }
            Variable::Scalar(value) => self.scalar_form(*value),
            Variable::Subtree(position) if self.opaque_factors[*position].1.is_empty() => {
                self.scalar_form(self.opaque_factors[*position].0)
            }
            Variable::Dot(_, vectors) => {
                // Vector admission already proved intrinsic heads, compatible
                // atomic dimensions and scalar metadata. Require the stronger
                // identity-free certificate before skipping future family work.
                for &vector in vectors {
                    let source = self.vectors[vector];
                    for argument in source.iter().take(source.get_nargs() - 1) {
                        self.scalar_form(argument)?;
                    }
                }
                Some(true)
            }
            _ => None,
        }
    }

    fn register_factor(
        &mut self,
        scope: &'a ContractionScope<'_>,
        node: NodeIndex,
        positions: &mut AHashMap<(usize, NodeIndex), usize>,
    ) -> Option<usize> {
        let key = (scope as *const _ as usize, node);
        if let Some(&position) = positions.get(&key) {
            return Some(position);
        }
        let position = self.opaque_factors.len();
        let (ports, interface) = &scope.ports[&node];
        self.scalar_interfaces
            .insert(scope.literal(node)?, ports.is_empty());
        self.opaque_factors.push((
            scope.literal(node)?,
            ports.iter().map(Atom::as_view).collect(),
            interface.clone(),
        ));
        if let Some(observed) =
            self.contractor.observations.as_ref().and_then(|observed| {
                observed.certify_terminal_region(scope.literal(node)?, interface)
            })
            && let Ok(slots) = interface.slots()
        {
            self.terminal_factors.insert(
                position,
                TerminalPorts {
                    slots,
                    needs_dots: observed.needs_dot_normalization(),
                },
            );
        }
        self.factor_roots.push(None);
        positions.insert(key, position);
        Some(position)
    }

    fn retain_factor_result(&mut self, position: usize, result: &'a FactorizedContraction) {
        self.opaque_factors[position].0 = result.root.as_view();
        self.terminal_factors.remove(&position);
        if let Some(observed) = result
            .observations
            .as_ref()
            .filter(|observed| observed.is_terminal_polynomial())
            && let Ok(slots) = self.opaque_factors[position].2.slots()
        {
            self.terminal_factors.insert(
                position,
                TerminalPorts {
                    slots,
                    needs_dots: observed.needs_dot_normalization(),
                },
            );
        }
    }

    fn opaque_factor(&mut self, position: usize, arguments: &[Argument<'a>]) -> Option<()> {
        let tensor = self.tensors.len();
        self.tensors
            .push((TensorSource::Factor(position), arguments.to_vec()));
        for (port, &argument) in arguments.iter().enumerate() {
            if let Argument::Original(slot) = argument {
                if !self.permits_slot(slot) {
                    continue;
                }
                let (space, index) = self.resolve_endpoint(slot)?;
                let node = self.endpoint(slot, space, index)?;
                self.nodes[node]
                    .terminals
                    .push(Terminal::Tensor(tensor, port));
            }
        }
        Some(())
    }

    fn residual_factor(&mut self, variable: usize, exponent: u32) -> Option<()> {
        match self.variables[variable].clone() {
            Variable::Vector(head, slot) => {
                let (space, index) = self.resolve_endpoint(slot)?;
                if exponent > 1 {
                    self.dot(space, head, head, exponent / 2);
                }
                if exponent % 2 == 1 {
                    let node = self.endpoint(slot, space, index)?;
                    self.nodes[node].terminals.push(Terminal::Vector(head));
                }
            }
            Variable::Metric([first, second]) => {
                self.metric(
                    Argument::Original(first),
                    Argument::Original(second),
                    exponent,
                )?;
            }
            Variable::Tensor(source, arguments) => {
                if exponent != 1 {
                    return None;
                }
                self.tensor_with_arguments(source, arguments)?;
            }
            value => self.variable(value, exponent),
        }
        Some(())
    }

    fn reset_components(&mut self, coefficient: &Rational) {
        self.nodes.clear();
        self.occurrences.clear();
        self.factors.clear();
        self.tensors.clear();
        self.metrics = 0;
        self.coefficient = coefficient.clone();
    }

    fn factor_terms(&mut self, position: usize) -> Option<Vec<FactorTerm>> {
        let node = self.factor_roots.get(position).copied().flatten()?;
        self.input.clear_terms();
        self.input
            .distribute(&mut vec![(node, 1)], &Rational::one())?;
        Some(self.input.take_terms())
    }

    fn contract_states(
        &mut self,
        remaining: Remaining<'a>,
        order: &[usize],
        compile: &mut impl FnMut(&mut Self, usize) -> Option<usize>,
    ) -> Option<FactorizedContraction> {
        #[cfg(feature = "reference-cases")]
        let _states_timing = timing::scope(Phase::ContractionStates);
        let mut states = vec![(remaining, ScalarPolynomial::new_one(&Q))];
        let mut emitted: AHashMap<usize, Atom> = AHashMap::new();
        #[cfg(feature = "reference-cases")]
        let mut emitted_bytes = 0usize;
        let mut scalar_variables: AHashMap<Atom, usize> = AHashMap::new();
        let mut shared = SharedScalars::new(self.contractor);
        let mut scalar_form = Some(false);
        #[cfg(feature = "reference-cases")]
        let mut generated = 0usize;
        for &selected in order {
            if self.factor_roots[selected].is_none() {
                self.factor_roots[selected] = Some(compile(self, selected)?);
            }
            let terms = self.factor_terms(selected)?;
            // Measure generated coefficient storage without placing an
            // automatic work limit on the selected contraction.
            #[cfg(feature = "reference-cases")]
            let state_bytes = states
                .iter()
                .map(|(remaining, coefficient)| {
                    Self::coefficient_bytes(coefficient) + remaining.storage_bytes()
                })
                .fold(0usize, usize::saturating_add);
            #[cfg(feature = "reference-cases")]
            let count = states
                .iter()
                .map(|(_, weight)| weight.nterms())
                .fold(0usize, usize::saturating_add)
                .saturating_mul(terms.len());
            #[cfg(feature = "reference-cases")]
            let predicted = state_bytes.saturating_add(emitted_bytes.saturating_mul(2));
            #[cfg(feature = "reference-cases")]
            timing::frontier(Frontier {
                selected,
                states: states.len(),
                local_terms: terms.len(),
                generated,
                predicted_bytes: predicted,
                coefficient_bytes: states
                    .iter()
                    .map(|(_, coefficient)| Self::coefficient_bytes(coefficient))
                    .sum(),
                max_coefficient_bytes: states
                    .iter()
                    .map(|(_, coefficient)| Self::coefficient_bytes(coefficient))
                    .max()
                    .unwrap_or(0),
            });
            #[cfg(feature = "reference-cases")]
            {
                generated = generated.saturating_add(count);
            }
            let mut positions: AHashMap<Remaining<'a>, usize> = AHashMap::new();
            // Every target state collects the scaled, shifted large source
            // states that reach it and the sum of its small products; one
            // merge per target then builds its coefficient.
            let mut next: Vec<(Remaining<'a>, Vec<Contribution>, Option<CoefficientSum>)> =
                Vec::new();
            let mut uses = vec![0usize; states.len()];
            for (origin, (remaining, weight)) in states.iter().enumerate() {
                let degrees = CoefficientSum::max_degrees(weight);
                let source = remaining
                    .factors
                    .iter()
                    .find(|(position, _)| *position == selected)?;
                self.overrides = self.opaque_factors[selected]
                    .1
                    .iter()
                    .copied()
                    .zip(source.1.iter().copied())
                    .collect();
                for (factors, coefficient) in &terms {
                    self.reset_components(coefficient);
                    for &(factor, exponent) in factors {
                        let leaf = self.input.leaves()[factor];
                        match leaf {
                            InputLeaf::Scalar(value) => {
                                self.variable(Variable::Scalar(value), exponent);
                            }
                            InputLeaf::Atom(value) if self.scalar_spectator(value) => {
                                self.variable(Variable::Scalar(value), exponent);
                            }
                            InputLeaf::Atom(value) => self.factor(value, exponent)?,
                            InputLeaf::Subtree(position) => {
                                let ports = &self.opaque_factors[position].1;
                                if ports.is_empty() {
                                    self.variable(Variable::Subtree(position), exponent);
                                } else {
                                    if exponent != 1 {
                                        return None;
                                    }
                                    let arguments = ports
                                        .iter()
                                        .map(|&slot| self.argument(slot))
                                        .collect::<Vec<_>>();
                                    self.opaque_factor(position, &arguments)?;
                                }
                            }
                        }
                    }
                    for &(variable, exponent) in &remaining.residual {
                        self.residual_factor(variable, exponent)?;
                    }
                    for (position, arguments) in &remaining.factors {
                        if *position != selected {
                            self.opaque_factor(*position, arguments)?;
                        }
                    }
                    self.reduce_components()?;
                    let mut key = Remaining {
                        factors: Vec::new(),
                        residual: Vec::new(),
                    };
                    let mut variables = Vec::new();
                    let mut exponents: Vec<u32> = Vec::new();
                    for factor in 0..self.factors.len() {
                        let (variable, exponent) = self.factors[factor];
                        match self.variables[variable].clone() {
                            Variable::Tensor(TensorSource::Factor(position), arguments) => {
                                key.factors.push((position, arguments));
                            }
                            value @ (Variable::Scalar(_)
                            | Variable::Subtree(_)
                            | Variable::Dot(_, _)) => {
                                let atom = if let Some(value) = emitted.get(&variable) {
                                    value.clone()
                                } else {
                                    // Dot-notation emission discharges dot work.
                                    let dot_notation = self.contractor.dot_notation;
                                    scalar_form = scalar_form
                                        .zip(self.variable_scalar_form(&value))
                                        .map(|(first, second)| (first, second && !dot_notation))
                                        .map(|(first, second)| first || second);
                                    let atom = self.emit_variable(&value)?;
                                    #[cfg(feature = "reference-cases")]
                                    {
                                        emitted_bytes = emitted_bytes
                                            .saturating_add(atom.as_view().get_byte_size());
                                    }
                                    emitted.insert(variable, atom.clone());
                                    atom
                                };
                                let index = *scalar_variables.entry(atom).or_insert(variable);
                                let variable = PolyVariable::Temporary(index);
                                if let Some(position) =
                                    variables.iter().position(|entry| entry == &variable)
                                {
                                    exponents[position] =
                                        exponents[position].checked_add(exponent)?;
                                } else {
                                    variables.push(variable);
                                    exponents.push(exponent);
                                }
                            }
                            _ => key.residual.push((variable, exponent)),
                        }
                    }
                    key.factors.sort_unstable_by_key(|(position, _)| *position);
                    key.residual.sort_unstable();
                    let mut shared_exponents = Vec::with_capacity(variables.len());
                    let mut composites = Vec::new();
                    for (variable, exponent) in variables.iter().zip(&exponents) {
                        let PolyVariable::Temporary(index) = variable else {
                            unreachable!("frontier variables are registered by temporary position")
                        };
                        let atom = emitted.get(index)?;
                        let slots = &mut *self.slots;
                        let composite = match shared.expansion(*index, atom, |atom| {
                            SharedScalars::index_free(slots, atom.as_view())
                        }) {
                            Expansion::Base(position) => Err(*position),
                            Expansion::Composite(expanded) => {
                                Ok(SharedScalars::small_power(expanded.nterms(), *exponent))
                            }
                        };
                        match composite {
                            Err(position) => shared_exponents.push((position, *exponent)),
                            Ok(true) => composites.push((*index, *exponent)),
                            // A large power of an expanded composite stays one
                            // opaque variable in this term.
                            Ok(false) => {
                                shared_exponents.push((shared.position(atom.clone()), *exponent));
                            }
                        }
                    }
                    // Distinct frontier variables can share one position (a
                    // nested `dot(p,q)` and this frontier's `g(p,q)`), so
                    // their exponents add.
                    let mut full = vec![0u32; shared.map.len()];
                    for (position, exponent) in shared_exponents {
                        full[position] = full[position].checked_add(exponent)?;
                    }
                    for (position, &exponent) in full.iter().enumerate().take(weight.nvars()) {
                        // Machine exponent capacity is not a work budget.
                        // A decline preserves the caller's exact source.
                        if exponent != 0 {
                            degrees[position].checked_add(exponent)?;
                        }
                    }
                    let mut monomial = ScalarPolynomial::new(&Q, Some(1), Arc::clone(&shared.map))
                        .monomial(self.coefficient.clone(), full);
                    for (index, exponent) in composites {
                        let Expansion::Composite(expanded) = &shared.expansions[&index] else {
                            unreachable!("only cached composite expansions enter this list")
                        };
                        for _ in 0..exponent {
                            CoefficientSum::check_product_degrees(
                                &CoefficientSum::max_degrees(&monomial),
                                expanded,
                            )?;
                            monomial = &monomial * expanded;
                        }
                    }
                    // Composite expansion adds degrees after the initial
                    // monomial was assembled. Validate the final contribution
                    // before either eager multiplication or lazy shifting.
                    CoefficientSum::check_product_degrees(&degrees, &monomial)?;
                    if monomial.is_zero() || weight.is_zero() {
                        continue;
                    }
                    // Contributions address positions of the shared map.
                    if monomial.get_vars_ref() != &shared.map[..monomial.nvars()] {
                        return None;
                    }
                    let target = match positions.get(&key) {
                        Some(&target) => target,
                        None => {
                            positions.insert(key.clone(), next.len());
                            next.push((key, Vec::new(), None));
                            next.len() - 1
                        }
                    };
                    if weight.nterms() < LAZY_SOURCE_TERMS {
                        // A small product is cheaper to form and fold now
                        // than to describe for the merge.
                        let product = weight * &monomial;
                        match &mut next[target].2 {
                            Some(sum) => sum.push(
                                product,
                                #[cfg(feature = "reference-cases")]
                                state_bytes.saturating_add(emitted_bytes.saturating_mul(2)),
                            ),
                            empty => *empty = Some(CoefficientSum::One(product)),
                        }
                        continue;
                    }
                    uses[origin] += monomial.nterms();
                    for term in 0..monomial.nterms() {
                        next[target].1.push(Contribution {
                            source: origin,
                            coefficient: monomial.coefficients[term].clone(),
                            shift: monomial.exponents(term).to_vec(),
                        });
                    }
                }
            }
            self.overrides.clear();
            let mut sources: Vec<Option<ScalarPolynomial>> = mem::take(&mut states)
                .into_iter()
                .map(|(_, weight)| Some(weight))
                .collect();
            #[cfg(feature = "reference-cases")]
            let mut next_bytes = 0usize;
            for (remaining, mut contributions, small) in next {
                // Released sources only lower the retained input storage, so
                // the step's input bytes bound it.
                #[cfg(feature = "reference-cases")]
                let surrounding_bytes = state_bytes
                    .saturating_add(next_bytes)
                    .saturating_add(emitted_bytes.saturating_mul(2));
                let small = small.map(|sum| {
                    sum.finish(
                        #[cfg(feature = "reference-cases")]
                        surrounding_bytes,
                    )
                });
                let weight = match small {
                    Some(small) if contributions.is_empty() => small,
                    small => {
                        if let Some(small) = small {
                            // The folded small products are one more stream.
                            contributions.push(Contribution {
                                source: sources.len(),
                                coefficient: Rational::one(),
                                shift: Vec::new(),
                            });
                            sources.push(Some(small));
                            uses.push(1);
                        }
                        #[cfg(feature = "reference-cases")]
                        timing::frontier_workspace(surrounding_bytes.saturating_add(
                            Self::merge_bytes(&sources, &contributions, shared.map.len()),
                        ));
                        CoefficientSum::frontier_sum(
                            &mut sources,
                            &mut uses,
                            &contributions,
                            &shared.map,
                        )?
                    }
                };
                if !weight.is_zero() {
                    #[cfg(feature = "reference-cases")]
                    {
                        next_bytes = next_bytes.saturating_add(
                            Self::coefficient_bytes(&weight)
                                .saturating_add(remaining.storage_bytes()),
                        );
                    }
                    states.push((remaining, weight));
                }
            }
        }
        self.overrides.clear();
        #[cfg(feature = "reference-cases")]
        let _emission_timing = timing::scope(Phase::ContractionEmission);
        let mut terms = Vec::with_capacity(states.len());
        // Open slots of vectors and metrics left in every state: the same
        // sorted explicit slots in each state, or no terminal certificate.
        let mut open_slots = None;
        let mut any_open = false;
        let mut terminal = true;
        for (remaining, weight) in states {
            let mut factors = vec![shared.emit(weight)?];
            let mut state_slots = Vec::new();
            for (position, arguments) in remaining.factors {
                let ports = &self.opaque_factors[position].1;
                if let Some(facts) = self.terminal_factors.get(&position).filter(|_| {
                    ports.len() == arguments.len()
                        && ports
                            .iter()
                            .zip(&arguments)
                            .all(|(&port, &argument)| argument == Argument::Original(port))
                }) {
                    state_slots.extend_from_slice(&facts.slots);
                    scalar_form = scalar_form.map(|needs_dots| needs_dots || facts.needs_dots);
                } else {
                    scalar_form = scalar_form
                        .zip(self.scalar_form(self.opaque_factors[position].0))
                        .filter(|_| ports.is_empty())
                        .map(|(first, second)| first || second);
                }
                factors.push(self.emit_factor(position, &arguments)?);
            }
            for (variable, exponent) in remaining.residual {
                if let Some(ports) = self.terminal_variable(variable) {
                    // Each power repeats the slots, so a repeated (internal)
                    // index never matches the distinct root ports.
                    scalar_form = scalar_form.map(|needs_dots| needs_dots || ports.needs_dots);
                    if exponent == 1 {
                        state_slots.extend_from_slice(&ports.slots);
                    } else {
                        // A repeated port cannot certify an open interface. Do
                        // not allocate one slot per power merely to reject it.
                        terminal = false;
                        scalar_form = None;
                    }
                    let atom = match emitted.get(&variable) {
                        Some(value) => value.clone(),
                        None => {
                            let atom = self.emit_variable(&self.variables[variable])?;
                            emitted.insert(variable, atom.clone());
                            atom
                        }
                    };
                    factors.push(atom.pow(exponent));
                    continue;
                }
                let atom = if let Some(value) = emitted.get(&variable) {
                    value.clone()
                } else {
                    scalar_form = scalar_form
                        .zip(self.variable_scalar_form(&self.variables[variable].clone()))
                        .map(|(first, second)| first || second);
                    let atom = self.emit_variable(&self.variables[variable])?;
                    emitted.insert(variable, atom.clone());
                    atom
                };
                factors.push(atom.pow(exponent));
            }
            state_slots.sort();
            any_open |= !state_slots.is_empty();
            match &open_slots {
                None => open_slots = Some(state_slots),
                Some(slots) => terminal &= *slots == state_slots,
            }
            terms.push(Atom::mul_many(factors));
        }
        self.terminal_slots = match (open_slots, scalar_form) {
            (Some(slots), Some(needs_dots)) if any_open && terminal => {
                Some(TerminalPorts { slots, needs_dots })
            }
            _ => None,
        };
        Some(FactorizedContraction {
            proofs: Default::default(),
            root: Atom::add_many(terms),
            status: ReductionStatus::Complete,
            // Open slots leave a tensor, never a closed scalar.
            observations: scalar_form
                .filter(|_| !any_open)
                .map(DomainObservations::closed_scalar),
        })
    }

    pub(super) fn argument(&self, value: AtomView<'a>) -> Argument<'a> {
        self.overrides
            .iter()
            .find_map(|&(source, target)| (source == value).then_some(target))
            .unwrap_or(Argument::Original(value))
    }

    pub(super) fn emit_argument(&self, argument: Argument<'a>) -> Atom {
        match argument {
            Argument::Original(value) => value.to_owned(),
            Argument::Vector(space, head) => {
                let (representation, dimension) = self.spaces[space];
                self.emit_vector(head, representation.to_symbolic([dimension]).as_view())
            }
        }
    }

    pub(super) fn emit_factor(&self, position: usize, arguments: &[Argument<'a>]) -> Option<Atom> {
        let (source, ports, interface) = &self.opaque_factors[position];
        let replacements = ports
            .iter()
            .zip(arguments)
            .enumerate()
            .filter(|(_, (port, argument))| **argument != Argument::Original(**port))
            .map(|(position, (_, &argument))| (position, self.emit_argument(argument)))
            .collect::<HashMap<_, _>>();
        if replacements.is_empty() {
            return Some(source.to_owned());
        }
        let source = self
            .contractor
            .proofs
            .as_ref()?
            .region_with_interface(*source, interface.clone())?;
        let (expression, produced) = composition::rewrite_interface_parts(
            &source,
            &replacements,
            Some(self.contractor.dummy_state()),
        )
        .ok()?;
        self.contractor.generated.borrow_mut().push(produced);
        Some(DotNormalizer::with_matcher(
            expression.as_view(),
            ContractSettings {
                metrics: self.contractor.metrics,
                rank_one: self.rank_one,
                representations: self.contractor.representations.as_deref(),
                ..Default::default()
            },
            self.contractor.normalization_matcher(),
        ))
    }

    pub(super) fn metric(
        &mut self,
        first: Argument<'a>,
        second: Argument<'a>,
        exponent: u32,
    ) -> Option<()> {
        match (first, second) {
            (Argument::Vector(space, a), Argument::Vector(other, b)) => {
                if space != other {
                    return None;
                }
                self.dot(space, a, b, exponent);
            }
            (Argument::Vector(space, head), Argument::Original(slot))
            | (Argument::Original(slot), Argument::Vector(space, head)) => {
                let (other, index) = self.resolve_endpoint(slot)?;
                if space != other {
                    return None;
                }
                if exponent > 1 {
                    self.dot(space, head, head, exponent / 2);
                    self.contracted = true;
                }
                if exponent % 2 == 1 {
                    let node = self.endpoint(slot, space, index)?;
                    self.nodes[node].terminals.push(Terminal::Vector(head));
                }
            }
            (Argument::Original(first), Argument::Original(second)) => {
                let (space, first_index) = self.resolve_endpoint(first)?;
                let (other, second_index) = self.resolve_endpoint(second)?;
                if space != other {
                    return None;
                }
                if exponent > 1 {
                    if first_index == second_index {
                        return None;
                    }
                    self.factor(self.spaces[space].1, exponent / 2)?;
                    self.contracted = true;
                }
                if exponent % 2 == 1 {
                    let a = self.endpoint(first, space, first_index)?;
                    let b = self.endpoint(second, space, second_index)?;
                    let a = self.root(a);
                    let b = self.root(b);
                    self.metrics += 1;
                    self.nodes[a].metrics += 1;
                    self.nodes[b].parent = a;
                }
            }
        }
        Some(())
    }
}

#[cfg(test)]
mod tests {
    use super::super::tests::{input, setup};
    use super::*;
    use crate::shorthands::schoonschip::Schoonschip;
    use crate::tensor::SymbolicNetParse;
    use crate::tensor::inference::tests::{INFERENCE_CALLS, INFERENCE_VISITS};
    use crate::test_support;
    use crate::{color_f, epsilon};
    use spenso::{g, mink, p, q, vector_symbol};
    use std::cell::Cell;
    use symbolica::{function, parse_lit};

    thread_local! {
        pub(super) static TERMINAL_VARIABLE_OBSERVATIONS: Cell<usize> = const { Cell::new(0) };
    }

    #[test]
    fn expanded_composite_degrees_decline_before_frontier_overflow() {
        test_support::test_initialize();
        let x = parse_lit!(frontier_overflow_x);
        let y = parse_lit!(frontier_overflow_y);
        let contractor = SlotContraction::new().expanding_scalars(true);
        let mut shared = SharedScalars::new(&contractor);
        shared.position(x.clone());
        shared.position(y);
        let mut source =
            ScalarPolynomial::new(&Q, Some(LAZY_SOURCE_TERMS), Arc::clone(&shared.map));
        for degree in 0..LAZY_SOURCE_TERMS {
            source.append_monomial_back(Rational::one(), &[u32::MAX - 1, degree as u32]);
        }
        let composite = x.pow(2) + Atom::one();
        let Expansion::Composite(expanded) = shared.expansion(7, &composite, |_| true) else {
            panic!("the small index-free spectator must expand");
        };
        let degrees = CoefficientSum::max_degrees(&source);
        assert!(CoefficientSum::check_product_degrees(&degrees, expanded).is_none());

        let mut sources = vec![Some(source)];
        let contribution = Contribution {
            source: 0,
            coefficient: Rational::one(),
            shift: vec![2, 0],
        };
        assert!(
            CoefficientSum::frontier_sum(&mut sources, &mut [1], &[contribution], &shared.map)
                .is_none()
        );
        assert!(
            sources[0].is_some(),
            "overflow must not consume the exact source"
        );
        let contribution = Contribution {
            source: 0,
            coefficient: Rational::one(),
            shift: vec![1, 0],
        };
        let shifted =
            CoefficientSum::frontier_sum(&mut sources, &mut [1], &[contribution], &shared.map)
                .unwrap();
        assert_eq!(shifted.nterms(), LAZY_SOURCE_TERMS);
        assert_eq!(CoefficientSum::max_degrees(&shifted)[0], u32::MAX);
    }

    #[test]
    fn frontier_retains_composites_with_unwriteable_or_overflowing_powers() {
        let contractor = setup().expanding_scalars(true);
        for spectator in [
            "(x^1073741824+y)^2",
            "(x^1073741824+y)^4",
            "(x^4294967296+y)^2",
        ] {
            let source = input(&format!(
                "(({spectator})*p(mink(4,a))+q(mink(4,a)))*q(mink(4,a))"
            ));
            let result = contractor
                .contract_factorized(source.as_view(), None, true)
                .unwrap();
            let expected = input(&format!(
                "({spectator})*g(p(mink(4)),q(mink(4)))+g(q(mink(4)),q(mink(4)))"
            ));
            assert_eq!(result.root, expected, "{source}");
        }
        let source = input("(x^2147483645+y)*(x^2+z)*p(mink(4,a))*q(mink(4,a))");
        let result = contractor
            .contract_factorized(source.as_view(), None, true)
            .unwrap();
        let expected = input("(x^2147483647+x^2147483645*z+y*x^2+y*z)*g(p(mink(4)),q(mink(4)))");
        assert_eq!(result.root.expand(), expected.expand());
    }

    #[test]
    fn terminal_vector_metadata_is_certified_once_per_variable() {
        let contractor = setup();
        let reps = test_support::test_initialize();
        let slot = reps.mink4.to_symbolic([Atom::num(991127)]);
        let known = parse_lit!(terminal_vector_known(x));
        let unknown = parse_lit!(terminal_vector_unknown(x));
        let head = vector_symbol!("terminal_variable_metadata");
        let known_vector = function!(head, &known, &slot);
        let unknown_vector = function!(head, &unknown, &slot);
        let mut matcher = SlotMatcher::default();
        let mut sum = ComponentSum::new(&contractor, Intake::Contraction, &mut matcher);
        let AtomView::Fun(known_fun) = known_vector.as_view() else {
            unreachable!()
        };
        let AtomView::Fun(unknown_fun) = unknown_vector.as_view() else {
            unreachable!()
        };
        let first = sum.vector(known_fun).unwrap();
        let second = sum.vector(unknown_fun).unwrap();
        sum.variable(Variable::Vector(first, slot.as_view()), 1);
        sum.variable(Variable::Vector(second, slot.as_view()), 1);
        // An admitted scalar metadata certificate can retain pending notation
        // work. Plain inert metadata alone cannot certify identity freedom.
        sum.scalar_forms.insert(known.as_view(), true);
        TERMINAL_VARIABLE_OBSERVATIONS.with(|count| count.set(0));
        for _ in 0..20 {
            let certificate = sum.terminal_variable(0).unwrap();
            assert_eq!(certificate.slots.len(), 1);
            assert!(certificate.needs_dots);
            assert!(sum.terminal_variable(1).is_none());
        }
        assert_eq!(TERMINAL_VARIABLE_OBSERVATIONS.with(Cell::get), 2);
    }

    #[test]
    fn frontier_sum_merges_shifted_sources_over_prefix_maps() {
        let map: Arc<Vec<PolyVariable>> = Arc::new((0..3).map(PolyVariable::Temporary).collect());
        let polynomial = |nvars: usize, terms: &[(i64, &[u32])]| {
            ScalarPolynomial::from_coefficient_list(
                terms
                    .iter()
                    .map(|(coefficient, _)| Rational::from(*coefficient))
                    .collect(),
                terms
                    .iter()
                    .flat_map(|(_, exponents)| exponents.iter().copied())
                    .collect(),
                Arc::new(map[..nvars].to_vec()),
                &Q,
            )
        };
        let monomial = |coefficient: i64, shift: &[u32]| {
            ScalarPolynomial::new(&Q, Some(1), map.clone())
                .monomial(Rational::from(coefficient), shift.to_vec())
        };
        // Sources over prefixes of the shared map, with overlapping shifts,
        // a non-unit scale and an exact cancellation between two streams.
        let first = polynomial(2, &[(1, &[0, 0]), (2, &[1, 0]), (3, &[0, 2])]);
        let second = polynomial(3, &[(1, &[1, 0, 0]), (-1, &[0, 1, 1])]);
        let third = polynomial(1, &[(4, &[2])]);
        let contributions = [
            (0, 1, [0, 0, 1]),
            (0, -2, [1, 1, 0]),
            (1, 3, [0, 0, 0]),
            (2, 1, [0, 0, 0]),
            // -4*x0^2 cancels the term 4*x0^2 of the previous stream.
            (0, -4, [2, 0, 0]),
        ];
        let expected = contributions
            .iter()
            .map(|&(source, coefficient, shift)| {
                [&first, &second, &third][source] * &monomial(coefficient, &shift)
            })
            .fold(ScalarPolynomial::new(&Q, None, map.clone()), |sum, term| {
                &sum + &term
            });
        let mut sources = vec![Some(first), Some(second), Some(third)];
        let mut uses = vec![3, 1, 1];
        let records = contributions
            .iter()
            .map(|&(source, coefficient, shift)| Contribution {
                source,
                coefficient: Rational::from(coefficient),
                shift: shift.to_vec(),
            })
            .collect::<Vec<_>>();
        let sum = CoefficientSum::frontier_sum(&mut sources, &mut uses, &records, &map).unwrap();
        assert_eq!(sum, expected);
        assert_eq!(sum.get_vars_ref(), map.as_ref());
        // Every source is used up and released.
        assert!(sources.iter().all(Option::is_none));

        // A sole unscaled, unshifted use over the full map moves the source.
        let full = polynomial(3, &[(5, &[0, 1, 2])]);
        let mut sources = vec![Some(full.clone())];
        let identity = [Contribution {
            source: 0,
            coefficient: Rational::one(),
            shift: vec![0; 3],
        }];
        assert_eq!(
            CoefficientSum::frontier_sum(&mut sources, &mut [1], &identity, &map).unwrap(),
            full
        );
        assert!(sources[0].is_none());

        // A source whose variables are not a prefix of the map is declined.
        let foreign =
            ScalarPolynomial::new(&Q, Some(1), Arc::new(vec![PolyVariable::Temporary(7)]))
                .monomial(Rational::one(), vec![1]);
        let mut sources = vec![Some(foreign)];
        assert!(CoefficientSum::frontier_sum(&mut sources, &mut [2], &identity, &map).is_none());
    }

    #[test]
    fn frontier_coefficient_sum_cancels_across_sizes_and_variable_orders() {
        let polynomial = |variables: &[usize], terms: &[(i64, [u32; 2])]| {
            ScalarPolynomial::from_coefficient_list(
                terms
                    .iter()
                    .map(|(coefficient, _)| Rational::from(*coefficient))
                    .collect(),
                terms.iter().flat_map(|(_, exponents)| *exponents).collect(),
                Arc::new(
                    variables
                        .iter()
                        .copied()
                        .map(PolyVariable::Temporary)
                        .collect(),
                ),
                &Q,
            )
        };
        let first = polynomial(
            &[0, 1],
            &[
                (1, [1, 0]),
                (2, [0, 1]),
                (3, [1, 1]),
                (-3, [0, 0]),
                (7, [0, u32::MAX]),
            ],
        );
        let opposite = polynomial(
            &[1, 0],
            &[(-2, [1, 0]), (-1, [0, 1]), (-3, [1, 1]), (3, [0, 0])],
        );
        // Addition must preserve machine-width exponents without multiplying
        // them, and cancellation must work across the two map orderings.
        let large_power = polynomial(&[1, 0], &[(-7, [u32::MAX, 0])]);
        let mut sum = CoefficientSum::One(first);
        // The five- and four-term inputs cancel down into the occupied
        // one-term bin, where the remaining large power must cancel as well.
        for coefficient in [large_power, opposite] {
            sum.push(
                coefficient,
                #[cfg(feature = "reference-cases")]
                0,
            );
        }
        assert!(
            sum.finish(
                #[cfg(feature = "reference-cases")]
                0
            )
            .is_zero()
        );

        let retained = polynomial(&[3, 2], &[(5, [0, 1])]);
        let mut sum = CoefficientSum::One(retained.clone());
        for coefficient in [
            -retained.clone(),
            ScalarPolynomial::new_zero(&Q),
            retained.clone(),
        ] {
            sum.push(
                coefficient,
                #[cfg(feature = "reference-cases")]
                0,
            );
        }
        assert_eq!(
            sum.finish(
                #[cfg(feature = "reference-cases")]
                0
            ),
            retained
        );

        // Finalizing small bins first must not impose a small partial's map
        // order on a large partial. Swapping these variables reverses the
        // large polynomial's monomial order, making insertion quadratic.
        let terms = (0..=128u32)
            .map(|power| (1, [power, 128 - power]))
            .collect::<Vec<_>>();
        let large = polynomial(&[0, 1], &terms);
        let small = polynomial(&[1, 0], &[(1, [129, 0])]);
        let large_variables = large.variables().clone();
        let expected = &large + &small;
        let mut sum = CoefficientSum::One(large);
        sum.push(
            small,
            #[cfg(feature = "reference-cases")]
            0,
        );
        let result = sum.finish(
            #[cfg(feature = "reference-cases")]
            0,
        );
        assert_eq!(result.variables(), &large_variables);
        assert_eq!(result, expected);

        #[cfg(feature = "reference-cases")]
        {
            let surrounding_bytes = 1 << 20;
            let ((result, retained_bytes), clocks) =
                timing::measure("coefficient workspace", || {
                    timing::frontier(Frontier {
                        selected: 0,
                        states: 1,
                        local_terms: 4,
                        generated: 0,
                        predicted_bytes: 0,
                        coefficient_bytes: 0,
                        max_coefficient_bytes: 0,
                    });
                    let mut sum = CoefficientSum::One(polynomial(&[0, 1], &[(1, [0, 0])]));
                    let mut retained_bytes = sum.storage_bytes();
                    for coefficient in [
                        polynomial(&[1, 0], &[(1, [1, 0]), (1, [2, 0])]),
                        polynomial(
                            &[2, 0],
                            &[(1, [1, 0]), (1, [2, 0]), (1, [3, 0]), (1, [4, 0])],
                        ),
                        polynomial(&[0, 1], &[(1, [1, 0])]),
                    ] {
                        sum.push(coefficient, surrounding_bytes);
                        retained_bytes = retained_bytes.max(sum.storage_bytes());
                    }
                    (sum.finish(surrounding_bytes), retained_bytes)
                });
            assert_eq!(result.nterms(), 8);
            let report = clocks.to_string();
            let predicted: usize = report
                .split_once("\"predicted_bytes\":")
                .unwrap()
                .1
                .split(',')
                .next()
                .unwrap()
                .parse()
                .unwrap();
            assert!(
                predicted >= surrounding_bytes + retained_bytes,
                "workspace includes all retained partial buffers: {report}"
            );
        }
    }

    #[test]
    fn frontier_coefficient_sum_bounds_partials_and_matches_bulk_collection() {
        let variables = Arc::new(vec![PolyVariable::Temporary(0)]);
        let mut coefficients = Vec::new();
        let mut exponents = Vec::new();
        let mut sum = CoefficientSum::One(ScalarPolynomial::new_zero(&Q));
        // A skewed arrival order includes both very small and large exponents;
        // no contiguous degree box or input-order assumption is permitted.
        for index in 0..1025u32 {
            let exponent = if index % 2 == 0 {
                index / 2
            } else {
                u32::MAX - index / 2
            };
            let coefficient = Rational::from(if index % 3 == 0 { -1 } else { 1 });
            coefficients.push(coefficient.clone());
            exponents.push(exponent);
            let term = ScalarPolynomial::new(&Q, Some(1), variables.clone())
                .monomial(coefficient, vec![exponent]);
            sum.push(
                term,
                #[cfg(feature = "reference-cases")]
                0,
            );
            if let CoefficientSum::Binned(bins) = &sum {
                assert!(bins.len() <= usize::BITS as usize);
                for (bin, partial) in bins.iter().enumerate() {
                    if let Some(partial) = partial {
                        assert_eq!(partial.nterms().ilog2() as usize, bin);
                    }
                }
            }
        }
        let expected =
            ScalarPolynomial::from_coefficient_list(coefficients, exponents, variables, &Q);
        assert_eq!(
            sum.finish(
                #[cfg(feature = "reference-cases")]
                0
            ),
            expected
        );
    }

    #[test]
    fn frontier_coefficient_sum_matches_expanded_contraction_with_cancellations() {
        let contractor = setup();
        let source = input(concat!(
            "(x*p(mink(4,a))+y*q(mink(4,a))+z*r(mink(4,a))+w*s(mink(4,a)))",
            "*(x*p(mink(4,b))-y*q(mink(4,b))+z*r(mink(4,b))-w*s(mink(4,b)))",
            "*g(mink(4,a),mink(4,b))",
        ));
        let expected = input(concat!(
            "x^2*g(p(mink(4)),p(mink(4)))+2*x*z*g(p(mink(4)),r(mink(4)))",
            "+z^2*g(r(mink(4)),r(mink(4)))-y^2*g(q(mink(4)),q(mink(4)))",
            "-2*y*w*g(q(mink(4)),s(mink(4)))-w^2*g(s(mink(4)),s(mink(4)))",
        ));
        assert_eq!(source.expand().schoonschip().expand(), expected);
        for order in [None, Some([2, 1, 0].as_slice())] {
            let result = contractor
                .contract_factorized(source.as_view(), order, true)
                .unwrap();
            assert_eq!(result.status, ReductionStatus::Complete);
            assert_eq!(result.root.expand(), expected);
        }
    }

    #[test]
    fn factorized_builtin_tensor_sums_preserve_antisymmetric_signs_and_typed_zero() {
        test_support::test_initialize();
        let [mu, nu, rho, sigma, replacement] =
            [97901, 97902, 97903, 97904, 97905].map(|index| mink!(4, Atom::num(index)));
        let volume = epsilon!(&mu, &nu, &rho, &sigma);
        let spectator = input("(x+y)^11");
        let compact_p = p!(mink!(4));
        let compact_q = q!(mink!(4));
        let cases = [
            (
                (p!(&mu) + q!(&mu)) * p!(&nu) * &volume,
                epsilon!(&compact_q, &compact_p, &rho, &sigma),
            ),
            (
                (p!(&mu) + q!(&mu)) * (p!(&nu) + q!(&nu)) * &volume,
                Atom::Zero,
            ),
            (
                g!(&mu, &replacement) * &volume,
                epsilon!(&replacement, &nu, &rho, &sigma),
            ),
        ];
        for (expression, expected) in cases {
            let source = SymbolicTensor::infer(&spectator * expression).unwrap();
            let contractor = SlotContraction::new().admitted(&source);
            INFERENCE_CALLS.with(|count| count.set(0));
            INFERENCE_VISITS.with(|count| count.set(0));
            let result = contractor
                .contract_factorized(source.expression.as_view(), None, true)
                .expect("intrinsic antisymmetric leaves admit factorized contraction");
            assert_eq!(result.status, ReductionStatus::Complete);
            assert_eq!(result.root, &spectator * &expected);
            let contracted = source.contract(Default::default()).unwrap();
            assert_eq!(contracted.expression(), &result.root);
            assert_eq!(contracted.structure, source.structure);
            if !contracted.is_zero() {
                let network = contracted
                    .parse_to_symbolic_net::<AbstractIndex>(&ParseSettings::default())
                    .unwrap();
                assert_eq!(network.graph.dangling_indices().len(), contracted.rank());
            }
            assert_eq!(
                contracted
                    .contract(Default::default())
                    .unwrap()
                    .expression(),
                &result.root
            );
            assert_eq!(INFERENCE_CALLS.with(|count| count.get()), 0);
            assert_eq!(INFERENCE_VISITS.with(|count| count.get()), 0);
        }
    }

    #[test]
    fn factorized_color_tensor_sums_preserve_dummy_orientation() {
        let reps = test_support::test_initialize();
        let [a, b, c, d, e, renamed] = [97911, 97912, 97913, 97914, 97915, 97916]
            .map(|index| reps.coad_da.to_symbolic([Atom::num(index)]));
        let mu = mink!(4, 97917);
        let nu = mink!(4, 97918);
        let first = color_f!(&a, &b, &e) * color_f!(&c, &d, &e);
        let second = color_f!(&a, &b, &renamed) * color_f!(&c, &d, &renamed);
        // The color factors share a sum with Lorentz sources, as in the
        // quartic vertex. They must not refuse the whole connected frontier.
        let source =
            SymbolicTensor::infer(&first * p!(&mu) * q!(&mu) + &second * p!(&nu) * q!(&nu))
                .unwrap();
        let contractor = SlotContraction::new().admitted(&source);
        INFERENCE_CALLS.with(|count| count.set(0));
        INFERENCE_VISITS.with(|count| count.set(0));
        let result = contractor
            .contract_factorized(source.expression.as_view(), None, true)
            .expect("explicit color tensor ports admit the factorized frontier");
        assert_eq!(result.status, ReductionStatus::Complete);
        let coefficient = Atom::num(2) * g!(p!(mink!(4)), q!(mink!(4)));
        assert!(result.root == &coefficient * first || result.root == coefficient * second);
        let contracted = result.into_tensor(&source).unwrap().root;
        assert_eq!(contracted.structure, source.structure);
        let network = contracted
            .parse_to_symbolic_net::<AbstractIndex>(&ParseSettings::default())
            .unwrap();
        assert_eq!(network.graph.dangling_indices().len(), contracted.rank());
        assert_eq!(INFERENCE_CALLS.with(|count| count.get()), 0);
        assert_eq!(INFERENCE_VISITS.with(|count| count.get()), 0);
    }

    #[test]
    fn source_free_color_scopes_stay_opaque_beside_real_contractions() {
        let reps = test_support::test_initialize();
        let [a, b, c, d, e, f] = [97801, 97802, 97803, 97804, 97805, 97806]
            .map(|index| reps.coad_da.to_symbolic([Atom::num(index)]));
        let rep = reps.cof_nc.to_symbolic([] as [Atom; 0]);
        let symmetric = |x: &Atom| {
            spenso::trace_sym!(
                &rep,
                crate::color_t!(x),
                crate::color_t!(&c),
                crate::color_t!(&d)
            )
        };
        let invariant = crate::color_f!(&a, &b, &e)
            * symmetric(&e)
            * crate::color_f!(&a, &b, &f)
            * symmetric(&f);
        let spectator = Atom::mul_many((0..32).map(|index| {
            let coefficient =
                Atom::var(symbolica::symbol!(&format!("opaque_color_scope::x{index}")));
            (coefficient + &invariant).pow(2)
        }));
        let mu = spenso::mink!(4, 97807);
        let nu = spenso::mink!(4, 97808);
        let source =
            SymbolicTensor::infer(&spectator * spenso::g!(&mu, &nu) * spenso::p!(&nu)).unwrap();
        let contractor = SlotContraction::new().expanding(false).admitted(&source);
        crate::tensor::SELECTED_OPENS.with(|count| count.set(0));
        let result = contractor
            .contract_factorized(source.expression.as_view(), None, true)
            .unwrap();
        assert_eq!(result.root, spectator * spenso::p!(&mu));
        assert_eq!(result.status, ReductionStatus::Complete);
        assert!(
            result.observations.is_none(),
            "remaining colour identities are not terminal scalar algebra"
        );
        assert_eq!(crate::tensor::SELECTED_OPENS.with(|count| count.get()), 0);
    }

    #[test]
    fn source_free_open_color_sum_keeps_positional_metric_bindings() {
        let reps = test_support::test_initialize();
        let [a, b, c, d, e, z] = [97811, 97812, 97813, 97814, 97815, 97816]
            .map(|index| reps.coad_da.to_symbolic([Atom::num(index)]));
        let rep = reps.cof_nc.to_symbolic([] as [Atom; 0]);
        let sum = |port: &Atom| {
            crate::color_f!(port, &b, &e)
                * spenso::trace_sym!(
                    &rep,
                    crate::color_t!(&e),
                    crate::color_t!(&c),
                    crate::color_t!(&d)
                )
                + crate::color_f!(port, &c, &e)
                    * spenso::trace_sym!(
                        &rep,
                        crate::color_t!(&e),
                        crate::color_t!(&b),
                        crate::color_t!(&d)
                    )
        };
        let source = SymbolicTensor::infer(spenso::g!(&a, &z) * sum(&a)).unwrap();
        let contractor = SlotContraction::new().expanding(false).admitted(&source);
        crate::tensor::SELECTED_OPENS.with(|count| count.set(0));
        let result = contractor
            .contract_factorized(source.expression.as_view(), None, true)
            .unwrap();
        assert_eq!(result.root, sum(&z));
        assert_eq!(result.status, ReductionStatus::Complete);
        assert!(result.observations.is_none());
        assert_eq!(crate::tensor::SELECTED_OPENS.with(|count| count.get()), 0);
        assert_eq!(
            SymbolicTensor::infer(result.root).unwrap().structure,
            source.structure
        );
    }

    #[test]
    fn source_free_callback_region_still_requires_checked_contraction() {
        test_support::test_initialize();
        let a = spenso::mink!(4, 97821);
        let b = spenso::mink!(4, 97822);
        let target = b.clone();
        let callback = spenso::tensor_symbol!(
            "source_proof_callback",
            norm = move |value, output| {
                if let AtomView::Fun(function) = value
                    && function.iter().any(|argument| argument == target.as_view())
                {
                    **output = Atom::one();
                }
            }
        );
        let sibling = spenso::tensor_symbol!("source_proof_callback_sibling");
        let source = SymbolicTensor::infer(
            spenso::g!(&a, &b)
                * (symbolica::function!(callback, &a) + symbolica::function!(sibling, &a)),
        )
        .unwrap();
        let contractor = SlotContraction::new().admitted(&source);
        assert!(
            contractor
                .contract_factorized(source.expression.as_view(), None, true)
                .is_none()
        );
        assert!(source.contract(Default::default()).is_err());
    }

    #[test]
    fn negative_power_contracts_its_existing_graph_scope() {
        let contractor = setup();
        for (source, expected) in [
            (
                "spenso::bracket(p(spenso::mink(4,a))*q(spenso::mink(4,a)))^-1",
                "spenso::g(p(spenso::mink(4)),q(spenso::mink(4)))^-1",
            ),
            (
                "x+spenso::bracket(p(spenso::mink(4,a))*q(spenso::mink(4,a)))^-2",
                "x+spenso::g(p(spenso::mink(4)),q(spenso::mink(4)))^-2",
            ),
        ] {
            let source = input(source);
            let result = contractor
                .contract_factorized(source.as_view(), None, true)
                .unwrap();
            assert_eq!(result.root, input(expected), "{source}");
        }
    }

    #[test]
    fn factorized_states_retain_exact_complex_coefficients() {
        let contractor = setup();
        let source = input("1i*p(mink(4,a))*q(mink(4,a))+2i*r(mink(4,b))*s(mink(4,b))");
        let result = contractor
            .contract_factorized(source.as_view(), None, true)
            .unwrap();
        assert!(result.status == ReductionStatus::Complete);
        assert_eq!(
            result.root,
            input("1i*g(p(mink(4)),q(mink(4)))+2i*g(r(mink(4)),s(mink(4)))")
        );
        let inexact = input("1.25*p(mink(4,a))*q(mink(4,a))+2.5*r(mink(4,b))*s(mink(4,b))");
        assert!(
            contractor
                .contract_factorized(inexact.as_view(), None, true)
                .is_none()
        );
    }

    #[test]
    fn factorized_states_use_the_component_reducer_in_both_orders() {
        let contractor = setup();
        for source in [
            "g(mink(4,a),mink(4,b))*(p(mink(4,a))+q(mink(4,a)))*r(mink(4,b))",
            "(g(mink(4,a),mink(4,b))*p(mink(4,c))+g(mink(4,a),mink(4,c))*q(mink(4,b)))*(p(mink(4,a))*q(mink(4,b))*r(mink(4,c))+q(mink(4,a))*r(mink(4,b))*s(mink(4,c)))",
            "(g(mink(4,a),mink(4,b))+p(mink(4,a))*q(mink(4,b)))*(g(mink(4,a),mink(4,b))+r(mink(4,a))*s(mink(4,b)))",
            "(g(mink(4,a),mink(4,b))+t(mink(4,a),mink(4,b)))*p(mink(4,a))*q(mink(4,b))",
        ] {
            let source = input(source);
            let expected = source.expand().schoonschip().expand();
            let count = match source.as_view() {
                AtomView::Mul(product) => product.iter().len(),
                _ => 1,
            };
            let reversed = (0..count).rev().collect::<Vec<_>>();
            for order in [None, Some(reversed.as_slice())] {
                let contracted = contractor
                    .contract_factorized(source.as_view(), order, true)
                    .unwrap();
                SymbolicTensor::infer(contracted.root.clone()).unwrap();
                assert_eq!(contracted.root.expand(), expected, "{source}");
            }
        }
    }

    #[test]
    fn factorized_states_keep_foreign_and_disconnected_sums() {
        let contractor = setup();
        let spectator = input("(1+x)*(t(mink(4,c))+routing(z)*u(mink(4,c)))");
        let first = input("(p(mink(4,a))+q(mink(4,a)))*r(mink(4,a))");
        let second = input("(p(mink(4,b))+s(mink(4,b)))*q(mink(4,b))");
        let source = Atom::mul_many([&spectator, &first, &second]);
        let contracted = contractor
            .contract_factorized(source.as_view(), None, true)
            .unwrap();
        let result = contracted.root;
        let expected = Atom::mul_many([
            spectator,
            first.expand().schoonschip(),
            second.expand().schoonschip(),
        ]);
        assert_eq!(result, expected);
    }

    #[test]
    fn factorized_port_bindings_remain_occurrence_local() {
        let contractor = setup();
        let source =
            input("(p(mink(4,a))+q(mink(4,a)))*(p(mink(4,b))+q(mink(4,b)))*t(mink(4,a),mink(4,b))");
        let contracted = contractor
            .contract_factorized(source.as_view(), None, true)
            .unwrap();
        assert_eq!(
            contracted.root.expand(),
            source.expand().schoonschip().expand()
        );
        assert!(
            contractor
                .contract_factorized(source.as_view(), Some(&[0, 0, 1]), true)
                .is_none()
        );
    }

    #[test]
    fn polynomial_exponents_beyond_u16_remain_exact_and_complete() {
        let contractor = setup();
        let x = input("x");
        let spectator = input("(y+z)^17");
        let first = &spectator * x.pow(u16::MAX);
        let second = x.clone();
        let mut slots = SlotMatcher::default();
        let mut sum = ComponentSum::new(&contractor, Intake::Contraction, &mut slots);
        let leaf = sum.input.leaf(InputLeaf::Scalar(x.as_view())).unwrap();
        let power = sum.input.power(leaf, u32::from(u16::MAX)).unwrap();
        let scalar = sum
            .input
            .leaf(InputLeaf::Scalar(spectator.as_view()))
            .unwrap();
        let weighted = sum.input.group([power, scalar].into_iter(), false).unwrap();
        sum.opaque_factors = vec![
            (
                first.as_view(),
                Vec::new(),
                PartialStructure::from_logical_slots([]),
            ),
            (
                second.as_view(),
                Vec::new(),
                PartialStructure::from_logical_slots([]),
            ),
        ];
        sum.factor_roots = vec![Some(weighted.0), Some(leaf.0)];
        let result = sum
            .contract_states(
                Remaining {
                    factors: vec![(0, Vec::new()), (1, Vec::new())],
                    residual: Vec::new(),
                },
                &[0, 1],
                &mut |_, _| unreachable!("the two selected factors are already compiled"),
            )
            .unwrap();
        assert_eq!(result.status, ReductionStatus::Complete);
        assert_eq!(result.root, &spectator * x.pow(u32::from(u16::MAX) + 1));
        assert!(result.observations.is_some());
        assert!(sum.overrides.is_empty());
    }

    #[test]
    fn skewed_frontier_counts_each_opaque_coefficient_once() {
        let contractor = setup();
        let scalar = Atom::add_many(
            (0..512).map(|i| Atom::var(symbolica::symbol!(&format!("skewed_frontier::s{i}")))),
        );
        let source = (&scalar * input("p(mink(4,a))") + input("q(mink(4,a))"))
            * input(
                "(p(mink(4,b))+q(mink(4,b)))*(p(mink(4,c))+q(mink(4,c)))*t(mink(4,a),mink(4,b),mink(4,c))",
            );
        let result = contractor
            .contract_factorized(source.as_view(), None, true)
            .unwrap();
        assert_eq!(result.status, ReductionStatus::Complete);
        let mut expected = Vec::new();
        for a in ["p", "q"] {
            for b in ["p", "q"] {
                for c in ["p", "q"] {
                    let tensor = input(&format!("t({a}(mink(4)),{b}(mink(4)),{c}(mink(4)))"));
                    expected.push(if a == "p" { &scalar * tensor } else { tensor });
                }
            }
        }
        assert_eq!(result.root, Atom::add_many(expected));
    }

    #[test]
    fn factorized_frontier_completes_all_selected_contractions() {
        let contractor = setup();
        let slots = (0..9).map(|i| format!("mink(4,a{i})")).collect::<Vec<_>>();
        let source = input(&format!(
            "t({})*{}",
            slots.join(","),
            slots
                .iter()
                .map(|slot| format!("(x*p({slot})+y*q({slot}))"))
                .collect::<Vec<_>>()
                .join("*")
        ));
        let contracted = contractor
            .contract_factorized(source.as_view(), None, true)
            .unwrap();
        assert!(!contracted.root.is_zero());
        assert_eq!(contracted.status, ReductionStatus::Complete);
        assert!(
            !contracted.root.has_repeated_explicit_indices(),
            "all selected connections must be contracted"
        );
        let result = contracted.root;
        assert_eq!(
            result.expand().schoonschip().expand(),
            source.expand().schoonschip().expand()
        );
    }
    #[test]
    fn scalar_power_without_graph_incidence_keeps_its_original_source() {
        test_support::test_initialize();
        let x = Atom::var(symbolica::symbol!("no_work_power_x"));
        let y = Atom::var(symbolica::symbol!("no_work_power_y"));
        let z = Atom::var(symbolica::symbol!("no_work_power_z"));
        let source = spenso::bracket!(x + y, z).pow(3);
        let result = SlotContraction::new()
            .contract_factorized(source.as_view(), None, true)
            .unwrap();
        assert!(result.status == ReductionStatus::Complete);
        assert_eq!(result.root, source);
    }

    #[test]
    fn scalar_rational_sums_keep_their_factorization() {
        let contractor = setup();
        for scalar in ["z", "q(0,spenso::cind(0))", "g(p(mink(4)),q(mink(4)))"] {
            let source = input(&format!(
                "n*(1/({scalar}+x)+1/({scalar}+y))+m*(1/({scalar}-x)+1/({scalar}-y))"
            ));
            for rank_one in [false, true] {
                let result = contractor
                    .contract_factorized(source.as_view(), None, rank_one)
                    .unwrap();
                assert!(result.status == ReductionStatus::Complete);
                // An exact Atom comparison protects the factored scalar input;
                // expanding both sides would hide this regression.
                assert_eq!(result.root, source, "{scalar}, rank_one={rank_one}");

                let tensor = SymbolicTensor::infer(source.clone()).unwrap();
                let settings = ContractSettings::default();
                let settings = if rank_one {
                    settings
                } else {
                    settings.without_rank_one_tensors()
                };
                for settings in [settings, settings.with_order(&[0])] {
                    let result = tensor.contract(settings).unwrap();
                    assert_eq!(result, tensor);
                    assert_eq!(result.expression, source);
                }
            }
        }
    }
}

#[cfg(test)]
mod lazy_tests {
    use super::super::tests as component_tests;
    use super::*;
    use crate::test_support;
    use spenso::network::tags::SPENSO_TAG;
    use std::ptr;
    use std::sync::atomic::{AtomicUsize, Ordering};
    use symbolica::domains::float::Float;
    use symbolica::parser::ParseSettings as SymbolicaParseSettings;
    use symbolica::symbol;

    #[test]
    fn partial_arithmetic_scopes_have_only_one_level() {
        let _contractor = component_tests::setup();
        for (source, nodes) in [
            ("p(spenso::mink(4,a))*q(spenso::mink(4,a))+x*y", 3),
            ("(p(spenso::mink(4,a))*q(spenso::mink(4,a))+x)^2", 2),
            (
                "spenso::bracket(p(spenso::mink(4,a))*q(spenso::mink(4,a)))^-1",
                2,
            ),
        ] {
            let source = component_tests::input(source);
            let scope = ContractionScope::parse(source.as_view()).unwrap();
            assert_eq!(scope.ports.len(), nodes, "{source}");
            assert_eq!(scope.network.graph.graph.n_nodes(), nodes, "{source}");
            assert!(scope.opened.values().all(|entry| entry.get().is_none()));
        }
    }

    #[test]
    fn planned_scope_borrows_graph_and_reopens_composite_leaves() {
        let setup = component_tests::setup();
        let product = component_tests::input("p(mink(4,a))*q(mink(4,a))");
        let source = &product + component_tests::input("x*y");
        let tensor = SymbolicTensor::infer(source.clone()).unwrap();
        let planned = tensor.shallow_graph().unwrap();
        let contractor = setup.admitted(&tensor).planned(
            source.as_view(),
            Arc::clone(&planned),
            Some(tensor.proofs.parser_leaf_interfaces()),
        );
        assert!(contractor.graph.as_ref().unwrap().fallback_root.is_none());
        let scope = contractor.scope(source.as_view()).unwrap();
        assert!(ptr::eq(scope.network.graph.as_ref(), &planned.graph));
        let child = scope
            .tree
            .iter_children(scope.root, &planned.graph.graph)
            .find(|&node| scope.literal(node) == Some(product.as_view()))
            .unwrap();
        let opened = scope.open(child).unwrap();
        assert!(Arc::ptr_eq(
            opened.network.admission.as_ref().unwrap(),
            planned.admission.as_ref().unwrap(),
        ));
        assert!(matches!(
            opened.network.graph.graph[opened.root],
            NetworkNode::Op(NetworkOp::Product)
        ));
    }

    #[test]
    fn selected_planned_scope_matches_its_payload_instead_of_its_enclosing_admission() {
        let setup = component_tests::setup();
        let selected = component_tests::input("p(mink(4,a))*q(mink(4,a))");
        let enclosing = SymbolicTensor::infer(&selected + component_tests::input("x*y")).unwrap();
        let mut local = SymbolicTensor::infer(selected.clone()).unwrap();
        local.replace_admission(enclosing.proofs.clone());
        let graph = local.shallow_graph().unwrap();
        let contractor =
            setup
                .admitted(&local)
                .planned(selected.as_view(), Arc::clone(&graph), None);
        assert!(
            contractor
                .graph
                .as_ref()
                .unwrap()
                .matches(selected.as_view())
        );
        assert!(
            !contractor
                .graph
                .as_ref()
                .unwrap()
                .matches(enclosing.expression.as_view())
        );
        assert!(ptr::eq(
            contractor
                .scope(selected.as_view())
                .unwrap()
                .network
                .graph
                .as_ref(),
            &graph.graph
        ));
        assert!(!ptr::eq(
            contractor
                .scope(enclosing.expression.as_view())
                .unwrap()
                .network
                .graph
                .as_ref(),
            &graph.graph
        ));
    }

    #[test]
    fn selected_arithmetic_reuses_admission_without_nested_graphs() {
        let setup = component_tests::setup();
        let input = component_tests::input;
        let scalar =
            Atom::add_many((0..512).map(|i| Atom::var(symbol!(&format!("selected_tape::s{i}")))));
        let source = (&scalar * input("p(mink(4,a))") + input("q(mink(4,a))"))
            * input("t(mink(4,a),mink(4,b))");
        let tensor = SymbolicTensor::infer(source.clone()).unwrap();
        let contractor = setup.admitted(&tensor);
        let scope = ContractionScope::parse(source.as_view()).unwrap();
        let result = contractor
            .contract_graph(&scope, scope.root, None, true)
            .unwrap();
        assert_eq!(result.status, ReductionStatus::Complete);
        assert_eq!(
            result.root,
            &scalar * input("t(p(mink(4)),mink(4,b))") + input("t(q(mink(4)),mink(4,b))")
        );
        assert!(
            scope.opened.values().all(|entry| entry.get().is_none()),
            "admitted arithmetic needs a tape traversal, not new port graphs"
        );
    }

    #[test]
    fn direct_selected_arithmetic_keeps_independent_dummy_scopes() {
        let setup = component_tests::setup();
        let input = component_tests::input;
        let source = input(
            "(p(mink(4,a))*q(mink(4,a))*r(mink(4,b))+s(mink(4,b)))*(p(mink(4,a))*q(mink(4,a))*t(mink(4,b),mink(4,c))+u(mink(4,b),mink(4,c)))",
        );
        let tensor = SymbolicTensor::infer(source.clone()).unwrap();
        let contractor = setup.admitted(&tensor);
        let scope = ContractionScope::parse(source.as_view()).unwrap();
        let result = contractor
            .contract_graph(&scope, scope.root, None, true)
            .unwrap();
        let dot = input("g(p(mink(4)),q(mink(4)))");
        let expected = dot.pow(2) * input("t(r(mink(4)),mink(4,c))")
            + &dot * input("t(s(mink(4)),mink(4,c))")
            + &dot * input("u(r(mink(4)),mink(4,c))")
            + input("u(s(mink(4)),mink(4,c))");
        assert_eq!(result.status, ReductionStatus::Complete);
        assert_eq!(result.root.expand(), expected.expand());
        assert!(scope.opened.values().all(|entry| entry.get().is_none()));
    }

    #[test]
    fn unobserved_powers_require_intrinsic_scalar_certificates() {
        let contractor = component_tests::setup();
        let opaque = symbol!("unobserved_power::opaque"; Scalar);
        for value in [
            component_tests::input("m^2"),
            component_tests::input("(m+n)^-2"),
            component_tests::input("(m+n)^(1/3)"),
            opaque.call(component_tests::input("(m+n)^11")).pow(2),
        ] {
            let mut slots = SlotMatcher::default();
            let mut sum = ComponentSum::new(&contractor, Intake::Contraction, &mut slots);
            let mut tape = TermTape::default();
            contractor
                .compile_arithmetic(value.as_view(), &mut tape, &mut sum)
                .unwrap();
            assert_eq!(tape.leaves().len(), 1);
            assert!(
                matches!(tape.leaves()[0], InputLeaf::Scalar(original) if original == value.as_view())
            );
        }
        let calls = Arc::new(AtomicUsize::new(0));
        let observed = Arc::clone(&calls);
        let callback = symbol!("unobserved_power::callback"; Scalar; norm=move|_,_| {
            observed.fetch_add(1, Ordering::Relaxed);
        });
        let rounded = Atom::num(Float::parse("0.1", Some(11)).unwrap());
        for value in [
            component_tests::input("p(mink(4,a))^2"),
            component_tests::input("(p(mink(4,a))*q(mink(4,a))+m)^2"),
            callback.call(component_tests::input("m")).pow(2),
            (rounded + component_tests::input("m")).pow(2),
        ] {
            calls.store(0, Ordering::Relaxed);
            let mut slots = SlotMatcher::default();
            let mut sum = ComponentSum::new(&contractor, Intake::Contraction, &mut slots);
            let mut tape = TermTape::default();
            assert!(
                contractor
                    .compile_arithmetic(value.as_view(), &mut tape, &mut sum)
                    .is_none()
            );
            assert_eq!(calls.load(Ordering::Relaxed), 0);
        }
    }

    #[test]
    fn scalar_powers_in_selected_sums_do_not_rebuild_graphs() {
        let setup = component_tests::setup();
        let input = component_tests::input;
        let source = input(
            "(m^4*p(mink(4,a))+m^2*g(q(mink(4)),s(mink(4)))^2*r(mink(4,a)))*t(mink(4,a),mink(4,b))",
        );
        let tensor = SymbolicTensor::infer(source.clone()).unwrap();
        let contractor = setup.admitted(&tensor);
        let scope = ContractionScope::parse(source.as_view()).unwrap();
        let result = contractor
            .contract_graph(&scope, scope.root, None, true)
            .unwrap();
        assert_eq!(result.status, ReductionStatus::Complete);
        assert_eq!(
            result.root,
            input(
                "m^4*t(p(mink(4)),mink(4,b))+m^2*g(q(mink(4)),s(mink(4)))^2*t(r(mink(4)),mink(4,b))",
            )
        );
        assert!(scope.opened.values().all(|entry| entry.get().is_none()));
    }

    #[test]
    fn compact_tensor_sums_do_not_enter_the_scalar_coefficient_tape() {
        let setup = component_tests::setup();
        let input = component_tests::input;
        let source = input("p(mink(4))+q(mink(4))");
        let tensor = SymbolicTensor::infer(source.clone()).unwrap();
        assert_eq!(tensor.structure.logical_slots().len(), 1);
        let contractor = setup.admitted(&tensor);
        let mut slots = SlotMatcher::default();
        let mut sum = ComponentSum::new(&contractor, Intake::Contraction, &mut slots);
        let mut tape = TermTape::default();
        contractor
            .compile_arithmetic(source.as_view(), &mut tape, &mut sum)
            .unwrap();
        // The rank-one alternatives keep their tensor ports. Zero written
        // indices does not establish the scalar interface required by a
        // coefficient, even before unresolved ports have been materialized.
        assert!(sum.scalar_interfaces.is_empty());
        assert_eq!(tape.leaves().len(), 2);
        assert!(
            tape.leaves()
                .iter()
                .all(|leaf| matches!(leaf, InputLeaf::Atom(_)))
        );
    }

    #[test]
    fn indexed_powers_in_selected_sums_keep_the_scoped_graph() {
        let setup = component_tests::setup();
        let input = component_tests::input;
        let source = input(
            "((p(mink(4,a))*q(mink(4,a))+x)^2*r(mink(4,b))+s(mink(4,b)))*t(mink(4,b),mink(4,c))",
        );
        let tensor = SymbolicTensor::infer(source.clone()).unwrap();
        let contractor = setup.admitted(&tensor);
        let scope = ContractionScope::parse(source.as_view()).unwrap();
        let result = contractor
            .contract_graph(&scope, scope.root, None, true)
            .unwrap();
        assert_eq!(result.status, ReductionStatus::Complete);
        assert_eq!(
            result.root.expand(),
            input(
                "(g(p(mink(4)),q(mink(4)))+x)^2*t(r(mink(4)),mink(4,c))+t(s(mink(4)),mink(4,c))",
            )
            .expand()
        );
        assert!(scope.opened.values().any(|entry| entry.get().is_some()));
    }

    #[test]
    fn fixed_component_scalar_sums_stay_closed() {
        let contractor = component_tests::setup();
        let source =
            component_tests::input("n*(1/(q(0,spenso::cind(0))+x)+1/(q(0,spenso::cind(0))+y))");
        let scope = ContractionScope::parse(source.as_view()).unwrap();
        let result = contractor
            .contract_graph(&scope, scope.root, None, true)
            .unwrap();
        assert!(result.status == ReductionStatus::Complete);
        assert_eq!(result.root, source);
        assert!(scope.opened.values().all(|entry| entry.get().is_none()));
    }

    #[test]
    fn graph_order_starts_with_a_small_boundary_and_accepts_an_override() {
        let contractor = component_tests::setup();
        let source = component_tests::input(
            "spenso::g(spenso::mink(4,a),spenso::mink(4,b))*p(spenso::mink(4,a))*q(spenso::mink(4,b))",
        );
        let scope = ContractionScope::parse(source.as_view()).unwrap();
        let nodes = scope.factors(scope.root);
        let selected = (0..nodes.len()).collect::<Vec<_>>();
        let order = contractor
            .factor_order(&scope, &nodes, selected.clone(), None)
            .unwrap();
        assert_eq!(scope.ports[&nodes[order[0]]].0.len(), 1);
        let explicit = selected.iter().rev().copied().collect::<Vec<_>>();
        assert_eq!(
            contractor
                .factor_order(&scope, &nodes, selected, Some(&explicit))
                .unwrap(),
            explicit
        );
        let expected = component_tests::input("spenso::g(p(spenso::mink(4)),q(spenso::mink(4)))");
        for order in [None, Some(explicit.as_slice())] {
            let result = contractor
                .contract_graph(&scope, scope.root, order, true)
                .unwrap();
            assert!(result.status == ReductionStatus::Complete);
            assert_eq!(result.root, expected);
        }
    }

    #[test]
    fn unselected_thousand_term_factors_stay_closed() {
        test_support::test_initialize();
        SPENSO_TAG.rank_one_tensor_symbol("lazy_contraction::p");
        SPENSO_TAG.tensor_symbol("lazy_contraction::t");
        SPENSO_TAG.tensor_symbol("lazy_contraction::u");
        let parse = |source: &str| {
            Atom::parse(
                source,
                "lazy_contraction",
                SymbolicaParseSettings::symbolica(),
            )
            .unwrap()
        };
        for foreign in [false, true] {
            let first = Atom::add_many((0..1000).map(|i| {
                parse(&if foreign {
                    format!("t({i},spenso::mink(4,c))")
                } else {
                    format!("x({i})")
                })
            }));
            let second = Atom::add_many((0..1000).map(|i| {
                parse(&if foreign {
                    format!("u({i},spenso::mink(4,d))")
                } else {
                    format!("y({i})")
                })
            }));
            let core = parse("spenso::g(spenso::mink(4,a),spenso::mink(4,b))*p(spenso::mink(4,a))");
            let source = Atom::mul_many([core, first.clone(), second.clone()]);
            let scope = ContractionScope::parse(source.as_view()).unwrap();
            assert_eq!(scope.ports.len(), 5, "one product and four factor leaves");
            assert_eq!(scope.network.graph.graph.n_nodes(), 5);
            let protected = scope
                .factors(scope.root)
                .into_iter()
                .filter(|&node| {
                    let value = scope.literal(node).unwrap();
                    value == first.as_view() || value == second.as_view()
                })
                .collect::<Vec<_>>();
            assert_eq!(protected.len(), 2);
            let result = SlotContraction::new()
                .contract_graph(&scope, scope.root, None, true)
                .unwrap();
            assert!(result.status == ReductionStatus::Complete);
            assert_eq!(
                result.root,
                Atom::mul_many([parse("p(spenso::mink(4,b))"), first, second])
            );
            for node in protected {
                assert!(
                    scope.opened[&node].get().is_none(),
                    "unselected factor was opened"
                );
            }
            assert!(scope.opened.values().all(|opened| opened.get().is_none()));
        }
    }
}
