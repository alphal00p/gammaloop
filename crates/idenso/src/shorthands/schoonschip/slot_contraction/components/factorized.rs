//! State merging on the collector's existing port graph.
//!
//! Remaining factors are opaque terminals of the same component reducer used
//! for a completed monomial. Only the selected factor contributes alternatives;
//! substitutions remain port bindings until a distinct output variable is emitted.

use super::*;
use crate::tensor::{
    ContractSettings, SymbolicNet, simplification::observation::DomainObservations,
};
use ahash::AHashSet;
use linnet::{
    half_edge::{
        NodeIndex,
        subgraph::{Inclusion, SuBitGraph},
        tree::SimpleTraversalTree,
    },
    tree::child_vec::ChildVecStore,
};
use spenso::{
    network::{
        graph::{NetworkEdge, NetworkLeaf, NetworkNode, NetworkOp},
        library::DummyLibrary,
        parsing::{LeafInterfaceCache, ParseSettings, ShorthandParsing},
        store::{TensorScalarStore, TensorScalarStoreMapping},
    },
    structure::{
        OrderedStructure, TensorStructure,
        abstract_index::AbstractIndex,
        partial::{PartialIndex, PartialStructureExt},
        slot::IsAbstractSlot,
    },
};
use std::sync::Arc;
use symbolica::{atom::AtomOrView, coefficient::CoefficientView};

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
        std::mem::size_of::<Self>()
            + self.factors.capacity() * std::mem::size_of::<(usize, Vec<Argument<'_>>)>()
            + self
                .factors
                .iter()
                .map(|(_, arguments)| arguments.capacity() * std::mem::size_of::<Argument<'_>>())
                .sum::<usize>()
            + self.residual.capacity() * std::mem::size_of::<(usize, u32)>()
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
            let old = std::mem::replace(old, ScalarPolynomial::new_zero(&Q));
            if first == second {
                *self = Self::One(Self::merge(
                    old,
                    coefficient,
                    #[cfg(feature = "reference-cases")]
                    surrounding_bytes.saturating_add(std::mem::size_of::<Self>()),
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
            unreachable!()
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
                        - std::mem::size_of::<ScalarPolynomial>();
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
        let _merge_timing = crate::reference_cases::timing::scope(
            crate::reference_cases::timing::Phase::CoefficientMerge,
        );
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
                    std::mem::size_of::<Rational>()
                        .saturating_add(variables.saturating_mul(std::mem::size_of::<u32>())),
                );
            let maps = variables.saturating_mul(std::mem::size_of::<PolyVariable>());
            let retained = ComponentSum::coefficient_bytes(&first)
                .saturating_add(ComponentSum::coefficient_bytes(&second));
            crate::reference_cases::timing::frontier_workspace(
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
        std::mem::size_of::<Self>()
            .saturating_add(
                bins.capacity()
                    .saturating_mul(std::mem::size_of::<Option<ScalarPolynomial>>()),
            )
            .saturating_add(
                bins.iter()
                    .flatten()
                    .map(|coefficient| {
                        ComponentSum::coefficient_bytes(coefficient)
                            - std::mem::size_of::<ScalarPolynomial>()
                    })
                    .fold(0usize, usize::saturating_add),
            )
    }

    #[cfg(feature = "reference-cases")]
    // Unused since 75d852458 (frontier merge); kept for the reference-cases accounting.
    #[allow(dead_code)]
    fn storage_bytes(&self) -> usize {
        match self {
            Self::One(coefficient) => {
                std::mem::size_of::<Self>()
                    .saturating_add(ComponentSum::coefficient_bytes(coefficient))
                    - std::mem::size_of::<ScalarPolynomial>()
            }
            Self::Binned(bins) => Self::binned_bytes(bins),
        }
    }
}

/// Facts that certify scalar syntax: everything a contraction reports except
/// a terminal polynomial with open ports, which certifies only its own root.
fn scalar_facts(observed: &DomainObservations) -> bool {
    observed.is_closed_scalar() || !observed.is_terminal_polynomial()
}

/// Sources with fewer terms are multiplied out and summed at once; larger
/// ones are merged lazily (see [`frontier_sum`]).
const LAZY_SOURCE_TERMS: usize = 32;

/// `IDENSO_FRONTIER_TRACE=1`: print one line per frontier step and per emission
/// (states, terms, coefficient storage, time, process RSS) to stderr.
static FRONTIER_TRACE: std::sync::LazyLock<bool> =
    std::sync::LazyLock::new(|| std::env::var_os("IDENSO_FRONTIER_TRACE").is_some());

fn trace_rss_mib() -> u64 {
    std::fs::read_to_string("/proc/self/statm")
        .ok()
        .and_then(|s| {
            s.split_whitespace()
                .nth(1)
                .and_then(|v| v.parse::<u64>().ok())
        })
        .map_or(0, |pages| pages * 4096 / (1 << 20))
}

/// Exponent and inline coefficient storage of one frontier coefficient.
fn trace_poly_bytes(polynomial: &ScalarPolynomial) -> usize {
    polynomial.exponents.len() * std::mem::size_of::<u32>()
        + polynomial.coefficients.len() * std::mem::size_of::<Rational>()
}

/// One contribution to a frontier target: a source state times one monomial,
/// `coefficient` times the variables of the shared map raised to `shift`.
struct Contribution {
    source: usize,
    coefficient: Rational,
    shift: Vec<u32>,
}

/// A cooperative resource-budget checkpoint of the frontier (see the per-term
/// checkpoint in `contract_states`). The frontier cannot decline once
/// started, so an exhausted budget aborts the operation.
fn budget_checkpoint() {
    if let Err(error) = symbolica::utils::resources::checkpoint() {
        panic!("resource budget exhausted in a factorized contraction: {error}");
    }
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
    use symbolica::domains::{Ring, RingOps};
    let nvars = map.len();
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
        let mut sum = ScalarPolynomial::new(&Q, Some(capacity), map.clone());
        let mut cursors = vec![0usize; streams.len()];
        let mut current = vec![0u32; streams.len() * nvars];
        let load = |stream: usize, cursor: usize, current: &mut [u32]| {
            let (source, contribution) = streams[stream];
            let width = source.nvars();
            let exponents = &source.exponents[cursor * width..(cursor + 1) * width];
            let target = &mut current[stream * nvars..(stream + 1) * nvars];
            for (variable, entry) in target.iter_mut().enumerate() {
                *entry = exponents.get(variable).copied().unwrap_or(0)
                    + contribution.shift.get(variable).copied().unwrap_or(0);
            }
        };
        let mut heap = Vec::with_capacity(streams.len());
        for (stream, (source, _)) in streams.iter().enumerate() {
            if !source.is_zero() {
                load(stream, 0, &mut current);
                heap.push(stream);
            }
        }
        for index in (0..heap.len() / 2).rev() {
            sift_down(&mut heap, index, &current, nvars);
        }
        let mut appended = 0usize;
        while let Some(&top) = heap.first() {
            // The merge writes the next frontier: budgets must see it grow.
            appended = appended.wrapping_add(1);
            if appended.is_multiple_of(4096) {
                budget_checkpoint();
            }
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
                load(top, cursors[top], &mut current);
            } else {
                heap.swap_remove(0);
            }
            sift_down(&mut heap, 0, &current, nvars);
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

/// The scalar variables of one frontier, in one append-only polynomial
/// variable map: products and merges of frontier coefficients then never
/// reorder and rebuild a polynomial. When enabled, small composite scalar
/// spectators (such as `2-dim` or a short closed trace) are expanded over their
/// base variables, so equal monomials merge instead of carrying an opaque copy
/// of the sum.
struct SharedScalars {
    atoms: Vec<Atom>,
    positions: AHashMap<Atom, usize>,
    dot_positions: AHashMap<Atom, usize>,
    /// Dot-notation spelling of every position.
    dots: Vec<Atom>,
    map: Arc<Vec<PolyVariable>>,
    expansions: AHashMap<usize, Expansion>,
    expand: bool,
    dot_notation: bool,
}

/// Composite scalars with at most this many polynomial terms are expanded
/// over their base variables; larger ones stay opaque variables.
const EXPANSION_TERM_LIMIT: usize = 64;

/// Composites whose structural term bound (products of sum lengths, powers
/// expanded) exceeds this are not converted at all, so that a large spectator
/// such as `(x+y)^1000` is never expanded only to be discarded.
const EXPANSION_WORK_LIMIT: usize = 1024;

/// Whether `terms^exponent` stays within the expansion limit.
fn small_power(terms: usize, exponent: u32) -> bool {
    terms
        .checked_pow(exponent)
        .is_some_and(|terms| terms <= EXPANSION_TERM_LIMIT)
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

enum Expansion {
    Base(usize),
    Composite(ScalarPolynomial),
}

impl SharedScalars {
    fn new(expand: bool, dot_notation: bool) -> Self {
        Self {
            atoms: Vec::new(),
            positions: AHashMap::new(),
            dot_positions: AHashMap::new(),
            dots: Vec::new(),
            dot_notation,
            map: Arc::new(Vec::new()),
            expansions: AHashMap::new(),
            expand,
        }
    }

    /// One position per scalar up to dot spelling: `g(p,q)` from this frontier
    /// and `dot(p,q)` from a nested contraction are one variable, so their
    /// exponents add instead of forming distinct monomials.
    fn position(&mut self, atom: Atom) -> usize {
        if let Some(&position) = self.positions.get(&atom) {
            return position;
        }
        let dots = crate::shorthands::schoonschip::DotNormalizer::notation(atom.as_view());
        if let Some(&position) = self.dot_positions.get(&dots) {
            self.positions.insert(atom, position);
            return position;
        }
        let position = self.atoms.len();
        self.atoms.push(atom.clone());
        self.dots.push(dots.clone());
        self.dot_positions.insert(dots, position);
        self.positions.insert(atom, position);
        let mut variables = self.map.as_ref().clone();
        variables.push(PolyVariable::Temporary(position));
        self.map = Arc::new(variables);
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
        fn composite(view: AtomView<'_>) -> bool {
            match view {
                AtomView::Add(_) => true,
                AtomView::Mul(product) => product.iter().any(composite),
                AtomView::Pow(power) => {
                    let (base, exponent) = power.get_base_exp();
                    matches!(exponent, AtomView::Num(_)) && composite(base)
                }
                _ => false,
            }
        }
        // An upper bound of the expanded term count, from the structure only.
        fn bound(view: AtomView<'_>) -> Option<usize> {
            match view {
                AtomView::Add(sum) => sum
                    .iter()
                    .try_fold(0usize, |total, term| total.checked_add(bound(term)?)),
                AtomView::Mul(product) => product
                    .iter()
                    .try_fold(1usize, |total, factor| total.checked_mul(bound(factor)?)),
                AtomView::Pow(power) => {
                    let (base, exponent) = power.get_base_exp();
                    let exponent = u32::try_from(i64::try_from(exponent).ok()?).ok()?;
                    bound(base)?.checked_pow(exponent)
                }
                _ => Some(1),
            }
        }
        if self.expand
            && composite(atom.as_view())
            && bound(atom.as_view()).is_some_and(|terms| terms <= EXPANSION_WORK_LIMIT)
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
                ScalarPolynomial::new(&Q, Some(polynomial.nterms()), self.map.clone());
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

    /// Frontier closing: scalar values (closed trace words evaluated by the
    /// gamma kernel) as polynomials over this frontier's variable map; their
    /// products of sums are expanded by polynomial arithmetic, never as atoms.
    /// Every variable is registered before any polynomial is built, so all of
    /// them share the final map (later arithmetic never re-orders a grown
    /// map). The kernel writes a value in factored form, the same subword
    /// result repeated wherever the recurrence reused it: each value is
    /// converted level by level (sums merged, products multiplied, like terms
    /// combined at every level), each repeated sum or product once.
    fn polynomials_of(&mut self, values: &[Atom]) -> Vec<Option<ScalarPolynomial>> {
        fn leaves<'v>(view: AtomView<'v>, seen: &mut AHashSet<AtomView<'v>>, out: &mut Vec<Atom>) {
            match view {
                AtomView::Num(_) => {}
                AtomView::Add(_) | AtomView::Mul(_) => {
                    if !seen.insert(view) {
                        return;
                    }
                    let children: Vec<AtomView<'v>> = match view {
                        AtomView::Add(add) => add.iter().collect(),
                        AtomView::Mul(product) => product.iter().collect(),
                        _ => unreachable!(),
                    };
                    for child in children {
                        leaves(child, seen, out);
                    }
                }
                AtomView::Pow(power) if natural_exponent(view).is_some() => {
                    leaves(power.get_base_exp().0, seen, out);
                }
                _ => {
                    if seen.insert(view) {
                        out.push(view.to_owned());
                    }
                }
            }
        }
        fn natural_exponent(view: AtomView<'_>) -> Option<usize> {
            let AtomView::Pow(power) = view else {
                return None;
            };
            let exponent = i64::try_from(power.get_base_exp().1).ok()?;
            usize::try_from(exponent).ok()
        }
        fn convert<'v>(
            view: AtomView<'v>,
            zero: &ScalarPolynomial,
            positions: &AHashMap<Atom, usize>,
            cache: &mut AHashMap<AtomView<'v>, ScalarPolynomial>,
        ) -> Option<ScalarPolynomial> {
            match view {
                AtomView::Num(_) => Some(zero.constant(Rational::try_from(view).ok()?)),
                AtomView::Add(_) | AtomView::Mul(_) => {
                    if let Some(known) = cache.get(&view) {
                        return Some(known.clone());
                    }
                    let polynomial = match view {
                        AtomView::Add(add) => {
                            let mut sum = CoefficientSum::One(zero.clone());
                            for term in add.iter() {
                                sum.push(
                                    convert(term, zero, positions, cache)?,
                                    #[cfg(feature = "reference-cases")]
                                    0,
                                );
                            }
                            sum.finish(
                                #[cfg(feature = "reference-cases")]
                                0,
                            )
                        }
                        AtomView::Mul(product) => {
                            let mut out = zero.one();
                            for factor in product.iter() {
                                out = &out * &convert(factor, zero, positions, cache)?;
                            }
                            out
                        }
                        _ => unreachable!(),
                    };
                    cache.insert(view, polynomial.clone());
                    Some(polynomial)
                }
                AtomView::Pow(power) if natural_exponent(view).is_some() => {
                    let base = convert(power.get_base_exp().0, zero, positions, cache)?;
                    Some(base.pow(natural_exponent(view)?))
                }
                _ => {
                    let mut exponents = vec![0u32; zero.nvars()];
                    exponents[*positions.get(&view.to_owned())?] = 1;
                    Some(zero.monomial(Rational::one(), exponents))
                }
            }
        }
        let mut found = Vec::new();
        {
            let mut seen = AHashSet::new();
            for value in values {
                leaves(value.as_view(), &mut seen, &mut found);
            }
        }
        for leaf in found {
            self.position(leaf);
        }
        let zero = ScalarPolynomial::new(&Q, None, self.map.clone());
        values
            .iter()
            .map(|value| {
                budget_checkpoint();
                // one value's repeated subwords; values do not share them
                let mut cache = AHashMap::new();
                convert(value.as_view(), &zero, &self.positions, &mut cache)
            })
            .collect()
    }

    /// Frontier closing: spell the routed-sum stand-ins of closed values as
    /// their sums, in every variable that holds one (once per variable,
    /// instead of rewriting the emitted expression). The values stay compact
    /// over the stand-ins: a scalar product with a routed sum is emitted as one
    /// factor (`dot(k1 + k2, p)` normalised by its linearity), not expanded
    /// into the frontier's own dots. Measured (2026-11-19, L4 B->B q1 and the
    /// L5 sample): spelling the sums at the leaves BEFORE the conversion gives
    /// a fully expanded output, but costs more in the contraction, the
    /// hand-over and the lowering (which reads products of sums sparsely) than
    /// the compact form; values are equal either way.
    fn substitute_spellings(&mut self, stand_ins: &[(Atom, Atom)]) {
        if stand_ins.is_empty() {
            return;
        }
        for spelling in self.atoms.iter_mut().chain(self.dots.iter_mut()) {
            let substituted =
                crate::dirac::closing::substitute_stand_ins(spelling.clone(), stand_ins);
            if substituted != *spelling {
                *spelling = substituted;
            }
        }
    }

    /// Emit each variable in the spelling it first had in this frontier, or
    /// in dot notation when the request ends there. The polynomial is
    /// relabelled with those atoms and written into one expression buffer,
    /// which is normalized once.
    fn emit(&self, mut polynomial: ScalarPolynomial) -> Option<Atom> {
        let spellings = if self.dot_notation {
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
        labelled.coefficients = std::mem::take(&mut polynomial.coefficients);
        labelled.exponents = std::mem::take(&mut polynomial.exponents);
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
}

// A parse belongs to one arithmetic scope. Opening an occurrence is cached on
// that scope, so borrowed payloads and its port atoms outlive every DP state.
// This is parser scratch, not another symbolic tensor representation.
struct ContractionScope<'src> {
    source: AtomView<'src>,
    network: SymbolicNet<AbstractIndex, AtomOrView<'src>>,
    tree: SimpleTraversalTree<ChildVecStore<()>>,
    root: NodeIndex,
    ports: AHashMap<NodeIndex, (Vec<Atom>, PartialStructure)>,
    opened: AHashMap<NodeIndex, std::cell::OnceCell<Option<Box<Self>>>>,
    reductions: AHashMap<NodeIndex, std::cell::OnceCell<Option<FactorizedContraction>>>,
    // Each scope serves one fixed contraction or admission policy.
    work: AHashMap<NodeIndex, std::cell::OnceCell<bool>>,
    interfaces: Option<LeafInterfaceCache>,
}

impl<'src> ContractionScope<'src> {
    fn parse(source: AtomView<'src>) -> Option<Self> {
        let network = SlotContraction::shallow_network(source, None).ok()?;
        Self::from_network(source, network, None)
    }

    fn from_network(
        source: AtomView<'src>,
        network: SymbolicNet<AbstractIndex, AtomOrView<'src>>,
        interfaces: Option<LeafInterfaceCache>,
    ) -> Option<Self> {
        let tree = network.graph.expr_tree().cast::<ChildVecStore<()>>();
        let root = network.graph.graph.node_id(network.graph.head());
        let nodes = tree
            .iter_preorder_tree_nodes(&network.graph.graph, root)
            .collect::<Vec<_>>();
        let external: SuBitGraph = network.graph.graph.external_filter();
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
            for child in tree.iter_children(node, &network.graph.graph) {
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
                if let NetworkNode::Leaf(NetworkLeaf::LocalTensor(index)) =
                    &network.graph.graph[node]
                    && let Some((_, ports)) = interfaces.get(
                        network.store.tensors[*index]
                            .expression
                            .as_atom_view()
                            .get_data(),
                    )
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
                for hedge in network.graph.graph.iter_crown(member) {
                    let NetworkEdge::Slot(slot) = network.graph.graph[[&hedge]] else {
                        continue;
                    };
                    if !external.includes(&hedge)
                        && positions
                            .get(&network.graph.graph.node_id(network.graph.graph.inv(hedge)))
                            .is_some_and(|&other| begin <= other && other < end)
                    {
                        continue;
                    }
                    // Sewing shares one edge label; the occurrence's storage
                    // axis still owns this endpoint's variance.
                    let slot = match &network.graph.graph[member] {
                        NetworkNode::Leaf(NetworkLeaf::LocalTensor(index)) => {
                            *network.store.tensors[*index]
                                .structure
                                .external_structure()
                                .get(usize::from(network.graph.slot_order[hedge.0]))?
                        }
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

    fn payload(&self, node: NodeIndex) -> Option<&AtomOrView<'src>> {
        match self.network.graph.graph[node] {
            NetworkNode::Leaf(NetworkLeaf::LocalTensor(index)) => {
                Some(&self.network.store.tensors[index].expression)
            }
            NetworkNode::Leaf(NetworkLeaf::Scalar(index)) => {
                Some(self.network.store.get_scalar_ref(index))
            }
            _ => None,
        }
    }

    fn literal(&self, node: NodeIndex) -> Option<AtomView<'_>> {
        self.payload(node)
            .map(AtomOrView::as_view)
            .or_else(|| (node == self.root).then_some(self.source))
    }

    fn open(&self, node: NodeIndex) -> Option<&Self> {
        self.opened[&node]
            .get_or_init(|| {
                // Composite arithmetic leaves retain their original view. The
                // parser's owned scalar coefficients contain no tensor interior.
                let AtomOrView::View(value) = self.payload(node)? else {
                    return None;
                };
                #[cfg(test)]
                crate::tensor::SELECTED_OPENS.with(|count| count.set(count.get() + 1));
                let network =
                    SlotContraction::shallow_network(*value, self.interfaces.as_ref()).ok()?;
                let opened = Self::from_network(*value, network, self.interfaces.clone())?;
                if matches!(
                    opened.network.graph.graph[opened.root],
                    NetworkNode::Leaf(_)
                ) && opened.literal(opened.root)? == *value
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
        settings: crate::tensor::ContractSettings<'_>,
    ) -> Option<bool> {
        let scope = self.scope(source)?;
        self.minimal_scope_work(&scope, scope.root, settings)
    }

    fn minimal_scope_work(
        &self,
        scope: &ContractionScope<'_>,
        root: NodeIndex,
        settings: crate::tensor::ContractSettings<'_>,
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
        if let Some((planned, graph)) = &self.graph
            && planned.as_view() == source
        {
            let network = graph.map_ref(
                |scalar| AtomOrView::View(scalar.as_view()),
                |tensor| SymbolicTensor {
                    expression: AtomOrView::View(tensor.expression.as_view()),
                    structure: tensor.structure.clone(),
                    is_metric: tensor.is_metric,
                    is_composite: tensor.is_composite,
                    proofs: tensor.proofs.clone(),
                },
            );
            ContractionScope::from_network(source, network, self.interfaces.clone())
        } else {
            ContractionScope::parse(source)
        }
    }

    /// The shared depth-one structural view keeps arithmetic leaves borrowed.
    /// Family collection opens selected leaves through its existing tape reader.
    pub(crate) fn shallow_network<'src>(
        source: AtomView<'src>,
        interfaces: Option<&LeafInterfaceCache>,
    ) -> Result<
        SymbolicNet<AbstractIndex, AtomOrView<'src>>,
        crate::tensor::inference::TensorInferenceError,
    > {
        #[cfg(feature = "reference-cases")]
        let _phase = crate::reference_cases::timing::scope(
            crate::reference_cases::timing::Phase::ShallowGraph,
        );
        let settings = ParseSettings {
            precontract_scalars: false,
            depth_limit: Some(1),
            shorthand_parsing: ShorthandParsing::Opaque,
            parse_composite_scalars_as_tensors: true,
            ..ParseSettings::default()
        };
        type Tensor<'src> = SymbolicTensor<OrderedStructure, AtomOrView<'src>>;
        let library = DummyLibrary::<Tensor<'src>>::new();
        if let Some(interfaces) = interfaces {
            SymbolicNet::<AbstractIndex, AtomOrView<'src>>::try_from_admitted_view::<
                OrderedStructure,
                _,
            >(source, &library, &settings, Arc::clone(interfaces))
        } else {
            SymbolicNet::<AbstractIndex, AtomOrView<'src>>::try_from_view::<OrderedStructure, _>(
                source, &library, &settings,
            )
        }
        .map_err(|error| crate::tensor::inference::TensorInferenceError::Invalid(error.to_string()))
    }

    /// Collect notation on connected scopes of the same shallow graph. Scalar
    /// spectators and components without a selected representation are reused.
    pub(crate) fn collect_matrix_notation(
        &self,
        source: AtomView<'_>,
        representations: &[LibraryRep],
        collect: bool,
        traces: bool,
    ) -> Option<Atom> {
        use crate::shorthands::chain::Chain;
        let scope = self.scope(source)?;
        let nodes = scope.factors(scope.root);
        let components = scope.network.graph.slot_components(&scope.tree, &nodes);
        let external: SuBitGraph = scope.network.graph.graph.external_filter();
        let mut output = Vec::new();
        let mut changed = false;
        for component in components {
            let factors = component
                .iter()
                .map(|&node| scope.literal(node))
                .collect::<Option<Vec<_>>>()?;
            let mut expression = if factors.len() == 1 {
                factors[0].to_owned()
            } else {
                Atom::mul_many(factors.iter())
            };
            for &representation in representations {
                let connection = component.len() > 1
                    // An opaque arithmetic leaf can own connections that the
                    // outer graph intentionally does not expose. Regional
                    // representation facts below select its local notation.
                    || factors.iter().any(|factor| {
                        matches!(factor, AtomView::Add(_) | AtomView::Mul(_) | AtomView::Pow(_))
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
                    let next = expression.collect_chains(representation, collect, traces, false);
                    changed |= next != expression;
                    expression = next;
                }
            }
            output.push(expression);
        }
        Some(if changed {
            Atom::mul_many(output)
        } else {
            source.to_owned()
        })
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

    /// Plan an already admitted interface; the typed caller owns source validation.
    pub(crate) fn contract_factorized(
        &self,
        source: AtomView<'_>,
        order: Option<&[usize]>,
        rank_one: bool,
    ) -> Option<FactorizedContraction> {
        if rank_one && self.observations.is_none() && source.contains_symbol(self.tags.chain) {
            self.dummy_state().reserve_indices(source);
        }
        let intrinsic = self
            .graph
            .as_ref()
            .filter(|(planned, _)| planned.as_view() == source)
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
        self.contract_graph(&scope, scope.root, order, rank_one)
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
        use spenso::network::parsing::AtomStructureExt;
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
                    if matches!(scope.literal(child)?, AtomView::Fun(function)
                        if function.get_symbol() == self.tags.bracket)
                        && scope.contains_work(child, self, self.contraction_settings(rank_one))
                    {
                        let opened = scope.open(child)?;
                        let mut result =
                            self.contract_graph(opened, opened.root, None, rank_one)?;
                        result.root = result.root.pow(power);
                        result.observations = result.observations.filter(scalar_facts);
                        return Some(result);
                    }
                    return Some(FactorizedContraction {
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
                result.observations = result.observations.filter(scalar_facts);
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
                    .collect::<std::collections::BTreeSet<_>>()
            })
            .collect::<Vec<_>>();
        let mut boundary = std::collections::BTreeSet::new();
        let mut planned = Vec::with_capacity(selected.len());
        while !selected.is_empty() {
            let best = (0..selected.len()).min_by_key(|&i| {
                let position = selected[i];
                (
                    boundary.symmetric_difference(&edges[position]).count(),
                    std::cmp::Reverse(boundary.intersection(&edges[position]).count()),
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
                        sum.opaque_factors[position].0 = result.root.as_view();
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
                            use crate::reference_cases::timing::{ScopeOpenReason, scope_open};
                            let reason = match region.as_ref() {
                                None => ScopeOpenReason::Unobserved,
                                Some(region) if region.contraction_sources(self.contraction_settings(sum.rank_one)) => ScopeOpenReason::MetricOrVector,
                                Some(_) => ScopeOpenReason::OpaqueOrBracket,
                            };
                            scope_open(reason, region.as_ref().is_some_and(|region| region.has_indexed_powers()));
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
                            if let Some(observed) = result.observations.as_ref().filter(|observed| scalar_facts(observed)) {
                                sum.scalar_forms.insert(result.root.as_view(), observed.needs_dot_normalization());
                            }
                            return Some(if scope.ports[&node].0.is_empty() {
                                TermLeaf::Value(InputLeaf::Scalar(result.root.as_view()))
                            } else {
                                let position = sum.register_factor(scope, node, positions)?;
                                sum.opaque_factors[position].0 = result.root.as_view();
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
                root: source.to_owned(),
                status: ReductionStatus::Complete,
                observations: self
                    .observations
                    .as_ref()
                    .and_then(|observed| observed.certify_scalar_region(source)),
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
                    .map(|(needs_dots, observed)| needs_dots || observed.needs_dot_normalization());
                terms.push(result.root);
            }
            return Some(FactorizedContraction {
                root: Atom::add_many(terms),
                status,
                observations: scalar_form
                    .filter(|_| scalar_boundary)
                    .map(DomainObservations::closed_scalar),
            });
        }
        let nodes = scope.factors(root);
        let selected = nodes
            .iter()
            .enumerate()
            .filter_map(|(i, &node)| {
                // Frontier closing (opt-in, inert when unset): a trace word is
                // never selected; it rides along with its bound ports and is
                // closed per state at emission.
                let word = crate::dirac::closing::frontier_active()
                    && scope
                        .literal(node)
                        .is_some_and(crate::dirac::closing::is_word);
                (!word && scope.contains_work(node, self, self.contraction_settings(rank_one)))
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
        let mut terminal = Some((Vec::new(), false));
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
                (Some((slots, needs_dots)), Some(observed), None) => {
                    Some((slots, needs_dots || observed.needs_dot_normalization()))
                }
                (Some((mut slots, needs_dots)), None, Some((open, open_dots))) => {
                    slots.extend(open);
                    Some((slots, needs_dots || open_dots))
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
        let root = if status == ReductionStatus::Complete && !changed {
            let source = scope.literal(root)?;
            // The bracket normalizer owns whether a product can be exposed:
            // closed scalar wrappers disappear, while AUTO operand order stays.
            if matches!(source, AtomView::Fun(function)
                if function.get_symbol() == self.tags.bracket)
            {
                crate::shorthands::bracket::BracketNormalizer::normalize(source)
            } else {
                source.to_owned()
            }
        } else {
            Atom::mul_many(roots)
        };
        let observations = if scalar_boundary {
            scalar_form
                .filter(|_| status == ReductionStatus::Complete)
                .map(DomainObservations::closed_scalar)
        } else {
            // Every component is a certified scalar polynomial, and the open
            // vectors and metrics sit exactly on the root's ports: no internal
            // connection and no identity candidate remains. Certify the open
            // output as terminal, so that its replacement need not be rescanned.
            terminal
                .filter(|_| status == ReductionStatus::Complete && changed)
                .and_then(|(mut slots, needs_dots)| {
                    slots.sort();
                    let mut ports = scope.ports[&root_node].0.clone();
                    ports.sort();
                    (slots == ports).then(|| {
                        DomainObservations::terminal_polynomial(
                            &scope.ports[&root_node].1,
                            needs_dots,
                        )
                    })
                })
        };
        Some(FactorizedContraction {
            root,
            status,
            observations,
        })
    }

    /// Execute one connected component on the admitted graph. Scratch bindings,
    /// compiled factors and tentative proofs cannot escape a refused component.
    #[allow(clippy::type_complexity)]
    fn contract_component<'a>(
        &self,
        scope: &'a ContractionScope<'_>,
        nodes: &[NodeIndex],
        selected: &[usize],
        exposed_sums: &AHashSet<NodeIndex>,
        rank_one: bool,
    ) -> Option<(FactorizedContraction, bool, Option<(Vec<Atom>, bool)>)> {
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
            let mut factors = Vec::with_capacity(remaining.len());
            for (position, arguments) in remaining {
                scalar_form = scalar_form
                    .zip(sum.scalar_form(sum.opaque_factors[position].0))
                    .filter(|_| sum.opaque_factors[position].1.is_empty())
                    .map(|(first, second)| first || second);
                factors.push(sum.emit_factor(position, &arguments)?);
            }
            return Some((
                FactorizedContraction {
                    root: Atom::mul_many(factors),
                    status: ReductionStatus::Complete,
                    observations: scalar_form.map(DomainObservations::closed_scalar),
                },
                false,
                None,
            ));
        }
        let mut scoped = AHashMap::new();
        let mut compile = |sum: &mut ComponentSum<'a, '_>, position: usize| {
            let mut tape = std::mem::take(&mut sum.input);
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
        let changed = sum.contracted || !scoped.is_empty();
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
        use crate::shorthands::schoonschip::{SchoonschipSettings, SchoonschipWithSettings};
        let mut settings = SchoonschipSettings::default().with_chain_like_functions();
        settings.schoonschip_rank1_tensors = rank_one;
        settings.metrics = self.metrics;
        settings.representations = self.representations.clone();
        let rewrite = SchoonschipWithSettings {
            settings: &settings,
        };
        match source {
            AtomView::Add(sum) => Atom::add_many(
                sum.iter()
                    .map(|term| rewrite.run_observed(term, observed, self)),
            ),
            _ => rewrite.run_observed(source, observed, self),
        }
    }
}

impl<'a> ComponentSum<'a, '_> {
    /// Storage of one [`frontier_sum`]: its contribution records and the
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
                std::mem::size_of::<Contribution>().saturating_add(
                    contribution
                        .shift
                        .capacity()
                        .saturating_mul(std::mem::size_of::<u32>()),
                )
            })
            .fold(0usize, usize::saturating_add);
        terms
            .saturating_mul(
                std::mem::size_of::<Rational>()
                    .saturating_add(nvars.saturating_mul(std::mem::size_of::<u32>())),
            )
            .saturating_add(records)
    }

    /// Reuse the admitted inventory; never inspect a completed coefficient to
    /// rediscover whether its source contains an identity or an index.
    /// The explicit slots of a vector or metric variable, whose only open
    /// indices are these slots; `None` for any other variable.
    fn open_slots(&self, variable: &Variable<'a>) -> Option<Vec<Atom>> {
        let explicit = |slot: AtomView<'_>| {
            spenso::structure::slot::Slot::<
                spenso::structure::representation::LibraryRep,
                spenso::structure::abstract_index::AbstractIndex,
            >::try_from(slot)
            .is_ok()
            .then(|| slot.to_owned())
        };
        match variable {
            Variable::Vector(_, slot) => Some(vec![explicit(*slot)?]),
            Variable::Metric([first, second]) => Some(vec![explicit(*first)?, explicit(*second)?]),
            _ => None,
        }
    }

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
        self.factor_roots.push(None);
        positions.insert(key, position);
        Some(position)
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
        let _states_timing = crate::reference_cases::timing::scope(
            crate::reference_cases::timing::Phase::ContractionStates,
        );
        let mut states = vec![(remaining, ScalarPolynomial::new_one(&Q))];
        let mut budget_checks = 0usize;
        let mut emitted: AHashMap<usize, Atom> = AHashMap::new();
        #[cfg(feature = "reference-cases")]
        let mut emitted_bytes = 0usize;
        let mut scalar_variables: AHashMap<Atom, usize> = AHashMap::new();
        let mut shared =
            SharedScalars::new(self.contractor.expand_scalars, self.contractor.dot_notation);
        let mut scalar_form = Some(false);
        #[cfg(feature = "reference-cases")]
        let mut generated = 0usize;
        let trace = *FRONTIER_TRACE;
        let trace_start = std::time::Instant::now();
        if trace {
            eprintln!(
                "FRONTIER begin order={} factors={} rss_mib={}",
                order.len(),
                states[0].0.factors.len(),
                trace_rss_mib()
            );
        }
        for (trace_step, &selected) in order.iter().enumerate() {
            let step_start = std::time::Instant::now();
            if self.factor_roots[selected].is_none() {
                self.factor_roots[selected] = Some(compile(self, selected)?);
            }
            let terms = self.factor_terms(selected)?;
            let trace_in = trace.then(|| {
                (
                    states.len(),
                    states.iter().map(|(_, w)| w.nterms()).sum::<usize>(),
                )
            });
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
            crate::reference_cases::timing::frontier(crate::reference_cases::timing::Frontier {
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
                let degrees = max_degrees(weight);
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
                    // Cooperative resource budgets (none by default) see a
                    // growing frontier here; a product of many index-sharing
                    // sums can otherwise grow without bound before any other
                    // checkpoint runs. The frontier cannot decline once
                    // started, so an exhausted budget aborts the operation.
                    budget_checks = budget_checks.wrapping_add(1);
                    if budget_checks.is_multiple_of(4096) {
                        budget_checkpoint();
                    }
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
                            unreachable!()
                        };
                        let atom = emitted.get(index)?;
                        let slots = &mut *self.slots;
                        let composite = match shared
                            .expansion(*index, atom, |atom| index_free(slots, atom.as_view()))
                        {
                            Expansion::Base(position) => Err(*position),
                            Expansion::Composite(expanded) => {
                                Ok(small_power(expanded.nterms(), *exponent))
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
                    let mut monomial = ScalarPolynomial::new(&Q, Some(1), shared.map.clone())
                        .monomial(self.coefficient.clone(), full);
                    for (index, exponent) in composites {
                        let Expansion::Composite(expanded) = &shared.expansions[&index] else {
                            unreachable!()
                        };
                        for _ in 0..exponent {
                            monomial = &monomial * expanded;
                        }
                    }
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
            let trace_contributions = if trace {
                next.iter().map(|(_, c, _)| c.len()).sum::<usize>()
            } else {
                0
            };
            let trace_generate = step_start.elapsed();
            let mut sources: Vec<Option<ScalarPolynomial>> = std::mem::take(&mut states)
                .into_iter()
                .map(|(_, weight)| Some(weight))
                .collect();
            #[cfg(feature = "reference-cases")]
            let mut next_bytes = 0usize;
            for (remaining, mut contributions, small) in next {
                budget_checkpoint();
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
                        crate::reference_cases::timing::frontier_workspace(
                            surrounding_bytes.saturating_add(Self::merge_bytes(
                                &sources,
                                &contributions,
                                shared.map.len(),
                            )),
                        );
                        frontier_sum(&mut sources, &mut uses, &contributions, &shared.map)?
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
            if let Some((in_states, in_terms)) = trace_in {
                let terms_out = states.iter().map(|(_, w)| w.nterms()).collect::<Vec<_>>();
                let bytes = states
                    .iter()
                    .map(|(_, w)| trace_poly_bytes(w))
                    .sum::<usize>();
                eprintln!(
                    "FRONTIER step={trace_step} selected={selected} local_terms={} in_states={in_states} in_terms={in_terms} out_states={} out_terms={} max_state_terms={} nvars={} coeff_bytes={bytes} contributions={trace_contributions} gen_s={:.3} merge_s={:.3} rss_mib={}",
                    terms.len(),
                    states.len(),
                    terms_out.iter().sum::<usize>(),
                    terms_out.iter().max().copied().unwrap_or(0),
                    shared.map.len(),
                    trace_generate.as_secs_f64(),
                    (step_start.elapsed() - trace_generate).as_secs_f64(),
                    trace_rss_mib(),
                );
            }
        }
        self.overrides.clear();
        let emit_start = std::time::Instant::now();
        let mut trace_emit_poly = std::time::Duration::ZERO;
        if trace {
            eprintln!(
                "FRONTIER states_done states={} terms={} nvars={} states_s={:.3} rss_mib={}",
                states.len(),
                states.iter().map(|(_, w)| w.nterms()).sum::<usize>(),
                shared.map.len(),
                trace_start.elapsed().as_secs_f64(),
                trace_rss_mib()
            );
            if std::env::var_os("IDENSO_FRONTIER_TRACE_VARS").is_some() {
                for (i, a) in shared.atoms.iter().enumerate() {
                    eprintln!("FRONTIER var {i} {}", a.to_canonical_string());
                }
            }
        }
        #[cfg(feature = "reference-cases")]
        let _emission_timing = crate::reference_cases::timing::scope(
            crate::reference_cases::timing::Phase::ContractionEmission,
        );
        let mut terms = Vec::with_capacity(states.len());
        // Open slots of vectors and metrics left in every state: the same
        // sorted explicit slots in each state, or no terminal certificate.
        let mut open_slots: Option<Vec<Atom>> = None;
        let mut any_open = false;
        let mut terminal = true;
        // Frontier closing (opt-in): every distinct binding of a word is
        // evaluated once by the gamma kernel. The states whose words all
        // closed and that keep no tensor are grouped by their closed values:
        // one group writes (sum of its weights) x (its values), after the
        // variable map is complete. The routed sums stay stand-ins until the
        // variables are spelled.
        let closing = crate::dirac::closing::frontier_active();
        let mut words: AHashMap<(usize, Vec<Argument<'a>>), Result<usize, Atom>> = AHashMap::new();
        let mut values: Vec<Atom> = Vec::new();
        let mut groups: Vec<(Vec<usize>, CoefficientSum)> = Vec::new();
        let mut group_of: AHashMap<Vec<usize>, usize> = AHashMap::new();
        let mut stand_ins: Vec<(Atom, Atom)> = Vec::new();
        // a closed value written beside a tensor: its stand-ins are
        // substituted in the root
        let mut loose_values = false;
        for (remaining, weight) in states {
            let remaining_count = remaining.factors.len();
            let mut closed = Vec::new();
            let mut extra = Vec::new();
            let mut state_slots = Vec::new();
            for (position, arguments) in remaining.factors {
                // A word closed by this state's bound vectors is a scalar
                // (conservatively in need of dot normalisation). A word still
                // open is emitted as today (no certificate).
                if closing && crate::dirac::closing::is_word(self.opaque_factors[position].0) {
                    let key = (position, arguments);
                    let evaluated = match words.get(&key) {
                        Some(evaluated) => evaluated.clone(),
                        None => {
                            let word = self.emit_factor(position, &key.1)?;
                            budget_checkpoint();
                            let evaluated =
                                match crate::dirac::closing::evaluate_closed_word(word.as_view()) {
                                    Some((value, opaque)) => {
                                        for pair in opaque {
                                            // two routed sums under one stand-in
                                            // name (a 64-bit hash collision):
                                            // decline, the source is kept
                                            if stand_ins.iter().any(|(stand_in, sum)| {
                                                *stand_in == pair.0 && *sum != pair.1
                                            }) {
                                                return None;
                                            }
                                            if !stand_ins.contains(&pair) {
                                                stand_ins.push(pair);
                                            }
                                        }
                                        values.push(value);
                                        Ok(values.len() - 1)
                                    }
                                    None => Err(word),
                                };
                            words.insert(key, evaluated.clone());
                            evaluated
                        }
                    };
                    match evaluated {
                        Ok(value) => closed.push(value),
                        Err(word) => {
                            scalar_form = None;
                            extra.push(word);
                        }
                    }
                    continue;
                }
                scalar_form = scalar_form
                    .zip(self.scalar_form(self.opaque_factors[position].0))
                    .filter(|_| self.opaque_factors[position].1.is_empty())
                    .map(|(first, second)| first || second);
                extra.push(self.emit_factor(position, &arguments)?);
            }
            if !closed.is_empty() {
                scalar_form = scalar_form.map(|_| true);
                if extra.is_empty() && remaining.residual.is_empty() {
                    closed.sort_unstable();
                    let group = match group_of.get(&closed) {
                        Some(&group) => group,
                        None => {
                            group_of.insert(closed.clone(), groups.len());
                            groups.push((
                                closed,
                                CoefficientSum::One(ScalarPolynomial::new_zero(&Q)),
                            ));
                            groups.len() - 1
                        }
                    };
                    groups[group].1.push(
                        weight,
                        #[cfg(feature = "reference-cases")]
                        0,
                    );
                    // no open slots: a closed state of the terminal bookkeeping
                    match &open_slots {
                        None => open_slots = Some(Vec::new()),
                        Some(slots) => terminal &= slots.is_empty(),
                    }
                    continue;
                }
                loose_values = true;
                extra.extend(closed.into_iter().map(|value| values[value].clone()));
            }
            let poly_start = std::time::Instant::now();
            let trace_terms = weight.nterms();
            let mut factors = vec![shared.emit(weight)?];
            trace_emit_poly += poly_start.elapsed();
            if trace {
                eprintln!(
                    "FRONTIER emit_state terms={trace_terms} remaining={} residual={} atom_bytes={} emit_s={:.3} rss_mib={}",
                    remaining_count,
                    remaining.residual.len(),
                    factors[0].as_view().get_byte_size(),
                    poly_start.elapsed().as_secs_f64(),
                    trace_rss_mib()
                );
            }
            factors.extend(extra);
            for (variable, exponent) in remaining.residual {
                if let Some(slots) = self.open_slots(&self.variables[variable]) {
                    // Each power repeats the slots, so a repeated (internal)
                    // index never matches the distinct root ports.
                    for _ in 0..exponent {
                        state_slots.extend(slots.iter().cloned());
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
        if !groups.is_empty() {
            let polynomials = shared.polynomials_of(&values);
            shared.substitute_spellings(&stand_ins);
            let mut total = CoefficientSum::One(ScalarPolynomial::new_zero(&Q));
            for (closed, weights) in groups {
                budget_checkpoint();
                let mut weight = weights.finish(
                    #[cfg(feature = "reference-cases")]
                    0,
                );
                if closed.iter().all(|&value| polynomials[value].is_some()) {
                    // one expanded polynomial of all closed states
                    for &value in &closed {
                        let mut value = polynomials[value].clone()?;
                        weight.unify_variables(&mut value);
                        weight = &weight * &value;
                    }
                    total.push(
                        weight,
                        #[cfg(feature = "reference-cases")]
                        0,
                    );
                    continue;
                }
                // a value the polynomial conversion declined: weights x value
                let mut factors = vec![shared.emit(weight)?];
                for value in closed {
                    factors.push(match &polynomials[value] {
                        Some(polynomial) => shared.emit(polynomial.clone())?,
                        None => crate::dirac::closing::substitute_stand_ins(
                            values[value].clone(),
                            &stand_ins,
                        ),
                    });
                }
                terms.push(Atom::mul_many(factors));
            }
            let total = total.finish(
                #[cfg(feature = "reference-cases")]
                0,
            );
            if !total.is_zero() {
                // A small closed sum is expanded once its routed sums are
                // spelled, so that an exact cancellation shows as zero, as on
                // the unclosed path; a large one stays compact (products of
                // sums, read sparsely downstream).
                let small = total.nterms() <= EXPANSION_TERM_LIMIT;
                let emitted = shared.emit(total)?;
                terms.push(if small && !stand_ins.is_empty() {
                    emitted.expand()
                } else {
                    emitted
                });
            }
        }
        self.terminal_slots = match (open_slots, scalar_form) {
            (Some(slots), Some(needs_dots)) if any_open && terminal => Some((slots, needs_dots)),
            _ => None,
        };
        let add_start = std::time::Instant::now();
        let root = if loose_values {
            crate::dirac::closing::substitute_stand_ins(Atom::add_many(terms), &stand_ins)
        } else {
            Atom::add_many(terms)
        };
        if trace {
            eprintln!(
                "FRONTIER emitted bytes={} emit_poly_s={:.3} add_many_s={:.3} emission_s={:.3} total_s={:.3} rss_mib={}",
                root.as_view().get_byte_size(),
                trace_emit_poly.as_secs_f64(),
                add_start.elapsed().as_secs_f64(),
                emit_start.elapsed().as_secs_f64(),
                trace_start.elapsed().as_secs_f64(),
                trace_rss_mib()
            );
        }
        Some(FactorizedContraction {
            root,
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

    fn emit_argument(&self, argument: Argument<'a>) -> Atom {
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
        let source = source.to_owned();
        let replacements = ports
            .iter()
            .zip(arguments)
            .enumerate()
            .filter(|(_, (port, argument))| **argument != Argument::Original(**port))
            .map(|(position, (_, &argument))| (position, self.emit_argument(argument)))
            .collect::<std::collections::HashMap<_, _>>();
        if replacements.is_empty() {
            return Some(source);
        }
        let source = SymbolicTensor::from_normalized_parts(source, interface.clone());
        let expression = crate::tensor::composition::rewrite_interface_ports(
            &source,
            &replacements,
            Some(self.contractor.dummy_state()),
        )
        .ok()?;
        Some(
            crate::shorthands::schoonschip::DotNormalizer::with_settings(
                expression.as_view(),
                crate::tensor::ContractSettings {
                    metrics: self.contractor.metrics,
                    rank_one: self.rank_one,
                    representations: self.contractor.representations.as_deref(),
                    ..Default::default()
                },
            ),
        )
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
        let sum = frontier_sum(&mut sources, &mut uses, &records, &map).unwrap();
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
            frontier_sum(&mut sources, &mut [1], &identity, &map).unwrap(),
            full
        );
        assert!(sources[0].is_none());

        // A source whose variables are not a prefix of the map is declined.
        let foreign =
            ScalarPolynomial::new(&Q, Some(1), Arc::new(vec![PolyVariable::Temporary(7)]))
                .monomial(Rational::one(), vec![1]);
        let mut sources = vec![Some(foreign)];
        assert!(frontier_sum(&mut sources, &mut [2], &identity, &map).is_none());
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
            use crate::reference_cases::timing::{self, Frontier};
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
    fn source_free_color_scopes_stay_opaque_beside_real_contractions() {
        let reps = crate::test_support::test_initialize();
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
        let contractor = SlotContraction::new()
            .expanding(false)
            .observed(source.reduction_observations());
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
        let reps = crate::test_support::test_initialize();
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
        let contractor = SlotContraction::new()
            .expanding(false)
            .observed(source.reduction_observations());
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
        crate::test_support::test_initialize();
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
        let contractor = SlotContraction::new().observed(source.reduction_observations());
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
        use spenso::network::parsing::AtomStructureExt;
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
        crate::test_support::test_initialize();
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
                let settings = crate::tensor::ContractSettings::default();
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
    use super::*;
    use spenso::network::tags::SPENSO_TAG;

    #[test]
    fn partial_arithmetic_scopes_have_only_one_level() {
        let _contractor = super::super::tests::setup();
        for (source, nodes) in [
            ("p(spenso::mink(4,a))*q(spenso::mink(4,a))+x*y", 3),
            ("(p(spenso::mink(4,a))*q(spenso::mink(4,a))+x)^2", 2),
            (
                "spenso::bracket(p(spenso::mink(4,a))*q(spenso::mink(4,a)))^-1",
                2,
            ),
        ] {
            let source = super::super::tests::input(source);
            let scope = ContractionScope::parse(source.as_view()).unwrap();
            assert_eq!(scope.ports.len(), nodes, "{source}");
            assert_eq!(scope.network.graph.graph.n_nodes(), nodes, "{source}");
            assert!(scope.opened.values().all(|entry| entry.get().is_none()));
        }
    }

    #[test]
    fn selected_arithmetic_reuses_admission_without_nested_graphs() {
        let setup = super::super::tests::setup();
        let input = super::super::tests::input;
        let scalar = Atom::add_many(
            (0..512).map(|i| Atom::var(symbolica::symbol!(&format!("selected_tape::s{i}")))),
        );
        let source = (&scalar * input("p(mink(4,a))") + input("q(mink(4,a))"))
            * input("t(mink(4,a),mink(4,b))");
        let tensor = SymbolicTensor::infer(source.clone()).unwrap();
        let contractor = setup.observed(tensor.reduction_observations());
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
        let setup = super::super::tests::setup();
        let input = super::super::tests::input;
        let source = input(
            "(p(mink(4,a))*q(mink(4,a))*r(mink(4,b))+s(mink(4,b)))*(p(mink(4,a))*q(mink(4,a))*t(mink(4,b),mink(4,c))+u(mink(4,b),mink(4,c)))",
        );
        let tensor = SymbolicTensor::infer(source.clone()).unwrap();
        let contractor = setup.observed(tensor.reduction_observations());
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
        use super::super::tests::{input, setup};
        let contractor = setup();
        let opaque = symbolica::symbol!("unobserved_power::opaque"; Scalar);
        for value in [
            input("m^2"),
            input("(m+n)^-2"),
            input("(m+n)^(1/3)"),
            opaque.call(input("(m+n)^11")).pow(2),
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
        let calls = std::sync::Arc::new(std::sync::atomic::AtomicUsize::new(0));
        let observed = std::sync::Arc::clone(&calls);
        let callback = symbolica::symbol!("unobserved_power::callback"; Scalar; norm=move|_,_| {
            observed.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
        });
        let rounded = Atom::num(symbolica::domains::float::Float::parse("0.1", Some(11)).unwrap());
        for value in [
            input("p(mink(4,a))^2"),
            input("(p(mink(4,a))*q(mink(4,a))+m)^2"),
            callback.call(input("m")).pow(2),
            (rounded + input("m")).pow(2),
        ] {
            calls.store(0, std::sync::atomic::Ordering::Relaxed);
            let mut slots = SlotMatcher::default();
            let mut sum = ComponentSum::new(&contractor, Intake::Contraction, &mut slots);
            let mut tape = TermTape::default();
            assert!(
                contractor
                    .compile_arithmetic(value.as_view(), &mut tape, &mut sum)
                    .is_none()
            );
            assert_eq!(calls.load(std::sync::atomic::Ordering::Relaxed), 0);
        }
    }

    #[test]
    fn scalar_powers_in_selected_sums_do_not_rebuild_graphs() {
        let setup = super::super::tests::setup();
        let input = super::super::tests::input;
        let source = input(
            "(m^4*p(mink(4,a))+m^2*g(q(mink(4)),s(mink(4)))^2*r(mink(4,a)))*t(mink(4,a),mink(4,b))",
        );
        let tensor = SymbolicTensor::infer(source.clone()).unwrap();
        let contractor = setup.observed(tensor.reduction_observations());
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
        let setup = super::super::tests::setup();
        let input = super::super::tests::input;
        let source = input("p(mink(4))+q(mink(4))");
        let tensor = SymbolicTensor::infer(source.clone()).unwrap();
        assert_eq!(tensor.structure.logical_slots().len(), 1);
        let contractor = setup.observed(tensor.reduction_observations());
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
        let setup = super::super::tests::setup();
        let input = super::super::tests::input;
        let source = input(
            "((p(mink(4,a))*q(mink(4,a))+x)^2*r(mink(4,b))+s(mink(4,b)))*t(mink(4,b),mink(4,c))",
        );
        let tensor = SymbolicTensor::infer(source.clone()).unwrap();
        let contractor = setup.observed(tensor.reduction_observations());
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
        let contractor = super::super::tests::setup();
        let source =
            super::super::tests::input("n*(1/(q(0,spenso::cind(0))+x)+1/(q(0,spenso::cind(0))+y))");
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
        let contractor = super::super::tests::setup();
        let source = super::super::tests::input(
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
        let expected =
            super::super::tests::input("spenso::g(p(spenso::mink(4)),q(spenso::mink(4)))");
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
        crate::test_support::test_initialize();
        SPENSO_TAG.rank_one_tensor_symbol("lazy_contraction::p");
        SPENSO_TAG.tensor_symbol("lazy_contraction::t");
        SPENSO_TAG.tensor_symbol("lazy_contraction::u");
        let parse = |source: &str| {
            Atom::parse(
                source,
                "lazy_contraction",
                symbolica::parser::ParseSettings::symbolica(),
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
