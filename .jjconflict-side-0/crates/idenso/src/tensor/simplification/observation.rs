//! Candidate facts retained by an immutable domain. The initial
//! inventory and regional observations come from the same syntactic traversal.

use super::super::SymbolicTensor;
use crate::{
    dirac::AGS,
    representations::{Bispinor, ColorAdjoint, ColorFundamental, ColorSextet},
    shorthands::schoonschip::SimplificationCandidates,
};
use spenso::structure::partial::{
    PartialIndex, PartialSlot, PartialStructure, PartialStructureExt,
};
use spenso::{
    network::{library::symbolic::ETS, tags::SPENSO_TAG},
    structure::{
        abstract_index::AbstractIndex,
        representation::LibraryRep,
        slot::{SlotMatch, SlotMatcher},
    },
};
use std::{
    collections::{HashMap, HashSet},
    sync::{Arc, LazyLock},
};
use symbolica::atom::{Atom, AtomCore, AtomView, Symbol};

pub(crate) const HEAD_COUNT: usize = 24;
pub(crate) static DOMAIN_HEADS: LazyLock<[Symbol; HEAD_COUNT]> = LazyLock::new(|| {
    [
        LibraryRep::from(Bispinor {}).symbol(),
        LibraryRep::from(ColorFundamental {}).symbol(),
        LibraryRep::from(ColorAdjoint {}).symbol(),
        LibraryRep::from(ColorSextet {}).symbol(),
        *crate::epsilon::EPSILON_SYMBOL,
        SPENSO_TAG.chain,
        SPENSO_TAG.trace,
        SPENSO_TAG.bracket,
        AGS.gamma,
        AGS.gamma0,
        AGS.gamma5,
        AGS.projm,
        AGS.projp,
        AGS.sigma,
        AGS.charge_conjugation,
        AGS.gammaconj,
        AGS.gammaadj,
        ETS.metric,
        crate::color::CS.t,
        crate::color::CS.f,
        crate::color::CS.d,
        crate::color::CS.gram,
        crate::color::CS.cas,
        crate::color::CS.idx,
    ]
});

#[cfg(test)]
thread_local! {
    pub(crate) static INITIAL_SCANS: std::cell::Cell<usize> = const { std::cell::Cell::new(0) };
    pub(crate) static INITIAL_DOMAINS: std::cell::RefCell<Vec<Atom>> = const { std::cell::RefCell::new(Vec::new()) };
    pub(crate) static REPLACEMENT_SCANS: std::cell::Cell<usize> = const { std::cell::Cell::new(0) };
    pub(crate) static REUSED_KERNEL_REGIONS: std::cell::Cell<usize> = const { std::cell::Cell::new(0) };
    pub(crate) static REUSED_REGIONS: std::cell::Cell<usize> = const { std::cell::Cell::new(0) };
}

#[derive(Clone, Debug)]
pub(crate) struct RegionObservation {
    pub(crate) counts: [usize; HEAD_COUNT],
    pub(crate) representations: HashSet<LibraryRep>,
    pub(crate) leaves: HashSet<Atom>,
    reserved_indices: HashSet<AbstractIndex>,
    indices: usize,
    max_indices: usize,
    metric_sources: HashSet<LibraryRep>,
    vector_sources: HashSet<LibraryRep>,
    powered: bool,
    epsilon_degree: u8,
    intrinsic: bool,
    opaque: bool,
    dots: bool,
}

impl Default for RegionObservation {
    fn default() -> Self {
        Self {
            counts: [0; HEAD_COUNT],
            representations: HashSet::new(),
            leaves: HashSet::new(),
            reserved_indices: HashSet::new(),
            indices: 0,
            max_indices: 0,
            metric_sources: HashSet::new(),
            vector_sources: HashSet::new(),
            powered: false,
            epsilon_degree: 0,
            intrinsic: true,
            opaque: false,
            dots: false,
        }
    }
}

impl RegionObservation {
    fn head(&mut self, head: Symbol) {
        for (count, candidate) in self.counts.iter_mut().zip(DOMAIN_HEADS.iter()) {
            *count += usize::from(head == *candidate);
        }
    }

    fn observe(&mut self, node: AtomView<'_>, slot: &SlotMatch<'_>, slots: &mut SlotMatcher) {
        match node {
            AtomView::Var(variable) => {
                self.head(variable.get_symbol());
            }
            AtomView::Fun(function) => {
                let head = function.get_symbol();
                self.head(head);
                self.epsilon_degree = u8::from(head == *crate::epsilon::EPSILON_SYMBOL);
                // Reuse these observed function boundaries for arbitrary
                // collector predicates as well as built-in identity heads.
                // Opaque slot metadata remains outside this traversal.
                self.leaves.insert(node.to_owned());
                self.intrinsic &=
                    super::super::inference::InterfaceInference::intrinsic_normalization_head(head);
                self.dots |= head.has_tag(&SPENSO_TAG.rank1);
                // Compact vectors have no explicit endpoint to substitute.
                // Retain actual source ports during the existing observation,
                // rather than treating every rank-one head as remaining work.
                if head == ETS.metric {
                    for argument in function.iter() {
                        if let SlotMatch::Explicit(port) = slots.classify(argument)
                            && let Ok(rep) = slots.representation(port)
                        {
                            self.metric_sources.insert(rep);
                        }
                    }
                } else if head.has_tag(&SPENSO_TAG.rank1)
                    && let Some(argument) = slots.vector_argument(function)
                    && let SlotMatch::Explicit(port) = slots.classify(argument)
                    && let Ok(rep) = slots.representation(port)
                {
                    self.vector_sources.insert(rep);
                }

                if let SlotMatch::Explicit(slot) = slot {
                    self.indices += 1;
                    self.max_indices += 1;
                    if let Ok(index) = AbstractIndex::try_from(slot.index()) {
                        self.reserved_indices.insert(index);
                    }
                    // The candidate walker keeps slot payloads opaque. Reserve
                    // nested written names there without treating them as heads.
                    for payload in [slot.dimension(), slot.index()] {
                        if !matches!(payload, AtomView::Var(_) | AtomView::Num(_)) {
                            self.reserve_metadata(payload, slots);
                        }
                    }
                    if let Ok(rep) = slots.representation(*slot) {
                        self.representations.insert(rep);
                    }
                    if slot.representation().wrapper().is_some() {
                        self.head(slot.representation().head());
                    }
                    for payload in [slot.dimension(), slot.index()] {
                        match payload {
                            AtomView::Var(variable) => self.head(variable.get_symbol()),
                            AtomView::Num(_) => {}
                            _ => self.opaque = true,
                        }
                    }
                } else if let Ok(rep) = slots.parse_representation::<LibraryRep>(node) {
                    self.representations.insert(rep.rep);
                } else if matches!(slot, SlotMatch::Opaque) {
                    self.reserve_metadata(node, slots);
                    self.opaque |= function
                        .iter()
                        .any(|arg| !matches!(arg, AtomView::Var(_) | AtomView::Num(_)));
                }
            }
            AtomView::Pow(power) => {
                self.powered = true;
                self.dots |= matches!(power.get_base(), AtomView::Fun(base) if base.get_symbol() == ETS.metric);
            }
            _ => {}
        }
    }

    fn reserve_metadata(&mut self, value: AtomView<'_>, slots: &mut SlotMatcher) {
        value.visitor(&mut |node| {
            if let Ok(slot) = slots.parse::<LibraryRep, AbstractIndex>(node) {
                self.reserved_indices.insert(slot.aind);
            }
            true
        });
    }

    pub(crate) fn identity_candidates(&self) -> u8 {
        Self::identity_heads(&self.counts)
    }

    fn identity_heads(counts: &[usize; HEAD_COUNT]) -> u8 {
        let gamma =
            counts[8..17].iter().any(|&count| count > 0) || (counts[0] > 0 && counts[6] > 0);
        let color = counts[18..].iter().any(|&count| count > 0)
            || (counts[1..4].iter().any(|&count| count > 0) && counts[6] > 0);
        (u8::from(gamma) * super::GAMMA)
            | (u8::from(color) * super::COLOR)
            | (u8::from(counts[4] > 0) * super::EPSILON)
    }

    fn merge(&mut self, other: &Self) {
        for (count, additional) in self.counts.iter_mut().zip(other.counts) {
            *count += additional;
        }
        self.representations
            .extend(other.representations.iter().copied());
        self.leaves.extend(other.leaves.iter().cloned());
        self.reserved_indices
            .extend(other.reserved_indices.iter().copied());
        self.indices += other.indices;
        self.max_indices += other.max_indices;
        self.metric_sources
            .extend(other.metric_sources.iter().copied());
        self.vector_sources
            .extend(other.vector_sources.iter().copied());
        self.powered |= other.powered;
        self.epsilon_degree = (self.epsilon_degree + other.epsilon_degree).min(2);
        self.intrinsic &= other.intrinsic;
        self.opaque |= other.opaque;
        self.dots |= other.dots;
    }

    pub(crate) fn has_indices(&self) -> bool {
        self.indices > 0
    }

    pub(crate) fn structural_work(&self) -> bool {
        self.indices > 0 || self.dots || self.counts[7] > 0 || self.counts[17] > 0
    }

    fn certifies_scalar_syntax(&self) -> bool {
        self.intrinsic
            && !self.opaque
            && self.identity_candidates() == 0
            && self.counts[5..8].iter().all(|&count| count == 0)
            && !self.has_indices()
            && self.reserved_indices.is_empty()
    }
}

#[derive(Clone, Debug)]
pub(crate) struct DomainObservations {
    pub(crate) candidates: SimplificationCandidates<HEAD_COUNT>,
    regions: HashMap<Atom, Arc<RegionObservation>>,
    scopes: HashMap<Atom, Arc<RegionObservation>>,
    // A producer can establish scalar closure without collecting syntactic
    // regions. Missing regions remain unknown to arbitrary collector predicates.
    scalar_form: Option<ScalarPolynomialForm>,
    pub(crate) representations: HashSet<LibraryRep>,
    pub(crate) epsilon_degree: u8,
    max_indices: usize,
    pub(crate) reserved_indices: HashSet<AbstractIndex>,
}

#[derive(Clone, Copy, Debug)]
enum ScalarPolynomialForm {
    Metrics,
    Dots,
}

impl DomainObservations {
    /// Certify intrinsic scalar syntax with no local identity or contraction work.
    /// The producer owns this proof, including its compact vector metadata.
    pub(crate) fn closed_scalar(needs_dots: bool) -> Self {
        Self {
            candidates: SimplificationCandidates {
                repeated_indices: false,
                brackets: false,
                dots: needs_dots,
                symbols: [false; HEAD_COUNT],
                counts: [0; HEAD_COUNT],
                complete: false,
                traversal_complete: false,
                intrinsic: true,
            },
            regions: HashMap::new(),
            scopes: HashMap::new(),
            scalar_form: Some(if needs_dots {
                ScalarPolynomialForm::Metrics
            } else {
                ScalarPolynomialForm::Dots
            }),
            representations: HashSet::new(),
            epsilon_degree: 0,
            max_indices: 0,
            reserved_indices: HashSet::new(),
        }
    }

    pub(crate) fn is_closed_scalar(&self) -> bool {
        self.scalar_form.is_some()
    }

    pub(crate) fn needs_dot_normalization(&self) -> bool {
        match self.scalar_form {
            Some(ScalarPolynomialForm::Metrics) => true,
            Some(ScalarPolynomialForm::Dots) => false,
            None => self.candidates.dots || self.candidates.symbols[17],
        }
    }

    /// Establish local scalar closure from observations and a caller-proven scalar
    /// interface.
    pub(crate) fn can_certify_scalar(&self) -> bool {
        self.is_closed_scalar()
            || (self.candidates.complete
                && self.candidates.traversal_complete
                && self.candidates.intrinsic
                && self
                    .regions
                    .values()
                    .all(|region| region.certifies_scalar_syntax()))
    }

    /// Reuse an observed subtree when composition has separately established its
    /// scalar interface. Unobserved or unrelated regions are never certified.
    pub(crate) fn certify_scalar_region(&self, expression: AtomView<'_>) -> Option<Self> {
        let region = self.region(expression)?;
        region
            .certifies_scalar_syntax()
            .then(|| Self::closed_scalar(region.dots || region.counts[17] > 0))
    }

    fn regions(expression: AtomView<'_>) -> Vec<AtomView<'_>> {
        match expression {
            AtomView::Mul(product) => product.iter().collect(),
            AtomView::Add(sum) => sum.iter().collect(),
            _ => vec![expression],
        }
    }

    fn initial(expression: AtomView<'_>) -> Self {
        #[cfg(test)]
        INITIAL_SCANS.with(|count| count.set(count.get() + 1));
        #[cfg(test)]
        INITIAL_DOMAINS.with(|domains| domains.borrow_mut().push(expression.to_owned()));
        Self::scan_regions(expression)
    }

    fn scan_regions(expression: AtomView<'_>) -> Self {
        let wanted = Self::regions(expression)
            .into_iter()
            .map(|region| region.to_owned())
            .collect::<HashSet<_>>();
        let mut regions = HashMap::new();
        let mut slots = SlotMatcher::default();
        let mut scopes = HashMap::new();
        let mut stack = Vec::<(usize, usize, Atom, RegionObservation)>::new();
        let candidates = SimplificationCandidates::scan_observing(
            expression,
            *DOMAIN_HEADS,
            || true,
            |node, slot| {
                // Contiguous subtree bytes provide exit events during the same
                // traversal. The stack retains no full arithmetic tree.
                let begin = node.get_data().as_ptr() as usize;
                let end = begin + node.get_data().len();
                while stack
                    .last()
                    .is_some_and(|(start, stop, _, _)| begin < *start || end > *stop)
                {
                    Self::finish_observation(&mut stack, &wanted, &mut regions, &mut scopes);
                }
                let mut observation = RegionObservation::default();
                observation.observe(node, slot, &mut slots);
                stack.push((begin, end, node.to_owned(), observation));
            },
        );
        while !stack.is_empty() {
            Self::finish_observation(&mut stack, &wanted, &mut regions, &mut scopes);
        }
        Self::from_regions(
            candidates,
            regions,
            scopes,
            matches!(expression, AtomView::Add(_)),
        )
    }

    fn finish_observation(
        stack: &mut Vec<(usize, usize, Atom, RegionObservation)>,
        wanted: &HashSet<Atom>,
        regions: &mut HashMap<Atom, Arc<RegionObservation>>,
        scopes: &mut HashMap<Atom, Arc<RegionObservation>>,
    ) {
        let (_, _, value, mut observed) = stack.pop().unwrap();
        // A power can contribute multiple epsilon copies. Retaining a candidate
        // for unsupported exponents is conservative; the kernel decides legality.
        if matches!(value.as_view(), AtomView::Pow(_)) && observed.epsilon_degree > 0 {
            observed.epsilon_degree = 2;
        }
        if let Some((_, _, parent, aggregate)) = stack.last_mut() {
            let degree = if matches!(parent.as_view(), AtomView::Add(_)) {
                aggregate.epsilon_degree.max(observed.epsilon_degree)
            } else {
                (aggregate.epsilon_degree + observed.epsilon_degree).min(2)
            };
            let max_indices = if matches!(parent.as_view(), AtomView::Add(_)) {
                aggregate.max_indices.max(observed.max_indices)
            } else {
                aggregate.max_indices + observed.max_indices
            };
            aggregate.merge(&observed);
            aggregate.epsilon_degree = degree;
            aggregate.max_indices = max_indices;
        }
        let observed = Arc::new(observed);
        if wanted.contains(&value) {
            regions.insert(value.clone(), Arc::clone(&observed));
        }
        if matches!(
            value.as_view(),
            AtomView::Add(_) | AtomView::Mul(_) | AtomView::Pow(_)
        ) || matches!(value.as_view(), AtomView::Fun(function) if function.get_symbol() == SPENSO_TAG.bracket)
        {
            scopes.insert(value, observed);
        }
    }

    fn from_regions(
        candidates: SimplificationCandidates<HEAD_COUNT>,
        regions: HashMap<Atom, Arc<RegionObservation>>,
        scopes: HashMap<Atom, Arc<RegionObservation>>,
        additive: bool,
    ) -> Self {
        let epsilon_degree = if additive {
            regions
                .values()
                .map(|region| region.epsilon_degree)
                .max()
                .unwrap_or(0)
        } else {
            regions.values().fold(0u8, |degree, region| {
                (degree + region.epsilon_degree).min(2)
            })
        };
        let max_indices = if additive {
            regions
                .values()
                .map(|region| region.max_indices)
                .max()
                .unwrap_or(0)
        } else {
            regions.values().map(|region| region.max_indices).sum()
        };
        let reserved_indices = regions
            .values()
            .flat_map(|region| region.reserved_indices.iter().copied())
            .collect();
        let representations = regions
            .values()
            .flat_map(|region| region.representations.iter().copied())
            .collect();
        Self {
            candidates,
            regions,
            scopes,
            scalar_form: None,
            representations,
            epsilon_degree,
            max_indices,
            reserved_indices,
        }
    }

    /// A kernel supplies the changed replacement. Reuse each surviving region;
    /// only newly written regions lose their syntactic and interface facts.
    pub(crate) fn updated(&self, expression: AtomView<'_>) -> Self {
        let mut regions = HashMap::new();
        let mut scopes = self.scopes.clone();
        for view in Self::regions(expression) {
            let region = if let Some(region) = self.regions.get(view.get_data()) {
                #[cfg(test)]
                REUSED_REGIONS.with(|count| count.set(count.get() + 1));
                Arc::clone(region)
            } else {
                #[cfg(test)]
                REPLACEMENT_SCANS.with(|count| count.set(count.get() + 1));
                let replacement = Self::scan_regions(view);
                let mut region = RegionObservation::default();
                for part in replacement.regions.values() {
                    region.merge(part);
                }
                region.epsilon_degree = replacement.epsilon_degree;
                region.max_indices = replacement.max_indices;
                scopes.extend(replacement.scopes);
                Arc::new(region)
            };
            regions.insert(view.to_owned(), region);
        }
        let mut candidates = self.candidates.clone();
        candidates.counts =
            std::array::from_fn(|i| regions.values().map(|region| region.counts[i]).sum());
        candidates.symbols = candidates.counts.map(|count| count > 0);
        candidates.complete = regions.values().all(|region| !region.opaque);
        candidates.intrinsic = regions
            .values()
            .all(|region| region.intrinsic && !region.opaque);
        candidates.traversal_complete = true;
        // Index incidence is a graph fact. Conservatively invalidate it after a
        // replacement; the contraction graph decides compatibility and scope.
        candidates.repeated_indices =
            regions.values().map(|region| region.indices).sum::<usize>() > 1;
        candidates.brackets = candidates.symbols[7];
        candidates.dots = regions.values().any(|region| region.dots);
        Self::from_regions(
            candidates,
            regions,
            scopes,
            matches!(expression, AtomView::Add(_)),
        )
    }

    /// No permitted metric/vector endpoint remains in this observed domain.
    pub(crate) fn excludes_contraction_sources(
        &self,
        settings: crate::tensor::ContractSettings<'_>,
    ) -> bool {
        if self.is_closed_scalar() {
            return true;
        }
        self.candidates.complete
            && self.candidates.traversal_complete
            && self.regions.values().all(|region| {
                !(settings.metrics
                    && region
                        .metric_sources
                        .iter()
                        .any(|rep| settings.permits(*rep)))
                    && !(settings.rank_one
                        && region
                            .vector_sources
                            .iter()
                            .any(|rep| settings.permits(*rep)))
            })
    }

    pub(crate) fn identity_candidates(&self) -> u8 {
        if self.is_closed_scalar() {
            0
        } else {
            RegionObservation::identity_heads(&self.candidates.counts)
        }
    }

    pub(crate) fn region(&self, expression: AtomView<'_>) -> Option<&RegionObservation> {
        self.regions
            .get(expression.get_data())
            .or_else(|| self.scopes.get(expression.get_data()))
            .map(AsRef::as_ref)
    }
}

impl SymbolicTensor<PartialStructure> {
    pub(crate) fn reduction_observations(&self) -> &Arc<DomainObservations> {
        self.proofs.observations.get_or_init(|| {
            let mut observed = DomainObservations::initial(self.expression.as_view());
            self.exclude_external_connections(&mut observed);
            Arc::new(observed)
        })
    }

    fn exclude_external_connections(&self, observed: &mut DomainObservations) {
        // Admission/checked rewriting proves the same external boundary for
        // every sum branch. If even the branch with most explicit occurrences
        // has only those external ports, none contains a dummy connection.
        // Repeated spelling across alternatives is not a contraction. Powers
        // can encode additional copies, so retain their candidate.
        if !observed.candidates.complete
            || self.proofs.validated.get() != Some(&true)
            || observed.regions.values().any(|region| region.powered)
        {
            return;
        }
        let external = self
            .structure
            .logical_slots()
            .iter()
            .filter(|slot| matches!(slot.aind, PartialIndex::Explicit(_)))
            .count();
        if observed.max_indices == external {
            observed.candidates.repeated_indices = false;
        }
    }

    pub(crate) fn observe_replacement(&self, previous: &Self) {
        // Exact unchanged opaque boundaries retain their logical interfaces.
        // A callback can invalidate these facts, so only intrinsic, admitted
        // transformations share the cache with a new payload.
        if self.proofs.validated.get() == Some(&true)
            && let Some(interfaces) = previous.proofs.leaf_interfaces.get()
            && previous.normalization_is_intrinsic()
        {
            let _ = self.proofs.leaf_interfaces.set(Arc::clone(interfaces));
        }
        if self.proofs.observations.get().is_some() {
            return;
        }
        let previous_observations = previous.reduction_observations();
        let observed = if self.expression == previous.expression {
            Arc::clone(previous_observations)
        } else {
            let mut observed = previous_observations.updated(self.expression.as_view());
            self.exclude_external_connections(&mut observed);
            Arc::new(observed)
        };
        let _ = self.proofs.observations.set(observed);
    }
}

/// Completed, unchanged intrinsic kernel inputs within one domain and family.
/// Rewritten outputs are never cached: their fresh bindings belong to one use.
#[derive(Default)]
pub(crate) struct SettledRegions {
    inputs: HashSet<(Atom, Vec<PartialSlot>)>,
    bytes: usize,
}

impl SettledRegions {
    pub(crate) fn run(
        &mut self,
        source: SymbolicTensor<PartialStructure>,
        complete: bool,
        operation: impl FnOnce(
            SymbolicTensor<PartialStructure>,
        )
            -> Result<SymbolicTensor<PartialStructure>, super::TensorInferenceError>,
    ) -> Result<SymbolicTensor<PartialStructure>, super::TensorInferenceError> {
        let key = (source.expression.clone(), source.structure.logical_slots());
        if complete && self.inputs.contains(&key) {
            #[cfg(test)]
            REUSED_KERNEL_REGIONS.with(|count| count.set(count.get() + 1));
            return Ok(source);
        }
        let intrinsic = complete && source.normalization_is_intrinsic();
        let result = operation(source)?;
        if intrinsic
            && result.proofs.frontier == super::ReductionStatus::Complete
            && result.expression == key.0
            && result.structure.logical_slots() == key.1
        {
            let bytes = key.0.as_view().get_byte_size().saturating_add(
                key.1
                    .len()
                    .saturating_mul(std::mem::size_of::<PartialSlot>()),
            );
            if self.bytes.saturating_add(bytes)
                <= spenso::shadowing::TermTape::<usize>::MAX_GENERATED_FACTOR_BYTES
                && self.inputs.insert(key)
            {
                self.bytes += bytes;
            }
        }
        Ok(result)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn scalar_certificate_excludes_work_without_fabricating_regional_inventory() {
        crate::test_support::test_initialize();
        let observed = DomainObservations::closed_scalar(true);
        assert!(observed.is_closed_scalar());
        assert!(observed.needs_dot_normalization());
        assert_eq!(observed.identity_candidates(), 0);
        assert!(observed.excludes_contraction_sources(Default::default()));
        assert!(!observed.candidates.complete);
        assert!(!observed.candidates.traversal_complete);

        let replacement = spenso::p!(spenso::mink!(4, 98925));
        let changed = observed.updated(replacement.as_view());
        assert!(!changed.is_closed_scalar());
        assert!(!changed.can_certify_scalar());
        assert!(!changed.excludes_contraction_sources(Default::default()));
        assert!(changed.region(replacement.as_view()).is_some());
    }

    #[test]
    fn scalar_region_certification_reuses_only_the_selected_spectator() {
        crate::test_support::test_initialize();
        let x = Atom::var(symbolica::symbol!("scalar_region::x"));
        let spectator = &x + Atom::one();
        let tensor = spenso::p!(spenso::mink!(4, 98926));
        let expression = &spectator * &tensor;
        let observed = DomainObservations::initial(expression.as_view());
        let scans = INITIAL_SCANS.with(|count| count.get());
        let certified = observed.certify_scalar_region(spectator.as_view()).unwrap();
        assert!(certified.is_closed_scalar());
        assert!(!certified.needs_dot_normalization());
        assert_eq!(INITIAL_SCANS.with(|count| count.get()), scans);
        assert!(observed.certify_scalar_region(tensor.as_view()).is_none());
        assert!(observed.certify_scalar_region(x.as_view()).is_none());
    }

    #[test]
    fn changed_scalar_factor_reuses_sum_incidence_without_an_initial_rescan() {
        crate::test_support::test_initialize();
        let a = spenso::mink!(4, 98921);
        let b = spenso::mink!(4, 98922);
        let c = spenso::mink!(4, 98923);
        let d = spenso::mink!(4, 98924);
        let x = Atom::var(symbolica::symbol!("incidence_reuse::x"));
        let y = Atom::var(symbolica::symbol!("incidence_reuse::y"));
        let free_sum =
            spenso::g!(&a, &b) * spenso::g!(&c, &d) + spenso::g!(&a, &c) * spenso::g!(&b, &d);
        let source = SymbolicTensor::infer(&x * &free_sum).unwrap();
        let before = source.reduction_observations();
        let leaf = Arc::clone(before.regions.get(free_sum.as_view().get_data()).unwrap());
        let initial_scans = INITIAL_SCANS.with(|count| count.get());
        let changed = source.with_identity_result(&y * &free_sum, None).unwrap();
        changed.observe_replacement(&source);
        let after = changed.reduction_observations();
        assert_eq!(INITIAL_SCANS.with(|count| count.get()), initial_scans);
        assert!(Arc::ptr_eq(
            &leaf,
            after.regions.get(free_sum.as_view().get_data()).unwrap()
        ));
        assert_eq!(after.max_indices, 4);
        assert!(!after.candidates.repeated_indices);
        assert!(
            after
                .regions
                .values()
                .map(|region| region.indices)
                .sum::<usize>()
                > 4
        );

        // A scalar-interface leaf may still contain a genuine internal pair.
        // Taking the maximum over alternatives must retain that branch's work.
        let open_sum = spenso::g!(&c, &d) + spenso::p!(&a) * spenso::q!(&a) * spenso::g!(&c, &d);
        let source = SymbolicTensor::infer(&x * &open_sum).unwrap();
        source.reduction_observations();
        let changed = source.with_identity_result(&y * &open_sum, None).unwrap();
        changed.observe_replacement(&source);
        assert_eq!(changed.reduction_observations().max_indices, 4);
        assert!(changed.reduction_observations().candidates.repeated_indices);
    }
}
