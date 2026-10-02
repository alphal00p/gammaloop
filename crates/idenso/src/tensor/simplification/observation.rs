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
        representation::{LibraryRep, RepName},
        slot::{IsAbstractSlot, SlotMatch, SlotMatcher},
    },
};
use std::{
    collections::{HashMap, HashSet},
    sync::{Arc, LazyLock, RwLock},
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
    static OBSERVED_NODES: std::cell::Cell<usize> = const { std::cell::Cell::new(0) };
    static REUSED_SCOPES: std::cell::Cell<usize> = const { std::cell::Cell::new(0) };
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct RegionObservation {
    pub(crate) counts: [usize; HEAD_COUNT],
    pub(crate) representations: HashSet<LibraryRep>,
    // Possible matrix channels, including existing chain endpoints. Merely
    // mentioning a representation in a vector or scalar dot is insufficient.
    chain_channels: HashSet<LibraryRep>,
    // Only functions observed directly in this scope. Descendant inventories
    // stay in their shared regional observations rather than being cloned and
    // rehashed into every arithmetic ancestor.
    leaves: HashSet<Atom>,
    leaf_children: Vec<Arc<RegionObservation>>,
    leaf_inventory_complete: bool,
    reserved_indices: HashSet<AbstractIndex>,
    indices: usize,
    max_indices: usize,
    metric_sources: HashSet<LibraryRep>,
    vector_sources: HashSet<LibraryRep>,
    indexed_power: bool,
    scalar_interface: bool,
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
            chain_channels: HashSet::new(),
            leaves: HashSet::new(),
            leaf_children: Vec::new(),
            leaf_inventory_complete: true,
            reserved_indices: HashSet::new(),
            indices: 0,
            max_indices: 0,
            metric_sources: HashSet::new(),
            vector_sources: HashSet::new(),
            indexed_power: false,
            scalar_interface: true,
            epsilon_degree: 0,
            intrinsic: true,
            opaque: false,
            dots: false,
        }
    }
}

impl RegionObservation {
    pub(crate) fn any_leaf(&self, mut selected: impl FnMut(AtomView<'_>) -> bool) -> Option<bool> {
        if !self.leaf_inventory_complete {
            return None;
        }
        let mut pending = vec![self];
        let mut regions = HashSet::new();
        let mut leaves = HashSet::new();
        while let Some(region) = pending.pop() {
            if !regions.insert(std::ptr::from_ref(region)) {
                continue;
            }
            for leaf in &region.leaves {
                // Arbitrary predicates retain the previous set semantics even
                // when an unchanged child occurs at several call sites.
                if leaves.insert(leaf.as_view()) && selected(leaf.as_view()) {
                    return Some(true);
                }
            }
            pending.extend(region.leaf_children.iter().map(AsRef::as_ref));
        }
        Some(false)
    }

    pub(crate) fn has_chain_channel(&self, representation: LibraryRep) -> bool {
        self.chain_channels.contains(&representation.base())
    }

    pub(crate) fn contraction_sources(
        &self,
        settings: crate::tensor::ContractSettings<'_>,
    ) -> bool {
        (settings.metrics && self.metric_sources.iter().any(|rep| settings.permits(*rep)))
            || (settings.rank_one && self.vector_sources.iter().any(|rep| settings.permits(*rep)))
    }

    /// A complete regional inventory can rule out metric/vector work without
    /// opening its arithmetic. Brackets can expose additional endpoints when
    /// unfolded, so their existing contraction path remains authoritative.
    pub(crate) fn excludes_contraction_sources(
        &self,
        settings: crate::tensor::ContractSettings<'_>,
    ) -> bool {
        // A terminal producer can certify sources without enumerating leaves.
        // Its surviving ports conservatively populate both source inventories.
        !self.opaque && self.counts[7] == 0 && !self.contraction_sources(settings)
    }

    fn head(&mut self, head: Symbol) {
        for (count, candidate) in self.counts.iter_mut().zip(DOMAIN_HEADS.iter()) {
            *count += usize::from(head == *candidate);
        }
    }

    fn observe(&mut self, node: AtomView<'_>, slot: &SlotMatch<'_>, slots: &mut SlotMatcher) {
        #[cfg(test)]
        OBSERVED_NODES.with(|count| count.set(count.get() + 1));
        match node {
            AtomView::Var(variable) => {
                let head = variable.get_symbol();
                self.head(head);
                self.scalar_interface &=
                    !head.has_tag(&SPENSO_TAG.tensor) && !head.has_tag(&SPENSO_TAG.representation);
            }
            AtomView::Fun(function) => {
                let head = function.get_symbol();
                self.head(head);
                self.scalar_interface &= head.is_scalar()
                    || !(head.has_tag(&SPENSO_TAG.tensor)
                        || head.has_tag(&SPENSO_TAG.rank1)
                        || head.has_tag(&SPENSO_TAG.representation)
                        || [
                            SPENSO_TAG.chain,
                            SPENSO_TAG.trace,
                            SPENSO_TAG.bracket,
                            SPENSO_TAG.dot,
                            ETS.metric,
                        ]
                        .contains(&head));
                self.epsilon_degree = (self.epsilon_degree
                    + u8::from(head == *crate::epsilon::EPSILON_SYMBOL))
                .min(2);
                // Reuse these observed function boundaries for arbitrary
                // collector predicates as well as built-in identity heads.
                // Opaque slot metadata remains outside this traversal.
                self.leaves.insert(node.to_owned());
                self.intrinsic &=
                    super::super::inference::InterfaceInference::intrinsic_normalization_head(head);
                let rank_one = head.has_tag(&SPENSO_TAG.rank1);
                self.dots |= rank_one;
                let vector_argument = rank_one.then(|| slots.vector_argument(function)).flatten();
                if head != ETS.metric
                    && vector_argument.is_none()
                    && function.get_nargs() >= 2
                    && matches!(slot, SlotMatch::Other)
                {
                    // Observe direct ports while this function is already being
                    // visited. Index payloads remain opaque; only representation,
                    // variance and dimension constrain a possible matrix channel.
                    // Compact ports can become explicit when the shared notation
                    // owner opens a composite dot, so retain those candidates.
                    let mut ports = Vec::new();
                    for argument in function.iter() {
                        // Supplied rank-one tensors still occupy matrix endpoints.
                        let endpoint = slots.port_representation(argument);
                        let Some((representation, dimension)) = endpoint else {
                            continue;
                        };
                        for &(previous, previous_dimension) in &ports {
                            if dimension == previous_dimension && representation.matches(&previous)
                            {
                                self.chain_channels.insert(representation.base());
                            }
                        }
                        ports.push((representation, dimension));
                    }
                }
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
                } else if let Some(argument) = vector_argument
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
                    // The complete index is reserved above. Its internal spelling
                    // is neither another index nor a tensor candidate.
                    if !matches!(slot.dimension(), AtomView::Var(_) | AtomView::Num(_)) {
                        self.reserve_metadata(slot.dimension(), slots);
                    }
                    if let Ok(rep) = slots.representation(*slot) {
                        self.representations.insert(rep);
                    }
                    if slot.representation().wrapper().is_some() {
                        let head = slot.representation().head();
                        self.head(head);
                        self.intrinsic &=
                            super::super::inference::InterfaceInference::intrinsic_normalization_head(head);
                    }
                    self.opaque |= !matches!(slot.dimension(), AtomView::Var(_) | AtomView::Num(_));
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

    fn merge(&mut self, other: &Arc<Self>) {
        for (count, additional) in self.counts.iter_mut().zip(other.counts) {
            *count += additional;
        }
        self.representations
            .extend(other.representations.iter().copied());
        self.chain_channels
            .extend(other.chain_channels.iter().copied());
        if !other.leaves.is_empty() || !other.leaf_children.is_empty() {
            self.leaf_children.push(Arc::clone(other));
        }
        self.leaf_inventory_complete &= other.leaf_inventory_complete;
        self.reserved_indices
            .extend(other.reserved_indices.iter().copied());
        self.indices += other.indices;
        self.max_indices += other.max_indices;
        self.metric_sources
            .extend(other.metric_sources.iter().copied());
        self.vector_sources
            .extend(other.vector_sources.iter().copied());
        self.indexed_power |= other.indexed_power;
        self.scalar_interface &= other.scalar_interface;
        self.epsilon_degree = (self.epsilon_degree + other.epsilon_degree).min(2);
        self.intrinsic &= other.intrinsic;
        self.opaque |= other.opaque;
        self.dots |= other.dots;
    }

    pub(crate) fn has_indices(&self) -> bool {
        self.indices > 0
    }

    pub(crate) fn has_indexed_powers(&self) -> bool {
        self.indexed_power
    }

    pub(crate) fn structural_work(&self) -> bool {
        self.indices > 0 || self.dots || self.counts[7] > 0 || self.counts[17] > 0
    }

    pub(crate) fn certifies_scalar_interface_syntax(&self) -> bool {
        self.intrinsic
            && !self.opaque
            && self.scalar_interface
            && self.counts[5..8].iter().all(|&count| count == 0)
            && !self.has_indices()
            && self.reserved_indices.is_empty()
    }

    fn certifies_scalar_syntax(&self) -> bool {
        self.certifies_scalar_interface_syntax() && self.identity_candidates() == 0
    }
}

#[derive(Clone, Debug)]
pub(crate) struct DomainObservations {
    pub(crate) candidates: SimplificationCandidates<HEAD_COUNT>,
    regions: HashMap<Atom, Arc<RegionObservation>>,
    // Exact syntax facts survive replacements independently of live candidate
    // multiplicities. Share this memo instead of cloning every unchanged scope.
    scopes: Arc<RwLock<HashMap<Atom, Arc<RegionObservation>>>>,
    // A producer can establish terminal polynomial form without collecting syntactic
    // regions. Missing regions remain unknown to arbitrary collector predicates.
    polynomial_form: Option<TerminalPolynomialForm>,
    // Candidate eligibility can be complete even when a kernel certified a
    // terminal region without enumerating its syntactic leaves.
    identity_inventory_complete: bool,
    // The admitted interface accounts for every explicit incidence. Unlike
    // terminal polynomial form, this says nothing about remaining identities.
    excludes_internal_connections: bool,
    pub(crate) representations: HashSet<LibraryRep>,
    pub(crate) epsilon_degree: u8,
    max_indices: usize,
    pub(crate) reserved_indices: HashSet<AbstractIndex>,
}

#[derive(Clone, Copy, Debug)]
enum TerminalPolynomialForm {
    Metrics,
    Dots,
}

impl DomainObservations {
    /// The regional inventory was collected with the initial candidate walk or
    /// replacement observation. Unknown producer facts must remain conservative.
    pub(crate) fn has_chain_channel(&self, representation: LibraryRep) -> bool {
        !self.is_terminal_polynomial()
            && (!self.candidates.traversal_complete
                || self
                    .regions
                    .values()
                    .any(|region| region.has_chain_channel(representation)))
    }

    /// Certify intrinsic scalar syntax with no local identity or contraction work.
    /// The producer owns this proof, including its compact vector metadata.
    pub(crate) fn closed_scalar(needs_dots: bool) -> Self {
        Self::terminal_polynomial(&PartialStructure::from_logical_slots([]), needs_dots)
    }

    /// A kernel has removed all internal connections and identity candidates.
    /// The remaining free ports are still needed for subsequent composition and
    /// dummy reservation. This proof cannot be propagated through a composition
    /// that connects those ports, or through an uncontrolled rewrite.
    ///
    /// Unlike an initial observation, this does not claim a complete syntactic
    /// inventory. Arbitrary collector predicates must inspect unknown regions.
    pub(crate) fn terminal_polynomial(structure: &PartialStructure, needs_dots: bool) -> Self {
        let slots = structure.logical_slots();
        let reserved_indices = slots
            .iter()
            .filter_map(|slot| match slot.aind {
                PartialIndex::Explicit(index) => Some(index),
                _ => None,
            })
            .collect::<HashSet<_>>();
        let max_indices = slots
            .iter()
            .filter(|slot| matches!(slot.aind, PartialIndex::Explicit(_)))
            .count();
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
            scopes: Arc::default(),
            polynomial_form: Some(if needs_dots {
                TerminalPolynomialForm::Metrics
            } else {
                TerminalPolynomialForm::Dots
            }),
            identity_inventory_complete: true,
            excludes_internal_connections: true,
            representations: slots.iter().map(IsAbstractSlot::rep_name).collect(),
            epsilon_degree: 0,
            max_indices,
            reserved_indices,
        }
    }

    pub(crate) fn is_terminal_polynomial(&self) -> bool {
        self.polynomial_form.is_some()
    }

    fn terminal_region(&self) -> Arc<RegionObservation> {
        debug_assert!(self.is_terminal_polynomial());
        Arc::new(RegionObservation {
            representations: self.representations.clone(),
            reserved_indices: self.reserved_indices.clone(),
            indices: self.max_indices,
            max_indices: self.max_indices,
            // Free polynomial ports can acquire partners during composition.
            // The producer excludes internal connections, not those new paths.
            metric_sources: self.representations.clone(),
            vector_sources: self.representations.clone(),
            scalar_interface: self.is_closed_scalar(),
            dots: self.needs_dot_normalization(),
            leaf_inventory_complete: false,
            ..Default::default()
        })
    }

    /// Retain a produced terminal region beside exactly observed spectators.
    /// Normalization must preserve its opaque boundary; newly merged factors
    /// decline this shortcut. Arbitrary predicates still inspect the terminal
    /// region because its leaf inventory has deliberately not been collected.
    pub(crate) fn compose_terminal_product(
        expression: AtomView<'_>,
        structure: &PartialStructure,
        terminal_expression: AtomView<'_>,
        terminal: &Self,
        spectators: &[SymbolicTensor<PartialStructure>],
    ) -> Option<Self> {
        if !terminal.is_terminal_polynomial()
            || spectators
                .iter()
                .any(|value| !value.normalization_is_intrinsic())
        {
            return None;
        }
        let observed = spectators
            .iter()
            .map(|value| (value.expression.as_view(), value.reduction_observations()))
            .collect::<Vec<_>>();
        let factors = match expression {
            AtomView::Mul(product) => product.iter().collect::<Vec<_>>(),
            _ => vec![expression],
        };
        let mut regions = HashMap::new();
        let mut retained = false;
        for factor in factors {
            let region = if factor == terminal_expression {
                retained = true;
                terminal.terminal_region()
            } else {
                observed.iter().find_map(|(body, facts)| {
                    facts.region(factor).or_else(|| {
                        (factor == *body && facts.is_terminal_polynomial())
                            .then(|| facts.terminal_region())
                    })
                })?
            };
            regions.insert(factor.to_owned(), region);
        }
        if !retained {
            return None;
        }
        let counts = std::array::from_fn(|i| regions.values().map(|region| region.counts[i]).sum());
        let candidates = SimplificationCandidates {
            repeated_indices: regions.values().map(|region| region.indices).sum::<usize>() > 1,
            brackets: counts[7] > 0,
            dots: regions.values().any(|region| region.dots),
            symbols: counts.map(|count| count > 0),
            counts,
            complete: false,
            traversal_complete: false,
            intrinsic: true,
        };
        // Keep nested scalar-interface proofs when a later contraction opens
        // one of these unchanged spectators. The usual same-domain case shares
        // its cache outright; independently admitted spectators contribute only
        // their existing exact-syntax facts, without another candidate walk.
        let scopes = observed
            .first()
            .map(|(_, facts)| Arc::clone(&facts.scopes))
            .unwrap_or_default();
        for (_, facts) in observed.iter().skip(1) {
            if !Arc::ptr_eq(&scopes, &facts.scopes) {
                let inherited = facts
                    .scopes
                    .read()
                    .unwrap()
                    .iter()
                    .map(|(expression, region)| (expression.clone(), Arc::clone(region)))
                    .collect::<Vec<_>>();
                scopes.write().unwrap().extend(inherited);
            }
        }
        let mut result = Self::from_regions(candidates, regions, scopes, false);
        result.identity_inventory_complete = terminal.identity_inventory_complete
            && observed
                .iter()
                .all(|(_, facts)| facts.identity_inventory_complete);
        // The collector has established the combined logical interface. A rank
        // drop signals new cross-boundary connections and invalidates closure.
        let ports = terminal.max_indices
            + spectators
                .iter()
                .map(|value| value.structure.logical_slots().len())
                .sum::<usize>();
        result.excludes_internal_connections = ports == structure.logical_slots().len()
            && observed
                .iter()
                .all(|(_, facts)| facts.excludes_internal_connections());
        if result.excludes_internal_connections {
            result.candidates.repeated_indices = false;
        }
        Some(result)
    }

    pub(crate) fn excludes_internal_connections(&self) -> bool {
        self.excludes_internal_connections
    }

    pub(crate) fn is_closed_scalar(&self) -> bool {
        self.is_terminal_polynomial() && self.representations.is_empty()
    }

    pub(crate) fn needs_dot_normalization(&self) -> bool {
        match self.polynomial_form {
            Some(TerminalPolynomialForm::Metrics) => true,
            Some(TerminalPolynomialForm::Dots) => false,
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
        Self::scan_regions(expression, None)
    }

    fn is_scope(value: AtomView<'_>) -> bool {
        matches!(
            value,
            AtomView::Add(_) | AtomView::Mul(_) | AtomView::Pow(_)
        ) || matches!(value, AtomView::Fun(function)
            if function.get_symbol().is_scalar()
                || function.get_symbol().has_tag(&SPENSO_TAG.tensor)
                || function.get_symbol().has_tag(&SPENSO_TAG.rank1)
                || [SPENSO_TAG.chain, SPENSO_TAG.trace, SPENSO_TAG.bracket,
                    SPENSO_TAG.scalar, SPENSO_TAG.pure_scalar, SPENSO_TAG.dot, ETS.metric,
                    *spenso::shadowing::SYM, *spenso::shadowing::ANTISYM, *spenso::shadowing::CYCLIC]
                    .contains(&function.get_symbol()))
    }

    fn scan_regions(
        expression: AtomView<'_>,
        cached: Option<&HashMap<Atom, Arc<RegionObservation>>>,
    ) -> Self {
        let wanted = Self::regions(expression)
            .into_iter()
            .map(|region| region.get_data().as_ptr() as usize)
            .collect::<HashSet<_>>();
        let mut regions = HashMap::new();
        let mut slots = SlotMatcher::default();
        let mut scopes = HashMap::new();
        let mut stack = Vec::<(usize, usize, Atom, Arc<RegionObservation>)>::new();
        let candidates = SimplificationCandidates::scan_observing(
            expression,
            *DOMAIN_HEADS,
            || true,
            |node, slot| {
                // Contiguous subtree bytes provide exit events during the same
                // traversal. Only retained scopes need aggregate frames; ordinary
                // function arguments and slot metadata update their owner directly.
                let begin = node.get_data().as_ptr() as usize;
                let end = begin + node.get_data().len();
                while stack
                    .last()
                    .is_some_and(|(start, stop, _, _)| begin < *start || end > *stop)
                {
                    Self::finish_observation(&mut stack, &wanted, &mut regions, &mut scopes);
                }
                // Each additive branch needs a separate maximum for incidence
                // and epsilon degree, including branches that are simple leaves.
                let retain = stack.is_empty()
                    || wanted.contains(&begin)
                    || Self::is_scope(node)
                    || stack.last().is_some_and(|(_, _, parent, _)| {
                        matches!(parent.as_view(), AtomView::Add(_))
                    });
                if retain {
                    if let Some(observation) = cached.and_then(|cache| cache.get(node.get_data())) {
                        #[cfg(test)]
                        REUSED_SCOPES.with(|count| count.set(count.get() + 1));
                        stack.push((begin, end, node.to_owned(), Arc::clone(observation)));
                        // Stop the shared candidate/index traversal itself, not
                        // just regional bookkeeping. Replacement candidates are
                        // reconstructed from these complete regional summaries.
                        return false;
                    }
                    stack.push((begin, end, node.to_owned(), Arc::default()));
                }
                Arc::make_mut(&mut stack.last_mut().unwrap().3).observe(node, slot, &mut slots);
                true
            },
        );
        while !stack.is_empty() {
            Self::finish_observation(&mut stack, &wanted, &mut regions, &mut scopes);
        }
        Self::from_regions(
            candidates,
            regions,
            Arc::new(RwLock::new(scopes)),
            matches!(expression, AtomView::Add(_)),
        )
    }

    fn finish_observation(
        stack: &mut Vec<(usize, usize, Atom, Arc<RegionObservation>)>,
        wanted: &HashSet<usize>,
        regions: &mut HashMap<Atom, Arc<RegionObservation>>,
        scopes: &mut HashMap<Atom, Arc<RegionObservation>>,
    ) {
        let (begin, _, value, mut observed) = stack.pop().unwrap();
        if let AtomView::Fun(function) = value.as_view() {
            let head = function.get_symbol();
            // Scalar owners hide representation metadata from their interface.
            // A metric consumes the sole compact port of each vector in a dot.
            // Generic compact tensor leaves keep their unresolved interfaces.
            let scalar = head.is_scalar()
                || head == SPENSO_TAG.scalar
                || head == SPENSO_TAG.pure_scalar
                || ([ETS.metric, SPENSO_TAG.dot].contains(&head) && function.get_nargs() == 2 && {
                    let mut slots = SlotMatcher::default();
                    function.iter().all(|argument| {
                        matches!(argument, AtomView::Fun(vector)
                            if vector.get_symbol().has_tag(&SPENSO_TAG.rank1)
                                && slots.vector_argument(vector)
                                    .and_then(|port| slots.compact_representation(port))
                                    .is_some())
                    })
                });
            if scalar && !observed.scalar_interface {
                Arc::make_mut(&mut observed).scalar_interface = true;
            }
        }
        if matches!(value.as_view(), AtomView::Pow(_)) {
            // Scalar coefficient powers stay opaque to contraction. Anything
            // without that certificate retains the existing per-copy scope,
            // including indices hidden in metadata and tensor shorthands.
            let indexed_power = !observed.certifies_scalar_interface_syntax();
            if (indexed_power && !observed.indexed_power) || observed.epsilon_degree == 1 {
                let observed = Arc::make_mut(&mut observed);
                observed.indexed_power |= indexed_power;
                // A power can contribute multiple epsilon copies. Retaining a candidate
                // for unsupported exponents is conservative; the kernel decides legality.
                if observed.epsilon_degree == 1 {
                    observed.epsilon_degree = 2;
                }
            }
        }
        if let Some((_, _, parent, aggregate)) = stack.last_mut() {
            let aggregate = Arc::make_mut(aggregate);
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
        if wanted.contains(&begin) {
            regions.insert(value.clone(), Arc::clone(&observed));
        }
        if Self::is_scope(value.as_view()) {
            scopes.insert(value, observed);
        }
    }

    fn from_regions(
        candidates: SimplificationCandidates<HEAD_COUNT>,
        regions: HashMap<Atom, Arc<RegionObservation>>,
        scopes: Arc<RwLock<HashMap<Atom, Arc<RegionObservation>>>>,
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
            identity_inventory_complete: candidates.traversal_complete,
            candidates,
            regions,
            scopes,
            polynomial_form: None,
            excludes_internal_connections: false,
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
        for view in Self::regions(expression) {
            let region = if let Some(region) = self.region(view) {
                #[cfg(test)]
                REUSED_REGIONS.with(|count| count.set(count.get() + 1));
                region
            } else {
                #[cfg(test)]
                REPLACEMENT_SCANS.with(|count| count.set(count.get() + 1));
                let replacement = {
                    let cached = self.scopes.read().unwrap();
                    Self::scan_regions(view, Some(&cached))
                };
                let region = replacement
                    .region(view)
                    .expect("the replacement root is observed");
                self.scopes
                    .write()
                    .unwrap()
                    .extend(std::mem::take(&mut *replacement.scopes.write().unwrap()));
                region
            };
            regions.insert(view.to_owned(), region);
        }
        let mut candidates = self.candidates.clone();
        candidates.counts =
            std::array::from_fn(|i| regions.values().map(|region| region.counts[i]).sum());
        candidates.symbols = candidates.counts.map(|count| count > 0);
        candidates.complete = regions
            .values()
            .all(|region| !region.opaque && region.leaf_inventory_complete);
        candidates.intrinsic = regions
            .values()
            .all(|region| region.intrinsic && !region.opaque);
        candidates.traversal_complete = regions
            .values()
            .all(|region| region.leaf_inventory_complete);
        // Index incidence is a graph fact. Conservatively invalidate it after a
        // replacement; the contraction graph decides compatibility and scope.
        candidates.repeated_indices =
            regions.values().map(|region| region.indices).sum::<usize>() > 1;
        candidates.brackets = candidates.symbols[7];
        candidates.dots = regions.values().any(|region| region.dots);
        let mut result = Self::from_regions(
            candidates,
            regions,
            Arc::clone(&self.scopes),
            matches!(expression, AtomView::Add(_)),
        );
        // Each replacement root is either observed now or reuses this domain's
        // exact regional inventory, including certified identity-free regions.
        result.identity_inventory_complete = self.identity_inventory_complete;
        result
    }

    /// No permitted metric/vector endpoint remains in this observed domain.
    pub(crate) fn excludes_contraction_sources(
        &self,
        settings: crate::tensor::ContractSettings<'_>,
    ) -> bool {
        if self.is_terminal_polynomial() {
            return true;
        }
        // Composed terminal outputs retain conservative source inventories even
        // when arbitrary leaf enumeration is unavailable. Require complete
        // regional coverage, rather than inventing a syntactic traversal proof.
        self.identity_inventory_complete
            && !self.regions.is_empty()
            && self
                .regions
                .values()
                .all(|region| region.excludes_contraction_sources(settings))
    }

    pub(crate) fn identity_candidates(&self) -> u8 {
        if self.is_terminal_polynomial() {
            0
        } else {
            RegionObservation::identity_heads(&self.candidates.counts)
        }
    }

    /// Estimate local distribution from already observed factor boundaries.
    /// Reducing an independent family first can avoid carrying a large selected
    /// sum through its kernels. This estimate affects scheduling, never scope.
    pub(crate) fn identity_branch_cost(&self, family: u8) -> usize {
        self.regions
            .iter()
            .fold(1usize, |cost, (expression, region)| {
                if region.identity_candidates() & family != 0
                    && let AtomView::Add(sum) = expression.as_view()
                {
                    cost.saturating_add(sum.get_nargs().saturating_sub(1))
                } else {
                    cost
                }
            })
    }

    /// Invalidate only families whose observed bodies or incident representations
    /// changed. The caller separately checks the arithmetic boundary and callback
    /// proofs. Opaque index labels do not invalidate a completed candidate walk.
    pub(crate) fn changed_identity_families(&self, previous: &Self) -> Option<u8> {
        if !self.identity_inventory_complete
            || !previous.identity_inventory_complete
            || self.regions.is_empty()
            || previous.regions.is_empty()
        {
            return None;
        }
        let changed = self
            .regions
            .iter()
            .filter(|(body, _)| !previous.regions.contains_key(*body))
            .chain(
                previous
                    .regions
                    .iter()
                    .filter(|(body, _)| !self.regions.contains_key(*body)),
            )
            .map(|(_, region)| region)
            .collect::<Vec<_>>();
        let mut dirty = changed.iter().fold(0, |families, region| {
            families | region.identity_candidates()
        });
        let representations = changed
            .iter()
            .flat_map(|region| region.representations.iter().map(|rep| rep.base()))
            .collect::<HashSet<_>>();
        // A changed metric/vector path can expose an identity in an otherwise
        // unchanged body. Representation intersection is conservative: no new
        // graph is needed merely to decide whether to revisit that family.
        for region in self.regions.values().chain(previous.regions.values()) {
            if region
                .representations
                .iter()
                .any(|rep| representations.contains(&rep.base()))
            {
                dirty |= region.identity_candidates();
            }
        }
        Some(dirty)
    }

    pub(crate) fn region(&self, expression: AtomView<'_>) -> Option<Arc<RegionObservation>> {
        self.regions
            .get(expression.get_data())
            .cloned()
            .or_else(|| {
                self.scopes
                    .read()
                    .unwrap()
                    .get(expression.get_data())
                    .cloned()
            })
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
        // Compound slot labels are opaque metadata, not hidden incidences;
        // their conservative candidate flags do not invalidate this proof.
        // Indexed powers can encode extra copies, while scalar powers cannot.
        if !observed.candidates.traversal_complete
            || self.proofs.validated.get() != Some(&true)
            || observed
                .regions
                .values()
                .any(|region| region.counts[7] > 0 || region.has_indexed_powers())
        {
            return;
        }
        let slots = self.structure.logical_slots();
        if observed.max_indices != slots.len()
            || slots
                .iter()
                .any(|slot| !matches!(slot.aind, PartialIndex::Explicit(_)))
            || observed.reserved_indices.iter().any(|&index| {
                let mut index = index;
                while let AbstractIndex::Scoped(scope) = index {
                    index = scope.index();
                }
                matches!(index, AbstractIndex::Open { .. })
            })
            || !self.normalization_is_intrinsic()
        {
            return;
        }
        observed.excludes_internal_connections = true;
        observed.candidates.repeated_indices = false;
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
            // This bounds only memoization. Larger unchanged inputs still run
            // normally; omitting a cache entry never stops reduction.
            const MAX_CACHED_BYTES: usize = 64 * 1024 * 1024;
            if self.bytes.saturating_add(bytes) <= MAX_CACHED_BYTES && self.inputs.insert(key) {
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
    fn leaf_inventories_share_descendants_instead_of_copying_each_ancestor() {
        crate::test_support::test_initialize();
        let leaf = symbolica::symbol!("local_inventory::leaf");
        let left = symbolica::symbol!("local_inventory::left");
        let right = symbolica::symbol!("local_inventory::right");
        for depth in [4, 16, 48] {
            let mut nested = symbolica::function!(leaf, 0);
            for number in 1..=depth {
                nested = (nested + symbolica::function!(leaf, number)).pow(2);
            }
            let source = symbolica::function!(left, &nested) + symbolica::function!(right, &nested);
            let observed = DomainObservations::initial(source.as_view());
            let root = observed.region(source.as_view()).unwrap();
            let mut pending = vec![root.as_ref()];
            let mut visited = HashSet::new();
            let mut retained_leaves = 0;
            while let Some(region) = pending.pop() {
                if !visited.insert(std::ptr::from_ref(region)) {
                    continue;
                }
                retained_leaves += region.leaves.len();
                pending.extend(region.leaf_children.iter().map(AsRef::as_ref));
            }
            // Every function occurrence is stored once, even though each one
            // lies below many retained arithmetic scopes. Flattening child sets
            // at every ancestor made this quadratic in the nesting depth.
            assert_eq!(retained_leaves, 2 * (depth + 1) + 2);

            let mut selected = HashSet::new();
            assert_eq!(
                root.any_leaf(|value| {
                    assert!(
                        selected.insert(value.to_owned()),
                        "predicate repeated for {value}"
                    );
                    false
                }),
                Some(false)
            );
            assert_eq!(selected.len(), depth + 3);
            assert_eq!(
                root.any_leaf(|value| value == symbolica::function!(leaf, depth).as_view()),
                Some(true)
            );
        }
    }

    #[test]
    fn replacement_reuses_tensor_and_symmetric_trace_boundaries() {
        use crate::shorthands::schoonschip::CANDIDATE_NODE_VISITS;
        let reps = crate::test_support::test_initialize();
        let label = spenso::index_symbol!("terminal_inventory::hedge");
        let ports = [0, 1, 2].map(|axis| {
            reps.coad_da
                .to_symbolic([symbolica::function!(label, 81, axis)])
        });
        let structure = crate::color_f!(&ports[0], &ports[1], &ports[2]);
        let trace = spenso::trace_sym!(
            reps.cof_nc.to_symbolic([] as [Atom; 0]);
            ports.iter().map(|port| crate::color_t!(port))
        );
        let source = SymbolicTensor::infer(&structure * &trace).unwrap();
        let initial = source.reduction_observations();
        assert!(initial.scopes.read().unwrap().contains_key(&structure));
        assert!(initial.scopes.read().unwrap().contains_key(&trace));
        let AtomView::Fun(trace_view) = trace.as_view() else {
            panic!("expected trace")
        };
        let symmetric = trace_view.iter().nth(1).unwrap();
        assert!(
            initial
                .scopes
                .read()
                .unwrap()
                .contains_key(symmetric.get_data())
        );

        let changed = Atom::one()
            + Atom::var(symbolica::symbol!("terminal_inventory::weight")) * &source.expression;
        CANDIDATE_NODE_VISITS.with(|count| count.set(0));
        let updated = initial.updated(changed.as_view());
        // One scalar summand, then Mul, coefficient, f and trace. Their tensor
        // arguments, projector words and compound slot payloads remain unopened.
        assert_eq!(CANDIDATE_NODE_VISITS.with(|count| count.get()), 5);
        let fresh = DomainObservations::initial(changed.as_view());
        assert_eq!(updated.regions, fresh.regions);
        assert_eq!(updated.candidates.counts, fresh.candidates.counts);
        assert_eq!(updated.reserved_indices, fresh.reserved_indices);
    }

    #[test]
    fn replacement_reuses_nested_scope_facts_without_losing_occurrences() {
        use crate::shorthands::schoonschip::CANDIDATE_NODE_VISITS;
        crate::test_support::test_initialize();
        let epsilon = *crate::epsilon::EPSILON_SYMBOL;
        let x = Atom::var(symbolica::symbol!("nested_observation::x"));
        let y = Atom::var(symbolica::symbol!("nested_observation::y"));
        let a = symbolica::symbol!("nested_observation::a");
        let b = symbolica::symbol!("nested_observation::b");
        let nested = symbolica::function!(epsilon, &x) + symbolica::function!(epsilon, &y);
        let source = &nested * &x;
        let observed = DomainObservations::initial(source.as_view());
        let changed = symbolica::function!(a, &nested) + symbolica::function!(b, &nested);
        REUSED_SCOPES.with(|count| count.set(0));
        OBSERVED_NODES.with(|count| count.set(0));
        CANDIDATE_NODE_VISITS.with(|count| count.set(0));
        let updated = observed.updated(changed.as_view());
        let visited = OBSERVED_NODES.with(|count| count.get());
        let scanned = CANDIDATE_NODE_VISITS.with(|count| count.get());
        // Two changed wrapper nodes and two reused sum roots. The candidate
        // scanner itself must not revisit either sum's children.
        assert_eq!(scanned, 4);
        assert_eq!(REUSED_SCOPES.with(|count| count.get()), 2);
        assert!(Arc::ptr_eq(&observed.scopes, &updated.scopes));
        OBSERVED_NODES.with(|count| count.set(0));
        CANDIDATE_NODE_VISITS.with(|count| count.set(0));
        let fresh = DomainObservations::initial(changed.as_view());
        assert!(visited < OBSERVED_NODES.with(|count| count.get()));
        assert!(scanned < CANDIDATE_NODE_VISITS.with(|count| count.get()));
        assert_eq!(updated.regions, fresh.regions);
        assert_eq!(updated.candidates.counts, fresh.candidates.counts);
        assert_eq!(updated.candidates.counts[4], 4);
        assert_eq!(updated.epsilon_degree, 1);

        let surviving = updated.updated(symbolica::function!(b, &nested).as_view());
        assert_eq!(surviving.candidates.counts[4], 2);
        assert_ne!(surviving.identity_candidates() & super::super::EPSILON, 0);
        let removed = surviving.updated(x.as_view());
        assert_eq!(removed.candidates.counts[4], 0);
        assert_eq!(removed.identity_candidates(), 0);
    }

    #[test]
    fn replacement_skip_retains_compound_indices_callbacks_and_opaque_metadata() {
        use crate::shorthands::schoonschip::CANDIDATE_NODE_VISITS;
        use spenso::structure::slot::ParseableAind;
        crate::test_support::test_initialize();
        let label = spenso::index_symbol!("cached_observation::hedge");
        let scope = symbolica::symbol!("cached_observation::scope_");
        let indices = std::array::from_fn::<_, 4, _>(|axis| {
            AbstractIndex::Named(label.into(), 37, axis).scoped(scope)
        });
        let ports = indices.map(|index| spenso::mink!(4, index.to_atom()));
        let epsilon = symbolica::function!(
            *crate::epsilon::EPSILON_SYMBOL,
            &ports[0],
            &ports[1],
            &ports[2],
            &ports[3]
        );
        let callback = spenso::tensor_symbol!("cached_observation::callback", norm = |_, _| {});
        let tensor = spenso::tensor_symbol!("cached_observation::T");
        let dimension =
            symbolica::function!(symbolica::symbol!("cached_observation::dimension"), 7);
        let opaque_port = spenso::structure::representation::Minkowski {}
            .to_symbolic([dimension, indices[0].to_atom()]);
        let unchanged = epsilon + symbolica::function!(callback, &ports[0]);
        let first = symbolica::symbol!("cached_observation::first");
        let second = symbolica::symbol!("cached_observation::second");
        let original = symbolica::function!(first, &unchanged);
        let initial = DomainObservations::initial(original.as_view());
        let changed = symbolica::function!(first, &unchanged)
            + symbolica::function!(
                second,
                &unchanged,
                symbolica::function!(tensor, opaque_port)
            );
        CANDIDATE_NODE_VISITS.with(|count| count.set(0));
        let updated = initial.updated(changed.as_view());
        let visits = CANDIDATE_NODE_VISITS.with(|count| count.get());
        CANDIDATE_NODE_VISITS.with(|count| count.set(0));
        let fresh = DomainObservations::initial(changed.as_view());
        assert!(visits < CANDIDATE_NODE_VISITS.with(|count| count.get()));
        assert_eq!(updated.regions, fresh.regions);
        assert_eq!(updated.candidates.counts, fresh.candidates.counts);
        assert_eq!(updated.reserved_indices, fresh.reserved_indices);
        assert_eq!(updated.reserved_indices, HashSet::from(indices));
        assert!(!updated.candidates.complete && !updated.candidates.intrinsic);
        assert_eq!(updated.candidates.counts[4], 2);
        assert!(updated.candidates.repeated_indices);

        let survivor = updated.updated(symbolica::function!(second, &unchanged).as_view());
        assert_eq!(survivor.candidates.counts[4], 1);
        assert!(!survivor.candidates.intrinsic);
        assert!(!survivor.region(unchanged.as_view()).unwrap().intrinsic);
    }

    #[test]
    fn chain_admission_observes_dimensions_variance_and_surviving_regions() {
        crate::test_support::test_initialize();
        let rep = LibraryRep::new_dual("chain_observation::R").unwrap();
        let r = rep.new_rep(3);
        let r4 = rep.new_rep(4);
        let a = SPENSO_TAG.tensor_symbol("chain_observation::A");
        let b = SPENSO_TAG.tensor_symbol("chain_observation::B");
        let first = r.to_symbolic([Atom::num(1)]);
        let compatible = r.dual().to_symbolic([Atom::num(2)]);
        let same_variance = r.to_symbolic([Atom::num(2)]);
        let other_dimension = r4.dual().to_symbolic([Atom::num(2)]);
        let matrix = |head, right: &Atom| symbolica::function!(head, &first, right);
        for right in [&same_variance, &other_dimension] {
            let expression = matrix(a, right);
            assert!(!DomainObservations::initial(expression.as_view()).has_chain_channel(rep));
        }
        let metric = symbolica::function!(ETS.metric, &first, &compatible);
        let first = matrix(a, &compatible);
        let second = matrix(b, &compatible);
        let expression = &first + &second;
        let observed = DomainObservations::initial(expression.as_view());
        let initial = INITIAL_SCANS.with(|count| count.get());
        assert!(observed.has_chain_channel(rep));
        let surviving = observed.updated((&metric + second).as_view());
        assert!(surviving.has_chain_channel(rep));
        let finished = surviving.updated((&metric * Atom::num(2)).as_view());
        assert!(!finished.has_chain_channel(rep));
        assert_eq!(INITIAL_SCANS.with(|count| count.get()), initial);
    }

    #[test]
    fn terminal_polynomial_preserves_open_ports_without_fabricating_scalar_closure() {
        crate::test_support::test_initialize();
        let expression = spenso::g!(spenso::mink!(4, 98928), spenso::mink!(4, 98929));
        let tensor = SymbolicTensor::infer(expression).unwrap();
        let observed = DomainObservations::terminal_polynomial(&tensor.structure, true);
        assert!(observed.is_terminal_polynomial());
        assert!(observed.excludes_internal_connections());
        assert!(!observed.is_closed_scalar());
        assert!(!observed.can_certify_scalar());
        assert_eq!(observed.max_indices, 2);
        assert_eq!(observed.reserved_indices.len(), 2);
        assert_eq!(observed.representations.len(), 1);
        assert_eq!(observed.identity_candidates(), 0);
        assert!(observed.needs_dot_normalization());
        assert!(observed.excludes_contraction_sources(Default::default()));
        assert!(!observed.candidates.complete);
        assert!(!observed.candidates.traversal_complete);
        assert!(observed.region(tensor.expression.as_view()).is_none());

        let epsilon = symbolica::function!(*crate::epsilon::EPSILON_SYMBOL, Atom::num(1));
        let changed = observed.updated(epsilon.as_view());
        assert!(!changed.is_terminal_polynomial());
        assert!(!changed.excludes_internal_connections());
        assert_ne!(changed.identity_candidates() & super::super::EPSILON, 0);
    }

    #[test]
    fn terminal_product_reuses_spectators_without_inventing_leaf_inventory() {
        crate::test_support::test_initialize();
        let ports = [98941, 98942, 98943, 98944].map(|index| spenso::mink!(4, Atom::num(index)));
        assert_eq!(ports.iter().collect::<HashSet<_>>().len(), 4);
        let [a, b, c, d] = ports;
        let polynomial =
            spenso::g!(&a, &b) * spenso::g!(&c, &d) + spenso::g!(&a, &c) * spenso::g!(&b, &d);
        let selected = SymbolicTensor::infer(polynomial).unwrap();
        let terminal = DomainObservations::terminal_polynomial(&selected.structure, false);
        let first =
            SymbolicTensor::infer(crate::color_t!([8, 98945], [3, 98946], [3, 98947])).unwrap();
        let second =
            SymbolicTensor::infer(crate::color_t!([8, 98948], [3, 98949], [3, 98950])).unwrap();
        let spectators = [first, second];
        let first_region = spectators[0]
            .reduction_observations()
            .region(spectators[0].expression.as_view())
            .unwrap();
        spectators[1].reduction_observations();
        let expression =
            &selected.expression * &spectators[0].expression * &spectators[1].expression;
        let combined = SymbolicTensor::infer(expression).unwrap();
        let scans = (
            INITIAL_SCANS.with(|count| count.get()),
            REPLACEMENT_SCANS.with(|count| count.get()),
        );
        let observed = DomainObservations::compose_terminal_product(
            combined.expression.as_view(),
            &combined.structure,
            selected.expression.as_view(),
            &terminal,
            &spectators,
        )
        .unwrap();
        assert_eq!(
            scans,
            (
                INITIAL_SCANS.with(|count| count.get()),
                REPLACEMENT_SCANS.with(|count| count.get())
            )
        );
        assert_eq!(observed.candidates.counts[18], 2);
        assert_eq!(observed.identity_candidates(), super::super::COLOR);
        assert!(observed.excludes_internal_connections());
        assert!(Arc::ptr_eq(
            &first_region,
            &observed.region(spectators[0].expression.as_view()).unwrap()
        ));
        assert!(
            observed
                .region(selected.expression.as_view())
                .unwrap()
                .any_leaf(|_| false)
                .is_none()
        );
        assert!(
            observed
                .region(spectators[0].expression.as_view())
                .unwrap()
                .any_leaf(|_| false)
                .is_some()
        );

        let survivor = &selected.expression * &spectators[1].expression;
        let updated = observed.updated(survivor.as_view());
        assert_eq!(updated.candidates.counts[18], 1);
        assert_eq!(updated.identity_candidates(), super::super::COLOR);
        assert_eq!(
            updated.changed_identity_families(&observed),
            Some(super::super::COLOR)
        );
        assert!(!updated.excludes_internal_connections());
        assert_eq!(
            scans,
            (
                INITIAL_SCANS.with(|count| count.get()),
                REPLACEMENT_SCANS.with(|count| count.get())
            )
        );
    }

    #[test]
    fn terminal_scalar_product_keeps_structural_completion_beside_color_invariants() {
        use crate::tensor::{AlgebraContraction, AlgebraSettings, ReductionStatus};
        let reps = crate::test_support::test_initialize();
        let selected = SymbolicTensor::infer(
            spenso::dot!(spenso::p!(spenso::mink!(4)), spenso::q!(spenso::mink!(4)))
                + spenso::dot!(spenso::p!(spenso::mink!(4)), spenso::p!(spenso::mink!(4))),
        )
        .unwrap();
        let ports = [98961, 98962, 98963].map(|index| reps.coad_da.to_symbolic([Atom::num(index)]));
        let color = SymbolicTensor::infer(spenso::trace_sym!(
            reps.cof_nc.to_symbolic([] as [Atom; 0]);
            ports.iter().chain(&ports).map(|port| crate::color_t!(port))
        ))
        .unwrap();
        let combined = SymbolicTensor::infer(&selected.expression * &color.expression).unwrap();
        let facts = DomainObservations::compose_terminal_product(
            combined.expression.as_view(),
            &combined.structure,
            selected.expression.as_view(),
            &DomainObservations::closed_scalar(false),
            std::slice::from_ref(&color),
        )
        .unwrap();
        assert!(!facts.candidates.complete && !facts.candidates.traversal_complete);
        assert!(
            facts
                .region(selected.expression.as_view())
                .unwrap()
                .any_leaf(|_| false)
                .is_none()
        );
        assert_eq!(facts.identity_candidates(), super::super::COLOR);
        assert!(facts.excludes_contraction_sources(Default::default()));
        // A changed spectator must reuse the scalar certificate without making
        // its unenumerated leaves look like a newly scanned expression.
        let changed_expression =
            &combined.expression * Atom::var(symbolica::symbol!("terminal_color_weight"));
        let changed = facts.updated(changed_expression.as_view());
        assert!(changed.excludes_contraction_sources(Default::default()));
        assert!(
            changed
                .region(selected.expression.as_view())
                .unwrap()
                .any_leaf(|_| false)
                .is_none()
        );

        let cached = SymbolicTensor::from_validated_parts(
            combined.expression.clone(),
            combined.structure.clone(),
        );
        cached.proofs.observations.set(Arc::new(facts)).unwrap();
        let settings = AlgebraSettings {
            contract: AlgebraContraction::Dots,
            ..Default::default()
        };
        for result in [
            cached.contract(Default::default()).unwrap(),
            cached.simplify_algebra(&settings).unwrap(),
        ] {
            assert_eq!(result.expression, combined.expression);
            assert_eq!(result.structure, combined.structure);
            assert_eq!(result.reduction_status(), ReductionStatus::Complete);
            let repeated = result.simplify_algebra(&settings).unwrap();
            assert_eq!(repeated, result);
            assert_eq!(repeated.reduction_status(), ReductionStatus::Complete);
            let admitted = SymbolicTensor::infer(result.expression.clone())
                .unwrap()
                .simplify_algebra(&settings)
                .unwrap();
            assert_eq!(admitted, result);
            assert_eq!(admitted.reduction_status(), ReductionStatus::Complete);
        }
    }

    #[test]
    fn terminal_product_retains_nested_scalar_scopes_without_rescanning() {
        crate::test_support::test_initialize();
        let [a, b, c, d] =
            [99101, 99102, 99103, 99104].map(|index| spenso::mink!(4, Atom::num(index)));
        let first_weight = symbolica::parse_lit!((scope_weight_x + scope_weight_y) ^ 10);
        let second_weight = symbolica::parse_lit!((scope_weight_z + scope_weight_t) ^ 9);
        let spectators = [
            SymbolicTensor::infer(&first_weight * spenso::p!(&a) + spenso::q!(&a)).unwrap(),
            SymbolicTensor::infer(&second_weight * spenso::p!(&c) + spenso::q!(&c)).unwrap(),
        ];
        let selected = SymbolicTensor::infer(
            spenso::g!(&a, &b) * spenso::g!(&c, &d) + spenso::g!(&a, &c) * spenso::g!(&b, &d),
        )
        .unwrap();
        let terminal = DomainObservations::terminal_polynomial(&selected.structure, false);
        let combined = SymbolicTensor::infer(
            &selected.expression * &spectators[0].expression * &spectators[1].expression,
        )
        .unwrap();
        let source_scopes = Arc::clone(&spectators[0].reduction_observations().scopes);
        let weights = [first_weight, second_weight];
        let source_regions = std::array::from_fn::<_, 2, _>(|i| {
            spectators[i]
                .reduction_observations()
                .region(weights[i].as_view())
                .unwrap()
        });
        let scans = (
            INITIAL_SCANS.with(|count| count.get()),
            REPLACEMENT_SCANS.with(|count| count.get()),
        );
        let observed = DomainObservations::compose_terminal_product(
            combined.expression.as_view(),
            &combined.structure,
            selected.expression.as_view(),
            &terminal,
            &spectators,
        )
        .unwrap();
        assert!(Arc::ptr_eq(&source_scopes, &observed.scopes));
        assert!(!observed.excludes_internal_connections());
        for (weight, source) in weights.iter().zip(source_regions) {
            let retained = observed.region(weight.as_view()).unwrap();
            assert!(retained.certifies_scalar_interface_syntax());
            assert!(Arc::ptr_eq(&source, &retained));
        }
        assert_eq!(
            scans,
            (
                INITIAL_SCANS.with(|count| count.get()),
                REPLACEMENT_SCANS.with(|count| count.get())
            )
        );
    }

    #[test]
    fn terminal_product_does_not_hide_new_connections_or_normalized_factors() {
        crate::test_support::test_initialize();
        let ports = [98951, 98952, 98953, 98954].map(|index| spenso::mink!(4, Atom::num(index)));
        assert_eq!(ports.iter().collect::<HashSet<_>>().len(), 4);
        let [a, b, c, d] = ports;
        let polynomial =
            spenso::g!(&a, &b) * spenso::g!(&c, &d) + spenso::g!(&a, &c) * spenso::g!(&b, &d);
        let selected = SymbolicTensor::infer(polynomial).unwrap();
        let terminal = DomainObservations::terminal_polynomial(&selected.structure, false);
        let vector = SymbolicTensor::infer(spenso::p!(&a)).unwrap();
        let combined = SymbolicTensor::infer(&selected.expression * &vector.expression).unwrap();
        let observed = DomainObservations::compose_terminal_product(
            combined.expression.as_view(),
            &combined.structure,
            selected.expression.as_view(),
            &terminal,
            std::slice::from_ref(&vector),
        )
        .unwrap();
        assert!(!observed.excludes_internal_connections());
        assert!(!observed.excludes_contraction_sources(Default::default()));
        assert!(observed.candidates.repeated_indices);
        assert!(
            observed
                .region(selected.expression.as_view())
                .unwrap()
                .contraction_sources(Default::default())
        );
        assert!(
            DomainObservations::compose_terminal_product(
                selected.expression.pow(2).as_view(),
                &selected.structure,
                selected.expression.as_view(),
                &terminal,
                &[],
            )
            .is_none()
        );
    }

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
    fn scalar_identity_powers_do_not_hide_compact_tensor_interfaces() {
        let reps = crate::test_support::test_initialize();
        let casimir = crate::color::CS.cas(Atom::num(2), reps.coad_da.to_symbolic([]));
        let x = Atom::var(symbolica::symbol!("scalar_power_observation::x"));
        for expression in [casimir.pow(2), (&x + &casimir).pow(3)] {
            let tensor = SymbolicTensor::infer(expression.clone()).unwrap();
            let observed = tensor.reduction_observations();
            assert!(
                !observed
                    .region(expression.as_view())
                    .unwrap()
                    .has_indexed_powers()
            );
            assert!(observed.excludes_internal_connections());
            assert_ne!(observed.identity_candidates() & super::super::COLOR, 0);
            assert!(!observed.can_certify_scalar());
        }

        let compact = reps.mink_d.to_symbolic([]);
        let tensor = spenso::tensor_symbol!("scalar_power_observation::T");
        let vector = spenso::p!(&compact);
        // Inspect the written form before intake lowers unresolved powers.
        for expression in [symbolica::function!(tensor, &compact).pow(2), vector.pow(2)] {
            let observed = DomainObservations::initial(expression.as_view());
            assert!(
                observed
                    .region(expression.as_view())
                    .unwrap()
                    .has_indexed_powers()
            );
            assert!(!observed.can_certify_scalar());
        }
        let dot = spenso::g!(spenso::p!(&compact), spenso::q!(&compact)).pow(2);
        let observed = DomainObservations::initial(dot.as_view());
        assert!(!observed.region(dot.as_view()).unwrap().has_indexed_powers());
        assert!(observed.can_certify_scalar());
    }

    #[test]
    fn admitted_index_identity_does_not_authorize_algebra_or_invalidate_observations() {
        use spenso::structure::slot::ParseableAind;
        crate::test_support::test_initialize();
        let named = spenso::index_symbol!("index_observation::hedge");
        let scope = symbolica::symbol!("index_observation::scope_");
        let tensor = spenso::tensor_symbol!("index_observation::T");
        for index in [
            AbstractIndex::Named(named.into(), 13, 1),
            AbstractIndex::Named(named.into(), 13, 1).scoped(scope),
            AbstractIndex::Symbol(AGS.gamma.into()),
            AbstractIndex::Symbol(crate::color::CS.f.into()),
            AbstractIndex::Symbol((*crate::epsilon::EPSILON_SYMBOL).into()),
            AbstractIndex::Symbol(SPENSO_TAG.trace.into()),
        ] {
            let port = spenso::mink!(4, index.to_atom());
            let source = SymbolicTensor::infer(symbolica::function!(tensor, &port)).unwrap();
            let observed = source.reduction_observations();
            assert!(observed.candidates.complete && observed.candidates.intrinsic);
            assert_eq!(observed.identity_candidates(), 0);
            assert_eq!(observed.reserved_indices, HashSet::from([index]));
            assert!(observed.excludes_contraction_sources(Default::default()));
            let result = source
                .simplify_algebra(&super::super::AlgebraSettings::hep())
                .unwrap();
            assert_eq!(result.expression, source.expression);
            assert_eq!(result.structure, source.structure);
            assert_eq!(
                result.reduction_status(),
                super::super::ReductionStatus::Complete
            );
            let again = result
                .simplify_algebra(&super::super::AlgebraSettings::hep())
                .unwrap();
            assert_eq!(again.expression, source.expression);
            assert_eq!(
                again.reduction_status(),
                super::super::ReductionStatus::Complete
            );
        }
    }

    #[test]
    fn external_incidence_proof_accepts_compound_labels_and_scalar_powers() {
        let reps = crate::test_support::test_initialize();
        let label = spenso::index_symbol!("external_incidence::hedge");
        let [a, b, c, d] = std::array::from_fn(|axis| {
            reps.mink_d
                .to_symbolic([symbolica::function!(label, 91, axis)])
        });
        let x = Atom::var(symbolica::symbol!("external_incidence::x"));
        let scalar = (&x + Atom::one()).pow(3);
        let polynomial =
            spenso::g!(&a, &b) * spenso::g!(&c, &d) + spenso::g!(&a, &c) * spenso::g!(&b, &d);
        let tensor = SymbolicTensor::infer(&scalar * &polynomial).unwrap();
        let observed = tensor.reduction_observations();
        assert!(observed.candidates.complete && observed.candidates.intrinsic);
        assert!(observed.excludes_internal_connections());
        assert!(!observed.is_terminal_polynomial());
        assert!(!observed.candidates.repeated_indices);
        assert_eq!(observed.max_indices, 4);
        super::super::super::collection::SHALLOW_BUILDS.with(|count| count.set(0));
        let contracted = tensor.contract(Default::default()).unwrap();
        assert_eq!(contracted.expression, tensor.expression);
        assert_eq!(contracted.structure, tensor.structure);
        assert_eq!(
            contracted.reduction_status(),
            super::super::ReductionStatus::Complete
        );
        super::super::super::collection::SHALLOW_BUILDS.with(|count| assert_eq!(count.get(), 0));

        let changed = tensor
            .with_identity_result((&x + Atom::num(2)).pow(4) * polynomial, None)
            .unwrap();
        assert!(
            changed
                .reduction_observations()
                .excludes_internal_connections()
        );

        // This proof excludes connections, not gamma identities.
        let left = reps.bis4.to_symbolic([Atom::num(98931)]);
        let right = reps.bis4.to_symbolic([Atom::num(98932)]);
        let gamma = SymbolicTensor::infer(crate::gamma!(&left, &right, &a)).unwrap();
        let observed = gamma.reduction_observations();
        assert!(observed.excludes_internal_connections());
        assert_ne!(observed.identity_candidates() & super::super::GAMMA, 0);
    }

    #[test]
    fn external_incidence_proof_requires_admission_and_no_hidden_connections() {
        crate::test_support::test_initialize();
        let a = spenso::mink!(4, 98933);
        let b = spenso::mink!(4, 98934);
        let c = spenso::mink!(4, 98935);
        let tensor = SymbolicTensor::infer(spenso::g!(&b, &c)).unwrap();
        let unchecked = SymbolicTensor::from_normalized_parts(
            tensor.expression.clone(),
            tensor.structure.clone(),
        );
        assert!(
            !unchecked
                .reduction_observations()
                .excludes_internal_connections()
        );

        let internal = SymbolicTensor::infer(
            spenso::g!(&b, &c) + spenso::p!(&a) * spenso::q!(&a) * spenso::g!(&b, &c),
        )
        .unwrap();
        assert!(
            !internal
                .reduction_observations()
                .excludes_internal_connections()
        );
        assert!(
            internal
                .reduction_observations()
                .candidates
                .repeated_indices
        );
        let indexed_power = SymbolicTensor::infer(spenso::p!(&a).pow(2)).unwrap();
        assert!(
            !indexed_power
                .reduction_observations()
                .excludes_internal_connections()
        );
        let unresolved = SymbolicTensor::infer(spenso::p!(spenso::mink!(4))).unwrap();
        assert!(
            !unresolved
                .reduction_observations()
                .excludes_internal_connections()
        );
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
