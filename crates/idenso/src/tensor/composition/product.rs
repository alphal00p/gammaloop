//! Ordered products with one payload operation per exact named run.

use super::{
    ExplicitIndexOccurrences, PortPair, PortRewriteState, ProductPlan, TensorCompositionError,
};
use crate::tensor::SymbolicTensor;
use crate::tensor::composition;
use ahash::AHashMap;
use spenso::network::tags::SPENSO_TAG;
use spenso::structure::abstract_index::AbstractIndex;
use spenso::structure::dimension::Dimension;
use spenso::structure::partial::{PartialIndex, PartialSlot, PartialStructure};
use spenso::structure::representation::LibraryRep;
use spenso::structure::slot::{IsAbstractSlot, Slot};
use spenso::structure::{Canonicalized, OrderedStructure, TensorStructure};
#[cfg(test)]
use std::cell::Cell;
use std::collections::hash_map::Entry;
use std::mem;
use std::ops::Range;
use std::sync::Arc;
use symbolica::atom::{Atom, AtomView};

#[cfg(test)]
#[derive(Clone, Copy, Debug, Default)]
pub(super) struct NamedProductWork {
    pub(super) interface_builds: usize,
    pub(super) inventory_checks: usize,
    pub(super) cartesian_checks: usize,
}

#[cfg(test)]
thread_local! {
    pub(super) static NAMED_PRODUCT_WORK: Cell<NamedProductWork> = const {
        Cell::new(NamedProductWork {
            interface_builds: 0,
            inventory_checks: 0,
            cartesian_checks: 0,
        })
    };
}

/// Dimension comparison is consistent within an oriented representation when
/// all dimensions are concrete, or all comparisons of dimensions tie.
#[derive(Clone, Copy, Debug, Default)]
struct DimensionOrder {
    concrete: Option<usize>,
    distinct_concrete: bool,
    symbolic: bool,
}

impl DimensionOrder {
    fn observe(&mut self, dimension: Dimension) {
        match dimension {
            Dimension::Concrete(value) => {
                self.distinct_concrete |= self.concrete.is_some_and(|previous| previous != value);
                self.concrete = self.concrete.or(Some(value));
            }
            Dimension::Symbolic(_) => self.symbolic = true,
        }
    }

    fn combined(self, other: Self) -> Self {
        Self {
            concrete: self.concrete.or(other.concrete),
            distinct_concrete: self.distinct_concrete
                || other.distinct_concrete
                || matches!((self.concrete, other.concrete), (Some(left), Some(right)) if left != right),
            symbolic: self.symbolic || other.symbolic,
        }
    }

    fn is_consistent(self) -> bool {
        !self.symbolic || !self.distinct_concrete
    }
}

#[derive(Debug)]
struct IncomingPort {
    slot: PartialSlot,
    logical_position: usize,
    contracted: bool,
}

#[derive(Debug, Default)]
struct NamedFactor {
    ports: Vec<IncomingPort>,
    dimensions: Vec<(LibraryRep, DimensionOrder)>,
}

impl NamedFactor {
    /// Classification, source correspondence and the ordering certificate share
    /// one borrowed canonical traversal. No cross-factor pairing happens here.
    fn load(&mut self, interface: &PartialStructure) -> bool {
        self.ports.clear();
        self.dimensions.clear();
        self.ports.reserve(interface.canonical().order());
        for (position, slot) in interface.canonical().external_structure_iter().enumerate() {
            if !SymbolicTensor::index_is_named(slot.aind) {
                self.ports.clear();
                self.dimensions.clear();
                return false;
            }
            let representation = slot.rep();
            if let Some((_, order)) = self
                .dimensions
                .iter_mut()
                .find(|(rep, _)| *rep == representation.rep)
            {
                order.observe(representation.dim);
            } else {
                let mut order = DimensionOrder::default();
                order.observe(representation.dim);
                self.dimensions.push((representation.rep, order));
            }
            self.ports.push(IncomingPort {
                slot,
                logical_position: interface.layout().canonical_to_logical_position(position),
                contracted: false,
            });
        }
        true
    }
}

#[derive(Debug)]
struct NamedSlot {
    slot: PartialSlot,
    factor: usize,
    canonical_position: usize,
    live: bool,
}

#[derive(Debug)]
struct CanonicalPrefix {
    interface: PartialStructure,
    slots: Vec<usize>,
}

/// Most live keys retain at most two source occurrences. Larger transient sets
/// still have ordinary owned storage; this container is not a multiplicity proof.
#[derive(Debug, Default)]
struct LiveOccurrences {
    inline: [usize; 2],
    len: usize,
    spill: Vec<usize>,
}

impl LiveOccurrences {
    fn push(&mut self, id: usize) {
        if self.len < self.inline.len() {
            self.inline[self.len] = id;
        } else {
            self.spill.push(id);
        }
        self.len += 1;
    }

    fn iter(&self) -> impl Iterator<Item = &usize> {
        self.inline[..self.len.min(self.inline.len())]
            .iter()
            .chain(self.spill.iter())
    }

    fn remove(&mut self, id: usize) {
        let mut retained = 0;
        for position in 0..self.len.min(self.inline.len()) {
            let candidate = self.inline[position];
            if candidate != id {
                self.inline[retained] = candidate;
                retained += 1;
            }
        }
        self.spill.retain(|&candidate| candidate != id);
        while retained < self.inline.len() && !self.spill.is_empty() {
            self.inline[retained] = self.spill.remove(0);
            retained += 1;
        }
        self.len = retained + self.spill.len();
    }

    fn is_empty(&self) -> bool {
        self.len == 0
    }
}

/// The inventory owns cross-factor pairing, never occurrence history. Stable
/// record IDs retain source occurrences even when their indices coincide.
#[derive(Debug)]
struct PreparedNamedProduct {
    inventory: AHashMap<Slot<LibraryRep, AbstractIndex>, LiveOccurrences>,
    slots: Vec<NamedSlot>,
    dimensions: AHashMap<LibraryRep, DimensionOrder>,
    legacy: Option<CanonicalPrefix>,
    rank: usize,
    factors: usize,
}

impl ProductPlan {
    fn named_key(slot: PartialSlot) -> Slot<LibraryRep, AbstractIndex> {
        let PartialIndex::Explicit(index) = slot.aind else {
            unreachable!("named preparation admits only explicit named ports");
        };
        slot.rep().base().slot(index)
    }

    fn prepare_named(interface: &PartialStructure) -> Option<PreparedNamedProduct> {
        let mut factor = NamedFactor::default();
        if !factor.load(interface) {
            return None;
        }
        let mut result = PreparedNamedProduct {
            inventory: AHashMap::new(),
            slots: Vec::with_capacity(factor.ports.len()),
            dimensions: AHashMap::new(),
            legacy: None,
            rank: 0,
            factors: 0,
        };
        let _ = result.add_ports(&mut factor, &[]);
        let certified = result.include_dimensions(&factor);
        if !certified {
            // The seed is an established prefix. Do not sort it an extra time.
            result.legacy = Some(CanonicalPrefix {
                interface: Canonicalized::identity(interface.canonical().clone()),
                slots: (0..result.slots.len()).collect(),
            });
        }
        Some(result)
    }
}

impl PreparedNamedProduct {
    fn pairings(
        &mut self,
        factor: &NamedFactor,
        pairs: &mut Vec<PortPair>,
    ) -> Result<(), TensorCompositionError> {
        pairs.clear();
        // All incoming occurrences see the same old inventory. In particular,
        // duplicate partners cannot be greedily consumed before ambiguity checks.
        for (right, port) in factor.ports.iter().enumerate() {
            if let Some(candidates) = self.inventory.get(&ProductPlan::named_key(port.slot)) {
                for &left in candidates.iter() {
                    #[cfg(test)]
                    NAMED_PRODUCT_WORK.with(|work| {
                        let mut observed = work.get();
                        observed.inventory_checks += 1;
                        work.set(observed);
                    });
                    if composition::automatically_contractible(&self.slots[left].slot, &port.slot) {
                        pairs.push(PortPair { left, right });
                    }
                }
            }
        }
        if ProductPlan::named_pairs_are_ambiguous(pairs) {
            let prefix = self.legacy.take().unwrap_or_else(|| {
                let order = self.source_order();
                self.canonicalize(order)
            });
            let mut positions = vec![0; self.slots.len()];
            for (position, &id) in prefix.slots.iter().enumerate() {
                positions[id] = position;
            }
            for pair in pairs.iter_mut() {
                pair.left = positions[pair.left];
                pair.right = factor.ports[pair.right].logical_position;
            }
            pairs.sort_unstable_by_key(|pair| (pair.left, pair.right));
            return Err(TensorCompositionError::Ambiguous(mem::take(pairs)));
        }
        Ok(())
    }

    fn include_dimensions(&mut self, factor: &NamedFactor) -> bool {
        let mut certified = true;
        for &(rep, incoming) in &factor.dimensions {
            let order = self.dimensions.entry(rep).or_default();
            *order = order.combined(incoming);
            certified &= order.is_consistent();
        }
        certified
    }

    fn add_ports(&mut self, factor: &mut NamedFactor, pairs: &[PortPair]) -> Range<usize> {
        let start = self.slots.len();
        for pair in pairs {
            factor.ports[pair.right].contracted = true;
            let previous = &mut self.slots[pair.left];
            previous.live = false;
            let key = ProductPlan::named_key(previous.slot);
            match self.inventory.entry(key) {
                Entry::Occupied(mut entry) => {
                    entry.get_mut().remove(pair.left);
                    if entry.get().is_empty() {
                        let _ = entry.remove();
                    }
                }
                Entry::Vacant(_) => unreachable!("a selected port is live"),
            }
        }
        for (position, port) in factor.ports.iter().enumerate() {
            if port.contracted {
                continue;
            }
            let id = self.slots.len();
            self.slots.push(NamedSlot {
                slot: port.slot,
                factor: self.factors,
                canonical_position: position,
                live: true,
            });
            self.inventory
                .entry(ProductPlan::named_key(port.slot))
                .or_default()
                .push(id);
        }
        self.rank = self.rank + factor.ports.len() - 2 * pairs.len();
        self.factors += 1;
        start..self.slots.len()
    }

    fn source_order(&self) -> Vec<usize> {
        let mut order = self
            .inventory
            .values()
            .flat_map(|occurrences| occurrences.iter())
            .copied()
            .collect::<Vec<_>>();
        order.sort_unstable_by_key(|&id| {
            let slot = &self.slots[id];
            (slot.factor, slot.canonical_position)
        });
        order
    }

    fn canonicalize(&self, order: Vec<usize>) -> CanonicalPrefix {
        #[cfg(test)]
        NAMED_PRODUCT_WORK.with(|work| {
            let mut observed = work.get();
            observed.interface_builds += 1;
            work.set(observed);
        });
        // The existing two-stage usize sorter also supplies correspondence for
        // diagnostics and subsequent legacy prefixes, without matching again.
        let interface =
            OrderedStructure::new(order.iter().map(|&id| self.slots[id].slot).collect());
        let slots = (0..order.len())
            .map(|position| order[interface.layout().canonical_to_logical_position(position)])
            .collect();
        CanonicalPrefix {
            interface: Canonicalized::identity(interface.into_canonical()),
            slots,
        }
    }

    /// Commit only after source admission, pairing, eligibility and history have
    /// succeeded. Certificate failure changes ordering work, not dispatch mode.
    fn commit(&mut self, factor: &mut NamedFactor, pairs: &[PortPair]) {
        if !self.include_dimensions(factor) && self.legacy.is_none() {
            let order = self.source_order();
            self.legacy = Some(self.canonicalize(order));
        }
        let added = self.add_ports(factor, pairs);
        if let Some(prefix) = self.legacy.take() {
            let order = prefix
                .slots
                .into_iter()
                .filter(|&id| self.slots[id].live)
                .chain(added)
                .collect();
            self.legacy = Some(self.canonicalize(order));
        }
    }

    fn finish(mut self) -> PartialStructure {
        if let Some(prefix) = self.legacy.take() {
            return prefix.interface;
        }
        #[cfg(test)]
        NAMED_PRODUCT_WORK.with(|work| {
            let mut observed = work.get();
            observed.interface_builds += 1;
            work.set(observed);
        });
        let order = self.source_order();
        Canonicalized::identity(OrderedStructure::merge_slots(
            order.into_iter().map(|id| self.slots[id].slot),
        ))
    }
}

#[derive(Debug)]
struct ProductOccurrence {
    slot: Slot<LibraryRep, AbstractIndex>,
    total: usize,
    base: usize,
    dual: usize,
}

impl ExplicitIndexOccurrences {
    /// Scratch lookup for a run, using exact dimension equality rather than
    /// structural ordering. Keep the first representative's orientation.
    fn accumulate_product(
        &self,
        entries: &mut AHashMap<Slot<LibraryRep, AbstractIndex>, ProductOccurrence>,
    ) -> Result<(), TensorCompositionError> {
        self.validate()?;
        let mut overflow = None;
        for (slot, total, base, dual) in self.oriented_entries() {
            let key = slot.rep().base().slot(slot.aind());
            let entry = entries.entry(key).or_insert(ProductOccurrence {
                slot,
                total: 0,
                base: 0,
                dual: 0,
            });
            entry.total = entry.total.saturating_add(total);
            entry.base = entry.base.saturating_add(base);
            entry.dual = entry.dual.saturating_add(dual);
            if entry.total > 2
                && overflow.is_none_or(|previous| Self::compare(&slot, &previous).is_lt())
            {
                // Binary history merging reports the first collision in its
                // total key order, within the earliest failing prefix.
                overflow = Some(slot);
            }
        }
        if let Some(slot) = overflow {
            return Err(TensorCompositionError::InvalidExplicitMultiplicity {
                index: slot.aind(),
                occurrences: 3,
            });
        }
        Ok(())
    }

    fn finish_product(
        entries: AHashMap<Slot<LibraryRep, AbstractIndex>, ProductOccurrence>,
        selected_pair: bool,
    ) -> Self {
        let mut entries = entries.into_values().collect::<Vec<_>>();
        entries.sort_unstable_by(|left, right| Self::compare(&left.slot, &right.slot));
        let mut result = Self {
            selected_pair,
            ..Self::default()
        };
        for entry in entries {
            result.indices.push(entry.slot);
            result.dummies.push(entry.total == 2);
            result
                .orientations
                .push((entry.base as u8) | ((entry.dual as u8) << 2));
        }
        result
    }
}

impl SymbolicTensor<PartialStructure> {
    /// Multiply borrowed factors in input order, starting from scalar one.
    ///
    /// Iteration stops at its first `None`, including for non-fused sources.
    /// Empty input returns one. A singleton still passes through scalar
    /// multiplication, including canonical free-axis ordering and nominal zero
    /// ports. Pairing decisions and collision checks occur at each prefix, just
    /// as in a left fold of [`Self::multiply`]. Exact named runs share one
    /// factor-preserving payload multiplication; positional ports, compact dots,
    /// rounded arithmetic and callbacks retain sequential operation dispatch.
    /// Named pairing advances over live source occurrences. Canonical ordering
    /// is deferred only when dimension comparisons are consistent across the run.
    ///
    /// # Errors
    ///
    /// Returns the first prefix's error from [`Self::multiply`], including invalid
    /// source metadata, ambiguous pairings, unsupported channel composition or
    /// excessive explicit-index multiplicity. Later factors are not admitted
    /// after an earlier prefix fails.
    ///
    /// ```
    /// use idenso::tensor::SymbolicTensor;
    /// use symbolica::atom::Atom;
    ///
    /// let factors = [SymbolicTensor::infer(Atom::num(2))?,
    ///                SymbolicTensor::infer(Atom::num(3))?];
    /// let product = SymbolicTensor::multiply_many(&factors)?;
    /// assert_eq!(product.expression(), &Atom::num(6));
    /// # Ok::<_, Box<dyn std::error::Error>>(())
    /// ```
    pub fn multiply_many<'a>(
        factors: impl IntoIterator<Item = &'a Self>,
    ) -> Result<Self, TensorCompositionError> {
        let mut factors = factors.into_iter().fuse().peekable();
        let mut product = Self::infer(Atom::one())?;
        let mut incoming = NamedFactor::default();
        let mut pairs = Vec::new();
        while factors.peek().is_some() {
            if let Some(mut prepared) = ProductPlan::prepare_named(&product.structure)
                && product.borrowed_parts().product_is_normalized()
                && PortRewriteState::exact_arithmetic(product.expression.as_view())
            {
                let mut sources = vec![&product];
                let mut entries = AHashMap::new();
                let history = product.explicit_index_occurrences()?;
                history.accumulate_product(&mut entries)?;
                let mut selected_pair = history.selected_pair;
                while let Some(&factor) = factors.peek() {
                    if !incoming.load(&factor.structure) {
                        break;
                    }
                    // Admission precedes planning, which precedes combined
                    // history validation, exactly as for binary multiplication.
                    let history = factor.explicit_index_occurrences()?;
                    prepared.pairings(&incoming, &mut pairs)?;
                    let potential_dot = !pairs.is_empty()
                        && prepared.rank == 1
                        && incoming.ports.len() == 1
                        && matches!(factor.expression.as_view(), AtomView::Fun(vector)
                            if vector.get_symbol().has_tag(&SPENSO_TAG.rank1));
                    if potential_dot
                        || !factor.borrowed_parts().product_is_normalized()
                        || !PortRewriteState::exact_arithmetic(factor.expression.as_view())
                    {
                        break;
                    }
                    history.accumulate_product(&mut entries)?;
                    selected_pair |= history.selected_pair;
                    prepared.commit(&mut incoming, &pairs);
                    sources.push(factor);
                    factors.next();
                }
                if sources.len() > 1 {
                    let expression =
                        Atom::mul_many(sources.iter().map(|source| &source.expression));
                    let value = Self::from_validated_parts(expression, prepared.finish());
                    let history = ExplicitIndexOccurrences::finish_product(entries, selected_pair);
                    let _ = value.proofs.occurrences.set(Arc::new(history));
                    let _ = value.proofs.normalized_product.set(true);
                    product = value.validate_rewrite("multiply_many", &sources)?;
                }
            }
            // Flush before dispatching an operation whose spelling, rounding or
            // callbacks depend on the actual accumulated payload.
            if let Some(factor) = factors.next() {
                product = product.multiply(factor)?;
            }
        }
        Ok(product)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::representations::ColorFundamental;
    use crate::test_support;
    use spenso::structure::partial::PartialStructureExt;
    use spenso::structure::representation::{ExtendibleReps, RepName, Representation};
    use symbolica::symbol;

    fn rep() -> Representation<LibraryRep> {
        ExtendibleReps::EUCLIDEAN.new_rep(4).cast()
    }

    #[test]
    fn live_occurrences_preserve_order_across_inline_spill_removal_and_refill() {
        let mut occurrences = LiveOccurrences::default();
        assert!(occurrences.is_empty());
        occurrences.push(10);
        assert_eq!(occurrences.iter().copied().collect::<Vec<_>>(), [10]);
        occurrences.push(20);
        assert_eq!(occurrences.iter().copied().collect::<Vec<_>>(), [10, 20]);
        assert_eq!(occurrences.spill.capacity(), 0);
        occurrences.push(30);
        occurrences.push(40);
        assert_eq!(
            occurrences.iter().copied().collect::<Vec<_>>(),
            [10, 20, 30, 40]
        );
        occurrences.remove(20);
        assert_eq!(
            occurrences.iter().copied().collect::<Vec<_>>(),
            [10, 30, 40]
        );
        occurrences.push(50);
        occurrences.remove(40);
        assert_eq!(
            occurrences.iter().copied().collect::<Vec<_>>(),
            [10, 30, 50]
        );
        occurrences.remove(10);
        assert_eq!(occurrences.iter().copied().collect::<Vec<_>>(), [30, 50]);
        assert!(occurrences.spill.is_empty());
        for absent in [99, 10, 20, 40] {
            occurrences.remove(absent);
            assert_eq!(occurrences.iter().copied().collect::<Vec<_>>(), [30, 50]);
        }
        occurrences.remove(30);
        assert_eq!(occurrences.iter().copied().collect::<Vec<_>>(), [50]);
        occurrences.remove(50);
        assert!(occurrences.is_empty());
        assert_eq!(occurrences.iter().count(), 0);
        occurrences.push(60);
        occurrences.push(70);
        assert_eq!(occurrences.iter().copied().collect::<Vec<_>>(), [60, 70]);
        assert!(occurrences.spill.is_empty());
    }

    #[test]
    fn prepared_inventory_deletes_and_reinserts_an_empty_key_without_stale_pairings() {
        test_support::test_initialize();
        let slot = rep().slot(PartialIndex::Explicit(AbstractIndex::Normal(99901)));
        let interface = PartialStructure::from_logical_slots([slot]);
        let key = ProductPlan::named_key(slot);
        let mut prepared = ProductPlan::prepare_named(&interface).unwrap();
        let mut incoming = NamedFactor::default();
        let mut pairs = Vec::new();
        assert!(incoming.load(&interface));
        prepared.pairings(&incoming, &mut pairs).unwrap();
        assert_eq!(pairs, [PortPair { left: 0, right: 0 }]);
        prepared.commit(&mut incoming, &pairs);
        assert_eq!(prepared.rank, 0);
        assert!(!prepared.inventory.contains_key(&key));
        assert!(incoming.ports[0].contracted);

        // Inventory storage can refill independently of the public operation's
        // separate history check. This does not authorize a third index use.
        assert!(incoming.load(&interface));
        assert!(!incoming.ports[0].contracted);
        prepared.pairings(&incoming, &mut pairs).unwrap();
        assert!(pairs.is_empty());
        prepared.commit(&mut incoming, &pairs);
        assert_eq!(prepared.rank, 1);
        let live = prepared.inventory[&key].iter().copied().collect::<Vec<_>>();
        assert_eq!(live, [1]);
        assert!(prepared.slots[live[0]].live);
        assert_eq!(prepared.finish().logical_slots(), [slot]);
    }

    #[test]
    fn named_factor_scratch_resets_size_dimensions_marks_and_logical_correspondence() {
        test_support::test_initialize();
        let base: Representation<LibraryRep> = ColorFundamental {}.new_rep(3).cast();
        let named = |representation: Representation<LibraryRep>, index| {
            representation.slot(PartialIndex::Explicit(AbstractIndex::Normal(index)))
        };
        let symbolic = Representation {
            rep: rep().rep,
            dim: Dimension::Symbolic(symbol!("storage_scratch_dimension").into()),
        };
        let large = PartialStructure::from_logical_slots([
            named(base, 9),
            named(rep().rep.new_rep(2), 7),
            named(symbolic, 5),
            named(rep(), 3),
        ]);
        let small = PartialStructure::from_logical_slots([named(rep(), 1)]);
        let scalar = PartialStructure::from_logical_slots([]);
        let rejected = PartialStructure::from_logical_slots([
            named(rep(), 11),
            base.slot(PartialIndex::open(0)),
        ]);
        let mut scratch = NamedFactor::default();
        assert!(scratch.load(&small));
        assert!(scratch.load(&large));
        assert_eq!(scratch.ports.len(), 4);
        assert_eq!(scratch.dimensions.len(), 2);
        let dimension = scratch
            .dimensions
            .iter()
            .find(|(key, _)| *key == rep().rep)
            .unwrap()
            .1;
        assert!(dimension.symbolic && dimension.distinct_concrete);
        let logical = large.logical_slots();
        for port in &mut scratch.ports {
            assert_eq!(logical[port.logical_position], port.slot);
            port.contracted = true;
        }
        let port_capacity = scratch.ports.capacity();
        let dimension_capacity = scratch.dimensions.capacity();
        assert!(scratch.load(&small));
        assert_eq!(scratch.ports.len(), 1);
        assert!(!scratch.ports[0].contracted);
        assert_eq!(scratch.ports[0].logical_position, 0);
        assert_eq!(scratch.dimensions.len(), 1);
        let dimension = scratch.dimensions[0].1;
        assert_eq!(dimension.concrete, Some(4));
        assert!(!dimension.symbolic && !dimension.distinct_concrete);
        assert_eq!(scratch.ports.capacity(), port_capacity);
        assert_eq!(scratch.dimensions.capacity(), dimension_capacity);
        assert!(scratch.load(&scalar));
        assert!(scratch.ports.is_empty() && scratch.dimensions.is_empty());
        assert!(!scratch.load(&rejected));
        assert!(scratch.ports.is_empty() && scratch.dimensions.is_empty());

        let dual = PartialStructure::from_logical_slots([named(base.dual().rep.new_rep(5), 13)]);
        assert!(scratch.load(&dual));
        assert_eq!(scratch.ports.len(), 1);
        assert!(!scratch.ports[0].contracted);
        assert_eq!(scratch.dimensions.len(), 1);
        assert_eq!(scratch.dimensions[0].0, base.dual().rep);
        let dimension = scratch.dimensions[0].1;
        assert_eq!(dimension.concrete, Some(5));
        assert!(!dimension.symbolic && !dimension.distinct_concrete);
    }

    #[test]
    fn accumulated_history_keeps_first_orientation_bounds_scopes_and_selected_pair() {
        test_support::test_initialize();
        let base: Representation<LibraryRep> = ColorFundamental {}.new_rep(3).cast();
        let index = AbstractIndex::Normal(99911);
        let scoped = index.scoped(symbol!("storage_history_scope"));
        let symbolic = Representation {
            rep: base.rep,
            dim: Dimension::Symbolic(symbol!("storage_history_dimension").into()),
        };
        let mut dual = ExplicitIndexOccurrences::default();
        dual.add(base.dual().slot(index), 1);
        dual.add(symbolic.slot(index), 1);
        dual.selected_pair = true;
        let mut fundamental = ExplicitIndexOccurrences::default();
        fundamental.add(base.slot(index), 1);
        fundamental.add(base.slot(scoped), 1);
        for (first, second, representative) in [
            (&dual, &fundamental, base.dual().slot(index)),
            (&fundamental, &dual, base.slot(index)),
        ] {
            let mut expected = first.clone();
            expected.append(second);
            expected.validate().unwrap();
            let mut entries = AHashMap::new();
            first.accumulate_product(&mut entries).unwrap();
            second.accumulate_product(&mut entries).unwrap();
            let actual = ExplicitIndexOccurrences::finish_product(
                entries,
                first.selected_pair || second.selected_pair,
            );
            assert_eq!(
                actual.oriented_entries().collect::<Vec<_>>(),
                expected.oriented_entries().collect::<Vec<_>>()
            );
            assert_eq!(actual.entries().count(), 3);
            let shared = actual
                .oriented_entries()
                .find(|(slot, _, _, _)| {
                    slot.aind() == index && slot.rep().dim == Dimension::Concrete(3)
                })
                .unwrap();
            assert_eq!(shared, (representative, 2, 1, 1));
            assert!(actual.selected_pair);
            assert_eq!(actual.count(symbolic.slot(index)), 1);
            assert_eq!(actual.count(base.slot(scoped)), 1);
            let mut entries = AHashMap::new();
            actual.accumulate_product(&mut entries).unwrap();
            let mut extra = ExplicitIndexOccurrences::default();
            extra.add(base.slot(index), 1);
            assert!(matches!(extra.accumulate_product(&mut entries),
                Err(TensorCompositionError::InvalidExplicitMultiplicity { index: failed, occurrences: 3 })
                    if failed == index
            ));
        }
    }
}
