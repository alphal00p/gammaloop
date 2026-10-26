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
use std::collections::HashMap;
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
}

#[derive(Debug)]
struct NamedFactor {
    ports: Vec<IncomingPort>,
    dimensions: Vec<(LibraryRep, DimensionOrder)>,
}

impl NamedFactor {
    /// Classification, source correspondence and the ordering certificate share
    /// one borrowed canonical traversal. No cross-factor pairing happens here.
    fn from_structure(interface: &PartialStructure) -> Option<Self> {
        let mut result = Self {
            ports: Vec::with_capacity(interface.canonical().order()),
            dimensions: Vec::new(),
        };
        for (position, slot) in interface.canonical().external_structure_iter().enumerate() {
            if !SymbolicTensor::index_is_named(slot.aind) {
                return None;
            }
            let representation = slot.rep();
            if let Some((_, order)) = result
                .dimensions
                .iter_mut()
                .find(|(rep, _)| *rep == representation.rep)
            {
                order.observe(representation.dim);
            } else {
                let mut order = DimensionOrder::default();
                order.observe(representation.dim);
                result.dimensions.push((representation.rep, order));
            }
            result.ports.push(IncomingPort {
                slot,
                logical_position: interface.layout().canonical_to_logical_position(position),
            });
        }
        Some(result)
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

/// The inventory owns cross-factor pairing, never occurrence history. Stable
/// record IDs retain source occurrences even when their indices coincide.
#[derive(Debug)]
struct PreparedNamedProduct {
    inventory: AHashMap<Slot<LibraryRep, AbstractIndex>, Vec<usize>>,
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
        let factor = NamedFactor::from_structure(interface)?;
        let mut result = PreparedNamedProduct {
            inventory: AHashMap::new(),
            slots: Vec::with_capacity(factor.ports.len()),
            dimensions: AHashMap::new(),
            legacy: None,
            rank: 0,
            factors: 0,
        };
        let _ = result.add_ports(&factor, &[]);
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
    fn pairings(&mut self, factor: &NamedFactor) -> Result<Vec<PortPair>, TensorCompositionError> {
        let mut pairs = Vec::new();
        // All incoming occurrences see the same old inventory. In particular,
        // duplicate partners cannot be greedily consumed before ambiguity checks.
        for (right, port) in factor.ports.iter().enumerate() {
            if let Some(candidates) = self.inventory.get(&ProductPlan::named_key(port.slot)) {
                for &left in candidates {
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
        if ProductPlan::named_pairs_are_ambiguous(&pairs) {
            let prefix = self.legacy.take().unwrap_or_else(|| {
                let order = self.source_order();
                self.canonicalize(order)
            });
            let mut positions = vec![0; self.slots.len()];
            for (position, &id) in prefix.slots.iter().enumerate() {
                positions[id] = position;
            }
            for pair in &mut pairs {
                pair.left = positions[pair.left];
                pair.right = factor.ports[pair.right].logical_position;
            }
            pairs.sort_unstable_by_key(|pair| (pair.left, pair.right));
            return Err(TensorCompositionError::Ambiguous(pairs));
        }
        Ok(pairs)
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

    fn add_ports(&mut self, factor: &NamedFactor, pairs: &[PortPair]) -> Range<usize> {
        let start = self.slots.len();
        let mut contracted = vec![false; factor.ports.len()];
        for pair in pairs {
            contracted[pair.right] = true;
            let previous = &mut self.slots[pair.left];
            previous.live = false;
            let key = ProductPlan::named_key(previous.slot);
            let occurrences = self
                .inventory
                .get_mut(&key)
                .expect("a selected port is live");
            occurrences.retain(|&id| id != pair.left);
            if occurrences.is_empty() {
                self.inventory.remove(&key);
            }
        }
        for (position, port) in factor.ports.iter().enumerate() {
            if contracted[position] {
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
            .flatten()
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
    fn commit(&mut self, factor: &NamedFactor, pairs: &[PortPair]) {
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
        entries: &mut HashMap<Slot<LibraryRep, AbstractIndex>, ProductOccurrence>,
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
        entries: HashMap<Slot<LibraryRep, AbstractIndex>, ProductOccurrence>,
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
        while factors.peek().is_some() {
            if let Some(mut prepared) = ProductPlan::prepare_named(&product.structure)
                && product.borrowed_parts().product_is_normalized()
                && PortRewriteState::exact_arithmetic(product.expression.as_view())
            {
                let mut sources = vec![&product];
                let mut entries = HashMap::new();
                let history = product.explicit_index_occurrences()?;
                history.accumulate_product(&mut entries)?;
                let mut selected_pair = history.selected_pair;
                while let Some(&factor) = factors.peek() {
                    let Some(incoming) = NamedFactor::from_structure(&factor.structure) else {
                        break;
                    };
                    // Admission precedes planning, which precedes combined
                    // history validation, exactly as for binary multiplication.
                    let history = factor.explicit_index_occurrences()?;
                    let pairs = prepared.pairings(&incoming)?;
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
                    prepared.commit(&incoming, &pairs);
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
