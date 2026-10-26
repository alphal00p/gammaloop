//! Sorted metadata products, with expression assembly independent of index history.

use super::{ExplicitIndexOccurrences, PortPair, ProductPlan, TensorCompositionError};
use crate::tensor::SymbolicTensor;
use crate::tensor::composition;
use spenso::network::tags::SPENSO_TAG;
use spenso::structure::abstract_index::AbstractIndex;
use spenso::structure::partial::{PartialIndex, PartialSlot, PartialStructure};
use spenso::structure::representation::{LibraryRep, RepName};
use spenso::structure::slot::{IsAbstractSlot, Slot};
use spenso::structure::{Canonicalized, OrderedStructure, SlotIndex};
use std::cmp::Reverse;
use std::collections::BinaryHeap;
use std::iter;
use std::sync::Arc;
use symbolica::atom::{Atom, AtomView};

/// Splitting canonical ports by orientation makes every borrowed stream monotone
/// in exact base representation, dimension and index, including the dual stream.
#[derive(Clone, Copy, Debug)]
enum NamedStream<'a> {
    Ports {
        slots: &'a [PartialSlot],
        canonical_offset: usize,
    },
    History(&'a ExplicitIndexOccurrences),
}

impl NamedStream<'_> {
    fn key(self, position: usize) -> Option<Slot<LibraryRep, AbstractIndex>> {
        match self {
            Self::Ports { slots, .. } => slots.get(position).copied().map(ProductPlan::named_key),
            Self::History(history) => history
                .indices
                .get(position)
                .map(|slot| slot.rep().base().slot(slot.aind())),
        }
    }
}

/// One cursor per nonempty stream; no per-port inventory survives a key.
#[derive(Clone, Copy, Debug, Eq, Ord, PartialEq, PartialOrd)]
struct NamedCursor {
    key: Slot<LibraryRep, AbstractIndex>,
    factor: usize,
    stream: usize,
    position: usize,
}

#[derive(Clone, Copy, Debug)]
struct SourcePort {
    slot: PartialSlot,
    factor: usize,
    /// Logical source ports flattened in factor order, used only in diagnostics.
    position: usize,
}

#[derive(Clone, Copy, Debug, Default)]
struct NamedStep {
    contractions: usize,
    /// An actual rank-one contraction needs this port, not a prefix interface.
    left_slot: Option<PartialSlot>,
}

#[derive(Debug, Default)]
struct ProductOccurrence {
    representative: Option<Slot<LibraryRep, AbstractIndex>>,
    total: usize,
    base: usize,
    dual: usize,
}

impl ProductOccurrence {
    fn include(&mut self, history: &ExplicitIndexOccurrences, position: usize) {
        self.representative.get_or_insert(history.indices[position]);
        self.total = self
            .total
            .saturating_add(1 + usize::from(history.dummies[SlotIndex(position)]));
        let orientations = history.orientations[position];
        // Alternative branches can carry both bounds without carrying their sum.
        self.base = self.base.saturating_add(usize::from(orientations & 3));
        self.dual = self
            .dual
            .saturating_add(usize::from((orientations >> 2) & 3));
    }
}

#[derive(Debug)]
struct NamedProduct {
    interface: PartialStructure,
    occurrences: ExplicitIndexOccurrences,
    steps: Vec<NamedStep>,
}

impl ProductPlan {
    fn named_key(slot: PartialSlot) -> Slot<LibraryRep, AbstractIndex> {
        let PartialIndex::Explicit(index) = slot.aind else {
            unreachable!("a named stream contains only explicit indices");
        };
        slot.rep().base().slot(index)
    }

    /// Find a cross-factor ambiguity witness in linear work, even for a large
    /// invalid key. Matching policy remains the same as binary named planning.
    fn named_group_ambiguity(ports: &[SourcePort]) -> Option<Vec<PortPair>> {
        if ports.len() <= 2 {
            return None;
        }
        let mut totals = [0usize; 2];
        let mut representatives = [None; 2];
        for port in ports {
            let orientation = usize::from(!port.slot.rep().rep.is_base());
            totals[orientation] += 1;
            representatives[orientation].get_or_insert(port.slot);
        }
        let mut start = 0;
        while start < ports.len() {
            let factor = ports[start].factor;
            let end = start
                + ports[start..]
                    .iter()
                    .take_while(|port| port.factor == factor)
                    .count();
            let mut local = [0usize; 2];
            for port in &ports[start..end] {
                local[usize::from(!port.slot.rep().rep.is_base())] += 1;
            }
            for port in &ports[start..end] {
                let partners = representatives
                    .iter()
                    .enumerate()
                    .filter(|(_, representative)| {
                        representative.is_some_and(|slot| {
                            composition::automatically_contractible(&port.slot, &slot)
                        })
                    })
                    .map(|(orientation, _)| totals[orientation] - local[orientation])
                    .sum::<usize>();
                if partners > 1 {
                    return Some(
                        ports
                            .iter()
                            .filter(|other| {
                                other.factor != factor
                                    && composition::automatically_contractible(
                                        &port.slot,
                                        &other.slot,
                                    )
                            })
                            .take(2)
                            .map(|other| PortPair {
                                left: port.position,
                                right: other.position,
                            })
                            .collect(),
                    );
                }
            }
            start = end;
        }
        None
    }

    /// A valid admitted exact key has at most two public occurrences. Repeated
    /// ports of one factor remain distinct; only different factors can contract.
    fn named_group_survivors(
        ports: &[SourcePort],
        steps: &mut [NamedStep],
        survivors: &mut Vec<PartialSlot>,
        dual_survivors: &mut Vec<PartialSlot>,
    ) {
        if let [left, right] = ports
            && left.factor != right.factor
            && composition::automatically_contractible(&left.slot, &right.slot)
        {
            let step = &mut steps[right.factor];
            step.contractions += 1;
            step.left_slot = Some(left.slot);
            return;
        }
        for port in ports {
            if port.slot.rep().rep.is_base() {
                survivors.push(port.slot);
            } else {
                dual_survivors.push(port.slot);
            }
        }
    }

    fn merge_named(
        sources: &[&SymbolicTensor<PartialStructure>],
    ) -> Result<NamedProduct, TensorCompositionError> {
        let mut streams = Vec::with_capacity(4 * sources.len());
        let mut cursors = Vec::with_capacity(4 * sources.len());
        let mut offsets = Vec::with_capacity(sources.len());
        let mut source_ports = 0;
        let mut occurrences = ExplicitIndexOccurrences::default();
        for (factor, source) in sources.iter().enumerate() {
            offsets.push(source_ports);
            source_ports += source.rank();
            let canonical = source.structure.canonical();
            let mut canonical_offset = 0;
            for slots in [
                canonical.self_dual_slice(),
                canonical.base_slice(),
                canonical.dual_slice(),
            ] {
                let stream = NamedStream::Ports {
                    slots,
                    canonical_offset,
                };
                if let Some(key) = stream.key(0) {
                    cursors.push(Reverse(NamedCursor {
                        key,
                        factor,
                        stream: streams.len(),
                        position: 0,
                    }));
                    streams.push(stream);
                }
                canonical_offset += slots.len();
            }
            let history = source.explicit_index_occurrences()?;
            occurrences.selected_pair |= history.selected_pair;
            let stream = NamedStream::History(history);
            if let Some(key) = stream.key(0) {
                cursors.push(Reverse(NamedCursor {
                    key,
                    factor,
                    stream: streams.len(),
                    position: 0,
                }));
                streams.push(stream);
            }
        }

        let mut heap = BinaryHeap::from(cursors);
        let mut ports = Vec::with_capacity(2);
        let mut survivors = Vec::with_capacity(source_ports);
        let mut dual_survivors = Vec::new();
        let mut steps = vec![NamedStep::default(); sources.len()];
        while let Some(Reverse(first)) = heap.peek().copied() {
            ports.clear();
            let mut occurrence = ProductOccurrence::default();
            while heap.peek().is_some_and(|cursor| cursor.0.key == first.key) {
                let Reverse(mut cursor) = heap.pop().expect("the current key has a cursor");
                match streams[cursor.stream] {
                    NamedStream::Ports {
                        slots,
                        canonical_offset,
                    } => {
                        let source = sources[cursor.factor];
                        ports.push(SourcePort {
                            slot: slots[cursor.position],
                            factor: cursor.factor,
                            position: offsets[cursor.factor]
                                + source.structure.layout().canonical_to_logical_position(
                                    canonical_offset + cursor.position,
                                ),
                        });
                    }
                    NamedStream::History(history) => occurrence.include(history, cursor.position),
                }
                cursor.position += 1;
                if let Some(key) = streams[cursor.stream].key(cursor.position) {
                    cursor.key = key;
                    heap.push(Reverse(cursor));
                }
            }
            if let Some(pairs) = Self::named_group_ambiguity(&ports) {
                return Err(TensorCompositionError::Ambiguous(pairs));
            }
            if occurrence.total > 2 {
                return Err(TensorCompositionError::InvalidExplicitMultiplicity {
                    index: first.key.aind(),
                    occurrences: 3,
                });
            }
            Self::named_group_survivors(&ports, &mut steps, &mut survivors, &mut dual_survivors);
            occurrences.push_product_occurrence(occurrence);
        }
        let interface = Canonicalized::identity(OrderedStructure::from_sorted_slots(
            survivors.into_iter().chain(dual_survivors),
        ));
        Ok(NamedProduct {
            interface,
            occurrences,
            steps,
        })
    }
}

impl ExplicitIndexOccurrences {
    /// Groups arrive in total key order, with the first factor's representative.
    /// Hidden pairs and cancelled indices need no public port to reach this output.
    fn push_product_occurrence(&mut self, occurrence: ProductOccurrence) {
        self.indices.push(
            occurrence
                .representative
                .expect("admitted public indices have retained occurrence metadata"),
        );
        self.dummies.push(occurrence.total == 2);
        self.orientations
            .push((occurrence.base as u8) | ((occurrence.dual as u8) << 2));
    }
}

impl SymbolicTensor<PartialStructure> {
    /// Flush arithmetic only where an actual rank-one contraction can choose
    /// the operation owner's compact dot spelling.
    fn named_product_expression(
        sources: &[&Self],
        steps: &[NamedStep],
    ) -> Result<Atom, TensorCompositionError> {
        let mut expression = Atom::one();
        let mut pending = 0;
        let mut rank = 0;
        for (factor, (source, step)) in sources.iter().zip(steps).enumerate() {
            if rank == 1
                && source.rank() == 1
                && step.contractions == 1
                && matches!(source.expression.as_view(), AtomView::Fun(vector)
                    if vector.get_symbol().has_tag(&SPENSO_TAG.rank1))
            {
                expression = Atom::mul_many(
                    iter::once(&expression).chain(
                        sources[pending..factor]
                            .iter()
                            .map(|source| &source.expression),
                    ),
                );
                pending = factor;
                if matches!(expression.as_view(), AtomView::Fun(vector)
                    if vector.get_symbol().has_tag(&SPENSO_TAG.rank1))
                {
                    expression = Self::compact_dot_parts(
                        &expression,
                        step.left_slot
                            .expect("a selected contraction has a left port"),
                        source,
                    )?;
                    pending += 1;
                }
            }
            rank = rank + source.rank() - 2 * step.contractions;
        }
        Ok(Atom::mul_many(iter::once(&expression).chain(
            sources[pending..].iter().map(|source| &source.expression),
        )))
    }

    /// Multiply borrowed tensors, using one sorted metadata merge per named run.
    ///
    /// Inputs own complete public interfaces and scoped explicit-index histories.
    /// Named runs merge those borrowed streams without inspecting payloads for
    /// admission, reconstructing prefixes or sorting output. Duplicates in one
    /// factor remain separate ports; contractions join different factors. Open
    /// ports and channels use the existing multiplication operation. Expression
    /// assembly uses Symbolica arithmetic and the existing compact-dot owner.
    ///
    /// Empty input returns scalar one. A singleton receives canonical named-axis
    /// order, including the nominal ports of zero. Iteration stops at its first
    /// `None`, including for non-fused sources.
    ///
    /// # Errors
    ///
    /// Returns multiplication errors for newly ambiguous pairings, excessive
    /// explicit multiplicity or unsupported channel composition. Named errors are
    /// selected in global exact-key order. Ambiguity witnesses identify logical
    /// source ports flattened in factor order; no earliest-prefix error is promised.
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
        let factors = factors.into_iter().fuse().collect::<Vec<_>>();
        let mut product = Self::from_validated_parts(
            Atom::one(),
            Canonicalized::identity(OrderedStructure::default()),
            Arc::new(ExplicitIndexOccurrences::default()),
        );
        let mut position = 0;
        while position < factors.len() {
            if Self::interface_is_named(&product.structure)
                && Self::interface_is_named(&factors[position].structure)
            {
                let start = position;
                position += 1;
                while position < factors.len()
                    && Self::interface_is_named(&factors[position].structure)
                {
                    position += 1;
                }
                let sources = iter::once(&product)
                    .chain(factors[start..position].iter().copied())
                    .collect::<Vec<_>>();
                let NamedProduct {
                    interface,
                    occurrences,
                    steps,
                } = ProductPlan::merge_named(&sources)?;
                let expression = Self::named_product_expression(&sources, &steps)?;
                product = Self::from_validated_parts(expression, interface, Arc::new(occurrences));
            } else {
                product = product.multiply(factors[position])?;
                position += 1;
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
    use spenso::structure::dimension::Dimension;
    use spenso::structure::partial::PartialStructureExt;
    use spenso::structure::representation::{ExtendibleReps, Representation};
    use symbolica::symbol;

    fn rep() -> Representation<LibraryRep> {
        ExtendibleReps::EUCLIDEAN.new_rep(4).cast()
    }

    fn admitted_zero(
        slots: impl IntoIterator<Item = PartialSlot>,
    ) -> SymbolicTensor<PartialStructure> {
        let interface = PartialStructure::from_logical_slots(slots);
        let occurrences = ExplicitIndexOccurrences::from_interface(&interface);
        occurrences.validate().unwrap();
        SymbolicTensor::from_validated_parts(Atom::Zero, interface, Arc::new(occurrences))
    }

    #[test]
    fn merge_keeps_same_factor_duplicates_and_only_contracts_different_factors() {
        test_support::test_initialize();
        let base: Representation<LibraryRep> = ColorFundamental {}.new_rep(3).cast();
        let [i, j, k] = [99901, 99903, 99905].map(AbstractIndex::Normal);
        let slot = |representation: Representation<LibraryRep>, index| {
            representation.slot(PartialIndex::Explicit(index))
        };
        let left = admitted_zero([
            slot(rep(), i),
            slot(rep(), i),
            slot(base, j),
            slot(base.dual(), j),
            slot(base.dual(), k),
        ]);
        let right = admitted_zero([slot(base, k)]);
        let merged = ProductPlan::merge_named(&[&left, &right]).unwrap();
        assert_eq!(merged.steps[0].contractions, 0);
        assert_eq!(merged.steps[1].contractions, 1);
        assert_eq!(merged.steps[1].left_slot, Some(slot(base.dual(), k)));
        assert_eq!(
            merged.interface.logical_slots(),
            [
                slot(rep(), i),
                slot(rep(), i),
                slot(base, j),
                slot(base.dual(), j)
            ]
        );
        for (representation, index) in [(rep(), i), (base, j), (base, k)] {
            assert_eq!(merged.occurrences.count(representation.slot(index)), 2);
        }
    }

    #[test]
    fn merge_keeps_equal_orientation_ports_from_different_factors() {
        test_support::test_initialize();
        let base: Representation<LibraryRep> = ColorFundamental {}.new_rep(3).cast();
        let index = AbstractIndex::Normal(99907);
        for representation in [base, base.dual()] {
            let slot = representation.slot(PartialIndex::Explicit(index));
            let left = admitted_zero([slot]);
            let right = admitted_zero([slot]);
            let merged = ProductPlan::merge_named(&[&left, &right]).unwrap();
            assert_eq!(merged.interface.logical_slots(), [slot, slot]);
            assert_eq!(merged.steps[1].contractions, 0);
            assert_eq!(merged.occurrences.count(base.slot(index)), 2);
        }
    }

    #[test]
    fn duplicate_partners_report_a_cross_factor_ambiguity_witness() {
        test_support::test_initialize();
        let slot = rep().slot(PartialIndex::Explicit(AbstractIndex::Normal(99909)));
        let repeated = admitted_zero([slot, slot]);
        let single = admitted_zero([slot]);
        for sources in [[&repeated, &single], [&single, &repeated]] {
            let Err(TensorCompositionError::Ambiguous(pairs)) = ProductPlan::merge_named(&sources)
            else {
                panic!("two cross-factor partners must be ambiguous");
            };
            assert_eq!(pairs.len(), 2);
            assert_eq!(pairs[0].left, pairs[1].left);
            assert_ne!(pairs[0].right, pairs[1].right);
            let boundary = sources[0].rank();
            assert!(pairs.iter().all(|pair| {
                (pair.left < boundary && pair.right >= boundary)
                    || (pair.right < boundary && pair.left >= boundary)
            }));
        }
    }

    #[test]
    fn history_preserves_first_orientation_independent_bounds_scopes_and_selected_pair() {
        test_support::test_initialize();
        let base: Representation<LibraryRep> = ColorFundamental {}.new_rep(3).cast();
        let index = AbstractIndex::Normal(99911);
        let scoped = index.scoped(symbol!("sorted_product_history_scope"));
        let symbolic = Representation {
            rep: base.rep,
            dim: Dimension::Symbolic(symbol!("sorted_product_history_dimension").into()),
        };
        let mut alternatives = ExplicitIndexOccurrences::default();
        alternatives.add_oriented(base.dual().slot(index), 1, 1, 1);
        alternatives.add(symbolic.slot(index), 1);
        alternatives.selected_pair = true;
        let mut fundamental = ExplicitIndexOccurrences::default();
        fundamental.add(base.slot(index), 1);
        fundamental.add(base.slot(scoped), 1);
        let left = SymbolicTensor::from_validated_parts(
            Atom::Zero,
            Canonicalized::identity(OrderedStructure::default()),
            Arc::new(alternatives),
        );
        let right = SymbolicTensor::from_validated_parts(
            Atom::Zero,
            Canonicalized::identity(OrderedStructure::default()),
            Arc::new(fundamental),
        );
        for (sources, representative) in [
            ([&left, &right], base.dual().slot(index)),
            ([&right, &left], base.slot(index)),
        ] {
            let merged = ProductPlan::merge_named(&sources).unwrap();
            assert!(merged.interface.logical_slots().is_empty());
            assert!(merged.occurrences.selected_pair);
            assert_eq!(merged.occurrences.entries().count(), 3);
            assert_eq!(merged.occurrences.count(symbolic.slot(index)), 1);
            assert_eq!(merged.occurrences.count(base.slot(scoped)), 1);
            assert!(
                merged
                    .occurrences
                    .oriented_entries()
                    .any(|entry| { entry == (representative, 2, 2, 1) })
            );
        }
    }

    #[test]
    fn hidden_history_collisions_follow_global_key_order() {
        test_support::test_initialize();
        let [low, high] = [99913, 99915].map(AbstractIndex::Normal);
        let hidden = |index, count| {
            let mut history = ExplicitIndexOccurrences::default();
            history.add(rep().slot(index), count);
            SymbolicTensor::from_validated_parts(
                Atom::Zero,
                Canonicalized::identity(OrderedStructure::default()),
                Arc::new(history),
            )
        };
        let first = hidden(high, 2);
        let second = hidden(high, 1);
        let third = hidden(low, 2);
        let fourth = hidden(low, 1);
        assert!(matches!(
            ProductPlan::merge_named(&[&first, &second, &third, &fourth]),
            Err(TensorCompositionError::InvalidExplicitMultiplicity { index, occurrences: 3 })
                if index == low
        ));
    }

    #[test]
    fn sorted_output_distinguishes_dimensions_and_appends_dual_survivors() {
        test_support::test_initialize();
        let base: Representation<LibraryRep> = ColorFundamental {}.new_rep(3).cast();
        let index = AbstractIndex::Normal(99917);
        let dimensions = [
            Dimension::Concrete(2),
            Dimension::Concrete(5),
            Dimension::Symbolic(symbol!("sorted_product_dim_a").into()),
            Dimension::Symbolic(symbol!("sorted_product_dim_b").into()),
        ];
        let slot = |representation, dim| {
            Representation {
                rep: representation,
                dim,
            }
            .slot(PartialIndex::Explicit(index))
        };
        let paired =
            dimensions.map(|dim| admitted_zero([slot(base.rep, dim), slot(base.dual().rep, dim)]));
        let independent = dimensions
            .into_iter()
            .rev()
            .flat_map(|dim| {
                [rep().rep, LibraryRep::Dummy]
                    .map(|representation| admitted_zero([slot(representation, dim)]))
            })
            .collect::<Vec<_>>();
        let borrowed = independent.iter().chain(paired.iter()).collect::<Vec<_>>();
        let merged = ProductPlan::merge_named(&borrowed).unwrap();
        let mut expected = dimensions
            .into_iter()
            .flat_map(|dim| {
                [rep().rep, base.rep, base.dual().rep, LibraryRep::Dummy]
                    .map(|representation| slot(representation, dim))
            })
            .collect::<Vec<_>>();
        expected.sort_unstable();
        assert_eq!(merged.interface.logical_slots(), expected);
        assert_eq!(merged.occurrences.entries().count(), 12);
        let history_keys = merged
            .occurrences
            .indices
            .iter()
            .map(|slot| slot.rep().base().slot::<AbstractIndex, _>(slot.aind()))
            .collect::<Vec<_>>();
        assert!(history_keys.windows(2).all(|pair| pair[0] < pair[1]));
    }
}
