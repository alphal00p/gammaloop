//! Ordered products with one payload operation per exact named run.

use super::{ExplicitIndexOccurrences, PortRewriteState, ProductPlan, TensorCompositionError};
use crate::tensor::SymbolicTensor;
use spenso::network::tags::SPENSO_TAG;
use spenso::structure::TensorStructure;
use spenso::structure::abstract_index::AbstractIndex;
use spenso::structure::partial::PartialStructure;
use spenso::structure::representation::LibraryRep;
use spenso::structure::slot::{IsAbstractSlot, Slot};
use std::collections::HashMap;
use std::sync::Arc;
use symbolica::atom::{Atom, AtomView};

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
            if Self::interface_is_named(&product.structure)
                && product.borrowed_parts().product_is_normalized()
                && PortRewriteState::exact_arithmetic(product.expression.as_view())
            {
                let mut interface = product.structure.clone();
                let mut sources = vec![&product];
                let mut entries = HashMap::new();
                let history = product.explicit_index_occurrences()?;
                history.accumulate_product(&mut entries)?;
                let mut selected_pair = history.selected_pair;
                while let Some(&factor) = factors.peek() {
                    if !Self::interface_is_named(&factor.structure) {
                        break;
                    }
                    // Admission precedes planning, which precedes combined
                    // history validation, exactly as for binary multiplication.
                    let history = factor.explicit_index_occurrences()?;
                    let plan = ProductPlan::named(&interface, &factor.structure)?;
                    let potential_dot = matches!(&plan, ProductPlan::Contract(_))
                        && interface.canonical().order() == 1
                        && factor.rank() == 1
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
                    let pairs = match &plan {
                        ProductPlan::Scalar | ProductPlan::Outer => &[][..],
                        ProductPlan::Contract(pairs) => pairs,
                        ProductPlan::Compose(..) => unreachable!("named ports choose no channels"),
                    };
                    // Preserve each prefix's canonical merge: Dimension's
                    // structural comparator does not justify one final sort.
                    interface = ProductPlan::free_interface(&interface, &factor.structure, pairs);
                    sources.push(factor);
                    factors.next();
                }
                if sources.len() > 1 {
                    let expression =
                        Atom::mul_many(sources.iter().map(|source| &source.expression));
                    let value = Self::from_validated_parts(expression, interface);
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
