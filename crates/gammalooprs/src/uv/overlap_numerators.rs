//! Offline, finite component contraction of already bound numerator coefficients.
//! This module is included only by the test-only overlap certificate. Neither
//! scalar nor tensor-valued graph numerators are polynomial-expanded here.

use std::{collections::HashMap, ops::Deref};

use color_eyre::{Result, eyre::ensure};
use spenso::{
    iterators::IteratableTensor,
    network::{ExecutionResult, MinResultRank, SequentialRef},
    structure::{TensorStructure, slot::IsAbstractSlot},
};
use symbolica::{coefficient::Coefficient, prelude::*};

use crate::{
    cff::expression::OrientationID,
    numerator::symbolica_ext::NumeratorAtomExt,
    utils::{FUN_LIB, GS, TENSORLIB},
};

const MAX_BODY_BYTES: usize = 2 * 1024 * 1024;
const MAX_COMPONENTS: usize = 4096;
const MAX_RESULT_BYTES: usize = 8 * 1024 * 1024;
const MAX_CACHE_BYTES: usize = 64 * 1024 * 1024;

#[derive(Debug, thiserror::Error)]
pub(crate) enum NumeratorJetFailure {
    #[error("unproven numerator coefficient: unsupported {0}")]
    Unsupported(&'static str),
    #[error("unproven numerator coefficient: budget exceeded ({0})")]
    Budget(&'static str),
}

/// A content key is an opaque numerator owner, including for a nonconstant
/// rank-zero tensor. Only a literal exact numeric factor leaves that owner.
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct CanonicalTensor {
    pub(crate) key: Atom,
    pub(crate) scalar_prefactor: Atom,
    pub(crate) is_zero: bool,
}

impl CanonicalTensor {
    fn number(value: Atom) -> Self {
        Self {
            key: Atom::num(1),
            is_zero: value.is_zero(),
            scalar_prefactor: value,
        }
    }
}

/// Cache by the complete physically bound coefficient, not by a family name.
/// The caller owns binding, routing, and proving that this coefficient has no
/// remaining Laurent variable. Unequal keys do not establish inequivalence.
pub(crate) struct NumeratorJetCache {
    values: HashMap<Atom, CanonicalTensor>,
    bytes: usize,
    max_bytes: usize,
}

impl Default for NumeratorJetCache {
    fn default() -> Self {
        Self {
            values: HashMap::new(),
            bytes: 0,
            max_bytes: MAX_CACHE_BYTES,
        }
    }
}

impl NumeratorJetCache {
    pub(crate) fn canonicalize(&mut self, body: &Atom) -> Result<CanonicalTensor> {
        if let Some(value) = self.values.get(body) {
            return Ok(value.clone());
        }
        ensure!(
            body.as_view().get_byte_size() <= MAX_BODY_BYTES,
            NumeratorJetFailure::Budget("input bytes")
        );
        for head in [
            Symbol::IF,
            OrientationID::symbol(),
            GS.theta,
            GS.orientation_delta,
            symbol!("gammalooprs::uv::numerator_family"),
        ] {
            ensure!(
                !body.contains_symbol(head),
                NumeratorJetFailure::Unsupported("guard or unresolved numerator family")
            );
        }
        let mut exact = true;
        body.visitor(&mut |view| {
            if let AtomView::Num(number) = view {
                exact &= matches!(number.get_coeff_view().to_owned(), Coefficient::Complex(_));
            }
            exact
        });
        ensure!(
            exact,
            NumeratorJetFailure::Unsupported("inexact or undefined number")
        );

        let (prefactor, body_without_number) = Self::split_number(body.clone());
        let mut result = if body_without_number.is_one() {
            CanonicalTensor::number(Atom::num(1))
        } else {
            Self::contract(&body_without_number)?
        };
        result.scalar_prefactor *= prefactor;
        result.is_zero |= result.scalar_prefactor.is_zero();
        let bytes = body.as_view().get_byte_size()
            + result.key.as_view().get_byte_size()
            + result.scalar_prefactor.as_view().get_byte_size();
        // Memoization cannot turn a completed exact contraction into an
        // unproven coefficient. Evict old results when full; a result larger
        // than the cache remains usable without caching. The individual input,
        // component-count and output-byte bounds above still apply.
        if bytes <= self.max_bytes {
            if self.bytes.saturating_add(bytes) > self.max_bytes {
                self.values.clear();
                self.bytes = 0;
            }
            self.bytes += bytes;
            self.values.insert(body.clone(), result.clone());
        }
        Ok(result)
    }

    fn split_number(atom: Atom) -> (Atom, Atom) {
        match atom.as_view() {
            AtomView::Num(_) => (atom, Atom::num(1)),
            AtomView::Mul(product) => {
                let mut number = Atom::num(1);
                let mut factors = Vec::new();
                for factor in product.iter() {
                    if matches!(factor, AtomView::Num(_)) {
                        number *= factor;
                    } else {
                        factors.push(factor.to_owned());
                    }
                }
                (number, Atom::mul_many(&factors))
            }
            _ => (Atom::num(1), atom),
        }
    }

    fn contract(body: &Atom) -> Result<CanonicalTensor> {
        let mut network = body.parse_into_net()?;
        let slots = network.graph.dangling_indices();
        let mut capacity = 1usize;
        for slot in slots {
            let dimension = usize::try_from(slot.dim())?;
            capacity = capacity.saturating_mul(dimension);
            ensure!(
                capacity <= MAX_COMPONENTS,
                NumeratorJetFailure::Budget("external tensor components")
            );
        }
        // Only finite tensor execution: this closes sum boundaries without
        // distributing scalar spectators, powers, or the surrounding graph.
        network.graph.contract_ready_sum_boundaries();
        network.execute::<SequentialRef, MinResultRank, _, _, _>(
            TENSORLIB.read().unwrap().deref(),
            FUN_LIB.deref(),
        )?;
        let tensor = match network.result_tensor(TENSORLIB.read().unwrap().deref())? {
            ExecutionResult::One => return Ok(CanonicalTensor::number(Atom::num(1))),
            ExecutionResult::Zero => return Ok(CanonicalTensor::number(Atom::Zero)),
            ExecutionResult::Val(tensor) => tensor.into_owned(),
        };
        ensure!(
            tensor.size()? <= MAX_COMPONENTS,
            NumeratorJetFailure::Budget("result tensor components")
        );
        // Canonicalize the external coordinate order together with each data
        // index. Keep the free labels and their variance; only contracted dummy
        // names have disappeared. Ignore storage names and sparse/dense layout.
        let original_slots = tensor
            .external_structure_iter()
            .map(|slot| slot.to_atom())
            .collect::<Vec<_>>();
        let mut order = (0..original_slots.len()).collect::<Vec<_>>();
        order.sort_by(|left, right| original_slots[*left].cmp(&original_slots[*right]));
        let slots = FunctionBuilder::new(symbol!("gammalooprs::uv::overlap_slots"))
            .add_args(order.iter().map(|index| &original_slots[*index]))
            .finish();
        let mut components = Vec::new();
        let mut bytes = slots.as_view().get_byte_size();
        let mut scalar_prefactor = Atom::num(1);
        for (index, value) in tensor.iter_flat() {
            let mut value = value.collect_compact_factors();
            if value.is_zero() {
                continue;
            }
            if original_slots.is_empty() {
                let (number, residual) = Self::split_number(value);
                scalar_prefactor = number;
                value = residual;
                if value.is_one() {
                    return Ok(CanonicalTensor::number(scalar_prefactor));
                }
            }
            let coordinates = tensor.expanded_index(index)?;
            let coordinate = FunctionBuilder::new(symbol!("gammalooprs::uv::overlap_component"))
                .add_args(
                    order
                        .iter()
                        .map(|position| Atom::num(coordinates.indices[*position])),
                )
                .finish();
            bytes = bytes.saturating_add(value.as_view().get_byte_size());
            ensure!(
                bytes <= MAX_RESULT_BYTES,
                NumeratorJetFailure::Budget("contracted component bytes")
            );
            components.push((coordinate, value));
        }
        if components.is_empty() {
            return Ok(CanonicalTensor::number(Atom::Zero));
        }
        components.sort_by(|left, right| left.0.cmp(&right.0));
        let entries = components.into_iter().map(|(index, value)| {
            function!(symbol!("gammalooprs::uv::overlap_entry"), index, value)
        });
        let key = FunctionBuilder::new(symbol!("gammalooprs::uv::overlap_numerator_tensor"))
            .add_arg(&slots)
            .add_args(entries)
            .finish();
        ensure!(
            key.as_view().get_byte_size() <= MAX_RESULT_BYTES,
            NumeratorJetFailure::Budget("content key bytes")
        );
        Ok(CanonicalTensor {
            key,
            scalar_prefactor,
            is_zero: false,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::initialisation::test_initialise;
    use linnet::half_edge::involution::EdgeIndex;

    #[test]
    fn open_tensor_keys_ignore_contracted_dummy_names_and_factor_order() {
        test_initialise().unwrap();
        let mut cache = NumeratorJetCache::default();
        // Q has the registered rank-one tensor owner; an arbitrary p(index)
        // would be an opaque scalar even if its argument resembles a slot.
        let a = parse_lit!((x + y) * spenso::g(spenso::mink(4, 10), spenso::mink(4, 20)))
            * GS.emr_mom(EdgeIndex(1), parse_lit!(spenso::mink(4, 20)));
        let b = GS.emr_mom(EdgeIndex(1), parse_lit!(spenso::mink(4, 30)))
            * parse_lit!(spenso::g(spenso::mink(4, 10), spenso::mink(4, 30)) * (y + x));
        assert_eq!(
            cache.canonicalize(&a).unwrap(),
            cache.canonicalize(&b).unwrap()
        );
        let other = parse_lit!(x + y) * GS.emr_mom(EdgeIndex(2), parse_lit!(spenso::mink(4, 10)));
        assert_ne!(
            cache.canonicalize(&a).unwrap().key,
            cache.canonicalize(&other).unwrap().key
        );
    }

    #[test]
    fn closed_numerator_is_opaque_but_literal_numbers_and_zero_are_not() {
        test_initialise().unwrap();
        let mut cache = NumeratorJetCache::default();
        let body = parse_lit!((x + y) * (z + w));
        let value = cache.canonicalize(&body).unwrap();
        assert!(!value.key.is_one());
        assert!(value.scalar_prefactor.is_one());
        let doubled = cache.canonicalize(&(Atom::num(2) * body)).unwrap();
        assert_eq!(value.key, doubled.key);
        assert_eq!(doubled.scalar_prefactor, Atom::num(2));
        let trace = parse_lit!(spenso::g(spenso::bis(4, 9), spenso::bis(4, 9)));
        let constant = cache.canonicalize(&trace).unwrap();
        assert!(constant.key.is_one());
        assert_eq!(constant.scalar_prefactor, Atom::num(4));
        assert!(cache.canonicalize(&Atom::Zero).unwrap().is_zero);
    }

    #[test]
    fn bounded_cache_eviction_and_recomputation_preserve_exact_keys() {
        test_initialise().unwrap();
        let bodies = [1, 2].map(|edge| {
            parse_lit!(x + y) * GS.emr_mom(EdgeIndex(edge), parse_lit!(spenso::mink(4, 10)))
        });
        let mut reference = NumeratorJetCache::default();
        let expected = bodies
            .each_ref()
            .map(|body| reference.canonicalize(body).unwrap());
        assert_ne!(expected[0].key, expected[1].key);
        let max_bytes = bodies
            .iter()
            .zip(&expected)
            .map(|(body, result)| {
                body.as_view().get_byte_size()
                    + result.key.as_view().get_byte_size()
                    + result.scalar_prefactor.as_view().get_byte_size()
            })
            .max()
            .unwrap();
        let mut bounded = NumeratorJetCache {
            max_bytes,
            ..Default::default()
        };
        for index in [0, 1, 0, 1] {
            assert_eq!(
                bounded.canonicalize(&bodies[index]).unwrap(),
                expected[index]
            );
            assert_eq!(bounded.values.len(), 1);
            assert!(bounded.values.contains_key(&bodies[index]));
            assert!(!bounded.values.contains_key(&bodies[1 - index]));
            assert!(bounded.bytes <= max_bytes);
        }
        let mut uncached = NumeratorJetCache {
            max_bytes: 0,
            ..Default::default()
        };
        for (body, expected) in bodies.iter().zip(expected) {
            assert_eq!(uncached.canonicalize(body).unwrap(), expected);
            assert!(uncached.values.is_empty());
            assert_eq!(uncached.bytes, 0);
        }
    }

    #[test]
    fn unresolved_families_and_large_tensors_remain_unproven() {
        test_initialise().unwrap();
        let mut cache = NumeratorJetCache::default();
        let unresolved = function!(symbol!("gammalooprs::uv::numerator_family"), 0, 1);
        assert!(cache.canonicalize(&unresolved).is_err());
        let wide = parse_lit!(
            spenso::g(spenso::mink(4, 1), spenso::mink(4, 2))
                * spenso::g(spenso::mink(4, 3), spenso::mink(4, 4))
                * spenso::g(spenso::mink(4, 5), spenso::mink(4, 6))
                * spenso::g(spenso::mink(4, 7), spenso::mink(4, 8))
        );
        let error = cache.canonicalize(&wide).unwrap_err();
        assert!(matches!(
            error.downcast_ref::<NumeratorJetFailure>(),
            Some(NumeratorJetFailure::Budget("external tensor components"))
        ));
    }
}
