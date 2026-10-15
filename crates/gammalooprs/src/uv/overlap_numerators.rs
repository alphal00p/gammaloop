//! Offline, finite component contraction of already bound numerator coefficients.
//! Shared by the overlap and on-shell test certificates. Neither
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
            // Numeric signs and weights belong to the shared compact owner.
            // In particular, do not run an unprotected numeric expansion here:
            // it would enter fractional powers and denominator bodies.
            let mut value = value.collect_compact_factors();
            if !value.is_zero() {
                // One bounded collection in repeated literal scalar factors
                // can close a finite component identity without opening those
                // factors. This is a single deterministic decomposition, not
                // a search over alternative groupings or polynomial variables.
                let mut occurrences = HashMap::<AtomView<'_>, usize>::new();
                value.visitor(&mut |part| {
                    match part {
                        AtomView::Fun(_) => {
                            *occurrences.entry(part).or_default() += 1;
                            return false;
                        }
                        AtomView::Pow(power) => {
                            if i64::try_from(power.get_exp()).is_ok_and(|power| power > 0)
                                && matches!(power.get_base(), AtomView::Add(_) | AtomView::Fun(_))
                            {
                                *occurrences.entry(power.get_base()).or_default() += 1;
                            }
                            return false;
                        }
                        AtomView::Add(_) if part != value.as_view() => {
                            *occurrences.entry(part).or_default() += 1;
                        }
                        _ => {}
                    }
                    true
                });
                let mut factors = occurrences
                    .into_iter()
                    .filter_map(|(factor, count)| (count > 1).then_some(factor))
                    .collect::<Vec<_>>();
                factors.sort_by(|left, right| {
                    right
                        .get_byte_size()
                        .cmp(&left.get_byte_size())
                        .then_with(|| left.cmp(right))
                });
                let mut keys = Vec::<Atom>::new();
                for factor in factors {
                    if keys
                        .iter()
                        .all(|key| !factor.contains(key) && !key.contains(factor))
                    {
                        keys.push(factor.to_owned());
                    }
                    if keys.len() > 128 {
                        break;
                    }
                }
                if !keys.is_empty() && keys.len() <= 128 {
                    // Large unrelated graph factors remain outside the row
                    // representation, including in this exact diagnostic.
                    let (spectator, selected) = if let AtomView::Mul(product) = value.as_view() {
                        let (selected, spectator): (Vec<_>, Vec<_>) = product
                            .iter()
                            .partition(|factor| keys.iter().any(|key| factor.contains(key)));
                        (Atom::mul_many(spectator), Atom::mul_many(selected))
                    } else {
                        (Atom::one(), value.clone())
                    };
                    if let Some(rows) = selected.coefficient_list_exact(&keys) {
                        let collected = spectator
                            * Atom::add_many(rows.into_iter().map(|(monomial, coefficient)| {
                                monomial * coefficient.collect_compact_factors()
                            }));
                        if collected.as_view().get_byte_size() < value.as_view().get_byte_size() {
                            value = collected;
                        }
                    }
                }
            }
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
    fn finite_numeric_coefficients_cancel_without_distributing_symbolic_factors() {
        test_initialise().unwrap();
        let mut cache = NumeratorJetCache::default();
        let p = parse_lit!(p);
        let k = parse_lit!(k);
        // Existing scalar names and tagged Casimirs remain literal coefficients
        // throughout finite component normalization.
        let reserved = Atom::var(symbol!("gammalooprs::uv::scalar_component_0"));
        for casimir in [Atom::one(), idenso::color_cas!(2, idenso::cof!(3))] {
            let shared = parse_lit!((a + b) * (c + d)) * casimir * &reserved;
            let difference = -(&shared * (-Atom::num(2) * &k)
                + &shared * (-Atom::num(2) * (&p - &k)))
                + &shared * (-Atom::num(2) * &p);
            assert!(cache.canonicalize(&difference).unwrap().is_zero);
            let retained = cache.canonicalize(&shared).unwrap();
            assert!(!retained.is_zero);
            assert!(
                retained
                    .key
                    .replace(shared.to_pattern())
                    .match_iter()
                    .next()
                    .is_some()
            );
        }
    }

    #[test]
    fn opposite_finite_component_factors_cancel_without_expanding_spectators() {
        test_initialise().unwrap();
        let mut cache = NumeratorJetCache::default();
        // Reduced (00) entry of the actual massive-quark bubble's U-US
        // comparison. y and w stand for the imaginary spatial components.
        let difference = parse_lit!(
            (-k0 + k3) * (2 * p0 - 2 * p3)
                + (-2 * p0 + 2 * p3) * (-k0 + k3)
                + (2 * x + 2 * y) * (x - y)
                + (-z + w) * (x + y)
                + (-z - w) * (x - y)
                + (-x + y) * (z + w)
                + (-x + y) * (2 * x + 2 * y)
                + (-x - y) * (z - w)
                + (-2 * x + 2 * y) * (x + y)
                + (-2 * x + 2 * y) * (-z - w)
                + (-2 * x + 2 * y) * (-x - y)
                + (-2 * x - 2 * y) * (-z + w)
        );
        let spectator =
            parse_lit!((a + b) * (c + d) * (a - b) ^ 11) * idenso::color_cas!(2, idenso::cof!(3));
        assert!(cache.canonicalize(&difference).unwrap().is_zero);
        assert!(
            cache
                .canonicalize(&(&spectator * &difference))
                .unwrap()
                .is_zero
        );
        let retained = cache.canonicalize(&spectator).unwrap();
        assert!(!retained.is_zero);
        assert_eq!(
            cache.canonicalize(&(&difference + &spectator)).unwrap(),
            retained
        );
        assert!(
            retained.key.contains(&parse_lit!((a - b) ^ 11))
                || retained.key.contains(&parse_lit!((-a + b) ^ 11))
        );
        assert!(retained.key.contains(&parse_lit!(a + b)));
        assert!(retained.key.contains(&parse_lit!(c + d)));
    }

    #[test]
    fn nested_parent_scalar_entry_cancels_by_bounded_factor_regrouping() {
        test_initialise().unwrap();
        // Reduced first scalar entry of the two-loop OS parent's vacuum
        // denominator class. All four binomial products remain factorized.
        let difference = parse_lit!(
            2 * (p1 - p2) * (2 * (p1 + p2) * (k0 + k3) - 2 * (k1 + k2) * (p0 + p3))
                + 2 * (p1 + p2) * (2 * (k1 - k2) * (p0 + p3) - 2 * (p1 - p2) * (k0 + k3))
                + 4 * (p0 + p3) * ((-k0 + k3) * (p0 + p3) + (p1 - p2) * (k1 + k2))
                - 4 * (p0 + p3) * ((-k0 + k3) * (p0 + p3) + (p1 + p2) * (k1 - k2))
        );
        let mut cache = NumeratorJetCache::default();
        assert!(cache.canonicalize(&difference).unwrap().is_zero);
        let spectator = parse_lit!((a + b) ^ 11) * idenso::color_cas!(2, idenso::cof!(3));
        assert!(
            cache
                .canonicalize(&(&spectator * &difference))
                .unwrap()
                .is_zero
        );
        let retained = parse_lit!((p1 - p2) * (a + b) ^ 11);
        assert!(
            !cache
                .canonicalize(&(&difference + &retained))
                .unwrap()
                .is_zero
        );
        let normalized = cache.canonicalize(&retained).unwrap().key;
        assert!(normalized.contains(&parse_lit!((a + b) ^ 11)));
        assert!(
            normalized.contains(&parse_lit!(p1 - p2)) || normalized.contains(&parse_lit!(-p1 + p2))
        );
    }

    #[test]
    fn opaque_common_factors_cancel_without_opening_their_bodies() {
        test_initialise().unwrap();
        let [a, b, c, d, e, f] = ["a", "b", "c", "d", "e", "f"].map(|name| {
            function!(
                symbol!(&format!("overlap_factor_test::{name}")),
                parse_lit!((x + y) ^ 11)
            )
        });
        let difference =
            &a * (&b * &c - &d * &e) + &b * (&f * &e - &a * &c) + &e * (&a * &d - &b * &f);
        let mut cache = NumeratorJetCache::default();
        assert!(cache.canonicalize(&difference).unwrap().is_zero);
        assert!(!cache.canonicalize(&(difference + &a)).unwrap().is_zero);
        assert!(cache.canonicalize(&a).unwrap().key.contains(&a));
    }

    #[test]
    fn repeated_scalar_products_keep_their_integer_power_relation() {
        test_initialise().unwrap();
        let contraction = parse_lit!(p0 * k0 - p1 * k1 - p2 * k2 - p3 * k3);
        let opposite = (-&contraction).expand_num();
        let spectator = parse_lit!((a + b) * (c + d) ^ 11);
        let difference = contraction.pow(2) * opposite + contraction.pow(3);
        let mut cache = NumeratorJetCache::default();
        assert!(cache.canonicalize(&difference).unwrap().is_zero);
        assert!(
            cache
                .canonicalize(&(&spectator * difference))
                .unwrap()
                .is_zero
        );
        let retained = cache
            .canonicalize(&(&spectator * contraction.pow(3)))
            .unwrap();
        assert!(!retained.is_zero);
        let retained = retained.key;
        assert!(retained.contains(&parse_lit!((c + d) ^ 11)));
        assert!(
            retained.contains(contraction.pow(3))
                || retained.contains((-&contraction).expand_num().pow(3))
        );
    }

    #[test]
    fn signed_integer_power_leaves_share_an_owner_without_expanding() {
        test_initialise().unwrap();
        let mut cache = NumeratorJetCache::default();
        let base = parse_lit!(a - b);
        let opposite = parse_lit!(-a + b);
        for exponent in [-2, -1, 2, 3, 11] {
            let left = base.pow(exponent);
            let right = Atom::num(if exponent % 2 == 0 { 1 } else { -1 }) * opposite.pow(exponent);
            assert_eq!(
                cache.canonicalize(&left).unwrap(),
                cache.canonicalize(&right).unwrap(),
            );
            assert!(cache.canonicalize(&(left - right)).unwrap().is_zero);
        }
        let fractional = parse_lit!((a - b) ^ (1 / 2));
        assert!(
            cache
                .canonicalize(&fractional)
                .unwrap()
                .key
                .contains(&fractional)
        );
    }

    #[test]
    fn finite_component_certificates_preserve_opaque_numeric_sum_payloads() {
        test_initialise().unwrap();
        let body = parse_lit!(c + 2 * (a + b));
        let fractional = body.pow(parse_lit!(1 / 2));
        let inverse = body.pow(-1);
        let function = function!(symbol!("finite_component_opaque_payload"), &body);
        let mut cache = NumeratorJetCache::default();
        for value in [fractional, inverse, function] {
            let before = value.clone();
            let certificate = cache.canonicalize(&value).unwrap();
            assert!(!certificate.is_zero);
            assert!(certificate.key.contains(&value));
            assert_eq!(value, before);
        }
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
