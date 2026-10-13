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

    fn normalize_scalar_component(value: AtomView<'_>) -> Atom {
        // This value is already a finite scalar tensor entry. A scalar
        // Casimir may still contain representation metadata in its arguments;
        // whole-graph tensor inference must not classify it as an open tensor.
        // Freeze function and branch-sensitive power leaves into an injective
        // neutral alphabet. Integer powers stay compressed algebraic powers,
        // so x*x^2=x^3 remains exact without distributing their bases.
        let normalize_sum = |sum: AtomView<'_>| {
            let AtomView::Add(sum) = sum else {
                unreachable!("only additive factors have a numeric pivot");
            };
            let (number, _) = sum
                .iter()
                .map(|term| Self::split_number(term.to_owned()))
                .min_by(|(_, left), (_, right)| left.cmp(right))
                .expect("an additive factor is nonempty");
            let body = (Atom::from(sum.to_owned()) / &number).expand_num();
            (number, body)
        };
        let mut used = value.get_all_symbols(true);
        let mut aliases = HashMap::<Atom, Atom>::new();
        let mut serial = 0usize;
        let protected = value.replace_map(|view, _, output| {
            let opaque = match view {
                AtomView::Fun(_) => true,
                AtomView::Pow(power) => i64::try_from(power.get_exp()).is_err(),
                _ => false,
            };
            if opaque {
                let alias = aliases.entry(view.to_owned()).or_insert_with(|| {
                    loop {
                        let symbol =
                            symbol!(&format!("gammalooprs::uv::scalar_component_{serial}"));
                        serial += 1;
                        if used.insert(symbol) {
                            break Atom::var(symbol);
                        }
                    }
                });
                **output = alias.clone();
            }
        });
        let mut normalized = protected;
        // Finite contractions can leave opposite additive factors whose
        // common factor appears only after a numeric sign is distributed.
        // Iterate exact factor collection to a fixed point. Function leaves
        // and compressed integer powers retain every graph factor.
        for _ in 0..32 {
            let signed = normalized.replace_map_bottom_up(|view, context, output| {
                if matches!(view, AtomView::Add(_))
                    && context.parent_type != Some(symbolica::atom::AtomType::Pow)
                {
                    let (number, body) = normalize_sum(view);
                    **output = number * body;
                    return;
                }
                if let AtomView::Pow(power) = view
                    && let Ok(exponent) = i64::try_from(power.get_exp())
                    && matches!(power.get_base(), AtomView::Add(_))
                {
                    let (number, body) = normalize_sum(power.get_base());
                    **output = number.pow(exponent) * body.pow(exponent);
                }
            });
            let next = signed.collect_factors().expand_num();
            if next == normalized {
                break;
            }
            normalized = next;
        }
        normalized.replace_multiple(
            aliases
                .iter()
                .map(|(source, alias)| Replacement::new(alias.to_pattern(), source.to_pattern())),
        )
    }

    /// Test certificate by alternative common-factor groupings. Failure to
    /// find a proof is not evidence of a nonzero polynomial.
    #[cfg(test)]
    pub(crate) fn factorized_is_zero(atom: Atom) -> bool {
        // Split along one complete common factor only, including an opaque
        // function or compressed power as an indivisible factor. Other products and
        // opaque functions/powers stay intact; this never builds monomials of
        // a graph numerator or distributes all of its additive factors.
        fn split(atom: AtomView<'_>, factor: &Atom, remaining: &mut usize) -> Option<(Atom, Atom)> {
            *remaining = remaining.checked_sub(1)?;
            if atom == factor.as_view() {
                return Some((Atom::one(), Atom::Zero));
            }
            match atom {
                AtomView::Add(sum) => {
                    let mut coefficient = Atom::Zero;
                    let mut residual = Atom::Zero;
                    for term in sum.iter() {
                        let (part, remainder) = split(term, factor, remaining)?;
                        coefficient += part;
                        residual += remainder;
                    }
                    Some((coefficient, residual))
                }
                AtomView::Mul(product) => {
                    let terms = product.iter().collect::<Vec<_>>();
                    for (index, term) in terms.iter().enumerate() {
                        let (coefficient, residual) = split(*term, factor, remaining)?;
                        if coefficient.is_zero() && residual.as_view() == *term {
                            continue;
                        }
                        let other =
                            Atom::mul_many(terms.iter().enumerate().filter_map(
                                |(position, term)| (position != index).then_some(*term),
                            ));
                        return Some((&other * coefficient, other * residual));
                    }
                    Some((Atom::Zero, atom.to_owned()))
                }
                _ => Some((Atom::Zero, atom.to_owned())),
            }
        }

        fn search(atom: Atom, visited: &mut std::collections::HashSet<Atom>) -> bool {
            if atom.as_view().get_byte_size() > MAX_BODY_BYTES {
                return false;
            }
            let atom = NumeratorJetCache::normalize_scalar_component(atom.as_view());
            if atom.is_zero() {
                return true;
            }
            if visited.len() >= 256 || !visited.insert(atom.clone()) {
                return false;
            }
            // A complete zero factor proves a scalar product zero. Testing
            // that smaller factor avoids retaining irrelevant opaque coupling,
            // color and denominator-jet spectators throughout the search.
            if let AtomView::Mul(product) = atom.as_view() {
                for factor in product.iter() {
                    if matches!(factor, AtomView::Add(_)) && search(factor.to_owned(), visited) {
                        return true;
                    }
                }
            }
            let mut factors = std::collections::BTreeSet::new();
            atom.visitor(&mut |view| {
                if let AtomView::Mul(product) = view {
                    factors.extend(
                        product
                            .iter()
                            .filter(|factor| {
                                matches!(
                                    factor,
                                    AtomView::Add(_) | AtomView::Fun(_) | AtomView::Pow(_)
                                )
                            })
                            .map(|factor| factor.to_owned()),
                    );
                }
                !matches!(view, AtomView::Fun(_) | AtomView::Pow(_))
            });
            for factor in factors {
                let Some((coefficient, residual)) = split(atom.as_view(), &factor, &mut 4096)
                else {
                    continue;
                };
                if coefficient.as_view().get_byte_size() + residual.as_view().get_byte_size()
                    > MAX_BODY_BYTES
                {
                    continue;
                }
                let coefficient =
                    NumeratorJetCache::normalize_scalar_component(coefficient.as_view());
                let residual = NumeratorJetCache::normalize_scalar_component(residual.as_view());
                let candidate = if coefficient.is_zero() {
                    residual
                } else {
                    factor * coefficient + residual
                };
                let candidate = NumeratorJetCache::normalize_scalar_component(candidate.as_view());
                if candidate.as_view().get_byte_size() < atom.as_view().get_byte_size()
                    && search(candidate, visited)
                {
                    return true;
                }
            }
            let mut sums = Vec::new();
            atom.visitor(&mut |view| {
                if matches!(view, AtomView::Fun(_) | AtomView::Pow(_)) {
                    return false;
                }
                if matches!(view, AtomView::Add(_)) {
                    sums.push(view.to_owned());
                }
                true
            });
            for sum in sums {
                let AtomView::Add(sum_view) = sum.as_view() else {
                    unreachable!();
                };
                let terms = sum_view.iter().collect::<Vec<_>>();
                for (i, left) in terms.iter().enumerate() {
                    for (j, right) in terms.iter().enumerate().skip(i + 1) {
                        let original_pair = *left + *right;
                        let pair =
                            NumeratorJetCache::normalize_scalar_component(original_pair.as_view());
                        if pair == original_pair {
                            continue;
                        }
                        let regrouped = Atom::add_many(
                            terms
                                .iter()
                                .enumerate()
                                .filter_map(|(k, term)| (k != i && k != j).then_some(*term))
                                .chain([pair.as_view()]),
                        );
                        let candidate = atom.replace_map(|view, _, output| {
                            if matches!(view, AtomView::Fun(_) | AtomView::Pow(_)) {
                                **output = view.to_owned();
                            } else if view == sum.as_view() {
                                **output = regrouped.clone();
                            }
                        });
                        let candidate =
                            NumeratorJetCache::normalize_scalar_component(candidate.as_view());
                        if candidate.as_view().get_byte_size() < atom.as_view().get_byte_size()
                            && search(candidate, visited)
                        {
                            return true;
                        }
                    }
                }
            }
            false
        }
        search(atom, &mut Default::default())
    }

    fn contract(body: &Atom) -> Result<CanonicalTensor> {
        let started = std::time::Instant::now();
        crate::debug_tags!(#generation, #profile, #uv, #numerator, #summary;
            stage = "numerator_certificate_parse_start",
            input_bytes = body.as_view().get_byte_size(),
            "Parsing the factorized numerator certificate"
        );
        let mut network = body.parse_into_net()?;
        crate::debug_tags!(#generation, #profile, #uv, #numerator, #summary;
            stage = "numerator_certificate_parse_done",
            elapsed_ms = started.elapsed().as_secs_f64() * 1000.0,
            "Parsed the factorized numerator certificate"
        );
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
        let started = std::time::Instant::now();
        crate::debug_tags!(#generation, #profile, #uv, #numerator, #summary;
            stage = "numerator_certificate_sum_boundaries_start",
            external_capacity = capacity,
            "Contracting ready numerator sum boundaries"
        );
        network.graph.contract_ready_sum_boundaries();
        crate::debug_tags!(#generation, #profile, #uv, #numerator, #summary;
            stage = "numerator_certificate_sum_boundaries_done",
            elapsed_ms = started.elapsed().as_secs_f64() * 1000.0,
            "Contracted ready numerator sum boundaries"
        );
        let started = std::time::Instant::now();
        crate::debug_tags!(#generation, #profile, #uv, #numerator, #summary;
            stage = "numerator_certificate_execute_start",
            strategy = "SequentialRef", external_capacity = capacity,
            "Executing the finite numerator certificate network"
        );
        network.execute::<SequentialRef, MinResultRank, _, _, _>(
            TENSORLIB.read().unwrap().deref(),
            FUN_LIB.deref(),
        )?;
        crate::debug_tags!(#generation, #profile, #uv, #numerator, #summary;
            stage = "numerator_certificate_execute_done",
            strategy = "SequentialRef",
            elapsed_ms = started.elapsed().as_secs_f64() * 1000.0,
            "Executed the finite numerator certificate network"
        );
        let tensor = match network.result_tensor(TENSORLIB.read().unwrap().deref())? {
            ExecutionResult::One => return Ok(CanonicalTensor::number(Atom::num(1))),
            ExecutionResult::Zero => return Ok(CanonicalTensor::number(Atom::Zero)),
            ExecutionResult::Val(tensor) => tensor.into_owned(),
        };
        let tensor_size = tensor.size()?;
        ensure!(
            tensor_size <= MAX_COMPONENTS,
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
        let started = std::time::Instant::now();
        crate::debug_tags!(#generation, #profile, #uv, #numerator, #summary;
            stage = "numerator_certificate_normalize_start",
            tensor_size,
            "Normalizing finite tensor certificate components"
        );
        for (index, value) in tensor.iter_flat() {
            // Finite tensor execution can leave nested numeric signs and
            // rational weights, such as -(C*(-2*k)+C*(-2*(p-k)))+C*(-2*p).
            // Distribute only numbers; symbolic products and function bodies
            // remain factorized, including the original numerator spectators.
            let mut value = Self::normalize_scalar_component(value);
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
        crate::debug_tags!(#generation, #profile, #uv, #numerator, #summary;
            stage = "numerator_certificate_normalize_done",
            elapsed_ms = started.elapsed().as_secs_f64() * 1000.0,
            nonzero_components = components.len(), component_bytes = bytes,
            "Normalized finite tensor certificate components"
        );
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
        // Reserving every source symbol must also protect deliberate overlap
        // with the temporary scalar alphabet used by this exact proof.
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
        assert!(NumeratorJetCache::factorized_is_zero(difference.clone()));
        let spectator = parse_lit!((a + b) ^ 11) * idenso::color_cas!(2, idenso::cof!(3));
        assert!(NumeratorJetCache::factorized_is_zero(
            &spectator * &difference
        ));
        let retained = parse_lit!((p1 - p2) * (a + b) ^ 11);
        assert!(!NumeratorJetCache::factorized_is_zero(
            &difference + &retained
        ));
        let normalized = NumeratorJetCache::normalize_scalar_component(retained.as_view());
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
        assert!(NumeratorJetCache::factorized_is_zero(difference.clone()));
        assert!(!NumeratorJetCache::factorized_is_zero(difference + &a));
        assert_eq!(
            NumeratorJetCache::normalize_scalar_component(a.as_view()),
            a,
        );
    }

    #[test]
    fn repeated_scalar_products_keep_their_integer_power_relation() {
        test_initialise().unwrap();
        let contraction = parse_lit!(p0 * k0 - p1 * k1 - p2 * k2 - p3 * k3);
        let opposite = (-&contraction).expand_num();
        let spectator = parse_lit!((a + b) * (c + d) ^ 11);
        let difference = contraction.pow(2) * opposite + contraction.pow(3);
        assert!(NumeratorJetCache::normalize_scalar_component(difference.as_view()).is_zero());
        assert!(NumeratorJetCache::factorized_is_zero(
            &spectator * difference
        ));
        let retained = NumeratorJetCache::normalize_scalar_component(
            (&spectator * contraction.pow(3)).as_view(),
        );
        assert!(!retained.is_zero());
        assert!(retained.contains(&parse_lit!((c + d) ^ 11)));
        assert!(
            retained.contains(contraction.pow(3))
                || retained.contains((-&contraction).expand_num().pow(3))
        );
    }

    #[test]
    fn signed_integer_power_leaves_share_an_owner_without_expanding() {
        test_initialise().unwrap();
        let base = parse_lit!(a - b);
        let opposite = parse_lit!(-a + b);
        for exponent in [-2, -1, 2, 3, 11] {
            let left = base.pow(exponent);
            let right = Atom::num(if exponent % 2 == 0 { 1 } else { -1 }) * opposite.pow(exponent);
            assert_eq!(
                NumeratorJetCache::normalize_scalar_component(left.as_view()),
                NumeratorJetCache::normalize_scalar_component(right.as_view()),
            );
            assert!(
                NumeratorJetCache::normalize_scalar_component((left - right).as_view()).is_zero()
            );
        }
        let fractional = parse_lit!((a - b) ^ (1 / 2));
        assert_eq!(
            NumeratorJetCache::normalize_scalar_component(fractional.as_view()),
            fractional,
        );
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
