//! Symbolic Lorentz contractions on factorized expressions.

use std::collections::BTreeMap;

use feynkit_tensor::{OrthogonalWeingarten, TensorReductionError};
use spenso::network::tags::SPENSO_TAG;
use symbolica::atom::{Atom, AtomCore, AtomView, FunctionBuilder};

use crate::{VakintError, VakintTerm, symbols::S};

/// Original multiplicities survive contraction: a scalar power cannot hide an
/// index from a consumer in its enclosing product.
pub(crate) struct LorentzTensor {
    pub(crate) expression: Atom,
    indices: BTreeMap<Atom, u8>,
    native: bool,
}

type SourceSlots = BTreeMap<Atom, Option<Atom>>;

impl LorentzTensor {
    pub(crate) fn parse(view: AtomView<'_>, dimension: &Atom) -> Result<Self, VakintError> {
        Self::parse_with_slots(view, dimension, &Self::source_slots(view))
    }

    fn source_slots(view: AtomView<'_>) -> SourceSlots {
        let mut slots = SourceSlots::new();
        view.visitor(&mut |value| {
            if matches!(value, AtomView::Fun(function) if function.get_symbol() == SPENSO_TAG.pure_scalar)
            {
                return false;
            }
            if let Some((_, index)) = VakintTerm::spenso_minkowski_slot_parts(value) {
                slots
                    .entry(index.to_owned())
                    .and_modify(|slot| {
                        if slot.as_ref().is_some_and(|slot| slot.as_view() != value) {
                            *slot = None;
                        }
                    })
                    .or_insert_with(|| Some(value.to_owned()));
                return false;
            }
            if let Some(slot) = Self::expanded_native_slot(value) {
                let (_, index) = VakintTerm::spenso_minkowski_slot_parts(slot.as_view())
                    .expect("expanded vector slot has its declared representation");
                slots
                    .entry(index.to_owned())
                    .and_modify(|original| {
                        if original.as_ref().is_some_and(|original| original != &slot) {
                            *original = None;
                        }
                    })
                    .or_insert(Some(slot));
            }
            true
        });
        slots
    }

    fn expanded_native_slot(view: AtomView<'_>) -> Option<Atom> {
        let AtomView::Fun(function) = view else {
            return None;
        };
        let arguments = function.iter().collect::<Vec<_>>();
        if arguments.len() < 2
            || !VakintTerm::is_compact_spenso_minkowski_representation(
                arguments[arguments.len() - 2],
            )
            || !VakintTerm::is_dot_dummy_index(arguments[arguments.len() - 1])
        {
            return None;
        }
        let AtomView::Fun(representation) = arguments[arguments.len() - 2] else {
            return None;
        };
        Some(
            FunctionBuilder::new(representation.get_symbol())
                .add_args(representation.iter())
                .add_arg(arguments[arguments.len() - 1])
                .finish(),
        )
    }

    fn index_key(index: AtomView<'_>, slots: &SourceSlots) -> Result<Atom, VakintError> {
        if VakintTerm::spenso_minkowski_slot_parts(index).is_some() {
            return Ok(index.to_owned());
        }
        match slots.get(&index.to_owned()) {
            Some(Some(slot)) => Ok(slot.clone()),
            Some(None) => Err(VakintError::InvalidNumerator(format!(
                "unrepresented Lorentz index {index} has multiple source dimensions"
            ))),
            None => Ok(index.to_owned()),
        }
    }

    fn parse_with_slots(
        view: AtomView<'_>,
        dimension: &Atom,
        slots: &SourceSlots,
    ) -> Result<Self, VakintError> {
        match view {
            AtomView::Add(sum) => Self::sum(
                sum.iter()
                    .map(|term| {
                        let mut branch_slots = slots.clone();
                        branch_slots.extend(Self::source_slots(term));
                        Self::parse_with_slots(term, dimension, &branch_slots)
                    })
                    .collect::<Result<Vec<_>, _>>()?,
            ),
            AtomView::Mul(product) => Self::product(
                product
                    .iter()
                    .map(|factor| Self::parse_with_slots(factor, dimension, slots))
                    .collect::<Result<Vec<_>, _>>()?,
                dimension,
                slots,
            ),
            AtomView::Pow(power) => {
                let (base, exponent) = power.get_base_exp();
                let original_base = base;
                let base = Self::parse_with_slots(base, dimension, slots)?;
                let exponent = Self::parse_with_slots(exponent, dimension, slots)?;
                if !exponent.is_scalar() {
                    return Err(VakintError::InvalidNumerator(format!(
                        "open Lorentz indices inside a scalar exponent: {view}"
                    )));
                }
                if base.is_scalar() {
                    return Ok(Self {
                        expression: base.expression.pow(exponent.expression),
                        indices: base.indices,
                        native: base.native,
                    });
                }
                // Preserve Vakint's existing identical-integrated-power policy.
                // The bound applies to this indexed vector, never to a power
                // of a scalar contraction that was already validated above.
                if matches!(original_base, AtomView::Fun(vector)
                    if vector.get_symbol() == S.k && base.indices.len() == 1)
                    && let Ok(count) = usize::try_from(exponent.expression.as_view())
                {
                    if count > OrthogonalWeingarten::MAX_RANK {
                        return Err(TensorReductionError::UnsupportedRank {
                            rank: count,
                            maximum: OrthogonalWeingarten::MAX_RANK,
                        }
                        .into());
                    }
                    if count != 0 && count.is_multiple_of(2) {
                        let mut copy = base;
                        copy.indices.retain(|_, count| *count == 1);
                        let other = Self {
                            expression: copy.expression.clone(),
                            indices: copy.indices.clone(),
                            native: copy.native,
                        };
                        let mut closed = Self::product(vec![copy, other], dimension, slots)?;
                        closed.expression = if closed.native {
                            view.to_owned()
                        } else {
                            closed.expression.pow(Atom::num(count / 2))
                        };
                        return Ok(closed);
                    }
                }
                if !matches!(i64::try_from(exponent.expression.as_view()), Ok(2)) {
                    return Err(VakintError::InvalidNumerator(format!(
                        "an indexed tensor can only have power two: {view}"
                    )));
                }
                let mut copy = Self::parse_with_slots(base.expression.as_view(), dimension, slots)?;
                copy.indices.retain(|_, count| *count == 1);
                let other = Self {
                    expression: copy.expression.clone(),
                    indices: copy.indices.clone(),
                    native: copy.native,
                };
                Self::product(vec![copy, other], dimension, slots)
            }
            AtomView::Fun(function)
                if [S.dot, SPENSO_TAG.dot, SPENSO_TAG.pure_scalar]
                    .contains(&function.get_symbol()) =>
            {
                // These are already closed notation or an explicitly opaque
                // scalar payload. The checked reducer still validates dots.
                Ok(Self::scalar(view.to_owned()))
            }
            AtomView::Fun(function) if function.get_symbol() == S.g => {
                if function.get_nargs() != 2 {
                    return Err(VakintError::InvalidNumerator(format!(
                        "a Lorentz metric requires two indices: {view}"
                    )));
                }
                let left = function.get(0);
                let right = function.get(1);
                let mut indices = BTreeMap::new();
                *indices.entry(Self::index_key(left, slots)?).or_insert(0) += 1;
                *indices.entry(Self::index_key(right, slots)?).or_insert(0) += 1;
                let native = [left, right]
                    .iter()
                    .any(|slot| VakintTerm::spenso_minkowski_slot_parts(*slot).is_some());
                Ok(Self {
                    expression: if !native && left == right {
                        dimension.clone()
                    } else {
                        view.to_owned()
                    },
                    indices,
                    native,
                })
            }
            AtomView::Fun(function)
                if [S.k, S.p].contains(&function.get_symbol()) && function.get_nargs() == 2 =>
            {
                let index = function.get(1);
                if VakintTerm::is_compact_spenso_minkowski_representation(index) {
                    return Ok(Self::scalar(view.to_owned()));
                }
                Ok(Self {
                    expression: view.to_owned(),
                    indices: BTreeMap::from([(Self::index_key(index, slots)?, 1)]),
                    native: VakintTerm::spenso_minkowski_slot_parts(index).is_some(),
                })
            }
            AtomView::Fun(function) => {
                let mut indices = BTreeMap::new();
                if let Some(slot) = Self::expanded_native_slot(view) {
                    indices.insert(slot, 1_u8);
                }
                for argument in function.iter() {
                    if VakintTerm::spenso_minkowski_slot_parts(argument).is_some() {
                        let count = indices.entry(argument.to_owned()).or_insert(0_u8);
                        *count += 1;
                        if *count > 2 {
                            return Err(VakintError::InvalidNumerator(format!(
                                "Lorentz index {argument} is contracted more than once"
                            )));
                        }
                    }
                }
                if !indices.is_empty() {
                    // Preserve complete native leaves and their metadata. Their
                    // contractions belong to the existing checked tensor owner.
                    return Ok(Self {
                        expression: view.to_owned(),
                        indices,
                        native: true,
                    });
                }
                let mut result = FunctionBuilder::new(function.get_symbol());
                for argument in function {
                    let argument = Self::parse_with_slots(argument, dimension, slots)?;
                    if !argument.is_scalar() {
                        return Err(VakintError::InvalidNumerator(format!(
                            "open Lorentz indices inside a scalar function: {view}"
                        )));
                    }
                    result = result.add_arg(argument.expression);
                }
                Ok(Self::scalar(result.finish()))
            }
            _ => Ok(Self::scalar(view.to_owned())),
        }
    }

    fn scalar(expression: Atom) -> Self {
        Self {
            expression,
            indices: BTreeMap::new(),
            native: false,
        }
    }

    pub(crate) fn is_scalar(&self) -> bool {
        self.indices.values().all(|count| *count == 2)
    }

    fn shared_indices(&self, other: &Self) -> Vec<Atom> {
        self.indices
            .iter()
            .filter(|(_, count)| **count == 1)
            .filter(|(index, _)| other.indices.get(*index) == Some(&1))
            .map(|(index, _)| index.clone())
            .collect()
    }

    fn sum(terms: Vec<Self>) -> Result<Self, VakintError> {
        let mut indices = BTreeMap::new();
        let mut external = None;
        let mut expression = Atom::zero();
        let mut native = false;
        let mut mixed_interface = false;
        for term in terms {
            let free = term
                .indices
                .iter()
                .filter(|(_, count)| **count == 1)
                .map(|(index, _)| index.clone())
                .collect::<Vec<_>>();
            if external.as_ref().is_some_and(|expected| expected != &free) {
                // Vakint also accepts formal sums of distinct tensor ranks.
                // Their projection remains with the checked native reducer;
                // do not distribute or contract this mixed component here.
                mixed_interface = true;
            }
            external = Some(free);
            for (index, count) in term.indices {
                let original = indices.entry(index).or_insert(0_u8);
                *original = (*original).max(count);
            }
            expression += term.expression;
            native |= term.native;
        }
        Ok(Self {
            expression,
            indices,
            native: native || mixed_interface,
        })
    }

    fn product(
        mut factors: Vec<Self>,
        dimension: &Atom,
        slots: &SourceSlots,
    ) -> Result<Self, VakintError> {
        let mut indices = BTreeMap::new();
        for factor in &factors {
            for (index, count) in &factor.indices {
                let total = indices.entry(index.clone()).or_insert(0_u8);
                *total += count;
                if *total > 2 {
                    return Err(VakintError::InvalidNumerator(format!(
                        "Lorentz index {index} is contracted more than once"
                    )));
                }
            }
        }

        // Only connected index components enter bilinear contraction. Closed
        // scalar factors and products of their sums remain factorized.
        loop {
            let pair = (0..factors.len())
                .flat_map(|left| (left + 1..factors.len()).map(move |right| (left, right)))
                .map(|(left, right)| {
                    (
                        factors[left].shared_indices(&factors[right]).len(),
                        left,
                        right,
                    )
                })
                .filter(|(shared, _, _)| *shared != 0)
                .max();
            let Some((_, left, right)) = pair else { break };
            let right = factors.remove(right);
            let left = factors.remove(left);
            factors.push(left.contract(right, dimension, slots)?);
        }
        let native = factors.iter().any(|factor| factor.native);
        Ok(Self {
            expression: factors
                .into_iter()
                .fold(Atom::one(), |product, factor| product * factor.expression),
            indices,
            native,
        })
    }

    fn contract(
        self,
        other: Self,
        dimension: &Atom,
        slots: &SourceSlots,
    ) -> Result<Self, VakintError> {
        if self.native || other.native {
            // Propagate opacity through the entire connected component. In
            // particular, a second attachment to T(mu,nu)*k(mu) must not split
            // that product and recursively rediscover the first attachment.
            let mut indices = self.indices;
            for (index, count) in other.indices {
                *indices.entry(index).or_insert(0) += count;
            }
            return Ok(Self {
                expression: self.expression * other.expression,
                indices,
                native: true,
            });
        }
        if self.expression.is_zero() || other.expression.is_zero() {
            let mut indices = BTreeMap::new();
            for tensor in [&self, &other] {
                for (index, count) in &tensor.indices {
                    if *count == 1 {
                        *indices.entry(index.clone()).or_insert(0) += count;
                    }
                }
            }
            return Ok(Self {
                expression: Atom::zero(),
                indices,
                native: false,
            });
        }
        for (left, right) in [(&self, &other), (&other, &self)] {
            match left.expression.as_view() {
                AtomView::Mul(product) => {
                    let mut factors = product
                        .iter()
                        .map(|factor| Self::parse_with_slots(factor, dimension, slots))
                        .collect::<Result<Vec<_>, _>>()?;
                    factors.push(Self::parse_with_slots(
                        right.expression.as_view(),
                        dimension,
                        slots,
                    )?);
                    return Self::product(factors, dimension, slots);
                }
                AtomView::Add(sum) => {
                    return Self::sum(
                        sum.iter()
                            .map(|term| {
                                Self::parse_with_slots(term, dimension, slots)?.contract(
                                    Self::parse_with_slots(
                                        right.expression.as_view(),
                                        dimension,
                                        slots,
                                    )?,
                                    dimension,
                                    slots,
                                )
                            })
                            .collect::<Result<Vec<_>, _>>()?,
                    );
                }
                _ => {}
            }
        }
        let shared = self.shared_indices(&other);
        let (AtomView::Fun(left), AtomView::Fun(right)) =
            (self.expression.as_view(), other.expression.as_view())
        else {
            return Err(VakintError::InvalidNumerator(format!(
                "cannot contract {} with {}",
                self.expression, other.expression
            )));
        };
        let nonshared = |function: symbolica::atom::representation::FunView<'_>| {
            function
                .iter()
                .find(|index| Self::index_key(*index, slots).is_ok_and(|index| index != shared[0]))
                .map(|index| index.to_owned())
        };
        let expression = if left.get_symbol() == S.g && right.get_symbol() == S.g {
            if shared.len() == 2 {
                dimension.clone()
            } else {
                FunctionBuilder::new(S.g)
                    .add_arg(nonshared(left).expect("metric has one uncontracted index"))
                    .add_arg(nonshared(right).expect("metric has one uncontracted index"))
                    .finish()
            }
        } else if left.get_symbol() == S.g || right.get_symbol() == S.g {
            let (metric, vector) = if left.get_symbol() == S.g {
                (left, right)
            } else {
                (right, left)
            };
            FunctionBuilder::new(vector.get_symbol())
                .add_arg(vector.get(0))
                .add_arg(nonshared(metric).expect("metric has one uncontracted index"))
                .finish()
        } else {
            S.dot(
                FunctionBuilder::new(left.get_symbol())
                    .add_arg(left.get(0))
                    .finish(),
                FunctionBuilder::new(right.get_symbol())
                    .add_arg(right.get(0))
                    .finish(),
            )
        };
        Self::parse_with_slots(expression.as_view(), dimension, slots)
    }
}

#[cfg(test)]
mod tests {
    use spenso::{
        network::library::symbolic::ETS,
        structure::representation::{Minkowski, RepName},
    };
    use symbolica::{function, symbol};

    use crate::{Vakint, utils::vakint_macros::vk_parse};

    use super::*;

    fn dimension() -> Atom {
        Vakint::initialize_vakint_symbols();
        vk_parse!("4-2*ε").unwrap()
    }

    #[test]
    fn closed_scalar_bases_keep_powers_and_spectators_factorized() {
        let dimension = dimension();
        let source =
            vk_parse!("(alpha+beta)^7*(k(1,mu)^2-2*k(1,nu)*k(2,nu)+k(2,rho)^2-1)^11").unwrap();
        let expected =
            vk_parse!("(alpha+beta)^7*(dot(k(1),k(1))-2*dot(k(1),k(2))+dot(k(2),k(2))-1)^11")
                .unwrap();
        let parsed = LorentzTensor::parse(source.as_view(), &dimension).unwrap();
        assert!(parsed.is_scalar());
        assert_eq!(parsed.expression, expected);
    }

    #[test]
    fn original_index_multiplicities_survive_closed_powers() {
        let dimension = dimension();
        for source in [
            "k(1,mu)^2*k(2,mu)^2",
            "k(1,mu)^2*p(1,mu)",
            "(k(1,mu)^2+1)^3*k(2,mu)",
            "g(mu,mu)*k(1,mu)",
            "(k(1,mu)^2+k(2,mu))*p(1,mu)",
            "(k(1,mu)+k(2,mu)^2)*p(1,mu)",
        ] {
            assert!(
                LorentzTensor::parse(vk_parse!(source).unwrap().as_view(), &dimension).is_err(),
                "{source}"
            );
        }
    }

    #[test]
    fn mixed_tensor_interfaces_are_deferred_without_losing_branch_multiplicity() {
        let dimension = dimension();
        let source = vk_parse!("(k(1,mu)+1)^2").unwrap();
        assert_eq!(
            LorentzTensor::parse(source.as_view(), &dimension)
                .unwrap()
                .expression,
            source
        );
        let source = vk_parse!("k(1,mu)^2+k(2,nu)*k(2,rho)").unwrap();
        let expected = vk_parse!("dot(k(1),k(1))+k(2,nu)*k(2,rho)").unwrap();
        assert_eq!(
            LorentzTensor::parse(source.as_view(), &dimension)
                .unwrap()
                .expression,
            expected
        );
    }

    #[test]
    fn legacy_identical_integrated_powers_keep_the_existing_bound() {
        let dimension = dimension();
        for count in [2, 4, 6, OrthogonalWeingarten::MAX_RANK] {
            let source = vk_parse!(&format!("k(1,mu)^{count}")).unwrap();
            let expected = vk_parse!(&format!("dot(k(1),k(1))^{}", count / 2)).unwrap();
            assert_eq!(
                LorentzTensor::parse(source.as_view(), &dimension)
                    .unwrap()
                    .expression,
                expected
            );
        }
        assert!(matches!(
            LorentzTensor::parse(vk_parse!("k(1,mu)^22").unwrap().as_view(), &dimension),
            Err(VakintError::FeynKitTensorReduction(
                TensorReductionError::UnsupportedRank {
                    rank: 22,
                    maximum: 20
                }
            ))
        ));
        for source in ["k(1,mu)^3", "p(1,mu)^4"] {
            assert!(
                LorentzTensor::parse(vk_parse!(source).unwrap().as_view(), &dimension).is_err()
            );
        }
    }

    #[test]
    fn native_two_attachment_components_stay_opaque_and_terminate() {
        let dimension = dimension();
        let mu = vk_parse!("mu").unwrap();
        let nu = vk_parse!("nu").unwrap();
        let mink = Minkowski {}.new_rep(6);
        let tensor = symbol!(
            "lorentz_scope_test::two_ports",
            tags = [SPENSO_TAG.tensor.clone()]
        );
        let body = vk_parse!("payload(mu)").unwrap();
        let native = function!(
            tensor,
            &body,
            mink.to_symbolic([mu.clone()]),
            mink.to_symbolic([nu.clone()])
        );
        let source = &native * function!(S.k, 1, &mu) * function!(S.k, 2, &nu);
        let parsed = LorentzTensor::parse(source.as_view(), &dimension).unwrap();
        assert!(parsed.is_scalar());
        assert_eq!(parsed.expression, source);

        let overcontracted = &native * function!(S.k, 1, &mu).pow(2);
        assert!(LorentzTensor::parse(overcontracted.as_view(), &dimension).is_err());
    }

    #[test]
    fn represented_metrics_and_distinct_dimensions_keep_their_metadata() {
        let dimension = dimension();
        let mu = vk_parse!("mu").unwrap();
        let six = Minkowski {}.new_rep(6).to_symbolic([mu.clone()]);
        let eight = Minkowski {}.new_rep(8).to_symbolic([mu]);
        for source in [
            function!(ETS.metric, &six, &six),
            function!(S.g, &six, &six),
            function!(S.p, 1, &six) * function!(S.p, 2, &eight),
        ] {
            assert_eq!(
                LorentzTensor::parse(source.as_view(), &dimension)
                    .unwrap()
                    .expression,
                source
            );
        }
    }

    #[test]
    fn expanded_native_vector_ports_keep_their_source_scope() {
        let dimension = dimension();
        let index = S.dot_dummy_ind(91);
        let vector = symbol!(
            "lorentz_scope_test::expanded_vector",
            tags = [SPENSO_TAG.rank1.clone(), SPENSO_TAG.tensor.clone()]
        );
        let native = function!(vector, Minkowski {}.new_rep(6).to_symbolic([]), &index);
        let source = &native * function!(S.k, 1, &index);
        let parsed = LorentzTensor::parse(source.as_view(), &dimension).unwrap();
        assert!(parsed.is_scalar());
        assert_eq!(parsed.expression, source);

        let overcontracted = &native * function!(S.k, 1, index).pow(2);
        assert!(LorentzTensor::parse(overcontracted.as_view(), &dimension).is_err());
    }

    #[test]
    fn metrics_close_only_their_connected_raw_component() {
        let dimension = dimension();
        let source = vk_parse!("(alpha+beta)^7*g(mu,nu)*k(1,mu)*k(2,nu)*(k(3,rho)^2+1)^5").unwrap();
        let expected = vk_parse!("(alpha+beta)^7*dot(k(1),k(2))*(dot(k(3),k(3))+1)^5").unwrap();
        assert_eq!(
            LorentzTensor::parse(source.as_view(), &dimension)
                .unwrap()
                .expression,
            expected
        );
    }
}
