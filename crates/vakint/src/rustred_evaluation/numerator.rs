//! Transport scalarized numerators with an already validated parent witness.

use symbolica::atom::{Atom, AtomCore, AtomView};
use symbolica::id::Replacement;

use rustred::family::IntegralKey;

use crate::symbols::S;
use crate::utils::vakint_macros::vk_symbol;
use crate::{
    ReplacementRules, VakintSettings, VakintTerm, get_integer_from_atom, get_prop_with_id,
};

use super::RustRedEvaluationError;

#[derive(Clone, Copy)]
pub(super) enum ScalarNumeratorForm {
    /// Existing unweighted scalar services consume expanded polynomial syntax.
    ExpandedPolynomial,
    /// Weighted lowering substitutes the complete denominator basis first.
    PreserveProducts,
}

/// Scalar input in the retained parent's loop basis, with explicit artifact slots.
pub(super) struct MatchedScalarInput {
    pub target: IntegralKey,
    pub mass_squared: Atom,
    pub numerator: Atom,
}

impl MatchedScalarInput {
    pub fn prepare(
        numerator: AtomView,
        matched: &ReplacementRules,
        parent_momenta: &[Atom],
        physical_slots: &[usize],
        index_count: usize,
        settings: &VakintSettings,
        form: ScalarNumeratorForm,
    ) -> Result<Self, RustRedEvaluationError> {
        let invalid = |detail: &str| RustRedEvaluationError::InvalidMatchedFamily {
            detail: detail.into(),
        };
        let integral = matched.canonical_topology.get_integral();
        integral.validate_parent_routing(parent_momenta)?;
        let slots = physical_slots
            .iter()
            .copied()
            .collect::<std::collections::HashSet<_>>();
        if physical_slots.len() != integral.n_props
            || slots.len() != physical_slots.len()
            || physical_slots.iter().any(|slot| *slot >= index_count)
        {
            return Err(invalid(
                "physical-to-family slot map has invalid arity, duplicate or out-of-range slots",
            ));
        }
        let loop_momenta = (1..=integral.n_loops)
            .map(|axis| S.k.call(axis))
            .collect::<Vec<_>>();
        let routed = route_to_parent(
            numerator,
            &loop_momenta,
            integral
                .parent_routing
                .as_ref()
                .map(|(_, coordinates)| coordinates.as_ref()),
            settings,
            form,
        )?;
        let numerator = routed;
        let expression = integral
            .canonical_expression
            .as_ref()
            .ok_or_else(|| invalid("matched canonical expression is absent"))?;
        let mut powers = vec![0; index_count];
        let mut mass = None;
        for (physical, slot) in physical_slots.iter().enumerate() {
            let Some(properties) = get_prop_with_id(expression.as_view(), physical + 1) else {
                continue;
            };
            powers[*slot] = properties
                .get(&vk_symbol!("pow_"))
                .and_then(|power| get_integer_from_atom(power.as_view()))
                .ok_or_else(|| invalid("noninteger matched propagator power"))?;
            let current = properties
                .get(&vk_symbol!("mUVsq_"))
                .ok_or_else(|| invalid("matched propagator mass is absent"))?;
            if current.is_zero()
                || current.contains_symbol(S.k)
                || mass.as_ref().is_some_and(|mass| mass != current)
            {
                return Err(invalid(
                    "saved application requires one nonzero common mass",
                ));
            }
            mass.get_or_insert_with(|| current.clone());
        }
        Ok(Self {
            target: IntegralKey::try_new(powers).map_err(|error| {
                RustRedEvaluationError::IntegralKey {
                    detail: error.to_string(),
                }
            })?,
            mass_squared: mass.ok_or_else(|| invalid("matched scalar input has no mass"))?,
            numerator,
        })
    }
}

pub(super) fn route_to_parent(
    numerator: AtomView,
    loop_momenta: &[Atom],
    coordinates: Option<&[Atom]>,
    settings: &VakintSettings,
    form: ScalarNumeratorForm,
) -> Result<Atom, RustRedEvaluationError> {
    let scalar = VakintTerm::scalar_product_form(settings, numerator).map_err(|error| {
        RustRedEvaluationError::ScalarNumerator {
            detail: error.to_string(),
        }
    })?;
    let Some(coordinates) = coordinates else {
        return Ok(scalar);
    };
    // Restore native closed-contraction scope before expanding. Opening a
    // scalar power first can incorrectly pair repeated dummy indices, changing
    // (k1.k2)^2 to k1²*k2². Symmetric linear dot heads permit direct simultaneous
    // loop substitution without opening indices at all.
    let routing = loop_momenta
        .iter()
        .zip(coordinates)
        .map(|(source, target)| Replacement::new(source.clone(), target.clone()))
        .collect::<Vec<_>>();
    let routed = scalar.replace_multiple(&routing);
    Ok(match form {
        ScalarNumeratorForm::ExpandedPolynomial => routed.expand(),
        ScalarNumeratorForm::PreserveProducts => routed,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Vakint;
    use crate::utils::vakint_macros::vk_parse;
    use rustred::algebra::CoefficientContext;
    use rustred::family::{AffineDenominator, IntegralFamily};
    use rustred::scalar_numerator::FamilyScalarNumeratorService;

    fn indexed_variants(scalar: &Atom) -> Vec<Atom> {
        let indexed = Vakint::convert_from_dot_notation(scalar.as_view());
        let renamed = indexed
            .replace(vk_parse!("dot_dummy_ind(i_)").unwrap().to_pattern())
            .with(vk_parse!("dot_dummy_ind(i_+137)").unwrap().to_pattern());
        vec![scalar.clone(), indexed, renamed]
    }

    #[test]
    fn weighted_parent_handoff_preserves_products_and_legacy_expansion_is_explicit() {
        Vakint::initialize_vakint_symbols();
        let settings = VakintSettings::default();
        let momenta = vec![S.k.call(1), S.k.call(2)];
        let coordinates = vec![&momenta[0] + &momenta[1], momenta[1].clone()];
        let input = vk_parse!("(dot(k(1),k(2))+1)^11").unwrap();
        let preserved = route_to_parent(
            input.as_view(),
            &momenta,
            Some(&coordinates),
            &settings,
            ScalarNumeratorForm::PreserveProducts,
        )
        .unwrap();
        let expanded = route_to_parent(
            input.as_view(),
            &momenta,
            Some(&coordinates),
            &settings,
            ScalarNumeratorForm::ExpandedPolynomial,
        )
        .unwrap();
        assert!(matches!(preserved.as_view(), AtomView::Pow(_)));
        assert!(matches!(expanded.as_view(), AtomView::Add(_)));
        assert_ne!(preserved, expanded);
        assert_eq!(preserved.expand(), expanded);
        // Keep the pre-existing no-routing behavior exactly: legacy callers
        // did not expand when no retained parent coordinates were supplied.
        for form in [
            ScalarNumeratorForm::ExpandedPolynomial,
            ScalarNumeratorForm::PreserveProducts,
        ] {
            assert_eq!(
                route_to_parent(input.as_view(), &momenta, None, &settings, form).unwrap(),
                input
            );
        }
    }

    #[test]
    fn closed_cross_products_keep_scope_during_parent_transport() {
        Vakint::initialize_vakint_symbols();
        let settings = VakintSettings::default();
        let momenta = vec![S.k.call(1), S.k.call(2)];
        let coordinates = vec![&momenta[0] + &momenta[1], -&momenta[1]];
        let substitution = momenta
            .iter()
            .zip(&coordinates)
            .map(|(source, target)| Replacement::new(source.clone(), target.clone()))
            .collect::<Vec<_>>();
        for expression in [
            "dot(k(1),k(2))^2",
            "dot(k(1),k(2))^3",
            "(dot(k(1),k(2))+dot(k(1),k(1)))^2*(dot(k(1),k(2))-1)^3",
            "g(71,72)*dot(k(1),k(2))^2",
        ] {
            let scalar = vk_parse!(expression).unwrap();
            let expected = scalar.replace_multiple(&substitution).expand();
            for input in indexed_variants(&scalar) {
                let actual = route_to_parent(
                    input.as_view(),
                    &momenta,
                    Some(&coordinates),
                    &settings,
                    ScalarNumeratorForm::PreserveProducts,
                )
                .unwrap();
                assert!(
                    (actual - &expected).expand().together().is_zero(),
                    "{expression}"
                );
            }
        }
        assert!(
            VakintTerm::scalar_product_form(
                &settings,
                vk_parse!("k(1,71)*k(1,72)").unwrap().as_view()
            )
            .is_err()
        );
    }

    #[test]
    fn five_loop_internal_prepass_keeps_degree_eleven_compact_with_indexed_public_setting() {
        let vakint = Vakint::new().unwrap();
        vakint
            .topologies
            .prepare_five_loop_catalogue_for_tests()
            .unwrap();
        let settings = VakintSettings {
            rustred_artifact_path: Some("/prepass-does-not-load-this-artifact".into()),
            use_dot_product_notation: false,
            evaluation_order: crate::EvaluationOrder(vec![crate::EvaluationMethod::RustRed(
                crate::RustRedEvaluationOptions {
                    substitute_masters: false,
                },
            )]),
            ..Default::default()
        };
        let input = vk_parse!("(dot(k(1),k(1))-1)^11*topo(I5L_32745(1,0,0,0,0,1,0,0,1,1,0,1,16))")
            .unwrap()
            .expand();
        let mut expression = crate::VakintExpression::try_from(input).unwrap();
        expression
            .canonicalize(&settings, &vakint.topologies, false)
            .unwrap();
        let expected = expression.0[0].numerator.clone();
        expression.0[0]
            .tensor_reduce_for_evaluation(&vakint, &settings)
            .unwrap();
        let scalar =
            VakintTerm::scalar_product_form(&settings, expression.0[0].numerator.as_view())
                .unwrap();
        assert!((scalar - expected).expand().together().is_zero());
    }

    #[test]
    fn expanded_closed_scalar_degree_eleven_uses_native_invariant_dots() {
        Vakint::initialize_vakint_symbols();
        let settings = VakintSettings {
            use_dot_product_notation: true,
            ..Default::default()
        };
        let scalar = vk_parse!("(dot(k(1),k(1))-1)^11").unwrap().expand();
        for (input, expected) in [
            (scalar.clone(), scalar.clone()),
            (
                &scalar * vk_parse!("k(2,71)*k(2,72)").unwrap(),
                &scalar * vk_parse!("g(71,72)*dot(k(2),k(2))/(4-2*ε)").unwrap(),
            ),
            (
                &scalar * vk_parse!("dot(k(2),p(1))^2").unwrap(),
                &scalar * vk_parse!("dot(k(2),k(2))*dot(p(1),p(1))/(4-2*ε)").unwrap(),
            ),
        ] {
            let mut term = VakintTerm {
                numerator: input,
                integral: Atom::one(),
                vectors: Vec::new(),
            };
            term.tensor_reduce_feynkit(&settings, false).unwrap();
            assert!((term.numerator - expected).expand().together().is_zero());
        }
        assert!(
            (VakintTerm::scalar_product_form(&settings, scalar.as_view()).unwrap() - scalar)
                .expand()
                .together()
                .is_zero()
        );
    }

    #[test]
    fn matcher_scalar_substitution_seals_scope_before_changing_loop_basis() {
        Vakint::initialize_vakint_symbols();
        let settings = VakintSettings::default();
        let mut matched = ReplacementRules::default();
        matched.numerator_substitutions = [
            ("k(1,idx___)", "k(1,idx___)-k(3,idx___)"),
            ("k(2,idx___)", "k(2,idx___)-k(3,idx___)"),
            ("k(3,idx___)", "-k(3,idx___)"),
        ]
        .into_iter()
        .map(|(source, target)| (vk_parse!(source).unwrap(), vk_parse!(target).unwrap()))
        .collect();
        for expression in [
            "(dot(k(1),k(1))-1)^5",
            "(dot(k(1),k(1))-1)^11",
            "dot(k(1),k(2))^2",
            "(dot(k(1),k(1))-2*dot(k(1),k(3))+dot(k(3),k(3))-1)^5",
        ] {
            let scalar = vk_parse!(expression).unwrap();
            let mut expected = VakintTerm {
                numerator: scalar.clone(),
                integral: Atom::one(),
                vectors: Vec::new(),
            };
            expected
                .apply_numerator_replacement_rules(&matched, &settings)
                .unwrap();
            for numerator in indexed_variants(&scalar) {
                let mut actual = VakintTerm {
                    numerator,
                    integral: Atom::one(),
                    vectors: Vec::new(),
                };
                actual
                    .apply_scalar_numerator_replacement_rules(&matched, &settings)
                    .unwrap();
                assert!(
                    (actual.numerator - &expected.numerator)
                        .expand()
                        .together()
                        .is_zero(),
                    "{expression}"
                );
            }
        }
    }

    #[test]
    fn scoped_denominator_powers_lower_to_one_parent_key() {
        Vakint::initialize_vakint_symbols();
        let settings = VakintSettings::default();
        let context = CoefficientContext::try_new(["d"]).unwrap();
        let family = IntegralFamily::new(
            "scalar-scope-fixture",
            vec!["q1".into(), "q2".into()],
            vec![],
            context.clone(),
            context.parameter("d").unwrap(),
            [[1, 0, 0], [0, 0, 1], [1, -2, 1]]
                .into_iter()
                .map(|row| {
                    AffineDenominator::new(
                        context.integer(-1),
                        row.into_iter()
                            .map(|value| context.integer(value))
                            .collect(),
                    )
                })
                .collect(),
            vec![],
            vec![context.zero(); 3],
        )
        .unwrap();
        let momenta = vec![S.k.call(1), S.k.call(2)];
        let coordinates = vec![&momenta[0] + &momenta[1], momenta[1].clone()];
        let lowerer = FamilyScalarNumeratorService::try_new(
            &family,
            S.dot,
            momenta.clone(),
            Default::default(),
        )
        .unwrap();
        let denominator = vk_parse!("dot(k(1),k(1))-2*dot(k(1),k(2))+dot(k(2),k(2))-1").unwrap();
        for degree in [2, 5, 11] {
            for input in indexed_variants(&denominator.clone().pow(degree)) {
                let routed = route_to_parent(
                    input.as_view(),
                    &momenta,
                    Some(&coordinates),
                    &settings,
                    ScalarNumeratorForm::PreserveProducts,
                )
                .unwrap();
                let lowered = lowerer
                    .lower_weighted(
                        &routed,
                        &IntegralKey::try_new([degree + 1, 1, 1]).unwrap(),
                        &Atom::one(),
                    )
                    .unwrap();
                assert_eq!(lowered.terms().len(), 1, "degree {degree}");
                let term = &lowered.terms()[0];
                assert_eq!(term.integral().powers(), &[1, 1, 1]);
                assert_eq!(term.coefficient(), &context.one());
                assert_eq!(term.scalar_spectator(), &Atom::one());
            }
        }
    }
}
