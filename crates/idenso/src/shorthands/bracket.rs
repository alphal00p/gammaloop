use spenso::{
    bracket,
    network::{parsing::AtomStructureExt, tags::SPENSO_TAG},
    structure::{
        OrderedStructure,
        abstract_index::AbstractIndex,
        representation::{LibraryRep, Representation},
        slot::Slot,
    },
};
use symbolica::atom::{Atom, AtomCore, AtomView};

/// Expose a bracket's factors to ordinary product identities without erasing
/// the scope of an unresolved indexed product (in particular under a power).
pub(crate) struct BracketNormalizer;

impl BracketNormalizer {
    pub(crate) fn normalize(expression: AtomView<'_>) -> Atom {
        // Every transformation below opens a bracket at this node or in one
        // of its factors. Ordinary tensor expressions need no rebuilding.
        if !expression.contains_symbol(SPENSO_TAG.bracket) {
            return expression.to_owned();
        }

        expression.replace_map_bottom_up(|node, _context, out| {
            if let AtomView::Mul(factors) = node {
                if factors
                    .iter()
                    .any(|factor| Self::open_payload(factor).is_some())
                {
                    let product = factors.iter().fold(Atom::one(), |product, factor| {
                        product * Self::open_payload(factor).unwrap_or(factor)
                    });
                    **out = Self::scoped_product(product);
                }
                return;
            }
            let AtomView::Fun(function) = node else {
                return;
            };
            if function.get_symbol() != SPENSO_TAG.bracket || function.get_nargs() == 0 {
                return;
            }
            if function.iter().any(|factor| factor.is_zero()) {
                **out = Atom::Zero;
                return;
            }
            // Unnamed ports still depend on operand order. Only explicitly
            // wired operands have an equivalent commutative Einstein product.
            if function.iter().any(Self::has_implicit_ports) {
                return;
            }
            let product = function.iter().fold(Atom::one(), |product, factor| {
                product * Self::open_payload(factor).unwrap_or(factor)
            });
            **out = Self::scoped_product(product);
        })
    }

    fn open_payload(expression: AtomView<'_>) -> Option<AtomView<'_>> {
        let AtomView::Fun(function) = expression else {
            return None;
        };
        if function.get_symbol() != SPENSO_TAG.bracket || function.get_nargs() != 1 {
            return None;
        }
        let payload = function.iter().next()?;
        if Self::has_implicit_ports(payload) {
            return None;
        }
        // Fast inference returns EmptyStructure for a scalar. Its contracted
        // indices belong to a closed scope and must not join the outer product.
        payload.infer_structure::<OrderedStructure>().ok()?;
        Some(payload)
    }

    fn scoped_product(product: Atom) -> Atom {
        // A single tensor function is already atomic. Scalar coefficients can
        // leave the bracket, but two unresolved indexed factors must stay in
        // one scope: (sum_i A_i B_i)^n is not sum_i A_i^n B_i^n.
        let mut coefficient = Atom::one();
        let mut tensors = Vec::new();
        let mut classify = |factor: AtomView<'_>| {
            if Self::contains_indices(factor) {
                tensors.push(factor.to_owned());
            } else {
                coefficient *= factor;
            }
        };
        if let AtomView::Mul(factors) = product.as_view() {
            for factor in factors {
                classify(factor);
            }
        } else {
            classify(product.as_view());
        }
        let payload = tensors
            .iter()
            .fold(Atom::one(), |product, factor| product * factor);
        if tensors.len() > 1 || matches!(payload.as_view(), AtomView::Add(_) | AtomView::Pow(_)) {
            coefficient * bracket!(payload)
        } else {
            coefficient * payload
        }
    }

    fn has_implicit_ports(expression: AtomView<'_>) -> bool {
        match expression {
            AtomView::Fun(function) => {
                let symbol = function.get_symbol();
                if symbol.is_scalar() {
                    return false;
                }
                if symbol == SPENSO_TAG.bracket {
                    return function.iter().any(Self::has_implicit_ports);
                }
                function
                    .iter()
                    .skip(usize::from(symbol == SPENSO_TAG.trace))
                    .any(
                        |argument| match Slot::<LibraryRep, AbstractIndex>::try_from(argument) {
                            Ok(slot) => matches!(slot.aind, AbstractIndex::Open { .. }),
                            Err(_) => Representation::<LibraryRep>::try_from(argument).is_ok(),
                        },
                    )
            }
            AtomView::Mul(product) => product.iter().any(Self::has_implicit_ports),
            AtomView::Add(sum) => sum.iter().any(Self::has_implicit_ports),
            AtomView::Pow(power) => Self::has_implicit_ports(power.get_base_exp().0),
            _ => false,
        }
    }

    fn contains_indices(expression: AtomView<'_>) -> bool {
        let mut found = false;
        expression.visitor(&mut |node| {
            if let AtomView::Fun(function) = node {
                if function.get_symbol().is_scalar() {
                    return false;
                }
                found |= function.get_symbol().has_attributes_of(SPENSO_TAG.rep_)
                    && function.get_nargs() == 2;
            }
            !found
        });
        found
    }
}

#[cfg(test)]
mod tests {
    use spenso::{bracket, g, mink, slot, tensor, vector};
    use symbolica::{
        atom::{Atom, AtomCore},
        parse_lit,
    };

    use crate::{
        color::ColorSimplifySettings,
        color_f, epsilon,
        epsilon::EpsilonSimplifier,
        gamma, gamma0,
        shorthands::schoonschip::{Schoonschip, SchoonschipSettings},
        test_support::test_initialize,
    };

    #[test]
    fn encoded_open_slots_are_unresolved_but_normal_slots_are_explicit() {
        use spenso::structure::{
            abstract_index::AbstractIndex,
            representation::{LibraryRep, RepName},
            slot::IsAbstractSlot,
        };
        use symbolica::atom::FunctionBuilder;

        test_initialize();
        let rep = LibraryRep::from(spenso::structure::representation::Minkowski {}).new_rep(4);
        let head = spenso::tensor_symbol!("bracket_open_slot_predicate");
        let owner = AbstractIndex::fresh_open_owner();
        for (index, expected) in [
            (AbstractIndex::Open { owner, axis: 0 }, true),
            (AbstractIndex::Normal(931), false),
        ] {
            let value = FunctionBuilder::new(head)
                .add_arg(rep.slot::<AbstractIndex, _>(index).to_atom())
                .finish();
            assert_eq!(
                super::BracketNormalizer::has_implicit_ports(value.as_view()),
                expected
            );
        }
    }

    #[test]
    fn bracket_color_projector_contracts_across_arguments() {
        let r = test_initialize();
        let a = slot!(r.coad_da, a);
        let b = slot!(r.coad_da, b);
        let c = slot!(r.coad_da, c);
        let d = slot!(r.coad_da, d);
        let metric = g!(a, b);
        for color in [metric.clone(), color_f!(a, c, d) * color_f!(b, c, d)] {
            let expected =
                crate::tensor::SymbolicTensor::infer((&metric * &color).as_atom_view().to_owned())
                    .unwrap()
                    .simplify_color(crate::color::ColorSimplifySettings::default())
                    .unwrap()
                    .resolved()
                    .unwrap()
                    .into_expression();
            assert_eq!(
                crate::tensor::SymbolicTensor::infer(
                    (bracket!(&metric, &color)).as_atom_view().to_owned()
                )
                .unwrap()
                .simplify_color(crate::color::ColorSimplifySettings::default())
                .unwrap()
                .resolved()
                .unwrap()
                .into_expression(),
                expected
            );
            let settings = ColorSimplifySettings {
                simplify_non_color: false,
                ..Default::default()
            };
            assert_eq!(
                crate::tensor::SymbolicTensor::infer(
                    (bracket!(&metric, &color)).as_atom_view().to_owned()
                )
                .unwrap()
                .simplify_color(settings)
                .unwrap()
                .resolved()
                .unwrap()
                .into_expression(),
                expected
            );
        }
    }

    #[test]
    fn bracket_color_zero_annihilates_the_whole_product() {
        let r = test_initialize();
        let a = slot!(r.coad_da, a);
        let b = slot!(r.coad_da, b);
        let c = slot!(r.coad_da, c);
        let d = slot!(r.coad_da, d);
        let e = slot!(r.coad_da, e);
        let color = color_f!(a, b, c) * color_f!(c, d, e) * g!(d, e);
        assert!(
            crate::tensor::SymbolicTensor::infer(
                (bracket!(g!(a, b), color)).as_atom_view().to_owned()
            )
            .unwrap()
            .simplify_color(crate::color::ColorSimplifySettings::default())
            .unwrap()
            .resolved()
            .unwrap()
            .into_expression()
            .is_zero()
        );
    }

    #[test]
    fn bracket_gamma_trace_contracts_across_arguments() {
        let r = test_initialize();
        let left = gamma!(a, b, slot!(r.mink_d, mu));
        let right = gamma!(b, a, slot!(r.mink_d, nu));
        let expected = Atom::num(4) * g!(slot!(r.mink_d, mu), slot!(r.mink_d, nu));
        assert_eq!(
            crate::tensor::SymbolicTensor::infer((bracket!(left, right)).as_atom_view().to_owned())
                .unwrap()
                .simplify_gamma(crate::dirac::GammaSimplifySettings::default())
                .unwrap()
                .resolved()
                .unwrap()
                .into_expression(),
            expected
        );
    }

    #[test]
    fn bracket_gamma0_pair_contracts_across_arguments() {
        let r = test_initialize();
        let left = gamma0!(a, b);
        let right = gamma0!(b, c);
        assert_eq!(
            crate::dirac::DiracSimplifier::factor_gamma_zero(bracket!(left, right).as_view()),
            g!(slot!(r.bis4, a), slot!(r.bis4, c))
        );
    }

    #[test]
    fn bracket_epsilon_pair_contracts_across_arguments() {
        test_initialize();
        let eps = epsilon!(mink!(4, a), mink!(4, b), mink!(4, c), mink!(4, d));
        assert_eq!(bracket!(&eps, &eps).simplify_epsilon(), Atom::num(24));
    }

    #[test]
    fn bracket_nested_metrics_preserve_the_free_indices() {
        test_initialize();
        let expression = bracket!(
            bracket!(g!(mink!(4, a), mink!(4, b)), g!(mink!(4, b), mink!(4, c))),
            g!(mink!(4, c), mink!(4, d))
        );
        let expected = g!(mink!(4, a), mink!(4, d));
        assert_eq!(
            expression
                .schoonschip_with_settings(&SchoonschipSettings::default().without_rank1_tensors()),
            expected
        );
    }

    #[test]
    fn bracket_open_product_contracts_with_an_external_metric() {
        test_initialize();
        let left = tensor!(bracket_A, mink!(4, a), mink!(4, b));
        let right = tensor!(bracket_B, mink!(4, b), mink!(4, c));
        let expression = bracket!(&left, right) * g!(mink!(4, c), mink!(4, d));
        let expected = bracket!(left * tensor!(bracket_B, mink!(4, b), mink!(4, d)));
        assert_eq!(
            expression
                .schoonschip_with_settings(&SchoonschipSettings::default().without_rank1_tensors()),
            expected
        );
    }

    #[test]
    fn bracket_simplification_preserves_scalar_factorization() {
        test_initialize();
        let coefficient = parse_lit!((x + y) ^ 12 * (z + w) ^ 12);
        let metric = g!(mink!(4, a), mink!(4, b));
        let expression = bracket!(&coefficient * &metric, &metric);
        assert_eq!(
            expression
                .schoonschip_with_settings(&SchoonschipSettings::default().without_rank1_tensors()),
            Atom::num(4) * coefficient
        );
    }

    #[test]
    fn bracket_unnamed_ports_keep_their_order() {
        use spenso::structure::{
            abstract_index::AbstractIndex, representation::RepName, slot::IsAbstractSlot,
        };
        use std::collections::HashMap;
        use symbolica::atom::{AtomView, FunctionBuilder};

        test_initialize();
        let left = vector!(bracket_z, mink!(4));
        let right = vector!(bracket_a, mink!(4));
        let expression = bracket!(&left, &right);
        assert_eq!(
            expression
                .schoonschip_with_settings(&SchoonschipSettings::default().without_rank1_tensors()),
            expression
        );
        let source = crate::tensor::SymbolicTensor::infer(expression).unwrap();
        let result = source
            .simplify_color(crate::color::ColorSimplifySettings::default())
            .unwrap()
            .resolved()
            .unwrap();
        // The shared constructor can remove the bracket when normal product
        // order already preserves these two anonymous ports.
        assert!(matches!(source.expression().as_view(), AtomView::Mul(_)));
        assert_eq!(result.structure(), source.structure());
        let AtomView::Fun(left) = left.as_view() else {
            unreachable!()
        };
        let AtomView::Fun(right) = right.as_view() else {
            unreachable!()
        };
        let indices = [74881, 74882].map(AbstractIndex::Normal);
        let indexed = |head, index| {
            FunctionBuilder::new(head)
                .add_arg(
                    spenso::structure::representation::Minkowski {}
                        .new_rep(4)
                        .slot::<AbstractIndex, _>(index)
                        .to_atom(),
                )
                .finish()
        };
        let components = [
            (indexed(left.get_symbol(), indices[0]), 2),
            (indexed(left.get_symbol(), indices[1]), 3),
            (indexed(right.get_symbol(), indices[0]), 5),
            (indexed(right.get_symbol(), indices[1]), 7),
        ];
        for (positions, expected_component) in [(indices, 14), ([indices[1], indices[0]], 15)] {
            let bindings = HashMap::from([(0, positions[0]), (1, positions[1])]);
            let expected = indexed(left.get_symbol(), positions[0])
                * indexed(right.get_symbol(), positions[1]);
            for value in [&source, &result] {
                let materialized = value.materialize_interface_ports(&bindings).unwrap();
                assert_eq!(materialized, expected);
                let component = materialized.replace_map_bottom_up(|node, _, out| {
                    if let Some((_, number)) =
                        components.iter().find(|(atom, _)| atom.as_view() == node)
                    {
                        **out = Atom::num(*number);
                    }
                });
                assert_eq!(component, Atom::num(expected_component));
            }
        }
        assert!(
            bracket!(vector!(bracket_z, mink!(4)), Atom::Zero)
                .schoonschip_with_settings(&SchoonschipSettings::default().without_rank1_tensors())
                .is_zero()
        );
    }

    #[test]
    fn bracket_independent_scalar_contractions_keep_their_scopes() {
        test_initialize();
        let left = bracket!(vector!(bracket_p, mink!(4, a)) * vector!(bracket_q, mink!(4, a)));
        let right = bracket!(vector!(bracket_r, mink!(4, a)) * vector!(bracket_s, mink!(4, a)));
        let expression = &left * &right;
        assert_eq!(
            expression
                .schoonschip_with_settings(&SchoonschipSettings::default().without_rank1_tensors()),
            expression
        );
        assert_eq!(
            expression.schoonschip(),
            left.schoonschip() * right.schoonschip()
        );
    }

    #[test]
    fn bracket_inverse_power_keeps_the_contracted_sum_atomic() {
        test_initialize();
        let product = vector!(bracket_p, mink!(4, a)) * vector!(bracket_q, mink!(4, a));
        let denominator = bracket!(&product);
        for exponent in [-1, -2, 2] {
            let expression = denominator.clone().pow(exponent);
            assert_eq!(
                expression.schoonschip_with_settings(
                    &SchoonschipSettings::default().without_rank1_tensors()
                ),
                expression
            );
            assert_eq!(
                expression.schoonschip(),
                product.schoonschip().pow(exponent)
            );
        }
    }
}
