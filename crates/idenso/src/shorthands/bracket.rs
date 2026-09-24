use spenso::{
    bracket,
    network::{
        parsing::{AtomStructureExt, StructureInferenceMode},
        tags::SPENSO_TAG,
    },
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
        payload
            .infer_structure::<OrderedStructure>(StructureInferenceMode::Fast)
            .ok()?;
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
                    .any(|argument| {
                        Representation::<LibraryRep>::try_from(argument).is_ok()
                            && Slot::<LibraryRep, AbstractIndex>::try_from(argument).is_err()
                    })
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
        color::{ColorSimplifier, ColorSimplifySettings},
        color_f,
        dirac::GammaSimplifier,
        epsilon,
        epsilon::EpsilonSimplifier,
        gamma, gamma0,
        shorthands::{metric::MetricSimplifier, schoonschip::Schoonschip},
        test_support::test_initialize,
    };

    #[test]
    fn bracket_color_projector_contracts_across_arguments() {
        let r = test_initialize();
        let a = slot!(r.coad_da, a);
        let b = slot!(r.coad_da, b);
        let c = slot!(r.coad_da, c);
        let d = slot!(r.coad_da, d);
        let metric = g!(a, b);
        for color in [metric.clone(), color_f!(a, c, d) * color_f!(b, c, d)] {
            let expected = (&metric * &color).simplify_color();
            assert_eq!(bracket!(&metric, &color).simplify_color(), expected);
            let settings = ColorSimplifySettings {
                simplify_non_color: false,
                ..Default::default()
            };
            assert_eq!(
                bracket!(&metric, &color).simplify_color_with(settings),
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
        assert!(bracket!(g!(a, b), color).simplify_color().is_zero());
    }

    #[test]
    fn bracket_gamma_trace_contracts_across_arguments() {
        let r = test_initialize();
        let left = gamma!(a, b, slot!(r.mink_d, mu));
        let right = gamma!(b, a, slot!(r.mink_d, nu));
        let expected = Atom::num(4) * g!(slot!(r.mink_d, mu), slot!(r.mink_d, nu));
        assert_eq!(bracket!(left, right).simplify_gamma(), expected);
    }

    #[test]
    fn bracket_gamma0_pair_contracts_across_arguments() {
        let r = test_initialize();
        let left = gamma0!(a, b);
        let right = gamma0!(b, c);
        assert_eq!(
            bracket!(left, right).simplify_gamma0(),
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
        assert_eq!(expression.simplify_metrics(), expected);
    }

    #[test]
    fn bracket_open_product_contracts_with_an_external_metric() {
        test_initialize();
        let left = tensor!(bracket_A, mink!(4, a), mink!(4, b));
        let right = tensor!(bracket_B, mink!(4, b), mink!(4, c));
        let expression = bracket!(&left, right) * g!(mink!(4, c), mink!(4, d));
        let expected = bracket!(left * tensor!(bracket_B, mink!(4, b), mink!(4, d)));
        assert_eq!(expression.simplify_metrics(), expected);
    }

    #[test]
    fn bracket_simplification_preserves_scalar_factorization() {
        test_initialize();
        let coefficient = parse_lit!((x + y) ^ 12 * (z + w) ^ 12);
        let metric = g!(mink!(4, a), mink!(4, b));
        let expression = bracket!(&coefficient * &metric, &metric);
        assert_eq!(expression.simplify_metrics(), Atom::num(4) * coefficient);
    }

    #[test]
    fn bracket_unnamed_ports_keep_their_order() {
        test_initialize();
        let left = vector!(bracket_z, mink!(4));
        let right = vector!(bracket_a, mink!(4));
        let expression = bracket!(&left, right);
        assert_eq!(expression.simplify_metrics(), expression);
        assert_eq!(expression.simplify_color(), expression);
        assert!(bracket!(left, Atom::Zero).simplify_metrics().is_zero());
    }

    #[test]
    fn bracket_independent_scalar_contractions_keep_their_scopes() {
        test_initialize();
        let left = bracket!(vector!(bracket_p, mink!(4, a)) * vector!(bracket_q, mink!(4, a)));
        let right = bracket!(vector!(bracket_r, mink!(4, a)) * vector!(bracket_s, mink!(4, a)));
        let expression = &left * &right;
        assert_eq!(expression.simplify_metrics(), expression);
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
            assert_eq!(expression.simplify_metrics(), expression);
            assert_eq!(
                expression.schoonschip(),
                product.schoonschip().pow(exponent)
            );
        }
    }
}
