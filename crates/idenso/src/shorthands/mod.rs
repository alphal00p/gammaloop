use spenso::structure::{abstract_index::AbstractIndex, partial::PartialStructure};
use spenso::{
    network::parsing::{ParseSettings, SchoonschipExpansionMode, ShorthandParsing},
    structure::slot::{AbsInd, DummyAind, ParseableAind},
};
use std::sync::Arc;
use symbolica::{
    atom::{Atom, AtomView},
    id::AliasedAtom,
};

use crate::{
    NetworkToolingError,
    tensor::{
        SymbolicNetExt, SymbolicNetParse, SymbolicTensor, aliases::AliasInterfaces,
        inference::TensorInferenceError,
    },
};

pub(crate) mod bracket;
pub mod chain;
pub mod metric;
pub mod schoonschip;

impl SymbolicTensor<PartialStructure> {
    /// Render compact metric products as dots without contracting indexed factors.
    pub fn to_dots(&self) -> Result<Self, TensorInferenceError> {
        self.with_rewritten_expression(schoonschip::DotNormalizer::metric_shorthand_to_dot(
            self.expression().as_view(),
        ))
    }

    /// Open dots into symbolic indexed contractions, retaining other shorthands.
    /// This does not expand a numerator or evaluate finite tensor components.
    pub fn undo_dots(&self) -> Result<Self, TensorInferenceError> {
        let expression = self.expression().undo_dots::<AbstractIndex>()?;
        self.with_rewritten_expression(expression)
    }
}

impl SymbolicTensor<AliasInterfaces, AliasedAtom> {
    /// Render compact metric products in the root and reachable definitions.
    pub fn to_dots(self: &Arc<Self>) -> Result<Arc<Self>, TensorInferenceError> {
        self.map_domains(|domain, _| Ok((domain.to_dots()?, Vec::new())))
    }

    /// Open dots in the root and reachable definitions without resolving aliases.
    pub fn undo_dots(self: &Arc<Self>) -> Result<Arc<Self>, TensorInferenceError> {
        self.map_domains(|domain, _| Ok((domain.undo_dots()?, Vec::new())))
    }
}

pub trait UndoShorthands {
    fn undo_schoonschip<Aind: AbsInd + DummyAind + ParseableAind>(
        &self,
    ) -> Result<Atom, NetworkToolingError>;

    fn undo_all<Aind: AbsInd + DummyAind + ParseableAind>(
        &self,
    ) -> Result<Atom, NetworkToolingError>;
    fn undo_dots<Aind: AbsInd + DummyAind + ParseableAind>(
        &self,
    ) -> Result<Atom, NetworkToolingError>;
    fn undo_chain<Aind: AbsInd + DummyAind + ParseableAind>(
        &self,
    ) -> Result<Atom, NetworkToolingError>;
    fn undo_trace<Aind: AbsInd + DummyAind + ParseableAind>(
        &self,
    ) -> Result<Atom, NetworkToolingError>;
}

impl UndoShorthands for Atom {
    fn undo_schoonschip<Aind: AbsInd + DummyAind + ParseableAind>(
        &self,
    ) -> Result<Atom, NetworkToolingError> {
        self.as_view().undo_schoonschip::<Aind>()
    }

    fn undo_all<Aind: AbsInd + DummyAind + ParseableAind>(
        &self,
    ) -> Result<Atom, NetworkToolingError> {
        self.as_view().undo_all::<Aind>()
    }

    fn undo_dots<Aind: AbsInd + DummyAind + ParseableAind>(
        &self,
    ) -> Result<Atom, NetworkToolingError> {
        self.as_view().undo_dots::<Aind>()
    }

    fn undo_chain<Aind: AbsInd + DummyAind + ParseableAind>(
        &self,
    ) -> Result<Atom, NetworkToolingError> {
        self.as_view().undo_chain::<Aind>()
    }

    fn undo_trace<Aind: AbsInd + DummyAind + ParseableAind>(
        &self,
    ) -> Result<Atom, NetworkToolingError> {
        self.as_view().undo_trace::<Aind>()
    }
}

impl<'a> UndoShorthands for AtomView<'a> {
    fn undo_all<Aind: AbsInd + DummyAind + ParseableAind>(
        &self,
    ) -> Result<Atom, NetworkToolingError> {
        let net = self
            .parse_to_symbolic_net::<Aind>(&ParseSettings {
                shorthand_parsing: ShorthandParsing::Expand {
                    schoonschip: SchoonschipExpansionMode::full(),
                    trace: true,
                    chain: true,
                },
                ..Default::default()
            })
            .map_err(|error| NetworkToolingError::Parse {
                reason: error.to_string(),
            })?;

        net.simple_execute::<()>()
    }
    fn undo_chain<Aind: AbsInd + DummyAind + ParseableAind>(
        &self,
    ) -> Result<Atom, NetworkToolingError> {
        let net = self
            .parse_to_symbolic_net::<Aind>(&ParseSettings {
                shorthand_parsing: ShorthandParsing::Expand {
                    schoonschip: SchoonschipExpansionMode::none(),
                    trace: false,
                    chain: true,
                },
                ..Default::default()
            })
            .map_err(|error| NetworkToolingError::Parse {
                reason: error.to_string(),
            })?;

        net.simple_execute::<()>()
    }

    fn undo_dots<Aind: AbsInd + DummyAind + ParseableAind>(
        &self,
    ) -> Result<Atom, NetworkToolingError> {
        let net = self
            .parse_to_symbolic_net::<Aind>(&ParseSettings {
                shorthand_parsing: ShorthandParsing::Expand {
                    schoonschip: SchoonschipExpansionMode {
                        inner_products: true,
                        expand_inside_chains: false,
                        expand_schoonship: false,
                    },
                    trace: false,
                    chain: false,
                },
                ..Default::default()
            })
            .map_err(|error| NetworkToolingError::Parse {
                reason: error.to_string(),
            })?;

        net.simple_execute::<()>()
    }

    fn undo_schoonschip<Aind: AbsInd + DummyAind + ParseableAind>(
        &self,
    ) -> Result<Atom, NetworkToolingError> {
        let net = self
            .parse_to_symbolic_net::<Aind>(&ParseSettings {
                shorthand_parsing: ShorthandParsing::Expand {
                    schoonschip: SchoonschipExpansionMode {
                        inner_products: false,
                        expand_inside_chains: true,
                        expand_schoonship: true,
                    },
                    trace: false,
                    chain: false,
                },
                ..Default::default()
            })
            .map_err(|error| NetworkToolingError::Parse {
                reason: error.to_string(),
            })?;

        net.simple_execute::<()>()
    }

    fn undo_trace<Aind: AbsInd + DummyAind + ParseableAind>(
        &self,
    ) -> Result<Atom, NetworkToolingError> {
        let net = self
            .parse_to_symbolic_net::<Aind>(&ParseSettings {
                shorthand_parsing: ShorthandParsing::Expand {
                    schoonschip: SchoonschipExpansionMode::none(),
                    trace: true,
                    chain: false,
                },
                ..Default::default()
            })
            .map_err(|error| NetworkToolingError::Parse {
                reason: error.to_string(),
            })?;

        net.simple_execute::<()>()
    }
}

#[cfg(test)]
mod tests {
    use spenso::{
        chain, chain_factor, dot, mink, network::tags::SPENSO_TAG,
        structure::abstract_index::AbstractIndex, vector,
    };
    use symbolica::{
        atom::{Atom, AtomCore, FunctionBuilder},
        symbol,
    };
    use symbolica_utils::AtomPrintExt;

    use crate::shorthands::{bracket::BracketNormalizer, schoonschip::Schoonschip};
    use crate::test_support::test_initialize;

    use super::*;

    #[test]
    fn undo_dots_positive_integer_powers_roundtrip() {
        let _ = test_initialize();
        let product = dot!(
            vector!(powered_dot_p, mink!(4)),
            vector!(powered_dot_q, mink!(4))
        );

        for exponent in [2, 3] {
            let expression = product.clone().pow(Atom::num(exponent));
            let expanded = expression.undo_dots::<AbstractIndex>().unwrap();

            assert_eq!(
                crate::shorthands::schoonschip::DotNormalizer::metric_shorthand_to_dot(
                    (crate::test_support::contracted_atom(
                        BracketNormalizer::normalize(expanded.as_view())
                            .normalize_dots()
                            .as_view()
                    )
                    .unwrap())
                    .as_view()
                ),
                expression,
                "independent contractions must survive power {exponent}"
            );
        }
    }

    #[test]
    fn undo_dots_powers_preserve_factorized_sums() {
        let _ = test_initialize();
        let product = dot!(vector!(sum_p, mink!(4)), vector!(sum_q, mink!(4)));
        let coefficient = (Atom::var(symbol!("dot_coefficient")) + Atom::num(1)).pow(5);
        let expression = coefficient * (product + Atom::num(1)).pow(2);
        let expanded = expression.undo_dots::<AbstractIndex>().unwrap();

        let AtomView::Mul(factors) = expanded.as_view() else {
            panic!("the product of sums must stay factorized: {expanded}");
        };
        assert_eq!(
            factors
                .iter()
                .filter(|factor| matches!(factor, AtomView::Add(_)))
                .count(),
            2,
        );
        assert_eq!(
            crate::shorthands::schoonschip::DotNormalizer::metric_shorthand_to_dot(
                (crate::test_support::contracted_atom(
                    BracketNormalizer::normalize(expanded.as_view())
                        .normalize_dots()
                        .as_view()
                )
                .unwrap())
                .as_view()
            ),
            expression,
        );
    }

    #[test]
    fn undo_dots_powers_reserve_existing_indices() {
        let _ = test_initialize();
        let spectator = vector!(
            reserved_dot_r,
            mink!(4, AbstractIndex::new_dummy_at(1_000_000).to_atom())
        );
        let product = dot!(
            vector!(reserved_dot_p, mink!(4)),
            vector!(reserved_dot_q, mink!(4))
        );
        let expression = spectator * product.pow(2);
        let expanded = expression.undo_dots::<AbstractIndex>().unwrap();

        assert_eq!(
            crate::shorthands::schoonschip::DotNormalizer::metric_shorthand_to_dot(
                (crate::test_support::contracted_atom(
                    BracketNormalizer::normalize(expanded.as_view())
                        .normalize_dots()
                        .as_view()
                )
                .unwrap())
                .as_view()
            ),
            expression,
        );
    }

    #[test]
    fn undo_dots_powers_preserve_self_dual_external_ports() {
        let _ = test_initialize();
        let product = dot!(vector!(open_dot_p, mink!(4)), vector!(open_dot_q, mink!(4)));
        let external = vector!(open_dot_r, mink!(4, mu));
        let norm = dot!(vector!(open_dot_r, mink!(4)), vector!(open_dot_r, mink!(4)));

        for exponent in [2, 3] {
            // Keep the external port and lowered dot in the same powered parser group.
            let expression = spenso::bracket!(&external, &product).pow(exponent);
            let expanded = expression.undo_dots::<AbstractIndex>().unwrap();
            let expected =
                product.pow(exponent) * norm.pow(exponent / 2) * external.pow(exponent % 2);
            assert_eq!(
                crate::shorthands::schoonschip::DotNormalizer::metric_shorthand_to_dot(
                    (crate::test_support::contracted_atom(
                        BracketNormalizer::normalize(expanded.as_view())
                            .normalize_dots()
                            .as_view()
                    )
                    .unwrap())
                    .as_view()
                ),
                expected,
            );
        }
    }

    #[test]
    fn undo_dots_negative_integer_powers_roundtrip() {
        let _ = test_initialize();
        let p = vector!(inverse_dot_p, mink!(4));
        let q = vector!(inverse_dot_q, mink!(4));

        for product in [dot!(p.clone(), q), dot!(p.clone(), p)] {
            for exponent in [-1, -2, -3] {
                let expression = product.clone().pow(Atom::num(exponent));
                let expanded = expression.undo_dots::<AbstractIndex>().unwrap();

                assert!(expanded.contains_symbol(SPENSO_TAG.bracket));
                assert_eq!(expanded.undo_dots::<AbstractIndex>().unwrap(), expanded);
                let normalized = BracketNormalizer::normalize(expanded.as_view()).normalize_dots();
                let planned = crate::shorthands::schoonschip::SlotContraction::new()
                    .contract_factorized(normalized.as_view(), None, true)
                    .map(|plan| {
                        (
                            plan.root,
                            plan.status
                                == crate::shorthands::schoonschip::ContractionStatus::Complete,
                            plan.aliases,
                        )
                    });
                assert_eq!(
                    crate::shorthands::schoonschip::DotNormalizer::metric_shorthand_to_dot(
                        (crate::test_support::contracted_atom(normalized.as_view()).unwrap())
                            .as_view()
                    ),
                    expression,
                    "the contracted denominator must survive power {exponent}; undone={expanded}; normalized={normalized}; plan={planned:?}"
                );
            }
        }
    }

    #[test]
    fn undo_schoonschip_across_chain() {
        let _ = test_initialize();
        let expr = chain!(
            mink!(4, i),
            mink!(4, j),
            chain_factor!(
                schoonschip_only_factor,
                in,
                out,
                (vector!(compact_p, mink!(4)))
            )
        );

        insta::assert_snapshot!( expr.undo_schoonschip::<AbstractIndex>().unwrap().to_bare_ordered_string(), @"chain(mink(4,i),mink(4,j),schoonschip_only_factor(in,out,mink(4,d_1000000)))*compact_p(mink(4,d_1000000))");
    }

    #[test]
    fn network_tooling_malformed_shorthand_returns_a_parse_error() {
        let _ = test_initialize();
        let expression = FunctionBuilder::new(SPENSO_TAG.dot)
            .add_arg(Atom::var(symbol!("malformed_dot_operand")))
            .finish();

        assert!(matches!(
            expression
                .undo_dots::<AbstractIndex>()
                .expect_err("malformed dot notation should return an error"),
            NetworkToolingError::Parse { reason } if reason.contains("Invalid dot function")
        ));
    }
}

#[cfg(test)]
mod typed_dot_tests {
    use super::*;
    use crate::tensor::SymbolicTensor;
    use spenso::{
        dot, g, mink,
        network::tags::SPENSO_TAG,
        structure::partial::{PartialStructure, PartialStructureExt},
        tensor, vector,
    };
    use std::sync::{
        Arc,
        atomic::{AtomicUsize, Ordering},
    };
    use symbolica::{atom::AtomCore, function, symbol};

    #[test]
    fn typed_dot_formatting_does_not_contract_explicit_factors() {
        crate::test_support::test_initialize();
        let p = vector!(typed_dot_p, mink!(4));
        let q = vector!(typed_dot_q, mink!(4));
        let compact = SymbolicTensor::infer(g!(&p, &q)).unwrap();
        let expected = dot!(&p, &q);
        let formatted = compact.to_dots().unwrap();
        assert_eq!(formatted.expression(), &expected);
        assert_eq!(formatted.to_dots().unwrap(), formatted);

        let explicit =
            vector!(typed_dot_p, mink!(4, 97001)) * vector!(typed_dot_q, mink!(4, 97001));
        let indexed = SymbolicTensor::infer(explicit.clone()).unwrap();
        assert_eq!(indexed.to_dots().unwrap().expression(), &explicit);
        assert_eq!(indexed.expression(), &explicit);
        let contracted = indexed
            .contract(Default::default())
            .unwrap()
            .resolved()
            .unwrap();
        assert_eq!(contracted.to_dots().unwrap().expression(), &expected);
    }

    #[test]
    fn typed_dot_opening_preserves_logical_order_zero_and_factorization() {
        crate::test_support::test_initialize();
        let product = dot!(
            vector!(typed_open_p, mink!(4)),
            vector!(typed_open_q, mink!(4))
        );
        let scalar = (Atom::var(symbol!("typed_dot_coefficient")) + Atom::one()).pow(5);
        let open = tensor!(typed_dot_open, mink!(4, 97011), mink!(4, 97012));
        let expression = &open * &scalar * (product + Atom::one()).pow(2);
        let inferred = SymbolicTensor::infer(expression.clone()).unwrap();
        let mut slots = inferred.structure().logical_slots();
        slots.reverse();
        let source = SymbolicTensor::checked_parts(
            expression.clone(),
            PartialStructure::from_logical_slots(slots),
        )
        .unwrap();
        let opened = source.undo_dots().unwrap();
        assert_eq!(opened.structure(), source.structure());
        assert!(!opened.expression().contains_symbol(SPENSO_TAG.dot));
        assert!(
            matches!(opened.expression().as_view(), symbolica::atom::AtomView::Mul(factors)
            if factors.iter().any(|factor| factor == scalar.as_view()))
        );
        assert_eq!(source.expression(), &expression);
        assert_eq!(opened.undo_dots().unwrap(), opened);
        let zero = source.with_rewritten_expression(Atom::Zero).unwrap();
        assert_eq!(zero.to_dots().unwrap(), zero);
        assert_eq!(zero.undo_dots().unwrap(), zero);
    }

    #[test]
    fn typed_dot_formatting_checks_callback_rank_loss_once() {
        crate::test_support::test_initialize();
        let calls = Arc::new(AtomicUsize::new(0));
        let seen = Arc::clone(&calls);
        let dot_head = SPENSO_TAG.dot;
        let head = spenso::tensor_symbol!(
            "typed_dot_rank_loss",
            norm = move |value, out| {
                if value.contains_symbol(dot_head) {
                    seen.fetch_add(1, Ordering::Relaxed);
                    **out = Atom::one();
                }
            }
        );
        let metadata = symbol!("typed_dot_metadata"; Scalar);
        let expression = function!(
            head,
            function!(
                metadata,
                g!(
                    vector!(typed_callback_p, mink!(4)),
                    vector!(typed_callback_q, mink!(4))
                )
            ),
            mink!(4, 97021)
        );
        let source = SymbolicTensor::infer(expression).unwrap();
        calls.store(0, Ordering::Relaxed);
        assert!(source.to_dots().is_err());
        assert_eq!(calls.load(Ordering::Relaxed), 1);
    }

    #[test]
    fn typed_dot_aliases_keep_registry_and_reuse_identity() {
        crate::test_support::test_initialize();
        let body = SymbolicTensor::infer(g!(
            vector!(typed_alias_p, mink!(4)),
            vector!(typed_alias_q, mink!(4))
        ))
        .unwrap();
        let disconnected = SymbolicTensor::infer(g!(
            vector!(typed_unused_p, mink!(4)),
            vector!(typed_unused_q, mink!(4))
        ))
        .unwrap();
        let handle = body.alias_handle().unwrap();
        let unused = disconnected.alias_handle().unwrap();
        let source = Arc::new(
            handle
                .clone()
                .with_aliases([
                    (handle.clone(), body.clone()),
                    (unused.clone(), disconnected.clone()),
                ])
                .unwrap(),
        );
        let formatted = source.to_dots().unwrap();
        assert_eq!(formatted.root(), source.root());
        assert_eq!(formatted.aliases().unwrap().len(), 2);
        let definitions = formatted.aliases().unwrap();
        assert!(definitions.contains(&(unused, disconnected)));
        assert!(Arc::ptr_eq(&formatted, &formatted.to_dots().unwrap()));
        let opened = formatted.undo_dots().unwrap();
        assert_eq!(opened.root(), source.root());
        assert_eq!(opened.aliases().unwrap().len(), 2);
        assert_eq!(
            source
                .aliases()
                .unwrap()
                .iter()
                .find(|(h, _)| h == &handle)
                .unwrap()
                .1,
            body
        );
    }
}
