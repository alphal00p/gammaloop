use spenso::{
    network::parsing::{ParseSettings, SchoonschipExpansionMode, ShorthandParsing},
    structure::slot::{AbsInd, DummyAind, ParseableAind},
};
use symbolica::atom::{Atom, AtomView};

use crate::{
    NetworkToolingError,
    tensor::{SymbolicNetExt, SymbolicNetParse},
};

pub mod chain;
pub mod metric;
pub mod schoonschip;

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

    use crate::shorthands::{metric::MetricSimplifier, schoonschip::Schoonschip};
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
                expanded
                    .schoonschip_net::<AbstractIndex>()
                    .unwrap()
                    .metric_shorthand_to_dot(),
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
            expanded
                .schoonschip_net::<AbstractIndex>()
                .unwrap()
                .metric_shorthand_to_dot(),
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
            expanded
                .schoonschip_net::<AbstractIndex>()
                .unwrap()
                .metric_shorthand_to_dot(),
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
                expanded
                    .schoonschip_net::<AbstractIndex>()
                    .unwrap()
                    .metric_shorthand_to_dot(),
                expected,
            );
        }
    }

    #[test]
    fn undo_shorthands_preserves_radial_powers() {
        use symbolica::atom::AtomCore;

        test_initialize();
        for representation in [mink!(4), mink!(dim)] {
            let momentum = vector!(Q, 9, representation);
            let radial = spenso::dot!(&momentum, &momentum);
            for power in [-9, -5, -3, -2, -1, 0, 1, 2, 3, 5, 9] {
                let original = radial.clone().pow(power);
                // Integrated UV conversion runs these three parsing/execution
                // passes in sequence. A cube must stay a cube at every boundary.
                let schoonschip = original.undo_schoonschip::<AbstractIndex>().unwrap();
                assert_eq!(schoonschip, original, "schoonschip, power {power}");
                let chain = schoonschip.undo_chain::<AbstractIndex>().unwrap();
                assert_eq!(chain, original, "chain, power {power}");
                let trace = chain.undo_trace::<AbstractIndex>().unwrap();
                assert_eq!(trace, original, "trace, power {power}");
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
    use crate::tensor::composition::TensorCompositionError;
    use crate::test_support;
    use spenso::structure::representation::LibraryRep;
    use spenso::structure::slot::{IsAbstractSlot, SlotMatcher};
    use spenso::tensor_symbol;
    use symbolica::atom::FunctionBuilder;

    use super::*;
    use crate::tensor::SymbolicTensor;
    use spenso::network::tags::SPENSO_TAG;
    use spenso::structure::partial::{PartialStructure, PartialStructureExt};
    use spenso::{chain, chain_factor, dot, g, mink, tensor, trace, vector};
    use std::sync::Arc;
    use std::sync::atomic::{AtomicUsize, Ordering};
    use symbolica::atom::AtomCore;
    use symbolica::{function, symbol};

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
        let contracted = indexed.contract(Default::default()).unwrap();
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
        let source = SymbolicTensor::from_produced_parts(
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
    fn trace_unfolding_preserves_order_and_does_not_evaluate() {
        crate::test_support::test_initialize();
        let trace = trace!(
            mink!(4),
            chain_factor!(notation_A, in, out),
            chain_factor!(notation_B, in, out)
        );
        let source = SymbolicTensor::infer(trace).unwrap();
        let opened = source.undo_trace().unwrap();
        assert_eq!(opened.structure(), source.structure());
        assert!(opened.expression().contains_symbol(SPENSO_TAG.chain));
        assert!(!opened.expression().contains_symbol(SPENSO_TAG.trace));
        let again = source.undo_trace().unwrap();
        assert_ne!(
            opened.expression(),
            again.expression(),
            "independent unfoldings use fresh dummies"
        );
        let explicit = opened.undo_chain().unwrap();
        assert!(!explicit.expression().contains_symbol(SPENSO_TAG.chain));
        SymbolicTensor::validate_atom(explicit.expression()).unwrap();
        let joined = opened.contract(Default::default()).unwrap();
        assert_eq!(joined.expression(), source.expression());
        let empty = SymbolicTensor::infer(trace!(mink!(4)))
            .unwrap()
            .undo_trace()
            .unwrap();
        assert!(empty.expression().contains_symbol(SPENSO_TAG.chain));
        assert_ne!(empty.expression(), &Atom::num(4));
    }

    #[test]
    fn chain_unfolding_materializes_each_projector_factor() {
        use spenso::shadowing::{FactorProjector, ProjectorExpander};

        crate::test_support::test_initialize();
        let a = chain_factor!(projector_unfold_A, in, out);
        let b = chain_factor!(projector_unfold_B, in, out);
        let c = chain_factor!(projector_unfold_C, in, out);
        for projector in [
            FactorProjector::Symmetric,
            FactorProjector::Antisymmetric,
            FactorProjector::Cyclic,
        ] {
            let group = projector.apply([a.clone(), b.clone(), c.clone()]);
            for expression in [
                chain!(mink!(4, 97101), mink!(4, 97102), &group),
                trace!(mink!(4), &group),
            ] {
                let source = SymbolicTensor::infer(expression.clone()).unwrap();
                let opened = source.undo_trace().unwrap().undo_chain().unwrap();
                assert_eq!(opened.structure(), source.structure());
                assert!(!opened.expression().contains_symbol(projector.symbol()));
                assert!(!opened.expression().contains_symbol(SPENSO_TAG.chain));
                SymbolicTensor::validate_atom(opened.expression()).unwrap();
                let expected = SymbolicTensor::infer(expression.expand_projectors())
                    .unwrap()
                    .contract(Default::default())
                    .unwrap();
                let joined = opened.contract(Default::default()).unwrap();
                assert_eq!(joined.expression(), expected.expression());
            }
        }
    }

    #[test]
    fn chain_unfolding_rejects_unexpandable_projectors() {
        use spenso::shadowing::FactorProjector;

        crate::test_support::test_initialize();
        // The existing normalized projector operation bounds its factorial
        // denominator to i64. Do not silently wire an unsupported group as one
        // matrix if that operation cannot expand it.
        let factor = chain_factor!(projector_unfold_large_A, in, out);
        let group = FactorProjector::Symmetric.apply(vec![factor; 21]);
        let expression = chain!(mink!(4, 97111), mink!(4, 97112), group);
        assert!(expression.undo_chain::<AbstractIndex>().is_err());
    }

    #[test]
    fn projector_chain_unfolding_preserves_powers_spectators_and_fresh_copies() {
        use spenso::shadowing::{FactorProjector, ProjectorExpander};

        crate::test_support::test_initialize();
        let trace = trace!(
            mink!(4),
            FactorProjector::Symmetric.apply([
                chain_factor!(projector_power_A, in, out),
                chain_factor!(projector_power_B, in, out),
                chain_factor!(projector_power_C, in, out),
            ])
        );
        let coefficient = (Atom::var(symbol!("projector_spectator")) + Atom::one()).pow(40);
        let source = SymbolicTensor::infer(&coefficient * trace.pow(2)).unwrap();
        let opened = source.undo_trace().unwrap().undo_chain().unwrap();
        assert!(
            matches!(opened.expression().as_view(), AtomView::Mul(product)
            if product.iter().any(|factor| factor == coefficient.as_view()))
        );
        let expected = SymbolicTensor::infer(trace.expand_projectors().pow(2))
            .unwrap()
            .contract(Default::default())
            .unwrap();
        assert_eq!(
            opened.contract(Default::default()).unwrap().expression(),
            &(&coefficient * expected.expression())
        );

        let single = SymbolicTensor::infer(trace).unwrap();
        let first = single.undo_trace().unwrap().undo_chain().unwrap();
        let second = single.undo_trace().unwrap().undo_chain().unwrap();
        assert_ne!(first.expression(), second.expression());
        let product = first.expression() * second.expression();
        SymbolicTensor::validate_atom(&product).unwrap();
        assert_eq!(
            SymbolicTensor::infer(product)
                .unwrap()
                .contract(Default::default())
                .unwrap()
                .expression(),
            expected.expression()
        );
    }

    #[test]
    fn chain_unfolding_keeps_factor_ports_and_spectator_factorization() {
        crate::test_support::test_initialize();
        let left = mink!(4, 97071);
        let right = mink!(4, 97072);
        let coefficient = (Atom::var(symbolica::symbol!("notation_weight")) + Atom::one()).pow(40);
        let expression = &coefficient
            * chain!(
                &left,
                &right,
                chain_factor!(notation_order_A, in, out),
                chain_factor!(notation_order_B, in, out)
            );
        let source = SymbolicTensor::infer(expression).unwrap();
        let opened = source.undo_chain().unwrap();
        assert_eq!(opened.structure(), source.structure());
        assert!(
            matches!(opened.expression().as_view(), AtomView::Mul(product)
            if product.iter().any(|factor| factor == coefficient.as_view()))
        );
        let first = spenso::tensor_symbol!(notation_order_A);
        let second = spenso::tensor_symbol!(notation_order_B);
        let mut ports = std::collections::HashMap::new();
        opened.expression().visitor(&mut |value| {
            if let AtomView::Fun(function) = value
                && (function.get_symbol() == first || function.get_symbol() == second)
            {
                ports.insert(
                    function.get_symbol(),
                    function.iter().map(|a| a.to_owned()).collect::<Vec<_>>(),
                );
            }
            true
        });
        assert_eq!(ports[&first][0], left);
        assert_eq!(ports[&first][1], ports[&second][0]);
        assert_eq!(ports[&second][1], right);
    }

    #[test]
    fn independently_unfolded_dots_use_disjoint_dummy_names() {
        crate::test_support::test_initialize();
        let source = SymbolicTensor::infer(spenso::dot!(
            vector!(notation_copy_p, mink!(4)),
            vector!(notation_copy_q, mink!(4))
        ))
        .unwrap();
        let first = source.undo_dots().unwrap();
        let second = source.undo_dots().unwrap();
        assert_ne!(first.expression(), second.expression());
        let product = first.expression() * second.expression();
        SymbolicTensor::validate_atom(&product).unwrap();
        let contracted = SymbolicTensor::infer(product)
            .unwrap()
            .contract(Default::default())
            .unwrap()
            .to_dots()
            .unwrap();
        assert_eq!(contracted.expression(), &source.expression().pow(2));
    }

    #[test]
    fn dot_notation_normalizes_compact_scalars_without_index_contractions() {
        crate::test_support::test_initialize();
        let compact = spenso::g!(vector!(notation_p, mink!(4)), vector!(notation_q, mink!(4)));
        let source = SymbolicTensor::infer(compact.clone()).unwrap();
        assert_eq!(
            source.expression(),
            &compact,
            "construction does not canonicalize this spelling"
        );
        let dotted = source.to_dots().unwrap();
        assert_ne!(dotted.expression(), &compact);
        assert!(dotted.expression().contains_symbol(SPENSO_TAG.dot));
        let indexed = vector!(notation_explicit_p, mink!(4, 97081)).pow(2);
        let source = SymbolicTensor::infer(indexed.clone()).unwrap();
        assert_eq!(source.to_dots().unwrap().expression(), &indexed);
    }

    #[test]
    fn opening_dot_notation_records_its_allocated_dummy_pair() {
        test_support::test_initialize();
        let source = SymbolicTensor::infer(dot!(
            vector!(allocated_dot_p, mink!(4)),
            vector!(allocated_dot_q, mink!(4))
        ))
        .unwrap();
        let opened = source.undo_dots().unwrap();
        assert!(opened.is_scalar());
        let mut slots = SlotMatcher::default();
        let mut allocated = None;
        opened.expression().visitor(&mut |node| {
            if let Ok(slot) = slots.parse::<LibraryRep, AbstractIndex>(node) {
                allocated.get_or_insert(slot);
                return false;
            }
            true
        });
        let slot = allocated.expect("opening a dot allocates its contraction index");
        let extra = SymbolicTensor::infer(
            FunctionBuilder::new(tensor_symbol!("allocated_dot_third_use"))
                .add_arg(slot.to_atom())
                .finish(),
        )
        .unwrap();
        assert!(matches!(opened.multiply(&extra),
            Err(TensorCompositionError::InvalidExplicitMultiplicity { index, occurrences: 3 })
                if index == slot.aind()
        ));
    }
}
