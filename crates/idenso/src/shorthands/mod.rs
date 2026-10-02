use spenso::structure::{abstract_index::AbstractIndex, partial::PartialStructure};
use spenso::{
    network::parsing::{ParseSettings, SchoonschipExpansionMode, ShorthandParsing},
    structure::slot::{AbsInd, DummyAind, ParseableAind},
};
use symbolica::atom::{Atom, AtomCore, AtomView};

use crate::{
    NetworkToolingError,
    tensor::{SymbolicNetExt, SymbolicNetParse, SymbolicTensor, inference::TensorInferenceError},
};

pub(crate) mod bracket;
pub mod chain;
pub mod metric;
pub mod schoonschip;

impl SymbolicTensor<PartialStructure> {
    /// Normalize index-free scalar-product notation without contracting indexed factors.
    pub fn to_dots(&self) -> Result<Self, TensorInferenceError> {
        let observed = self.proofs.observations.get();
        if observed
            .is_some_and(|facts| facts.is_closed_scalar() && !facts.needs_dot_normalization())
        {
            return Ok(self.clone());
        }
        let facts = observed.filter(|facts| facts.is_closed_scalar()).map(|_| {
            crate::tensor::simplification::observation::DomainObservations::closed_scalar(false)
        });
        self.with_identity_result(
            schoonschip::DotNormalizer::notation(self.expression().as_view()),
            facts,
        )
    }

    /// Open dots into symbolic indexed contractions, retaining other shorthands.
    /// This does not expand a numerator or evaluate finite tensor components.
    pub fn undo_dots(&self) -> Result<Self, TensorInferenceError> {
        let expression = self.expression().undo_dots::<AbstractIndex>()?;
        self.with_identity_result(expression, None)
    }

    /// Expose chain factors with fresh compatible dummy indices and unchanged endpoints.
    pub fn undo_chain(&self) -> Result<Self, TensorInferenceError> {
        self.with_identity_result(self.expression().undo_chain::<AbstractIndex>()?, None)
    }

    /// Replace trace notation with a closed chain, without evaluating it.
    pub fn undo_trace(&self) -> Result<Self, TensorInferenceError> {
        self.with_identity_result(self.expression().undo_trace::<AbstractIndex>()?, None)
    }
}

pub trait UndoShorthands {
    #[cfg(test)]
    fn undo_schoonschip<Aind: AbsInd + DummyAind + ParseableAind>(
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
    #[cfg(test)]
    fn undo_schoonschip<Aind: AbsInd + DummyAind + ParseableAind>(
        &self,
    ) -> Result<Atom, NetworkToolingError> {
        self.as_view().undo_schoonschip::<Aind>()
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
    fn undo_chain<Aind: AbsInd + DummyAind + ParseableAind>(
        &self,
    ) -> Result<Atom, NetworkToolingError> {
        use spenso::network::{parsing::ParseState, tags::SPENSO_TAG};
        let state = ParseState::<Aind>::default();
        state.reserve_indices(*self);
        let mut error = None;
        let expression = self.replace_map(|value, context, output| {
            let AtomView::Fun(chain) = value else { return };
            if chain.get_symbol() != SPENSO_TAG.chain {
                return;
            }
            match state.materialize_indexed_chain(chain) {
                Ok(expression) => {
                    **output = if context.parent_type == Some(symbolica::atom::AtomType::Pow) {
                        symbolica::atom::FunctionBuilder::new(SPENSO_TAG.bracket)
                            .add_arg(expression)
                            .finish()
                    } else {
                        expression
                    };
                }
                Err(reason) => error = Some(reason),
            }
        });
        error.map_or(Ok(expression), |reason| {
            Err(NetworkToolingError::Parse { reason })
        })
    }

    fn undo_dots<Aind: AbsInd + DummyAind + ParseableAind>(
        &self,
    ) -> Result<Atom, NetworkToolingError> {
        use spenso::{
            network::{library::symbolic::ETS, parsing::ParseState, tags::SPENSO_TAG},
            structure::{
                representation::LibraryRep,
                slot::{IsAbstractSlot, SlotMatcher},
            },
        };
        let contains_dots = |value: AtomView<'_>| {
            value.contains_symbol(SPENSO_TAG.dot) || value.contains_symbol(ETS.metric)
        };
        if !contains_dots(*self) {
            return Ok(self.to_owned());
        }
        let state = ParseState::<Aind>::default();
        state.reserve_indices(*self);
        let unfold = |scope: AtomView<'_>| -> Result<Atom, NetworkToolingError> {
            // Preserve the parser's established per-copy scopes when notation
            // occurs inside a power. Ordinary arithmetic outside this selected
            // scope never enters the parser.
            let net = scope
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
            let mut slots = SlotMatcher::default();
            let mut written = std::collections::HashSet::new();
            scope.visitor(&mut |value| {
                if let Ok(slot) = slots.parse::<LibraryRep, Aind>(value) {
                    written.insert(slot.aind);
                    return false;
                }
                true
            });
            let expression = net.simple_execute::<()>()?;
            let mut fresh = std::collections::HashMap::new();
            // Only names allocated while opening this notation are freshened.
            // Existing bindings and external port identities remain unchanged.
            Ok(expression.replace_map(|value, _, output| {
                if let Ok(slot) = slots.parse::<LibraryRep, Aind>(value)
                    && (slot.aind.is_dummy()
                        || slot
                            .aind
                            .to_atom()
                            .as_view()
                            .get_symbol()
                            .is_some_and(|symbol| {
                                symbol.has_tag(spenso::structure::abstract_index::DUMMY_INDEX_TAG)
                            }))
                    && !written.contains(&slot.aind)
                {
                    let index = *fresh
                        .entry(slot.aind)
                        .or_insert_with(|| state.fresh_index());
                    **output = slot.reindex(index).to_atom();
                }
            }))
        };
        let mut error = None;
        let expression = self.replace_map(|value, _, output| {
            let selected = match value {
                AtomView::Pow(_) => contains_dots(value),
                AtomView::Fun(fun) => {
                    fun.get_symbol() == SPENSO_TAG.dot || fun.get_symbol() == ETS.metric
                }
                _ => false,
            };
            if selected {
                match unfold(value) {
                    Ok(expression) if expression.as_view() != value => **output = expression,
                    Ok(_) => {}
                    Err(reason) => error = Some(reason),
                }
            }
        });
        error.map_or(Ok(expression), Err)
    }

    #[cfg(test)]
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
        use spenso::{
            network::{parsing::ParseState, tags::SPENSO_TAG},
            shadowing,
            structure::{
                representation::{LibraryRep, Representation},
                slot::{DualSlotTo, IsAbstractSlot},
            },
        };
        use symbolica::atom::FunctionBuilder;
        let state = ParseState::<Aind>::default();
        state.reserve_indices(*self);
        let mut error = None;
        let expression = self.replace_map(|value, _, output| {
            let AtomView::Fun(trace) = value else { return };
            if trace.get_symbol() != SPENSO_TAG.trace {
                return;
            }
            let mut arguments = trace.iter();
            let Some(representation) = arguments.next() else {
                error = Some("trace requires a representation".to_owned());
                return;
            };
            let representation = match Representation::<LibraryRep>::try_from(representation) {
                Ok(rep) => rep,
                Err(err) => {
                    error = Some(err.to_string());
                    return;
                }
            };
            let start = representation.slot::<Aind, _>(state.fresh_index());
            let arguments = arguments.collect::<Vec<_>>();
            **output = FunctionBuilder::new(SPENSO_TAG.chain)
                .add_arg(start.to_atom())
                .add_arg(start.dual().to_atom())
                .add_args(shadowing::trace_factor_views(&arguments))
                .finish();
        });
        error.map_or(Ok(expression), |reason| {
            Err(NetworkToolingError::Parse { reason })
        })
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
                            plan.status == crate::tensor::ReductionStatus::Complete,
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
        chain, chain_factor, dot, g, mink,
        network::tags::SPENSO_TAG,
        structure::partial::{PartialStructure, PartialStructureExt},
        tensor, trace, vector,
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
}
