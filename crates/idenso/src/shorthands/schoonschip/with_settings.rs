use spenso::network::parsing::AtomStructureExt;
use symbolica::atom::{Atom, AtomView};

use crate::shorthands::bracket::BracketNormalizer;

use super::{
    DotNormalizer, SimplificationCandidates, settings::SchoonschipSettings,
    slot_contraction::SlotContraction,
};

pub(crate) struct SchoonschipWithSettings<'a> {
    pub(crate) settings: &'a SchoonschipSettings,
}

impl SchoonschipWithSettings<'_> {
    pub(crate) fn run<const OBSERVE_REMAINDER: bool>(&self, view: AtomView<'_>) -> Atom {
        let mut current = view.to_owned();
        let mut pending_candidates = None;
        loop {
            let mut components = None;
            let mut attempted_collection = false;
            let mut candidates = pending_candidates.take().unwrap_or_else(|| {
                SimplificationCandidates::scan(current.as_view(), [], || {
                    if !OBSERVE_REMAINDER
                        && self.settings.schoonschip_rank1_tensors
                        && SlotContraction::component_sum_candidate(
                            current.as_view(),
                            self.settings.expand_contracted_sums,
                        )
                    {
                        attempted_collection = true;
                        components = SlotContraction::component_sum(
                            current.as_view(),
                            self.settings.expand_contracted_sums,
                        );
                    }
                    // Failed admission resumes the existing observer-only
                    // suffix before any callback-sensitive cleanup decisions.
                    components.is_none()
                })
            });
            if let Some(components) = components {
                return components;
            }
            if candidates.normalized() {
                return current;
            }
            // Opt-in distribution feeds admitted factors directly to the
            // collector. Unsupported forms keep the ordinary traversal below.
            if !attempted_collection
                && self.settings.schoonschip_rank1_tensors
                && candidates.repeated_indices
                && let Some(components) = SlotContraction::component_sum(
                    current.as_view(),
                    self.settings.expand_contracted_sums,
                )
            {
                return components;
            }
            // Without callbacks or bracket scopes, dots only remove indices.
            // Carry the input's contraction candidate through normalization;
            // no new repeated-index or head scan is needed on its larger result.
            let reusable = candidates.intrinsic && !candidates.brackets;
            let bracketed = if candidates.brackets {
                BracketNormalizer::normalize(current.as_view())
            } else {
                current.clone()
            };
            let normalized = if candidates.dots || bracketed != current {
                DotNormalizer::run(bracketed.as_view())
            } else {
                bracketed
            };
            // Dot normalization stays outside this guard: compact vector
            // rewrites can require no repeated explicit index at all.
            let repeated = if reusable || normalized == current {
                candidates.repeated_indices
            } else {
                normalized.has_repeated_explicit_indices()
            };
            if repeated {
                // The initial expanded collector already declined this exact
                // input. Factored admission and changed normalization remain
                // separate opportunities for the expanded collector.
                let declined_expanded = self.settings.schoonschip_rank1_tensors
                    && !self.settings.expand_contracted_sums
                    && SlotContraction::component_sum_candidate(current.as_view(), false)
                    && normalized == current;
                let contracted = SlotContraction::run(
                    normalized.as_view(),
                    self.settings.simplify_chain_like_functions,
                    self.settings.schoonschip_rank1_tensors,
                    !declined_expanded,
                );
                if contracted != normalized {
                    let dotted = DotNormalizer::normalize(contracted.as_view());
                    let next = if reusable {
                        dotted.expression
                    } else {
                        BracketNormalizer::normalize(dotted.expression.as_view())
                    };
                    // The contractor reached its fixed point. Only an odd
                    // power can expose a new tensor factor during intrinsic
                    // cleanup; keep that work pending instead of rescanning.
                    if reusable && !dotted.exposes_product {
                        return self.finish_changed_sum(next, view);
                    }
                    // The contractor already reaches its fixed point. If
                    // cleanup changes nothing, no new contraction or
                    // normalization can be exposed by another full pass.
                    if next == contracted || next == current {
                        return self.finish_changed_sum(next, view);
                    }
                    current = next;
                    if reusable {
                        candidates.dots = false;
                        candidates.repeated_indices = dotted.exposes_product;
                        pending_candidates = Some(candidates);
                    }
                    continue;
                }
            }
            // Initial normalization can expose another bracket identity even
            // without a contraction, such as bracket(p(mu)^2) becoming scalar.
            if reusable || normalized == current {
                return self.finish_changed_sum(normalized, view);
            }
            current = normalized;
        }
    }

    fn finish_changed_sum(&self, expression: Atom, original: AtomView<'_>) -> Atom {
        // Fallback contractions and dot cleanup can make a previously declined
        // sum eligible (or small enough) for alpha collection. Complete that
        // work now rather than exposing a different result on a second call.
        // Successful initial collection, unchanged inputs, and scalar results
        // retain their existing return paths.
        if self.settings.schoonschip_rank1_tensors
            && expression.as_view() != original
            && (matches!(expression.as_view(), AtomView::Add(_))
                || self.settings.expand_contracted_sums)
            && expression.has_repeated_explicit_indices()
            && let Some(result) = SlotContraction::component_sum(
                expression.as_view(),
                self.settings.expand_contracted_sums,
            )
        {
            return result;
        }
        expression
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::shorthands::schoonschip::Schoonschip;
    use symbolica::{atom::AtomCore, parser::ParseSettings};

    #[test]
    fn early_collection_matches_full_observation_across_settings_and_fallbacks() {
        crate::test_support::test_initialize();
        let _ = spenso::p!(spenso::mink!(4));
        let _ = spenso::q!(spenso::mink!(4));
        let sources = [
            "g(mink(4,a),mink(4,b))*T(mink(4,b))+g(mink(4,a),mink(4,c))*U(mink(4,c))",
            "p(mink(4,a))*(T(mink(4,a))+U(mink(4,a)))",
            "(p(mink(4,a))+q(mink(4,a)))^2",
            "(x+y)^3*p(mink(4,a))*(T(mink(4,a))+U(mink(4,a)))",
            "scope(mink(4,a),mink(4,a),bracket(p(mink(4,b))^2))+x",
            "scope(mink(4,a),mink(4,a),p(mink(4,b))^3)+x",
            "scope(mink(4,a),mink(4,a),T(mink(4,hidden(b))))+x",
            "p(mink(4,a))*(T(mink(4,a))+U(mink(4,b)))",
            "(p(mink(4,a))^(-2)+q(mink(4,a))^(-3))*p(mink(4,a))",
            "x+y",
        ];
        for expand in [false, true] {
            for rank_one in [false, true] {
                for chain_like in [false, true] {
                    let settings = SchoonschipSettings {
                        expand_contracted_sums: expand,
                        schoonschip_rank1_tensors: rank_one,
                        simplify_chain_like_functions: chain_like,
                        ..Default::default()
                    };
                    let owner = SchoonschipWithSettings {
                        settings: &settings,
                    };
                    for source in sources {
                        let expression =
                            Atom::parse(source, "spenso", ParseSettings::symbolica()).unwrap();
                        let full = owner.run::<true>(expression.as_view());
                        let early = owner.run::<false>(expression.as_view());
                        assert_eq!(
                            early, full,
                            "{source}; expand={expand}, rank_one={rank_one}, chain_like={chain_like}"
                        );
                        assert_eq!(
                            owner.run::<false>(early.as_view()),
                            early,
                            "fixed point: {source}"
                        );
                    }
                }
            }
        }
    }

    #[test]
    fn no_repeat_and_admitted_collection_need_only_one_observation() {
        use super::super::analysis::SCAN_COUNTS;
        crate::test_support::test_initialize();
        for head in ["spenso::early_scan_tensor_t", "spenso::early_scan_tensor_u"] {
            spenso::network::tags::SPENSO_TAG.tensor_symbol(head);
        }
        let settings = SchoonschipSettings::default();
        let owner = SchoonschipWithSettings {
            settings: &settings,
        };
        for source in [
            "x+y",
            "g(mink(4,a),mink(4,b))*early_scan_tensor_t(mink(4,b))+g(mink(4,a),mink(4,c))*early_scan_tensor_u(mink(4,c))",
        ] {
            let expression = Atom::parse(source, "spenso", ParseSettings::symbolica()).unwrap();
            SCAN_COUNTS.with(|count| count.set(0));
            let early = owner.run::<false>(expression.as_view());
            assert_eq!(SCAN_COUNTS.with(|count| count.get()), 1, "{source}");
            SCAN_COUNTS.with(|count| count.set(0));
            let full = owner.run::<true>(expression.as_view());
            assert_eq!(SCAN_COUNTS.with(|count| count.get()), 1, "{source}");
            assert_eq!(early, full);
            if expression.has_repeated_explicit_indices() {
                assert_ne!(
                    early, expression,
                    "the tensor fixture must be admitted and contracted"
                );
            } else {
                assert_eq!(early, expression);
            }
        }
    }

    #[test]
    fn unchanged_expanded_decline_is_not_retried_but_factored_fallback_remains() {
        use super::super::{analysis::SCAN_COUNTS, slot_contraction::COLLECTION_COUNTS};
        crate::test_support::test_initialize();
        spenso::network::tags::SPENSO_TAG.tensor_symbol("spenso::retry_tensor");
        let expression = Atom::parse(
            "x+g(mink(4,a),mink(4,b))*retry_tensor(mink(4,b))^2",
            "spenso",
            ParseSettings::symbolica(),
        )
        .unwrap();
        for expand in [false, true] {
            let settings = SchoonschipSettings {
                expand_contracted_sums: expand,
                ..Default::default()
            };
            COLLECTION_COUNTS.with(|counts| counts.set((0, 0)));
            SCAN_COUNTS.with(|count| count.set(0));
            let actual = SchoonschipWithSettings {
                settings: &settings,
            }
            .run::<false>(expression.as_view());
            assert_eq!(
                actual, expression,
                "unsupported opaque powers remain unchanged"
            );
            assert_eq!(
                SCAN_COUNTS.with(|count| count.get()),
                1,
                "refusal must continue the same walk"
            );
            assert_eq!(
                COLLECTION_COUNTS.with(|counts| counts.get()),
                if expand { (1, 1) } else { (1, 0) }
            );
            // Direct ordinary contraction retains its independent expanded
            // attempt; the optimization applies only to the unchanged caller.
            assert_eq!(
                SlotContraction::run(expression.as_view(), false, true, true),
                actual
            );
        }
    }

    #[test]
    fn declined_early_collection_completes_callback_and_dot_tails() {
        use super::super::analysis::SCAN_COUNTS;
        use std::sync::{Arc, Mutex};
        use symbolica::atom::FunctionBuilder;
        crate::test_support::test_initialize();
        let calls = Arc::new(Mutex::new(Vec::new()));
        let observed = Arc::clone(&calls);
        let callback = spenso::tensor_symbol!(
            "early_tail_callback",
            norm = move |value, _| {
                observed.lock().unwrap().push(value.to_owned());
            }
        );
        let vector = spenso::network::tags::SPENSO_TAG.rank_one_tensor_symbol("early_tail_vector");
        let a = spenso::mink!(4, 83801);
        let b = spenso::mink!(4, 83803);
        let tail = FunctionBuilder::new(callback)
            .add_arg(FunctionBuilder::new(vector).add_arg(&b).finish().pow(2))
            .finish();
        let scope = FunctionBuilder::new(symbolica::symbol!("early_tail_scope"))
            .add_arg(&a)
            .add_arg(&a)
            .add_arg(tail)
            .finish();
        let expression = scope + Atom::var(symbolica::symbol!("early_tail_scalar"));
        for expand in [false, true] {
            let settings = SchoonschipSettings {
                expand_contracted_sums: expand,
                ..Default::default()
            };
            let owner = SchoonschipWithSettings {
                settings: &settings,
            };
            calls.lock().unwrap().clear();
            SCAN_COUNTS.with(|count| count.set(0));
            let full = owner.run::<true>(expression.as_view());
            let full_scans = SCAN_COUNTS.with(|count| count.get());
            let transcript = calls.lock().unwrap().clone();
            assert!(
                !transcript.is_empty(),
                "the late power must trigger callback-sensitive cleanup"
            );
            calls.lock().unwrap().clear();
            SCAN_COUNTS.with(|count| count.set(0));
            let early = owner.run::<false>(expression.as_view());
            assert_eq!(early, full);
            assert_eq!(*calls.lock().unwrap(), transcript);
            assert_eq!(
                SCAN_COUNTS.with(|count| count.get()),
                full_scans,
                "refusal resumes the first walk; later changed expressions may need their own observation"
            );
        }
    }

    #[test]
    fn unchanged_contractions_preserve_normalization_fixed_points() {
        crate::test_support::test_initialize();
        let _ = spenso::p!(spenso::mink!(4));
        let _ = spenso::q!(spenso::mink!(4));
        let sources = [
            "g(mink(4,a),mink(4,b))*epsilon(mink(4,c),mink(4,d),mink(4,e),mink(4,f))",
            "bracket(p(mink(4,a))^2)",
            "bracket(g(mink(4,a),mink(4,b))^2)",
            "bracket(p(mink(4,a)),p(mink(4,a)))",
            "bracket(bracket(g(mink(4,a),mink(4,b))^2),p(mink(4,c)))",
            "bracket(p(mink(4,a))^2+q(mink(4,a))^2)^2",
            "g(mink(4,a),mink(4,b))*unknown(mink(4,b))",
            "p(mink(4,a))*q(mink(4,a))",
            "unknown(mink(4,a))*another(mink(4,a))",
            "bracket(g(mink(4,a),mink(4,b)),unknown(mink(4,b)))",
            "(x+y)^6*g(mink(4,a),mink(4,b))*(unknown(mink(4,b))+another(mink(4,b)))",
            "(g(mink(4,a),mink(4,b))+g(mink(4,a),mink(4,c)))*epsilon(mink(4,a),mink(4,d),mink(4,e),mink(4,f))",
            "p(mink(4,a))*(p(mink(4,a))+q(mink(4,a)))",
            "bracket(g(mink(4,a),mink(4,b))*p(mink(4,a))*q(mink(4,b)))",
            "p(mink(4,a))^3*q(mink(4,a))",
            "g(mink(4,a),mink(4,b))^3*p(mink(4,a))*q(mink(4,b))",
            "(p(mink(4,a))^3+q(mink(4,a))^3)*p(mink(4,a))",
            "x*(p(mink(4,a))*q(mink(4,a))+g(mink(4,b),mink(4,c))^4)",
            "(p(mink(4,a))^(-2)+q(mink(4,a))^(-3))*p(mink(4,a))",
        ];
        for (mode, settings) in [
            SchoonschipSettings::default(),
            SchoonschipSettings::default().with_chain_like_functions(),
            SchoonschipSettings::default().without_rank1_tensors(),
        ]
        .into_iter()
        .enumerate()
        {
            for source in sources {
                let expression = Atom::parse(source, "spenso", ParseSettings::symbolica()).unwrap();
                let mut expected = expression.clone();
                loop {
                    // Preserve the previous eager cleanup schedule as an oracle
                    // for scalar powers and unresolved indexed bracket scopes.
                    let normalized =
                        BracketNormalizer::normalize(expected.as_view()).normalize_dots();
                    let contracted = if normalized.has_repeated_explicit_indices() {
                        SlotContraction::run(
                            normalized.as_view(),
                            settings.simplify_chain_like_functions,
                            settings.schoonschip_rank1_tensors,
                            true,
                        )
                    } else {
                        normalized
                    };
                    let next = BracketNormalizer::normalize(contracted.normalize_dots().as_view());
                    if next == expected {
                        break;
                    }
                    expected = next;
                }
                let actual = expression.schoonschip_with_settings(&settings);
                assert_eq!(actual, expected, "{source}, mode {mode}");
                assert_eq!(
                    actual.schoonschip_with_settings(&settings),
                    actual,
                    "{source}"
                );
            }
        }
    }

    #[test]
    fn shared_state_keeps_callback_discovery_and_rechecks_introduced_syntax() {
        use spenso::network::library::symbolic::ETS;
        use std::sync::{Arc, Mutex};
        use symbolica::atom::FunctionBuilder;

        crate::test_support::test_initialize();
        let calls = Arc::new(Mutex::new(Vec::new()));
        let observed = Arc::clone(&calls);
        let leaf = spenso::tensor_symbol!(
            "state_callback_leaf",
            norm = move |value, _| observed.lock().unwrap().push(value.to_owned())
        );
        let observed = Arc::clone(&calls);
        let wrapper = spenso::tensor_symbol!(
            "state_callback_wrapper",
            norm = move |value, _| observed.lock().unwrap().push(value.to_owned())
        );
        let a = spenso::mink!(4, 73231);
        let b = spenso::mink!(4, 73237);
        let input = ETS.metric(&a, &b) * FunctionBuilder::new(leaf).add_arg(&a).finish();
        let output = FunctionBuilder::new(leaf).add_arg(&b).finish();
        let nested_input = FunctionBuilder::new(wrapper).add_arg(input).finish();
        let nested_output = FunctionBuilder::new(wrapper).add_arg(&output).finish();
        calls.lock().unwrap().clear();
        assert_eq!(nested_input.schoonschip(), nested_output);
        assert_eq!(
            *calls.lock().unwrap(),
            [output.clone(), output, nested_output]
        );

        // A callback can introduce a bracket and a contraction after a power
        // rewrite. These flags must be rediscovered from the actual result.
        let emitted =
            spenso::bracket!(ETS.metric(&a, &b) * FunctionBuilder::new(leaf).add_arg(&a).finish());
        let expected = emitted.schoonschip();
        let trigger = spenso::tensor_symbol!(
            "state_callback_rewrite",
            norm = move |value, out| {
                if let AtomView::Fun(function) = value
                    && matches!(function.iter().next(), Some(AtomView::Fun(inner))
                    if inner.get_symbol() == ETS.metric)
                {
                    **out = emitted.clone();
                }
            }
        );
        let vector =
            spenso::network::tags::SPENSO_TAG.rank_one_tensor_symbol("state_callback_vector");
        let input = FunctionBuilder::new(trigger)
            .add_arg(FunctionBuilder::new(vector).add_arg(&a).finish().pow(2))
            .finish();
        assert_eq!(input.schoonschip(), expected);
    }
}

#[cfg(test)]
mod fused_expansion_tests {
    // Proposed child test module for with_settings.rs. No production implementation.
    use crate::{
        shorthands::schoonschip::{Schoonschip, SchoonschipSettings},
        tensor::SymbolicTensor,
    };
    use spenso::{
        network::{library::symbolic::ETS, tags::SPENSO_TAG},
        structure::partial::{PartialStructure, PartialStructureExt},
    };
    use std::sync::{Arc, Mutex};
    use symbolica::{
        atom::{Atom, AtomCore, AtomView, FunctionBuilder},
        parser::ParseSettings,
    };

    fn setup() {
        crate::test_support::test_initialize();
        for head in ["p", "q"] {
            SPENSO_TAG.rank_one_tensor_symbol(&format!("fused_raw_routing::{head}"));
        }
        for head in ["t", "u"] {
            SPENSO_TAG.tensor_symbol(&format!("fused_raw_routing::{head}"));
        }
    }

    fn input(source: &str) -> Atom {
        Atom::parse(source, "fused_raw_routing", ParseSettings::symbolica()).unwrap()
    }

    #[test]
    fn fused_raw_setting_contracts_opaque_sum_products_only_when_enabled() {
        setup();
        let source = input(
            "(t(spenso::mink(4,a),spenso::mink(4,b))+u(spenso::mink(4,a),spenso::mink(4,b)))*(p(spenso::mink(4,a))+q(spenso::mink(4,a)))",
        );
        let expected = input(
            "t(p(spenso::mink(4)),spenso::mink(4,b))+t(q(spenso::mink(4)),spenso::mink(4,b))+u(p(spenso::mink(4)),spenso::mink(4,b))+u(q(spenso::mink(4)),spenso::mink(4,b))",
        );
        let settings = SchoonschipSettings::default().with_expanded_contracted_sums();
        assert_eq!(
            source.schoonschip(),
            source,
            "default keeps the existing factorized route"
        );
        assert_eq!(source.schoonschip_with_settings(&settings), expected);
        assert_ne!(
            source, expected,
            "success must perform contraction, not only agree after expansion"
        );
        assert_eq!(expected.schoonschip_with_settings(&settings), expected);
        assert_eq!(source.expand().schoonschip(), expected);
        let disabled = settings.without_rank1_tensors();
        assert_eq!(source.schoonschip_with_settings(&disabled), source);
    }

    #[test]
    fn fused_raw_setting_preserves_common_scalar_spectator_factorization() {
        setup();
        let spectator = input(
            "(x+y)^3*(spenso::g(p(spenso::mink(4)),p(spenso::mink(4)))+spenso::g(q(spenso::mink(4)),q(spenso::mink(4))))",
        );
        let tensor = input(
            "(t(spenso::mink(4,a),spenso::mink(4,b))+u(spenso::mink(4,a),spenso::mink(4,b)))*(p(spenso::mink(4,a))+q(spenso::mink(4,a)))",
        );
        let contracted = input(
            "t(p(spenso::mink(4)),spenso::mink(4,b))+t(q(spenso::mink(4)),spenso::mink(4,b))+u(p(spenso::mink(4)),spenso::mink(4,b))+u(q(spenso::mink(4)),spenso::mink(4,b))",
        );
        let source = &spectator * tensor;
        let expected = &spectator * contracted;
        let settings = SchoonschipSettings::default().with_expanded_contracted_sums();
        let actual = source.schoonschip_with_settings(&settings);
        assert_eq!(
            actual, expected,
            "symbolic scalar factors stay outside the tensor sum"
        );
        assert_ne!(
            actual,
            actual.expand(),
            "this fixture distinguishes factor protection from global expansion"
        );
        assert_eq!(actual.expand(), source.expand().schoonschip().expand());
        assert_eq!(actual.schoonschip_with_settings(&settings), actual);
    }

    #[test]
    fn fused_raw_setting_preserves_typed_logical_order_and_zero() {
        setup();
        let product = input(
            "(t(spenso::mink(4,a),spenso::mink(4,b),spenso::mink(6,c))+u(spenso::mink(4,a),spenso::mink(4,b),spenso::mink(6,c)))*(p(spenso::mink(4,a))+q(spenso::mink(4,a)))",
        );
        let expected = input(
            "t(p(spenso::mink(4)),spenso::mink(4,b),spenso::mink(6,c))+t(q(spenso::mink(4)),spenso::mink(4,b),spenso::mink(6,c))+u(p(spenso::mink(4)),spenso::mink(4,b),spenso::mink(6,c))+u(q(spenso::mink(4)),spenso::mink(4,b),spenso::mink(6,c))",
        );
        let settings = SchoonschipSettings::default().with_expanded_contracted_sums();
        for (source, expected) in [
            (product.clone(), expected.clone()),
            (product - expected, Atom::Zero),
        ] {
            let inferred = SymbolicTensor::<PartialStructure>::infer(source).unwrap();
            let mut tensor =
                SymbolicTensor::checked_parts(inferred.expression, inferred.structure).unwrap();
            tensor.structure = PartialStructure::from_logical_slots(
                tensor.structure.logical_slots().into_iter().rev(),
            );
            let original_ports = tensor.structure.logical_slots();
            assert_eq!(original_ports.len(), 2);
            let output = tensor.expression.schoonschip_with_settings(&settings);
            assert_eq!(output, expected);
            let result = tensor.with_rewritten_expression(output).unwrap();
            assert_eq!(result.structure.logical_slots(), original_ports);
            assert_eq!(
                result.expression.schoonschip_with_settings(&settings),
                result.expression
            );
        }
    }

    #[test]
    fn fused_raw_setting_coalesces_vector_powers_before_counting_ports() {
        setup();
        let settings = SchoonschipSettings::default().with_expanded_contracted_sums();
        for (source, expected) in [
            (
                "(p(spenso::mink(4,a))+q(spenso::mink(4,a)))^2",
                "spenso::g(p(spenso::mink(4)),p(spenso::mink(4)))+2*spenso::g(p(spenso::mink(4)),q(spenso::mink(4)))+spenso::g(q(spenso::mink(4)),q(spenso::mink(4)))",
            ),
            (
                "p(spenso::mink(4,a))*(p(spenso::mink(4,a))+q(spenso::mink(4,a)))*(p(spenso::mink(4,a))-q(spenso::mink(4,a)))",
                "p(spenso::mink(4,a))*(spenso::g(p(spenso::mink(4)),p(spenso::mink(4)))-spenso::g(q(spenso::mink(4)),q(spenso::mink(4))))",
            ),
        ] {
            let source = input(source);
            let expected = input(expected).expand();
            let actual = source.schoonschip_with_settings(&settings);
            assert_eq!(actual, expected, "{source}");
            assert_eq!(
                actual,
                source.expand().schoonschip(),
                "normalized duplicate powers retain their local-pair rule"
            );
            assert_eq!(actual.schoonschip_with_settings(&settings), actual);
        }
    }

    #[test]
    fn fused_raw_setting_does_not_relax_missing_sum_branch_interfaces() {
        setup();
        let source = input("p(spenso::mink(4,a))*(t(spenso::mink(4,a))+u(spenso::mink(4,b)))");
        assert!(SymbolicTensor::<PartialStructure>::infer(source.clone()).is_err());
        assert!(SymbolicTensor::<PartialStructure>::infer(source.expand()).is_err());
        let settings = SchoonschipSettings::default().with_expanded_contracted_sums();
        let expected = input("t(p(spenso::mink(4)))+p(spenso::mink(4,a))*u(spenso::mink(4,b))");
        assert_eq!(
            source.schoonschip_with_settings(&settings),
            expected,
            "raw monomial contraction does not invent a missing branch port"
        );
        assert!(SymbolicTensor::<PartialStructure>::infer(expected).is_err());
    }

    #[test]
    fn fused_raw_setting_preserves_callback_fallback_and_hidden_metadata_schedule() {
        setup();
        let calls = Arc::new(Mutex::new(Vec::new()));
        let observed = Arc::clone(&calls);
        let head = spenso::tensor_symbol!(
            "fused_raw_callback_leaf",
            norm = move |value, _| observed.lock().unwrap().push(value.to_owned())
        );
        let a = spenso::mink!(4, 78311);
        let b = spenso::mink!(4, 78317);
        let leaf = FunctionBuilder::new(head).add_arg(&a).finish();
        let product = ETS.metric(&a, &b) * leaf;
        let scalar = symbolica::symbol!("fused_raw_routing::opaque_metadata"; Scalar);
        let metadata = FunctionBuilder::new(scalar).add_arg(&product).finish();
        let outside = input("(p(spenso::mink(4,c))+q(spenso::mink(4,c)))*t(spenso::mink(4,c))");
        let settings = SchoonschipSettings::default().with_expanded_contracted_sums();
        for source in [product, metadata * outside] {
            calls.lock().unwrap().clear();
            let expected = source.schoonschip();
            let transcript = calls.lock().unwrap().clone();
            assert!(!transcript.is_empty());
            calls.lock().unwrap().clear();
            assert_eq!(source.schoonschip_with_settings(&settings), expected);
            assert_eq!(
                *calls.lock().unwrap(),
                transcript,
                "declined planning must not replay normalization"
            );
            calls.lock().unwrap().clear();
            assert_eq!(expected.schoonschip_with_settings(&settings), expected);
            assert!(calls.lock().unwrap().is_empty());
        }
    }

    #[test]
    fn fused_raw_setting_keeps_checked_metric_callback_rank_loss_rejection() {
        setup();
        let calls = Arc::new(Mutex::new(Vec::new()));
        let observed = Arc::clone(&calls);
        let a = spenso::mink!(4, 78401);
        let b = spenso::mink!(4, 78403);
        let target = b.clone();
        let head = spenso::tensor_symbol!(
            "fused_raw_metric_rank_loss",
            norm = move |value, output| {
                observed.lock().unwrap().push(value.to_owned());
                if let AtomView::Fun(function) = value
                    && function.iter().next() == Some(target.as_view())
                {
                    **output = Atom::num(7);
                }
            }
        );
        let source = ETS.metric(&a, &b) * FunctionBuilder::new(head).add_arg(&a).finish();
        let inferred = SymbolicTensor::<PartialStructure>::infer(source.clone()).unwrap();
        let tensor =
            SymbolicTensor::checked_parts(inferred.expression, inferred.structure).unwrap();
        calls.lock().unwrap().clear();
        let expected = source.schoonschip();
        let transcript = calls.lock().unwrap().clone();
        assert_eq!(expected, Atom::num(7));
        assert!(!transcript.is_empty());
        calls.lock().unwrap().clear();
        let output = source.schoonschip_with_settings(
            &SchoonschipSettings::default().with_expanded_contracted_sums(),
        );
        assert_eq!(output, expected);
        assert_eq!(*calls.lock().unwrap(), transcript);
        calls.lock().unwrap().clear();
        assert!(tensor.with_rewritten_expression(output).is_err());
        assert!(
            calls.lock().unwrap().is_empty(),
            "observing rank loss must not invoke the callback again"
        );
    }
}
