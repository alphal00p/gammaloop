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
    pub(crate) fn run(
        &self,
        view: AtomView<'_>,
        literal_relabellings: &mut Vec<(Atom, Atom)>,
    ) -> Atom {
        let mut current = view.to_owned();
        let mut pending_candidates = None;
        loop {
            let mut candidates = pending_candidates
                .take()
                .unwrap_or_else(|| SimplificationCandidates::scan(current.as_view(), [], || true));
            if candidates.normalized() {
                return current;
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
                let contracted = SlotContraction::run(
                    normalized.as_view(),
                    self.settings.simplify_chain_like_functions,
                    self.settings.schoonschip_rank1_tensors,
                    literal_relabellings,
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
                        return next;
                    }
                    // The contractor already reaches its fixed point. If
                    // cleanup changes nothing, no new contraction or
                    // normalization can be exposed by another full pass.
                    if next == contracted || next == current {
                        return next;
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
                return normalized;
            }
            current = normalized;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::shorthands::schoonschip::Schoonschip;
    use symbolica::{atom::AtomCore, parser::ParseSettings};

    #[test]
    fn ordered_fallback_preserves_fixedpoints_across_settings_and_scopes() {
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
        for rank_one in [false, true] {
            for chain_like in [false, true] {
                let settings = SchoonschipSettings {
                    schoonschip_rank1_tensors: rank_one,
                    simplify_chain_like_functions: chain_like,
                };
                let owner = SchoonschipWithSettings {
                    settings: &settings,
                };
                for source in sources {
                    let expression =
                        Atom::parse(source, "spenso", ParseSettings::symbolica()).unwrap();
                    let result = owner.run(expression.as_view(), &mut Vec::new());
                    assert_eq!(
                        owner.run(result.as_view(), &mut Vec::new()),
                        result,
                        "{source}; rank_one={rank_one}, chain_like={chain_like}"
                    );
                }
            }
        }
    }

    #[test]
    fn scalar_and_local_contraction_need_only_one_observation() {
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
            let result = owner.run(expression.as_view(), &mut Vec::new());
            assert_eq!(SCAN_COUNTS.with(|count| count.get()), 1, "{source}");
            if expression.has_repeated_explicit_indices() {
                assert_ne!(result, expression, "the tensor fixture must contract");
            } else {
                assert_eq!(result, expression);
            }
        }
    }

    #[test]
    fn unsupported_opaque_power_is_unchanged_after_one_observation() {
        use super::super::analysis::SCAN_COUNTS;
        crate::test_support::test_initialize();
        spenso::network::tags::SPENSO_TAG.tensor_symbol("spenso::retry_tensor");
        let expression = Atom::parse(
            "x+g(mink(4,a),mink(4,b))*retry_tensor(mink(4,b))^2",
            "spenso",
            ParseSettings::symbolica(),
        )
        .unwrap();
        SCAN_COUNTS.with(|count| count.set(0));
        let actual = expression.schoonschip();
        assert_eq!(actual, expression);
        assert_eq!(SCAN_COUNTS.with(|count| count.get()), 1);
        assert_eq!(
            SlotContraction::run(expression.as_view(), false, true, &mut Vec::new()),
            actual
        );
    }

    #[test]
    fn repeated_indices_do_not_hide_callback_and_dot_tails() {
        use std::sync::{Arc, Mutex};
        use symbolica::atom::FunctionBuilder;
        crate::test_support::test_initialize();
        let calls = Arc::new(Mutex::new(Vec::new()));
        let observed = Arc::clone(&calls);
        let callback = spenso::tensor_symbol!(
            "early_tail_callback",
            norm = move |value, _| observed.lock().unwrap().push(value.to_owned())
        );
        let vector = spenso::network::tags::SPENSO_TAG.rank_one_tensor_symbol("early_tail_vector");
        let a = spenso::mink!(4, 83801);
        let b = spenso::mink!(4, 83803);
        let vector = FunctionBuilder::new(vector).add_arg(&b).finish();
        let tail = FunctionBuilder::new(callback)
            .add_arg(vector.pow(2))
            .finish();
        let scope_head = symbolica::symbol!("early_tail_scope");
        let scalar = Atom::var(symbolica::symbol!("early_tail_scalar"));
        let expression = FunctionBuilder::new(scope_head)
            .add_arg(&a)
            .add_arg(&a)
            .add_arg(tail)
            .finish()
            + &scalar;
        let expected_tail = FunctionBuilder::new(callback)
            .add_arg(vector.pow(2).normalize_dots())
            .finish();
        let expected = FunctionBuilder::new(scope_head)
            .add_arg(&a)
            .add_arg(&a)
            .add_arg(&expected_tail)
            .finish()
            + scalar;
        calls.lock().unwrap().clear();
        let result = expression.schoonschip();
        assert_eq!(result, expected);
        assert!(calls.lock().unwrap().contains(&expected_tail));
        calls.lock().unwrap().clear();
        assert_eq!(result.schoonschip(), result);
        assert!(calls.lock().unwrap().is_empty());
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
                            &mut Vec::new(),
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
