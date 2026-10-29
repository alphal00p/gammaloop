//! Whole-expression handoff to RustRed's ordered scheduler.
//!
//! Matching, scalar lowering and normalization retain their existing owners.
//! Only compatible unit-mass coefficient channels share a native frontier.

use rustred::solver::CandidateRoutedApplicationObserver;
use symbolica::atom::Atom;

use super::{MatchedContribution, RustRedEvaluationError, four_loop, materialize};
use crate::{Vakint, VakintError, VakintSettings};

impl Vakint {
    pub(crate) fn rustred_reduce_ordered(
        &self,
        settings: &VakintSettings,
        contributions: Vec<MatchedContribution>,
        observation: Option<&mut CandidateRoutedApplicationObserver<'_>>,
        progress: Option<&super::monitor::CallerProgress>,
    ) -> Result<Atom, VakintError> {
        let Some(first) = contributions.first() else {
            return Ok(Atom::Zero);
        };
        let loop_count = first.matched.canonical_topology.get_integral().n_loops;
        if contributions.iter().any(|contribution| {
            contribution
                .matched
                .canonical_topology
                .get_integral()
                .n_loops
                != loop_count
        }) {
            return Err(RustRedEvaluationError::InvalidMatchedFamily {
                detail: "one ordered family batch contains different loop counts".into(),
            }
            .into());
        }
        let mut scalar_settings = settings.clone();
        scalar_settings.use_dot_product_notation = true;
        crate::citations::USED_CITATIONS.record(crate::citations::CitationSource::RustRed);
        log::debug!(target: "vakint::ordered_reduction",
            "ordered reduction: {} matched contributions; DFS result caches inactive",
            contributions.len());
        match loop_count {
            1..=3 => materialize::evaluate_ordered_batch(
                self,
                &scalar_settings,
                &contributions,
                observation,
                progress,
            )
            .map_err(Into::into),
            4 => four_loop::evaluate_ordered_batch(
                self,
                &scalar_settings,
                &contributions,
                observation,
                progress,
            ),
            5 => self
                .saved_application(&scalar_settings)?
                .evaluate_ordered_batch(&scalar_settings, &contributions, observation, progress)
                .map(|result| result.value),
            _ => Err(RustRedEvaluationError::UnsupportedMatchedFamily {
                detail: format!("no loaded adapter for {loop_count} loops"),
            }
            .into()),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        EvaluationMethod, EvaluationOrder, LoopNormalizationFactor, RustRedEvaluationOptions,
        RustRedMasterMode, RustRedReductionStrategy, vakint_parse,
    };
    use rustred::solver::{CandidateRoutedObservationOptions, CandidateRoutedObserverControl};
    use symbolica::atom::AtomCore;

    fn settings(strategy: RustRedReductionStrategy) -> VakintSettings {
        let mut settings = VakintSettings {
            evaluation_order: EvaluationOrder(vec![EvaluationMethod::RustRed(
                RustRedEvaluationOptions {
                    master_mode: RustRedMasterMode::Retain,
                },
            )]),
            integral_normalization_factor: LoopNormalizationFactor::FMFTandMATAD,
            use_dot_product_notation: true,
            form_exe_path: "/ordered-reduction-must-not-call-form".into(),
            ..Default::default()
        };
        settings.rustred_application_options.reduction_strategy = strategy;
        settings
    }

    fn assert_parity(inputs: &[&str]) {
        let engine = Vakint::new().unwrap();
        for input in inputs {
            let input = vakint_parse!(input).unwrap();
            let depth_first = engine
                .evaluate_integral(
                    &settings(RustRedReductionStrategy::DepthFirst),
                    input.as_view(),
                )
                .unwrap();
            let ordered = engine
                .evaluate_integral(
                    &settings(RustRedReductionStrategy::Ordered),
                    input.as_view(),
                )
                .unwrap();
            assert!(
                (depth_first - ordered).together().expand().is_zero(),
                "{input}"
            );
        }
    }

    #[test]
    fn one_loop_ordered_and_depth_first_agree() {
        assert_parity(&["3*topo(I1L(muvsq,7))-2*topo(I1L(muvsq,8))+topo(I1L(muvsq,9))"]);
    }

    #[test]
    fn ordered_msbar_retained_normalization_matches_depth_first() {
        let engine = Vakint::new().unwrap();
        let euler_gamma = crate::vakint_symbol!("EulerGamma");
        for input in [
            "topo(I1L(1,3))+2*topo(I1L(1,4))",
            "topo(I2L(7/3,2,1,1))+2*topo(I2L(7/3,1,2,1))",
        ] {
            let input = vakint_parse!(input).unwrap();
            let mut results = Vec::new();
            for strategy in [
                RustRedReductionStrategy::DepthFirst,
                RustRedReductionStrategy::Ordered,
            ] {
                let mut settings = settings(strategy);
                settings.integral_normalization_factor = LoopNormalizationFactor::MSbar;
                let result = engine
                    .evaluate_integral(&settings, input.as_view())
                    .unwrap();
                assert!(
                    !result.contains_symbol(euler_gamma),
                    "MSbar must cancel normalization exponentials before retaining masters: {result}"
                );
                results.push(result);
            }
            assert!((&results[0] - &results[1]).together().expand().is_zero());
        }
    }

    #[test]
    fn ordered_msbar_normalization_is_shared_by_direct_entry_points() {
        let engine = Vakint::new().unwrap();
        let mut settings = settings(RustRedReductionStrategy::Ordered);
        settings.integral_normalization_factor = LoopNormalizationFactor::MSbar;
        for integral in ["topo(I1L(1,3))", "topo(I2L(1,2,1,1))"] {
            let integral = vakint_parse!(integral).unwrap();
            let batched = engine
                .evaluate_integral(&settings, integral.as_view())
                .unwrap();
            let mut term = crate::VakintTerm {
                integral: integral.clone(),
                numerator: Atom::one(),
                vectors: Vec::new(),
            };
            let matched = term
                .prepare_integral_specs(&engine, &settings, None)
                .unwrap();
            let direct = settings.evaluation_order.0[0]
                .evaluate_integral(&engine, &settings, term.numerator.as_view(), &matched)
                .unwrap();
            let mut term = crate::VakintTerm {
                integral,
                numerator: Atom::one(),
                vectors: Vec::new(),
            };
            term.evaluate_integral(&engine, &settings).unwrap();
            assert_eq!(direct, batched);
            assert_eq!(Atom::from(term), batched);
        }
    }

    #[test]
    fn ordered_msbar_substitution_preserves_numerical_precision() {
        use ahash::HashMap;

        let engine = Vakint::new().unwrap();
        for (loop_count, input) in [(1, "topo(I1L(1,1))"), (2, "topo(I2L(1,1,1,1))")] {
            let input = vakint_parse!(input).unwrap();
            let mut values = Vec::new();
            for mode in [RustRedMasterMode::Symbolic, RustRedMasterMode::Numerical] {
                for strategy in [
                    RustRedReductionStrategy::DepthFirst,
                    RustRedReductionStrategy::Ordered,
                ] {
                    let mut settings = settings(strategy);
                    settings.integral_normalization_factor = LoopNormalizationFactor::MSbar;
                    settings.number_of_terms_in_epsilon_expansion = loop_count + 1;
                    settings.run_time_decimal_precision = 40;
                    settings.evaluation_order = EvaluationOrder(vec![EvaluationMethod::RustRed(
                        RustRedEvaluationOptions { master_mode: mode },
                    )]);
                    let expression = engine
                        .evaluate_integral(&settings, input.as_view())
                        .unwrap()
                        .replace(
                            Atom::var(crate::vakint_symbol!(&settings.mu_r_sq_symbol)).to_pattern(),
                        )
                        .with(Atom::one());
                    let (value, statistical_error) = engine
                        .numerical_evaluation(
                            &settings,
                            expression.as_view(),
                            &HashMap::default(),
                            &HashMap::default(),
                            None,
                        )
                        .unwrap();
                    assert!(statistical_error.is_none());
                    for coefficient in value.coefficients().values() {
                        for part in [
                            coefficient.precision().real(),
                            coefficient.precision().imaginary(),
                        ] {
                            assert!(
                                part.exact_zero()
                                    || part.significant_digits().is_some_and(|digits| digits >= 37)
                            );
                        }
                    }
                    values.push(value);
                }
            }
            for actual in &values[1..] {
                let (equal, diagnostic) = actual.does_approx_match(&values[0], None, 1e-30, 5.0);
                assert!(equal, "{diagnostic}");
            }
        }
    }

    #[test]
    fn two_loop_ordered_and_depth_first_agree() {
        assert_parity(&[
            "topo(I2L(muvsq,2,1,1))+2*topo(I2L(muvsq,1,2,1))",
            "dot(k(1),k(2))*topo(I2L(muvsq,1,1,1))",
        ]);
    }

    #[test]
    fn two_loop_complex_coefficients_agree_with_depth_first() {
        assert_parity(&[
            "-8*𝑖/3*topo(I2L(muvsq,2,1,1))",
            "(2/3+3*𝑖/5)*topo(I2L(muvsq,2,1,1))+(-2/3+7*𝑖/5)*topo(I2L(muvsq,3,1,1))",
        ]);
    }

    #[test]
    fn complex_dimension_numerators_share_channels_across_regulator_aliases() {
        let engine = Vakint::new().unwrap();
        let mut ordered_settings = settings(RustRedReductionStrategy::Ordered);
        ordered_settings.epsilon_symbol = "ordered_carrier_control::epsilon".into();
        let mut depth_first_settings = ordered_settings.clone();
        depth_first_settings
            .rustred_application_options
            .reduction_strategy = RustRedReductionStrategy::DepthFirst;
        let directory = tempfile::tempdir().unwrap();
        for (alias_index, epsilon) in [ordered_settings.epsilon_symbol.as_str(), "vakint::ep"]
            .into_iter()
            .enumerate()
        {
            let d = format!("(4-2*{epsilon})");
            let inputs = [
                format!("({d}+𝑖)/({d}+1)*topo(I2L(1,2,1,1))"),
                format!("{d}/({d}+1)*topo(I2L(1,2,1,1))+𝑖/({d}+1)*topo(I2L(1,2,1,1))"),
            ];
            let mut values = Vec::new();
            for (form_index, input) in inputs.iter().enumerate() {
                let input = vakint_parse!(input).unwrap();
                let snapshot = directory
                    .path()
                    .join(format!("{alias_index}-{form_index}.json"));
                let ordered = engine
                    .evaluate_with_monitor(
                        &ordered_settings,
                        input.as_view(),
                        &crate::RustRedMonitorOptions {
                            progress_snapshot: Some(snapshot.clone()),
                            interval: std::time::Duration::ZERO,
                            ..Default::default()
                        },
                    )
                    .unwrap();
                let depth_first = engine
                    .evaluate_integral(&depth_first_settings, input.as_view())
                    .unwrap();
                assert!((&ordered - depth_first).together().expand().is_zero());
                let document: serde_json::Value =
                    serde_json::from_reader(std::fs::File::open(snapshot).unwrap()).unwrap();
                assert_eq!(document["monitor"]["reduction_cohorts"]["total"], 2);
                values.push(ordered);
            }
            assert!((&values[0] - &values[1]).together().expand().is_zero());
            let original_pole = vakint_parse!(&format!(
                "(rustred::d+𝑖)*topo(I2L(1,2,1,1))/(rustred::d+2*{epsilon}-4)"
            ))
            .unwrap();
            let error = engine
                .evaluate_integral(&ordered_settings, original_pole.as_view())
                .unwrap_err();
            assert!(matches!(error, VakintError::RustRedEvaluation(
                RustRedEvaluationError::InvalidMatchedFamily { ref detail }
            ) if detail.contains("original scalar pole")));
        }
    }

    #[test]
    fn three_loop_ordered_and_depth_first_agree() {
        assert_parity(&["topo(I3L(muvsq,2,1,1,1,1,1))"]);
    }

    #[test]
    fn ordered_combines_distinct_roots_before_full_reduction() {
        let engine = Vakint::new().unwrap();
        let settings = settings(RustRedReductionStrategy::Ordered);
        for mass in ["1", "7/3", "muvsq"] {
            for epsilon in [settings.epsilon_symbol.as_str(), "vakint::ep"] {
                for phase in ["1", "𝑖", "-8*𝑖/3", "2/3+3*𝑖/5"] {
                    let input = vakint_parse!(&format!(
                        "({phase})*2*({mass})*topo(I1L({mass},8))+({phase})*(10+2*{})/7*topo(I1L({mass},7))",
                        epsilon,
                    ))
                    .unwrap();
                    let mut callback =
                        |_: &rustred::solver::CandidateRoutedApplicationEnvelope<'_>| {
                            Ok(CandidateRoutedObserverControl::Continue)
                        };
                    let mut observer = CandidateRoutedApplicationObserver::new(
                        CandidateRoutedObservationOptions::default(),
                        &mut callback,
                    );
                    let result = engine
                        .evaluate_observed(&settings, input.as_view(), Some(&mut observer), None)
                        .unwrap();
                    assert!(result.is_zero(), "mass {mass}: {result}");
                    let progress = observer.last_ordered_progress().unwrap();
                    assert_eq!(progress.processed_keys, 1, "mass {mass}");
                    assert!(progress.cancelled_terms > 0, "mass {mass}");
                    assert_eq!(progress.completed_phases, progress.total_phases);
                }
            }
        }
    }

    #[test]
    fn strategy_switches_keep_depth_first_default() {
        assert_eq!(
            VakintSettings::default()
                .rustred_application_options
                .reduction_strategy,
            RustRedReductionStrategy::DepthFirst
        );
    }

    #[test]
    fn reduction_math_capture_observes_both_engines_without_changing_results() {
        use crate::rustred_evaluation::experimental::timing::capture_reduction_timing;

        let engine = Vakint::new().unwrap();
        let input = vakint_parse!("topo(I1L(1,4))+2*topo(I1L(1,5))").unwrap();
        for strategy in [
            RustRedReductionStrategy::DepthFirst,
            RustRedReductionStrategy::Ordered,
        ] {
            let settings = settings(strategy);
            let expected = engine
                .evaluate_integral(&settings, input.as_view())
                .unwrap();
            let captured =
                capture_reduction_timing(|| engine.evaluate_integral(&settings, input.as_view()));
            let actual = captured.result.unwrap();
            assert!((actual - expected).together().expand().is_zero());
            assert!(captured.timing.complete, "{:?}", captured.timing);
            assert_eq!(captured.timing.scope_id, "vakint.rustred.reduction-math.v1");
            assert!(captured.timing.interval_count > 0);
            assert_eq!(
                captured.timing.measured_owner_preparation,
                std::time::Duration::ZERO
            );
            assert_eq!(captured.timing.reduction_math, Some(captured.timing.gross));
        }
    }

    #[test]
    fn ordered_normal_monitor_reports_all_loop_batches_and_mass_cohorts() {
        let engine = Vakint::new().unwrap();
        let settings = settings(RustRedReductionStrategy::Ordered);
        let input =
            vakint_parse!("topo(I1L(1,2))+topo(I1L(2,3))+topo(I2L(1,1,1,1))+topo(I2L(2,1,1,1))")
                .unwrap();
        let directory = tempfile::tempdir().unwrap();
        let snapshot = directory.path().join("ordered-progress.json");
        let result = engine
            .evaluate_with_monitor(
                &settings,
                input.as_view(),
                &crate::RustRedMonitorOptions {
                    progress_snapshot: Some(snapshot.clone()),
                    interval: std::time::Duration::ZERO,
                    ..Default::default()
                },
            )
            .unwrap();
        assert!(!result.is_zero());
        let document: serde_json::Value =
            serde_json::from_reader(std::fs::File::open(snapshot).unwrap()).unwrap();
        assert_eq!(document["monitor"]["ordered_batches"]["completed"], 2);
        assert_eq!(document["monitor"]["ordered_batches"]["total"], 2);
        assert_eq!(document["monitor"]["reduction_cohorts"]["completed"], 2);
        assert_eq!(document["monitor"]["reduction_cohorts"]["total"], 2);
        assert_eq!(
            document["monitor"]["reduction_cohorts"]["fraction_is_time_fraction"],
            false
        );
    }

    #[test]
    fn ordered_rejects_original_poles_that_vanish_under_dimension_identification() {
        let engine = Vakint::new().unwrap();
        let settings = settings(RustRedReductionStrategy::Ordered);
        for epsilon in [settings.epsilon_symbol.as_str(), "vakint::ep"] {
            let input =
                vakint_parse!(&format!("topo(I1L(1,1))/(rustred::d+2*{epsilon}-4)")).unwrap();
            let error = engine
                .evaluate_integral(&settings, input.as_view())
                .unwrap_err();
            assert!(
                matches!(error, VakintError::RustRedEvaluation(
                RustRedEvaluationError::InvalidMatchedFamily { ref detail }
            ) if detail.contains("original scalar pole")),
                "{error}"
            );

            // A visible pole that remains nonzero for generic epsilon may be
            // removable after substitution; this is valid on its original domain.
            let removable =
                vakint_parse!(&format!("(rustred::d-4)/({epsilon})*topo(I1L(1,1))")).unwrap();
            let actual = engine
                .evaluate_integral(&settings, removable.as_view())
                .unwrap();
            let master = engine
                .evaluate_integral(
                    &settings,
                    vakint_parse!("topo(I1L(1,1))").unwrap().as_view(),
                )
                .unwrap();
            assert!(
                (actual + Atom::num(2) * master)
                    .together()
                    .expand()
                    .is_zero()
            );
        }
    }

    #[test]
    fn ordered_requests_leave_attached_depth_first_cache_untouched() {
        use crate::{RustRedCachePreparationOptions, RustRedCacheRoot, RustRedCacheWorkload};
        use rustred::family::IntegralKey;
        use std::sync::atomic::AtomicBool;

        let engine = Vakint::new().unwrap();
        let directory = tempfile::tempdir().unwrap();
        let path = directory.path().join("depth-first-cache");
        let options = RustRedCachePreparationOptions {
            workload: RustRedCacheWorkload::Explicit(vec![RustRedCacheRoot {
                namespace: 0,
                integral: IntegralKey::try_new([7]).unwrap(),
            }]),
            ..Default::default()
        };
        engine
            .prepare_cache(1, &path, &options, &AtomicBool::new(false), |_| Ok(()))
            .unwrap();
        engine.load_cache(1, &path, &Default::default()).unwrap();
        let before = engine.cache_statistics(1).unwrap();
        let input = vakint_parse!("topo(I1L(1,7))").unwrap();
        let ordered = engine
            .evaluate_integral(
                &settings(RustRedReductionStrategy::Ordered),
                input.as_view(),
            )
            .unwrap();
        assert_eq!(engine.cache_statistics(1).unwrap(), before);
        let depth_first = engine
            .evaluate_integral(
                &settings(RustRedReductionStrategy::DepthFirst),
                input.as_view(),
            )
            .unwrap();
        assert!((ordered - depth_first).together().expand().is_zero());
        assert_eq!(
            engine.cache_statistics(1).unwrap().native.backing_hits,
            before.native.backing_hits + 1
        );
    }
}
