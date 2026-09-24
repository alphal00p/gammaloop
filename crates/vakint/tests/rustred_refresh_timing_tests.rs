//! Matched public scalar timings for packaged one-, two-, and three-loop rules.
//!
//! Invoke one exact ignored test per fresh process. The first native call
//! includes lazy loading and application; subsequent calls retain public caches.
//! No public cache-reset API exists, so repeats are not cold-point measurements.
mod test_utils;

#[path = "experimental_rustred_4l/timing_metrics.rs"]
mod metrics;

use std::collections::HashMap;
use symbolica::atom::Atom;
use vakint::{
    EvaluationOrder, LoopNormalizationFactor, TensorReductionMethod, Vakint, VakintSettings,
    vakint_parse,
};

const WARM_REPEATS: usize = 31;

fn scalar_timings(loops: i64, inputs: [(&'static str, &'static str); 2]) {
    test_utils::run_multi_lane_acceptance(move || {
        let settings = VakintSettings {
            integral_normalization_factor: LoopNormalizationFactor::FMFTandMATAD,
            tensor_reduction_method: TensorReductionMethod::FeynKit,
            number_of_terms_in_epsilon_expansion: loops + 1,
            run_time_decimal_precision: 32,
            use_dot_product_notation: true,
            form_exe_path: "/rustred-refresh-timing-must-not-use-form".into(),
            evaluation_order: EvaluationOrder::rustred_only(),
            ..VakintSettings::default()
        };
        let mut matad = settings.clone();
        matad.form_exe_path = std::env::var("VAKINT_REFRESH_ORACLE_FORM_PATH")
            .expect("explicit MATAD FORM executable required");
        matad.evaluation_order = EvaluationOrder::matad_only(None);
        let vakint = Vakint::new().unwrap();
        let masses = vakint.params_from_f64(
            &settings,
            &HashMap::from([("muvsq".into(), 1.0), ("mursq".into(), 1.0)]),
        );
        let ticks = metrics::clock_ticks_per_second();
        println!(
            "PUBLIC_SCALAR_TIMING_META\tloops={loops}\twarm_repeats={WARM_REPEATS}\tclock_ticks={ticks:?}\tCPU=self_plus_waited_children\tRSS=parent_snapshots_not_phase_peaks\tfirst_use=lazy_load_plus_application\tcold_point_cache=unavailable_public_API\tsetup_tensor_numerical_conversion=excluded\tprogram_generation=never\tprecision=32\tterms={}\tsymmetric_relative_threshold=1e-20\terror_estimate=none",
            loops + 1
        );
        for (index, (name, expression)) in inputs.into_iter().enumerate() {
            let input = vakint_parse!(expression).unwrap();
            let canonical = vakint
                .to_canonical(&settings, input.as_view(), false)
                .unwrap();
            let scalar = vakint
                .tensor_reduce(&settings, canonical.as_view())
                .unwrap();
            // Keep the prepared scalar Atom identical. Native first-use runs
            // before MATAD or numerical conversion can register shared state.
            let first = metrics::Observation::measure(&vakint, &settings, &scalar, ticks);
            let reference = metrics::Observation::measure(&vakint, &matad, &scalar, ticks);
            let numerical = |value: &Atom| {
                let result = Vakint::full_numerical_evaluation_without_error(
                    &settings,
                    value.as_view(),
                    &masses,
                    &HashMap::default(),
                    None,
                )
                .unwrap();
                assert!(!result.0.is_empty(), "{name}: empty numerical output");
                assert!(
                    result
                        .0
                        .iter()
                        .all(|(_, value)| value.re.is_finite() && value.im.is_finite()),
                    "{name}: non-finite numerical output"
                );
                result
            };
            let oracle = numerical(&reference.value);
            let check = |value: &Atom| {
                let (equal, detail) = numerical(value).does_approx_match(&oracle, None, 1e-20, 0.0);
                assert!(equal, "{name}: MATAD/RustRed numerical mismatch: {detail}");
            };
            check(&first.value);
            first.report(
                name,
                if index == 0 {
                    "first-family-use"
                } else {
                    "first-input-call"
                },
                "rustred",
                0,
            );
            reference.report(name, "initial", "matad", 0);
            for repeat in 0..WARM_REPEATS {
                let repeated = metrics::Observation::measure(&vakint, &settings, &scalar, ticks);
                check(&repeated.value);
                repeated.report(name, "repeated-public-call", "rustred", repeat);
            }
            println!("PUBLIC_SCALAR_PASS\t{name}");
        }
    });
}

#[test]
#[ignore = "public 1L first-use/warm timings and MATAD parity; fresh process and explicit FORM required"]
fn one_loop_public_scalar_timings() {
    scalar_timings(
        1,
        [
            ("1L/D4", "topo(I1L(muvsq,4))"),
            ("1L/D6", "topo(I1L(muvsq,6))"),
        ],
    );
}

#[test]
#[ignore = "public 2L first-use/warm timings and MATAD parity; fresh process and explicit FORM required"]
fn two_loop_public_scalar_timings() {
    scalar_timings(
        2,
        [
            ("2L/D1_cubed", "topo(I2L(muvsq,3,1,1))"),
            ("2L/numerator", "topo(I2L(muvsq,-1,2,2))"),
        ],
    );
}

#[test]
#[ignore = "public 3L first-use/warm timings and MATAD parity; fresh process and explicit FORM required"]
fn three_loop_public_scalar_timings() {
    scalar_timings(
        3,
        [
            ("3L/D1_squared", "topo(I3L(muvsq,2,1,1,1,1,1))"),
            ("3L/pinch6", "topo(I3L_pinch_6(muvsq,1,1,1,1,1,0))"),
        ],
    );
}
