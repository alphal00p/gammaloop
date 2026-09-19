//! Matched-input public-backend timings; no generated or injected reducers.
//!
//! Run this exact ignored test in a fresh process. First family use includes
//! lazy program loading and first point application, not an isolated load
//! timer. Public cache clearing/counters are unavailable: repeats retain the
//! process-local cache and must not be advertised as cold point timings.

use std::collections::HashMap;
use std::time::Instant;

use rustred::persistence::{BinaryIoLimits, ExactTerminalCatalog};
use symbolica::atom::Atom;
use vakint::{
    EvaluationMethod, EvaluationOrder, FMFTOptions, LoopNormalizationFactor, TensorReductionMethod,
    Vakint, VakintSettings,
};

use super::{acceptance_inputs, input::ParentInput, retained_parent_descriptor, test_utils};
use acceptance_inputs::ParentDescriptor;
#[path = "timing_metrics.rs"]
mod metrics;

const PAIRED_SAMPLES: usize = 5;

struct Workload {
    name: String,
    parent: ParentDescriptor,
    input: Atom,
    first_parent_use: bool,
}

impl Workload {
    fn matrix() -> Vec<Self> {
        let mut matrix = Vec::new();
        for parent in [
            ParentDescriptor::H,
            ParentDescriptor::Fg,
            ParentDescriptor::Bmw,
            ParentDescriptor::X,
        ] {
            let descriptor = ParentInput::from_csv(parent.csv());
            let mut powers = (0..10)
                .map(|i| i64::from(i < descriptor.physical_momenta.len()))
                .collect::<Vec<_>>();
            powers[0] = 3;
            matrix.push(Self {
                name: format!("{parent:?}/D1_cubed"),
                parent,
                input: descriptor.integral(&powers),
                first_parent_use: true,
            });
            matrix.push(Self {
                name: format!("{parent:?}/expanded_D7"),
                parent,
                input: descriptor
                    .numerator_equals_propagators(&[6], "D7", "P7")
                    .with_numerator,
                first_parent_use: false,
            });
        }
        let clover = acceptance_inputs::cases()
            .into_iter()
            .find(|case| case.name == "test_integrate_4l_clover")
            .unwrap();
        matrix.push(Self {
            name: "factorized/clover".into(),
            parent: retained_parent_descriptor(&clover.input),
            input: clover.input,
            first_parent_use: false,
        });
        matrix
    }

    fn assert_nonterminal_probe(&self) {
        if !self.first_parent_use {
            return;
        }
        let descriptor = ParentInput::from_csv(self.parent.csv());
        let mut signature = (0..10)
            .map(|i| i64::from(i < descriptor.physical_momenta.len()))
            .collect::<Vec<_>>();
        signature[0] = 3;
        signature.sort_unstable();
        let bytes = match self.parent {
            ParentDescriptor::H => {
                include_bytes!("../../data/rustred/four_loop/h.rrcat.bin").as_slice()
            }
            ParentDescriptor::Fg => {
                include_bytes!("../../data/rustred/four_loop/fg.rrcat.bin").as_slice()
            }
            ParentDescriptor::Bmw => {
                include_bytes!("../../data/rustred/four_loop/bmw.rrcat.bin").as_slice()
            }
            ParentDescriptor::X => {
                include_bytes!("../../data/rustred/four_loop/x.rrcat.bin").as_slice()
            }
        };
        let catalog = ExactTerminalCatalog::decode_generated(
            bytes,
            descriptor.family.fingerprint(),
            10,
            BinaryIoLimits::default(),
        )
        .unwrap();
        // No routing permutation can turn this probe into a declared terminal.
        // Inspect after timing so native import cannot prewarm first-use state.
        for key in catalog.terms().keys() {
            let mut terminal = key.powers().to_vec();
            terminal.sort_unstable();
            assert_ne!(
                terminal, signature,
                "{:?} dotted timing probe is a declared terminal",
                self.parent
            );
        }
    }
}

struct Observation {
    value: Atom,
    wall_ns: u128,
    cpu_seconds: Option<f64>,
    rss_before_kib: Option<u64>,
    rss_after_kib: Option<u64>,
}

impl Observation {
    fn measure(
        vakint: &Vakint,
        settings: &VakintSettings,
        scalar: &Atom,
        ticks: Option<f64>,
    ) -> Self {
        let rss_before_kib = metrics::resident_kib();
        let cpu_before = metrics::cpu_ticks();
        let started = Instant::now();
        let value = vakint
            .evaluate_integral(settings, scalar.as_view())
            .unwrap();
        let wall_ns = started.elapsed().as_nanos();
        let cpu_after = metrics::cpu_ticks();
        let rss_after_kib = metrics::resident_kib();
        let cpu_seconds =
            cpu_before
                .zip(cpu_after)
                .zip(ticks)
                .and_then(|((before, after), ticks)| {
                    after.checked_sub(before).map(|delta| delta as f64 / ticks)
                });
        Self {
            value,
            wall_ns,
            cpu_seconds,
            rss_before_kib,
            rss_after_kib,
        }
    }

    fn report(&self, case: &str, phase: &str, lane: &str, repeat: usize) {
        println!(
            "PUBLIC_SCALAR_TIMING\t{case}\t{phase}\t{lane}\t{repeat}\t{}\t{:?}\t{:?}\t{:?}",
            self.wall_ns, self.cpu_seconds, self.rss_before_kib, self.rss_after_kib
        );
    }
}

#[test]
#[ignore = "five paired public scalar timings; explicit FMFT executable required; fresh process"]
fn representative_public_scalar_timings() {
    test_utils::run_multi_lane_acceptance(|| {
        let oracle = std::env::var("VAKINT_4L_CANDIDATE_ORACLE_FORM_PATH").unwrap();
        let settings = VakintSettings {
            integral_normalization_factor: LoopNormalizationFactor::FMFTandMATAD,
            tensor_reduction_method: TensorReductionMethod::FeynKit,
            number_of_terms_in_epsilon_expansion: 5,
            run_time_decimal_precision: 32,
            use_dot_product_notation: true,
            form_exe_path: "/public-rustred-timing-must-not-use-form".into(),
            evaluation_order: EvaluationOrder::rustred_only(),
            ..VakintSettings::default()
        };
        let mut fmft = settings.clone();
        fmft.form_exe_path = oracle;
        fmft.evaluation_order =
            EvaluationOrder(vec![EvaluationMethod::FMFT(FMFTOptions::default())]);
        let vakint = Vakint::new().unwrap();
        let masses = vakint.params_from_f64(
            &settings,
            &HashMap::from([("muvsq".into(), 1.0), ("mursq".into(), 1.0)]),
        );
        let ticks = metrics::clock_ticks_per_second();
        println!(
            "PUBLIC_SCALAR_TIMING_META\tpaired_samples={PAIRED_SAMPLES}\tclock_ticks={ticks:?}\tCPU=self_plus_waited_children\tRSS=parent_snapshots_not_phase_peaks\tfirst_use=lazy_load_plus_application\tcold_point_cache=unavailable_public_API\tsetup_tensor_numerical_conversion=excluded\tFORM=startup_and_IO_included\tprogram_generation=never"
        );
        for workload in Workload::matrix() {
            assert_eq!(retained_parent_descriptor(&workload.input), workload.parent);
            let canonical = vakint
                .to_canonical(&settings, workload.input.as_view(), false)
                .unwrap();
            let scalar = vakint
                .tensor_reduce(&settings, canonical.as_view())
                .unwrap();
            println!(
                "PUBLIC_SCALAR_START\t{}\t{:?}",
                workload.name, workload.parent
            );
            // The prepared Atom is identical for both methods. No symbolic
            // numerical conversion or comparison belongs to a timed interval.
            let check = |left: &Atom, right: &Atom| {
                let values = [left, right].map(|value| {
                    Vakint::full_numerical_evaluation_without_error(
                        &settings,
                        value.as_view(),
                        &masses,
                        &HashMap::default(),
                        None,
                    )
                    .unwrap()
                });
                let (equal, detail) = values[0].does_approx_match(&values[1], None, 1e-20, 0.0);
                assert!(
                    equal,
                    "{}: numerical FMFT/RustRed mismatch: {detail}",
                    workload.name
                );
            };
            let first = Observation::measure(&vakint, &settings, &scalar, ticks);
            let reference = Observation::measure(&vakint, &fmft, &scalar, ticks);
            check(&first.value, &reference.value);
            let phase = if workload.first_parent_use {
                "first-parent-use"
            } else {
                "first-input-call"
            };
            first.report(&workload.name, phase, "rustred", 0);
            reference.report(&workload.name, "initial", "fmft", 0);
            for repeat in 0..PAIRED_SAMPLES {
                let (native, legacy) = if repeat % 2 == 0 {
                    let native = Observation::measure(&vakint, &settings, &scalar, ticks);
                    let legacy = Observation::measure(&vakint, &fmft, &scalar, ticks);
                    (native, legacy)
                } else {
                    let legacy = Observation::measure(&vakint, &fmft, &scalar, ticks);
                    let native = Observation::measure(&vakint, &settings, &scalar, ticks);
                    (native, legacy)
                };
                check(&native.value, &legacy.value);
                native.report(&workload.name, "repeated-public-call", "rustred", repeat);
                legacy.report(&workload.name, "repeated-public-call", "fmft", repeat);
            }
            workload.assert_nonterminal_probe();
            println!("PUBLIC_SCALAR_PASS\t{}", workload.name);
        }
    });
}
