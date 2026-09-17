//! Opt-in matched finite-target timing after numerical acceptance has passed.
//!
//! Generation and offline catalogs are outside every measured interval. RSS
//! values are parent-process snapshots, not per-lane peaks or FORM child memory.

use std::collections::HashMap;
use std::time::Instant;

use symbolica::atom::Atom;
use vakint::rustred_evaluation::experimental::ExperimentalRustRed;
use vakint::rustred_evaluation::experimental::native::NativeCandidate;
use vakint::{
    EvaluationMethod, EvaluationOrder, FMFTOptions, TensorReductionMethod, Vakint, VakintSettings,
};

fn cpu_ticks() -> Option<u64> {
    let status = std::fs::read_to_string("/proc/self/stat").ok()?;
    // Field two (comm) may contain spaces and parentheses.
    let fields = status
        .rsplit_once(')')?
        .1
        .split_whitespace()
        .collect::<Vec<_>>();
    [11, 12, 13, 14].into_iter().try_fold(0u64, |sum, field| {
        sum.checked_add(fields.get(field)?.parse::<u64>().ok()?)
    })
}

fn resident_kib() -> Option<u64> {
    std::fs::read_to_string("/proc/self/status")
        .ok()?
        .lines()
        .find_map(|line| line.strip_prefix("VmRSS:"))?
        .split_whitespace()
        .next()?
        .parse()
        .ok()
}

fn scalar_input(vakint: &Vakint, settings: &VakintSettings, input: &Atom) -> Atom {
    let canonical = vakint
        .to_canonical(settings, input.as_view(), false)
        .unwrap();
    vakint.tensor_reduce(settings, canonical.as_view()).unwrap()
}

pub fn run<const N: usize>(
    vakint: &Vakint,
    settings: &VakintSettings,
    native: &NativeCandidate<N>,
    evaluator: &ExperimentalRustRed,
    inputs: &[(&str, Atom)],
    repeats: usize,
) {
    assert!(
        (1..=9).contains(&repeats),
        "bounded timing repeats are 1..=9"
    );
    let clock_ticks = std::process::Command::new("getconf")
        .arg("CLK_TCK")
        .output()
        .ok()
        .filter(|output| output.status.success())
        .and_then(|output| String::from_utf8(output.stdout).ok())
        .and_then(|output| output.trim().parse::<f64>().ok())
        .filter(|ticks| ticks.is_finite() && *ticks > 0.0);
    let mut legacy_settings = settings.clone();
    legacy_settings.evaluation_order =
        EvaluationOrder(vec![EvaluationMethod::FMFT(FMFTOptions::default())]);
    let mut native_settings = settings.clone();
    native_settings.form_exe_path = "/timed-candidate-scalar-must-not-use-form".into();
    let mut prepass_settings = native_settings.clone();
    prepass_settings.tensor_reduction_method = TensorReductionMethod::FeynKit;
    let masses = vakint.params_from_f64(
        settings,
        &HashMap::from([("mursq".into(), 1.0), ("muvsq".into(), 1.0)]),
    );
    println!(
        "EVALUATION_TIMING_META\trepeats={repeats}\tclock_ticks={clock_ticks:?}\tRSS=parent_snapshots_not_phase_peaks\tFORM=startup_and_IO_included\tsetup=excluded"
    );
    let mut measured_cold_applications = 0usize;
    for (name, input) in inputs {
        let scalar = scalar_input(vakint, &prepass_settings, input);
        for phase in ["scalar", "full-feynkit"] {
            for repeat in 0..repeats {
                // Alternate which implementation runs first. The native warm
                // observation immediately follows its cold observation.
                let lanes = if repeat % 2 == 0 {
                    ["fmft", "rustred-cold", "rustred-warm"]
                } else {
                    ["rustred-cold", "rustred-warm", "fmft"]
                };
                let mut results = Vec::new();
                for lane in lanes {
                    if lane == "rustred-cold" {
                        native.clear_cache().unwrap();
                    }
                    let before_rules = native.rule_applications().unwrap();
                    let before_cpu = cpu_ticks();
                    let before_rss = resident_kib();
                    let started = Instant::now();
                    let prepared = if phase == "full-feynkit" {
                        scalar_input(vakint, &prepass_settings, input)
                    } else {
                        scalar.clone()
                    };
                    let value = if lane == "fmft" {
                        vakint
                            .evaluate_integral(&legacy_settings, prepared.as_view())
                            .unwrap()
                    } else {
                        evaluator
                            .evaluate_integral(vakint, &native_settings, prepared.as_view())
                            .unwrap()
                            .value
                    };
                    let elapsed = started.elapsed();
                    let after_cpu = cpu_ticks();
                    let after_rss = resident_kib();
                    let applications = native.rule_applications().unwrap() - before_rules;
                    if lane == "rustred-warm" {
                        assert_eq!(
                            applications, 0,
                            "warm observation must reuse the point cache"
                        );
                    } else if lane == "rustred-cold" {
                        measured_cold_applications += applications;
                    }
                    let cpu_seconds = before_cpu.zip(after_cpu).zip(clock_ticks).and_then(
                        |((before, after), ticks)| {
                            after.checked_sub(before).map(|delta| delta as f64 / ticks)
                        },
                    );
                    results.push((
                        lane,
                        value,
                        elapsed.as_nanos(),
                        cpu_seconds,
                        before_rss,
                        after_rss,
                        applications,
                    ));
                }
                // Numerical conversion and the unchanged strict comparison
                // are outside all timers. No measured result is trusted just
                // because an earlier correctness run passed.
                let numerical = results
                    .iter()
                    .map(|(_, value, ..)| {
                        Vakint::full_numerical_evaluation_without_error(
                            settings,
                            value.as_view(),
                            &masses,
                            &HashMap::default(),
                            None,
                        )
                        .unwrap()
                    })
                    .collect::<Vec<_>>();
                let reference = results
                    .iter()
                    .position(|(lane, ..)| *lane == "fmft")
                    .unwrap();
                for value in &numerical {
                    let (matches, detail) =
                        numerical[reference].does_approx_match(value, None, 1e-20, 10.0);
                    assert!(matches, "timed result differs from FMFT: {detail}");
                }
                for (lane, _, wall_ns, cpu, before_rss, after_rss, applications) in results {
                    println!(
                        "EVALUATION_TIMING\t{name}\t{phase}\t{lane}\t{repeat}\t{wall_ns}\t{cpu:?}\t{before_rss:?}\t{after_rss:?}\t{applications}"
                    );
                }
            }
        }
    }
    assert!(
        measured_cold_applications > 0,
        "timing must include genuine cold recurrence application"
    );
}
