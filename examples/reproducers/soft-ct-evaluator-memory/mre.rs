#!/usr/bin/env rust-script
//! ```cargo
//! [package]
//! edition = "2024"
//! [dependencies]
//! symbolica = { git = "https://github.com/symbolica-dev/symbolica", rev = "578dfcb55fb0662d871456ea8ecc57c1390caa77", default-features = false, features = ["float-mpfr", "integer-gmp"] }
//! [patch.crates-io]
//! numerica = { git = "https://github.com/symbolica-dev/symbolica", rev = "578dfcb55fb0662d871456ea8ecc57c1390caa77" }
//! graphica = { git = "https://github.com/symbolica-dev/symbolica", rev = "578dfcb55fb0662d871456ea8ecc57c1390caa77" }
//! ```
// Usage: rust-script mre.rs [inline|shared|raw] [calls=512] [depth=128] [--verbose]
// Run modes in separate processes when comparing peak memory. Requires a normal
// Symbolica user license; no GammaLoop files, captured expressions, or input data.
//
// H_0(x)=x; H_j(x)=1/(1+x*H_{j-1}(x)); F(u,x)=1/(u+H_depth(x)).
// inline: sum_i F(i,x); shared: sum_i 1/(i+H_depth(x)), i=1..calls.
// Both use the same compact alias chain. Only the ordinary function boundary
// differs. Its fresh per-call translation cache duplicates the shared chain.
//
// raw measures translation before CSE/CPE/stack compaction. At the pinned
// revision, root(0) produces an empty Horner scheme. Explicit coefficient
// collection reproduces its preprocessing before setting horner_iterations(0).
// Numeric checks subsequently compact the raw stack, outside the build timer.
use std::{error::Error, time::Instant};
use symbolica::{
    evaluate::{FunctionMap, FunctionRegistrationOptions, InliningPolicy, OptimizationSettings},
    prelude::*,
};

fn peak_rss_kib() -> Option<u64> {
    std::fs::read_to_string("/proc/self/status")
        .ok()?
        .lines()
        .find_map(|line| {
            line.strip_prefix("VmHWM:")?
                .split_whitespace()
                .next()?
                .parse()
                .ok()
        })
}

fn main() -> Result<(), Box<dyn Error>> {
    let mut args = std::env::args().skip(1).collect::<Vec<_>>();
    let verbose = args.last().is_some_and(|arg| arg == "--verbose");
    if verbose {
        args.pop();
    }
    let mode = args.first().map_or("inline", String::as_str);
    if args.len() > 3 || !["inline", "shared", "raw"].contains(&mode) {
        return Err("usage: mre.rs [inline|shared|raw] [calls=512] [depth=128] [--verbose]".into());
    }
    let calls = args.get(1).map_or(Ok(512), |arg| arg.parse::<usize>())?;
    let depth = args.get(2).map_or(Ok(128), |arg| arg.parse::<usize>())?;
    if calls == 0 || depth == 0 {
        return Err("calls and depth must be positive".into());
    }
    let raw = mode == "raw";
    let preprocess = |atom: Atom| {
        if raw {
            atom.collect_by_coefficient()
        } else {
            atom
        }
    };
    let mut aliases = vec![(parse!("h0"), parse!("x"))];
    for j in 1..=depth {
        aliases.push((
            parse!(&format!("h{j}")),
            preprocess(parse!(&format!("1/(1+x*h{})", j - 1))),
        ));
    }
    let body = preprocess(parse!(&format!("1/(u+h{depth})")));
    let sum = (1..=calls)
        .map(|i| {
            if mode == "shared" {
                format!("1/({i}+h{depth})")
            } else {
                format!("f({i},x)")
            }
        })
        .collect::<Vec<_>>()
        .join("+");
    let sum = preprocess(parse!(&sum));
    let source_bytes = body.as_view().get_byte_size()
        + sum.as_view().get_byte_size()
        + aliases
            .iter()
            .map(|(lhs, rhs)| lhs.as_view().get_byte_size() + rhs.as_view().get_byte_size())
            .sum::<usize>();
    let root = parse!("root(0)");
    aliases.push((root.clone(), sum));
    let mut functions = FunctionMap::new();
    functions.add_function_with_options(
        symbol!("f"),
        vec![symbol!("u"), symbol!("x")],
        body,
        FunctionRegistrationOptions::new().inlining(InliningPolicy::Always),
    )?;
    let settings = OptimizationSettings::new()
        .direct_translation(true)
        .cores(1)
        .horner_iterations(if raw { 0 } else { 1 })
        .cpe_iterations(Some(0))
        .verbose(verbose);
    println!("start mode={mode} calls={calls} depth={depth} source_bytes={source_bytes}");
    let params = vec![parse!("x")];
    let start = Instant::now();
    let mut evaluator = root
        .evaluator(&params)
        .function_map(functions)
        .add_aliases(aliases)?
        .optimization_settings(settings)
        .build()?;
    let seconds = start.elapsed().as_secs_f64();
    let ops = evaluator.count_operations();
    println!(
        "built seconds={seconds:.6} peak_rss_kib={:?} additions={} multiplications={} inversions={} function_calls={}",
        peak_rss_kib(),
        ops.additions,
        ops.multiplications,
        ops.inversions,
        ops.function_calls,
    );
    if raw {
        evaluator.optimize_stack();
    }
    let mut evaluator = evaluator.map_coeff(&|c| c.re.to_f64());
    for x in [0.1_f64, 0.2, 0.3] {
        let mut h = x;
        for _ in 0..depth {
            h = 1.0 / (1.0 + x * h);
        }
        let expected = (1..=calls).map(|i| 1.0 / (i as f64 + h)).sum::<f64>();
        let actual = evaluator.evaluate_single(&[x]);
        let error = (actual - expected).abs() / expected.abs();
        println!(
            "check x={x} actual={actual:.16e} expected={expected:.16e} relative_error={error:.3e}"
        );
        assert!(actual.is_finite() && error < 1e-12, "numeric oracle failed");
    }
    Ok(())
}
