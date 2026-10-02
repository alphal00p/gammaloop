#!/usr/bin/env rust-script
//! ```cargo
//! [dependencies]
//! symbolica = "=3.0.0"
//! ```
//!
//! Standalone Symbolica 3.0.0 collect_factors() memory reproducer.
//! Run with SYMBOLICA_LICENSE in the environment:
//!   rust-script collect_factors_mre.rs [TERMS=1024] [WIDTH=64]
//!   rust-script collect_factors_mre.rs 4096 256  # deliberately huge allocations
//!
//! E = sum_i c(i) * (sum_j a(i,j)). There are no common factors to extract.
//! Exactly one collect_factors() call; no expand(), GammaLoop, or fixed-point loop.
//! The outer sum has 2*TERMS distinct factor keys. In Symbolica 3.0.0,
//! inserting absent keys with exponent zero creates 2*TERMS^2 map entries.
//! Each owned inner-sum key is also cloned across the other summands, adding
//! O(TERMS^2 * WIDTH) bytes although the input is only O(TERMS * WIDTH).
//! Linux reports process RSS before/after and its peak; other platforms can
//! use an external memory profiler. Run each size in a fresh process.

use std::{env, error::Error, fs, time::Instant};
use symbolica::{
    atom::{Atom, AtomCore},
    function, symbol,
};

fn main() -> Result<(), Box<dyn Error>> {
    let args: Vec<_> = env::args().skip(1).collect();
    if args.iter().any(|arg| arg == "--help" || arg == "-h") {
        println!("Usage: collect_factors_mre.rs [TERMS=1024] [WIDTH=64]");
        println!("TERMS >= 2, WIDTH >= 2. Stress case: 4096 256.");
        return Ok(());
    }
    if args.len() > 2 {
        return Err("expected at most two arguments: TERMS WIDTH".into());
    }
    let terms = args
        .first()
        .map(|x| x.parse::<usize>())
        .transpose()?
        .unwrap_or(1024);
    let width = args
        .get(1)
        .map(|x| x.parse::<usize>())
        .transpose()?
        .unwrap_or(64);
    if terms < 2 || width < 2 {
        return Err("TERMS and WIDTH must both be at least 2".into());
    }
    let memory = || {
        fs::read_to_string("/proc/self/status")
            .map(|status| {
                status
                    .lines()
                    .filter(|line| line.starts_with("VmRSS:") || line.starts_with("VmHWM:"))
                    .collect::<Vec<_>>()
                    .join("; ")
            })
            .unwrap_or_else(|_| "RSS unavailable on this platform".into())
    };
    let a = symbol!("collect_factors_mre::a");
    let c = symbol!("collect_factors_mre::c");
    let input = Atom::add_many((0..terms).map(|i| {
        let inner = Atom::add_many((0..width).map(|j| function!(a, i, j)));
        function!(c, i) * inner
    }));
    eprintln!(
        "terms={terms} width={width} input_bytes={} dense_map_entries={}",
        input.as_view().get_byte_size(),
        2 * (terms as u128).pow(2),
    );
    eprintln!("before: {}", memory());

    let start = Instant::now();
    let result = input.collect_factors();
    let elapsed = start.elapsed();

    eprintln!("after:  {}", memory());
    println!(
        "collect_seconds={:.6} result_bytes={} unchanged={}",
        elapsed.as_secs_f64(),
        result.as_view().get_byte_size(),
        result == input,
    );
    assert!(result == input, "factor collection changed this expression");
    Ok(())
}
