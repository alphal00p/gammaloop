#!/usr/bin/env rust-script
//! Streaming expansion: visiting the terms a product of sums generates, without
//! building the expanded sum.
//!
//! `expand()` is the only way Symbolica offers to see the terms of a product of
//! sums, and it materializes every one of them as part of one normalized atom.
//! A contraction engine needs to *visit* each generated term (borrowed factor
//! views), contract it and merge it into a small table; the expanded atom is
//! never wanted. This script times `expand()` against a walk over the factor
//! views that touches every generated term without allocating an atom. The
//! walk is what idenso's collector re-implements (tape + distributor) because
//! Symbolica has no streaming `expand`/`replace`.
//!
//! Arguments: `[factors] [terms]`, default 8 factors of 6 terms (6^8 = 1,679,616
//! generated terms, none of which merge).
//!
//! ```cargo
//! [dependencies]
//! symbolica = { git = "https://github.com/symbolica-dev/symbolica", rev = "06906976bca24fefc5203aee699d90d62ebe08cd", default-features = false, features = ["float-mpfr", "integer-gmp", "native_code_generation"] }
//! ```
use std::time::Instant;

use symbolica::{
    atom::{Atom, AtomCore, AtomView},
    coefficient::CoefficientView,
    parser::ParseSettings,
};

fn visit<'a>(
    factors: &[Vec<AtomView<'a>>],
    depth: usize,
    coefficient: i128,
    count: &mut u64,
    checksum: &mut i128,
) {
    if depth == factors.len() {
        *count += 1;
        *checksum += coefficient;
        return;
    }
    for term in &factors[depth] {
        let mut next = coefficient;
        if let AtomView::Mul(product) = term {
            for factor in product.iter() {
                if let AtomView::Num(number) = factor {
                    if let CoefficientView::Natural(numerator, 1, 0, 1) = number.get_coeff_view() {
                        next *= numerator as i128;
                    }
                }
            }
        }
        visit(factors, depth + 1, next, count, checksum);
    }
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let factors: usize = args.get(1).and_then(|a| a.parse().ok()).unwrap_or(8);
    let terms: usize = args.get(2).and_then(|a| a.parse().ok()).unwrap_or(6);

    // Product of sums with distinct variables: every generated term survives.
    let sums: Vec<Atom> = (0..factors)
        .map(|f| {
            let body = (0..terms)
                .map(|t| format!("{}*x{f}_{t}", (t as i64 % 5) - 2 + if t % 5 == 2 { 3 } else { 0 }))
                .collect::<Vec<_>>()
                .join("+");
            Atom::parse(&body, "mre", ParseSettings::symbolica()).unwrap()
        })
        .collect();
    let product = sums.iter().fold(Atom::num(1), |acc, s| &acc * s);
    println!("product of {factors} sums with {terms} terms: {} generated terms", (terms as u64).pow(factors as u32));

    let start = Instant::now();
    let expanded = product.expand();
    let expand_time = start.elapsed();
    println!(
        "expand():   {:>8.3} s, {} terms, {:.1} MB",
        expand_time.as_secs_f64(),
        expanded.nterms(),
        expanded.as_view().get_byte_size() as f64 / 1e6
    );

    let views: Vec<Vec<AtomView<'_>>> = sums
        .iter()
        .map(|s| match s.as_view() {
            AtomView::Add(sum) => sum.iter().collect(),
            other => vec![other],
        })
        .collect();
    let start = Instant::now();
    let (mut count, mut checksum) = (0u64, 0i128);
    visit(&views, 0, 1, &mut count, &mut checksum);
    let walk_time = start.elapsed();
    println!(
        "walk views: {:>8.3} s, {} terms visited, coefficient checksum {}, 0 atoms built",
        walk_time.as_secs_f64(),
        count,
        checksum
    );
    println!(
        "materialization costs {:.0}x the traversal",
        expand_time.as_secs_f64() / walk_time.as_secs_f64().max(1e-9)
    );
    println!();
    println!("ask: an `expand`/`replace` that hands each generated term to a consumer as");
    println!("     borrowed factor views (with Symbolica's coefficient and power folding),");
    println!("     so an engine can contract and merge terms without the expanded atom.");
}
