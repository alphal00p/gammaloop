#!/usr/bin/env rust-script
//! Evaluators from aliased atoms: Rust has them, Python does not.
//!
//! A contraction result is a DAG: every node is an alias whose definition
//! refers to earlier aliases. `AliasedAtom::evaluator_multiple` builds an
//! evaluator straight from the alias table. The only alternative is to resolve
//! the aliases into one tree first, which duplicates every shared node; a chain
//! in which each definition uses its predecessor twice (`w/2 + w^2/4`) doubles
//! in size per level. Python exposes only the tree route (`Expression.evaluator` has no
//! alias input), which is what `gammalooprs` would have to use to consume a
//! contraction result without expanding it.
//!
//! Argument: chain depth, default 18.
//!
//! ```cargo
//! [dependencies]
//! symbolica = { git = "https://github.com/symbolica-dev/symbolica", rev = "06906976bca24fefc5203aee699d90d62ebe08cd", default-features = false, features = ["float-mpfr", "integer-gmp", "native_code_generation"] }
//! ```
use std::time::Instant;

use symbolica::{
    atom::{AliasedAtom, Atom, AtomCore},
    parser::ParseSettings,
};

fn atom(source: &str) -> Atom {
    Atom::parse(source, "mre", ParseSettings::symbolica()).unwrap()
}

fn main() {
    let depth: usize = std::env::args().nth(1).and_then(|a| a.parse().ok()).unwrap_or(18);
    let x = atom("x");

    // w(0) = x, w(i) = w(i-1)/2 + w(i-1)^2/4: each level uses the previous alias
    // twice, so the resolved tree doubles per level while the value stays finite.
    let mut aliased = AliasedAtom::from(atom(&format!("w({depth})")));
    aliased.register_alias(atom("w(0)"), x.clone());
    for i in 1..=depth {
        let handle = atom(&format!("w({i})"));
        let body = atom(&format!("w({})/2 + w({})^2/4", i - 1, i - 1));
        aliased.register_alias(handle, body);
    }
    println!("alias table: {} definitions, {} bytes", aliased.get_aliases().len(), aliased.get_byte_size());

    let start = Instant::now();
    let mut from_aliases = AliasedAtom::evaluator_multiple(std::slice::from_ref(&aliased), &[x.clone()])
        .unwrap()
        .build()
        .unwrap()
        .map_coeff(&|c| c.re.to_f64());
    let alias_build = start.elapsed();
    let value_from_aliases = from_aliases.evaluate_single(&[0.5]);
    println!(
        "evaluator from aliases:  built in {:.3} s, operations {:?}",
        alias_build.as_secs_f64(),
        from_aliases.count_operations()
    );

    let start = Instant::now();
    let tree = aliased.clone().into_inner();
    let resolve_time = start.elapsed();
    println!(
        "into_inner():            {:.3} s, {:.1} MB tree ({} bytes per definition otherwise)",
        resolve_time.as_secs_f64(),
        tree.as_view().get_byte_size() as f64 / 1e6,
        aliased.get_byte_size() / (depth + 1)
    );
    let start = Instant::now();
    let mut from_tree = tree
        .evaluator(&[x.clone()])
        .build()
        .unwrap()
        .map_coeff(&|c| c.re.to_f64());
    let tree_build = start.elapsed();
    let value_from_tree = from_tree.evaluate_single(&[0.5]);
    println!(
        "evaluator from the tree: built in {:.3} s, operations {:?}",
        tree_build.as_secs_f64(),
        from_tree.count_operations()
    );
    println!(
        "same value: {} ({value_from_aliases:e})",
        (value_from_aliases - value_from_tree).abs() <= 1e-9 * value_from_aliases.abs().max(1.0)
    );
    println!();
    println!("ask: expose the alias route in Python — an `AliasedExpression` (root + alias");
    println!("     pairs) with `evaluator(...)`, or an `aliases=` input on `Expression.evaluator`,");
    println!("     so results are consumed as DAGs without resolving them to trees.");
}
