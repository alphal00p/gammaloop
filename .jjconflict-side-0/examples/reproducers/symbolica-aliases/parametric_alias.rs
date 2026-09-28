#!/usr/bin/env rust-script
//! Parametric aliases: an alias handle that carries slot arguments must resolve
//! by pattern, not by literal atom equality.
//!
//! A tensor-valued alias `T(1, mink(4,mu), mink(4,nu))` stands for a definition
//! with two free Minkowski ports. Contraction rewrites the *ports* of such a
//! handle without opening it: `g(mink(4,rho), mink(4,mu)) * T(1, mink(4,mu), mink(4,nu))`
//! becomes `T(1, mink(4,rho), mink(4,nu))`. `AliasedAtom::into_inner` looks the
//! handle up as a literal atom, so the relabeled use stays unresolved, while a
//! one-line pattern replacement resolves every use. Registering one literal
//! alias per labelling works today but duplicates the definition per use.
//!
//! Run standalone with `rust-script parametric_alias.rs` or through the
//! reproducer package with `cargo run --release --bin parametric_alias`.
//!
//! ```cargo
//! [dependencies]
//! symbolica = { git = "https://github.com/symbolica-dev/symbolica", rev = "06906976bca24fefc5203aee699d90d62ebe08cd", default-features = false, features = ["float-mpfr", "integer-gmp", "native_code_generation"] }
//! ```
use symbolica::{
    atom::{AliasedAtom, Atom, AtomCore, AtomView},
    parser::ParseSettings,
};

const NAMESPACE: &str = "mre";

fn atom(source: &str) -> Atom {
    Atom::parse(source, NAMESPACE, ParseSettings::symbolica()).unwrap()
}

fn main() {
    // The definition, with the labelling the engine happened to produce first.
    let handle = atom("T(1, mink(4,mu), mink(4,nu))");
    let body = atom("p(mink(4,mu))*q(mink(4,nu)) + m*g(mink(4,mu), mink(4,nu))");
    let alias_head = match handle.as_view() {
        AtomView::Fun(function) => function.get_symbol(),
        _ => unreachable!("the handle is a function"),
    };
    // Two uses of the same alias: one literal, one after a port was relabeled
    // by a metric contraction. Both mean the definition with the ports renamed.
    let root = atom("T(1, mink(4,mu), mink(4,nu)) + T(1, mink(4,rho), mink(4,nu))");

    let mut aliased = AliasedAtom::from(root.clone());
    aliased.register_alias(handle.clone(), body.clone());
    let resolved = aliased.clone().into_inner();
    println!("root:                 {root}");
    println!("alias:                {handle} = {body}");
    println!("into_inner():         {resolved}");
    println!(
        "relabeled handle left unresolved by literal lookup: {}",
        resolved.contains_symbol(alias_head)
    );

    // What the engine needs: the same alias, resolved by pattern.
    let pattern = atom("T(1, mink(4,a_), mink(4,b_))").to_pattern();
    let rhs = atom("p(mink(4,a_))*q(mink(4,b_)) + m*g(mink(4,a_), mink(4,b_))").to_pattern();
    let by_pattern = root.replace(&pattern).with(&rhs);
    println!("pattern resolution:   {by_pattern}");
    println!(
        "handles left after pattern resolution: {}",
        by_pattern.contains_symbol(alias_head)
    );

    // The workaround available today: one literal alias per distinct labelling.
    // Every labelling repeats the whole definition.
    let labels = ["mu", "rho", "sigma", "tau", "alpha", "beta", "gamma", "delta"];
    let mut per_labelling = AliasedAtom::from(root.clone());
    let slot = atom("mink(4,mu)").to_pattern();
    for label in labels {
        let target = atom(&format!("mink(4,{label})")).to_pattern();
        let handle = handle.replace(&slot).with(&target);
        let body = body.replace(&slot).with(&target);
        per_labelling.register_alias(handle, body);
    }
    println!(
        "bytes: one pattern alias {} | {} literal aliases {} ({}x the definition)",
        aliased.get_byte_size(),
        labels.len(),
        per_labelling.get_byte_size(),
        per_labelling.get_byte_size() / aliased.get_byte_size().max(1)
    );
    println!(
        "per-labelling into_inner() resolves the relabeled use: {}",
        !per_labelling.into_inner().contains_symbol(alias_head)
    );
    println!();
    println!("ask: AliasedAtom::register_pattern_alias(pattern, body) used by into_inner,");
    println!("     get_aliases and EvaluatorBuilder::add_aliases, so a handle with slot");
    println!("     arguments resolves like a function of its ports.");
}
