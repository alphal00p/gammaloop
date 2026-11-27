//! Compact adjoint lines and their explicit structure constants must agree
//! under changes of dummy-label order.

use idenso::color::{CS, ColorSimplifySettings};
use idenso::dirac::GammaSimplifySettings;
use idenso::epsilon::EPSILON_SYMBOL;
use idenso::representations;
use idenso::tensor::{AlgebraSettings, SymbolicTensor};
use idenso::{CookMode, CookSettings};
use spenso::network::tags::SPENSO_TAG;
use std::collections::BTreeSet;
use symbolica::atom::{Atom, AtomCore};
use symbolica::parser::ParseSettings;

fn initialize() {
    representations::initialize();
    let _ = *EPSILON_SYMBOL;
}

fn parse(source: &str) -> Atom {
    Atom::parse(source, "spenso", ParseSettings::symbolica()).unwrap()
}

fn simplify(expression: &Atom) -> Atom {
    let cooking = CookSettings::indices()
        .with_mode(CookMode::ReversibleEncoding)
        .with_representation_payloads(true, true);
    let tensor = SymbolicTensor::infer_indexed(cooking.try_cook(expression.as_view()).unwrap())
        .unwrap()
        .simplify_algebra(&AlgebraSettings {
            gamma: Some(GammaSimplifySettings::canonical()),
            color: Some(ColorSimplifySettings::default()),
            epsilon: true,
            ..Default::default()
        })
        .unwrap()
        .to_dots()
        .unwrap()
        .expanded(None, false)
        .unwrap()
        .canonize()
        .unwrap();
    cooking.uncook(tensor.expression().as_view()).expand()
}

fn generator(label: &str) -> String {
    format!("f(coad(8,{label}),in,out)")
}

fn line(start: &str, end: &str, labels: &[&str]) -> Atom {
    let generators = labels.iter().map(|l| generator(l)).collect::<Vec<_>>();
    parse(&format!(
        "chain(coad(8,{start}),coad(8,{end}),{})",
        generators.join(",")
    ))
}

/// The structure constants of `line`, written out: `F^x_{ij} = f(x, i, j)`.
fn structure_constants(start: &str, end: &str, labels: &[&str], scope: &str) -> Atom {
    let mut product = Atom::one();
    let mut previous = start.to_string();
    for (k, label) in labels.iter().enumerate() {
        let next = if k + 1 == labels.len() {
            end.to_string()
        } else {
            format!("{scope}_{k}")
        };
        product *= parse(&format!(
            "f(coad(8,{label}),coad(8,{previous}),coad(8,{next}))"
        ));
        previous = next;
    }
    product
}

fn reduced(expression: &Atom) -> bool {
    [CS.f, SPENSO_TAG.chain, SPENSO_TAG.trace]
        .into_iter()
        .all(|head| !expression.contains_symbol(head))
}

/// All orderings of the multiset `labels`.
fn arrangements<'a>(labels: &[&'a str]) -> Vec<Vec<&'a str>> {
    let mut words = BTreeSet::from([Vec::new()]);
    for &label in labels {
        let mut next = BTreeSet::new();
        for word in words {
            for position in 0..=word.len() {
                let mut longer = word.clone();
                longer.insert(position, label);
                next.insert(longer);
            }
        }
        words = next;
    }
    words.into_iter().collect()
}

#[test]
fn four_generator_adjoint_lines_with_summed_labels_reduce_to_invariants() {
    initialize();
    let mut failures = Vec::new();
    // Label names chosen so that the summed labels sort before, between and
    // after the endpoints.
    for (a, b, x, y) in [
        ("a", "b", "x", "y"),
        ("x", "y", "a", "b"),
        ("a", "y", "b", "x"),
        ("m", "n", "c", "z"),
    ] {
        let mut cases = Vec::new();
        // Both summed labels on the line.
        for labels in arrangements(&[x, x, y, y]) {
            cases.push((a, b, labels));
        }
        // One summed label on the line and on an endpoint.
        for labels in arrangements(&[x, x, y, b]) {
            cases.push((a, y, labels.clone()));
            cases.push((y, a, labels));
        }
        for (start, end, labels) in cases {
            let chain = line(start, end, &labels);
            let explicit = structure_constants(start, end, &labels, "internal");
            let from_chain = simplify(&chain);
            let from_explicit = simplify(&explicit);
            if !reduced(&from_chain)
                || !reduced(&from_explicit)
                || !(&from_chain - &from_explicit).expand().is_zero()
            {
                failures.push(format!(
                    "{chain}\n  line -> {from_chain}\n  structure constants -> {from_explicit}"
                ));
            }
        }
    }
    assert!(
        failures.is_empty(),
        "{} unreduced or inconsistent lines:\n{}",
        failures.len(),
        failures.join("\n")
    );
}

#[test]
fn crossed_adjoint_lines_match_the_metric_with_the_expected_coefficient() {
    initialize();
    let ca = CS.cas(2, parse("coad(8)"));
    let metric = parse("g(coad(8,q),coad(8,p))") * ca.pow(Atom::num(2));
    for (lines, divisor) in [
        ([("a", "b", ["q", "c"]), ("c", "p", ["a", "b"])], 2),
        ([("a", "p", ["q", "e"]), ("c", "e", ["a", "c"])], 1),
    ] {
        let mut compact = Atom::one();
        let mut explicit = Atom::one();
        for (i, (start, end, labels)) in lines.into_iter().enumerate() {
            compact *= line(start, end, &labels);
            explicit *= structure_constants(start, end, &labels, &format!("internal_{i}"));
        }
        let expected = &metric / Atom::num(divisor);
        for input in [compact, explicit] {
            assert!((simplify(&input) - &expected).expand().is_zero(), "{input}");
        }
    }
}
