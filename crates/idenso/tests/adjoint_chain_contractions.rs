//! Adjoint generator lines with every colour index summed except the two
//! line endpoints reduce to a group invariant times the adjoint metric.
//!
//! The local UV parity tests in GammaLoop (`factorized_growth_matches_*`)
//! compare two forms of the same integrated poles; which form a run produces
//! depends on the symbol registration order of the process. These are the
//! colour structures of the form that failed to reduce.

use idenso::color::ColorSimplifySettings;
use idenso::dirac::GammaSimplifySettings;
use idenso::tensor::{AlgebraSettings, SymbolicTensor};
use idenso::{CookMode, CookSettings};
use symbolica::atom::{Atom, AtomCore};

fn initialize() {
    idenso::representations::initialize();
    let _ = *idenso::epsilon::EPSILON_SYMBOL;
}

fn parse(source: &str) -> Atom {
    Atom::parse(
        source,
        "spenso",
        symbolica::parser::ParseSettings::symbolica(),
    )
    .unwrap()
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
    let text = expression.to_string();
    !text.contains("f(") && !text.contains("chain(") && !text.contains("trace(")
}

/// All orderings of the multiset `labels`.
fn arrangements<'a>(labels: &[&'a str]) -> Vec<Vec<&'a str>> {
    let mut out = std::collections::BTreeSet::new();
    let n = labels.len();
    let mut indices: Vec<usize> = (0..n).collect();
    // Heap's algorithm over positions.
    fn heap<'a>(
        k: usize,
        indices: &mut Vec<usize>,
        labels: &[&'a str],
        out: &mut std::collections::BTreeSet<Vec<&'a str>>,
    ) {
        if k == 1 {
            out.insert(indices.iter().map(|&i| labels[i]).collect());
            return;
        }
        heap(k - 1, indices, labels, out);
        for i in 0..k - 1 {
            if k % 2 == 0 {
                indices.swap(i, k - 1);
            } else {
                indices.swap(0, k - 1);
            }
            heap(k - 1, indices, labels, out);
        }
    }
    heap(n, &mut indices, labels, &mut out);
    out.into_iter().collect()
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
