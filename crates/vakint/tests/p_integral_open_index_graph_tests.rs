//! D6 (helper VH) on GRAPHS: the background-field `B -> B` self-energy graphs
//! of one and two loops (the beta thread's prepared FeynKit DOT files) with
//! the two external-leg indices left OPEN instead of contracted with a
//! projector. The open indices are detected (the indices that occur once in
//! the product of the graph's numerators), the external momentum and the
//! routing come from the graph (`VakintExpression::from_diagram`), and the
//! open tensor result `Π^{μν}` of vakint + RustRed is closed AFTERWARDS with
//! the three projections of the reference path (lane refpath):
//!
//! * metric:        `g_{μν} Π^{μν} = (D-1) 𝑖 q² × refpath's metric value`,
//! * transverse:    `(g_{μν} - p_μ p_ν/q²) Π^{μν} = (D-1) 𝑖 q² × refpath's transverse value`,
//! * longitudinal:  `p_μ p_ν Π^{μν} = 𝑖 q⁴ × refpath's longitudinal value`,
//!
//! order by order in `ε` through the last order of the reference value (the
//! polynomial factors `D - 1`, `q²` do not mix truncated orders), exactly.
//!
//! Data: `VAKINT_TEST_FACTORIZED_DOT_ROOT` (default
//! `/common/ben/beta/results/prepared`); without it the test prints SKIP and
//! passes. The built-in 1L-3L p library; no FORM.
#![cfg(feature = "feynkit-ingress")]

#[path = "support/open_indices.rs"]
mod open_indices;

use std::collections::BTreeMap;
use std::path::{Path, PathBuf};
use std::sync::Arc;

use feynkit_graph::FeynmanDiagram;
use feynkit_kinematics::Kinematics;
use open_indices::*;
use symbolica::atom::{Atom, AtomCore, AtomView};
use vakint::diagram_integral::DiagramIntegralOptions;
use vakint::refpath::{self, Projection, RunOptions};
use vakint::{Vakint, VakintExpression, vakint_parse, vakint_symbol};

const DEFAULT_ROOT: &str = "/common/ben/beta/results/prepared";

fn root() -> Option<PathBuf> {
    let root = std::env::var("VAKINT_TEST_FACTORIZED_DOT_ROOT")
        .map(PathBuf::from)
        .unwrap_or_else(|_| PathBuf::from(DEFAULT_ROOT));
    if root.join("L1").join("model.json").is_file() {
        Some(root)
    } else {
        eprintln!("SKIP: no DOT data under {}", root.display());
        None
    }
}

fn dots(dir: &Path) -> Vec<PathBuf> {
    let mut files = std::fs::read_dir(dir)
        .unwrap()
        .map(|entry| entry.unwrap().path())
        .filter(|path| {
            path.file_name()
                .and_then(|name| name.to_str())
                .is_some_and(|name| name.starts_with("diagram-") && name.ends_with(".dot"))
        })
        .collect::<Vec<_>>();
    files.sort();
    files
}

fn graph_identity(engine: &Vakint, model: &Arc<feynkit_model::Model>, path: &Path) -> &'static str {
    let dot = std::fs::read_to_string(path).unwrap();
    let diagram = FeynmanDiagram::from_dot(Arc::clone(model), &dot).unwrap();
    if refpath::scaleless_contact(&diagram) {
        return "scaleless";
    }
    let (branches, _) = refpath::colour_branches(&diagram).unwrap();
    if branches.is_empty() {
        return "colour zero";
    }
    let loops = diagram.loop_count();
    let settings = refpath::vakint_settings(loops as i64 + 1);
    let epsilon = vakint_symbol!(settings.epsilon_symbol.as_str());
    let in_epsilon = Atom::num(4) - Atom::num(2) * Atom::var(epsilon);

    // The reference values (refpath, projector applied before vakint).
    let reference = |projection| {
        refpath::run_graph(
            engine,
            model,
            &dot,
            &RunOptions {
                projection,
                count_terms: false,
                ..RunOptions::default()
            },
        )
        .unwrap()
        .results
        .into_iter()
        .map(|result| (result.projection, result.value))
        .collect::<BTreeMap<_, _>>()
    };
    let mut references = reference(Projection::Metric);
    references.extend(reference(Projection::Both));

    // The open path: ports detected, left open, mapped to the indices 1, 2.
    let (open, ports) = refpath::contract_open(&diagram, &branches).unwrap();
    assert_eq!(ports.len(), 2, "{}: open ports {ports:?}", path.display());
    let mut open = open;
    for (port, label) in ports.iter().zip([1, 2]) {
        let AtomView::Fun(slot) = port.as_view() else {
            panic!("port {port}")
        };
        let index = slot.iter().nth(1).unwrap().to_owned();
        open = open.replace(index.to_pattern()).with(Atom::num(label).to_pattern());
    }
    let weight = feynkit_graph::expressions::evaluate_overall_factor(diagram.overall_factor().as_view())
        * Atom::var(refpath::SYMBOLS.nf).pow(Atom::num(refpath::closed_quark_loops(&diagram).unwrap() as i64));
    let dimension_symbol = refpath::SYMBOLS.d;
    let kinematics = Kinematics::in_dimension(&Atom::var(dimension_symbol)).unwrap();
    let family = diagram.propagator_family(&kinematics).unwrap();
    let input = Atom::from(
        VakintExpression::from_diagram(&diagram, &family, &(open * &weight), &DiagramIntegralOptions::default())
            .unwrap_or_else(|error| panic!("{}: {error}", path.display())),
    )
    .replace(Atom::var(dimension_symbol))
    .with(in_epsilon.clone());
    let value = engine
        .evaluate(&settings, input.as_view())
        .unwrap_or_else(|error| panic!("{}: {error}", path.display()));
    let rendered = value.to_canonical_string();
    assert!(!rendered.contains("topo("), "{}: unreduced {rendered}", path.display());

    let q2 = vakint_parse!("dot(p(1),p(1))").unwrap();
    let i = Atom::i();
    let d_minus_1 = &in_epsilon - Atom::num(1);
    let g = vakint_parse!("g(1,2)").unwrap();
    let pp = vakint_parse!("p(1,1)*p(1,2)").unwrap();
    for (name, closing, factor) in [
        ("metric", g.clone(), &d_minus_1 * &i * &q2),
        ("transverse", &g - &pp / &q2, &d_minus_1 * &i * &q2),
        ("longitudinal", pp.clone(), &i * &q2 * &q2),
    ] {
        let reference = references.get(name).cloned().unwrap_or(Atom::Zero);
        let closed = close_indices(&(&value * &closing), epsilon);
        assert_series_through(
            &closed,
            &(&factor * &reference),
            last_order(&reference, &closed, epsilon),
            &format!("{} {name}", path.display()),
            epsilon,
        );
    }
    "checked"
}

fn run_set(set: &str) {
    let Some(root) = root() else {
        return;
    };
    let directory = root.join(set);
    let model = refpath::load_model(&directory).unwrap();
    let engine = Vakint::new().unwrap();
    let mut counts = BTreeMap::<&str, usize>::new();
    for path in dots(&directory) {
        *counts.entry(graph_identity(&engine, &model, &path)).or_default() += 1;
    }
    eprintln!("D6 graph identity {set}: {counts:?}");
    assert!(counts.get("checked").copied().unwrap_or(0) > 0);
}

#[test]
fn open_index_self_energies_close_to_the_reference_at_one_loop() {
    run_set("L1");
}

#[test]
fn open_index_self_energies_close_to_the_reference_at_two_loops() {
    run_set("L2");
}
