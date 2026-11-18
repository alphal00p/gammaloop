//! The graph-level entry `evaluate_graph` (FeynKit graph -> value) on the beta
//! thread's one-loop background-field graphs (prepared DOT files): it derives
//! the propagator kind from the graph and equals the internal layer
//! (`from_feynkit` branches through `evaluate_factorized`, each branch's
//! prefactor once), exactly. Data: `VAKINT_TEST_FACTORIZED_DOT_ROOT` (default
//! /common/ben/beta/results/prepared); without it the test prints SKIP.
#![cfg(feature = "feynkit-ingress")]

use std::path::PathBuf;
use std::sync::Arc;

use feynkit_graph::FeynmanDiagram;
use symbolica::atom::{Atom, AtomCore};
use vakint::p_integral::FactorizedNumerator;
use vakint::p_integral::factorized::{
    FromFeynkitOptions, GraphIntegral, evaluate_factorized, evaluate_graph,
};
use vakint::{Vakint, refpath};

#[test]
fn the_graph_entry_derives_the_propagator_kind_and_equals_the_internal_layer() {
    let root = std::env::var("VAKINT_TEST_FACTORIZED_DOT_ROOT")
        .map(PathBuf::from)
        .unwrap_or_else(|_| PathBuf::from("/common/ben/beta/results/prepared"));
    let directory = root.join("L1");
    if !directory.join("model.json").is_file() {
        eprintln!("SKIP: no DOT data under {}", root.display());
        return;
    }
    let vakint = Vakint::new().unwrap();
    let settings = refpath::vakint_settings(2);
    let model = refpath::load_model(&directory).unwrap();
    let options = FromFeynkitOptions::default();
    let mut files = std::fs::read_dir(&directory)
        .unwrap()
        .map(|e| e.unwrap().path())
        .filter(|p| p.extension().is_some_and(|x| x == "dot"))
        .collect::<Vec<_>>();
    files.sort();
    let mut checked = 0;
    for path in files {
        let diagram =
            FeynmanDiagram::from_dot(Arc::clone(&model), &std::fs::read_to_string(&path).unwrap())
                .unwrap();
        let got = evaluate_graph(&vakint, &settings, &diagram, &model, &options).unwrap();
        assert_eq!(
            got.integral,
            GraphIntegral::Propagator,
            "{}",
            path.display()
        );
        let mut expected = Atom::Zero;
        for branch in FactorizedNumerator::from_feynkit(&diagram, &model, &options).unwrap() {
            expected += evaluate_factorized(&vakint, &settings, &branch, &branch.prefactor.total())
                .unwrap()
                .value;
        }
        assert!(
            (&got.value - &expected).expand().is_zero(),
            "{}: graph entry {} vs internal layer {}",
            path.display(),
            got.value,
            expected
        );
        checked += 1;
    }
    eprintln!("graph entry: {checked} L1 graphs equal to the internal layer");
    assert!(checked >= 3);
}
