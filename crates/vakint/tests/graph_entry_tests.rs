//! The graph-level entry `evaluate_graph` (FeynKit graph -> value):
//! - on the beta thread's one-loop background-field graphs (prepared DOT
//!   files): it derives the propagator kind from the graph and equals the
//!   internal layer (`from_feynkit` branches through `evaluate_factorized`,
//!   each branch's prefactor once), exactly. Data:
//!   `VAKINT_TEST_FACTORIZED_DOT_ROOT` (default
//!   /common/ben/beta/results/prepared); without it the test prints SKIP;
//! - on GENERATED scalar vacuum graphs with ONE massive line (feynkit-model's
//!   fixture scalars_2p_3p): it derives the vacuum kind and the massive
//!   closure, and the value equals the graph's constant numerator times the
//!   hand-cut closure, exactly.
#![cfg(feature = "feynkit-ingress")]

use std::path::PathBuf;
use std::sync::Arc;

use feynkit_generator::{GenerationOptions, Process};
use feynkit_graph::FeynmanDiagram;
use feynkit_model::Model;
use symbolica::atom::{Atom, AtomCore};
use vakint::p_integral::factorized::vacuum::contract_to_vakint_vacuum;
use vakint::p_integral::factorized::{
    FromFeynkitOptions, GraphIntegral, ProjectorChoice, evaluate_factorized, evaluate_graph,
};
use vakint::p_integral::{FactorizedNumerator, MassiveClosure};
use vakint::{Vakint, refpath, vakint_parse as vk_parse};

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

/// Generated scalar VACUUM graphs (no external legs) of feynkit-model's fixture
/// scalars_2p_3p (scalar_0 massless, scalar_1 of mass `mass_scalar_1`) at one
/// and two loops with at most two vertices, kept when exactly one line is
/// massive (scalar_1) and the others are massless (scalar_0): the one-loop
/// bubble at zero momentum and the two-loop sunrise with masses (M, 0, 0).
/// `evaluate_graph` alone derives the vacuum kind, the cut and the routing;
/// its value equals SUM over colour branches of the branch prefactor times the
/// graph's (momentum-free) numerator, times vakint's massive closure of the
/// hand-cut graph: `1/p^2` (one loop), the massless bubble (two loops).
#[test]
fn the_graph_entry_integrates_a_vacuum_graph_with_one_massive_line() {
    let model = Arc::new(
        Model::from_json(include_str!(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../feynkit-model/tests/fixtures/scalars_2p_3p.json"
        )))
        .unwrap(),
    );
    let vakint = Vakint::new().unwrap();
    let settings = refpath::vakint_settings(3);
    // a vacuum graph is never projected (evaluate_graph drops the projector itself)
    let options = FromFeynkitOptions {
        projector: ProjectorChoice::None,
        ..FromFeynkitOptions::default()
    };
    let parse = |text: &str| vk_parse!(text).unwrap();
    let mass_squared = parse("UFO::mass_scalar_1^2");
    let cuts = [
        (1usize, "dot(p(1),p(1))^(-1)"),
        (
            2,
            "topo(prop(1,edge(1,2),k(1),0,1)*prop(2,edge(2,1),k(1)-p(1),0,1))",
        ),
    ];
    let mut checked = [0usize; 2];
    for (loops, cut) in cuts {
        let generation = GenerationOptions::default()
            .threads(1)
            .max_vertices(2)
            .with_loop_count(loops, loops)
            .unwrap();
        let generated = Process::new(Vec::<&str>::new(), Vec::<&str>::new())
            .generate_diagrams(Arc::clone(&model), &generation)
            .unwrap();
        let hand = MassiveClosure {
            momentum: parse("p(1)"),
            mass_squared: mass_squared.clone(),
            power: 1,
        }
        .evaluate(&vakint, &settings, parse(cut).as_view())
        .unwrap();
        for graph in &generated.diagrams {
            let diagram =
                FeynmanDiagram::from_dot(Arc::clone(&model), &graph.to_dot().unwrap()).unwrap();
            let branches = FactorizedNumerator::from_feynkit(&diagram, &model, &options).unwrap();
            let internal = branches[0]
                .edges
                .iter()
                .filter(|edge| !edge.is_external())
                .collect::<Vec<_>>();
            let massive = internal
                .iter()
                .filter(|edge| !edge.mass_squared.is_zero())
                .collect::<Vec<_>>();
            if massive.len() != 1
                || massive[0].mass_squared != mass_squared
                || massive[0].power != 1
                || internal.len() != loops + 1
                || internal.iter().any(|edge| edge.power != 1)
            {
                continue;
            }
            let got = evaluate_graph(&vakint, &settings, &diagram, &model, &options).unwrap();
            assert_eq!(got.integral, GraphIntegral::Vacuum);
            let mut constant = Atom::Zero;
            for branch in &branches {
                let (scalar, _) = contract_to_vakint_vacuum(branch).unwrap();
                constant += branch.prefactor.total() * scalar;
            }
            let expected = &constant * &hand;
            assert!(
                (&got.value - &expected).expand().is_zero(),
                "{loops} loops: graph entry {}\nhand closure {}",
                got.value,
                expected
            );
            assert!(!expected.expand().is_zero());
            checked[loops - 1] += 1;
        }
    }
    eprintln!(
        "graph entry, one massive line: {} one-loop and {} two-loop generated graphs equal to the hand closure",
        checked[0], checked[1]
    );
    assert!(checked[0] >= 1 && checked[1] >= 1, "{checked:?}");
}
