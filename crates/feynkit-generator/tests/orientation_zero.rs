//! Opt-in closed fermion-loop orientation sums
//! ([`GenerationOptions::cancel_closed_loop_orientations`]).
//!
//! Fixtures in `tests/fixtures/orientation` (numerator attributes stripped,
//! since the test rebuilds every numerator from the model):
//! - `qed_massless.json`: the graph-census QED model (photon, massless electron).
//! - `qed_l5`, `qed_l6`: hepkit `up_to_sign` survivors of the QED photon
//!   self-energy (built-in `Model::qed()`, particle selection `[a, e-]`): the
//!   L5 nf3 pair D1286/D1293 and the 24 L6 nf3 graphs. Each one cancels
//!   against its own loop reversal through an automorphism.
//! - `qcd_l5`: the background-field QCD model of the five-loop beta run and 17
//!   raw L5 source graphs: 12 whose orientation-summed color vanishes (the 6
//!   C-merged zero orbits kept by `generate --reduce`) and 5 nf4 two-triangle
//!   controls whose f-part survives.

use std::{
    collections::BTreeMap,
    sync::{Arc, Mutex},
};

use feynkit_generator::{
    FeynmanDiagramRulesExt, GenerationControl, GenerationFilter, GenerationOptions,
    GenerationResult, GraphGroupingOptions, NumeratorGrouping, OrientationSumVerdict, Process,
};
use feynkit_graph::FeynmanDiagram;
use feynkit_model::Model;
use symbolica::graph::Graph;

const FIXTURES: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/fixtures/orientation");

fn fixture(path: &str) -> String {
    std::fs::read_to_string(format!("{FIXTURES}/{path}")).unwrap()
}

fn qed_massless() -> Model {
    Model::from_json(&fixture("qed_massless.json")).unwrap()
}

fn qcd_background() -> Model {
    Model::from_json(&fixture("qcd_l5/model.json")).unwrap()
}

fn photon_options(loops: usize) -> GenerationOptions {
    GenerationOptions::default()
        .with_loop_count(loops, loops)
        .unwrap()
        .allow_self_loops(true)
        .with_graph_filter(GenerationFilter::LoopOneParticleIrreducible)
}

fn photon_self_energy(model: Model, options: &GenerationOptions) -> GenerationResult {
    Process::new([22_i64], [22_i64])
        .generate_diagrams(model, options)
        .unwrap()
}

fn dot_fixtures(model: &Model, directory: &str) -> Vec<(String, FeynmanDiagram)> {
    let mut paths = std::fs::read_dir(format!("{FIXTURES}/{directory}"))
        .unwrap()
        .map(|entry| entry.unwrap().path())
        .filter(|path| path.extension().is_some_and(|extension| extension == "dot"))
        .collect::<Vec<_>>();
    paths.sort();
    paths
        .into_iter()
        .map(|path| {
            let name = path.file_stem().unwrap().to_string_lossy().into_owned();
            let text = std::fs::read_to_string(&path).unwrap();
            (
                name,
                FeynmanDiagram::from_dot(model.clone(), &text).unwrap(),
            )
        })
        .collect()
}

/// Closed non-ghost fermion loops of a finalized diagram.
fn closed_fermion_loops(diagram: &FeynmanDiagram) -> usize {
    let model = diagram.model();
    let vertices = diagram.vertices().count();
    let mut adjacency = vec![Vec::new(); vertices];
    for (_, ends, edge) in diagram.edges() {
        let particle = model.particle_by_id(edge.particle).unwrap();
        if !particle.is_fermion() || particle.is_ghost() {
            continue;
        }
        if let (Some(source), Some(target)) = (ends.source, ends.target) {
            adjacency[source.0].push(target.0);
            adjacency[target.0].push(source.0);
        }
    }
    let mut seen = vec![false; vertices];
    let mut loops = 0;
    for start in 0..vertices {
        if seen[start] || adjacency[start].is_empty() {
            continue;
        }
        let mut stack = vec![start];
        let mut closed = true;
        while let Some(vertex) = stack.pop() {
            if std::mem::replace(&mut seen[vertex], true) {
                continue;
            }
            closed &= adjacency[vertex].len() == 2;
            stack.extend(adjacency[vertex].iter().copied().filter(|n| !seen[*n]));
        }
        loops += usize::from(closed);
    }
    loops
}

fn report_json(result: &GenerationResult) -> serde_json::Value {
    serde_json::to_value(&result.report).unwrap()
}

fn diagrams_json(result: &GenerationResult) -> BTreeMap<String, String> {
    result
        .diagrams
        .iter()
        .map(|diagram| (diagram.name().to_owned(), diagram.to_json().unwrap()))
        .collect()
}

#[test]
fn orientation_zero_is_opt_in() {
    let stages = Arc::new(Mutex::new(Vec::new()));
    let observed = Arc::clone(&stages);
    let default = photon_self_energy(qed_massless(), &photon_options(3).threads(1));
    let off = photon_self_energy(
        qed_massless(),
        &photon_options(3)
            .threads(1)
            .cancel_closed_loop_orientations(false)
            .progress(move |progress| {
                observed.lock().unwrap().push(progress.stage);
                GenerationControl::Continue
            }),
    );
    assert_eq!(default.report, off.report);
    assert_eq!(report_json(&default), report_json(&off));
    assert_eq!(diagrams_json(&default), diagrams_json(&off));
    assert_eq!(format!("{:?}", default.groups), format!("{:?}", off.groups));
    assert!(!stages.lock().unwrap().contains(&"orientation_sum_zero"));
    let serialized = report_json(&off);
    assert!(serialized.get("orientation_zero_count").is_none());
    assert!(serialized.get("orientation_test_skipped").is_none());
    assert_eq!(off.report.retained_count, 20); // QGRAF onepi

    let stages = Arc::new(Mutex::new(Vec::new()));
    let observed = Arc::clone(&stages);
    let on = photon_self_energy(
        qed_massless(),
        &photon_options(3)
            .threads(1)
            .cancel_closed_loop_orientations(true)
            .progress(move |progress| {
                observed.lock().unwrap().push(progress.stage);
                GenerationControl::Continue
            }),
    );
    assert!(stages.lock().unwrap().contains(&"orientation_sum_zero"));
    assert_eq!(on.report.retained_count, 18); // QGRAF onepi, floop
    assert_eq!(on.report.orientation_zero_count, 2);
    assert_eq!(on.report.orientation_test_skipped, 0);
    assert_eq!(report_json(&on)["orientation_zero_count"], 2);
    // Survivors keep their names, numerators and overall factors.
    let off = diagrams_json(&off);
    for (name, json) in diagrams_json(&on) {
        assert_eq!(off.get(&name), Some(&json), "{name}");
    }
}

#[test]
fn furry_odd_loop_pairs_cancel_without_grouping() {
    let off = photon_self_energy(qed_massless(), &photon_options(4));
    let on = photon_self_energy(
        qed_massless(),
        &photon_options(4).cancel_closed_loop_orientations(true),
    );
    assert_eq!(off.report.retained_count, 189); // QGRAF onepi
    assert_eq!(on.report.retained_count, 153); // QGRAF onepi, floop
    assert_eq!(on.report.orientation_zero_count, 36);
    let per_nf = |result: &GenerationResult| {
        let mut counts = BTreeMap::<usize, usize>::new();
        for diagram in &result.diagrams {
            *counts.entry(closed_fermion_loops(diagram)).or_default() += 1;
        }
        counts
    };
    assert_eq!(per_nf(&on), BTreeMap::from([(1, 105), (2, 45), (3, 3)]));
    // Exactly the diagrams with an odd closed loop are removed. In QED that
    // is Furry's theorem; QGRAF's `floop` drops the same set.
    for diagram in &on.diagrams {
        assert_ne!(
            diagram
                .closed_loop_orientation_sum(&GenerationOptions::default())
                .unwrap(),
            OrientationSumVerdict::Zero
        );
    }
}

#[test]
fn self_conjugate_furry_graph_is_removed() {
    // D1286 and D1293 survive hepkit's up_to_sign grouping together, with
    // ratio +1: each cancels only against its own automorphic reversal.
    let model = Model::qed();
    let fixtures = dot_fixtures(&model, "qed_l5");
    assert_eq!(fixtures.len(), 2);
    for (name, diagram) in &fixtures {
        assert_eq!(closed_fermion_loops(diagram), 3, "{name}");
        assert_eq!(
            diagram
                .closed_loop_orientation_sum(&GenerationOptions::default())
                .unwrap(),
            OrientationSumVerdict::Zero,
            "{name}"
        );
    }
}

#[test]
fn self_conjugate_furry_graphs_at_six_loops_are_removed() {
    let model = Model::qed();
    let fixtures = dot_fixtures(&model, "qed_l6");
    assert_eq!(fixtures.len(), 24);
    for (name, diagram) in &fixtures {
        assert_eq!(
            diagram
                .closed_loop_orientation_sum(&GenerationOptions::default())
                .unwrap(),
            OrientationSumVerdict::Zero,
            "{name}"
        );
    }
}

/// Generated end to end: with `UpToSign`, the pair survives as one group
/// without the flag and disappears with it.
#[test]
fn self_conjugate_furry_graph_is_removed_with_up_to_sign_grouping() {
    let model = Model::qed();
    let ids = dot_fixtures(&model, "qed_l5")
        .into_iter()
        .map(|(_, diagram)| diagram.id())
        .collect::<Vec<_>>();
    let process = Process::new([22_i64], [22_i64])
        .with_particle_filters(None, Some(vec!["a".into(), "e-".into()]), None, vec![])
        .unwrap();
    let options = photon_options(5)
        .with_graph_filter(GenerationFilter::FermionLoopCountRange((3, 3)))
        .numerator_grouping(NumeratorGrouping::UpToSign(GraphGroupingOptions {
            number_of_numerical_samples: 3,
            fully_numerical_substitution: true,
            ..GraphGroupingOptions::default()
        }))
        .select_diagram_ids(ids);
    let off = process.generate_diagrams(model.clone(), &options).unwrap();
    assert_eq!(off.groups.len(), 1);
    assert_eq!(off.groups[0].members.len(), 2);
    let on = process
        .generate_diagrams(model, &options.cancel_closed_loop_orientations(true))
        .unwrap();
    assert!(on.diagrams.is_empty());
    assert_eq!(on.report.orientation_zero_count, 2);
}

#[test]
fn qcd_orientation_summed_colour_zero() {
    let model = qcd_background();
    let fixtures = dot_fixtures(&model, "qcd_l5");
    assert_eq!(fixtures.len(), 17);
    let zeros = [
        178837, 179017, 179305, 179700, 179342, 179737, 179390, 179794, 179396, 179800, 179405,
        179809,
    ];
    let controls = [178677, 179783, 179790, 179823, 180169];
    for (name, diagram) in &fixtures {
        let index: usize = name.trim_start_matches("diagram-").parse().unwrap();
        let verdict = diagram
            .closed_loop_orientation_sum(&GenerationOptions::default())
            .unwrap();
        if zeros.contains(&index) {
            assert_eq!(verdict, OrientationSumVerdict::Zero, "{name}");
        } else {
            assert!(controls.contains(&index), "{name}");
            // nf4 two-triangle graphs: the f-part survives, proportional to C_A.
            assert_eq!(closed_fermion_loops(diagram), 4, "{name}");
            assert_eq!(verdict, OrientationSumVerdict::NonZero, "{name}");
        }
    }
}

#[test]
fn even_loop_bubble_never_dropped() {
    // QED: every one- and two-loop photon self-energy has one even loop.
    for (loops, count) in [(1, 1), (2, 3)] {
        let on = photon_self_energy(
            qed_massless(),
            &photon_options(loops).cancel_closed_loop_orientations(true),
        );
        assert_eq!(on.report.retained_count, count);
        assert_eq!(on.report.orientation_zero_count, 0);
    }
    // Background-field QCD: the one-loop quark bubble has c_sym = c != 0.
    let options = GenerationOptions::default().with_loop_count(1, 1).unwrap();
    let process = Process::new(["B"], ["B"]);
    let off = process
        .generate_diagrams(qcd_background(), &options)
        .unwrap();
    let on = process
        .generate_diagrams(
            qcd_background(),
            &options.clone().cancel_closed_loop_orientations(true),
        )
        .unwrap();
    assert_eq!(off.report.retained_count, on.report.retained_count);
    assert_eq!(on.report.orientation_zero_count, 0);
    let bubbles = on
        .diagrams
        .iter()
        .filter(|diagram| closed_fermion_loops(diagram) == 1)
        .collect::<Vec<_>>();
    assert!(!bubbles.is_empty());
    for bubble in bubbles {
        assert_eq!(
            bubble
                .closed_loop_orientation_sum(&GenerationOptions::default())
                .unwrap(),
            OrientationSumVerdict::NonZero
        );
    }
}

#[test]
fn reversed_propagator_fallback_fails_closed() {
    let mut definition: serde_json::Value =
        serde_json::from_str(&fixture("qed_massless.json")).unwrap();
    definition["propagators"]
        .as_array_mut()
        .unwrap()
        .retain(|propagator| propagator["particle"] != "e+");
    let model = Model::from_json(&definition.to_string()).unwrap();
    let on = photon_self_energy(
        model.clone(),
        &photon_options(3).cancel_closed_loop_orientations(true),
    );
    assert_eq!(on.report.retained_count, 20);
    assert_eq!(on.report.orientation_zero_count, 0);
    assert_eq!(on.report.orientation_test_skipped, 20);
    assert!(matches!(
        on.diagrams[0]
            .closed_loop_orientation_sum(&GenerationOptions::default())
            .unwrap(),
        OrientationSumVerdict::Skipped(_)
    ));
}

/// Undirected skeleton with fixed externals: orientation families share it.
fn undirected_skeleton(diagram: &FeynmanDiagram) -> (String, u64) {
    let model = diagram.model();
    let mut graph = Graph::<i64, i64>::new();
    let vertices = diagram.vertices().count();
    for _ in 0..vertices {
        graph.add_node(0);
    }
    for (_, ends, edge) in diagram.edges() {
        let pdg = model.particle_by_id(edge.particle).unwrap().pdg_code.abs();
        match (ends.source, ends.target, &edge.external) {
            (Some(source), Some(target), None) => {
                graph.add_edge(source.0, target.0, false, pdg).unwrap();
            }
            (Some(vertex), None, Some(external)) | (None, Some(vertex), Some(external)) => {
                let node = graph.add_node(1 + external.index as i64);
                graph.add_edge(node, vertex.0, false, pdg).unwrap();
            }
            other => panic!("unexpected attachment {other:?}"),
        }
    }
    let canonical = graph.canonize();
    (
        format!("{:?}", canonical.graph),
        canonical
            .automorphism_group_size
            .to_string()
            .parse()
            .unwrap(),
    )
}

#[test]
fn weights_of_survivors_unchanged() {
    let off = photon_self_energy(qed_massless(), &photon_options(4));
    let on = photon_self_energy(
        qed_massless(),
        &photon_options(4).cancel_closed_loop_orientations(true),
    );
    let off_factors = off
        .diagrams
        .iter()
        .map(|diagram| (diagram.name().to_owned(), diagram.overall_factor().clone()))
        .collect::<BTreeMap<_, _>>();
    for diagram in &on.diagrams {
        assert_eq!(
            off_factors.get(diagram.name()),
            Some(diagram.overall_factor()),
            "{}",
            diagram.name()
        );
    }
    // The distinct oriented members C of one family satisfy
    // sum_C |Aut(U)|/|Aut(C)| = 2^n; their automorphism factors differ.
    let mut families = BTreeMap::<String, (u64, usize, Vec<u64>)>::new();
    for diagram in &off.diagrams {
        let (key, automorphisms) = undirected_skeleton(diagram);
        let family = families.entry(key).or_insert((
            automorphisms,
            closed_fermion_loops(diagram),
            Vec::new(),
        ));
        family.2.push(diagram.symmetry_factor());
    }
    for (automorphisms, loops, members) in families.values() {
        let total = members
            .iter()
            .map(|member| {
                assert_eq!(automorphisms % member, 0);
                automorphisms / member
            })
            .sum::<u64>();
        assert_eq!(total, 1 << loops);
    }
}

/// Every member of one orientation family (same undirected skeleton) must
/// get the same verdict: the sum from any member is the same family sum up
/// to a sign. A wrong spinor or colour transposition breaks this.
fn assert_family_verdicts_agree(result: &GenerationResult) -> (usize, usize, usize) {
    let mut skipped = 0;
    let mut families = BTreeMap::<String, Vec<(String, OrientationSumVerdict)>>::new();
    for diagram in &result.diagrams {
        let model = diagram.model();
        // Ghost loops are not varied, so ghost reversals are separate families.
        if diagram
            .edges()
            .any(|(_, _, edge)| model.particle_by_id(edge.particle).unwrap().is_ghost())
        {
            continue;
        }
        let verdict = diagram
            .closed_loop_orientation_sum(&GenerationOptions::default())
            .unwrap();
        if matches!(verdict, OrientationSumVerdict::Skipped(_)) {
            skipped += 1;
            continue;
        }
        families
            .entry(undirected_skeleton(diagram).0)
            .or_default()
            .push((diagram.name().to_owned(), verdict));
    }
    let mut zeros = 0;
    for members in families.values() {
        assert!(
            members.iter().all(|(_, verdict)| *verdict == members[0].1),
            "{members:?}"
        );
        zeros += members
            .iter()
            .filter(|(_, verdict)| *verdict == OrientationSumVerdict::Zero)
            .count();
    }
    (families.len(), zeros, skipped)
}

#[test]
fn orientation_family_members_share_verdicts() {
    let qed = photon_self_energy(qed_massless(), &photon_options(4));
    let (_, zeros, skipped) = assert_family_verdicts_agree(&qed);
    assert_eq!((zeros, skipped), (36, 0));
    let qcd = Process::new(["B"], ["B"])
        .generate_diagrams(
            qcd_background(),
            &GenerationOptions::default().with_loop_count(2, 2).unwrap(),
        )
        .unwrap();
    let (families, _, skipped) = assert_family_verdicts_agree(&qcd);
    eprintln!(
        "QCD B->B: {} diagrams, {families} families, {skipped} skipped",
        qcd.diagrams.len()
    );
    assert!(families > 0);
}
