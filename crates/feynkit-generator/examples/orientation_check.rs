//! Photon self-energy counts with and without the closed-loop orientation test.
//!
//! ```text
//! cargo run --release -p feynkit-generator --example orientation_check -- \
//!     MODEL.json LOOPS [none|sign] [on|off] [THREADS] [NF]
//! ```
//!
//! Prints one JSON object: retained diagrams per number of closed fermion
//! loops (nf), the report counters and the wall time. The model must contain
//! the photon (PDG 22). `MODEL.json` may be `builtin:qed`: the built-in QED
//! model with massive electrons, restricted to `[a, e-]` as in the hepkit runs.
//! With `NF`, only diagrams with exactly that many closed fermion loops are
//! generated. With `ORIENTATION_CHECK_GROUPS=PATH`, the groups are written to
//! PATH as `master<TAB>member<TAB>ratio` lines in result order.
//!
//! Further environment switches:
//! - `ORIENTATION_CHECK_PROCESS=NAME`: the self-energy of particle `NAME`
//!   (for example `B` in the background-field QCD model) instead of the photon;
//! - `ORIENTATION_CHECK_FILTER=none`: no one-particle-irreducibility filter;
//! - `ORIENTATION_CHECK_SELF_LOOPS=off`: no self-loops (on by default);
//! - `ORIENTATION_CHECK_FILTER_ZERO_COLOR=on`: also drop diagrams whose own
//!   colour vanishes (`GenerationOptions::filter_zero_color`);
//! - `ORIENTATION_CHECK_BETA=on`: beta's self-energy options (unit projector,
//!   no tadpoles, no internal background line);
//! - `ORIENTATION_CHECK_MATCH=DIR`: report, per DOT file in DIR, the retained
//!   diagram with the same `DiagramId` (on stderr);
//! - `ORIENTATION_CHECK_NAMES=PATH`: write the retained diagram names and ids;
//! - `ORIENTATION_CHECK_VERDICTS=PATH`: write `name<TAB>nf<TAB>verdict` for every
//!   retained diagram, the verdict recomputed per diagram (in parallel);
//!   with `ORIENTATION_CHECK_DOT_DIR=DIR` also the DOT of every `Zero` and
//!   `Skipped` diagram.

use std::{collections::BTreeMap, time::Instant};

use rayon::prelude::*;

use feynkit_generator::{
    FeynmanDiagramRulesExt, GenerationFilter, GenerationOptions, GraphGroupingOptions,
    NumeratorGrouping, OrientationSumVerdict, Process, TadpoleFilterOptions,
};
use feynkit_graph::FeynmanDiagram;
use feynkit_model::Model;
use symbolica::atom::Atom;

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

fn main() {
    let args = std::env::args().collect::<Vec<_>>();
    let (model, process) = if args[1] == "builtin:qed" {
        let process = Process::new([22_i64], [22_i64])
            .with_particle_filters(None, Some(vec!["a".into(), "e-".into()]), None, vec![])
            .unwrap();
        (Model::qed(), process)
    } else {
        let model = Model::from_json(&std::fs::read_to_string(&args[1]).unwrap()).unwrap();
        let process = match std::env::var("ORIENTATION_CHECK_PROCESS") {
            Ok(name) => Process::new([name.as_str()], [name.as_str()]),
            Err(_) => Process::new([22_i64], [22_i64]),
        };
        (model, process)
    };
    let loops: usize = args[2].parse().unwrap();
    let grouping = match args.get(3).map(String::as_str).unwrap_or("none") {
        "none" => NumeratorGrouping::None,
        "sign" => NumeratorGrouping::UpToSign(GraphGroupingOptions {
            number_of_numerical_samples: 3,
            fully_numerical_substitution: true,
            ..GraphGroupingOptions::default()
        }),
        other => panic!("unknown grouping {other}"),
    };
    let flag = args.get(4).map(String::as_str).unwrap_or("on") == "on";
    let threads: usize = args.get(5).map_or(8, |value| value.parse().unwrap());
    let mut options = GenerationOptions::default()
        .with_loop_count(loops, loops)
        .unwrap()
        .allow_self_loops(std::env::var("ORIENTATION_CHECK_SELF_LOOPS").as_deref() != Ok("off"))
        .threads(threads)
        .numerator_grouping(grouping)
        .cancel_closed_loop_orientations(flag);
    if std::env::var("ORIENTATION_CHECK_FILTER_ZERO_COLOR").as_deref() == Ok("on") {
        options = options.filter_zero_color(true);
    }
    if std::env::var("ORIENTATION_CHECK_FILTER").as_deref() != Ok("none") {
        options = options.with_graph_filter(GenerationFilter::LoopOneParticleIrreducible);
    }
    if let Some(nf) = args.get(6) {
        let nf: usize = nf.parse().unwrap();
        options = options.with_graph_filter(GenerationFilter::FermionLoopCountRange((nf, nf)));
    }
    if std::env::var("ORIENTATION_CHECK_BETA").as_deref() == Ok("on") {
        // beta's `self_energy_options`: unit projector, no tadpoles and no
        // background line between two interaction vertices.
        let background = model.particle_id("B").unwrap();
        options = options
            .projector(Atom::one())
            .with_graph_filter(GenerationFilter::Tadpoles(TadpoleFilterOptions::default()))
            .filter(move |graph, _| {
                graph.edges().iter().all(|edge| {
                    edge.data.particle != background
                        || graph.nodes()[edge.vertices.0].data.external.is_some()
                        || graph.nodes()[edge.vertices.1].data.external.is_some()
                })
            });
    }
    let fixtures = std::env::var_os("ORIENTATION_CHECK_MATCH").map(|directory| {
        let mut fixtures = std::fs::read_dir(directory)
            .unwrap()
            .map(|entry| entry.unwrap().path())
            .filter(|path| path.extension().is_some_and(|extension| extension == "dot"))
            .map(|path| {
                let diagram = FeynmanDiagram::from_dot(
                    model.clone(),
                    &std::fs::read_to_string(&path).unwrap(),
                )
                .unwrap();
                (
                    path.file_stem().unwrap().to_string_lossy().into_owned(),
                    diagram.id(),
                )
            })
            .collect::<Vec<_>>();
        fixtures.sort();
        fixtures
    });
    let start = Instant::now();
    let result = process.generate_diagrams(model, &options).unwrap();
    let seconds = start.elapsed().as_secs_f64();
    if let Some(fixtures) = &fixtures {
        let ids = result
            .diagrams
            .iter()
            .map(|diagram| (diagram.id(), diagram.name().to_owned()))
            .collect::<BTreeMap<_, _>>();
        for (fixture, id) in fixtures {
            eprintln!(
                "{{\"fixture\":\"{fixture}\",\"retained_as\":{}}}",
                ids.get(id)
                    .map_or("null".to_owned(), |name| format!("\"{name}\""))
            );
        }
    }
    if let Some(path) = std::env::var_os("ORIENTATION_CHECK_GROUPS") {
        let mut lines = String::new();
        for (diagram, group) in result.diagrams.iter().zip(&result.groups) {
            for member in &group.members {
                lines.push_str(&format!(
                    "{}\t{}\t{}\n",
                    diagram.name(),
                    member.source_name,
                    member.ratio
                ));
            }
        }
        std::fs::write(path, lines).unwrap();
    }
    if let Some(path) = std::env::var_os("ORIENTATION_CHECK_NAMES") {
        let names = result
            .diagrams
            .iter()
            .map(|diagram| format!("{}\t{:?}\n", diagram.name(), diagram.id()))
            .collect::<String>();
        std::fs::write(path, names).unwrap();
    }
    if let Some(path) = std::env::var_os("ORIENTATION_CHECK_VERDICTS") {
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .unwrap();
        let lines = pool.install(|| {
            result
                .diagrams
                .par_iter()
                .map(|diagram| {
                    let verdict = diagram
                        .closed_loop_orientation_sum(&GenerationOptions::default())
                        .unwrap();
                    if let Some(directory) = std::env::var_os("ORIENTATION_CHECK_DOT_DIR")
                        && !matches!(verdict, OrientationSumVerdict::NonZero)
                        && !matches!(verdict, OrientationSumVerdict::NoClosedLoops)
                    {
                        let path = std::path::Path::new(&directory)
                            .join(format!("{}.dot", diagram.name()));
                        std::fs::write(path, diagram.to_dot().unwrap()).unwrap();
                    }
                    format!(
                        "{}\t{}\t{verdict:?}\n",
                        diagram.name(),
                        closed_fermion_loops(diagram)
                    )
                })
                .collect::<Vec<_>>()
        });
        std::fs::write(path, lines.concat()).unwrap();
    }
    let mut per_nf = BTreeMap::<usize, usize>::new();
    let mut members_per_nf = BTreeMap::<usize, usize>::new();
    for (diagram, group) in result.diagrams.iter().zip(&result.groups) {
        let nf = closed_fermion_loops(diagram);
        *per_nf.entry(nf).or_default() += 1;
        *members_per_nf.entry(nf).or_default() += group.members.len();
    }
    let map = |counts: &BTreeMap<usize, usize>| {
        counts
            .iter()
            .map(|(nf, count)| format!("\"{nf}\":{count}"))
            .collect::<Vec<_>>()
            .join(",")
    };
    println!(
        "{{\"loops\":{loops},\"grouping\":\"{}\",\"flag\":{flag},\"threads\":{threads},\"retained\":{},\"retained_per_nf\":{{{}}},\"members_per_nf\":{{{}}},\"zero_numerator_count\":{},\"orientation_zero_count\":{},\"orientation_test_skipped\":{},\"completed\":{},\"seconds\":{seconds:.2}}}",
        args.get(3).map(String::as_str).unwrap_or("none"),
        result.report.retained_count,
        map(&per_nf),
        map(&members_per_nf),
        result.report.zero_numerator_count,
        result.report.orientation_zero_count,
        result.report.orientation_test_skipped,
        result.report.completed,
    );
}
