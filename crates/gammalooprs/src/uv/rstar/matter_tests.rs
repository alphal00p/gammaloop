use std::{collections::BTreeMap, sync::Arc, time::Instant};

use eyre::{Result, eyre};
use feynkit_generator::{
    GenerationFilter, GenerationOptions, NumeratorGrouping, Process, VertexSelector,
};
use linnet::half_edge::{
    NodeIndex,
    involution::{EdgeIndex, HedgePair},
};
use symbolica::atom::Atom;
use vakint::Vakint;

use crate::{
    debug_tags,
    feyngen::feynkit::FeynmanDiagramGammaLoopExt,
    graph::Graph,
    initialisation::test_initialise,
    model::Model,
    utils::load_generic_model,
    uv::{UVgenerationSettings, UltravioletGraph, uv_graph::UVE},
};

use super::{
    RStarPlan,
    tests::{assert_same_counterterm, canonical_tensor, massive_oracle, settings},
};

/// Select a ghost line through 2L vertices with L nested gluon chords.
/// The marks are the innermost and outermost gluon, respectively.
fn nested_ghost_marks(graph: &Graph, model: &Model) -> Option<(EdgeIndex, EdgeIndex)> {
    let mut external = Vec::new();
    let mut ghost_neighbors = BTreeMap::<NodeIndex, Vec<NodeIndex>>::new();
    let mut gluons = BTreeMap::new();
    let mut ghost_edges = 0;
    let key = |a: NodeIndex, b: NodeIndex| (a.min(b), a.max(b));
    for (pair, edge, data) in graph.iter_edges() {
        let pdg = data.data.particle_pdg_code(model)?.abs();
        match pair {
            HedgePair::Unpaired { hedge, .. } if pdg == 9_000_005 => {
                external.push(graph.node_id(hedge));
            }
            HedgePair::Paired { source, sink } => {
                let a = graph.node_id(source);
                let b = graph.node_id(sink);
                match pdg {
                    9_000_005 => {
                        ghost_edges += 1;
                        ghost_neighbors.entry(a).or_default().push(b);
                        ghost_neighbors.entry(b).or_default().push(a);
                    }
                    21 if gluons.insert(key(a, b), edge).is_none() => {}
                    _ => return None,
                }
            }
            _ => return None,
        }
    }
    let loops = gluons.len();
    let vertices = 2 * loops;
    if external.len() != 2 || loops == 0 || ghost_edges + 1 != vertices {
        return None;
    }
    let mut line = vec![external[0]];
    while line.len() < vertices {
        let current = *line.last()?;
        let next = ghost_neighbors
            .get(&current)?
            .iter()
            .filter(|next| !line.contains(next))
            .copied()
            .collect::<Vec<_>>();
        if next.len() != 1 {
            return None;
        }
        line.push(next[0]);
    }
    if line[vertices - 1] != external[1]
        || (0..loops)
            .any(|index| !gluons.contains_key(&key(line[index], line[vertices - 1 - index])))
    {
        return None;
    }
    Some((
        *gluons.get(&key(line[loops - 1], line[loops]))?,
        *gluons.get(&key(line[0], line[vertices - 1]))?,
    ))
}

fn ghost_self_energy(
    model: &Model,
    loops: usize,
    accept: impl Fn(&Graph) -> bool,
) -> Result<Graph> {
    // V_35 is the SM ghost-antighost-gluon rule, including its antighost
    // derivative, adjoint structure constant and coupling. No numerator
    // fragments or species are replaced by scalar stand-ins.
    let options = GenerationOptions::default()
        .with_loop_count(loops, loops)?
        .max_vertices(2 * loops)
        .allow_self_loops(false)
        .projector(Atom::one())
        .numerator_grouping(NumeratorGrouping::None)
        .with_graph_filter(GenerationFilter::VertexAllow(vec![VertexSelector::Name(
            "V_35".into(),
        )]))
        .with_graph_filter(GenerationFilter::LoopOneParticleIrreducible)
        .threads(1);
    let generated =
        Process::new(["ghG"], ["ghG"]).generate_diagrams(Arc::new(model.clone()), &options)?;
    for diagram in generated.diagrams {
        let graph = diagram.to_gamma_loop_graph(None, true)?;
        if accept(&graph) {
            return Ok(graph);
        }
    }
    Err(eyre!(
        "no matching {loops}-loop SM ghost self-energy was generated"
    ))
}

fn automatic_ghost_plan(
    graph: &Graph,
    model: &Model,
    settings: &UVgenerationSettings,
) -> Result<RStarPlan> {
    let started = Instant::now();
    let plan = RStarPlan::with_automatic_edge(graph, model, settings)?;
    debug_tags!(#uv, #integrated, #rstar;
        stage = "ghost_automatic_hard_line_selection",
        loops = graph.n_loops(&plan.support),
        marked_edge = usize::from(plan.marked_edge),
        elapsed_ms = started.elapsed().as_millis() as u64,
        "Selected the massive ghost self-energy line"
    );
    Ok(plan)
}

fn timed_ghost_counterterm(
    plan: &RStarPlan,
    graph: &Graph,
    model: &Model,
    vakint: &Vakint,
    settings: &UVgenerationSettings,
) -> Result<Atom> {
    let HedgePair::Paired { source, sink } = graph[&plan.marked_edge].1 else {
        unreachable!("the marked line is an internal propagator");
    };
    let started = Instant::now();
    let result = plan.counterterm(graph, model, vakint, settings)?;
    let elapsed_ms = started.elapsed().as_millis() as u64;
    let cost = plan.estimated_cost(graph);
    debug_tags!(#uv, #integrated, #rstar;
        stage = "ghost_hard_line_counterterm",
        loops = graph.n_loops(&plan.support),
        marked_edge = usize::from(plan.marked_edge),
        particle = graph[plan.marked_edge].particle_pdg_code(model),
        source = usize::from(graph.node_id(source)),
        sink = usize::from(graph.node_id(sink)),
        elapsed_ms,
        max_reduction_loops = cost.max_reduction_loops,
        ordinary_terms = cost.ordinary_terms,
        ir_terms = cost.ir_terms,
        quotient_calls = cost.quotient_calls,
        file.counterterm = %canonical_tensor(&result)?,
        "Evaluated the ghost self-energy counterterm for this massive line"
    );
    Ok(result)
}

#[test]
fn one_loop_ghost_taylor_keeps_the_external_vertex_momentum() -> Result<()> {
    test_initialise()?;
    let model = load_generic_model("sm");
    let graph = ghost_self_energy(&model, 1, |_| true)?;
    let marked = graph
        .iter_edges()
        .find_map(|(pair, edge, data)| {
            (pair.is_paired() && data.data.particle_pdg_code(&model) == Some(21)).then_some(edge)
        })
        .unwrap();
    let settings = settings();
    let vakint = Vakint::new()?;
    let plan = RStarPlan::new(&graph, &model, marked, &settings)?;
    // One antighost derivative already supplies an external momentum. The
    // loop UV degree is one, while Taylor expansion of the complete numerator
    // needs total degree two to retain its nonzero p^2 counterterm.
    assert_eq!(graph.compute_dod(&plan.support), 1);
    assert_eq!(plan.taylor_order, 2);
    let actual = plan.counterterm(&graph, &model, &vakint, &settings)?;
    let expected = massive_oracle(&plan, &graph, &model, &vakint, &settings)?;
    assert!(!canonical_tensor(&expected)?.is_zero());
    assert_same_counterterm(&actual, &expected)
}

#[test]
fn three_loop_nested_ghost_hard_taylor_matches_massive_r() -> Result<()> {
    test_initialise()?;
    let model = load_generic_model("sm");
    let graph = ghost_self_energy(&model, 3, |graph| {
        nested_ghost_marks(graph, &model).is_some()
    })?;
    let (inside, _) = nested_ghost_marks(&graph, &model).unwrap();
    let settings = settings();
    let vakint = Vakint::new()?;
    let inside = RStarPlan::new(&graph, &model, inside, &settings)?;
    assert_eq!(graph.n_loops(&inside.support), 3);
    assert_eq!(graph.iter_edges_of(&inside.support).count(), 8);
    assert!(inside.taylor_order > 0);
    assert!(
        inside
            .ir
            .iter()
            .any(|hard| { graph.n_loops(&hard.hard) == 2 && hard.taylor_order > 0 })
    );
    assert!(inside.ir.iter().any(|hard| {
        hard.taylor_order > 0 && hard.forest.iter().any(|forest| forest.steps.len() >= 2)
    }));

    let automatic = automatic_ghost_plan(&graph, &model, &settings)?;
    let automatic_result = timed_ghost_counterterm(&automatic, &graph, &model, &vakint, &settings)?;
    let inside_result = if inside.marked_edge == automatic.marked_edge {
        automatic_result.clone()
    } else {
        timed_ghost_counterterm(&inside, &graph, &model, &vakint, &settings)?
    };
    assert_same_counterterm(&automatic_result, &inside_result)?;
    // Vary every internal line, including the five directed ghost propagators
    // and all three gluon chords. Different species and positions produce
    // different hard preparations while the local UV counterterm is fixed.
    for (pair, marked, _) in graph.iter_edges() {
        if !pair.is_paired() || marked == inside.marked_edge || marked == automatic.marked_edge {
            continue;
        }
        let plan = RStarPlan::new(&graph, &model, marked, &settings)?;
        let result = timed_ghost_counterterm(&plan, &graph, &model, &vakint, &settings)?;
        assert_same_counterterm(&inside_result, &result)?;
    }

    // Repeat the selected line after the sweep so its diagnostic timing can
    // also be compared after the first-call initialization costs are paid.
    let repeated = timed_ghost_counterterm(&automatic, &graph, &model, &vakint, &settings)?;
    assert_same_counterterm(&automatic_result, &repeated)?;

    // Only the independent oracle integrates a three-loop massive vacuum.
    // All R* paths retain actual ghost/gluon Feynman numerators and use
    // massive R only for lower-loop UV operators.
    let expected = massive_oracle(&inside, &graph, &model, &vakint, &settings)?;
    assert!(!canonical_tensor(&expected)?.is_zero());
    debug_tags!(#uv, #integrated, #rstar;
        stage = "three_loop_ghost_massive_oracle",
        log.counterterm = canonical_tensor(&expected)?,
        "Three-loop nested ghost massive R reference"
    );
    assert_same_counterterm(&inside_result, &expected)
}

#[test]
fn four_loop_nested_ghost_hard_taylor_is_hard_line_independent() -> Result<()> {
    test_initialise()?;
    let model = load_generic_model("sm");
    // Seven ghost propagators form A-B-C-D-E-F-G-H. The four gluons join
    // AH, BG, CF and DE, with all eight vertices generated from V_35.
    let graph = ghost_self_energy(&model, 4, |graph| {
        nested_ghost_marks(graph, &model).is_some()
    })?;
    let (inside, _) = nested_ghost_marks(&graph, &model).unwrap();
    let settings = settings();
    let vakint = Vakint::new()?;
    let inside = RStarPlan::new(&graph, &model, inside, &settings)?;
    assert_eq!(graph.n_loops(&inside.support), 4);
    assert_eq!(graph.iter_edges_of(&inside.support).count(), 11);
    assert_eq!(inside.taylor_order, 2);
    assert!(
        inside
            .ir
            .iter()
            .any(|hard| { graph.n_loops(&hard.hard) == 3 && hard.taylor_order > 0 })
    );
    assert!(inside.ir.iter().any(|hard| {
        hard.taylor_order > 0 && hard.forest.iter().any(|forest| forest.steps.len() >= 3)
    }));

    let automatic = automatic_ghost_plan(&graph, &model, &settings)?;
    let automatic_result = timed_ghost_counterterm(&automatic, &graph, &model, &vakint, &settings)?;
    let inside_result = if inside.marked_edge == automatic.marked_edge {
        automatic_result.clone()
    } else {
        timed_ghost_counterterm(&inside, &graph, &model, &vakint, &settings)?
    };
    assert!(!canonical_tensor(&inside_result)?.is_zero());
    assert_same_counterterm(&automatic_result, &inside_result)?;
    // Vary all eleven internal lines, including every directed ghost and gluon.
    // Every production UV sub-operator has at most three loops; the root uses
    // Vakint's massless propagator reduction and the single-mass closure.
    for (pair, marked, _) in graph.iter_edges() {
        if !pair.is_paired() || marked == inside.marked_edge || marked == automatic.marked_edge {
            continue;
        }
        let plan = RStarPlan::new(&graph, &model, marked, &settings)?;
        let result = timed_ghost_counterterm(&plan, &graph, &model, &vakint, &settings)?;
        assert_same_counterterm(&inside_result, &result)?;
    }
    let repeated = timed_ghost_counterterm(&automatic, &graph, &model, &vakint, &settings)?;
    assert_same_counterterm(&automatic_result, &repeated)?;
    Ok(())
}
