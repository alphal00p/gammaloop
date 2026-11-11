//! Regressions for massive R child counterterms (fixed upstream in the UV child
//! boundary routing): the correct pole-part counterterm of the massive 3-loop
//! banana (four lines of mass `m`, momentum `p` through it) in the conventions
//! of `tests::settings()` is
//!
//! `-2 m^4/eps^3 + 10/3 m^4/eps^2 - 2/3 m^2 p^2/eps^2
//!  + 1/2 m^2 p^2/eps - 1/2 m^4/eps - 1/36 p^4/eps`,
//!
//! whether the banana is a two-point graph written with ported or unported
//! legs, or the in-graph child "lines 0-3" of the 5-line banana vacuum graph
//! whose boundary is both ends of line 4 (formerly `+1/12 m^2 p^2/eps^2` and
//! m_uv terms). Only Vakint's built-in 1L-3L data are needed.

use eyre::Result;
use linnet::half_edge::{
    involution::EdgeIndex,
    subgraph::{ModifySubSet, SuBitGraph},
};
use spenso::{
    network::tags::SPENSO_TAG,
    structure::representation::{Minkowski, RepName},
};
use symbolica::{
    atom::{Atom, AtomCore, Symbol},
    function, symbol,
};
use vakint::Vakint;

use crate::{
    graph::{Graph, parse::IntoFinalizedRuntimeGraph},
    initialisation::test_initialise,
    model::Model,
    utils::{GS, load_generic_model},
    uv::UVgenerationSettings,
};

use super::RStarPlan;

fn mass() -> Atom {
    Atom::var(symbol!("massive_r_regression_mass"; Scalar))
}

fn p2() -> Symbol {
    symbol!("massive_r_regression_p2")
}

fn expected_banana_counterterm() -> Atom {
    let m2 = mass().pow(Atom::num(2));
    let p2 = Atom::var(p2());
    let e = Atom::var(GS.dim_epsilon);
    (Atom::num(-2) * m2.pow(Atom::num(2)) / e.pow(Atom::num(3))
        + Atom::num(10) / Atom::num(3) * m2.pow(Atom::num(2)) / e.pow(Atom::num(2))
        - Atom::num(2) / Atom::num(3) * &m2 * &p2 / e.pow(Atom::num(2))
        + Atom::num(1) / Atom::num(2) * &m2 * &p2 / &e
        - Atom::num(1) / Atom::num(2) * m2.pow(Atom::num(2)) / &e
        - Atom::num(1) / Atom::num(36) * p2.pow(Atom::num(2)) / &e)
        .expand()
}

/// Replace `Q(edge).Q(edge)` by `p2` for each listed edge.
fn in_p2(counterterm: &Atom, edges: &[usize]) -> Atom {
    edges.iter().fold(counterterm.clone(), |value, edge| {
        let momentum = function!(
            GS.emr_mom,
            *edge,
            Minkowski {}.new_rep(GS.dim).to_symbolic([])
        );
        value
            .replace(function!(SPENSO_TAG.dot, &momentum, &momentum).to_pattern())
            .with(Atom::var(p2()).to_pattern())
            .expand()
    })
}

fn massive(mut graph: Graph) -> Graph {
    let edges: Vec<_> = graph.iter_edges().map(|(_, edge, _)| edge).collect();
    for edge in edges {
        let data = &mut graph.underlying[edge];
        data.particle = data.particle.clone().override_mass(Some(mass()));
    }
    graph
}

/// The full forest counterterm of a graph's root (the ordinary path).
fn root_counterterm(
    graph: &Graph,
    model: &Model,
    vakint: &Vakint,
    settings: &UVgenerationSettings,
) -> Result<Atom> {
    let plan = RStarPlan::new(graph, model, EdgeIndex::from(0), settings)?;
    let settings = UVgenerationSettings {
        generate_integrated: true,
        add_marker: false,
        ..settings.clone()
    };
    let mut forest = plan.uv_forest(graph, &settings);
    forest.integrate(graph, model, vakint, &settings)?;
    Ok(forest
        .counterterm_for_support(&plan.support)?
        .replace(GS.m_uv_expansion)
        .with(GS.m_uv_vacuum)
        .expand())
}

fn two_point_banana(model: &Model, ports: bool) -> Result<Graph> {
    let (a, b) = if ports { (":0", ":1") } else { ("", "") };
    let dot = format!(
        "digraph banana {{
        graph [projector=1]; node [num=1];
        edge [particle=\"scalar_1\" num=1]; e [style=invis];
        e -> A{a} [id=4 sink=\"{{ufo_order:0}}\"];
        B{b} -> e [id=5 source=\"{{ufo_order:0}}\"];
        A -> B [id=0 source=\"{{ufo_order:1}}\" sink=\"{{ufo_order:1}}\"];
        A -> B [id=1 lmb_id=0 source=\"{{ufo_order:2}}\" sink=\"{{ufo_order:2}}\"];
        A -> B [id=2 lmb_id=1 source=\"{{ufo_order:3}}\" sink=\"{{ufo_order:3}}\"];
        A -> B [id=3 lmb_id=2 source=\"{{ufo_order:4}}\" sink=\"{{ufo_order:4}}\"];
    }}"
    );
    Ok(massive(dot.into_finalized_runtime_graph(model)?))
}

fn assert_banana(actual: &Atom, what: &str) {
    let expected = expected_banana_counterterm();
    assert!(
        !actual.contains_symbol(GS.m_uv_vacuum),
        "{what}: the counterterm depends on m_uv: {actual}"
    );
    assert!(
        (actual - &expected).expand().is_zero(),
        "{what}: {actual}; expected {expected}"
    );
}

/// The banana two-point graph with ported legs.
#[test]
fn massive_banana_two_point_counterterm_with_ported_legs() -> Result<()> {
    test_initialise()?;
    let model = load_generic_model("scalars");
    let settings = super::tests::settings();
    let vakint = Vakint::new()?;
    let graph = two_point_banana(&model, true)?;
    let actual = in_p2(
        &root_counterterm(&graph, &model, &vakint, &settings)?,
        &[4, 5],
    );
    assert_banana(&actual, "ported two-point banana");
    Ok(())
}

/// The identical graph with unported legs (`e -> A`, `B -> e`); formerly
/// `13/3 m^4/eps^2` instead of `10/3` and `m_uv` terms.
#[test]
fn massive_banana_two_point_counterterm_with_unported_legs() -> Result<()> {
    test_initialise()?;
    let model = load_generic_model("scalars");
    let settings = super::tests::settings();
    let vakint = Vakint::new()?;
    let graph = two_point_banana(&model, false)?;
    let actual = in_p2(
        &root_counterterm(&graph, &model, &vakint, &settings)?,
        &[4, 5],
    );
    assert_banana(&actual, "unported two-point banana");
    Ok(())
}

/// The 3-loop banana child (lines 0-3) of the 5-line banana vacuum graph,
/// whose boundary is both ends of line 4, through the R* children path:
/// formerly `+1/12 m^2 p^2/eps^2` instead of `-2/3 m^2 p^2/eps^2`, and `m_uv` terms.
#[test]
fn massive_r_child_with_both_ends_of_one_boundary_line() -> Result<()> {
    test_initialise()?;
    let model = load_generic_model("scalars");
    let settings = super::tests::settings();
    let vakint = Vakint::new()?;
    let graph = massive(
        "digraph five_line_banana {
        graph [projector=1]; node [num=1];
        edge [particle=\"scalar_1\" num=1];
        B -> A [id=0 source=\"{ufo_order:0}\" sink=\"{ufo_order:0}\"];
        B -> A [id=1 lmb_id=0 source=\"{ufo_order:1}\" sink=\"{ufo_order:1}\"];
        A -> B [id=2 lmb_id=1 source=\"{ufo_order:2}\" sink=\"{ufo_order:2}\"];
        A -> B [id=3 lmb_id=2 source=\"{ufo_order:3}\" sink=\"{ufo_order:3}\"];
        B -> A [id=4 lmb_id=3 source=\"{ufo_order:4}\" sink=\"{ufo_order:4}\"];
    }"
        .into_finalized_runtime_graph(&model)?,
    );
    let plan = RStarPlan::new(&graph, &model, EdgeIndex::from(0), &settings)?;
    let mut child: SuBitGraph = graph.empty_subgraph();
    for edge in 0..4 {
        child.add(graph[&EdgeIndex::from(edge)].1);
    }
    let children = plan.evaluate_uv_children(&graph, &model, &mass(), &vakint, &settings)?;
    let actual = in_p2(
        &children[&child]
            .replace(GS.m_uv_expansion)
            .with(GS.m_uv_vacuum)
            .expand(),
        &[4],
    );
    assert_banana(&actual, "in-graph banana child");
    Ok(())
}
