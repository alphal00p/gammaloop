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
    uv::{
        UVgenerationSettings,
        approx::streamed::{CheckTally, with_algebraic_check},
    },
};

use super::{RStarPlan, UvChildrenRoute};

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
    // Massive R children (the R* route, the default, needs a massless graph).
    let plan = RStarPlan::new(&graph, &model, EdgeIndex::from(0), &settings)?
        .with_uv_children(UvChildrenRoute::Massive);
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

/// A numerator stress case (built-in 1L-3L data): name, vacuum lines
/// `(from, to, power)` and numerator monomials `[(line, power)]`.
type StressCase = (
    &'static str,
    &'static [(usize, usize, u32)],
    &'static [&'static [(usize, i64)]],
);

/// n1: bubble crossed by q = Q_a + Q_b, whose loop dependence cancels.
const N1: StressCase = (
    "n1_cancel_bubble",
    &[(0, 1, 1), (0, 1, 1), (1, 2, 1), (1, 2, 1), (2, 0, 1)],
    &[&[(0, 1), (1, 1)]],
);

/// n3c: dotted self-loop with a cancelling bubble numerator.
const N3C: StressCase = (
    "n3c_selfloop_dot_cancel",
    &[(0, 0, 2), (0, 1, 1), (0, 1, 1), (1, 0, 1)],
    &[&[(1, 1), (2, 1)]],
);

/// n4: one-vertex-reducible subgraph (two bubbles at a cut vertex) with a
/// numerator straddling both blocks.
const N4: StressCase = (
    "n4_1VR_straddle",
    &[(0, 1, 1), (0, 1, 1), (1, 2, 1), (1, 2, 1), (2, 0, 1)],
    &[&[(0, 1), (2, 1)]],
);

/// n5: disconnected spinney (two vertex-disjoint bubbles) with a numerator
/// straddling both components.
const N5: StressCase = (
    "n5_disconnected_straddle",
    &[
        (0, 1, 1),
        (0, 1, 1),
        (2, 3, 1),
        (2, 3, 1),
        (1, 2, 1),
        (3, 0, 1),
    ],
    &[&[(0, 1), (2, 1)]],
);

/// The numerator-monomial path on 3-loop stress cases (built-in 1L-3L data):
/// the identity `K G = -Delta(G) - K sum_S Delta(S) G/S` must equal Vakint's
/// direct evaluation with `Delta(G)` from the ordinary forest AND from R*.
///
/// The identity is about the massive R operation, not the momentum
/// linearisation: with `GAMMALOOP_LINEARIZE_MOMENTA` on (the default) the
/// factor-by-factor momentum routing of n5 is a different spelling than the
/// whole map's, which the test-only `report_check` rejects. The
/// linearisation is off here, as before it became the default; the
/// `_linearized` tests ([`assert_identity_with_numerators_linearized`]) run
/// the same identities with it on.
fn assert_identity_with_numerators(case: StressCase) -> Result<()> {
    super::linearize::with_linearize(false, || identity_with_numerators(case))
}

/// [`assert_identity_with_numerators`] with the momentum linearisation on
/// (the production default), through the same two pole assertions. The
/// streamed maps are compared with the whole maps algebraically (a
/// structural difference must have an expanded difference of zero), and the
/// streamed result goes on: the forest and R* poles are computed from the
/// linearised factor-by-factor routing. Returns the stream checks per stage.
fn assert_identity_with_numerators_linearized(case: StressCase) -> Result<CheckTally> {
    let (result, tally) = super::linearize::with_linearize(true, || {
        assert!(super::linearize::enabled());
        with_algebraic_check(|| identity_with_numerators(case))
    });
    result?;
    eprintln!(
        "{}: stream checks with the linearisation on, per stage (identical, algebraic only): {tally:?}",
        case.0
    );
    let (identical, algebraic) = tally.get("momentum_routing").copied().unwrap_or_default();
    assert!(
        identical + algebraic > 0,
        "{}: the momentum routing was never checked",
        case.0
    );
    Ok(tally)
}

fn identity_with_numerators((name, edges, numerators): StressCase) -> Result<()> {
    test_initialise()?;
    let model = load_generic_model("scalars");
    let settings = super::tests::settings();
    let vakint = Vakint::new()?;
    let topology = super::VacuumTopology {
        name: name.into(),
        edges: edges.to_vec(),
        numerators: numerators
            .iter()
            .map(|monomial| monomial.to_vec())
            .collect(),
    };
    let (mut poles, mut direct, mut delta_forest) = (Atom::Zero, Atom::Zero, Atom::Zero);
    for monomial in topology.numerator_monomials(&mass()) {
        let graph = topology.monomial_graph(&model, &mass(), &monomial.attachments)?;
        // Massive R children (the R* route, the default, needs a massless graph).
        let plan = RStarPlan::new(&graph, &model, EdgeIndex::from(0), &settings)?
            .with_uv_children(UvChildrenRoute::Massive);
        let sub = plan.massive_vacuum_subdivergence_poles(
            &graph,
            &model,
            &mass(),
            &vakint,
            &settings,
            0,
        )?;
        poles += &monomial.coefficient * &sub.poles;
        direct += &monomial.coefficient
            * super::massive_vacuum_direct(&graph, &model, &vakint, &settings, 0)?;
        delta_forest += &monomial.coefficient
            * plan.massive_vacuum_counterterm(&graph, &model, &vakint, &settings)?;
    }
    let direct = direct.series(GS.dim_epsilon, 0, -1)?.to_atom().expand();
    let poles = poles.series(GS.dim_epsilon, 0, -1)?.to_atom().expand();
    let (delta_rstar, _) =
        topology.overall_counterterm_rstar(&model, &mass(), &vakint, &settings)?;
    for (route, delta) in [("forest", delta_forest), ("R*", delta_rstar)] {
        assert!(
            !delta.contains_symbol(GS.m_uv_vacuum),
            "{name}: Delta ({route}) depends on m_uv"
        );
        let difference = (&poles - &delta - &direct).expand().together();
        assert!(
            difference.is_zero(),
            "{name}: identity with Delta ({route}) fails by {difference}"
        );
    }
    Ok(())
}

/// n1 ([`N1`]), linearisation off.
#[test]
fn numerator_cancelling_inside_a_bubble() -> Result<()> {
    assert_identity_with_numerators(N1)
}

/// n3c ([`N3C`]), linearisation off.
#[test]
fn numerator_on_a_bubble_next_to_a_dotted_self_loop() -> Result<()> {
    assert_identity_with_numerators(N3C)
}

/// n4 ([`N4`]), linearisation off.
#[test]
fn numerator_straddling_a_one_vertex_reducible_subgraph() -> Result<()> {
    assert_identity_with_numerators(N4)
}

/// n5 ([`N5`]), linearisation off.
#[test]
fn numerator_straddling_a_disconnected_spinney() -> Result<()> {
    assert_identity_with_numerators(N5)
}

/// n1 ([`N1`]), linearisation on.
#[test]
fn numerator_cancelling_inside_a_bubble_linearized() -> Result<()> {
    assert_identity_with_numerators_linearized(N1).map(drop)
}

/// n3c ([`N3C`]), linearisation on.
#[test]
fn numerator_on_a_bubble_next_to_a_dotted_self_loop_linearized() -> Result<()> {
    assert_identity_with_numerators_linearized(N3C).map(drop)
}

/// n4 ([`N4`]), linearisation on.
#[test]
fn numerator_straddling_a_one_vertex_reducible_subgraph_linearized() -> Result<()> {
    assert_identity_with_numerators_linearized(N4).map(drop)
}

/// n5 ([`N5`]), linearisation on. Its factor-by-factor momentum routing is
/// spelled differently from the whole map's (the stream check that stopped
/// the default-on run before both pole assertions); here that linearised
/// spelling is what the forest and R* poles are computed from.
#[test]
fn numerator_straddling_a_disconnected_spinney_linearized() -> Result<()> {
    let tally = assert_identity_with_numerators_linearized(N5)?;
    let (_, algebraic) = tally["momentum_routing"];
    assert!(
        algebraic > 0,
        "n5: the linearised momentum routing is now spelled as the whole map's; \
         this test no longer exercises a different spelling than the off-mode test"
    );
    Ok(())
}

/// `power_distributions` with a line of power 0: no underflow, the absent
/// propagator takes no extra power, and the weights are the binomials of the
/// other lines (`1/(k^2-m^2)^2 = sum_a (a+1) m^(2a) / (k^2)^(2+a)`).
#[test]
fn power_distributions_with_a_power_zero_line() {
    let topology = super::VacuumTopology {
        name: "power_zero_line".into(),
        edges: vec![(0, 1, 2), (0, 1, 0), (1, 0, 1)],
        numerators: Vec::new(),
    };
    let mut distributions = topology.power_distributions(1);
    distributions.sort();
    assert_eq!(distributions, vec![(vec![0, 0, 1], 1), (vec![1, 0, 0], 2)]);
    assert_eq!(topology.power_distributions(0), vec![(vec![0, 0, 0], 1)]);
}
