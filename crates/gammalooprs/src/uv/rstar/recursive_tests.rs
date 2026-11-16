//! Proper UV children computed recursively by R* (`GAMMALOOP_UV_CHILDREN=rstar`)
//! must equal the massive R operation's `Delta(S)` exactly.
//!
//! The inputs are production background-field QCD diagrams (the prepared
//! two- to four-loop inputs of the beta-function project and its five-loop
//! diagram 41295), loaded with the project's model. Both results are routed
//! in the parent graph's loop-momentum basis, which imposes momentum
//! conservation on the child's boundary momenta, before tensor canonization.
//!
//! Child comparisons use colour-factored diagrams, as production does (beta's
//! `color::factor_color` keeps one colour row of every local rule and moves
//! the closed colour network into a coefficient outside R*): the open-index
//! colour simplifier has no unique basis for children with several open
//! adjoint indices, so full-colour children are compared only where it has.
//!
//! The tests run with the production stream piece limit: they compare the two
//! routes, not the streamed maps. (With the test limit of one byte, the
//! streamed outer simplification of a child root with open Dirac chains is a
//! different spelling of the whole map's result, which `report_check` rejects.)

use std::{
    collections::BTreeSet,
    sync::Arc,
    time::{Duration, Instant},
};

use eyre::{Result, WrapErr, ensure, eyre};
use feynkit_graph::FeynmanDiagram;
use idenso::{
    CookMode, CookSettings,
    representations::{ColorAdjoint, ColorFundamental, ColorSextet},
    tensor::SymbolicTensor,
};
use spenso::shadowing::TensorCollectFilter;
use linnet::half_edge::{
    involution::EdgeIndex,
    subgraph::{Inclusion, ModifySubSet, SuBitGraph, SubSetLike},
};
use symbolica::atom::{Atom, AtomCore};
use vakint::Vakint;

use crate::{
    feyngen::feynkit::FeynmanDiagramGammaLoopExt,
    graph::{Graph, LMBext},
    initialisation::test_initialise,
    model::Model,
    utils::{GS, W_},
    uv::{
        ApproximationType, RenormalizationPrescriptionSettings, UVgenerationSettings,
        approx::streamed::{with_piece_limit, with_production_piece_limit},
        settings::VakintSettings as GammaVakintSettings, uv_graph::UVE,
    },
};

use super::{ChildMass, CountertermCache, RStarPlan, UvChildrenRoute, tests::canonical_tensor};

const MODEL: &str = include_str!("testdata/background_qcd.json");
const RSTAR: UvChildrenRoute = UvChildrenRoute::Rstar { min_loops: 1 };

fn load(dot: &str) -> Result<FeynmanDiagram> {
    let model = Model::from_json(MODEL).map_err(|error| eyre!("{error}"))?;
    FeynmanDiagram::from_dot(model, dot).map_err(|error| eyre!("{error}"))
}

/// One colour channel of the diagram with its colour factored out: every
/// vertex and propagator rule keeps the residual of its first colour row, as
/// in a branch of beta's `color::factor_color` (which also weights rows of
/// multi-row rules; the R* operations are linear in each local rule).
fn colour_factored(diagram: &FeynmanDiagram) -> Result<FeynmanDiagram> {
    let filter = || {
        TensorCollectFilter::Reps([
            ColorAdjoint {}.into(),
            ColorFundamental {}.into(),
            ColorSextet {}.into(),
        ])
    };
    let mut vertices = Vec::new();
    for (id, vertex) in diagram.vertices() {
        let rows = vertex.numerator.coefficient_list(filter())?;
        if let Some((_, residual)) = rows.into_iter().next() {
            vertices.push((id, residual));
        }
    }
    let mut edges = Vec::new();
    for (id, _, edge) in diagram.edges() {
        let rows = edge.numerator.coefficient_list(filter())?;
        if let Some((_, residual)) = rows.into_iter().next() {
            edges.push((id, residual));
        }
    }
    diagram
        .clone()
        .with_numerator_fragments(vertices, edges)
        .map_err(|error| eyre!("{error}"))
}

/// The beta engine's UV settings (MS-bar, RustRed before Forcer).
fn settings() -> UVgenerationSettings {
    UVgenerationSettings {
        renormalization_prescription: RenormalizationPrescriptionSettings {
            log_divergent: ApproximationType::PolePart,
            massive_power_divergent: ApproximationType::PolePart,
            massless_power_divergent: ApproximationType::PolePart,
            ..Default::default()
        },
        vakint: GammaVakintSettings {
            evaluation_methods: vec!["rustred".into(), "forcer".into()],
            normalization: "MSbar".into(),
            additional_normalization: "1".into(),
            ..Default::default()
        },
        ..Default::default()
    }
}

/// The settings of the massive route's child requests (`uv_family_requests`).
fn massive_settings(settings: &UVgenerationSettings) -> UVgenerationSettings {
    UVgenerationSettings {
        generate_integrated: true,
        add_marker: false,
        subtract_uv: true,
        renormalization_prescription: RenormalizationPrescriptionSettings {
            log_divergent: ApproximationType::PolePart,
            massive_power_divergent: ApproximationType::PolePart,
            massless_power_divergent: ApproximationType::PolePart,
            overrides: Vec::new(),
        },
        ..settings.clone()
    }
}

/// Route every edge momentum in the parent's basis: boundary momenta of a
/// child are then independent variables, whichever subset each route used.
fn parent_frame(graph: &Graph, expression: &Atom) -> Atom {
    let all: SuBitGraph = graph.underlying.full_filter();
    let rules =
        graph.normal_emr_replacement(&all, &graph.loop_momentum_basis, &[W_.x___], |_| true);
    expression.replace_multiple(&rules).expand()
}

fn assert_same_delta(graph: &Graph, label: &str, actual: &Atom, expected: &Atom) -> Result<()> {
    let cooking = CookSettings::indices()
        .with_mode(CookMode::ReversibleEncoding)
        .with_representation_payloads(true, true);
    let canonical = |value: &Atom| -> Result<Atom> {
        let simplified = canonical_tensor(&parent_frame(graph, value))?;
        let tensor =
            SymbolicTensor::infer_indexed(cooking.try_cook(simplified.as_view())?)?.canonize()?;
        Ok(cooking.uncook(tensor.expression().as_view()).expand())
    };
    let actual = canonical(actual)?;
    let expected = canonical(expected)?;
    if (&actual - &expected).together().is_zero() {
        return Ok(());
    }
    // Different canonical spellings of one tensor: reduce the difference.
    let difference = canonical_tensor(&(&actual - &expected))?;
    ensure!(
        difference.together().is_zero(),
        "{label}: recursive R*: {actual}\nmassive R: {expected}\ndifference: {difference}"
    );
    Ok(())
}

#[derive(Default, Debug)]
struct Coverage {
    children: usize,
    nonzero: usize,
    quark: usize,
    ghost: usize,
    degrees: BTreeSet<i32>,
    loops: BTreeSet<usize>,
    nested: usize,
    marked: usize,
    rstar_time: Duration,
    massive_time: Duration,
}

fn species(graph: &Graph, model: &Model, support: &SuBitGraph) -> (bool, bool) {
    let mut quark = false;
    let mut ghost = false;
    for (_, _, data) in graph.iter_edges_of(support) {
        match data.data.particle_pdg_code(model).map(isize::abs) {
            Some(1) => quark = true,
            Some(9000005) => ghost = true,
            _ => {}
        }
    }
    (quark, ghost)
}

/// Compare every proper UV child of the diagram's automatic R* plan, computed
/// on the massless graph by recursive R* and by the massive R operation.
fn compare_children(
    name: &str,
    dot: &str,
    vakint: &Vakint,
    max_loops: usize,
    full_colour: bool,
    coverage: &mut Coverage,
) -> Result<()> {
    let diagram = if full_colour {
        load(dot)?
    } else {
        colour_factored(&load(dot)?)?
    };
    let graph = diagram.to_gamma_loop_graph(None, true)?;
    let model = diagram.model();
    let settings = settings();
    let massive = massive_settings(&settings);
    let plan = RStarPlan::with_automatic_edge(&graph, model, &settings)?;
    let massive_cache = CountertermCache::default();
    let rstar_cache = Arc::new(CountertermCache::default());
    let marked = graph[&plan.marked_edge].1;
    for call in &plan.counterterms {
        let support = &call.argument.support;
        if call.argument.loops > max_loops {
            continue;
        }
        let started = Instant::now();
        let expected = massive_cache
            .counterterm(&graph, support, model, vakint, &massive)?
            .replace(GS.m_uv_expansion)
            .with(GS.m_uv_vacuum);
        coverage.massive_time += started.elapsed();
        let started = Instant::now();
        let actual = RStarPlan::rstar_child_counterterm(
            &rstar_cache,
            &graph,
            support,
            model,
            vakint,
            &settings,
            RSTAR,
        )
        .wrap_err_with(|| format!("{name}: R* child {}", support.string_label()))?;
        coverage.rstar_time += started.elapsed();
        ensure!(
            !actual.contains_symbol(GS.m_uv_vacuum)
                && !actual.contains_symbol(GS.m_uv_expansion)
                && !actual.contains_symbol(*super::REFERENCE_MASS_SQUARED),
            "{name}: an R* child retains an auxiliary or reference mass: {actual}"
        );
        let label = format!(
            "{name}{} child {} ({} loops, degree {})",
            if full_colour { " full colour" } else { "" },
            support.string_label(),
            call.argument.loops,
            call.taylor_order
        );
        assert_same_delta(&graph, &label, &actual, &expected)?;
        let (quark, ghost) = species(&graph, model, support);
        coverage.children += 1;
        coverage.nonzero += usize::from(!expected.expand().is_zero());
        coverage.quark += usize::from(quark);
        coverage.ghost += usize::from(ghost);
        coverage.degrees.insert(call.taylor_order);
        coverage.loops.insert(call.argument.loops);
        coverage.nested += usize::from(
            plan.counterterms
                .iter()
                .any(|other| other != call && support.includes(&other.argument.support)),
        );
        coverage.marked += usize::from(support.includes(&marked));
        eprintln!("RECURSIVE_CHILD_EQUAL {label} quark={quark} ghost={ghost}");
    }
    Ok(())
}

#[test]
fn background_two_and_three_loop_children_match_massive_r() -> Result<()> {
    with_production_piece_limit(background_two_and_three_loop_children_match_massive_r_body)
}

fn background_two_and_three_loop_children_match_massive_r_body() -> Result<()> {
    test_initialise()?;
    let vakint = Vakint::new()?;
    let mut coverage = Coverage::default();
    for (name, dot) in [
        ("L2/0", include_str!("testdata/l2-d0.dot")),
        ("L2/6", include_str!("testdata/l2-d6.dot")),
        ("L2/7", include_str!("testdata/l2-d7.dot")),
        ("L2/8", include_str!("testdata/l2-d8.dot")),
        ("L2/12", include_str!("testdata/l2-d12.dot")),
        ("L3/166", include_str!("testdata/l3-d166.dot")),
        ("L3/185", include_str!("testdata/l3-d185.dot")),
        ("L3/241", include_str!("testdata/l3-d241.dot")),
    ] {
        compare_children(name, dot, &vakint, usize::MAX, false, &mut coverage)?;
    }
    eprintln!("RECURSIVE_COVERAGE two/three loops: {coverage:?}");
    // Full colour, where every child has a unique open colour invariant.
    let mut full = Coverage::default();
    for (name, dot) in [
        ("L2/6", include_str!("testdata/l2-d6.dot")),
        ("L2/7", include_str!("testdata/l2-d7.dot")),
        ("L2/8", include_str!("testdata/l2-d8.dot")),
        ("L2/12", include_str!("testdata/l2-d12.dot")),
    ] {
        compare_children(name, dot, &vakint, usize::MAX, true, &mut full)?;
    }
    eprintln!("RECURSIVE_COVERAGE full colour: {full:?}");
    ensure!(coverage.quark > 0 && coverage.ghost > 0, "{coverage:?}");
    ensure!(coverage.loops.contains(&2) && coverage.nested > 0, "{coverage:?}");
    ensure!(coverage.degrees.is_superset(&BTreeSet::from([0, 1, 2])), "{coverage:?}");
    ensure!(coverage.nonzero > 0, "{coverage:?}");
    Ok(())
}

/// Four-loop graphs whose three-loop children are quark self-energies (the
/// one-loop-lower analogue of the four-loop insertion of diagram 41295), with
/// a gluon bubble or a ghost loop inside, and all of their nested children.
#[test]
fn background_four_loop_children_match_massive_r() -> Result<()> {
    with_production_piece_limit(background_four_loop_children_match_massive_r_body)
}

fn background_four_loop_children_match_massive_r_body() -> Result<()> {
    test_initialise()?;
    let vakint = Vakint::new()?;
    let mut coverage = Coverage::default();
    for (name, dot) in [
        ("L4/2774", include_str!("testdata/l4-d2774.dot")),
        ("L4/4613", include_str!("testdata/l4-d4613.dot")),
    ] {
        compare_children(name, dot, &vakint, 3, false, &mut coverage)?;
    }
    eprintln!("RECURSIVE_COVERAGE four loops: {coverage:?}");
    ensure!(coverage.loops.contains(&3), "{coverage:?}");
    ensure!(coverage.quark > 0 && coverage.ghost > 0, "{coverage:?}");
    Ok(())
}

/// The plan-level switch: children of either route, the root counterterm
/// built from them, and the extraction path's massless children.
#[test]
fn rstar_route_reproduces_plan_children_and_root_counterterms() -> Result<()> {
    with_production_piece_limit(rstar_route_reproduces_plan_children_and_root_counterterms_body)
}

fn rstar_route_reproduces_plan_children_and_root_counterterms_body() -> Result<()> {
    test_initialise()?;
    let vakint = Vakint::new()?;
    let settings = settings();
    let mass = Atom::var(*super::REFERENCE_MASS_SQUARED).sqrt();
    let mut fallbacks = 0;
    for (name, dot) in [
        ("L2/8", include_str!("testdata/l2-d8.dot")),
        ("L3/185", include_str!("testdata/l3-d185.dot")),
        ("L3/241", include_str!("testdata/l3-d241.dot")),
    ] {
        let diagram = load(dot)?;
        let graph = diagram.to_gamma_loop_graph(None, true)?;
        let model = diagram.model();
        let massive_plan = RStarPlan::with_automatic_edge(&graph, model, &settings)?
            .with_uv_children(UvChildrenRoute::Massive);
        let rstar_plan =
            RStarPlan::with_automatic_edge(&graph, model, &settings)?.with_uv_children(RSTAR);
        assert_eq!(massive_plan.marked_edge, rstar_plan.marked_edge);
        let marked = graph[&rstar_plan.marked_edge].1;
        let massive = massive_plan.evaluate_uv_children(&graph, model, &mass, &vakint, &settings)?;
        let rstar = rstar_plan.evaluate_uv_children(&graph, model, &mass, &vakint, &settings)?;
        let massless = rstar_plan.evaluate_uv_children_for(
            &graph,
            model,
            &mass,
            &vakint,
            &settings,
            ChildMass::Zero,
        )?;
        ensure!(massive.len() == rstar.len() && massive.len() == massless.len());
        for call in &rstar_plan.counterterms {
            let support = &call.argument.support;
            let label = format!("{name} plan child {}", support.string_label());
            if support.includes(&marked) && call.taylor_order >= 2 {
                fallbacks += 1;
            }
            assert_same_delta(&graph, &label, &rstar[support], &massive[support])?;
            let at_zero = massive[support]
                .replace(*super::REFERENCE_MASS_SQUARED)
                .with(Atom::Zero);
            assert_same_delta(&graph, &label, &massless[support], &at_zero)?;
        }
        let expected = massive_plan.counterterm(&graph, model, &vakint, &settings)?;
        let actual = rstar_plan.counterterm(&graph, model, &vakint, &settings)?;
        ensure!(!expected.expand().is_zero(), "{name}: zero root counterterm");
        assert_same_delta(&graph, &format!("{name} root"), &actual, &expected)?;
        // Extraction mode: bare poles plus the original forest, whose children
        // are consumed only at zero reference mass.
        let bare = massive_plan.pole_part(&graph, model, &vakint, &settings)?;
        let recovered =
            rstar_plan.counterterm_from_pole_part(&graph, model, &vakint, &settings, &bare)?;
        assert_same_delta(&graph, &format!("{name} from poles"), &recovered, &expected)?;
        eprintln!("RECURSIVE_ROOT_EQUAL {name}");
    }
    eprintln!("RECURSIVE_PLAN mass-dependent children kept massive: {fallbacks}");
    Ok(())
}

/// A marked line inside a quadratically divergent child: that child depends
/// on the reference mass, keeps the massive route for the full-graph root and
/// is computed by R* for the massless consumers. Both roots stay exact.
#[test]
fn mass_dependent_children_keep_the_massive_route_only_when_needed() -> Result<()> {
    with_production_piece_limit(mass_dependent_children_keep_the_massive_route_only_when_needed_body)
}

fn mass_dependent_children_keep_the_massive_route_only_when_needed_body() -> Result<()> {
    test_initialise()?;
    let vakint = Vakint::new()?;
    let settings = settings();
    let mass = Atom::var(*super::REFERENCE_MASS_SQUARED).sqrt();
    let mut fallbacks = 0;
    for (name, dot) in [
        ("L2/6", include_str!("testdata/l2-d6.dot")),
        ("L2/7", include_str!("testdata/l2-d7.dot")),
        ("L3/166", include_str!("testdata/l3-d166.dot")),
    ] {
        let diagram = colour_factored(&load(dot)?)?;
        let graph = diagram.to_gamma_loop_graph(None, true)?;
        let model = diagram.model();
        let automatic = RStarPlan::with_automatic_edge(&graph, model, &settings)?;
        let Some(quadratic) = automatic
            .counterterms
            .iter()
            .filter(|call| call.taylor_order >= 2)
            .max_by_key(|call| call.argument.loops)
        else {
            continue;
        };
        let marked = graph
            .iter_edges_of(&quadratic.argument.support)
            .find_map(|(pair, edge, _)| pair.is_paired().then_some(edge))
            .unwrap();
        let massive_plan = RStarPlan::new(&graph, model, marked, &settings)?
            .with_uv_children(UvChildrenRoute::Massive);
        let rstar_plan = RStarPlan::new(&graph, model, marked, &settings)?.with_uv_children(RSTAR);
        let massive = massive_plan.evaluate_uv_children(&graph, model, &mass, &vakint, &settings)?;
        let rstar = rstar_plan.evaluate_uv_children(&graph, model, &mass, &vakint, &settings)?;
        let massless = rstar_plan.evaluate_uv_children_for(
            &graph,
            model,
            &mass,
            &vakint,
            &settings,
            ChildMass::Zero,
        )?;
        let dependent = &massive[&quadratic.argument.support];
        ensure!(
            dependent.contains_symbol(*super::REFERENCE_MASS_SQUARED),
            "{name}: the marked quadratic child has no reference-mass term: {dependent}"
        );
        fallbacks += 1;
        for call in &rstar_plan.counterterms {
            let support = &call.argument.support;
            let label = format!("{name} marked-inside child {}", support.string_label());
            assert_same_delta(&graph, &label, &rstar[support], &massive[support])?;
            let at_zero = massive[support]
                .replace(*super::REFERENCE_MASS_SQUARED)
                .with(Atom::Zero);
            assert_same_delta(&graph, &label, &massless[support], &at_zero)?;
        }
        let expected = massive_plan.counterterm(&graph, model, &vakint, &settings)?;
        let actual = rstar_plan.counterterm(&graph, model, &vakint, &settings)?;
        assert_same_delta(&graph, &format!("{name} marked-inside root"), &actual, &expected)?;
        eprintln!(
            "RECURSIVE_MASS_DEPENDENT {name} marked={} child_loops={}",
            usize::from(marked),
            quadratic.argument.loops
        );
    }
    ensure!(fallbacks >= 2, "only {fallbacks} mass-dependent children exercised");
    Ok(())
}

#[test]
fn uv_children_route_parses_the_environment_syntax() -> Result<()> {
    assert_eq!(UvChildrenRoute::parse("")?, UvChildrenRoute::Massive);
    assert_eq!(UvChildrenRoute::parse("massive")?, UvChildrenRoute::Massive);
    assert_eq!(UvChildrenRoute::parse("rstar")?, RSTAR);
    assert_eq!(
        UvChildrenRoute::parse("rstar:3")?,
        UvChildrenRoute::Rstar { min_loops: 3 }
    );
    assert!(UvChildrenRoute::parse("rstar:x").is_err());
    assert!(UvChildrenRoute::parse("recursive").is_err());
    assert_eq!(RSTAR.to_string(), "rstar:1");
    Ok(())
}

fn edges(graph: &Graph, selected: &[usize]) -> SuBitGraph {
    let mut support = graph.empty_subgraph::<SuBitGraph>();
    for edge in selected {
        support.add(graph[&EdgeIndex::from(*edge)].1);
    }
    support
}

/// The four-loop quark self-energy insertion of the five-loop diagram 41295
/// (edges 5-15 of its saved DOT), whose massive-R child took most of that
/// diagram's 26,250 s in production. Recursive R* needs at most three-loop
/// p-integrals. Its value must not depend on the marked line: the automatic
/// line is compared with a second, explicitly marked quark line.
/// Heavy: `cargo test ... -- --ignored four_loop_quark_self_energy`.
#[test]
#[ignore = "computes a four-loop UV child of a five-loop diagram twice (minutes)"]
fn four_loop_quark_self_energy_of_41295_is_independent_of_the_marked_line() -> Result<()> {
    // Its outer Taylor polynomial exceeds the production piece limit: the
    // streamed simplification of an open Dirac chain is a different spelling
    // than the whole map's, which the test-only `report_check` rejects. Use
    // whole maps here; the streamed path is checked end to end by the zb of
    // diagram 41295 (IMPLEMENTATION.md).
    with_piece_limit(usize::MAX, four_loop_quark_self_energy_of_41295_body)
}

fn four_loop_quark_self_energy_of_41295_body() -> Result<()> {
    test_initialise()?;
    let vakint = Vakint::new()?;
    let settings = settings();
    let diagram = colour_factored(&load(include_str!("testdata/l5-d41295.dot"))?)?;
    let graph = diagram.to_gamma_loop_graph(None, true)?;
    let model = diagram.model();
    let plan = RStarPlan::with_automatic_edge(&graph, model, &settings)?;
    let insertion = edges(&graph, &[5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]);
    let call = plan
        .counterterms
        .iter()
        .find(|call| call.argument.support == insertion)
        .ok_or_else(|| eyre!("41295 has no four-loop quark self-energy child"))?;
    ensure!(call.argument.loops == 4 && call.taylor_order == 1);
    let cache = Arc::new(CountertermCache::default());
    let started = Instant::now();
    let automatic = RStarPlan::rstar_child_counterterm(
        &cache, &graph, &insertion, model, &vakint, &settings, RSTAR,
    )?;
    let automatic_seconds = started.elapsed().as_secs_f64();
    ensure!(!automatic.expand().is_zero(), "zero four-loop child");
    let child = RStarPlan::for_support(&graph, model, &insertion, &settings)?;
    let alternative = [5_usize, 6, 8, 9, 11]
        .map(EdgeIndex::from)
        .into_iter()
        .find(|edge| *edge != child.marked_edge)
        .unwrap();
    let started = Instant::now();
    let explicit =
        RStarPlan::for_support_with_edge(&graph, model, &insertion, Some(alternative), &settings)?
            .with_uv_children(RSTAR)
            .with_counterterm_cache(Arc::clone(&cache))
            .counterterm(&graph, model, &vakint, &settings)?;
    eprintln!(
        "RECURSIVE_41295 automatic_mark={} seconds={automatic_seconds:.1} alternative_mark={} seconds={:.1} cache={:?}",
        usize::from(child.marked_edge),
        usize::from(alternative),
        started.elapsed().as_secs_f64(),
        cache.stats()
    );
    assert_same_delta(&graph, "41295 four-loop child", &explicit, &automatic)
}
