//! Complete-delta reuse is checked against the independently unfolded forest.

use std::sync::Arc;

use eyre::Result;
use idenso::{CookMode, CookSettings, tensor::SymbolicTensor};
use linnet::half_edge::involution::EdgeIndex;
use spenso::{
    network::library::symbolic::ETS,
    structure::{
        representation::{Minkowski, RepName},
        slot::ParseableAind,
    },
};
use symbolica::{
    atom::{Atom, AtomCore},
    parse, symbol,
};
use vakint::Vakint;

use crate::{
    debug_tags, finalized_runtime_dot,
    graph::{Autogen, Graph, parse::IntoFinalizedRuntimeGraph},
    initialisation::test_initialise,
    numerator::aind::Aind,
    utils::{GS, load_generic_model},
    uv::{
        UltravioletGraph,
        approx::{
            UVCtx,
            integrated::Integrated,
            local_4d::{self, Full4dCts},
        },
        counterterm_identity::CountertermIdentity,
        hedge_poset::ForestNode,
    },
};

use super::{
    CountertermCache, MassiveUvEvaluator, RStarPlan,
    tests::{
        assert_same_counterterm, canonical_tensor, gluon_self_energy, parallel_pair, settings,
    },
};

fn massless(mut graph: Graph) -> Graph {
    let edges = graph
        .iter_edges()
        .map(|(_, edge, _)| edge)
        .collect::<Vec<_>>();
    for edge in edges {
        let data = &mut graph.underlying[edge];
        data.particle = data.particle.clone().override_mass(Some(Atom::Zero));
    }
    graph
}

fn nested_bubbles() -> Result<Graph> {
    Ok(massless(finalized_runtime_dot!(
        digraph G {
            graph [projector=1];
            node [num=1];
            edge [particle="scalar_1" num=1];
            e [style=invis];
            e -> A [id=6 sink="{ufo_order:0}"];
            C -> e [id=7 source="{ufo_order:0}"];
            A -> D [id=0 source="{ufo_order:1}" sink="{ufo_order:0}"];
            A -> D [id=1 lmb_id=0 source="{ufo_order:2}" sink="{ufo_order:1}"];
            D -> B [id=2 source="{ufo_order:2}" sink="{ufo_order:0}"];
            A -> B [id=3 lmb_id=1 source="{ufo_order:3}" sink="{ufo_order:1}"];
            A -> C [id=4 source="{ufo_order:4}" sink="{ufo_order:1}"];
            B -> C [id=5 lmb_id=2 source="{ufo_order:2}" sink="{ufo_order:2}"];
        }, "scalars"
    )?))
}

fn disjoint_bubbles() -> Result<Graph> {
    Ok(massless(finalized_runtime_dot!(
        digraph G {
            graph [projector=1];
            node [num=1];
            edge [particle="scalar_1" num=1];
            A -> B [id=0 source="{ufo_order:0}" sink="{ufo_order:0}"];
            A -> B [id=1 lmb_id=0 source="{ufo_order:1}" sink="{ufo_order:1}"];
            C -> D [id=2 source="{ufo_order:0}" sink="{ufo_order:0}"];
            C -> D [id=3 lmb_id=1 source="{ufo_order:1}" sink="{ufo_order:1}"];
            B -> C [id=4 source="{ufo_order:2}" sink="{ufo_order:2}"];
            D -> A [id=5 lmb_id=2 source="{ufo_order:2}" sink="{ufo_order:2}"];
        }, "scalars"
    )?))
}

fn assert_same_open_tensor(actual: &Atom, expected: &Atom) -> Result<()> {
    // Color reduction may introduce private sums with fresh names. Canonize
    // each result before subtraction: rerunning color reduction on an
    // uncollected difference can generate further alpha-equivalent summands.
    let cooking = CookSettings::indices()
        .with_mode(CookMode::ReversibleEncoding)
        .with_representation_payloads(true, true);
    let canonical = |value: &Atom| -> Result<Atom> {
        let simplified = canonical_tensor(value)?;
        let tensor =
            SymbolicTensor::infer_indexed(cooking.try_cook(simplified.as_view())?)?.canonize()?;
        Ok(cooking.uncook(tensor.expression().as_view()).expand())
    };
    let actual = canonical(actual)?;
    let expected = canonical(expected)?;
    assert!(
        (&actual - &expected).together().is_zero(),
        "cached: {actual}; forest: {expected}"
    );
    Ok(())
}

#[test]
fn cached_nested_delta_matches_unfolded_forest_and_reuses_completed_children() -> Result<()> {
    test_initialise()?;
    let graph = nested_bubbles()?;
    let model = load_generic_model("scalars");
    let settings = settings();
    let vakint = Vakint::new()?;
    let plan = RStarPlan::new(&graph, &model, EdgeIndex::from(4), &settings)?;
    let mut forest = plan.uv_forest(&graph, &settings);
    forest.integrate(&graph, &model, &vakint, &settings)?;
    let expected = forest.counterterm_for_support(&plan.support)?;
    let cache = CountertermCache::default();
    // Filling the two-loop child first makes the three-loop miss enter with a
    // completed nested delta, rather than recomputing that child's forest.
    let child = plan
        .counterterms
        .iter()
        .find(|call| call.argument.loops == 2)
        .unwrap();
    let _ = cache.counterterm(&graph, &child.argument.support, &model, &vakint, &settings)?;
    let before = cache.stats();
    let actual = cache.counterterm(&graph, &plan.support, &model, &vakint, &settings)?;
    assert_same_counterterm(&actual, &expected)?;
    let after = cache.stats();
    assert!(after.hits > before.hits);
    assert_eq!(after.misses, before.misses + 1);
    let repeated = cache.counterterm(&graph, &plan.support, &model, &vakint, &settings)?;
    assert_same_counterterm(&repeated, &expected)?;
    assert_eq!(cache.stats().misses, after.misses);
    Ok(())
}

#[test]
fn pole_depth_preserves_finite_quotient_terms_promoted_by_child_poles() -> Result<()> {
    test_initialise()?;
    let graph = nested_bubbles()?;
    let model = load_generic_model("scalars");
    let settings = settings();
    let vakint = Vakint::new()?;
    let plan = RStarPlan::new(&graph, &model, EdgeIndex::from(4), &settings)?;
    let catalogue =
        graph.classified_spinneys(&plan.support, &model, &settings, &graph.loop_momentum_basis);
    let current = ForestNode {
        spinney: catalogue
            .iter()
            .find(|spinney| spinney.filter() == &plan.support)
            .unwrap(),
        topo_order: 0,
    };
    let given = ForestNode {
        spinney: catalogue
            .iter()
            .find(|spinney| graph.n_loops(spinney.filter()) == 1)
            .unwrap(),
        topo_order: 0,
    };
    assert_eq!(graph.n_loops(&plan.support), 3);
    let epsilon = Atom::var(GS.dim_epsilon);
    let input = Full4dCts::from_counterterm_product(epsilon.pow(-1), 1);
    let ctx = UVCtx::new(&graph, &model, &settings);
    let local = local_4d::uv_limit(&input, &ctx, &current, &given, &current, &given)?;
    let mut legacy_settings = settings.vakint.true_settings();
    legacy_settings.number_of_terms_in_epsilon_expansion = 4;
    let expected = Integrated::new(&vakint, &legacy_settings)
        .run(&local, &ctx, &current, &given, &current, &given)?
        .physical_pole_atom();
    let mut pole_settings = legacy_settings.clone();
    // Two remaining loops mean N=2 requests the complete product through
    // epsilon^-1. The one-loop insertion still requires the vacuum's finite
    // coefficient; truncating that vacuum separately to poles would lose it.
    pole_settings.number_of_terms_in_epsilon_expansion = 2;
    let actual = Integrated::new(&vakint, &pole_settings)
        .run(&local, &ctx, &current, &given, &current, &given)?
        .physical_pole_atom();
    assert_same_counterterm(&actual, &expected)?;
    assert!(
        expected
            .coefficient_list::<i32>(std::slice::from_ref(&epsilon))
            .iter()
            .any(|(power, coefficient)| *power == epsilon.pow(-1) && !coefficient.is_zero()),
        "the child pole must promote a nonzero finite vacuum coefficient"
    );
    Ok(())
}

#[test]
fn pole_depth_combines_child_and_dimension_poles_with_open_tensor_ports() -> Result<()> {
    test_initialise()?;
    let graph = nested_bubbles()?;
    let model = load_generic_model("scalars");
    let settings = settings();
    let vakint = Vakint::new()?;
    let plan = RStarPlan::new(&graph, &model, EdgeIndex::from(4), &settings)?;
    let catalogue =
        graph.classified_spinneys(&plan.support, &model, &settings, &graph.loop_momentum_basis);
    let current = ForestNode {
        spinney: catalogue
            .iter()
            .find(|spinney| spinney.filter() == &plan.support)
            .unwrap(),
        topo_order: 0,
    };
    let given = ForestNode {
        spinney: catalogue
            .iter()
            .find(|spinney| graph.n_loops(spinney.filter()) == 2)
            .unwrap(),
        topo_order: 0,
    };
    assert_eq!(graph.n_loops(&plan.support), 3);
    let epsilon = Atom::var(GS.dim_epsilon);
    let mink = Minkowski {}.new_rep(GS.dim);
    let metric = ETS.metric(
        mink.to_symbolic([Aind::Hedge(6, 1).to_atom()]),
        mink.to_symbolic([Aind::Hedge(7, 1).to_atom()]),
    );
    let poles = metric * epsilon.pow(-2) / (Atom::var(GS.dim) - Atom::num(4));
    let input = Full4dCts::from_counterterm_product(poles, 2);
    let ctx = UVCtx::new(&graph, &model, &settings);
    let local = local_4d::uv_limit(&input, &ctx, &current, &given, &current, &given)?;
    let mut deep_settings = settings.vakint.true_settings();
    // N=3 retains the complete one-loop quotient through epsilon^1, two
    // orders beyond its poles. The combined prefactor is -1/(2 epsilon^3).
    deep_settings.number_of_terms_in_epsilon_expansion = 3;
    let expected = Integrated::new(&vakint, &deep_settings)
        .run(&local, &ctx, &current, &given, &current, &given)?
        .physical_pole_atom();
    let mut pole_settings = deep_settings.clone();
    pole_settings.number_of_terms_in_epsilon_expansion = 1;
    let actual = Integrated::new(&vakint, &pole_settings)
        .run(&local, &ctx, &current, &given, &current, &given)?
        .physical_pole_atom();
    assert_same_open_tensor(&actual, &expected)?;

    let cooking = CookSettings::indices()
        .with_mode(CookMode::ReversibleEncoding)
        .with_representation_payloads(true, true);
    let tensor = SymbolicTensor::infer_indexed(cooking.try_cook(actual.as_view())?)?;
    assert_eq!(tensor.rank(), 2, "both open Lorentz ports must survive");
    let coefficients = actual.coefficient_list::<i32>(std::slice::from_ref(&epsilon));
    let simple_pole = coefficients
        .iter()
        .find(|(power, _)| *power == epsilon.pow(-1))
        .expect("the combined poles promote a nonzero epsilon^2 master coefficient");
    // The one-loop master's first zeta(3) occurs at epsilon^2. Its presence
    // prevents this comparison from passing with only lower master orders.
    assert!(
        simple_pole
            .1
            .pattern_match(&parse!("symbolica::zeta(3)").to_pattern(), None, None)
            .next()
            .is_some(),
        "the simple pole must retain the promoted zeta(3) coefficient"
    );
    Ok(())
}

#[test]
fn cached_dimension_pole_retains_the_promoted_finite_master_coefficient() -> Result<()> {
    test_initialise()?;
    let mut graph = massless(finalized_runtime_dot!(
        digraph G {
            graph [projector=1];
            node [num=1];
            edge [particle="scalar_1" num=1];
            e [style=invis];
            e -> A [id=2 sink="{ufo_order:0}"];
            B -> e [id=3 source="{ufo_order:0}"];
            A -> B [id=0 source="{ufo_order:1}" sink="{ufo_order:1}"];
            A -> B [id=1 lmb_id=0 source="{ufo_order:2}" sink="{ufo_order:2}"];
        }, "scalars"
    )?);
    let model = load_generic_model("scalars");
    let settings = settings();
    let vakint = Vakint::new()?;
    let vertex = graph
        .iter_nodes()
        .find_map(|(node, _, data)| (data.name.value == "A").then_some(node))
        .unwrap();
    let epsilon = Atom::var(GS.dim_epsilon);
    // This independent expression already exposes its pole to Vakint, even
    // before the scalar-d regularization used by the cache was introduced.
    graph.underlying[vertex].num = Autogen::explicit(-Atom::one() / (Atom::num(2) * &epsilon));
    let explicit_plan = RStarPlan::new(&graph, &model, EdgeIndex::from(0), &settings)?;
    let mut explicit_forest = explicit_plan.uv_forest(&graph, &settings);
    explicit_forest.integrate(&graph, &model, &vakint, &settings)?;
    let expected = explicit_forest.counterterm_for_support(&explicit_plan.support)?;
    assert!(
        expected
            .coefficient_list::<i32>(std::slice::from_ref(&epsilon))
            .iter()
            .any(|(power, coefficient)| *power == epsilon.pow(-1) && !coefficient.is_zero()),
        "the finite master coefficient must contribute to the simple pole"
    );

    graph.underlying[vertex].num = Autogen::explicit((Atom::var(GS.dim) - Atom::num(4)).pow(-1));
    let plan = RStarPlan::new(&graph, &model, EdgeIndex::from(0), &settings)?;
    let actual = CountertermCache::default().counterterm(
        &graph,
        &plan.support,
        &model,
        &vakint,
        &settings,
    )?;
    assert_same_counterterm(&actual, &expected)?;
    let mut forest = plan.uv_forest(&graph, &settings);
    forest.integrate(&graph, &model, &vakint, &settings)?;
    assert_same_counterterm(&actual, &forest.counterterm_for_support(&plan.support)?)
}

#[test]
fn cached_disjoint_deltas_match_forest_and_share_particle_isomorphism() -> Result<()> {
    test_initialise()?;
    let graph = disjoint_bubbles()?;
    let model = load_generic_model("scalars");
    let settings = settings();
    let vakint = Vakint::new()?;
    let cache = Arc::new(CountertermCache::default());
    let plan = RStarPlan::new(&graph, &model, EdgeIndex::from(4), &settings)?;
    let mass = Atom::var(symbol!("cached_reference_mass"; Scalar));
    let expected = plan.evaluate_uv_children(&graph, &model, &mass, &vakint, &settings)?;
    let cached = plan.with_counterterm_cache(cache.clone());
    let actual = cached.evaluate_uv_children(&graph, &model, &mass, &vakint, &settings)?;
    assert_eq!(actual.len(), 2);
    for (support, value) in &expected {
        assert_same_counterterm(&actual[support], value)?;
    }
    assert_eq!(
        cache.stats().misses,
        1,
        "the two scalar bubbles are isomorphic"
    );
    assert_eq!(cache.stats().hits, 1);
    let mut forest = cached.uv_forest(&graph, &settings);
    forest.integrate(&graph, &model, &vakint, &settings)?;
    let actual = cache.counterterm(&graph, &cached.support, &model, &vakint, &settings)?;
    assert_same_counterterm(&actual, &forest.counterterm_for_support(&cached.support)?)
}

#[test]
fn cached_local_epsilon_coefficient_stays_inside_each_pole_projection() -> Result<()> {
    test_initialise()?;
    let mut graph = nested_bubbles()?;
    let model = load_generic_model("scalars");
    let settings = settings();
    let vakint = Vakint::new()?;
    let cache = CountertermCache::default();
    let plan = RStarPlan::new(&graph, &model, EdgeIndex::from(4), &settings)?;
    let plain = cache.counterterm(&graph, &plan.support, &model, &vakint, &settings)?;
    let vertex = graph
        .iter_nodes()
        .find_map(|(node, _, data)| (data.name.value == "D").then_some(node))
        .unwrap();
    graph.underlying[vertex].num = Autogen::explicit(Atom::num(1) + Atom::var(GS.dim_epsilon));
    let dressed_plan = RStarPlan::new(&graph, &model, EdgeIndex::from(4), &settings)?;
    let mut forest = dressed_plan.uv_forest(&graph, &settings);
    forest.integrate(&graph, &model, &vakint, &settings)?;
    let expected = forest.counterterm_for_support(&dressed_plan.support)?;
    let actual = cache.counterterm(&graph, &dressed_plan.support, &model, &vakint, &settings)?;
    assert_same_counterterm(&actual, &expected)?;
    assert!(!(&actual - &plain).expand().is_zero());
    Ok(())
}

#[test]
fn cached_decorated_quotients_distinguish_insertions_and_reference_masses() -> Result<()> {
    test_initialise()?;
    let graph = nested_bubbles()?;
    let model = load_generic_model("scalars");
    let settings = settings();
    let vakint = Vakint::new()?;
    let plan = RStarPlan::new(&graph, &model, EdgeIndex::from(0), &settings)?;
    let call = &plan.ir[0].forest[0].steps[0];
    assert_eq!(call.argument.insertions.len(), 1);
    let cache = Arc::new(CountertermCache::default());
    for coefficient in [Atom::num(1), Atom::num(1) + Atom::var(GS.dim_epsilon)] {
        for mass in [
            Atom::Zero,
            Atom::var(symbol!("cached_quotient_mass"; Scalar)),
        ] {
            let edge = graph
                .iter_edges_of(&call.argument.remaining)
                .next()
                .unwrap()
                .1;
            let ordinary = MassiveUvEvaluator::new(&graph, &model, &vakint, &settings)
                .with_reference_mass(edge, mass.clone());
            let expected = ordinary.counterterm(call, std::slice::from_ref(&coefficient))?;
            let cached = ordinary.with_counterterm_cache(cache.clone());
            let actual = cached.counterterm(call, std::slice::from_ref(&coefficient))?;
            assert_same_counterterm(&actual, &expected)?;
            let before = cache.stats();
            let repeated = cached.counterterm(call, std::slice::from_ref(&coefficient))?;
            assert_same_counterterm(&repeated, &expected)?;
            assert_eq!(cache.stats().misses, before.misses);
        }
    }
    assert!(cache.stats().misses >= 4);
    Ok(())
}

#[test]
fn cached_gluon_delta_preserves_open_tensor_and_external_momenta_between_parents() -> Result<()> {
    test_initialise()?;
    let model = load_generic_model("sm");
    let settings = settings();
    let vakint = Vakint::new()?;
    let cache = CountertermCache::default();
    let primitive = gluon_self_energy(&model, 1, 2, "V_36", |_| true)?;
    let primitive_mark = primitive
        .iter_edges()
        .find_map(|(pair, edge, _)| pair.is_paired().then_some(edge))
        .unwrap();
    let primitive_plan = RStarPlan::new(&primitive, &model, primitive_mark, &settings)?;
    let mut primitive_forest = primitive_plan.uv_forest(&primitive, &settings);
    primitive_forest.integrate(&primitive, &model, &vakint, &settings)?;
    let primitive_delta = cache.counterterm(
        &primitive,
        &primitive_plan.support,
        &model,
        &vakint,
        &settings,
    )?;
    assert_same_open_tensor(
        &primitive_delta,
        &primitive_forest.counterterm_for_support(&primitive_plan.support)?,
    )?;

    let parent = gluon_self_energy(&model, 2, 4, "V_36", |graph| parallel_pair(graph).is_some())?;
    let pair = parallel_pair(&parent).unwrap();
    let plan = RStarPlan::new(&parent, &model, pair[0], &settings)?;
    let mut forest = plan.uv_forest(&parent, &settings);
    forest.integrate(&parent, &model, &vakint, &settings)?;
    let before = cache.stats();
    for call in &plan.counterterms {
        let actual =
            cache.counterterm(&parent, &call.argument.support, &model, &vakint, &settings)?;
        assert_same_open_tensor(
            &actual,
            &forest.counterterm_for_support(&call.argument.support)?,
        )?;
    }
    assert!(
        cache.stats().hits > before.hits,
        "the gluon bubble must reuse the earlier parent with its open ports remapped"
    );
    Ok(())
}

#[test]
fn gluon_isomorphic_supports_have_identical_cache_keys() -> Result<()> {
    test_initialise()?;
    let model = load_generic_model("sm");
    let settings = settings();
    let primitive = gluon_self_energy(&model, 1, 2, "V_36", |_| true)?;
    let primitive_mark = primitive
        .iter_edges()
        .find_map(|(pair, edge, _)| pair.is_paired().then_some(edge))
        .unwrap();
    let primitive_plan = RStarPlan::new(&primitive, &model, primitive_mark, &settings)?;
    let primitive_key = CountertermIdentity::new(&primitive, &primitive_plan.support, &model)?.key;
    debug_tags!(#uv, #cache, #trace;
        stage = "primitive_gluon_cache_key",
        file.graph = %primitive.debug_dot(),
        "Primitive gluon bubble identity fixture"
    );
    let parent = gluon_self_energy(&model, 2, 4, "V_36", |graph| parallel_pair(graph).is_some())?;
    let plan = RStarPlan::new(
        &parent,
        &model,
        parallel_pair(&parent).unwrap()[0],
        &settings,
    )?;
    let mut matched = false;
    for call in &plan.counterterms {
        let key = CountertermIdentity::new(&parent, &call.argument.support, &model)?.key;
        debug_tags!(#uv, #cache, #trace;
            stage = "child_gluon_cache_key",
            file.graph = %parent.debug_dot(),
            support = ?call.argument.support,
            matches_primitive = key == primitive_key,
            "Parent's proper UV support identity fixture"
        );
        matched |= key == primitive_key;
    }
    assert!(
        matched,
        "the primitive gluon bubble must match its parent occurrence"
    );
    Ok(())
}

#[test]
fn static_delta_prewarm_deduplicates_and_completes_nested_levels() -> Result<()> {
    test_initialise()?;
    let graph = Arc::new(nested_bubbles()?);
    let model = load_generic_model("scalars");
    let settings = Arc::new(settings());
    let vakint = Vakint::new()?;
    let plan = RStarPlan::new(&graph, &model, EdgeIndex::from(4), &settings)?;
    let request = super::CountertermRequest::new(
        Arc::clone(&graph),
        plan.support.clone(),
        &model,
        Arc::clone(&settings),
    );
    let cache = CountertermCache::default().with_timing(true);
    let discovery = super::CountertermPrewarmPlan::default();
    cache.extend_prewarm_plan(&discovery, [request.clone()])?;
    let original_counts = discovery.unique_requests_by_loop();
    assert_eq!(discovery.retained_graph_families(), 1);
    // A duplicate discovered in a later batch must not retain a second family.
    let duplicate_family = Arc::new((*graph).clone());
    let dropped_duplicate = Arc::downgrade(&duplicate_family);
    let duplicate = super::CountertermRequest::new(
        duplicate_family,
        plan.support.clone(),
        &model,
        Arc::clone(&settings),
    );
    cache.extend_prewarm_plan(&discovery, [duplicate, request.clone()])?;
    assert!(dropped_duplicate.upgrade().is_none());
    assert_eq!(discovery.unique_requests_by_loop(), original_counts);
    assert_eq!(discovery.retained_graph_families(), 1);
    let levels = discovery.into_levels();
    let unique = levels.values().map(Vec::len).sum::<usize>();
    assert_eq!(levels.keys().copied().collect::<Vec<_>>(), vec![1, 2, 3]);
    assert_eq!(levels[&3].len(), 1);
    for (loops, requests) in levels {
        for request in requests {
            assert_eq!(request.loops(), loops);
            let _ = request.evaluate(&cache, &vakint)?;
        }
    }
    // A missing recursively discovered dependency would cause an extra miss.
    assert_eq!(cache.stats().misses, unique as u64);
    let actual = request.evaluate(&cache, &vakint)?;
    assert_eq!(cache.stats().misses, unique as u64);
    assert_eq!(cache.stats().waits, 0);
    assert!(cache.prewarm_levels([request])?.is_empty());
    let mut forest = plan.uv_forest(&graph, &settings);
    forest.integrate(&graph, &model, &vakint, &settings)?;
    assert_same_counterterm(&actual, &forest.counterterm_for_support(&plan.support)?)
}

#[test]
fn static_delta_prewarm_prunes_children_of_overall_convergent_requests() -> Result<()> {
    test_initialise()?;
    let model = load_generic_model("scalars");
    // Existing bubble-in-kite fixture: a convergent three-loop parent with
    // one logarithmically divergent child. Its complete local Delta is zero.
    let graph = Arc::new(massless(finalized_runtime_dot!(
        digraph G {
            graph [projector=1]; node [num=1]; edge [particle="scalar_1" num=1];
            e [style=invis];
            e -> A:0 [id=7 sink="{ufo_order:0}"];
            B:1 -> e [id=8 source="{ufo_order:0}"];
            C -> E [id=0 source="{ufo_order:0}" sink="{ufo_order:0}"];
            C -> E [id=1 lmb_id=0 source="{ufo_order:1}" sink="{ufo_order:1}"];
            A -> C [id=2 source="{ufo_order:1}" sink="{ufo_order:2}"];
            A -> D [id=3 source="{ufo_order:2}" sink="{ufo_order:0}"];
            B -> E [id=4 source="{ufo_order:1}" sink="{ufo_order:2}"];
            B -> D [id=5 lmb_id=1 source="{ufo_order:2}" sink="{ufo_order:1}"];
            C -> D [id=6 lmb_id=2 source="{ufo_order:3}" sink="{ufo_order:2}"];
        }, &model
    )?));
    let settings = Arc::new(settings());
    let plan = RStarPlan::new(&graph, &model, EdgeIndex::from(0), &settings)?;
    assert_eq!(graph.compute_dod(&plan.support), -2);
    assert_eq!(plan.counterterms.len(), 1);
    let request =
        super::CountertermRequest::new(Arc::clone(&graph), plan.support.clone(), &model, settings);
    let cache = CountertermCache::default();
    let levels = cache.prewarm_levels([request])?;
    assert_eq!(levels.keys().copied().collect::<Vec<_>>(), vec![3]);
    assert_eq!(levels[&3].len(), 1);
    assert!(levels[&3][0].evaluate(&cache, &Vakint::new()?)?.is_zero());
    assert_eq!(cache.stats().misses, 1);
    Ok(())
}

#[test]
fn prewarmed_original_uv_children_keep_the_exact_reference_mass_and_open_ports() -> Result<()> {
    test_initialise()?;
    let model = load_generic_model("sm");
    let graph = gluon_self_energy(&model, 2, 4, "V_36", |graph| parallel_pair(graph).is_some())?;
    let settings = settings();
    let vakint = Vakint::new()?;
    let mark = parallel_pair(&graph).unwrap()[0];
    let cache = Arc::new(CountertermCache::default());
    let plain = RStarPlan::new(&graph, &model, mark, &settings)?;
    let reference_mass = Atom::var(*super::REFERENCE_MASS_SQUARED).sqrt();
    let expected =
        plain.evaluate_uv_children(&graph, &model, &reference_mass, &vakint, &settings)?;
    let cached = plain.with_counterterm_cache(Arc::clone(&cache));
    for requests in cache
        .prewarm_levels(cached.uv_child_requests(&graph, &model, &settings))?
        .into_values()
    {
        for request in requests {
            let _ = request.evaluate(&cache, &vakint)?;
        }
    }
    let before = cache.stats();
    let actual =
        cached.evaluate_uv_children(&graph, &model, &reference_mass, &vakint, &settings)?;
    assert_eq!(cache.stats().misses, before.misses);
    assert_eq!(actual.len(), expected.len());
    for (support, tensor) in expected {
        assert_same_open_tensor(&actual[&support], &tensor)?;
    }
    Ok(())
}

#[test]
fn parallel_delta_discovery_matches_serial_closure_and_exact_results() -> Result<()> {
    test_initialise()?;
    let graph = Arc::new(nested_bubbles()?);
    let model = load_generic_model("scalars");
    let settings = Arc::new(settings());
    let vakint = Vakint::new()?;
    let plan = RStarPlan::new(&graph, &model, EdgeIndex::from(4), &settings)?;
    let child = plan
        .counterterms
        .iter()
        .find(|call| call.argument.loops == 2)
        .unwrap();
    let root = super::CountertermRequest::new(
        Arc::clone(&graph),
        plan.support.clone(),
        &model,
        Arc::clone(&settings),
    );
    let child_request = super::CountertermRequest::new(
        Arc::clone(&graph),
        child.argument.support.clone(),
        &model,
        Arc::clone(&settings),
    );
    let serial_cache = CountertermCache::default();
    let serial = super::CountertermPrewarmPlan::default();
    serial_cache.extend_prewarm_plan(&serial, [root.clone(), child_request])?;
    let expected_counts = serial.unique_requests_by_loop();
    let parallel_cache = CountertermCache::default();
    let parallel = super::CountertermPrewarmPlan::default();
    // Distinct parent allocations race for duplicate roots and their shared
    // child. No worker waits for another worker's recursive classification.
    let start = std::sync::Barrier::new(4);
    std::thread::scope(|scope| -> Result<()> {
        let mut handles = Vec::new();
        for index in 0..4 {
            let support = if index % 2 == 0 {
                &plan.support
            } else {
                &child.argument.support
            };
            let request = super::CountertermRequest::new(
                Arc::new((*graph).clone()),
                support.clone(),
                &model,
                Arc::clone(&settings),
            );
            let cache = &parallel_cache;
            let discovery = &parallel;
            let start = &start;
            handles.push(scope.spawn(move || {
                start.wait();
                cache.extend_prewarm_plan(discovery, [request])
            }));
        }
        for handle in handles {
            handle.join().expect("parallel discovery worker panicked")?;
        }
        Ok(())
    })?;
    assert_eq!(parallel.unique_requests_by_loop(), expected_counts);
    let unique = expected_counts.values().sum::<usize>();
    for requests in serial.into_levels().into_values() {
        for request in requests {
            let _ = request.evaluate(&serial_cache, &vakint)?;
        }
    }
    for requests in parallel.into_levels().into_values() {
        for request in requests {
            let _ = request.evaluate(&parallel_cache, &vakint)?;
        }
    }
    assert_eq!(serial_cache.stats().misses, unique as u64);
    assert_eq!(parallel_cache.stats().misses, unique as u64);
    assert_eq!(parallel_cache.stats().waits, 0);
    assert_same_counterterm(
        &root.evaluate(&parallel_cache, &vakint)?,
        &root.evaluate(&serial_cache, &vakint)?,
    )?;
    assert!(parallel_cache.prewarm_levels([root])?.is_empty());
    Ok(())
}
