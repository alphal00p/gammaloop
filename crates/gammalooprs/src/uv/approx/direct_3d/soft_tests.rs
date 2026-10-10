use super::super::branches::{DirectResidueBranches, DirectResidueKey};
use super::*;
use crate::cff::expression::OrientationID;
use crate::{
    dot,
    graph::parse::IntoGraph,
    graph::{FeynmanGraph, Graph, cuts::CutSet},
    initialisation::test_initialise,
    settings::global::OrientationPattern,
    utils::GS,
    uv::{
        Integrands, UVgenerationSettings,
        approx::{
            OrientationProjection, Rooted,
            direct_3d::{Direct3dCts, DirectSector},
            local_3d::{
                Local3DConceptualBranch, Local3DMaterialization, Local3DProjectionPath,
                Local3DSignedBranch, Localizer,
            },
        },
    },
};
use gammaloop_tracing_filter::LogMessage;
use idenso::shorthands::metric::MetricSimplifier;
use linnet::half_edge::{
    involution::EdgeIndex,
    subgraph::{ModifySubSet, SubSetOps},
};
use spenso::{
    network::tags::SPENSO_TAG,
    structure::representation::{Minkowski, RepName},
};
use symbolica::{function, symbol};

// #[debug_instrument(
//     current = %current.log_display(),
//     given = %given.log_display(),
//     reduced,
// )]
#[cfg(test)]
fn t_tilde<S: ForestNodeLike>(
    ctx: &UVCtx<'_>,
    current: &S,
    given: &S,
    started: &Atom,
    active_subgraph: Option<&SuBitGraph>,
    lmb: &LoopMomentumBasis,
) -> Result<Atom> {
    if current.dod() == 0 {
        return Ok(Atom::Zero);
    }
    Local3DLoopRescaling::FullSubgraph.project(
        Local3DDeformation::Soft,
        ctx,
        current,
        given,
        started,
        active_subgraph,
        lmb,
    )
}

struct TestNode {
    subgraph: SuBitGraph,
    lmb: LoopMomentumBasis,
    dod: i32,
    scheme: ApproximationType,
}

impl LogMessage for TestNode {
    fn log_display(&self) -> String {
        "local-3d test node".to_owned()
    }
}

impl ForestNodeLike for TestNode {
    fn subgraph(&self) -> &SuBitGraph {
        &self.subgraph
    }

    fn lmb(&self) -> &LoopMomentumBasis {
        &self.lmb
    }

    fn dod(&self) -> i32 {
        self.dod
    }

    fn renormalization_scheme(&self) -> ApproximationType {
        self.scheme
    }

    fn topo_order(&self) -> usize {
        0
    }

    fn reduced_subgraph(&self, given: &Self) -> SuBitGraph {
        self.subgraph.subtract(&given.subgraph)
    }
}

fn scalar_two_point_graph() -> Graph {
    dot!(
        digraph scalar_self_energy {
            edge [particle="H" num=1];
            node [num=1];
            ext [style=invis];
            ext -> A:0 [id=0];
            B:1 -> ext [id=1];
            A -> B [id=2];
            A -> B [id=3];
        }
    )
    .unwrap()
}

fn two_point_node(graph: &Graph, dod: i32) -> TestNode {
    let external: SuBitGraph = graph.external_filter();
    let subgraph = graph.full_filter().subtract(&external);
    TestNode {
        lmb: graph.lmb_of(&subgraph),
        subgraph,
        dod,
        scheme: ApproximationType::IR,
    }
}

fn root_node(graph: &Graph) -> TestNode {
    TestNode {
        subgraph: graph.empty_subgraph(),
        lmb: graph.empty_lmb(),
        dod: 0,
        scheme: ApproximationType::MUV,
    }
}

#[test]
fn zero_contributions_do_not_publish_projection_paths() -> Result<()> {
    test_initialise()?;
    let graph = scalar_two_point_graph();
    let nonzero = Integrands::root();
    let zero = nonzero.map(|_| Atom::Zero);
    let make_sector = |integrands: &Integrands| -> Result<DirectSector> {
        Ok(DirectSector {
            active_subgraph: graph.empty_subgraph(),
            coordinate_frames: Vec::new(),
            projection_paths: vec![Local3DProjectionPath::default()],
            active: DirectResidueBranches::production(OrientationID(0), integrands.clone())?,
            frozen_integrands: nonzero.clone(),
        })
    };
    let cts = Direct3dCts::from_sectors(vec![make_sector(&zero)?, make_sector(&nonzero)?])?;
    assert_eq!(cts.branches()?.materialize(false)?, nonzero);
    assert_eq!(
        cts.sectors()?.len(),
        2,
        "filtering provenance must not remove a numerical sector"
    );
    assert_eq!(
        cts.projection_paths(),
        vec![Local3DProjectionPath::default()],
        "only the nonzero sector contributes an exported projection path"
    );
    assert_eq!(
        cts.projection_paths().len(),
        1,
        "a zero integrated addend must not publish a no-op path"
    );
    let with_nonzero_path =
        Direct3dCts::from_sectors(vec![make_sector(&nonzero)?, make_sector(&nonzero)?])?;
    assert_eq!(
        with_nonzero_path.branches()?.materialize(false)?,
        nonzero.clone().zip_add([nonzero.clone()])?
    );
    assert_eq!(
        with_nonzero_path.projection_paths().len(),
        2,
        "a nonzero integrated addend must retain its projection path"
    );
    Ok(())
}

#[test]
fn orientation_term_keeps_external_selectors_and_adds_internal_ones() -> Result<()> {
    test_initialise()?;
    let reduced_expression = function!(GS.ose, 0);
    // Complete generalized residue-map selectors now own both the external
    // and integrated internal energy choices. Physical edge-sign selectors
    // cannot distinguish maps which share those directions.
    let selector = OrientationID(3);
    let integrands = Integrands::root().map(|_| reduced_expression.clone());
    let branches = DirectResidueBranches::production(selector, integrands)?;
    let localized = branches.materialize(true)?;
    let expected = reduced_expression * selector.atom();
    assert_eq!(localized.iter().next().unwrap().1, &expected);
    Ok(())
}

#[test]
fn one_orientation_cff_hard_charts_match_direct_external_jets() -> Result<()> {
    test_initialise().unwrap();
    let mut graph = scalar_two_point_graph();
    let current = two_point_node(&graph, 2);
    let root = root_node(&graph);
    let orientation_pattern = OrientationPattern::default();
    let cutset = CutSet::empty(1);

    let options = graph.denominator_only_cff_3d_expression_options();
    let canonization = graph.get_esurface_canonization(&graph.loop_momentum_basis);
    let mut production = graph.generate_3d_expression_for_integrand(
        &[],
        &canonization,
        &options,
        Some(&Atom::one()),
    )?;
    production.expression.orientations.truncate(1);
    let orientation =
        OrientationProjection::exact_expression(&production, &options, &orientation_pattern, false);
    let localizer = Localizer::new(&cutset, orientation);
    let input = Direct3dCts::root(&graph, localizer)?;
    let branches = input.branches()?;
    let mut selected_terms = branches.iter_keys();
    let (_residue_key, residue_integrands) = selected_terms
        .next()
        .expect("the selected scalar-bubble orientation has a CFF term");
    assert!(
        selected_terms.next().is_none(),
        "the admitted-orientation CFF must contain exactly one orientation term"
    );
    let oriented_cff = residue_integrands.iter().next().unwrap().1;
    let settings = UVgenerationSettings::default();
    let ctx = UVCtx::new(&graph, &settings);
    let active_started = start(
        &ctx,
        &current,
        oriented_cff,
        &graph
            .numerator(&current.reduced_subgraph(&root), root.subgraph())
            .get_single_atom()
            .unwrap(),
        None,
        current.lmb(),
    )?;
    assert!(
        active_started.contains_symbol(GS.on_shell_energy)
            && (0..=3).all(|component| {
                (0..graph.underlying.n_edges()).any(|edge| {
                    active_started.contains(GS.emr_mom(EdgeIndex(edge), GS.cind(component)))
                })
            }),
        "the active energies must retain explicit spatial components and external energies"
    );

    let loop_edge = *current
        .lmb
        .loop_edges
        .first()
        .expect("the scalar bubble has a loop carrier");
    let mass_scope = function!(*LOCAL_3D_MASS_SCOPE, usize::from(loop_edge) as i64);
    let physical_mass = graph
        .iter_edges_of(&current.subgraph)
        .find(|(pair, _, _)| pair.is_paired())
        .expect("the scalar bubble has an internal propagator")
        .2
        .data
        .mass_atom()
        * &mass_scope;
    let minkowski = Minkowski {}.new_rep(4).to_symbolic([]);
    let loop_momentum = GS.emr_vec(loop_edge, minkowski.clone());
    let homogeneous_numerator = physical_mass.pow(2)
        - function!(SPENSO_TAG.dot, loop_momentum.clone(), loop_momentum.clone());
    let tree_denominator = GS.wrap_tree_denoms(
        graph.denominator(&graph.tree_edges.subtract(&graph.initial_state_cut), |_| -1),
    );
    // Expand the bilinear dot shorthand before asking Symbolica for the
    // independent direct jet.  Otherwise it treats `dot` as an opaque
    // function and leaves formally equivalent derivatives in different
    // normal forms in the direct and hard charts.
    let active_with_numerator = (&active_started * &homogeneous_numerator).expand_dots()?;
    let started = &active_with_numerator * &tree_denominator;

    let reduced = current.reduced_subgraph(&root);
    let mut momentum_replacements =
        graph.uv_spatial_wrapped_replacement(&reduced, current.lmb(), &[W_.x___]);
    momentum_replacements.extend(graph.uv_wrapped_replacement(&reduced, current.lmb(), &[W_.x___]));
    let routed_active = active_with_numerator.replace_multiple(&momentum_replacements);
    let direct_parameter = symbol!("gammalooprs::local_3d_direct_jet");
    let direct_scale = Atom::var(direct_parameter);

    for (deformation, direct_order) in [
        (Local3DDeformation::Ordinary, 2),
        (Local3DDeformation::Soft, 1),
    ] {
        let mut direct = routed_active.clone();
        for external_edge in &current.lmb.ext_edges {
            let external_spatial = GS.emr_vec(*external_edge, minkowski.clone());
            direct = direct
                .replace(external_spatial.clone())
                .with(external_spatial * &direct_scale);
            let external_momentum =
                function!(GS.emr_mom, usize::from(*external_edge) as i64, W_.x___);
            direct = direct
                .replace(external_momentum.clone())
                .with(external_momentum * &direct_scale);
        }

        if matches!(deformation, Local3DDeformation::Ordinary) {
            direct = direct
                .replace(mass_scope.clone())
                .with(mass_scope.clone() * &direct_scale);
            let uv_expansion_mass_squared = (Atom::var(GS.m_uv_expansion) * &mass_scope).pow(2);
            let uv_vacuum_mass_squared = (Atom::var(GS.m_uv_vacuum) * &mass_scope).pow(2);
            let uv_completion =
                uv_vacuum_mass_squared - direct_scale.pow(2) * uv_expansion_mass_squared;
            for energy_edge in &current.lmb.loop_edges {
                direct = direct
                    .replace(function!(
                        GS.on_shell_energy,
                        usize::from(*energy_edge) as i64,
                        W_.prop_
                    ))
                    .with(function!(
                        GS.on_shell_energy,
                        usize::from(*energy_edge) as i64,
                        W_.prop_ + &uv_completion
                    ));
            }
        }

        let direct = direct
            .series(direct_parameter, Atom::Zero, direct_order)?
            .to_atom()
            .replace(direct_parameter)
            .with(Atom::one())
            * &tree_denominator;
        let hard = Local3DLoopRescaling::FullSubgraph.project(
            deformation,
            &ctx,
            &current,
            &root,
            &started,
            None,
            current.lmb(),
        )?;
        let difference = (hard - direct).together().cancel().expand();
        assert!(
            difference.is_zero(),
            "the {deformation:?} hard chart disagrees with its direct external jet: {difference}"
        );
    }

    Ok(())
}

#[test]
fn soft_hard_laurent_projection_matches_external_jets_and_keeps_mass() -> Result<()> {
    test_initialise().unwrap();
    let graph = scalar_two_point_graph();
    let root = root_node(&graph);
    let settings = UVgenerationSettings::default();
    let ctx = UVCtx::new(&graph, &settings);

    for dod in 0..=2 {
        let current = two_point_node(&graph, dod);
        let loop_edge = *current
            .lmb
            .loop_edges
            .first()
            .expect("the scalar bubble has a loop carrier");
        let external_edge = *current
            .lmb
            .ext_edges
            .first()
            .expect("the scalar bubble has an external carrier");
        let physical_mass = graph
            .iter_edges_of(&current.subgraph)
            .find(|(pair, _, _)| pair.is_paired())
            .expect("the scalar bubble has an internal propagator")
            .2
            .data
            .mass_atom();
        let minkowski = Minkowski {}.new_rep(4).to_symbolic([]);
        let mass_scope = function!(*LOCAL_3D_MASS_SCOPE, usize::from(loop_edge) as i64);
        let scoped_mass = &physical_mass * &mass_scope;
        let loop_momentum = GS.emr_vec(loop_edge, minkowski.clone());
        let loop_energy = GS.ose_full(loop_edge, loop_edge, scoped_mass.clone(), None);
        let q0 = GS.emr_mom(external_edge, GS.cind(0));
        let q_dot_k = function!(
            SPENSO_TAG.dot,
            GS.emr_vec(external_edge, minkowski.clone()),
            loop_momentum
        );
        let constant = &scoped_mass / loop_energy.pow(3);
        let (started, expected) = match dod {
            0 => (Atom::one() / loop_energy.pow(3), Atom::Zero),
            1 => (
                (&scoped_mass + &q0 + &q_dot_k / &loop_energy) / loop_energy.pow(3),
                constant.clone(),
            ),
            2 => {
                let expected = (&scoped_mass + &q0 + &q_dot_k / &loop_energy) / loop_energy.pow(2);
                (
                    &expected + q0.pow(2) / &scoped_mass / loop_energy.pow(2),
                    expected,
                )
            }
            _ => unreachable!(),
        };

        let actual = t_tilde(&ctx, &current, &root, &started, None, current.lmb())?;
        let difference = (actual.clone() - expected).expand();
        assert!(
            difference.is_zero(),
            "unexpected degree-{dod} soft projection: {difference}"
        );
        if dod == 2 {
            let ordinary = Local3DLoopRescaling::FullSubgraph.t(
                &ctx,
                &current,
                &root,
                &started,
                None,
                current.lmb(),
            )?;
            let overlap = Local3DLoopRescaling::FullSubgraph.t(
                &ctx,
                &current,
                &root,
                &actual,
                None,
                current.lmb(),
            )?;
            let soft_remainder = &started - &actual;
            let uv_remainder = Local3DLoopRescaling::FullSubgraph.t(
                &ctx,
                &current,
                &root,
                &soft_remainder,
                None,
                current.lmb(),
            )?;
            for (label, branch) in [
                ("S(X)", &actual),
                ("X-S(X)", &soft_remainder),
                ("U(S(X))", &overlap),
                ("U(X-S(X))", &uv_remainder),
            ] {
                let normalized = branch.clone().together().cancel().expand();
                assert!(
                    !normalized.is_zero(),
                    "the routed factorization oracle must have nonzero {label}"
                );
            }
            let expanded = ordinary + &actual - overlap;
            let factorized = &actual + uv_remainder;
            let difference = (factorized - expanded).together().cancel().expand();
            assert!(
                difference.is_zero(),
                "the routed nontrivial S(X)+U(X-S(X)) branch must equal U(X)+S(X)-U(S(X)): {difference}"
            );
        }
        assert!(!actual.contains_symbol(GS.rescale));
        assert!(!actual.contains_symbol(GS.m_uv_expansion));
        assert!(!actual.contains_symbol(GS.m_uv_vacuum));
        if dod > 0 {
            assert!(actual.contains_symbol(*LOCAL_3D_MASS_SCOPE));
        }
    }

    Ok::<_, eyre::Report>(())
}

#[test]
fn integrated_mass_vertex_keeps_its_soft_weight_and_scalar_logs() -> Result<()> {
    test_initialise()?;
    // Diagnostic effective-vertex fixture: the integrated tadpole supplies
    // the same mass/log structure as the finite massive-fermion self-energy.
    // Its enclosing scalar bubble isolates m/(k^2-m^2)^2, whose degree-one
    // soft completion is the complete physical-mass zero-external-momentum jet.
    let mut graph: Graph = dot!(
        digraph integrated_mass_vertex {
            edge [particle="H" num=1];
            node [num=1];
            ext [style=invis];
            ext -> A:0 [id=0];
            B:1 -> ext [id=1];
            A -> B [id=2 lmb_index=0];
            A -> B [id=3];
            A -> A [id=4 lmb_index=1];
        }
    )?;
    let current = two_point_node(&graph, 1);
    let mut child_subgraph = graph.empty_subgraph::<SuBitGraph>();
    child_subgraph.add(graph[&EdgeIndex(4)].1);
    let child = TestNode {
        lmb: graph.lmb_of(&child_subgraph),
        subgraph: child_subgraph,
        dod: 1,
        scheme: ApproximationType::IR,
    };
    let physical_mass = graph
        .iter_edges_of(child.subgraph())
        .next()
        .expect("the integrated source is the massive tadpole")
        .2
        .data
        .mass_atom();
    let coefficient = &physical_mass
        * (Atom::one() + Atom::num(2) * function!(Symbol::LOG, GS.mu_r_sq)
            - Atom::num(4) * function!(Symbol::LOG, physical_mass.clone()));
    let settings = UVgenerationSettings {
        add_marker: false,
        ..Default::default()
    };
    let orientation_pattern = OrientationPattern::default();
    let cutset = CutSet::empty(2);
    let options = graph.denominator_only_cff_3d_expression_options();
    let canonization = graph.get_esurface_canonization(&graph.loop_momentum_basis);
    let production = graph.generate_3d_expression_for_integrand(
        &[],
        &canonization,
        &options,
        Some(&Atom::one()),
    )?;
    let orientation =
        OrientationProjection::exact_expression(&production, &options, &orientation_pattern, false);
    let localizer = Localizer::new(&cutset, orientation);
    let unprojected = localizer.localize(&coefficient, &mut graph, &child)?;
    let active_subgraph = current.reduced_subgraph(&child);
    let expected = {
        let ctx = UVCtx::new(&graph, &settings);
        let (_, lmb) = coordinate_lmb(&ctx, &current, &child, None, &[], &active_subgraph)?;
        let started = DirectResidueBranches::from_transient(&unprojected.active)?.fallible_map(
            |_, atom| {
                start(
                    &ctx,
                    &current,
                    atom,
                    &Atom::one(),
                    Some(&active_subgraph),
                    &lmb,
                )
            },
        )?;
        // Express both physical propagators in the retained quotient frame
        // before taking the external jet: their edge momenta are not independent.
        let mut momentum_replacements =
            graph.uv_spatial_wrapped_replacement(&active_subgraph, &lmb, &[W_.x___]);
        momentum_replacements.extend(graph.uv_wrapped_replacement(
            &active_subgraph,
            &lmb,
            &[W_.x___],
        ));
        let mut expected = -started
            .factorized_sum()
            .replace_multiple(&momentum_replacements);
        // This independent degree-zero external jet leaves all scalar
        // logarithms literal and never acts on the child's frozen measure.
        for edge in &current.lmb.ext_edges {
            for momentum in [GS.emr_mom, GS.emr_vec] {
                expected = expected
                    .replace(function!(momentum, usize::from(*edge) as i64, W_.x___))
                    .with(Atom::Zero);
            }
        }
        expected
    };
    let actual = Direct3dApproximation::new(localizer, &mut graph, &settings).run_integrated(
        &[(&coefficient, &child)],
        &child,
        &current,
        &child,
        &current,
        &child,
    )?;
    let [sector] = actual.sectors()? else {
        panic!("one finite effective vertex must produce one active/frozen sector")
    };
    assert_eq!(sector.frozen_integrands, unprojected.frozen_integrands);
    let actual = sector.active.factorized_sum();
    assert!(
        !actual.is_zero(),
        "the physical-mass soft jet must be nonzero"
    );
    let residual = (actual - expected)
        .replace(function!(*LOCAL_3D_MASS_SCOPE, W_.a_))
        .with(Atom::one())
        .collect_factors()
        .expand_num()
        .together()
        .cancel();
    assert!(
        residual.is_zero(),
        "the finite mass vertex must retain its physical soft jet and unchanged mass/scale logarithms: {residual}"
    );
    Ok(())
}

#[test]
fn disconnected_finite_mass_vertex_keeps_its_sibling_mass_inert() -> Result<()> {
    test_initialise()?;
    // Diagnostic effective vertices: the first integrated tadpole lies inside
    // the bubble parent, while the second is outside it and has the same mass.
    // Giving the sibling a squared mass makes accidental co-scaling retain a
    // spurious quadratic external-soft term even after orientation summation.
    let mut graph: Graph = dot!(
        digraph disconnected_finite_mass_vertices {
            edge [particle="H" num=1];
            node [num=1];
            ext [style=invis];
            ext -> A:0 [id=0];
            C:1 -> ext [id=1];
            A -> B [id=2 lmb_index=0];
            A -> B [id=3];
            A -> A [id=4 lmb_index=1];
            C -> C [id=5 lmb_index=2];
            B -> C [id=6];
        }
    )?;
    let mut inside_subgraph = graph.empty_subgraph::<SuBitGraph>();
    inside_subgraph.add(graph[&EdgeIndex(4)].1);
    let inside = TestNode {
        lmb: graph.lmb_of(&inside_subgraph),
        subgraph: inside_subgraph,
        dod: 1,
        scheme: ApproximationType::IR,
    };
    let mut outside_subgraph = graph.empty_subgraph::<SuBitGraph>();
    outside_subgraph.add(graph[&EdgeIndex(5)].1);
    let outside = TestNode {
        lmb: graph.lmb_of(&outside_subgraph),
        subgraph: outside_subgraph,
        dod: 2,
        scheme: ApproximationType::IR,
    };
    let mut parent_subgraph = inside.subgraph.clone();
    for edge in [EdgeIndex(2), EdgeIndex(3)] {
        parent_subgraph.add(graph[&edge].1);
    }
    let current = TestNode {
        lmb: graph.lmb_of(&parent_subgraph),
        subgraph: parent_subgraph,
        dod: 1,
        scheme: ApproximationType::IR,
    };
    let prefix_subgraph = inside.subgraph.union(&outside.subgraph);
    let prefix = TestNode {
        lmb: graph.lmb_of(&prefix_subgraph),
        subgraph: prefix_subgraph,
        dod: 0,
        scheme: ApproximationType::MUV,
    };
    let mass = graph[&EdgeIndex(4)].0.mass_atom();
    assert_eq!(mass, graph[&EdgeIndex(5)].0.mass_atom());
    let sibling_coefficient = mass.pow(2);
    let one = Atom::one();
    let settings = UVgenerationSettings {
        add_marker: false,
        ..Default::default()
    };
    let orientation_pattern = OrientationPattern::default();
    let cutset = CutSet::empty(3);
    let options = graph.denominator_only_cff_3d_expression_options();
    let canonization = graph.get_esurface_canonization(&graph.loop_momentum_basis);
    // Match production: bridges carry fixed external momentum and belong to
    // the separately wrapped tree denominators, not the CFF loop source.
    let contract_edges = graph.paired_edges(&graph.tree_edges.subtract(&graph.initial_state_cut));
    let production = graph.generate_3d_expression_for_integrand(
        &contract_edges,
        &canonization,
        &options,
        Some(&one),
    )?;
    let orientation =
        OrientationProjection::exact_expression(&production, &options, &orientation_pattern, false);
    let localizer = Localizer::new(&cutset, orientation);
    let actual = Direct3dApproximation::new(localizer, &mut graph, &settings).run_integrated(
        &[(&mass, &inside), (&sibling_coefficient, &outside)],
        &prefix,
        &current,
        &inside,
        &current,
        &inside,
    )?;
    let control = Direct3dApproximation::new(localizer, &mut graph, &settings).run_integrated(
        &[(&mass, &inside), (&one, &outside)],
        &prefix,
        &current,
        &inside,
        &current,
        &inside,
    )?;
    let [actual] = actual.sectors()? else {
        panic!("one disconnected finite prefix must retain one active/frozen sector")
    };
    let [control] = control.sectors()? else {
        panic!("the scalar sibling control must retain the same sector")
    };
    assert_eq!(actual.frozen_integrands, control.frozen_integrands);
    let actual = actual.active.factorized_sum();
    assert!(!actual.is_zero());
    let sibling_scope = function!(*LOCAL_3D_MASS_SCOPE, usize::from(outside.lmb_id()) as i64);
    let expected = sibling_coefficient * sibling_scope.pow(2) * control.active.factorized_sum();
    let residual = (actual - expected).collect_factors().together().cancel();
    assert!(
        residual.is_zero(),
        "the disjoint finite vertex's equal physical mass must remain a scalar spectator: {residual}"
    );
    Ok(())
}

#[test]
fn nested_cff_toy_retains_four_signed_forest_families() {
    test_initialise().unwrap();
    let k1 = Atom::var(symbol!("nested_cff_toy_k1"));
    let p0 = Atom::var(symbol!("nested_cff_toy_p0"));
    let p1 = Atom::var(symbol!("nested_cff_toy_p1"));
    let q0 = Atom::var(symbol!("nested_cff_toy_q0"));
    let q1 = Atom::var(symbol!("nested_cff_toy_q1"));
    let mass = Atom::var(symbol!("nested_cff_toy_mass"));
    let m_uv = Atom::var(symbol!("nested_cff_toy_m_uv"));

    let energy =
        |momentum: &Atom, energy_mass: &Atom| (momentum.pow(2) + energy_mass.pow(2)).pow((1, 2));
    // One oriented two-propagator CFF residue after cancelling one
    // numerator energy against its residue measure.
    let cff_block = |loop_momentum: &Atom,
                     boundary_energy: &Atom,
                     boundary_spatial: &Atom,
                     block_mass: &Atom| {
        let loop_energy = energy(loop_momentum, block_mass);
        let shifted_energy = energy(&(loop_momentum - boundary_spatial), block_mass);
        Atom::one() / (&shifted_energy * (loop_energy + &shifted_energy + boundary_energy))
    };
    let physical_coefficients = |loop_momentum: &Atom,
                                 boundary_energy: &Atom,
                                 boundary_spatial: &Atom,
                                 block_mass: &Atom| {
        let loop_energy = energy(loop_momentum, block_mass);
        let shifted_linear = -(loop_momentum * boundary_spatial) / &loop_energy;
        let surface_constant = Atom::num(2) * &loop_energy;
        let surface_linear = boundary_energy + &shifted_linear;
        let denominator_constant = &loop_energy * &surface_constant;
        let denominator_linear = &loop_energy * surface_linear + shifted_linear * surface_constant;
        (denominator_constant, denominator_linear)
    };
    let ordinary_coefficients = |loop_momentum: &Atom,
                                 boundary_energy: &Atom,
                                 boundary_spatial: &Atom,
                                 block_mass: &Atom| {
        let loop_energy = energy(loop_momentum, &m_uv);
        let mass_gap = m_uv.pow(2) - block_mass.pow(2);
        let shifted_linear = -(loop_momentum * boundary_spatial) / &loop_energy;
        let shifted_quadratic = m_uv.pow(2) * boundary_spatial.pow(2)
            / (Atom::num(2) * loop_energy.pow(3))
            - &mass_gap / (Atom::num(2) * &loop_energy);
        let surface_constant = Atom::num(2) * &loop_energy;
        let surface_linear = boundary_energy + &shifted_linear;
        let surface_quadratic = m_uv.pow(2) * boundary_spatial.pow(2)
            / (Atom::num(2) * loop_energy.pow(3))
            - &mass_gap / &loop_energy;
        let denominator_constant = &loop_energy * &surface_constant;
        let denominator_linear =
            &loop_energy * &surface_linear + &shifted_linear * &surface_constant;
        let denominator_quadratic = &loop_energy * surface_quadratic
            + &shifted_linear * surface_linear
            + shifted_quadratic * surface_constant;
        (
            denominator_constant,
            denominator_linear,
            denominator_quadratic,
        )
    };
    let constant = |denominator_constant: &Atom| Atom::one() / denominator_constant;
    let linear = |denominator_constant: &Atom, denominator_linear: &Atom| {
        -denominator_linear / denominator_constant.pow(2)
    };
    let quadratic =
        |denominator_constant: &Atom, denominator_linear: &Atom, denominator_quadratic: &Atom| {
            denominator_linear.pow(2) / denominator_constant.pow(3)
                - denominator_quadratic / denominator_constant.pow(2)
        };

    let child_bare = cff_block(&k1, &p0, &p1, &mass);
    let child_bare_uv = cff_block(&k1, &p0, &p1, &m_uv);
    let (child_physical_constant, _) = physical_coefficients(&k1, &p0, &p1, &mass);
    let (child_uv_constant, child_uv_linear, _) = ordinary_coefficients(&k1, &p0, &p1, &mass);
    // H_1 keeps the physical constant jet and the ordinary linear jet.
    let completed_child =
        constant(&child_physical_constant) + linear(&child_uv_constant, &child_uv_linear);
    // A containing operation must re-expand this complete child.  Its
    // leading coefficient is not the bare child's UV completion.
    let completed_child_uv =
        constant(&child_uv_constant) + linear(&child_uv_constant, &child_uv_linear);

    let parent_bare = cff_block(&p1, &q0, &q1, &mass);
    let (parent_physical_constant, parent_physical_linear) =
        physical_coefficients(&p1, &q0, &q1, &mass);
    let (parent_uv_constant, parent_uv_linear, parent_uv_quadratic) =
        ordinary_coefficients(&p1, &q0, &q1, &mass);
    let physical_parent_jet = constant(&parent_physical_constant)
        + linear(&parent_physical_constant, &parent_physical_linear);
    let uv_parent_quadratic =
        quadratic(&parent_uv_constant, &parent_uv_linear, &parent_uv_quadratic);

    let f_empty = &child_bare * &parent_bare;
    let f_gamma = -(&completed_child * &parent_bare);
    let f_parent = -(&child_bare * &physical_parent_jet + &child_bare_uv * &uv_parent_quadratic);
    let f_gamma_gamma =
        &completed_child * &physical_parent_jet + &completed_child_uv * &uv_parent_quadratic;
    let regrouped = (&child_bare - &completed_child) * (&parent_bare - &physical_parent_jet)
        + (&completed_child_uv - &child_bare_uv) * &uv_parent_quadratic;
    let complete = &f_empty + &f_gamma + &f_parent + &f_gamma_gamma;

    // Reusing the bare child's leading coefficient in the nested family
    // is the precise stale-reexpansion error this oracle must reject.
    let stale_nested =
        &completed_child * &physical_parent_jet + &child_bare_uv * &uv_parent_quadratic;
    let stale_complete = &f_empty + &f_gamma + &f_parent + stale_nested;
    let fresh_child_uv = &completed_child_uv - &child_bare_uv;

    for (mass_label, mass_replacement) in [("generic mass", None), ("massless", Some(Atom::Zero))] {
        let specialize = |atom: &Atom| match &mass_replacement {
            Some(value) => atom.replace(mass.clone()).with(value.clone()),
            None => atom.clone(),
        };
        let f_empty = specialize(&f_empty);
        let f_gamma = specialize(&f_gamma);
        let f_parent = specialize(&f_parent);
        let f_gamma_gamma = specialize(&f_gamma_gamma);

        for (family_label, family) in [
            ("F_empty", &f_empty),
            ("F_gamma", &f_gamma),
            ("F_Gamma", &f_parent),
            ("F_gammaGamma", &f_gamma_gamma),
        ] {
            let normalized = family.clone().together().cancel().expand();
            assert!(
                !normalized.is_zero(),
                "{mass_label}: the direct-3D {family_label} family must be retained"
            );
        }

        let difference = (specialize(&complete) - specialize(&regrouped))
            .together()
            .cancel()
            .expand();
        assert!(
            difference.is_zero(),
            "{mass_label}: the four direct-3D forest families do not regroup: {difference}"
        );

        let fresh_child_uv = specialize(&fresh_child_uv).together().cancel().expand();
        assert!(
            !fresh_child_uv.is_zero(),
            "{mass_label}: the completed child must have a fresh outer-UV coefficient"
        );
        let stale_difference = (specialize(&stale_complete) - specialize(&regrouped))
            .together()
            .cancel()
            .expand();
        assert!(
            !stale_difference.is_zero(),
            "{mass_label}: reusing the bare child coefficient must fail the nested-family oracle"
        );
    }
}

#[test]
fn soft_hard_chart_coscales_free_and_terminal_uv_mass_from_inner_u() -> Result<()> {
    test_initialise().unwrap();
    let graph = scalar_two_point_graph();
    let root = root_node(&graph);
    let mut inner = two_point_node(&graph, 2);
    inner.scheme = ApproximationType::MUV;
    let mut outer = two_point_node(&graph, 2);
    outer.scheme = ApproximationType::IR;
    let settings = UVgenerationSettings::default();
    let ctx = UVCtx::new(&graph, &settings);

    let loop_edge = *inner
        .lmb
        .loop_edges
        .first()
        .expect("the scalar bubble has a loop carrier");
    let scope = function!(*LOCAL_3D_MASS_SCOPE, usize::from(loop_edge) as i64);
    let physical_mass = graph
        .iter_edges_of(&inner.subgraph)
        .find(|(pair, _, _)| pair.is_paired())
        .expect("the scalar bubble has an internal propagator")
        .2
        .data
        .mass_atom()
        * &scope;
    let uv_expansion_mass = Atom::var(GS.m_uv_expansion) * &scope;
    let uv_vacuum_mass = Atom::var(GS.m_uv_vacuum) * &scope;
    let physical_energy = GS.ose_full(loop_edge, loop_edge, physical_mass.clone(), None);

    // With one residue loop, U_2[1/E_m] retains the leading terminal-mUV
    // energy and the degree-two coefficient containing the free
    // (m^2-mUVexp^2) numerator factor.
    let inner_u = Local3DLoopRescaling::FullSubgraph.t(
        &ctx,
        &inner,
        &root,
        &(Atom::one() / physical_energy),
        None,
        inner.lmb(),
    )?;
    let terminal_energy = GS.ose_full(loop_edge, loop_edge, uv_vacuum_mass.clone(), None);
    let mass_difference = physical_mass.pow(2) - uv_expansion_mass.pow(2);
    let expected_inner =
        Atom::one() / &terminal_energy - &mass_difference / (Atom::num(2) * terminal_energy.pow(3));
    let difference = (inner_u.clone() - &expected_inner)
        .together()
        .cancel()
        .expand();
    assert!(
        difference.is_zero(),
        "the genuine inner U_2 projection lost its free/terminal UV-mass split: {difference}"
    );
    let without_terminal_argument = inner_u
        .replace(function!(GS.on_shell_energy, W_.a_, W_.prop_))
        .with(Atom::one());
    assert!(
        without_terminal_argument.contains_symbol(GS.m_uv_expansion),
        "the inner U oracle must retain a free mUVexp coefficient outside its terminal OSE"
    );
    assert!(
        !without_terminal_argument.contains_symbol(GS.m_uv_vacuum),
        "mUV must occur only in the terminal OSE"
    );
    assert!(inner_u.contains_symbol(GS.m_uv_vacuum));

    // Direct S holds the physical mass and both existing inner-U mass
    // roles fixed. Its dual hard chart must therefore co-scale every
    // scoped mass with the loop. The complete inner result is homogeneous
    // and S_1 reproduces it exactly.
    let outer_s = Local3DLoopRescaling::FullSubgraph.project(
        Local3DDeformation::Soft,
        &ctx,
        &outer,
        &root,
        &inner_u,
        None,
        outer.lmb(),
    )?;
    let difference = (outer_s.clone() - inner_u).together().cancel().expand();
    assert!(
        difference.is_zero(),
        "outer S did not co-scale the free and terminal UV masses: {difference}"
    );
    assert!(outer_s.contains_symbol(GS.m_uv_expansion));
    assert!(outer_s.contains_symbol(GS.m_uv_vacuum));
    assert!(outer_s.contains_symbol(*LOCAL_3D_MASS_SCOPE));
    assert!(!outer_s.contains_symbol(GS.rescale));

    Ok(())
}

#[test]
fn soft_laurent_projection_expands_materialized_internal_energy() -> Result<()> {
    test_initialise().unwrap();
    let graph = scalar_two_point_graph();
    let root = root_node(&graph);
    let settings = UVgenerationSettings::default();
    let ctx = UVCtx::new(&graph, &settings);
    let minkowski = Minkowski {}.new_rep(4).to_symbolic([]);

    for dod in 1..=2 {
        let current = two_point_node(&graph, dod);
        let loop_edge = *current
            .lmb
            .loop_edges
            .first()
            .expect("the scalar bubble has a loop carrier");
        let external_edge = *current
            .lmb
            .ext_edges
            .first()
            .expect("the scalar bubble has an external carrier");
        let physical_mass = graph
            .iter_edges_of(&current.subgraph)
            .find(|(pair, _, _)| pair.is_paired())
            .expect("the scalar bubble has an internal propagator")
            .2
            .data
            .mass_atom();
        let mass_scope = function!(*LOCAL_3D_MASS_SCOPE, usize::from(loop_edge) as i64);
        let mass = physical_mass * mass_scope;
        let loop_momentum = GS.emr_vec(loop_edge, minkowski.clone());
        let external_momentum = GS.emr_vec(external_edge, minkowski.clone());
        let loop_external = function!(SPENSO_TAG.dot, loop_momentum.clone(), external_momentum);
        let base = GS.ose_full(loop_edge, loop_edge, mass.clone(), None);
        let shifted_energy = GS
            .ose_full(loop_edge, loop_edge, mass.clone(), None)
            .replace(function!(GS.emr_mom, usize::from(loop_edge), W_.x___))
            .with(
                function!(GS.emr_mom, usize::from(loop_edge), W_.x___)
                    + function!(GS.emr_mom, usize::from(external_edge), W_.x___),
            );
        let started = if dod == 1 {
            Atom::one() / shifted_energy.pow(2)
        } else {
            Atom::one() / shifted_energy
        };

        let actual = t_tilde(&ctx, &current, &root, &started, None, current.lmb())?;
        let expected = if dod == 1 {
            Atom::one() / base.pow(2)
        } else {
            Atom::one() / &base + loop_external / base.pow(3)
        };
        let difference = (actual - expected).expand_dots()?.expand();
        assert!(
            difference.is_zero(),
            "unexpected degree-{dod} materialized internal-energy soft projection: {difference}"
        );
    }

    Ok(())
}

#[test]
fn soft_projection_holds_cograph_tree_denominators() -> Result<()> {
    test_initialise().unwrap();
    let graph = scalar_two_point_graph();
    let root = root_node(&graph);
    let current = two_point_node(&graph, 1);
    let settings = UVgenerationSettings::default();
    let ctx = UVCtx::new(&graph, &settings);
    let external_edge = *current
        .lmb
        .ext_edges
        .first()
        .expect("the scalar bubble has an external carrier");
    let external_energy = GS.emr_mom(external_edge, GS.cind(0));
    let mass = graph
        .iter_edges_of(&current.subgraph)
        .find(|(pair, _, _)| pair.is_paired())
        .expect("the scalar bubble has an internal propagator")
        .2
        .data
        .mass_atom();
    let loop_edge = *current
        .lmb
        .loop_edges
        .first()
        .expect("the scalar bubble has a loop carrier");
    let mass_scope = function!(*LOCAL_3D_MASS_SCOPE, usize::from(loop_edge) as i64);
    let mass = mass * mass_scope;
    let loop_energy = GS.ose_full(loop_edge, loop_edge, mass.clone(), None);
    let tree_denominator = GS.wrap_tree_denoms(Atom::one() / external_energy.clone().pow(2));
    let started = &tree_denominator * (&mass + external_energy) / loop_energy.pow(3);

    let actual = t_tilde(&ctx, &current, &root, &started, None, current.lmb())?;
    let difference = (actual - mass * tree_denominator / loop_energy.pow(3)).expand();
    assert!(
        difference.is_zero(),
        "a component soft jet must leave the co-graph bridge denominator literal: {difference}"
    );

    Ok(())
}

#[test]
fn nested_soft_projection_scales_only_the_component_boundary_after_routing() -> Result<()> {
    test_initialise().unwrap();
    let graph: Graph =
        include_str!("../../../../../../tests/resources/graphs/gamma_star_ddbar_top_bubble.dot")
            .into_graph(&crate::utils::load_generic_model("sm"))?;
    let mut child_subgraph = graph.empty_subgraph::<SuBitGraph>();
    for edge in [EdgeIndex(7), EdgeIndex(8)] {
        child_subgraph.add(graph[&edge].1);
    }
    let child = TestNode {
        lmb: graph.try_compatible_sub_lmb(
            &child_subgraph,
            graph.dummy_less_full_crown(&child_subgraph),
            &graph.loop_momentum_basis,
        )?,
        subgraph: child_subgraph,
        dod: 1,
        scheme: ApproximationType::IR,
    };
    let mut parent_subgraph = graph.empty_subgraph::<SuBitGraph>();
    for edge in [3, 4, 5, 6, 7, 8].map(EdgeIndex) {
        parent_subgraph.add(graph[&edge].1);
    }
    let parent_lmb = graph.try_compatible_sub_lmb(
        &parent_subgraph,
        graph.dummy_less_full_crown(&parent_subgraph),
        &graph.loop_momentum_basis,
    )?;
    let root = root_node(&graph);
    let settings = UVgenerationSettings::default();
    let ctx = UVCtx::new(&graph, &settings);

    assert_eq!(
        graph.loop_momentum_basis.loop_edges.raw,
        vec![EdgeIndex(6), EdgeIndex(8)]
    );
    assert_eq!(child.lmb.loop_edges.raw, vec![EdgeIndex(8)]);
    assert_eq!(parent_lmb.loop_edges.raw, vec![EdgeIndex(6), EdgeIndex(8)]);
    assert!(
        child
            .lmb
            .loop_edges
            .iter()
            .all(|edge| parent_lmb.loop_edges.contains(edge)),
        "a nested child's canonical loop carrier must remain a loop carrier of its containing parent"
    );
    assert!(child.lmb.ext_edges.contains(&EdgeIndex(0)));
    assert!(child.lmb.ext_edges.contains(&EdgeIndex(5)));

    let loop_edge = *child
        .lmb
        .loop_edges
        .first()
        .expect("the nested child has a loop carrier");
    let mass = graph[EdgeIndex(7)].mass_atom()
        * function!(*LOCAL_3D_MASS_SCOPE, usize::from(loop_edge) as i64);
    let graph_external_spectator = GS.emr_mom(EdgeIndex(0), GS.cind(0));
    let boundary_energy_shift = GS.emr_mom(EdgeIndex(5), GS.cind(0));
    let retained_parent_pole = GS.ose(EdgeIndex(5));
    let shifted_energy = GS.ose_full(EdgeIndex(7), loop_edge, mass.clone(), None);
    let base_energy = GS.ose_full(loop_edge, loop_edge, mass.clone(), None);
    let spectator = &graph_external_spectator / &retained_parent_pole;
    let started = &spectator * (&mass + boundary_energy_shift) / shifted_energy.pow(3);

    let actual = t_tilde(&ctx, &child, &root, &started, None, child.lmb())?;
    let expected = spectator * mass / base_energy.pow(3);
    let difference = (actual - expected).expand();
    assert!(
        difference.is_zero(),
        "the child soft jet must scale its routed e5 boundary, but neither the graph-external e0 spectator nor the retained parent OSE(e5): {difference}",
    );

    Ok(())
}

#[test]
fn appendix_b1_nested_routes_keep_the_graph_canonical_enclosing_chart() -> Result<()> {
    test_initialise()?;
    let graph: Graph = include_str!(
        "../../../../../../tests/resources/graphs/paper_appendix_b1_nested_gluon_self_energy.dot"
    )
    .into_graph(&crate::utils::load_generic_model("sm"))?;
    let mut parent_subgraph = graph.empty_subgraph::<SuBitGraph>();
    for edge in [2, 3, 4, 5, 6].map(EdgeIndex) {
        parent_subgraph.add(graph[&edge].1);
    }
    let parent = TestNode {
        lmb: graph.try_compatible_sub_lmb(
            &parent_subgraph,
            graph.dummy_less_full_crown(&parent_subgraph),
            &graph.loop_momentum_basis,
        )?,
        subgraph: parent_subgraph,
        dod: 2,
        scheme: ApproximationType::IR,
    };
    let settings = UVgenerationSettings::default();
    let ctx = UVCtx::new(&graph, &settings);

    for (child_edges, expected_child_route) in [
        (&[3, 4][..], vec![EdgeIndex(3)]),
        (&[2, 3, 5, 6][..], vec![EdgeIndex(2)]),
    ] {
        let mut child_subgraph = graph.empty_subgraph::<SuBitGraph>();
        for edge in child_edges.iter().copied().map(EdgeIndex) {
            child_subgraph.add(graph[&edge].1);
        }
        let child = TestNode {
            lmb: graph.try_compatible_sub_lmb(
                &child_subgraph,
                graph.dummy_less_full_crown(&child_subgraph),
                &graph.loop_momentum_basis,
            )?,
            subgraph: child_subgraph,
            dod: 1,
            scheme: ApproximationType::MUV,
        };
        assert_eq!(child.lmb.loop_edges.raw, expected_child_route);
        let prior_frame = DirectCoordinateFrame {
            active_subgraph: child.subgraph.clone(),
            lmb: child.lmb.clone(),
        };
        let (canonical, selected) = coordinate_lmb(
            &ctx,
            &parent,
            &child,
            Some(&child.subgraph),
            &[prior_frame],
            &parent.subgraph,
        )?;
        assert_eq!(canonical.loop_edges.raw, vec![EdgeIndex(2), EdgeIndex(3)]);
        assert_eq!(
            selected.loop_edges.raw, canonical.loop_edges.raw,
            "both Appendix-B.1 child histories must enter the outer Taylor operator in its graph-canonical chart"
        );
    }

    Ok(())
}

#[test]
fn nested_banana_retains_quotient_carriers_after_an_integrated_prefix() -> Result<()> {
    test_initialise()?;
    let graph: Graph = dot!(
        digraph banana {
            edge [particle=scalar_1 num=1]
            node [num=1]
            ext [style=invis]
            ext -> a:0 [id=5]
            b:1 -> ext [id=4]
            a -> b [id=1]
            a -> b [id=2]
            a -> b [id=3]
            a -> b [id=0]
        },
        "scalars"
    )?;
    let subgraph = |edges: &[usize]| {
        let mut subgraph = graph.empty_subgraph::<SuBitGraph>();
        for edge in edges.iter().copied().map(EdgeIndex) {
            subgraph.add(graph[&edge].1);
        }
        subgraph
    };
    let node = |edges: &[usize]| {
        let subgraph = subgraph(edges);
        TestNode {
            lmb: graph.lmb_of(&subgraph),
            subgraph,
            dod: 0,
            scheme: ApproximationType::MUV,
        }
    };
    let integrated = node(&[1, 2]);
    let child = node(&[0, 1, 2]);
    let parent = node(&[0, 1, 2, 3]);
    let settings = UVgenerationSettings::default();
    let ctx = UVCtx::new(&graph, &settings);
    let child_active = child.reduced_subgraph(&integrated);
    let (_, child_lmb) = coordinate_lmb(&ctx, &child, &integrated, None, &[], &child_active)?;
    assert_eq!(child_lmb.loop_edges.raw, vec![EdgeIndex(0)]);
    assert!(
        graph.loop_momentum_basis.edge_signatures[EdgeIndex(0)]
            .external
            .iter()
            .any(|sign| sign.is_sign()),
        "the original graph routes Q0 with a fixed external shift"
    );
    let parent_active = child_active.union(&parent.reduced_subgraph(&child));
    let prior_frame = DirectCoordinateFrame {
        active_subgraph: child_active.clone(),
        lmb: child_lmb,
    };
    let (_, selected) = coordinate_lmb(
        &ctx,
        &parent,
        &child,
        Some(&child_active),
        &[prior_frame],
        &parent_active,
    )?;
    // The original graph carrier Q3 keeps its first position; the retained
    // quotient carrier Q0 follows it in the compatible enclosing chart.
    assert_eq!(selected.loop_edges.raw, vec![EdgeIndex(3), EdgeIndex(0)]);
    // Contracting the integrated two-line bubble makes each remaining line
    // a self-loop. Certify its signed momentum, rather than only its loop id.
    for (index, edge) in selected.loop_edges.iter().enumerate() {
        let signature = &selected.edge_signatures[*edge];
        assert!(signature.external.iter().all(|sign| sign.is_zero()));
        for (loop_index, sign) in signature.internal.iter().enumerate() {
            assert_eq!(sign.is_sign(), loop_index == index);
            if loop_index == index {
                assert!(!sign.is_negative());
            }
        }
    }
    let q0 = GS.emr_mom(EdgeIndex(0), GS.cind(1));
    let q3 = GS.emr_mom(EdgeIndex(3), GS.cind(1));
    let fixed_external = GS.emr_mom(EdgeIndex(5), GS.cind(1));
    let input = q0 * q3 * fixed_external;
    let rescaled = Direct3dApproximation::t_rescale(
        Local3DDeformation::Ordinary,
        &ctx,
        &parent,
        &child,
        &input,
        Some(&parent_active),
        &selected,
        false,
    )?;
    assert!(
        (rescaled - input * Atom::var(GS.rescale).pow(-2))
            .expand()
            .is_zero(),
        "both quotient carriers must scale homogeneously while graph-external momentum stays fixed"
    );
    Ok(())
}

#[test]
fn nested_route_rejects_a_retained_affine_graph_external_carrier() -> Result<()> {
    test_initialise()?;
    let graph: Graph = dot!(
        digraph nested_soft_affine_routing {
            edge [particle=scalar_1 num=1]
            node [num=1]
            ext [style=invis]
            ext -> v1:0 [id=0]
            v4:1 -> ext [id=1]
            v1 -> v2 [id=2]
            v2 -> v3 [id=3 lmb_id=0]
            v2 -> v3 [id=4]
            v3 -> v4 [id=5]
            v1 -> v4 [id=6 lmb_id=1]
            v2 -> v2 [id=7 lmb_id=2]
        },
        "scalars"
    )?;
    let mut child_subgraph = graph.empty_subgraph::<SuBitGraph>();
    for edge in [2, 3, 5, 6].map(EdgeIndex) {
        child_subgraph.add(graph[&edge].1);
    }
    let affine_carrier = EdgeIndex(2);
    assert!(
        graph.loop_momentum_basis.edge_signatures[affine_carrier]
            .external
            .iter()
            .any(|sign| sign.is_sign()),
        "the retained test carrier must require a graph-external affine shift"
    );
    let affine_lmb = graph
        .generate_loop_momentum_bases_of(&child_subgraph)
        .into_iter()
        .find(|candidate| candidate.loop_edges.contains(&affine_carrier))
        .expect("the child cycle must admit its affine edge as a nominal loop carrier");
    let child = TestNode {
        lmb: affine_lmb.clone(),
        subgraph: child_subgraph,
        dod: 0,
        scheme: ApproximationType::MUV,
    };
    let external: SuBitGraph = graph.external_filter();
    let parent_subgraph = graph.full_filter().subtract(&external);
    let parent = TestNode {
        lmb: graph.lmb_of(&parent_subgraph),
        subgraph: parent_subgraph,
        dod: 2,
        scheme: ApproximationType::IR,
    };
    let prior_frame = DirectCoordinateFrame {
        active_subgraph: child.subgraph.clone(),
        lmb: affine_lmb,
    };
    let settings = UVgenerationSettings::default();
    let ctx = UVCtx::new(&graph, &settings);
    // Integrating an unrelated tadpole does not change this carrier's affine
    // routing. A nonempty inactive prefix must not switch off the guard.
    let mut tadpole = graph.empty_subgraph::<SuBitGraph>();
    tadpole.add(graph[&EdgeIndex(7)].1);
    for inactive in [graph.empty_subgraph(), tadpole] {
        let given = TestNode {
            subgraph: child.subgraph.union(&inactive),
            lmb: child.lmb.clone(),
            dod: child.dod,
            scheme: child.scheme,
        };
        let active = parent.subgraph.subtract(&inactive);
        let error = coordinate_lmb(
            &ctx,
            &parent,
            &given,
            Some(&child.subgraph),
            std::slice::from_ref(&prior_frame),
            &active,
        )
        .expect_err("an affine retained carrier must not become a fixed parent coordinate");
        let message = error.to_string();
        assert!(
            message.contains("affine graph-external momentum shift"),
            "unexpected routing error: {message}"
        );
        assert!(message.contains(&parent.subgraph.string_label()));
    }

    Ok(())
}

#[test]
fn production_ir_dispatch_composes_h_for_d0_d1_d2() -> Result<()> {
    test_initialise().unwrap();
    let mut graph = scalar_two_point_graph();
    let settings = UVgenerationSettings {
        generate_integrated: false,
        ..Default::default()
    };
    let options = graph.denominator_only_cff_3d_expression_options();
    let canonization = graph.get_esurface_canonization(&graph.loop_momentum_basis);
    let production = graph.generate_3d_expression_for_integrand(
        &[],
        &canonization,
        &options,
        Some(&Atom::one()),
    )?;
    let orientation_pattern = OrientationPattern::default();
    let orientation =
        OrientationProjection::exact_expression(&production, &options, &orientation_pattern, false);
    let key = DirectResidueKey::production(OrientationID(0));
    let ctx = UVCtx::new(&graph, &settings);

    for dod in 0..=2 {
        let ir = two_point_node(&graph, dod);
        let ordinary = TestNode {
            subgraph: ir.subgraph.clone(),
            lmb: ir.lmb.clone(),
            dod,
            scheme: ApproximationType::MUV,
        };
        let loop_edge = *ir
            .lmb
            .loop_edges
            .first()
            .expect("the scalar bubble has a loop carrier");
        let scope = function!(*LOCAL_3D_MASS_SCOPE, usize::from(loop_edge) as i64);
        let physical_mass = graph
            .iter_edges_of(&ir.subgraph)
            .find(|(pair, _, _)| pair.is_paired())
            .expect("the scalar bubble has an internal propagator")
            .2
            .data
            .mass_atom()
            * &scope;
        let physical_energy = GS.ose_full(loop_edge, loop_edge, physical_mass.clone(), None);
        let integrand = Atom::one() / physical_energy.pow(3 - dod as i64);

        // Taking `given=current` makes the supplied atom the completed
        // component expression, so this invokes the production dispatch
        // without multiplying another graph numerator into the oracle.
        let input =
            Integrands::from_iter([(crate::cff::CutCFFIndex::new_all_none(), integrand.clone())]);
        let input = DirectResidueBranches::from_keyed([(key.clone(), input)])?;
        // Resolve the retained definitions only for this scalar oracle.
        let ir_result = apply_taylor(&ctx, orientation, &ir, &ir, None, ir.lmb(), &input)?
            .materialize(false)?
            .resolved()?;
        let ir_step = Local3DLoopRescaling::FullSubgraph.projection_step(
            &ctx,
            &ir,
            &ir,
            ir.subgraph(),
            None,
            ir.lmb(),
            ir.lmb(),
        );
        if dod == 0 {
            assert_eq!(ir_step.materialization, Local3DMaterialization::DirectU);
            assert_eq!(
                ir_step.conceptual_branches,
                vec![Local3DSignedBranch {
                    branch: Local3DConceptualBranch::U,
                    coefficient: 1,
                }]
            );
        } else {
            assert_eq!(
                ir_step.materialization,
                Local3DMaterialization::FactorizedSoft
            );
            assert_eq!(
                ir_step.conceptual_branches,
                vec![
                    Local3DSignedBranch {
                        branch: Local3DConceptualBranch::U,
                        coefficient: 1,
                    },
                    Local3DSignedBranch {
                        branch: Local3DConceptualBranch::S,
                        coefficient: 1,
                    },
                    Local3DSignedBranch {
                        branch: Local3DConceptualBranch::US,
                        coefficient: -1,
                    },
                ]
            );
        }
        let ir_result = ir_result
            .iter()
            .next()
            .expect("the projection retains the oracle residue")
            .1
            .clone();
        let ordinary_result = apply_taylor(
            &ctx,
            orientation,
            &ordinary,
            &ordinary,
            None,
            ordinary.lmb(),
            &input,
        )?
        .materialize(false)?
        .resolved()?;
        let ordinary_result = ordinary_result
            .iter()
            .next()
            .expect("the projection retains the oracle residue")
            .1
            .clone();

        if dod > 0 {
            let started = start(
                &ctx,
                &ir,
                &integrand,
                &graph
                    .numerator(&ir.reduced_subgraph(&ir), ir.subgraph())
                    .get_single_atom()
                    .unwrap(),
                None,
                ir.lmb(),
            )?;
            let ordinary_branch =
                Local3DLoopRescaling::FullSubgraph.t(&ctx, &ir, &ir, &started, None, ir.lmb())?;
            let soft_branch = Local3DLoopRescaling::FullSubgraph.project(
                Local3DDeformation::Soft,
                &ctx,
                &ir,
                &ir,
                &started,
                None,
                ir.lmb(),
            )?;
            let normalized_soft = soft_branch.clone().together().cancel().expand();
            assert!(
                !normalized_soft.is_zero(),
                "the degree-{dod} factorization oracle must exercise a nonzero soft branch"
            );
            let overlap_branch = Local3DLoopRescaling::FullSubgraph.t(
                &ctx,
                &ir,
                &ir,
                &soft_branch,
                None,
                ir.lmb(),
            )?;
            let uv_remainder = Local3DLoopRescaling::FullSubgraph.t(
                &ctx,
                &ir,
                &ir,
                &(&started - &soft_branch),
                None,
                ir.lmb(),
            )?;
            let factorized = &soft_branch + uv_remainder;
            let expanded = ordinary_branch + soft_branch - overlap_branch;
            let difference = (factorized - expanded).together().cancel().expand();
            assert!(
                difference.is_zero(),
                "S(X)+U(X-S(X)) must equal U(X)+S(X)-U(S(X)) for nonzero degree-{dod} S: {difference}"
            );
        }

        assert!(!ir_result.contains_symbol(GS.rescale));
        if dod == 0 {
            assert_eq!(
                (ir_result - ordinary_result).expand(),
                Atom::Zero,
                "H_0 must be the ordinary U_0 projection"
            );
        } else {
            assert_eq!(
                (ir_result - &integrand).expand(),
                Atom::Zero,
                "S_{} contains the complete momentum-independent physical-mass term, so H_{dod} must reproduce it",
                dod - 1,
            );
            assert!(
                !(ordinary_result - &integrand).expand().is_zero(),
                "the degree-{dod} oracle must distinguish H from ordinary U"
            );
        }
    }

    Ok(())
}

#[test]
fn one_orientation_ir_dispatch_stores_one_combined_h_atom_per_residue() -> Result<()> {
    test_initialise().unwrap();
    let mut graph = scalar_two_point_graph();
    let current = two_point_node(&graph, 2);
    let root = root_node(&graph);
    let orientation_pattern = OrientationPattern::default();
    let cutset = CutSet::empty(1);

    let options = graph.denominator_only_cff_3d_expression_options();
    let canonization = graph.get_esurface_canonization(&graph.loop_momentum_basis);
    let mut production = graph.generate_3d_expression_for_integrand(
        &[],
        &canonization,
        &options,
        Some(&Atom::one()),
    )?;
    production.expression.orientations.truncate(1);
    let orientation =
        OrientationProjection::exact_expression(&production, &options, &orientation_pattern, false);
    let localizer = Localizer::new(&cutset, orientation);
    let input = Direct3dCts::root(&graph, localizer)?;
    let branches = input.branches()?;
    let mut selected_terms = branches.iter_keys();
    let (_residue_key, residue_integrands) = selected_terms
        .next()
        .expect("the selected scalar-bubble orientation has a CFF term");
    assert!(
        selected_terms.next().is_none(),
        "the admitted-orientation CFF must contain exactly one orientation term"
    );
    let oriented_cff = residue_integrands.iter().next().unwrap().1;
    let input_key = *residue_integrands.iter().next().unwrap().0;
    let input_atom = oriented_cff.clone();
    assert_eq!(
        residue_integrands.iter().count(),
        1,
        "the empty-cut one-orientation fixture must have exactly one residue"
    );
    let settings = UVgenerationSettings {
        generate_integrated: false,
        add_marker: false,
        ..Default::default()
    };
    let active_subgraph = current.reduced_subgraph(&root);
    let expected = {
        let ctx = UVCtx::new(&graph, &settings);
        let started = start(
            &ctx,
            &current,
            &input_atom,
            &graph
                .numerator(&current.reduced_subgraph(&root), root.subgraph())
                .get_single_atom()
                .unwrap(),
            Some(&active_subgraph),
            current.lmb(),
        )?;
        let ordinary = Local3DLoopRescaling::FullSubgraph.t(
            &ctx,
            &current,
            &root,
            &started,
            Some(&active_subgraph),
            current.lmb(),
        )?;
        let soft = Local3DLoopRescaling::FullSubgraph.project(
            Local3DDeformation::Soft,
            &ctx,
            &current,
            &root,
            &started,
            Some(&active_subgraph),
            current.lmb(),
        )?;
        let overlap = Local3DLoopRescaling::FullSubgraph.t(
            &ctx,
            &current,
            &root,
            &soft,
            Some(&active_subgraph),
            current.lmb(),
        )?;
        // This scalar one-orientation CFF has no strictly negative soft
        // coefficient after dressing, so S and US both vanish.  It is
        // nevertheless a useful dispatch oracle: production must store
        // the single completed H=U+S-US atom, not three branch atoms.
        assert!(soft.is_zero());
        assert!(overlap.is_zero());
        let expanded = &ordinary + &soft - &overlap;
        let factorized = &soft
            + Local3DLoopRescaling::FullSubgraph.t(
                &ctx,
                &current,
                &root,
                &(&started - &soft),
                Some(&active_subgraph),
                current.lmb(),
            )?;
        let difference = (&factorized - &expanded).together().cancel().expand();
        assert!(
            difference.is_zero(),
            "S(X)+U(X-S(X)) must retain the exact U(X)+S(X)-U(S(X)) operator: {difference}"
        );
        -expanded
    };

    let actual = Direct3dApproximation::new(localizer, &mut graph, &settings)
        .run_local(&input, &current, &root, &current, &root)?;
    let sectors = actual
        .sectors()
        .expect("a local operation retains its active component sector");
    assert_eq!(
        sectors.len(),
        1,
        "U, S, and US must not be exported as separate active sectors"
    );
    assert_eq!(
        sectors[0]
            .active
            .iter_keys()
            .flat_map(|(_, integrands)| integrands.iter())
            .count(),
        1,
        "the selected residue/orientation must contain one combined sector atom"
    );
    let paths = actual.projection_paths();
    let [path] = paths.as_slice() else {
        panic!("one connected local operation must retain one projection path")
    };
    let [step] = path.steps.as_slice() else {
        panic!("one connected local operation must retain one projection step")
    };
    assert_eq!(
        step.canonical_route_loop_edges,
        current
            .lmb
            .loop_edges
            .iter()
            .map(|edge| usize::from(*edge))
            .collect::<Vec<_>>()
    );
    assert_eq!(step.route_loop_edges, vec![usize::from(current.lmb_id())]);
    assert_eq!(step.materialization, Local3DMaterialization::FactorizedSoft);
    assert_eq!(
        step.conceptual_branches,
        vec![
            Local3DSignedBranch {
                branch: Local3DConceptualBranch::U,
                coefficient: 1,
            },
            Local3DSignedBranch {
                branch: Local3DConceptualBranch::S,
                coefficient: 1,
            },
            Local3DSignedBranch {
                branch: Local3DConceptualBranch::US,
                coefficient: -1,
            },
        ]
    );

    let materialized = actual.branches()?.materialize(false)?;
    let mut actual_entries = materialized.iter();
    let (actual_key, actual_atom) = actual_entries
        .next()
        .expect("the completed H counterterm must retain its residue atom");
    assert_eq!(*actual_key, input_key);
    assert!(
        actual_entries.next().is_none(),
        "the completed H counterterm must export exactly one atom for the selected residue/orientation"
    );
    let difference = (actual_atom - expected).together().cancel().expand();
    assert!(
        difference.is_zero(),
        "the sole exported atom must be the signed U+S-US composite: {difference}"
    );

    Ok(())
}

#[test]
#[should_panic(
    expected = "local on-shell counterterms are deferred until local counterterms can be derived from the 4D expanded representation"
)]
fn local_3d_os_dispatch_is_unconditionally_deferred() {
    test_initialise().unwrap();
    let mut graph = scalar_two_point_graph();
    let mut current = two_point_node(&graph, 1);
    current.scheme = ApproximationType::OS;
    let given = root_node(&graph);
    let orientation_pattern = OrientationPattern::default();
    let cutset = CutSet::empty(1);
    let localizer = Localizer::new(&cutset, OrientationProjection::four_d(&orientation_pattern));
    let settings = UVgenerationSettings {
        generate_integrated: false,
        ..Default::default()
    };

    Direct3dApproximation::new(localizer, &mut graph, &settings)
        .run_local(
            &Direct3dCts::Root(
                DirectResidueBranches::production(OrientationID(0), Integrands::root()).unwrap(),
            ),
            &current,
            &given,
            &current,
            &given,
        )
        .unwrap();
}
