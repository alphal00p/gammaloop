use super::*;

use color_eyre::Result;

use crate::graph::cuts::CutSet;
use crate::graph::parse::IntoGraph;
use crate::initialisation::test_initialise;
use crate::settings::global::GenerationSettings;
use crate::settings::global::MediumMode;
use crate::utils::load_generic_model;
use crate::uv::UVOrchestrator;
use crate::uv::approx::local_3d::{
    Local3DConceptualBranch, Local3DLoopRescaling, Local3DMaterialization, Local3DSignedBranch,
};
use crate::uv::approx::{CutStructure, ForestNodeLike, UVCtx};
use crate::uv::hedge_poset::{ForestNode, Wood};
use crate::uv::orchestrator::validate_local_counterterm_schemes;
use crate::uv::test_fixtures::{Ring, edge_set, medium, ring_uv_settings, thermal_medium};

impl Ring {
    /// The expected hosts of the host policy, keyed by their edges. The kernel tests read
    /// them on ring topologies with the same edges and test numerators.
    pub(in crate::uv) fn expected_regions(
        self,
        graph: &Graph,
        vacuum_subtraction: bool,
    ) -> BTreeMap<BTreeSet<EdgeIndex>, ThermalZeroModeRegion> {
        let union = [2, 3]
            .iter()
            .chain(self.second_loop_edges())
            .copied()
            .collect_vec();
        let mut ring_tuples = vec![tuple(&[4], ZeroModeSector::Complete, 1)];
        if vacuum_subtraction {
            ring_tuples.push(tuple(&[3], ZeroModeSector::HardVacuum, 0));
        }
        let inner = matches!(self, Self::Mixed | Self::Ghost);
        let union_region = ThermalZeroModeRegion {
            kind: ZeroModeHostKind::Union,
            sigma: edge_set(&[0, 1]),
            classes: vec![class(&[0, 1])],
            omega: -1,
            soft_loops: 1,
            tuples: ring_tuples,
            zero_bins: if inner {
                vec![graph.edge_subgraph(&[2, 3])]
            } else {
                Vec::new()
            },
            coordinate_loop_edges: edge_list(&[2, 4]),
        };
        let mut regions = BTreeMap::from([(edge_set(&union), union_region)]);
        if inner {
            regions.insert(
                edge_set(&[2, 3]),
                ThermalZeroModeRegion {
                    kind: ZeroModeHostKind::Connected,
                    sigma: edge_set(&[0, 1, 4, 5]),
                    // e5 = k - q is a static condition implied by q and k.
                    classes: vec![class(&[0, 1]), class(&[4])],
                    omega: 0,
                    soft_loops: 2,
                    tuples: vec![tuple(&[4, 2], ZeroModeSector::Complete, 0)],
                    zero_bins: Vec::new(),
                    coordinate_loop_edges: edge_list(&[2]),
                },
            );
        }
        regions
    }
}

fn edge_list(edges: &[usize]) -> Vec<EdgeIndex> {
    edges.iter().copied().map(EdgeIndex).collect()
}

fn class(edges: &[usize]) -> ZeroModeClassSpec {
    ZeroModeClassSpec::new(edge_set(edges)).unwrap()
}

fn tuple(orders: &[u32], sector: ZeroModeSector, taylor_order: i64) -> ZeroModeTuple {
    ZeroModeTuple {
        orders: orders.to_vec(),
        sector,
        taylor_order,
    }
}

/// `region` with only its static tuples, each through the relative Taylor
/// order `taylor_order` instead of the policy's: a host that keeps a single
/// Taylor term of its static mode, whose image under an enclosing U the
/// kernel's depth pins read on its own.
pub(in crate::uv) fn static_region_through(
    region: &ThermalZeroModeRegion,
    taylor_order: i64,
) -> ThermalZeroModeRegion {
    ThermalZeroModeRegion {
        tuples: region
            .tuples
            .iter()
            .filter(|tuple| tuple.sector == ZeroModeSector::Complete)
            .map(|tuple| ZeroModeTuple {
                taylor_order,
                ..tuple.clone()
            })
            .collect(),
        ..region.clone()
    }
}

/// A synthetic connected host region with the soft set `sigma`, which the
/// wood reads; it requests no orders.
pub(in crate::uv) fn connected_region(sigma: &[usize]) -> ThermalZeroModeRegion {
    ThermalZeroModeRegion {
        kind: ZeroModeHostKind::Connected,
        sigma: edge_set(sigma),
        classes: Vec::new(),
        omega: 0,
        soft_loops: 1,
        tuples: Vec::new(),
        zero_bins: Vec::new(),
        coordinate_loop_edges: Vec::new(),
    }
}

fn uncut(graph: &Graph) -> SuBitGraph {
    graph.full_filter().subtract(&graph.initial_state_cut.left)
}

fn classified(graph: &Graph, settings: &UVgenerationSettings) -> Result<Vec<Spinney>> {
    graph.classified_spinneys(&uncut(graph), settings, &graph.loop_momentum_basis)
}

/// The spinneys of the uncut graph after the zero-mode policy and the
/// full-observable restriction.
pub(in crate::uv) fn zero_mode_policy(
    graph: &Graph,
    settings: &UVgenerationSettings,
    medium: &MediumSettings,
) -> Result<Vec<Spinney>, ZeroModeHostError> {
    let mut spinneys = classified(graph, settings).unwrap();
    assign_zero_mode_hosts(
        graph,
        &uncut(graph),
        &mut spinneys,
        settings,
        medium,
        &CutStructure::empty(graph),
    )?;
    graph.restrict_full_observable_spinney(&mut spinneys, medium)?;
    Ok(spinneys)
}

fn edges_of(graph: &Graph, spinney: &Spinney) -> BTreeSet<EdgeIndex> {
    graph
        .underlying
        .paired_edges(spinney.filter())
        .into_iter()
        .collect()
}

fn regions(
    graph: &Graph,
    spinneys: &[Spinney],
) -> BTreeMap<BTreeSet<EdgeIndex>, ThermalZeroModeRegion> {
    spinneys
        .iter()
        .filter_map(|spinney| {
            spinney
                .thermal_zero_mode()
                .map(|region| (edges_of(graph, spinney), region.clone()))
        })
        .collect()
}

/// The spinneys as (edges, scheme, hosts a region), sorted.
fn classification(
    graph: &Graph,
    spinneys: &[Spinney],
) -> Vec<(BTreeSet<EdgeIndex>, ApproximationType, bool)> {
    spinneys
        .iter()
        .map(|spinney| {
            (
                edges_of(graph, spinney),
                spinney.renormalization_scheme,
                spinney.thermal_zero_mode().is_some(),
            )
        })
        .sorted_by(|a, b| (&a.0, a.2).cmp(&(&b.0, b.2)))
        .collect()
}

#[test]
fn ring_hosts_follow_the_eft_natural_policy() -> Result<()> {
    test_initialise()?;
    let settings = ring_uv_settings(ApproximationType::IR);
    for ring in Ring::ALL {
        let graph = ring.physical_graph()?;
        let all = edge_set(&[0, 1, 2, 3, 4, 5][..graph.underlying.n_edges()]);
        let before = classified(&graph, &settings)?;
        for vacuum_subtraction in [true, false] {
            let context = format!("{ring:?}, vacuum subtraction {vacuum_subtraction}");
            let spinneys =
                zero_mode_policy(&graph, &settings, &thermal_medium(vacuum_subtraction))?;
            let regions = regions(&graph, &spinneys);
            assert_eq!(
                regions,
                ring.expected_regions(&graph, vacuum_subtraction),
                "{context}: {:#?}",
                classification(&graph, &spinneys)
            );

            // The union host read through the kernel's accessors.
            let union = edge_set(&[[2, 3].as_slice(), ring.second_loop_edges()].concat());
            let union_region = &regions[&union];
            assert_eq!(union_region.omega, -1, "{context}");
            assert_eq!(union_region.classes()[0].edges(), &edge_set(&[0, 1]));
            assert_eq!(union_region.tuples()[0].orders(), [4]);
            assert_eq!(
                union_region.zero_bins.len(),
                usize::from(matches!(ring, Ring::Mixed | Ring::Ghost)),
                "{context}"
            );
            assert_eq!(union_region.coordinate_loop_edges(), edge_list(&[2, 4]));

            // Hosts keep IR; every other operation, including the plain union
            // of the union host, is a plain U.
            for (edges, scheme, host) in classification(&graph, &spinneys)
                .into_iter()
                .filter(|(edges, _, _)| !edges.is_empty())
            {
                assert_eq!(
                    scheme == ApproximationType::IR,
                    host,
                    "{context}: {edges:?} {scheme}"
                );
                if !host {
                    assert_eq!(scheme, ApproximationType::MUV, "{context}: {edges:?}");
                }
            }
            assert!(
                spinneys
                    .iter()
                    .any(|spinney| edges_of(&graph, spinney) == union
                        && spinney.thermal_zero_mode().is_none()),
                "{context}: the plain union must stay for the factorized product"
            );
            // The whole-graph U is the R* of the static cograph.
            assert!(
                spinneys
                    .iter()
                    .any(|spinney| edges_of(&graph, spinney) == all
                        && spinney.renormalization_scheme == ApproximationType::MUV),
                "{context}: the full graph must be retained as a plain U"
            );

            // Every divergent IR spinney that hosts nothing is demoted.
            let demoted = before
                .iter()
                .filter(|spinney| spinney.renormalization_scheme == ApproximationType::IR)
                .map(|spinney| edges_of(&graph, spinney))
                .filter(|edges| !edges.is_empty() && !regions.contains_key(edges))
                .collect::<BTreeSet<_>>();
            let mut expected = vec![all.clone(), edge_set(ring.second_loop_edges())];
            if matches!(ring, Ring::Quark | Ring::Tadpole) {
                // The quark bubble has a fermionic soft set; the tadpole's
                // host A is pruned.
                expected.push(edge_set(&[2, 3]));
            }
            assert!(
                expected.iter().all(|edges| demoted.contains(edges)),
                "{context}: demoted {demoted:?}"
            );
            for spinney in spinneys.iter().filter(|spinney| {
                demoted.contains(&edges_of(&graph, spinney))
                    && spinney.thermal_zero_mode().is_none()
            }) {
                assert_eq!(spinney.renormalization_scheme, ApproximationType::MUV);
            }
        }
    }
    Ok(())
}

#[test]
fn static_ir_degree_counts_only_vertices_inside_sigma() -> Result<()> {
    test_initialise()?;
    for ring in Ring::ALL {
        let graph = ring.physical_graph()?;
        // The ring gluons: no vertex lies inside Σ.
        assert_eq!(
            graph.static_ir_degree(&graph.edge_subgraph(&[0, 1]), 1),
            -1,
            "{ring:?}"
        );
    }
    let degree = |ring: Ring, sigma: &[usize], loops: usize| -> Result<i32> {
        let graph = ring.physical_graph()?;
        Ok(graph.static_ir_degree(&graph.edge_subgraph(sigma), loops))
    };
    // Host A: two three-gluon (ghost-gluon) vertices of dod +1 lie inside Σ.
    assert_eq!(degree(Ring::Mixed, &[0, 1, 4, 5], 2)?, 0);
    assert_eq!(degree(Ring::Ghost, &[0, 1, 4, 5], 2)?, 0);
    // The tadpole's four-gluon vertex has dod 0.
    assert_eq!(degree(Ring::Tadpole, &[0, 1, 4], 2)?, 0);
    // D1 = A + {e0, e1, e4}: one soft gluon e5 whose vertices are hard.
    assert_eq!(degree(Ring::Mixed, &[5], 1)?, 1);
    Ok(())
}

#[test]
fn ghost_lines_are_bosonic_zero_mode_classes() -> Result<()> {
    test_initialise()?;
    let graph = Ring::Ghost.physical_graph()?;
    for edge in [4, 5].map(EdgeIndex) {
        assert!(
            !graph[edge].is_fermion(),
            "the ghost {edge} has a zero mode"
        );
    }
    let spinneys = zero_mode_policy(
        &graph,
        &ring_uv_settings(ApproximationType::IR),
        &thermal_medium(true),
    )?;
    let regions = regions(&graph, &spinneys);
    let host = &regions[&edge_set(&[2, 3])];
    assert_eq!(host.kind(), ZeroModeHostKind::Connected);
    assert_eq!(host.classes(), [class(&[0, 1]), class(&[4])]);
    Ok(())
}

#[test]
fn tadpole_inner_host_is_pruned_and_certified_by_the_union_host() -> Result<()> {
    test_initialise()?;
    let graph = Ring::Tadpole.physical_graph()?;
    for vacuum_subtraction in [true, false] {
        let spinneys = zero_mode_policy(
            &graph,
            &ring_uv_settings(ApproximationType::IR),
            &thermal_medium(vacuum_subtraction),
        )?;
        let regions = regions(&graph, &spinneys);
        assert_eq!(
            regions.keys().cloned().collect_vec(),
            [edge_set(&[2, 3, 4])],
            "only the union host A+T remains"
        );
        assert!(spinneys.iter().any(|spinney| {
            edges_of(&graph, spinney) == edge_set(&[2, 3])
                && spinney.renormalization_scheme == ApproximationType::MUV
        }));

        // The joint soft set {e0, e1, e4} of the pruned host A is certified by
        // the union host: e4 = k is safe in its chart, and the remainder
        // converges by ω + J + 1 = 2.
        let union = &regions[&edge_set(&[2, 3, 4])];
        let omega = graph.static_ir_degree(&graph.edge_subgraph(&[0, 1, 4]), 2);
        assert_eq!(omega, 0);
        assert_eq!(i64::from(omega) + union.static_taylor_order() + 1, 2);
    }
    Ok(())
}

#[test]
fn natural_orders_satisfy_the_nested_ir_and_uv_conditions() -> Result<()> {
    test_initialise()?;
    for ring in [Ring::Mixed, Ring::Ghost] {
        let graph = ring.physical_graph()?;
        for vacuum_subtraction in [true, false] {
            let spinneys = zero_mode_policy(
                &graph,
                &ring_uv_settings(ApproximationType::IR),
                &thermal_medium(vacuum_subtraction),
            )?;
            let regions = regions(&graph, &spinneys);
            let inner = &regions[&edge_set(&[2, 3])];
            let outer = &regions[&edge_set(&[2, 3, 4, 5])];
            let static_order = |region: &ThermalZeroModeRegion| region.tuples()[0].taylor_order();
            let (j_in, j_out) = (static_order(inner), static_order(outer));
            // ω_Δ: the degree the extra soft lines Δ = Σ_in \ Σ_out add, from
            // the independently counted soft sets.
            let degree = |sigma: &BTreeSet<EdgeIndex>, loops: usize| {
                let sigma = sigma.iter().map(|edge| edge.0).collect_vec();
                i64::from(graph.static_ir_degree(&graph.edge_subgraph(&sigma), loops))
            };
            let omega_delta =
                degree(inner.sigma(), inner.soft_loops) - degree(outer.sigma(), outer.soft_loops);
            // The static B loop: a linear UV degree, +1.
            assert_eq!(omega_delta, 1, "{ring:?}");
            // IR: the inner order reaches the outer one less the extra
            // degree, J_in >= J_out - ω_Δ, with equality for J = -ω.
            assert_eq!(j_in, j_out - omega_delta, "{ring:?}");
            // UV: the outer operation expands the static cograph Δ = B far
            // enough, d_Δ <= J_out - J_in. Its UV degree d_Δ is counted
            // independently: the 4D degree of B, with its vertices, less the
            // energy integral of its loop.
            let static_cograph = graph.edge_subgraph(&[4, 5]);
            let d_delta = i64::from(graph.local_dod(&static_cograph)) - 1;
            assert_eq!(d_delta, 1, "{ring:?}");
            assert!(d_delta <= j_out - j_in, "{ring:?}");
        }
    }
    Ok(())
}

#[test]
fn zero_mode_tuples_follow_the_soft_degree() {
    // The ring: one class nu = 2 with ω = -1.
    assert_eq!(
        ZeroModeTuple::required(&[2], -1, true),
        Ok(Some(vec![
            tuple(&[4], ZeroModeSector::Complete, 1),
            tuple(&[3], ZeroModeSector::HardVacuum, 0),
        ]))
    );
    assert_eq!(
        ZeroModeTuple::required(&[2], -1, false),
        Ok(Some(vec![tuple(&[4], ZeroModeSector::Complete, 1)]))
    );
    // Host A: classes nu = 2 and nu = 1 with ω = 0. The vacuum pole (3, 1)
    // has J = -2, and mixed tuples vanish.
    for vacuum_subtraction in [true, false] {
        assert_eq!(
            ZeroModeTuple::required(&[2, 1], 0, vacuum_subtraction),
            Ok(Some(vec![tuple(&[4, 2], ZeroModeSector::Complete, 0)]))
        );
    }
    // A convergent static mode needs nothing.
    assert_eq!(ZeroModeTuple::required(&[1], 1, true), Ok(None));
    // A class nu = 1 has no order below its vacuum pole.
    assert_eq!(
        ZeroModeTuple::required(&[1], -3, true),
        Ok(Some(vec![
            tuple(&[2], ZeroModeSector::Complete, 3),
            tuple(&[1], ZeroModeSector::HardVacuum, 2),
        ]))
    );
    // Two classes nu = 1 without vacuum subtraction: only the static mode.
    assert_eq!(
        ZeroModeTuple::required(&[1, 1], -2, false),
        Ok(Some(vec![tuple(&[2, 2], ZeroModeSector::Complete, 2)]))
    );
}

#[test]
fn zero_mode_tuples_reject_orders_below_the_vacuum_pole() {
    for vacuum_subtraction in [true, false] {
        assert_eq!(
            ZeroModeTuple::required(&[2], -3, vacuum_subtraction),
            Err(UnsupportedOrders::BelowVacuumPole {
                orders: vec![2],
                taylor_order: 1,
            })
        );
        // Already the order two below the static one is needed at ω = -2.
        assert_eq!(
            ZeroModeTuple::required(&[2], -2, vacuum_subtraction),
            Err(UnsupportedOrders::BelowVacuumPole {
                orders: vec![2],
                taylor_order: 0,
            })
        );
    }
}

#[test]
fn zero_mode_tuples_reject_multi_class_vacuum_poles() {
    assert_eq!(
        ZeroModeTuple::required(&[1, 1], -2, true),
        Err(UnsupportedOrders::MultiClassVacuumPole {
            orders: vec![1, 1],
            taylor_order: 0,
        })
    );
}

#[test]
fn hostless_soft_sets_are_rejected() -> Result<()> {
    test_initialise()?;
    // A pure-glue ring: e3 = p - q and e5 = k - q are unsafe in each bubble,
    // so no bubble can host the soft set of the other, nor the ring.
    let graph: Graph = "digraph pure_glue_ring {
        a1:0 -> b0:1 [id=0 dir=none particle=g lmb_id=0];
        b1:2 -> a0:3 [id=1 dir=none particle=g];
        a0:4 -> a1:5 [id=2 dir=none particle=g lmb_id=1];
        a1:6 -> a0:7 [id=3 dir=none particle=g];
        b0:8 -> b1:9 [id=4 dir=none particle=g lmb_id=2];
        b1:10 -> b0:11 [id=5 dir=none particle=g];
    }"
    .into_graph(&load_generic_model("sm"))?;
    let error = zero_mode_policy(
        &graph,
        &ring_uv_settings(ApproximationType::IR),
        &thermal_medium(true),
    )
    .expect_err("the pure-glue ring has no zero-mode host");
    let ZeroModeHostError::Hostless { sigma, omega } = &error else {
        panic!("unexpected error: {error}");
    };
    assert_eq!((sigma, *omega), (&edge_set(&[0, 1]), -1));

    // With the default MUV prescription no quark-ring bubble is a candidate.
    let graph = Ring::Quark.physical_graph()?;
    let muv = ring_uv_settings(ApproximationType::MUV);
    let error = zero_mode_policy(&graph, &muv, &thermal_medium(true))
        .expect_err("the MUV prescription leaves the quark ring without a host");
    assert!(matches!(
        &error,
        ZeroModeHostError::Hostless { sigma, omega: -1 } if *sigma == edge_set(&[0, 1])
    ));
    let error = validate_local_counterterm_schemes(
        &graph,
        &CutStructure::empty(&graph),
        &muv,
        &thermal_medium(false),
    )
    .expect_err("the scheme validation runs the host policy");
    let message = format!("{error:#}");
    assert!(
        message.contains(&format!("graph '{}' cut 0", graph.name)),
        "{message}"
    );
    assert!(
        message.contains("renormalization_prescription.massless_power_divergent = IR"),
        "{message}"
    );
    Ok(())
}

/// A quark ring whose second bubble is a top triangle B = {e4, e5, e6},
/// attached through the zero-momentum Higgs bridge e7 to the top tadpole
/// T = {e8}.
fn bridged_ring() -> Result<Graph> {
    "digraph bridged_quark_ring {
        a1:0 -> b0:1 [id=0 dir=none particle=g lmb_id=0];
        b1:2 -> a0:3 [id=1 dir=none particle=g];
        a0:4 -> a1:5 [id=2 particle=d lmb_id=1];
        a1:6 -> a0:7 [id=3 particle=d];
        b0:8 -> b2:9 [id=4 particle=t lmb_id=2];
        b2:10 -> b1:11 [id=5 particle=t];
        b1:12 -> b0:13 [id=6 particle=t];
        b2:14 -> c0:15 [id=7 dir=none particle=H];
        c0:16 -> c0:17 [id=8 particle=t lmb_id=3];
    }"
    .into_graph(&load_generic_model("sm"))
}

#[test]
fn tree_lines_belong_to_no_soft_set() -> Result<()> {
    test_initialise()?;
    let graph = bridged_ring()?;
    let bridge = &graph.loop_momentum_basis.edge_signatures[EdgeIndex(7)];
    assert!(
        bridge.internal.iter().all(|sign| sign.is_zero()),
        "the Higgs bridge carries no loop momentum"
    );
    let settings = ring_uv_settings(ApproximationType::IR);
    for vacuum_subtraction in [true, false] {
        let spinneys = zero_mode_policy(&graph, &settings, &thermal_medium(vacuum_subtraction))?;
        let mut tuples = vec![tuple(&[4], ZeroModeSector::Complete, 1)];
        if vacuum_subtraction {
            tuples.push(tuple(&[3], ZeroModeSector::HardVacuum, 0));
        }
        // The two gapped bubbles and the tadpole host the ring gluons
        // jointly; the bridge lies in neither the host nor its soft set.
        let union_region = ThermalZeroModeRegion {
            kind: ZeroModeHostKind::Union,
            sigma: edge_set(&[0, 1]),
            classes: vec![class(&[0, 1])],
            omega: -1,
            soft_loops: 1,
            tuples,
            zero_bins: Vec::new(),
            coordinate_loop_edges: edge_list(&[2, 4, 8]),
        };
        assert_eq!(
            regions(&graph, &spinneys),
            BTreeMap::from([(edge_set(&[2, 3, 4, 5, 6, 8]), union_region)]),
            "vacuum subtraction {vacuum_subtraction}: {:#?}",
            classification(&graph, &spinneys)
        );
    }

    // Without a host, the divergent ring soft set is reported rather than
    // hidden behind the bridge.
    let error = zero_mode_policy(
        &graph,
        &ring_uv_settings(ApproximationType::MUV),
        &thermal_medium(true),
    )
    .expect_err("the MUV prescription leaves the bridged ring without a host");
    assert!(
        matches!(
            &error,
            ZeroModeHostError::Hostless { sigma, omega: -1 } if *sigma == edge_set(&[0, 1])
        ),
        "{error}"
    );
    Ok(())
}

#[test]
fn zero_mode_hosts_reject_unsupported_routes() -> Result<()> {
    test_initialise()?;
    let graph = Ring::Quark.physical_graph()?;
    let settings = ring_uv_settings(ApproximationType::IR);
    let empty = CutStructure::empty(&graph);
    validate_local_counterterm_schemes(&graph, &empty, &settings, &thermal_medium(true))?;

    // Spinneys are merged across cuts, so hosts need one uncut structure.
    let two_cuts = CutStructure {
        cuts: vec![
            CutSet::empty(graph.n_hedges()),
            CutSet::empty(graph.n_hedges()),
        ],
    };
    let error =
        validate_local_counterterm_schemes(&graph, &two_cuts, &settings, &thermal_medium(true))
            .expect_err("a host on a non-trivial cut structure");
    assert!(matches!(
        error.root_cause().downcast_ref::<ZeroModeHostError>(),
        Some(ZeroModeHostError::NonTrivialCut { host }) if *host == edge_set(&[2, 3, 4, 5])
    ));

    // Routes that do not freeze soft branches cannot subtract the zero mode.
    for (route, settings) in [
        (
            "legacy",
            UVgenerationSettings {
                orchestrator: UVOrchestrator::LegacyDagForest,
                ..settings.clone()
            },
        ),
        (
            "compare",
            UVgenerationSettings {
                orchestrator: UVOrchestrator::Compare,
                ..settings.clone()
            },
        ),
    ] {
        for vacuum_subtraction in [true, false] {
            let error = validate_local_counterterm_schemes(
                &graph,
                &empty,
                &settings,
                &thermal_medium(vacuum_subtraction),
            )
            .expect_err(route);
            let message = error.to_string();
            assert!(
                message.contains("thermal soft counterterm component")
                    && message.contains("uv.orchestrator=hedge_poset")
                    && !message.contains("generate_integrated"),
                "{route}: {message}"
            );
        }
    }

    // Integrated counterterms keep the route's frozen soft branches, so the
    // union host is assigned as without them.
    let integrated = UVgenerationSettings {
        generate_integrated: true,
        ..settings.clone()
    };
    for vacuum_subtraction in [true, false] {
        let medium = thermal_medium(vacuum_subtraction);
        validate_local_counterterm_schemes(&graph, &empty, &integrated, &medium)?;
        let mut spinneys = classified(&graph, &integrated)?;
        assign_zero_mode_hosts(
            &graph,
            &uncut(&graph),
            &mut spinneys,
            &integrated,
            &medium,
            &empty,
        )?;
        let hosts = spinneys
            .iter()
            .filter(|spinney| spinney.thermal_zero_mode().is_some())
            .map(|spinney| edges_of(&graph, spinney))
            .collect::<Vec<_>>();
        assert_eq!(hosts, [edge_set(&[2, 3, 4, 5])], "VS {vacuum_subtraction}");
    }

    // Without the full observable the static cograph of the host has no U.
    for vacuum_subtraction in [true, false] {
        let medium = thermal_medium(vacuum_subtraction);
        let mut spinneys = classified(&graph, &settings)?;
        assign_zero_mode_hosts(
            &graph,
            &uncut(&graph),
            &mut spinneys,
            &settings,
            &medium,
            &empty,
        )?;
        let full = edge_set(&[0, 1, 2, 3, 4, 5]);
        spinneys.retain(|spinney| edges_of(&graph, spinney) != full);
        let error = graph
            .restrict_full_observable_spinney(&mut spinneys, &medium)
            .expect_err("a host needs its enclosing whole-graph U");
        assert!(matches!(
            &error,
            ZeroModeHostError::MissingEnclosingSpinney { host, dod: 4 }
                if *host == edge_set(&[2, 3, 4, 5])
        ));
    }
    Ok(())
}

#[test]
fn zero_mode_policy_is_skipped_without_uv_subtraction() -> Result<()> {
    test_initialise()?;
    let settings = UVgenerationSettings {
        subtract_uv: false,
        ..ring_uv_settings(ApproximationType::IR)
    };
    for ring in Ring::ALL {
        let graph = ring.physical_graph()?;
        for vacuum_subtraction in [true, false] {
            let spinneys =
                zero_mode_policy(&graph, &settings, &thermal_medium(vacuum_subtraction))?;
            assert_eq!(spinneys.len(), 1, "{ring:?}");
            assert!(spinneys[0].filter().is_empty(), "{ring:?}");
            validate_local_counterterm_schemes(
                &graph,
                &CutStructure::empty(&graph),
                &settings,
                &thermal_medium(vacuum_subtraction),
            )?;
        }
    }
    Ok(())
}

#[test]
fn zero_mode_hosts_without_vacuum_subtraction_only_drop_the_vacuum_order() -> Result<()> {
    test_initialise()?;
    let settings = ring_uv_settings(ApproximationType::IR);
    for ring in Ring::ALL {
        let graph = ring.physical_graph()?;
        let with = zero_mode_policy(&graph, &settings, &thermal_medium(true))?;
        let without = zero_mode_policy(&graph, &settings, &thermal_medium(false))?;
        let static_only = |mut region: ThermalZeroModeRegion| {
            region
                .tuples
                .retain(|tuple| tuple.sector == ZeroModeSector::Complete);
            region
        };
        let with_regions = regions(&graph, &with)
            .into_iter()
            .map(|(edges, region)| (edges, static_only(region)))
            .collect::<BTreeMap<_, _>>();
        let without_regions = regions(&graph, &without);
        assert_eq!(with_regions, without_regions, "{ring:?}");
        assert!(
            without_regions.values().all(|region| {
                region
                    .tuples()
                    .iter()
                    .all(|tuple| tuple.sector() == ZeroModeSector::Complete)
            }),
            "{ring:?}"
        );
        // Without vacuum subtraction no spinney is restricted: the policy
        // only adds the union host copy.
        let union_hosts = without_regions
            .values()
            .filter(|region| region.kind() == ZeroModeHostKind::Union)
            .count();
        assert_eq!(
            without.len(),
            classified(&graph, &settings)?.len() + union_hosts,
            "{ring:?}"
        );
        assert_eq!(
            classification(&graph, &with),
            classification(&graph, &without)
        );
    }
    Ok(())
}

#[test]
fn zero_mode_policy_leaves_vacuum_and_zero_temperature_unchanged() -> Result<()> {
    test_initialise()?;
    let settings = ring_uv_settings(ApproximationType::IR);
    for ring in Ring::ALL {
        let graph = ring.physical_graph()?;
        let before = classified(&graph, &settings)?;
        for medium in [
            medium(MediumMode::Vacuum, false),
            medium(MediumMode::ZeroTemperatureEquilibrium, false),
            medium(MediumMode::ZeroTemperatureEquilibrium, true),
        ] {
            let mut spinneys = before.clone();
            assign_zero_mode_hosts(
                &graph,
                &uncut(&graph),
                &mut spinneys,
                &settings,
                &medium,
                &CutStructure::empty(&graph),
            )?;
            assert_eq!(
                classification(&graph, &spinneys),
                classification(&graph, &before),
                "{ring:?}: {medium:?}"
            );
        }
    }
    Ok(())
}

#[test]
fn zero_mode_hosts_label_their_projection_provenance() -> Result<()> {
    test_initialise()?;
    let graph = Ring::Mixed.physical_graph()?;
    let settings = ring_uv_settings(ApproximationType::IR);
    let ctx = UVCtx::new(&graph, &settings);
    let root = Spinney::empty(&graph);
    let provenance = |spinney: &Spinney| {
        let current = ForestNode {
            spinney,
            topo_order: 0,
        };
        let given = ForestNode {
            spinney: &root,
            topo_order: 0,
        };
        let step = Local3DLoopRescaling::FullSubgraph.projection_step(
            &ctx,
            &current,
            &given,
            current.subgraph(),
            None,
            current.lmb(),
            current.lmb(),
        );
        (step.conceptual_branches, step.materialization)
    };
    let signed = |branch, coefficient| Local3DSignedBranch {
        branch,
        coefficient,
    };
    let soft_and_hard = (
        vec![
            signed(Local3DConceptualBranch::U, 1),
            signed(Local3DConceptualBranch::S, 1),
            signed(Local3DConceptualBranch::US, -1),
        ],
        Local3DMaterialization::FactorizedSoft,
    );

    let spinneys = zero_mode_policy(&graph, &settings, &thermal_medium(true))?;
    let spinney = |edges: &[usize], kind: Option<ZeroModeHostKind>| {
        spinneys
            .iter()
            .find(|spinney| {
                edges_of(&graph, spinney) == edge_set(edges)
                    && spinney.thermal_zero_mode().map(ThermalZeroModeRegion::kind) == kind
            })
            .expect("the mixed ring has the spinney")
    };
    // The connected host A: S0 + U - U S0.
    assert_eq!(
        provenance(spinney(&[2, 3], Some(ZeroModeHostKind::Connected))),
        soft_and_hard
    );
    // The union host: S0 alone.
    assert_eq!(
        provenance(spinney(&[2, 3, 4, 5], Some(ZeroModeHostKind::Union))),
        (
            vec![signed(Local3DConceptualBranch::S, 1)],
            Local3DMaterialization::FactorizedSoft
        )
    );
    // Without a region the scheme decides: a plain U, and the vacuum soft
    // subtraction of an IR spinney of positive degree.
    assert_eq!(
        provenance(spinney(&[4, 5], None)),
        (
            vec![signed(Local3DConceptualBranch::U, 1)],
            Local3DMaterialization::DirectU
        )
    );
    let unassigned = classified(&graph, &settings)?
        .into_iter()
        .find(|spinney| edges_of(&graph, spinney) == edge_set(&[2, 3]))
        .expect("the mixed ring has the bubble A");
    assert_eq!(unassigned.renormalization_scheme, ApproximationType::IR);
    assert!(unassigned.thermal_zero_mode().is_none());
    assert_eq!(provenance(&unassigned), soft_and_hard);
    Ok(())
}

#[test]
fn ring_woods_host_a_zero_mode_only_when_subtracting_at_finite_temperature() -> Result<()> {
    test_initialise()?;
    for ring in Ring::ALL {
        let graph = ring.physical_graph()?;
        for subtract_uv in [true, false] {
            for medium in [
                thermal_medium(true),
                thermal_medium(false),
                medium(MediumMode::ZeroTemperatureEquilibrium, true),
            ] {
                let settings = GenerationSettings {
                    uv: UVgenerationSettings {
                        subtract_uv,
                        ..ring_uv_settings(ApproximationType::IR)
                    },
                    medium: medium.clone(),
                    ..Default::default()
                };
                let wood = Wood::new(CutStructure::empty(&graph), &graph, &settings)?;
                assert_eq!(
                    wood.hosts_zero_mode(),
                    subtract_uv && medium.mode.is_finite_temperature(),
                    "{ring:?}: subtract_uv={subtract_uv}, {medium:?}"
                );
            }
        }
    }
    Ok(())
}

/// An integrated prefix holds each frequency class of a host wholly or not
/// at all, or it contracts part of one: on the gluon ring, host A has the
/// classes `{e0, e1}` and `{e4}`; its own lines and the second loop B, which
/// holds the class `{e4}` and the implied line `e5`, contract none partly,
/// while a prefix with only one of the two gluons contracts part of
/// `{e0, e1}`.
#[test]
fn integrated_prefixes_contract_whole_frequency_classes() -> Result<()> {
    test_initialise()?;
    let graph = Ring::Mixed.topology()?;
    let regions = Ring::Mixed.expected_regions(&graph, true);
    let host = &regions[&edge_set(&[2, 3])];
    for whole in [&[2, 3][..], &[4, 5], &[0, 1, 4, 5], &[0, 1, 2, 3, 4]] {
        assert_eq!(
            host.partly_contracted_class(&edge_set(whole)),
            None,
            "{whole:?}"
        );
    }
    for partial in [&[0][..], &[1, 4, 5]] {
        assert_eq!(
            host.partly_contracted_class(&edge_set(partial)),
            Some(&class(&[0, 1])),
            "{partial:?}"
        );
    }
    Ok(())
}
