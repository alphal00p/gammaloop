use super::*;
use crate::{
    cff::{VertexSet, esurface::Esurface},
    dot,
    graph::parse::from_dot::IntoGraph,
    initialisation::test_initialise,
    momentum::{FourMomentum, Rotatable, RotationMethod},
};
use linnet::half_edge::involution::EdgeIndex;
use typed_index_collections::ti_vec;

#[test]
fn gl638_c7_keeps_a_certified_full_overlap_after_solver_failure() {
    test_initialise().unwrap();
    // GL638's directed topology with scalar numerators and common parent [4,10,7,12].
    // The paired initial-state halfedges are sewn exactly as in the imported graph.
    // Only routing and masses enter this test; no numerator or integrand is evaluated.
    let graph: Graph = dot!(digraph GL638_C7 {
        edge [num=1 mass=0]
        node [num=1]
        exte0 [style=invis]
        exte1 [style=invis]
        exte15 [style=invis]
        exte16 [style=invis]
        exte0 -> v9:0 [id=0 is_cut=28]
        exte1 -> v9:1 [id=1 is_cut=29]
        v0:2 -> v1:3 [id=2]
        v0:4 -> v6:5 [id=3]
        v7:6 -> v0:7 [id=4 lmb_id=0]
        v3:8 -> v1:9 [id=5]
        v1:10 -> v7:11 [id=6]
        v2:12 -> v5:13 [id=7 lmb_id=2]
        v6:14 -> v2:15 [id=8]
        v2:16 -> v9:17 [id=9]
        v4:18 -> v3:19 [id=10 lmb_id=1]
        v3:20 -> v8:21 [id=11]
        v5:22 -> v4:23 [id=12 lmb_id=3]
        v4:24 -> v7:25 [id=13]
        v5:26 -> v6:27 [id=14]
        v8:28 -> exte15 [id=15 is_cut=28]
        v8:29 -> exte16 [id=16 is_cut=29]
    })
    .unwrap();
    let lmbs = ti_vec![graph.loop_momentum_basis.clone()];
    assert_eq!(
        graph
            .loop_momentum_basis
            .loop_edges
            .iter()
            .map(|edge| usize::from(*edge))
            .collect_vec(),
        vec![4, 10, 7, 12]
    );
    let subspace = SubspaceData::new_from_parent_basis_edges(
        &[EdgeIndex::from(7)],
        &graph.full_filter(),
        LmbIndex::from(0),
        &graph,
        &lmbs,
    )
    .unwrap();
    assert_eq!(
        subspace.iter_lmb_indices().collect_vec(),
        vec![LoopIndex::from(2)]
    );
    // Exact fixed momenta from the live all-orientation point5586 SOCP capture.
    // Three last-bit differences from recomputing cut roots change Clarabel's status:
    // identity reported NumericalError and z InsufficientProgress, although both
    // candidates lay >70GeV inside every surface. Do not assert platform-specific status.
    let p3 = [-211.12977078180307, 54.95284617175551, -368.5885668818366];
    let p4 = [-177.57308684481387, 46.218714156193634, -310.0055920799352];
    let t4 = [6.095429701834923, -19.027527107331455, -28.559764597935715];
    let p0 = [-154.54706340567736, 40.22550192793898, -269.80695524633524];
    // Keep all seven instances, including the four identical A thresholds. Only
    // p=k12 and t=k10 enter the fixed complement; the other inactive coordinate is zero.
    let surfaces = [
        (vec![7, 8], p3, [0.0; 3]),
        (vec![8, 12, 14], p3, [0.0; 3]),
        (vec![7, 8], p4, t4),
        (vec![8, 10, 13, 14], p4, t4),
        (vec![7, 8], [0.0; 3], [0.0; 3]),
        (vec![7, 8], p0, [0.0; 3]),
        (vec![8, 12, 14], p0, [0.0; 3]),
    ];
    let thresholds: EsurfaceCollection = surfaces
        .iter()
        .map(|(edges, _, _)| Esurface {
            energies: edges.iter().copied().map(EdgeIndex::from).collect(),
            external_shift: vec![(EdgeIndex::from(0), -1), (EdgeIndex::from(1), -1)],
            vertex_set: VertexSet::dummy(),
        })
        .collect::<Vec<_>>()
        .into();
    let external: ExternalFourMomenta<_> = [
        FourMomentum::from_args(F(500.0), F(0.0), F(0.0), F(500.0)),
        FourMomentum::from_args(F(500.0), F(0.0), F(0.0), F(-500.0)),
    ]
    .into_iter()
    .collect();
    let masses = graph.underlying.new_edgevec(|_, edge, _| {
        F(if [3, 4, 5, 6, 7, 8, 10, 12].contains(&usize::from(edge)) {
            173.0
        } else {
            0.0
        })
    });
    let mut settings = RuntimeSettings::default();
    settings.kinematics.e_cm = 1000.0;
    let existing: ExistingThresholds = thresholds.keys().collect();
    let selected = (0..7).map(ExistingEsurfaceId::from).collect_vec();
    let mut identity_center: Option<LoopMomenta<F<f64>>> = None;
    for method in [
        RotationMethod::Identity,
        RotationMethod::Pi2Z,
        RotationMethod::EulerAngles(0.1, 0.2, 0.3),
    ] {
        let rotation = Rotation::new(method);
        let kinematics = surfaces
            .iter()
            .map(|(_, p, t)| {
                let loops: LoopMomenta<_> = [[0.0; 3], *t, [0.0; 3], *p]
                    .into_iter()
                    .map(|[x, y, z]| ThreeMomentum::new(F(x), F(y), F(z)))
                    .collect();
                OverlapKinematics {
                    loop_moms: loops.rotate(&rotation),
                    external_momenta: external
                        .iter()
                        .map(|p| FourMomentum {
                            temporal: p.temporal,
                            spatial: p.spatial.rotate(&rotation),
                        })
                        .collect(),
                    edge_masses: None,
                }
            })
            .collect_vec();
        let input = OverlapInput {
            graph: &graph,
            settings: &settings,
            subspace: &subspace,
            threshold_subspaces: None,
            lmbs: &lmbs,
            thresholds: &thresholds,
            edge_masses: masses.clone(),
            surface_kinematics: Some(&kinematics),
        };
        let first = &kinematics[0];
        let center = find_center(
            &input,
            &selected,
            &existing,
            &first.loop_moms,
            &first.external_momenta,
            false,
        )
        .unwrap()
        .expect("the full seven-surface overlap has a strictly interior witness");
        assert!(check_global_center(
            &input,
            &existing,
            &center,
            &first.loop_moms,
            &first.external_momenta
        ));
        if let Some(identity) = &identity_center {
            let expected = identity.rotate(&rotation);
            for (actual, expected) in center[LoopIndex::from(2)]
                .into_iter()
                .zip(&expected[LoopIndex::from(2)])
            {
                assert!(
                    (actual.0 - expected.0).abs() < 2e-6,
                    "center failed rotation covariance for {method}"
                );
            }
        } else {
            identity_center = Some(center.clone());
        }
        let overlap = find_maximal_overlap(
            &input,
            &existing,
            &first.loop_moms,
            &first.external_momenta,
            &rotation,
        )
        .unwrap();
        assert_eq!(
            overlap.overlap_groups.len(),
            1,
            "lost full overlap for {method}"
        );
        assert_eq!(overlap.overlap_groups[0].existing_esurfaces, selected);
        assert!(overlap.overlap_groups[0].complement.is_empty());
        assert!(check_global_center(
            &input,
            &existing,
            &overlap.overlap_groups[0].center,
            &first.loop_moms,
            &first.external_momenta
        ));
    }
}
