use super::*;
use crate::{
    cff::{VertexSet, esurface::RaisedEsurfaceGroup},
    dot,
    graph::{lmb::LMBext, parse::from_dot::IntoGraph},
    initialisation::test_initialise,
    momentum::{FourMomentum, RotationMethod, sample::BareMomentumSample},
    processes::{
        CrossSectionGraph, ResolvedCutGroupThresholdCounterterms,
        ResolvedThresholdCountertermVariant, ResolvedThresholdCounterterms,
        ThresholdCountertermVariantStatus,
    },
    subtraction::generate_rstar_t_dependence_evaluator,
};
use symbolica::atom::Atom;

#[test]
fn shared_group_weights_preserve_foreign_cut_data_and_radial_derivatives() {
    test_initialise().unwrap();
    let graph: Graph = dot!(digraph shared_weights {
        ext [style=invis]
        edge [num=1 mass=0]
        node [num=1]
        ext->a:0 [id=0]
        a->b [id=1]
        b->a [id=2]
        a->b [id=3]
        ext->b:1 [id=4]
    })
    .unwrap();
    let mut alternate = graph.loop_momentum_basis.clone();
    alternate.swap_loops(LoopIndex(0), LoopIndex(1));
    let lmbs = ti_vec![graph.loop_momentum_basis.clone(), alternate];
    let selected = lmbs[LmbIndex::from(0)].loop_edges[LoopIndex(0)];
    let support = graph
        .iter_loop_edges()
        .map(|(_, edge, _)| edge)
        .find(|edge| !graph.loop_momentum_basis.loop_edges.contains(edge))
        .unwrap();
    let external = lmbs[LmbIndex::from(0)].ext_edges[ExternalIndex(0)];
    let subspaces = [0, 1].map(|parent| {
        SubspaceData::new_from_parent_basis_edges(
            &[selected],
            &graph.full_filter(),
            LmbIndex::from(parent),
            &graph,
            &lmbs,
        )
        .unwrap()
    });
    assert_eq!(
        subspaces[0].solve_signature(&lmbs),
        subspaces[1].solve_signature(&lmbs)
    );
    let masses = graph.underlying.new_edgevec(|_, _, _| F(0.0));
    fn make_sample<T: FloatLike>(active: f64, fixed: f64, energy: f64) -> MomentumSample<T> {
        MomentumSample {
            sample: BareMomentumSample {
                loop_moms: LoopMomenta::from_iter([active, fixed].map(|px| {
                    ThreeMomentum::new(F::from_f64(px), F::from_f64(0.0), F::from_f64(0.0))
                })),
                dual_loop_moms: None,
                loop_mom_cache_id: 0,
                loop_mom_base_cache_id: 0,
                external_moms: [
                    FourMomentum::from_args(
                        F::from_f64(energy),
                        F::from_f64(0.0),
                        F::from_f64(0.0),
                        F::from_f64(0.0),
                    ),
                    FourMomentum::from_args(
                        F::from_f64(-energy),
                        F::from_f64(0.0),
                        F::from_f64(0.0),
                        F::from_f64(0.0),
                    ),
                ]
                .into_iter()
                .collect(),
                external_mom_cache_id: 0,
                external_mom_base_cache_id: 0,
                jacobian: F::from_f64(1.0),
                orientation: None,
                parameterization_branch: None,
            },
        }
    }
    fn make_point<T: FloatLike>(sample: MomentumSample<T>) -> LUCTKinematicPoint<T> {
        LUCTKinematicPoint {
            unrescaled_sample: sample.clone(),
            dualized_momentum_sample_cache: vec![sample],
            lu_cut_parameter_cache: vec![],
            lu_cut_esurface_values: vec![],
        }
    }
    let surface_a = Esurface {
        energies: vec![selected],
        external_shift: vec![(external, -1)],
        vertex_set: VertexSet::dummy(),
    };
    let surface_bc = Esurface {
        energies: vec![support],
        ..surface_a.clone()
    };
    let signature = &lmbs[LmbIndex::from(0)].edge_signatures[support].internal;
    let active_sign = signature[LoopIndex(0)] as i8 as f64;
    let fixed_sign = signature[LoopIndex(1)] as i8 as f64;
    assert_ne!(fixed_sign, 0.0);
    // At the common active coordinate x=4, B and C have eta=1 and eta=3.
    // Their own fixed coordinates differ from the target cut's fixed coordinate 17.
    let point = make_point(make_sample::<f64>(3.0, 17.0, 4.0));
    let mut overlap = OverlapStructure {
        existing_esurfaces: ti_vec![EsurfaceID(0), EsurfaceID(1), EsurfaceID(2)],
        overlap_groups: [2.0, 0.0]
            .map(|center| OverlapGroup {
                existing_esurfaces: vec![],
                complement: vec![],
                center: LoopMomenta::from_iter(
                    [center, 0.0].map(|px| ThreeMomentum::new(F(px), F(0.0), F(0.0))),
                ),
            })
            .to_vec(),
    };
    overlap.overlap_groups[0].existing_esurfaces =
        vec![ExistingEsurfaceId::from(0), ExistingEsurfaceId::from(1)];
    overlap.overlap_groups[1].existing_esurfaces =
        vec![ExistingEsurfaceId::from(0), ExistingEsurfaceId::from(2)];
    overlap.fill_in_complements();
    let group = Arc::new(ThresholdSolveGroup {
        thresholds: vec![surface_a.clone(), surface_bc.clone(), surface_bc].into(),
        subspace: subspaces[0].clone(),
        kinematics: vec![
            point.representative_sample().clone(),
            make_sample::<f64>(100.0, -2.0 * active_sign / fixed_sign, 1.0),
            make_sample::<f64>(-100.0, 0.0, 1.0),
        ],
        records: (0..3)
            .map(|index| GlobalOverlapRecord {
                cut_group_id: CutGroupId(index),
                side: ThresholdCountertermSide::Left,
                local_threshold_id: 0,
                variant_id: None,
            })
            .collect(),
        overlap,
        prefactor_power: 2,
    });
    // An independently labelled group can have identical coordinates. Its center must not
    // enter this group's denominator, even though both are dispatched by the same cut.
    let independent = Arc::new(ThresholdSolveGroup {
        thresholds: vec![surface_a.clone()].into(),
        subspace: subspaces[0].clone(),
        kinematics: vec![point.representative_sample().clone()],
        records: vec![GlobalOverlapRecord {
            cut_group_id: CutGroupId(0),
            side: ThresholdCountertermSide::Left,
            local_threshold_id: 1,
            variant_id: None,
        }],
        overlap: OverlapStructure {
            existing_esurfaces: ti_vec![EsurfaceID(0)],
            overlap_groups: vec![OverlapGroup {
                existing_esurfaces: vec![ExistingEsurfaceId::from(0)],
                complement: vec![],
                center: group.overlap.overlap_groups[1].center.clone(),
            }],
        },
        prefactor_power: 2,
    });
    // Local dispatch sees A and the independent variant. A's channel normalization must still
    // contain both foreign thresholds, although neither local complement contains them.
    let mut local = OverlapStructure {
        existing_esurfaces: ti_vec![EsurfaceID(0)],
        overlap_groups: group
            .overlap
            .overlap_groups
            .iter()
            .map(|center| OverlapGroup {
                existing_esurfaces: vec![ExistingEsurfaceId::from(0)],
                complement: vec![],
                center: center.center.clone(),
            })
            .collect(),
    };
    local.existing_esurfaces.push(EsurfaceID(1));
    local.overlap_groups.push(OverlapGroup {
        existing_esurfaces: vec![ExistingEsurfaceId::from(1)],
        complement: vec![],
        center: independent.overlap.overlap_groups[0].center.clone(),
    });
    let contexts = [
        (Arc::clone(&group), 0),
        (Arc::clone(&group), 1),
        (independent, 0),
    ];
    let settings = RuntimeSettings::default();
    let rotation = Rotation::new(RotationMethod::Identity);
    let thresholds: EsurfaceCollection = vec![surface_a.clone(), surface_a].into();
    let builder = CounterTermBuilder::new(
        &graph,
        &settings,
        &thresholds,
        point,
        &local,
        &masses,
        &lmbs,
        &subspaces[1],
        None,
        &rotation,
        CutGroupId(0),
        &contexts,
    );
    let root_sample = make_sample::<f64>(4.0, 17.0, 4.0)
        .lmb_transform(&lmbs[LmbIndex::from(0)], &lmbs[LmbIndex::from(1)]);
    let dual = |value, derivative| {
        HyperDual::from_values(simple_n_deriv_shape(1), vec![F(value), F(derivative)])
    };
    let mut dual_momenta = dualize_loop_momenta(&dual(0.0, 0.0), root_sample.loop_moms());
    dual_momenta[LoopIndex(1)].px = dual(4.0, 1.0);
    dual_momenta[LoopIndex(0)].px = dual(17.0, 7.0);
    let geometry = DualRstarGeometry {
        radius: dual(1.0, 0.0),
        radius_star: dual(2.0, 1.0),
        esurface_derivative: dual(1.0, 0.0),
        rstar_loop_momenta: dual_momenta,
        external_moms: dualize_external_momenta(&dual(0.0, 0.0), root_sample.external_moms()),
    };
    for (index, expected, derivative) in [(0, 0.9, -0.12), (1, 0.1, 0.12)] {
        let overlap_builder = builder.new_overlap_builder(index, EsurfaceID(0));
        assert_eq!(
            overlap_builder.center[LoopIndex(1)].px,
            F(if index == 0 { 2.0 } else { 0.0 })
        );
        assert_eq!(overlap_builder.center[LoopIndex(0)].px, F(0.0));
        let solution = RstarSolution {
            esurface_ct_builder: overlap_builder.new_esurface_builder(ExistingEsurfaceId::from(0)),
            solution: NewtonIterationResult {
                solution: F(2.0),
                derivative_at_solution: F(1.0),
                error_of_function: F(0.0),
                num_iterations_used: 0,
            },
            alpha: F(2.0) / overlap_builder.radius,
            alpha_t_dependence: None,
        };
        let weight = solution.non_dual_multichanneling_factor(&root_sample);
        assert!((weight.re.0 - expected).abs() < 1e-12);
        let weight = solution.dual_multichanneling_factor(&geometry);
        assert!((weight.values[0].0 - expected).abs() < 1e-12);
        assert!((weight.values[1].0 - derivative).abs() < 1e-12);
    }

    // A nontrivial center changes the physical radius even though the selected
    // energy is simply |q0|. This closed form independently fixes all t orders.
    {
        let full_subspace = SubspaceData::new_from_parent_basis_edges(
            &lmbs[LmbIndex::from(0)]
                .loop_edges
                .iter()
                .copied()
                .collect_vec(),
            &graph.full_filter(),
            LmbIndex::from(0),
            &graph,
            &lmbs,
        )
        .unwrap();
        let source = make_sample::<f64>(3.0, 0.0, 4.0);
        let center = make_sample::<f64>(0.0, 2.0, 4.0).loop_moms().clone();
        let shape = HyperDual::new(simple_n_deriv_shape(3));
        let t = shape.variable(0, F(1.0));
        let expected =
            ((t.clone() * &t) * F(9.0) + new_constant(&t, &F(4.0))).sqrt() * F(4.0 / 3.0) / &t;
        let alpha = generate_rstar_t_dependence_evaluator(3)
            .unwrap()
            .evaluate_alpha(RstarTDependenceInput {
                t_star: &F(1.0),
                alpha: &(F(4.0) / F(3.0)),
                overlap_center: &center,
                subspace: &full_subspace,
                unrescaled_momentum_sample: &source,
                representative_sample: &source,
                masses: &masses,
                threshold_esurface: &thresholds[EsurfaceID(0)],
                lmb: &graph.loop_momentum_basis,
                all_lmbs: &lmbs,
            });
        let actual = alpha * ((t.clone() * &t) * F(9.0) + new_constant(&t, &F(4.0))).sqrt();
        for (order, (actual, expected)) in actual.values.iter().zip(&expected.values).enumerate() {
            assert!(
                (actual - expected).abs() < F(1.0e-10),
                "physical rstar Taylor order {order}: actual={actual}, expected={expected}"
            );
        }
    }

    // A massive sphere with a nonzero center has alpha(t)=2/(3t-2).
    // The other loop is a null direction of this surface, but contributes to
    // the physical hyper-radius and all mixed threshold-r derivatives.
    fn check_alpha<T: FloatLike>(
        graph: &Graph,
        selected: EdgeIndex,
        support: EdgeIndex,
        surface: &Esurface,
    ) {
        let f = F::<T>::from_f64;
        let masses = graph.underlying.new_edgevec(|_, _, _| f(0.0));
        let settings = RuntimeSettings::default();
        let rotation = Rotation::new(RotationMethod::Identity);
        let alternate = graph
            .generate_loop_momentum_bases()
            .into_iter()
            .find(|basis| {
                basis.loop_edges.contains(&selected) && basis.loop_edges.contains(&support)
            })
            .unwrap();
        let mut projection_masses = masses.clone();
        projection_masses[selected] = f(3.0);
        let mut evaluator = generate_rstar_t_dependence_evaluator(3).unwrap();
        for native_lmb in [graph.loop_momentum_basis.clone(), alternate] {
            let native_lmbs = ti_vec![native_lmb];
            let native_lmb = &native_lmbs[LmbIndex::from(0)];
            let full_subspace = SubspaceData::new_from_parent_basis_edges(
                &native_lmb.loop_edges.iter().copied().collect_vec(),
                &graph.full_filter(),
                LmbIndex::from(0),
                graph,
                &native_lmbs,
            )
            .unwrap();
            let normal = native_lmb
                .loop_edges
                .iter_enumerated()
                .find_map(|(index, edge)| (*edge == selected).then_some(index))
                .unwrap();
            let tangent = LoopIndex(1 - normal.0);
            for p in [0.7, 1.0e8] {
                let tau = f(1.2);
                let mut raw = make_sample::<T>(3.0, p, 5.0);
                // Native source perturbation is introduced before root/frame work.
                raw.sample.loop_moms[LoopIndex(0)].px += f(1.0e-25);
                let velocity = raw.loop_moms()[LoopIndex(0)].px.clone();
                raw.sample.external_moms[ExternalIndex(0)].spatial.px = f(0.7);
                raw.sample.external_moms[ExternalIndex(1)].spatial.px = f(-0.7);
                if native_lmb != &graph.loop_momentum_basis {
                    let mut origin = raw.clone();
                    origin.sample.loop_moms = LoopMomenta::from_iter(
                        (0..2).map(|_| ThreeMomentum::new(f(0.0), f(0.0), f(0.0))),
                    );
                    let native_origin =
                        origin.lmb_transform(&graph.loop_momentum_basis, native_lmb);
                    assert_ne!(
                        native_origin.loop_moms()[tangent].px,
                        f(0.0),
                        "alternate parent must exercise an actual external affine shift"
                    );
                }
                let mut point = make_point(raw.rescaled_loop_momenta(&tau, None));
                point.unrescaled_sample = raw.clone();
                point.lu_cut_parameter_cache = vec![LUParams {
                    tstar: DualOrNot::NonDual(tau.clone()),
                    h_function: DualOrNot::NonDual(f(1.0)),
                }];
                for order in 1..=3 {
                    let shape = HyperDual::new(simple_n_deriv_shape(order));
                    let t = shape.variable(0, tau.clone());
                    let mut packet = raw.clone();
                    packet.sample.dual_loop_moms =
                        Some(raw.loop_moms().rescale_with_hyper_dual(&t, None));
                    point.dualized_momentum_sample_cache.push(packet);
                }
                let native_sample = point
                    .representative_sample()
                    .lmb_transform(&graph.loop_momentum_basis, native_lmb);
                let center = LoopMomenta::from_iter((0..2).map(|index| {
                    ThreeMomentum::new(
                        F(if index == normal.0 { 2.0 } else { 1.25 }),
                        F(0.0),
                        F(0.0),
                    )
                }));
                let overlap = OverlapStructure {
                    existing_esurfaces: ti_vec![EsurfaceID(0)],
                    overlap_groups: vec![OverlapGroup {
                        existing_esurfaces: vec![ExistingEsurfaceId::from(0)],
                        complement: vec![],
                        center,
                    }],
                };
                let shared = Arc::new(ThresholdSolveGroup {
                    thresholds: vec![surface.clone()].into(),
                    subspace: full_subspace.clone(),
                    kinematics: vec![native_sample],
                    records: vec![GlobalOverlapRecord {
                        cut_group_id: CutGroupId(0),
                        side: ThresholdCountertermSide::Left,
                        local_threshold_id: 0,
                        variant_id: None,
                    }],
                    overlap: overlap.clone(),
                    prefactor_power: 4,
                });
                let shared_context = [(Arc::clone(&shared), 0)];
                let builder = CounterTermBuilder::new(
                    graph,
                    &settings,
                    &shared.thresholds,
                    point,
                    &overlap,
                    &projection_masses,
                    &native_lmbs,
                    &full_subspace,
                    None,
                    &rotation,
                    CutGroupId(0),
                    &shared_context,
                );
                let overlap_builder = builder.new_overlap_builder(0, EsurfaceID(0));
                let solution = overlap_builder
                    .new_esurface_builder(ExistingEsurfaceId::from(0))
                    .solve_rstar(
                        &mut evaluator,
                        &RadialRootIdentity::new("analytic alpha projection".into()),
                        &mut RadialRootDiagnostics::default(),
                    )
                    .unwrap();
                let t = HyperDual::new(simple_n_deriv_shape(3)).variable(0, tau.clone());
                let expected_alpha =
                    new_constant(&t, &f(2.0)) / (t.clone() * &velocity - new_constant(&t, &f(2.0)));
                for (order, (actual, expected)) in solution
                    .alpha_t_dependence
                    .as_ref()
                    .unwrap()
                    .values
                    .iter()
                    .zip(&expected_alpha.values)
                    .enumerate()
                {
                    assert!(
                        (actual - expected).abs() < f(1.0e-9),
                        "alpha Taylor order {order}, p={p}: actual={actual}, expected={expected}"
                    );
                }
                if f(1.0e-25) > tau.epsilon() * f(1000.0) {
                    let unperturbed = f(2.0) / (f(3.0) * &tau - f(2.0));
                    assert!((&solution.alpha - unperturbed).abs() > f(1.0e-27));
                    assert!((&solution.alpha - &expected_alpha.values[0]).abs() < f(1.0e-28));
                    for (actual, expected) in solution
                        .alpha_t_dependence
                        .as_ref()
                        .unwrap()
                        .values
                        .iter()
                        .zip(&expected_alpha.values)
                    {
                        assert!((actual - expected).abs() < f(1.0e-27));
                    }
                }
                let base = solution.base_rstar_loop_momenta();
                assert!((&base[normal].px - f(4.0)).abs() < f(1.0e-11));
                assert!(
                    (&solution.solution.derivative_at_solution * &overlap_builder.radius
                        - f(0.8) * (&velocity * &tau - f(2.0)))
                    .abs()
                        < f(1.0e-11)
                );
                for order in 1..=3 {
                    let geometry = solution.dual_geometry_for_order(order);
                    for (momentum, base) in geometry.rstar_loop_momenta.iter().zip(base.iter()) {
                        assert_eq!(momentum.px.values[0], base.px);
                    }
                    for derivative in &geometry.rstar_loop_momenta[normal].px.values[1..] {
                        assert!(
                            derivative.abs() < f(1.0e-8),
                            "projected normal must remain on shell"
                        );
                    }
                    let index = CutCFFIndex {
                        lu_cut_order: Some(order + 1),
                        left_threshold_order: Some(3),
                        right_threshold_order: None,
                    };
                    let variable = variable_indices_from_cut_cff_index(&index).left_threshold;
                    let mixed = solution.dual_geometry_for_cut_cff_index(&index, variable);
                    let (orders, coefficients) =
                        extract_coefficient_t_duals(&mixed.rstar_loop_momenta[normal].px, 0);
                    let radial = orders.iter().position(|order| order == &[1]).unwrap();
                    let source = &overlap_builder
                        .transformed_kinematic_point
                        .sample_for_order(order);
                    let shifted = source
                        .sample
                        .dual_loop_moms
                        .as_ref()
                        .unwrap()
                        .iter()
                        .zip(overlap_builder.center.iter())
                        .map(|(momentum, center)| {
                            momentum.clone()
                                - center.map_ref(&|value| new_constant(&momentum.px, value))
                        })
                        .collect::<LoopMomenta<_>>();
                    let expected_direction =
                        shifted[normal].px.clone() / dual_shifted_radius(&shifted, &full_subspace);
                    for (actual, expected) in coefficients[radial]
                        .values
                        .iter()
                        .zip(&expected_direction.values)
                    {
                        assert!((actual - expected).abs() < f(1.0e-10));
                    }
                    // eta_r=E'(4+delta_r*D(t))*D(t), with E(q)=sqrt(q^2+9).
                    // Its first three radial Taylor coefficients are known exactly.
                    let (radial_orders, radial_coefficients) =
                        extract_coefficient_t_duals(&mixed.esurface_derivative, 0);
                    for (power, factor) in
                        [f(4.0) / f(5.0), f(9.0) / f(125.0), -f(108.0) / f(6250.0)]
                            .into_iter()
                            .enumerate()
                    {
                        let slot = radial_orders
                            .iter()
                            .position(|order| order == &[power])
                            .unwrap();
                        let expected = (0..power).fold(expected_direction.clone(), |value, _| {
                            value * &expected_direction
                        }) * factor;
                        for (actual, expected) in radial_coefficients[slot]
                            .values
                            .iter()
                            .zip(&expected.values)
                        {
                            assert!(
                                (actual - expected).abs() < f(1.0e-8),
                                "mixed eta_r order {power}: actual={actual}, expected={expected}"
                            );
                        }
                    }
                    assert_eq!(
                        mixed.rstar_loop_momenta[tangent].px.values[0],
                        base[tangent].px
                    );
                    assert_eq!(mixed.radius_star.values[0], geometry.radius_star.values[0]);
                }
            }
        }
    }

    check_alpha::<f64>(&graph, selected, support, &thresholds[EsurfaceID(0)]);
    check_alpha::<crate::utils::f128>(&graph, selected, support, &thresholds[EsurfaceID(0)]);
    check_alpha::<crate::utils::ArbPrec>(&graph, selected, support, &thresholds[EsurfaceID(0)]);

    // Exercise the complete raised iterated evaluation, including the recording branch: a
    // second multiplication by the channel weight after the helper changes this result.
    // This is an algebraic diagnostic with a constant residue, not a generated cross section.
    let right_edge = lmbs[LmbIndex::from(0)].loop_edges[LoopIndex(1)];
    let right_subspace = SubspaceData::new_from_parent_basis_edges(
        &[right_edge],
        &graph.full_filter(),
        LmbIndex::from(0),
        &graph,
        &lmbs,
    )
    .unwrap();
    let right_threshold = Esurface {
        energies: vec![right_edge],
        ..thresholds[EsurfaceID(0)].clone()
    };
    let mut runtime_point = make_point(make_sample::<f64>(3.0, 17.0, 4.0));
    runtime_point.lu_cut_parameter_cache = vec![LUParams {
        tstar: DualOrNot::NonDual(F(1.0)),
        h_function: DualOrNot::NonDual(F(1.0)),
    }];
    runtime_point.lu_cut_esurface_values = vec![DualOrNot::NonDual(F(1.0))];
    let raised_group = Arc::new(ThresholdSolveGroup {
        thresholds: group.thresholds.clone(),
        subspace: group.subspace.clone(),
        kinematics: group.kinematics.clone(),
        records: group.records.clone(),
        overlap: group.overlap.clone(),
        prefactor_power: 4,
    });
    let right_group = Arc::new(ThresholdSolveGroup {
        thresholds: vec![right_threshold.clone()].into(),
        subspace: right_subspace.clone(),
        kinematics: vec![runtime_point.representative_sample().clone()],
        records: vec![GlobalOverlapRecord {
            cut_group_id: CutGroupId(0),
            side: ThresholdCountertermSide::Right,
            local_threshold_id: 0,
            variant_id: None,
        }],
        overlap: OverlapStructure {
            existing_esurfaces: ti_vec![EsurfaceID(0)],
            overlap_groups: vec![OverlapGroup {
                existing_esurfaces: vec![ExistingEsurfaceId::from(0)],
                complement: vec![],
                center: group.overlap.overlap_groups[1].center.clone(),
            }],
        },
        prefactor_power: 2,
    });
    let shared = LUSharedOverlaps {
        left: OverlapStructure {
            existing_esurfaces: ti_vec![EsurfaceID(0)],
            overlap_groups: local.overlap_groups[..2].to_vec(),
        },
        right: right_group.overlap.clone(),
        left_groups: vec![(Arc::clone(&raised_group), 0), (raised_group, 1)],
        right_groups: vec![(right_group, 0)],
    };
    let index = CutCFFIndex {
        left_threshold_order: Some(2),
        right_threshold_order: Some(1),
        lu_cut_order: Some(1),
    };
    let resolved = ResolvedThresholdCounterterms {
        legacy_equivalent: true,
        variants: [
            (ThresholdCountertermSide::Left, &subspaces[1], 2),
            (ThresholdCountertermSide::Right, &right_subspace, 1),
        ]
        .into_iter()
        .map(
            |(side, subspace, order)| ResolvedThresholdCountertermVariant {
                name: format!("{side:?}"),
                group_id: None,
                cut_group_id: Some(CutGroupId(0)),
                associations: vec![],
                side,
                threshold_esurface_ids: vec![EsurfaceID(0)],
                raised_esurface_group: RaisedEsurfaceGroup {
                    esurface_ids: vec![EsurfaceID(0)],
                    max_occurence: order,
                },
                requested_subspace: None,
                requested_parent_lmb: None,
                resolved_parent_lmb: subspace.get_lmb(&lmbs).loop_edges.clone().into(),
                subspace: subspace.clone(),
                subspace_loop_count: 1,
                multiplier: None,
            },
        )
        .collect(),
        cross_section_cut_groups: ti_vec![ResolvedCutGroupThresholdCounterterms {
            left: vec![ThresholdCountertermVariantId(0)],
            right: vec![ThresholdCountertermVariantId(1)],
        }],
    };
    let metadata_registry = ThresholdCountertermMetadataRegistry::build(
        &graph.name,
        &resolved,
        &lmbs,
        &[ThresholdCountertermVariantStatus {
            generated: true,
            active: true,
        }; 2],
        vec![],
    )
    .unwrap();
    let mut param_builder = graph.param_builder.clone();
    param_builder.initialize_duals(2);
    let evaluator_settings = GlobalSettings::default().generation.evaluator;
    let orientation = graph.underlying.new_edgevec(|_, _, _| Orientation::Default);
    let (residue, _) = EvaluatorStack::new_with_timings(
        &[Atom::num(1)],
        &param_builder,
        std::slice::from_ref(&orientation),
        &[OrientationID(0)],
        shape_from_cut_cff_index(&index),
        &evaluator_settings,
    )
    .unwrap();
    let helper = CrossSectionGraph::new(graph.clone())
        .iterated_th_helper(
            2,
            1,
            1,
            1,
            true,
            LUThresholdHelperOutputs::LegacyAndPieces,
            None,
            evaluator_settings.optimization_settings(),
            &evaluator_settings,
        )
        .unwrap();
    let mut counterterm = LUCounterTerm {
        evaluators: ti_vec![LUCounterTermEvaluators {
            left_thresholds_evaluator: ti_vec![BTreeMap::new()],
            right_thresholds_evaluator: ti_vec![BTreeMap::new()],
            iterated_evaluator: IteratedCtCollection::new(
                vec![BTreeMap::from([(index, residue)])],
                1,
                1,
            ),
            threshold_helpers: LUThresholdHelperEvaluators {
                left_thresholds: ti_vec![BTreeMap::new()],
                right_thresholds: ti_vec![BTreeMap::new()],
                iterated: IteratedCtCollection::new(vec![BTreeMap::from([(index, helper)])], 1, 1,),
            },
            threshold_multipliers: None,
            residue_from_e_surface_evaluators: vec![build_derivative_structure(
                1,
                -1,
                &evaluator_settings,
            )],
        }],
        thresholds: ti_vec![(
            ti_vec![thresholds[EsurfaceID(0)].clone()],
            ti_vec![right_threshold]
        )],
        subspaces: ti_vec![(subspaces[1].clone(), right_subspace)],
        variant_subspaces: None,
        metadata_registry: Some(metadata_registry),
        active_cut_groups: ti_vec![true],
        active_left_thresholds: ti_vec![ti_vec![true]],
        active_right_thresholds: ti_vec![ti_vec![true]],
        active_iterated_thresholds: ti_vec![IteratedCtCollection::new(vec![true], 1, 1)],
        rstar_dependence_calculator: ti_vec![generate_rstar_t_dependence_evaluator(0).unwrap()],
    };
    let [unrecorded, recorded] = [false, true].map(|record_components| {
        counterterm
            .evaluate(
                &runtime_point,
                CutGroupId(0),
                &[],
                &lmbs,
                &graph,
                &masses,
                &rotation,
                &settings,
                &mut param_builder,
                SingleOrAllOrientations::Single {
                    orientation: &orientation,
                    id: OrientationID(0),
                },
                &mut EvaluationMetaData::new_empty(),
                false,
                record_components,
                Some(&shared),
            )
            .unwrap()
    });
    assert!(unrecorded.components.is_none());
    assert!(unrecorded.total.re.0.is_finite() && unrecorded.total.im.0.is_finite());
    let scale = unrecorded.total.norm_squared().0.sqrt();
    assert!(scale > 0.0);
    assert!((recorded.total - unrecorded.total).norm_squared().0.sqrt() < 1e-12 * scale);
    let components = recorded.components.unwrap();
    assert_eq!(components.len(), 8);
    assert!(components.iter().all(|component| matches!(
        &component.occurrence,
        ThresholdCountertermComponentOccurrence::LocalUnitarity {
            left_threshold_order: Some(2),
            right_threshold_order: Some(1),
            ..
        }
    )));
    let component_sum = components
        .iter()
        .fold(Complex::new_re(F(0.0)), |sum, component| {
            sum + component.weighted
        });
    assert!((component_sum - recorded.total).norm_squared().0.sqrt() < 1e-12 * scale);
}
