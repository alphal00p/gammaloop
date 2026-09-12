use super::*;
use crate::{
    cff::{VertexSet, esurface::RaisedEsurfaceGroup},
    dot,
    graph::parse::from_dot::IntoGraph,
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
    let make_sample = |active, fixed, energy| MomentumSample {
        sample: BareMomentumSample {
            loop_moms: LoopMomenta::from_iter(
                [active, fixed].map(|px| ThreeMomentum::new(F(px), F(0.0), F(0.0))),
            ),
            dual_loop_moms: None,
            loop_mom_cache_id: 0,
            loop_mom_base_cache_id: 0,
            external_moms: [
                FourMomentum::from_args(F(energy), F(0.0), F(0.0), F(0.0)),
                FourMomentum::from_args(F(-energy), F(0.0), F(0.0), F(0.0)),
            ]
            .into_iter()
            .collect(),
            external_mom_cache_id: 0,
            external_mom_base_cache_id: 0,
            jacobian: F(1.0),
            orientation: None,
            parameterization_branch: None,
        },
    };
    let make_point = |sample: MomentumSample<f64>| LUCTKinematicPoint {
        unrescaled_sample: sample.clone(),
        dualized_momentum_sample_cache: vec![sample],
        lu_cut_parameter_cache: vec![],
        lu_cut_esurface_values: vec![],
    };
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
    let point = make_point(make_sample(3.0, 17.0, 4.0));
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
            point.clone(),
            make_point(make_sample(100.0, -2.0 * active_sign / fixed_sign, 1.0)),
            make_point(make_sample(-100.0, 0.0, 1.0)),
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
        kinematics: vec![point.clone()],
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
    let root_sample = make_sample(4.0, 17.0, 4.0)
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
            t_dependent_solution: None,
        };
        let weight = solution.non_dual_multichanneling_factor(&root_sample);
        assert!((weight.re.0 - expected).abs() < 1e-12);
        let weight = solution.dual_multichanneling_factor(&geometry);
        assert!((weight.values[0].0 - expected).abs() < 1e-12);
        assert!((weight.values[1].0 - derivative).abs() < 1e-12);
    }

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
    let mut runtime_point = make_point(make_sample(3.0, 17.0, 4.0));
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
        kinematics: vec![runtime_point.clone()],
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
