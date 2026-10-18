//! Private native overlapping-forest consumers; public admission remains guarded.
use super::*;
use crate::{
    graph::feynman_graph::FeynmanGraph,
    initialisation::test_initialise,
    uv::{
        approx::integrated::native_mass_tests::QcdMassReference,
        orchestrator::validate_local_counterterm_schemes, settings::FinalIntegrandDimension,
    },
};

#[test]
fn native_d11_projected_wood_preserves_complete_local_histories() -> Result<()> {
    test_initialise()?;
    let fixture = QcdMassReference::new()?;
    // Reuse the two actual canonical-QCD child certificates and their retained
    // native graph. Their analytic checks are fixture prerequisites; this new
    // forest computation below deliberately disables integrated addbacks.
    let mut certified = fixture.certify_integrated_qg_children(&[6, 6, 21])?;
    assert_eq!(certified.len(), 2);
    assert!(
        certified
            .iter()
            .all(|record| record.diagram.name() == "D11")
    );
    let mut child_edges = certified
        .iter()
        .map(|record| {
            assert_eq!(record.child.dod(), 0);
            assert_eq!(
                record.child.renormalization_scheme(),
                ApproximationType::MUV
            );
            record
                .graph
                .paired_edges(record.child.subgraph())
                .into_iter()
                .map(|edge| edge.0)
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();
    child_edges.sort();
    assert_eq!(child_edges, [vec![2, 3, 5], vec![2, 4, 6]]);
    let a = child_edges[0].clone();
    let b = child_edges[1].clone();
    let a_set = a.iter().copied().collect::<BTreeSet<_>>();
    let b_set = b.iter().copied().collect::<BTreeSet<_>>();
    assert!(!a_set.is_disjoint(&b_set));
    assert!(!a_set.is_subset(&b_set) && !b_set.is_subset(&a_set));

    let retained = certified.remove(0);
    let mut graph = retained.graph;
    let model = graph.model.clone();
    let settings = GenerationSettings {
        explicit_orientation_sum_only: true,
        uv: UVgenerationSettings {
            generate_integrated: false,
            final_integrand: FinalIntegrandDimension::ThreeD,
            local_uv_cts_from_expanded_4d_integrands: true,
            orchestrator: crate::uv::UVOrchestrator::HedgePoset,
            add_marker: false,
            ..retained.settings
        },
        ..Default::default()
    };
    assert!(!settings.medium.vacuum_subtraction);
    assert_eq!(
        settings.uv.on_shell_velocity()?,
        [Atom::one(), Atom::Zero, Atom::Zero, Atom::Zero]
    );
    let cuts = CutStructure::empty(&graph);
    assert_eq!(cuts.cuts.len(), 1);
    let cut = cuts.cuts[0].clone();
    assert!(cut.union.is_empty());
    let rejection = validate_local_counterterm_schemes(&graph, &cuts, &settings.uv)
        .expect_err("public OS-over-MUV admission remains a separate guarded boundary");
    assert!(
        rejection
            .to_string()
            .contains("unsupported nested or disconnected scheme policy")
    );
    let mut forests = Wood::new(cuts, &graph, &model, &settings)?.unfold();
    let parent = graph
        .paired_edges(&graph.full_filter())
        .into_iter()
        .map(|edge| edge.0)
        .collect::<Vec<_>>();
    assert_eq!(parent, [2, 3, 4, 5, 6]);
    let expected = BTreeSet::from([
        vec![],
        vec![a.clone()],
        vec![b.clone()],
        vec![parent.clone()],
        vec![a, parent.clone()],
        vec![b, parent],
    ]);
    let mut histories = BTreeSet::new();
    for (node, _, operation) in forests.graph.iter_nodes() {
        assert!(
            !forests.graph.is_disjoint_union(node),
            "the children overlap"
        );
        let mut history = Vec::new();
        let mut cursor = node;
        while !forests.graph[cursor].key.is_empty() {
            let source = forests.source_spinney(cursor);
            history.push(
                graph
                    .paired_edges(source.filter())
                    .into_iter()
                    .map(|edge| edge.0)
                    .collect::<Vec<_>>(),
            );
            let (previous, _) = forests
                .graph
                .unique_parent(cursor)
                .expect("each connected ordered operation has one preceding prefix");
            assert_eq!(
                forests.graph[cursor].key.op_count(),
                forests.graph[previous].key.op_count() + 1
            );
            cursor = previous;
        }
        history.reverse();
        assert_eq!(operation.key.op_count(), history.len());
        assert!(histories.insert(history), "each actual history occurs once");
    }
    assert_eq!(
        histories, expected,
        "[], [A], [B], [P], [A,P], [B,P], with no [A,B]"
    );
    assert_eq!(forests.graph.n_nodes(), 6);
    assert!(forests.compute_store.entries.is_empty());
    debug_tags!(#uv, #local, #four_d, #inspect;
        stage = "native_d11_private_projected_wood_inventory",
        file.histories = ?histories, public_guard_rejected = true,
        "The existing classified wood owns all six overlapping-source histories"
    );

    // Same production CFF-source construction as full_forest_evaluator. No
    // manual forest sum, synthetic numerator or orientation list is substituted.
    let numerator = graph
        .numerator(&graph.full_filter(), &graph.empty_subgraph())?
        .state
        .expr;
    let options = graph.denominator_only_cff_3d_expression_options();
    let surfaces = graph.get_esurface_canonization(&graph.loop_momentum_basis);
    let contracted_edges = graph.paired_edges(&graph.tree_edges.subtract(&graph.initial_state_cut));
    let production = graph.generate_3d_expression_for_integrand(
        &contracted_edges,
        &surfaces,
        &options,
        Some(numerator.expression()),
    )?;
    let orientation = OrientationProjection::exact_expression(
        &production,
        &options,
        &settings.orientation_pattern,
        true,
    );
    assert!(!orientation.orientation_ids()?.is_empty());
    // This is the production private orchestrator sequence, below the public
    // scheme guard. Each projected node runs the existing canonical/full-pole
    // and immutable signed-CFF reconstruction certificates; first failure stops.
    // Retain the same production phases and actual projection context so
    // every constructed source passes its opt-in proof before CFF preparation.
    forests.integrate(&graph, &model, crate::utils::vakint()?, &settings.uv)?;
    let mut projection_context = Local4dProjectionContext::default();
    assert_eq!(projection_context.preparations.len(), 0);
    assert_eq!(projection_context.source_preparation_builds, 0);
    assert_eq!(
        projection_context.denominator_reconstruction_certificates,
        0
    );
    projection_context.certify_denominator_reconstruction = true;
    for (compatible_subset, cutset) in forests.cuts.clone() {
        let localizer = Localizer::new(&cutset, orientation);
        for node in forests.compatible_topological_order(&compatible_subset)? {
            let operation = forests.graph[node].clone();
            let computed = forests.local_3d_for_node(
                node,
                &mut graph,
                &cutset,
                localizer,
                &settings.uv,
                &mut projection_context,
            )?;
            forests
                .compute_store
                .entry(operation)
                .or_default()
                .cuts
                .insert(cutset.clone(), computed);
        }
    }
    assert!(projection_context.denominator_reconstruction_certificates > 0);
    assert_eq!(
        projection_context.denominator_reconstruction_certificates,
        projection_context.source_preparation_builds
    );
    assert!(
        projection_context.projected_component_requests
            >= projection_context.denominator_reconstruction_certificates
    );
    debug_tags!(#uv, #local, #four_d, #cff, #inspect;
        stage = "native_projected_wood_denominator_certificate_summary",
        certificates = projection_context.denominator_reconstruction_certificates,
        source_preparations = projection_context.source_preparation_builds,
        component_requests = projection_context.projected_component_requests,
        "Every actual source preparation in this fresh context has its exact denominator certificate"
    );
    assert_eq!(forests.compute_store.entries.len(), 6);
    for (node, _, operation) in forests.graph.iter_nodes() {
        let computed = forests.compute_store.require(operation)?;
        let _ = computed.local_4d(operation)?;
        let integrated = computed.integrated(operation)?;
        assert!(integrated.physical_pole()?.is_zero());
        assert!(integrated.physical_finite_counterterm()?.is_zero());
        let projected = computed.cut(operation, &cut)?;
        if node != forests.root {
            assert!(matches!(&projected.local_3d, Local3DCts::Projected4d(_)));
        }
    }
    let generated = forests.orientation_parametric_exprs(&graph, &settings.uv)?;
    assert_eq!(generated.len(), 1);
    assert_eq!(generated[0].cuts, cut);
    let integrands = &generated[0].integrands;
    assert!(integrands.iter().any(|(_, root)| !root.is_zero()));
    assert!(
        !integrands.numerators().is_empty(),
        "physical numerator families remain retained"
    );
    debug_tags!(#uv, #local, #four_d, #inspect;
        stage = "native_d11_private_projected_wood_complete",
        operations = 6, generate_integrated = false, public_admission_changed = false,
        numerical_evaluation_performed = false,
        file.roots = ?integrands.iter().collect::<Vec<_>>(),
        file.numerator_definitions = ?integrands.numerators(),
        "All six existing wood operations reconstructed and retained their factorized local CFF sum"
    );
    Ok(())
}
