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

impl QcdMassReference {
    fn certify_projected_wood_histories(
        &self,
        internal_pdgs: &[i64],
        generate_integrated: bool,
    ) -> Result<()> {
        let (source_name, expected_children) = match internal_pdgs {
            [6, 6, 21] => ("D11", [vec![2, 3, 5], vec![2, 4, 6]]),
            [6, 21, 21] => ("D02", [vec![2, 4, 5], vec![3, 4, 6]]),
            _ => return Err(eyre!("only the two certified native QCD vertex sources")),
        };
        // These independently certified child records select the canonical native
        // source and both attachments. Their separate graph frames and symbolic-
        // master payloads are not transplanted into this computation store.
        let mut certified = self.certify_integrated_qg_children(internal_pdgs)?;
        assert_eq!(certified.len(), 2);
        assert!(
            certified
                .iter()
                .all(|record| record.diagram.name() == source_name)
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
        assert_eq!(child_edges, expected_children);
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
                generate_integrated,
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
        let wood = Wood::new(cuts, &graph, &model, &settings)?;
        assert_eq!(wood.vakint_settings.number_of_terms_in_epsilon_expansion, 3);
        let mut forests = wood.unfold();
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
        let mut node_histories = BTreeMap::new();
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
            assert!(node_histories.insert(node, history.clone()).is_none());
            assert!(histories.insert(history), "each actual history occurs once");
        }
        assert_eq!(
            histories, expected,
            "[], [A], [B], [P], [A,P], [B,P], with no [A,B]"
        );
        assert_eq!(forests.graph.n_nodes(), 6);
        assert!(forests.compute_store.entries.is_empty());
        debug_tags!(#uv, #local, #four_d, #inspect;
            stage = "native_qg_private_projected_wood_inventory",
            source = source_name, generate_integrated,
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
        let contracted_edges =
            graph.paired_edges(&graph.tree_edges.subtract(&graph.initial_state_cut));
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
        let mut completed_children = BTreeSet::new();
        for (node, _, operation) in forests.graph.iter_nodes() {
            let history = &node_histories[&node];
            let source = forests.source_spinney(node);
            let computed = forests.compute_store.require(operation)?;
            let local = computed.local_4d(operation)?;
            let integrated = computed.integrated(operation)?;
            let finite = integrated.finite_counterterm()?;
            if !generate_integrated || node == forests.root {
                assert!(integrated.physical_pole()?.is_zero());
                assert!(integrated.physical_finite_counterterm()?.is_zero());
            }

            // Retain the complete actual Laurent payload and bookkeeping factor.
            // This equality uses only this node's graph, frame and computed values;
            // no raw coefficient from a separately imported child graph is used.
            let recursive = forests.recursion_input_4d(node)?;
            if finite.is_zero() {
                // A scalar zero completion does not replace the local ranked interface.
                assert!(recursive.atom().same_value(local.atom()));
            } else {
                assert!(recursive.atom().same_value(&local.atom().try_add(&finite)?));
            }
            let inherited = local
                .active_sectors()
                .iter()
                .chain(local.recursive_completion())
                .collect::<Vec<_>>();
            let sectors = recursive.sectors().collect::<Vec<_>>();
            assert_eq!(
                sectors.len(),
                inherited.len() + usize::from(!finite.is_zero())
            );
            assert_eq!(&sectors[..inherited.len()], inherited.as_slice());
            if !finite.is_zero() {
                let completion = sectors[inherited.len()];
                assert!(completion.atom.same_value(&finite));
                assert!(completion.active_components.is_empty());
                assert_eq!(completion.frozen_lmbs(), std::slice::from_ref(&source.lmb));
            }

            // A containing operation retains the nonzero completed child's
            // localization frame, never its overlapping sibling's frame. The
            // containing operation's own finite addback may independently vanish.
            let mut inherited_frozen_sectors = 0;
            for sector in local.active_sectors() {
                assert!(sector.frozen_lmbs().len() <= 1);
                if !sector.frozen_lmbs().is_empty() {
                    assert!(generate_integrated && history.len() == 2);
                    let (previous, _) = forests
                        .graph
                        .unique_parent(node)
                        .expect("a nested history has its actual child predecessor");
                    let child_lmb = &forests.source_spinney(previous).lmb;
                    assert_eq!(sector.frozen_lmbs(), std::slice::from_ref(child_lmb));
                    for (_, _, active_lmb) in &sector.active_components {
                        assert!(
                            active_lmb
                                .loop_edges
                                .iter()
                                .all(|edge| !child_lmb.loop_edges.contains(edge)),
                            "a frozen integration is not an active loop generator again"
                        );
                    }
                    inherited_frozen_sectors += 1;
                }
            }
            if generate_integrated && history.len() == 2 {
                assert!(
                    inherited_frozen_sectors > 0,
                    "the containing operation must retain its actual finite child"
                );
            }
            if history.len() == 1 && expected_children.contains(&history[0]) {
                assert_eq!(source.dod, 0);
                assert_eq!(source.renormalization_scheme, ApproximationType::MUV);
                if generate_integrated {
                    assert!(!integrated.physical_pole()?.is_zero());
                    assert!(!finite.is_zero());
                    let epsilon_one = finite
                        .expression()
                        .derivative(GS.dim_epsilon)
                        .replace(GS.dim_epsilon)
                        .with(0);
                    assert!(
                        !epsilon_one.is_zero(),
                        "retain the child's positive-epsilon coefficient"
                    );
                    assert!(local.recursive_completion().is_empty());
                    assert_eq!(inherited_frozen_sectors, 0);
                    assert!(completed_children.insert(history[0].clone()));
                    debug_tags!(#uv, #local, #integrated, #inspect;
                        stage = "native_qg_private_projected_wood_child_completion",
                        source = source_name, file.child_edges = ?history[0],
                        file.pole = %integrated.pole()?.expression().to_canonical_string(),
                        file.finite = %finite.expression().to_canonical_string(),
                        file.epsilon_one = %epsilon_one.to_canonical_string(),
                        file.frozen_lmb = ?source.lmb,
                        "The actual child store retains its full Laurent coefficient and one frozen completion"
                    );
                }
            }
            debug_tags!(#uv, #local, #integrated, #inspect;
                stage = "native_qg_private_projected_wood_recursion_ownership",
                source = source_name, generate_integrated, file.history = ?history,
                inherited_frozen_sectors,
                own_finite_is_zero = finite.is_zero(),
                retained_sectors = inherited.len(), recursive_sectors = sectors.len(),
                "Each recursive coefficient and frozen domain is checked against its own computation store"
            );
            let projected = computed.cut(operation, &cut)?;
            if node != forests.root {
                let Local3DCts::Projected4d(projected) = &projected.local_3d else {
                    panic!("every non-root operation uses the projected route");
                };
                // Projection groups equal active/frozen bindings. Compare its
                // actual grouped order, not a zip with ungrouped Taylor sectors.
                let grouped = local.projection_sectors()?;
                if grouped.is_empty() {
                    let [zero] = projected.sectors() else {
                        panic!("an empty local projection retains one typed zero");
                    };
                    assert!(zero.frozen_factor.is_one());
                    assert!(zero.coefficient.iter().all(|(_, root)| root.is_zero()));
                } else {
                    assert_eq!(projected.sectors().len(), grouped.len());
                    for (projected, raw) in projected.sectors().iter().zip(&grouped) {
                        let expected = raw.frozen_lmbs().iter().fold(Atom::one(), |factor, lmb| {
                            factor * GS.localizing_integrand(lmb)
                        });
                        assert_eq!(projected.frozen_factor, expected);
                    }
                }
            }
        }
        if generate_integrated {
            assert_eq!(completed_children, BTreeSet::from(expected_children));
        } else {
            assert!(completed_children.is_empty());
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
            stage = "native_qg_private_projected_wood_complete",
            source = source_name, operations = 6, generate_integrated,
            public_admission_changed = false,
            numerical_evaluation_performed = false,
            file.roots = ?integrands.iter().collect::<Vec<_>>(),
            file.numerator_definitions = ?integrands.numerators(),
            "All six wood operations reconstructed with retained coefficients, domains and factorized CFF definitions"
        );
        Ok(())
    }
}

#[test]
fn native_d11_projected_wood_preserves_complete_local_histories() -> Result<()> {
    test_initialise()?;
    QcdMassReference::new()?.certify_projected_wood_histories(&[6, 6, 21], false)
}

#[test]
fn native_d02_projected_wood_preserves_complete_local_histories() -> Result<()> {
    test_initialise()?;
    QcdMassReference::new()?.certify_projected_wood_histories(&[6, 21, 21], false)
}

#[test]
fn native_d11_projected_wood_preserves_integrated_prefix_histories() -> Result<()> {
    test_initialise()?;
    QcdMassReference::new()?.certify_projected_wood_histories(&[6, 6, 21], true)
}

#[test]
fn native_d02_projected_wood_preserves_integrated_prefix_histories() -> Result<()> {
    test_initialise()?;
    QcdMassReference::new()?.certify_projected_wood_histories(&[6, 21, 21], true)
}
