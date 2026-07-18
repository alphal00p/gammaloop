use super::utils::*;
use super::*;
use gammalooprs::{
    graph::GroupId,
    integrands::{HasIntegrand, process::ProcessIntegrand},
    integrate::{IntegrationWorkspaceManifest, workspace_manifest_path},
    settings::RuntimeSettings,
    utils::serde_utils::SmartSerde,
};
use std::collections::BTreeSet;
use symbolica::numerical_integration::{Grid, Sample};

fn effective_result(
    result: &gammalooprs::integrands::evaluation::EvaluationResult,
) -> Complex<f64> {
    let factor = result.parameterization_jacobian.unwrap_or(F(1.0)).0 * result.integrator_weight.0;
    Complex::new(
        result.integrand_result.re.0 * factor,
        result.integrand_result.im.0 * factor,
    )
}

fn assert_complex_close(actual: Complex<f64>, expected: Complex<f64>, context: &str) {
    let norm = |value: Complex<f64>| (value.re * value.re + value.im * value.im).sqrt();
    let scale = norm(actual).max(norm(expected)).max(1.0e-30);
    assert!(
        norm(actual - expected) <= 1.0e-10 * scale,
        "{context}: got {actual}, expected {expected}"
    );
}

fn discrete_sample(group_id: usize, weight: f64, point: &[f64]) -> Sample<F<f64>> {
    Sample::Discrete(
        F(weight),
        group_id,
        Some(Box::new(Sample::Continuous(
            F(weight),
            point.iter().copied().map(F).collect(),
        ))),
    )
}

#[test]
#[serial]
fn amplitude_runtime_graph_subset_preserves_complete_group_metadata() -> Result<()> {
    let test_name = "amplitude_runtime_graph_subset_preserves_complete_group_metadata";
    let mut cli = get_test_cli(
        Some("test_grouped_subtraction.toml".into()),
        get_tests_workspace_path().join(test_name),
        None,
        false,
    )?;
    let process_ref = ProcessRef::Id(0);
    let requested_integrand_name = "default".to_string();
    let (process_id, integrand_name) = cli
        .state
        .find_integrand_ref(Some(&process_ref), Some(&requested_integrand_name))?;
    let all_master_names = cli
        .state
        .process_list
        .get_integrand(process_id, &integrand_name)?
        .require_generated()?
        .graph_group_master_names()
        .into_iter()
        .map(str::to_string)
        .collect_vec();
    assert!(
        all_master_names.len() >= 2,
        "amplitude fixture must contain at least two graph groups"
    );
    let master_name = all_master_names
        .last()
        .expect("validated amplitude graph groups")
        .clone();
    cli.run_command(&format!(
        "set process -p 0 -i default string '\n[sampling]\ngraphs = \"monte_carlo\"\ngraph_names = [\"{master_name}\"]\n'"
    ))?;

    let source = cli
        .state
        .process_list
        .get_integrand(process_id, &integrand_name)?
        .require_generated()?
        .clone();
    let source_fingerprint = source.resume_fingerprint()?;
    let source_backend = source.active_f64_backend();
    let selected = source.clone_with_selected_graph_groups(std::slice::from_ref(&master_name))?;

    assert_eq!(source.resume_fingerprint()?, source_fingerprint);
    assert_eq!(selected.active_f64_backend(), source_backend);
    assert!(selected.graph_count() < source.graph_count());
    assert_eq!(selected.graph_group_master_names(), [master_name.as_str()]);
    assert_eq!(
        selected.get_settings().sampling.selected_graph_names(),
        std::slice::from_ref(&master_name)
    );

    let (ProcessIntegrand::Amplitude(source), ProcessIntegrand::Amplitude(selected)) =
        (&source, &selected)
    else {
        panic!("expected amplitude integrands")
    };
    let old_group_id = source
        .data
        .graph_group_structure
        .iter_enumerated()
        .find_map(|(group_id, group)| {
            group
                .into_iter()
                .find(|&graph_id| source.data.graph_terms[graph_id].graph.is_group_master)
                .filter(|&graph_id| source.data.graph_terms[graph_id].graph.name == master_name)
                .map(|_| group_id)
        })
        .expect("selected amplitude master must resolve");
    let expected_names = source.data.graph_group_structure[old_group_id]
        .into_iter()
        .map(|graph_id| source.data.graph_terms[graph_id].graph.name.clone())
        .collect_vec();
    let selected_names = selected.data.graph_group_structure[GroupId(0)]
        .into_iter()
        .map(|graph_id| selected.data.graph_terms[graph_id].graph.name.clone())
        .collect_vec();
    assert_eq!(selected_names, expected_names);
    assert_eq!(selected.data.group_derived_data.len(), 1);
    assert_eq!(
        selected.data.group_derived_data[GroupId(0)].esurface_map,
        source.data.group_derived_data[old_group_id].esurface_map
    );
    assert_eq!(
        selected.data.group_derived_data[GroupId(0)].esurface_atoms,
        source.data.group_derived_data[old_group_id].esurface_atoms
    );
    for selected_term in &selected.data.graph_terms {
        let source_term = source
            .data
            .graph_terms
            .iter()
            .find(|term| term.graph.name == selected_term.graph.name)
            .expect("selected amplitude term must come from source");
        assert_eq!(selected_term.lmbs.len(), source_term.lmbs.len());
        assert_eq!(
            selected_term.threshold_counterterm.generated_mask,
            source_term.threshold_counterterm.generated_mask
        );
        assert_eq!(
            selected_term.threshold_counterterm.active_mask,
            source_term.threshold_counterterm.active_mask
        );
    }

    clean_test(&cli.cli_settings.state.folder);
    Ok(())
}

#[test]
#[serial]
fn cross_section_runtime_graph_subset_is_compact_normalized_and_event_safe() -> Result<()> {
    let test_name = "cross_section_runtime_graph_subset_is_compact_normalized_and_event_safe";
    let mut cli = setup_sm_differential_lu_cli(test_name)?;
    set_process_incoming_helicities(&mut cli, "epem_ddxg", "default", "[1, -1]")?;
    cli.run_command(
        "set process kv general.generate_events=true general.store_additional_weights_in_event=true stability.rotation_axis=[]",
    )?;
    let (process_id, integrand_name) = cli.state.find_integrand_ref(None, None)?;
    let all_master_names = cli
        .state
        .process_list
        .get_integrand(process_id, &integrand_name)?
        .require_generated()?
        .graph_group_master_names()
        .into_iter()
        .map(str::to_string)
        .collect_vec();
    assert!(
        all_master_names.len() >= 2,
        "cross-section fixture must contain at least two graph groups"
    );
    let selected_in_input_order = vec![
        all_master_names[all_master_names.len() - 1].clone(),
        all_master_names[0].clone(),
    ];
    cli.run_command(&format!(
        "set process string '\n[sampling]\ngraphs = \"monte_carlo\"\ngraph_names = [\"{}\", \"{}\"]\n'",
        selected_in_input_order[0], selected_in_input_order[1]
    ))?;

    let source = cli
        .state
        .process_list
        .get_integrand(process_id, &integrand_name)?
        .require_generated()?
        .clone();
    let source_fingerprint = source.resume_fingerprint()?;
    let source_graph_count = source.graph_count();
    let expected_master_names = all_master_names
        .iter()
        .filter(|name| selected_in_input_order.contains(name))
        .cloned()
        .collect_vec();
    let selected = source.clone_with_selected_graph_groups(&selected_in_input_order)?;

    assert_eq!(source.resume_fingerprint()?, source_fingerprint);
    assert_eq!(source.graph_count(), source_graph_count);
    assert_eq!(selected.graph_group_master_names(), expected_master_names);
    assert_eq!(
        selected.get_settings().sampling.selected_graph_names(),
        expected_master_names
    );
    let single_reduced =
        source.clone_with_selected_graph_groups(std::slice::from_ref(&expected_master_names[0]))?;
    assert!(single_reduced.graph_count() < source.graph_count());
    assert_eq!(single_reduced.graph_group_count(), 1);
    let Grid::Discrete(selected_grid) = selected.create_grid() else {
        panic!("selected graph groups must produce a discrete graph grid")
    };
    assert_eq!(selected_grid.bins.len(), expected_master_names.len());

    let (ProcessIntegrand::CrossSection(source_xs), ProcessIntegrand::CrossSection(selected_xs)) =
        (&source, &selected)
    else {
        panic!("expected cross-section integrands")
    };
    let unknown_error = source
        .clone_with_selected_graph_groups(&["UNKNOWN_GRAPH".to_string()])
        .err()
        .expect("unknown graph name must be rejected");
    assert!(format!("{unknown_error:?}").contains("Unknown graph 'UNKNOWN_GRAPH'"));

    for (new_group_id, group) in selected_xs.data.graph_group_structure.iter_enumerated() {
        let selected_master_id = group
            .into_iter()
            .find(|&graph_id| selected_xs.data.graph_terms[graph_id].graph.is_group_master)
            .expect("selected group must contain its master");
        let selected_master = &selected_xs.data.graph_terms[selected_master_id];
        assert_eq!(selected_master.graph.group_id, Some(new_group_id));
        assert!(selected_master.graph.is_group_master);
        let old_group_id = source.resolve_group_id_by_master_name(&selected_master.graph.name)?;
        let expected_names = source_xs.data.graph_group_structure[old_group_id]
            .into_iter()
            .map(|graph_id| source_xs.data.graph_terms[graph_id].graph.name.clone())
            .collect_vec();
        let actual_names = group
            .into_iter()
            .map(|graph_id| selected_xs.data.graph_terms[graph_id].graph.name.clone())
            .collect_vec();
        assert_eq!(actual_names, expected_names);
    }
    assert_eq!(
        selected_xs.data.graph_to_group_id.len(),
        selected_xs.data.graph_terms.len()
    );
    for selected_term in &selected_xs.data.graph_terms {
        let source_term = source_xs
            .data
            .graph_terms
            .iter()
            .find(|term| term.graph.name == selected_term.graph.name)
            .expect("selected cross-section term must come from source");
        assert_eq!(selected_term.lmbs.len(), source_term.lmbs.len());
        assert_eq!(selected_term.cuts.len(), source_term.cuts.len());
        assert_eq!(
            selected_term.threshold_candidate_esurface_ids,
            source_term.threshold_candidate_esurface_ids
        );
        assert_eq!(
            selected_term.cut_threshold_associations.len(),
            source_term.cut_threshold_associations.len()
        );
        for (selected_cut, source_cut) in selected_term
            .cut_threshold_associations
            .iter()
            .zip(source_term.cut_threshold_associations.iter())
        {
            assert_eq!(selected_cut.left.len(), source_cut.left.len());
            assert_eq!(selected_cut.right.len(), source_cut.right.len());
        }
    }
    let selected_term_names = selected_xs
        .data
        .graph_terms
        .iter()
        .map(|term| term.graph.name.clone())
        .collect::<BTreeSet<_>>();

    let model = cli
        .state
        .resolve_model_for_integrand(process_id, &integrand_name)?;
    let point = gammaloop_integration_tests::default_xspace_point(&cli)?;
    let mut selected = selected;
    selected.warm_up(&model)?;
    let group_count = expected_master_names.len();
    let mut selected_channel_sum = Complex::new(0.0, 0.0);
    let mut individual_sum = Complex::new(0.0, 0.0);
    let mut individual_magnitude_sum = 0.0;
    for (group_id, master_name) in expected_master_names.iter().enumerate() {
        let selected_sample = discrete_sample(group_id, group_count as f64, &point);
        let selected_result = selected.evaluate_sample(
            &selected_sample,
            &model,
            F(group_count as f64),
            1,
            false,
            Complex::new(F(0.0), F(0.0)),
        )?;
        let selected_effective = effective_result(&selected_result);
        selected_channel_sum += selected_effective / group_count as f64;

        let mut individual =
            source.clone_with_selected_graph_groups(std::slice::from_ref(master_name))?;
        individual.warm_up(&model)?;
        let individual_sample = discrete_sample(0, 1.0, &point);
        let individual_result = individual.evaluate_sample(
            &individual_sample,
            &model,
            F(1.0),
            1,
            false,
            Complex::new(F(0.0), F(0.0)),
        )?;
        let individual_effective = effective_result(&individual_result);
        individual_magnitude_sum += individual_effective.re.abs() + individual_effective.im.abs();
        individual_sum += individual_effective;
        assert_complex_close(
            selected_effective,
            individual_effective * group_count as f64,
            "selected graph-grid inverse-PDF normalization",
        );

        let event_sum = selected_result
            .event_groups
            .iter()
            .flat_map(|group| group.iter())
            .fold(Complex::new(0.0, 0.0), |sum, event| {
                let graph_name = selected
                    .graph_name_by_id(event.cut_info.graph_id)
                    .expect("event graph id must resolve in reduced view");
                assert!(selected_term_names.contains(graph_name));
                sum + Complex::new(event.weight.re.0, event.weight.im.0)
            });
        assert_complex_close(event_sum, selected_effective, "selected event weights");
    }
    assert_complex_close(
        selected_channel_sum,
        individual_sum,
        "selected graph-grid expectation",
    );
    assert!(
        individual_magnitude_sum > 1.0e-30,
        "cross-section subset fixture must contain a nonzero selected contribution"
    );

    let persisted_master_name = &expected_master_names[0];
    cli.run_command(&format!(
        "set process string '\n[sampling]\ngraphs = \"monte_carlo\"\ngraph_names = [\"{persisted_master_name}\"]\n\n[integrator]\nn_bins = 4\nmin_samples_for_update = 1\nn_start = 8\nn_increase = 8\nn_max = 8\n'"
    ))?;
    let workspace = get_tests_workspace_path()
        .join(test_name)
        .join("subset_integration_workspace");
    let integrate = Integrate {
        integrand_name: vec![integrand_name.clone()],
        n_cores: Some(1),
        workspace_path: Some(workspace.clone()),
        restart: true,
        batch_size: Some(4),
        no_stream_updates: true,
        no_stream_iterations: true,
        ..Default::default()
    };
    integrate.run(&mut cli.state, &cli.cli_settings)?;

    let manifest: IntegrationWorkspaceManifest = IntegrationWorkspaceManifest::from_file(
        workspace_manifest_path(&workspace),
        "runtime graph-subset integration manifest",
    )?;
    assert_eq!(manifest.integrand_fingerprints, [source_fingerprint]);
    let slot_workspace = selected_slot_workspace(&cli, &workspace, None, Some(&integrand_name))?;
    let workspace_settings: RuntimeSettings = RuntimeSettings::from_file(
        slot_workspace.join("settings.toml"),
        "runtime graph-subset workspace settings",
    )?;
    assert_eq!(
        workspace_settings.sampling.selected_graph_names(),
        std::slice::from_ref(persisted_master_name)
    );

    let changed_master_name = &expected_master_names[1];
    cli.run_command(&format!(
        "set process string '\n[sampling]\ngraphs = \"monte_carlo\"\ngraph_names = [\"{changed_master_name}\"]\n'"
    ))?;
    let mut resume = integrate;
    resume.restart = false;
    resume.run(&mut cli.state, &cli.cli_settings)?;
    let resumed_source = cli
        .state
        .process_list
        .get_integrand(process_id, &integrand_name)?
        .require_generated()?;
    assert_eq!(
        resumed_source
            .get_settings()
            .sampling
            .selected_graph_names(),
        std::slice::from_ref(persisted_master_name)
    );
    assert_eq!(
        resumed_source.resume_fingerprint()?,
        manifest.integrand_fingerprints[0]
    );

    clean_test(&cli.cli_settings.state.folder);
    Ok(())
}
