use super::utils::*;
use super::*;
use std::fs;

use gammalooprs::integrands::process::ProcessIntegrand;
use gammalooprs::settings::runtime::{
    RotationSetting, StabilityLevelSetting, StabilityRecordingSettings,
};

#[test]
#[serial]
fn raised_cut_numerator_cancels_one_propagator_in_both_orientation_modes() -> Result<()> {
    for explicit_orientation_sum_only in [false, true] {
        let mode = if explicit_orientation_sum_only {
            "explicit"
        } else {
            "localized"
        };
        let test_root =
            get_tests_workspace_path().join(format!("raised_cut_numerator_cancellation_{mode}"));
        let mut cli = get_test_cli(None, &test_root, Some(mode.to_string()), true)?;
        run_commands(
            &mut cli,
            &[
                "import model ./assets/models/json/scalars/scalars_2p_3p.json",
                "import graphs ./tests/resources/graphs/raised_cut_numerator_cancellation.dot -p raised_cut_cancellation -i compare",
                &format!(
                    "set global kv global.generation.explicit_orientation_sum_only={explicit_orientation_sum_only} global.generation.tropical_subgraph_table.disable_tropical_generation=true global.generation.evaluator.iterative_orientation_optimization=false global.generation.evaluator.compile=false global.generation.evaluator.store_atom=true global.generation.threshold_subtraction.enable_thresholds=false"
                ),
                r#"set default-runtime string '
                    [general]
                    integral_unit = "none"
                    disable_flux_factor = true

                    [kinematics.externals]
                    type = "constant"

                    [kinematics.externals.data]
                    momenta = [[4.0, 0.0, 0.0, 0.0]]
                    helicities = [1]

                    [sampling]
                    graphs = "summed"
                    orientations = "summed"
                    lmb_multichanneling = true
                    lmb_channel_weight = "ose"
                    lmb_channels = "summed"

                    [subtraction]
                    disable_threshold_subtraction = true
                '"#,
                "generate existing -p raised_cut_cancellation -i compare",
                "select -p raised_cut_cancellation -i compare --with-only-graph-names powered_cancel --output_process raised_cut_powered --output_integrand compare",
                "select -p raised_cut_cancellation -i compare --with-only-graph-names lower_bubble --output_process raised_cut_lower --output_integrand compare",
                "generate existing -p raised_cut_powered -i compare",
                "generate existing -p raised_cut_lower -i compare",
            ],
        )?;

        for point in [[0.11, 0.23, 0.37], [0.71, 0.43, 0.19]] {
            let powered =
                inspect_xspace_process(&mut cli, "raised_cut_powered", "compare", &point)?;
            let lower = inspect_xspace_process(&mut cli, "raised_cut_lower", "compare", &point)?;
            let combined =
                inspect_xspace_process(&mut cli, "raised_cut_cancellation", "compare", &point)?;
            let scale = powered.re.hypot(powered.im).max(lower.re.hypot(lower.im));
            let direct_difference = powered + lower;

            assert!(
                scale > 0.0,
                "the raised-cut cancellation oracle is trivial at {point:?} in {mode} mode"
            );
            assert!(
                direct_difference.re.hypot(direct_difference.im) <= 1.0e-12 * scale,
                "the q^2-m^2 numerator did not cancel one raised propagator at {point:?} in {mode} mode: powered={powered:e}, lower={lower:e}"
            );
            assert!(
                combined.re.hypot(combined.im) <= 1.0e-12 * scale,
                "the summed LU graph did not preserve the raised-propagator cancellation at {point:?} in {mode} mode: combined={combined:e}, scale={scale:e}"
            );
        }

        clean_test(test_root);
    }
    Ok(())
}

#[test]
#[serial]
fn uv_profile_selects_lu_graph_and_cut() -> Result<()> {
    let test_root = get_tests_workspace_path().join("uv_profile_selects_lu_graph_and_cut");
    let mut cli = get_test_cli(None, &test_root, Some("uv_profile_cut".to_string()), true)?;
    run_commands(
        &mut cli,
        &[
            "import model ./assets/models/json/scalars/scalars_2p_3p.json",
            "import graphs ./tests/resources/graphs/mass_approach_scalar_self_energy.dot -p scalar_self_energy -i compare",
            "set global kv global.generation.tropical_subgraph_table.disable_tropical_generation=true global.generation.evaluator.iterative_orientation_optimization=false global.generation.evaluator.compile=false global.generation.evaluator.store_atom=true global.generation.threshold_subtraction.enable_thresholds=false",
            r#"set default-runtime string '
                [general]
                integral_unit = "none"
                disable_flux_factor = true

                [kinematics.externals]
                type = "constant"

                [kinematics.externals.data]
                momenta = [[4.0, 0.0, 0.0, 0.0]]
                helicities = [1]

                [sampling]
                graphs = "summed"
                orientations = "summed"
                lmb_multichanneling = true
                lmb_channel_weight = "ose"
                lmb_channels = "summed"

                [subtraction]
                disable_threshold_subtraction = true
            '"#,
            "generate existing -p scalar_self_energy -i compare",
        ],
    )?;

    let cut_edges = {
        let process = ProcessRef::Unqualified("scalar_self_energy".to_string());
        let integrand_name = "compare".to_string();
        let (process_id, integrand_name) = cli
            .state
            .find_integrand_ref(Some(&process), Some(&integrand_name))?;
        let ProcessIntegrand::CrossSection(integrand) = cli
            .state
            .process_list
            .get_integrand(process_id, &integrand_name)?
            .require_generated()?
        else {
            panic!("scalar self-energy fixture must generate a cross-section integrand")
        };
        let graph = integrand
            .data
            .graph_terms
            .iter()
            .find(|term| term.graph.name == "dotted_bubble")
            .expect("scalar self-energy fixture must contain dotted_bubble");
        let cut = graph
            .cuts
            .first()
            .expect("dotted_bubble must contain a physical Cutkosky cut");
        graph
            .graph
            .underlying
            .iter_edges_of(&cut.cut)
            .map(|(_, edge_id, _)| usize::from(edge_id))
            .sorted()
            .collect_vec()
    };
    assert!(!cut_edges.is_empty());

    let analysis = Profile::UltraViolet(UltraVioletProfile {
        process: Some(ProcessRef::Unqualified("scalar_self_energy".to_string())),
        integrand_name: Some("compare".to_string()),
        graph: Some("dotted_bubble".to_string()),
        cutkosky_cut: cut_edges.clone(),
        n_points: 6,
        ..Default::default()
    })
    .run(&mut cli.state, &cli.cli_settings)?
    .unwrap_uv();

    assert_eq!(analysis.graphs.len(), 1);
    let graph = &analysis.graphs[0];
    assert_eq!(graph.graph_name, "dotted_bubble");
    assert_eq!(
        graph
            .cutkosky_cut
            .as_ref()
            .expect("selected cut identity must be preserved")
            .iter()
            .copied()
            .map(usize::from)
            .collect_vec(),
        cut_edges
    );
    let fitted_limits = analysis
        .graphs
        .iter()
        .flat_map(|graph| &graph.lmbs)
        .flat_map(|lmb| &lmb.subsets)
        .filter(|subset| subset.estimated_dod().is_some())
        .count();
    let report = analysis.pass_fail(-0.9);
    assert!(
        report.total > 0,
        "the selected UV profile must be nonvacuous"
    );
    assert!(
        fitted_limits > 0,
        "the selected UV profile must fit at least one nonvanishing limit"
    );
    assert_eq!(report.failed, 0, "{report}");

    clean_test(test_root);
    Ok(())
}

#[test]
fn inspect_x_space_reports_invalid_coordinate_count_cleanly() -> Result<()> {
    let mut cli = get_test_cli(
        Some("z_decay_test.toml".into()),
        get_tests_workspace_path().join("z_decay_invalid_inspect_point"),
        None,
        true,
    )?;

    let error = Inspect {
        process: Some(ProcessRef::Id(0)),
        integrand_name: Some("default".to_string()),
        point: vec![0.1, 0.2],
        momentum_space: false,
        discrete_dim: vec![0, 0],
        ..Default::default()
    }
    .run(&mut cli)
    .expect_err("inspect should reject an invalid x-space point length");

    let rendered = format!("{error:#}");
    assert!(
        rendered.contains("Expected 3 x-space coordinates for this integrand selection, got 2."),
        "{rendered}"
    );

    clean_test(&cli.cli_settings.state.folder);
    Ok(())
}

#[test]
fn bench_cli_profiles_fixed_scalar_triangle_point_and_restores_settings() -> Result<()> {
    let mut cli = setup_scalar_topologies_cli("bench_cli_fixed_scalar_triangle")?;
    let point = default_xspace_point_for(&cli, "triangle", "scalar_tri")?;
    let point_arg = point.iter().map(|entry| format!("{entry:.17e}")).join(" ");
    let normal_json_path = cli.cli_settings.state.folder.join("inspect_normal.json");
    let bench_json_path = cli.cli_settings.state.folder.join("bench.json");

    cli.run_command(&format!(
        "inspect -p triangle -i scalar_tri -x {point_arg} --json-output {}",
        normal_json_path.display()
    ))?;
    let normal_json: serde_json::Value = serde_json::from_slice(&fs::read(&normal_json_path)?)?;
    assert!(normal_json.get("evaluation").is_some(), "{normal_json:#}");

    let process_ref = ProcessRef::Unqualified("triangle".to_string());
    let integrand_ref = "scalar_tri".to_string();
    let (process_id, integrand_name) =
        cli.find_integrand_ref(Some(&process_ref), Some(&integrand_ref))?;
    let original_settings = {
        let integrand = cli
            .process_list
            .get_integrand_mut(process_id, &integrand_name)?;
        let settings = integrand.get_mut_settings();
        settings.general.enable_cache = true;
        settings.general.debug_cache = true;
        settings.general.generate_events = true;
        settings.general.store_additional_weights_in_event = true;
        settings.stability.rotation_axis = vec![RotationSetting::Pi2X {}];
        settings.stability.levels = vec![
            StabilityLevelSetting::default_double(),
            StabilityLevelSetting::default_quad(),
            StabilityLevelSetting::default_arb(),
        ];
        settings.stability.recording = Some(StabilityRecordingSettings {
            record_rotated_results: true,
            record_all_stability_levels: true,
            record_loop_momenta_escalation: true,
        });
        settings.clone()
    };

    cli.run_command(&format!(
        "bench -p triangle -i scalar_tri -x {point_arg} --duration 0.001 --n-batches 2 --minimal-integrand --json-output {}",
        bench_json_path.display()
    ))?;
    let bench_json: serde_json::Value = serde_json::from_slice(&fs::read(&bench_json_path)?)?;
    assert_eq!(bench_json["n_batches"], 2);
    assert_eq!(bench_json["minimal_integrand"], true);
    let summary = bench_json["summary"]
        .as_array()
        .expect("bench JSON should contain summary rows");
    assert!(
        summary
            .iter()
            .any(|row| row["category"].as_str() == Some("Total")),
        "{bench_json:#}"
    );
    let restored_settings = cli
        .process_list
        .get_integrand_mut(process_id, &integrand_name)?
        .get_settings()
        .clone();
    assert_eq!(restored_settings, original_settings);

    let value = inspect_xspace_process(&mut cli, "triangle", "scalar_tri", &point)?;
    assert!(value.re.is_finite() && value.im.is_finite());

    clean_test(&cli.cli_settings.state.folder);
    Ok(())
}

#[test]
fn integrate_writes_numerical_stability_histograms_for_scalar_triangle() -> Result<()> {
    let mut cli = setup_scalar_topologies_cli("integrate_numerical_stability_histograms")?;
    cli.run_command(
        "set process -p triangle -i scalar_tri kv integrator.n_start=2 integrator.min_samples_for_update=2 integrator.n_max=2 integrator.n_increase=0",
    )?;

    let workspace = get_tests_workspace_path()
        .join("integrate_numerical_stability_histograms")
        .join("integration_workspace");
    Integrate {
        process: vec![ProcessRef::Unqualified("triangle".to_string())],
        integrand_name: vec!["scalar_tri".to_string()],
        workspace_path: Some(workspace.clone()),
        target: vec![],
        n_cores: Some(1),
        restart: true,
        ..Default::default()
    }
    .run(&mut cli.state, &cli.cli_settings)?;

    let stability_dir = workspace.join("numerical_stability");
    assert!(stability_dir.join("global.json").exists());
    assert!(stability_dir.join("global.hwu").exists());
    assert!(stability_dir.join("triangle@scalar_tri.json").exists());
    assert!(stability_dir.join("triangle@scalar_tri.hwu").exists());

    let global_bundle = gammalooprs::observables::ObservableSnapshotBundle::from_json_file(
        stability_dir.join("global.json"),
    )?;
    for key in ["Double", "Quad", "ArbPrec"] {
        assert!(
            global_bundle.histograms.contains_key(key),
            "missing numerical stability histogram {key}"
        );
    }

    clean_test(&cli.cli_settings.state.folder);
    Ok(())
}

#[test]
#[serial]
fn inspect_x_space_reports_missing_discrete_dimensions_cleanly() -> Result<()> {
    let mut cli =
        setup_gg_hhh_threshold_amplitude_cli("gg_hhh_inspect_missing_discrete_dimensions")?;

    let error = Inspect {
        process: Some(ProcessRef::Unqualified("gg_hhh".to_string())),
        integrand_name: Some("1L".to_string()),
        point: vec![0.1, 0.2, 0.3],
        momentum_space: false,
        ..Default::default()
    }
    .run(&mut cli)
    .expect_err("inspect should reject x-space points that omit required discrete dimensions");

    let rendered = format!("{error:#}");
    assert!(
        rendered.contains("This integrand uses discrete graph sampling")
            && rendered.contains("sample_orientations = true")
            && rendered.contains("sampling_type = discrete_multi_channeling")
            && rendered.contains(
                "requires 3 discrete dimensions [graph group, orientation, channel], but got 0."
            ),
        "{rendered}"
    );

    clean_test(&cli.cli_settings.state.folder);
    Ok(())
}

#[test]
fn inspect_momentum_space_graph_id_is_sampling_agnostic() -> Result<()> {
    let mut cli = setup_scalar_topologies_cli("scalar_box_graph_id_inspect")?;

    cli.run_command(
        r#"set process -p box -i scalar_box string '
[sampling]
graphs = "monte_carlo"
orientations = "summed"
lmb_multichanneling = false
lmb_channels = "summed"
'"#,
    )?;

    let (_, graph_selected) = Inspect {
        process: Some(ProcessRef::Unqualified("box".to_string())),
        integrand_name: Some("scalar_box".to_string()),
        point: vec![0.1, 0.2, 0.3],
        momentum_space: true,
        graph_id: Some(0),
        ..Default::default()
    }
    .run(&mut cli)?;

    let (_, discrete_selected) = Inspect {
        process: Some(ProcessRef::Unqualified("box".to_string())),
        integrand_name: Some("scalar_box".to_string()),
        point: vec![0.1, 0.2, 0.3],
        momentum_space: true,
        discrete_dim: vec![0],
        ..Default::default()
    }
    .run(&mut cli)?;

    assert_eq!(graph_selected, discrete_selected);

    clean_test(&cli.cli_settings.state.folder);
    Ok(())
}

#[test]
fn test_mass_approach_scalar_self_energy() -> Result<()> {
    let mut cli = get_test_cli(
        Some("mass_approach_scalar_self_energy.toml".into()),
        get_tests_workspace_path().join("mass_approach_scalar_self_energy"),
        None,
        true,
    )?;

    let mass_values = [1.1, 1.05, 1.01, 1.001, 1.0001];
    let mut inspect_magnitudes = Vec::with_capacity(mass_values.len());

    for mass in mass_values {
        cli.run_command(&format!("set model mass_scalar_2={mass}"))?;

        let (_, inspect) = Inspect {
            process: None,
            integrand_name: Some("default".to_string()),
            point: vec![0.1, 0.2, 0.3, 0.3, 0.4, 0.5],
            momentum_space: false,
            use_arb_prec: true,
            ..Default::default()
        }
        .run(&mut cli.state)?;

        let magnitude = (inspect.re * inspect.re + inspect.im * inspect.im).sqrt();
        inspect_magnitudes.push(magnitude);
    }

    assert!(
        inspect_magnitudes.windows(2).all(|pair| pair[1] <= pair[0]),
        "Inspect magnitude is not monotonically decreasing as mass_scalar_2 approaches 1: {inspect_magnitudes:?}"
    );
    assert!(
        inspect_magnitudes.last().copied().unwrap_or(f64::INFINITY) < 2.0e-10,
        "Inspect did not approach zero closely enough near mass_scalar_2=1: {inspect_magnitudes:?}"
    );

    Ok(())
}

#[test]
fn test_mass_approach_threshold_subtraction() -> Result<()> {
    let mut cli = get_test_cli(
        Some("generate_threshold_subtraction_mass_approach.toml".into()),
        get_tests_workspace_path().join("threshold_subtraction_mass_approach"),
        None,
        true,
    )?;

    let mass_values = [2.1, 2.05, 2.01, 2.001, 2.0001];
    let mut inspect_magnitudes = Vec::with_capacity(mass_values.len());

    for mass in mass_values {
        cli.run_command(&format!("set model mass_scalar_2={mass}"))?;

        let (_, inspect) = Inspect {
            process: None,
            integrand_name: Some("default".to_string()),
            point: vec![0.1, 0.2, 0.3, 0.3, 0.4, 0.5],
            momentum_space: false,
            use_arb_prec: true,
            ..Default::default()
        }
        .run(&mut cli.state)?;

        let magnitude = (inspect.re * inspect.re + inspect.im * inspect.im).sqrt();
        inspect_magnitudes.push(magnitude);
    }

    assert!(
        inspect_magnitudes.windows(2).all(|pair| pair[1] <= pair[0]),
        "Inspect magnitude is not monotonically decreasing as mass_scalar_2 approaches 2: {inspect_magnitudes:?}"
    );
    assert!(
        inspect_magnitudes.last().copied().unwrap_or(f64::INFINITY) < 4.0e-10,
        "Inspect did not approach zero closely enough near mass_scalar_2=2  : {inspect_magnitudes:?}"
    );

    Ok(())
}

#[test]
fn test_mass_approach_threshold_subtraction_reversed() -> Result<()> {
    let mut cli = get_test_cli(
        Some("generate_threshold_subtraction_mass_approach_rev.toml".into()),
        get_tests_workspace_path().join("threshold_subtraction_mass_approach_rev"),
        None,
        true,
    )?;

    let mass_values = [2.1, 2.05, 2.01, 2.001, 2.0001];
    let mut inspect_magnitudes = Vec::with_capacity(mass_values.len());

    for mass in mass_values {
        cli.run_command(&format!("set model mass_scalar_2={mass}"))?;

        let (_, inspect) = Inspect {
            process: None,
            integrand_name: Some("default".to_string()),
            point: vec![0.1, 0.2, 0.3, 0.3, 0.4, 0.5],
            momentum_space: false,
            use_arb_prec: true,
            ..Default::default()
        }
        .run(&mut cli)?;

        let magnitude = (inspect.re * inspect.re + inspect.im * inspect.im).sqrt();
        inspect_magnitudes.push(magnitude);
    }

    assert!(
        inspect_magnitudes.windows(2).all(|pair| pair[1] <= pair[0]),
        "Inspect magnitude is not monotonically decreasing as mass_scalar_2 approaches 2: {inspect_magnitudes:?}"
    );
    assert!(
        inspect_magnitudes.last().copied().unwrap_or(f64::INFINITY) < 4.0e-10,
        "Inspect did not approach zero closely enough near mass_scalar_2=2: {inspect_magnitudes:?}"
    );

    Ok(())
}

#[derive(Debug)]
struct OrientationMassApproachOutcome {
    orientation_id: usize,
    pattern: String,
    magnitudes: Vec<f64>,
    monotonic: bool,
    all_small: bool,
    final_small: bool,
}

impl OrientationMassApproachOutcome {
    fn passed(&self) -> bool {
        (self.monotonic || self.all_small) && self.final_small
    }
}

fn orientation_pattern_from_signature(signature: &[i8]) -> Result<String> {
    let entries = signature
        .iter()
        .map(|sign| match sign {
            1 => Ok("+"),
            -1 => Ok("-"),
            0 => Ok("0"),
            other => Err(eyre::eyre!(
                "Unexpected orientation signature entry {other}; expected one of -1, 0, 1"
            )),
        })
        .collect::<Result<Vec<_>>>()?;

    Ok(format!("({})", entries.join(",")))
}

fn inspect_magnitudes_for_mass_values(
    cli: &mut gammaloop_integration_tests::CLIState,
    mass_values: &[f64],
) -> Result<Vec<f64>> {
    let mut inspect_magnitudes = Vec::with_capacity(mass_values.len());

    for &mass in mass_values {
        cli.run_command(&format!("set model mass_scalar_2={mass}"))?;

        let (_, inspect) = Inspect {
            process: None,
            integrand_name: Some("default".to_string()),
            point: vec![0.1, 0.2, 0.3, 0.3, 0.4, 0.5],
            momentum_space: false,
            use_arb_prec: true,
            ..Default::default()
        }
        .run(&mut cli.state)?;

        let magnitude = (inspect.re * inspect.re + inspect.im * inspect.im).sqrt();
        inspect_magnitudes.push(magnitude);
    }

    Ok(inspect_magnitudes)
}

#[test]
fn test_mass_approach_threshold_subtraction_dotted() -> Result<()> {
    let mut cli = get_test_cli(
        Some("generate_threshold_subtraction_mass_approach_dotted.toml".into()),
        get_tests_workspace_path().join("threshold_subtraction_mass_approach_dotted"),
        None,
        true,
    )?;

    let mass_values = [2.1, 2.05, 2.01, 2.001, 2.0001];
    let info = cli.state.get_integrand_info(None, None)?;
    let canonical_group = info
        .graph_groups
        .first()
        .ok_or_else(|| eyre::eyre!("Generated integrand has no graph groups to classify"))?;
    let canonical_signature_set = canonical_group
        .orientations
        .iter()
        .map(|orientation| orientation.signature.clone())
        .collect::<std::collections::BTreeSet<_>>();

    for group in info.graph_groups.iter().skip(1) {
        let signature_set = group
            .orientations
            .iter()
            .map(|orientation| orientation.signature.clone())
            .collect::<std::collections::BTreeSet<_>>();
        if signature_set != canonical_signature_set {
            return Err(eyre::eyre!(
                "Per-orientation dotted threshold test requires all graph groups to share the same orientation signatures; group {} has {:?}, canonical group {} has {:?}",
                group.group_id,
                signature_set,
                canonical_group.group_id,
                canonical_signature_set,
            ));
        }
    }

    let mut outcomes = Vec::with_capacity(canonical_group.orientations.len());
    let small_threshold = 4.0e-10;
    for orientation in &canonical_group.orientations {
        let pattern = orientation_pattern_from_signature(&orientation.signature)?;
        cli.run_command(&format!(
            "set process -p {} -i {} string '\n[general.orientation_pat]\npat = \"{}\"\n'",
            info.process_id, info.integrand_name, pattern
        ))?;

        let magnitudes = inspect_magnitudes_for_mass_values(&mut cli, &mass_values)?;
        let monotonic = magnitudes.windows(2).all(|pair| pair[1] <= pair[0]);
        let all_small = magnitudes
            .iter()
            .all(|magnitude| *magnitude < small_threshold);
        let final_small = magnitudes.last().copied().unwrap_or(f64::INFINITY) < small_threshold;

        outcomes.push(OrientationMassApproachOutcome {
            orientation_id: orientation.orientation_id,
            pattern,
            magnitudes,
            monotonic,
            all_small,
            final_small,
        });
    }

    let passed = outcomes
        .iter()
        .filter(|outcome| outcome.passed())
        .map(|outcome| format!("#{} {}", outcome.orientation_id, outcome.pattern))
        .collect_vec();
    let failed = outcomes
        .iter()
        .filter(|outcome| !outcome.passed())
        .collect_vec();

    let mut report = String::new();
    writeln!(
        &mut report,
        "Per-orientation threshold-subtraction mass-approach summary for {}:{}",
        info.process_name, info.integrand_name
    )
    .unwrap();
    writeln!(
        &mut report,
        "Passed orientations: {}",
        if passed.is_empty() {
            "<none>".to_string()
        } else {
            passed.join(", ")
        }
    )
    .unwrap();
    writeln!(
        &mut report,
        "Failed orientations: {}",
        if failed.is_empty() {
            "<none>".to_string()
        } else {
            failed
                .iter()
                .map(|outcome| format!("#{} {}", outcome.orientation_id, outcome.pattern))
                .join(", ")
        }
    )
    .unwrap();
    for outcome in &failed {
        writeln!(
            &mut report,
            "  #{} {} monotonic={} all_small={} final_small={} magnitudes={:?}",
            outcome.orientation_id,
            outcome.pattern,
            outcome.monotonic,
            outcome.all_small,
            outcome.final_small,
            outcome.magnitudes,
        )
        .unwrap();
    }

    println!("{report}");
    assert!(failed.is_empty(), "{report}");

    Ok(())
}

mod failing {
    use super::*;

    #[test]
    fn scalar_sunrise_inspect() -> Result<()> {
        symbolica::GLOBAL_SETTINGS
            .initialize_tracing
            .store(false, std::sync::atomic::Ordering::Relaxed);
        let mut cli = get_test_cli(
            Some("scalar_sunrise.toml".into()),
            get_tests_workspace_path().join("scalar_sunrise"),
            Some("scalar_sunrise".to_string()),
            false,
        )?;

        let point = [1., 1., 1., 2., 3., 4.];

        let point = [1., 1., 1., -3., -4., -5.];

        let point = vec![2., 3., 4., 1., 1., 1.];
        let mut ins = Inspect {
            point: point.clone(),
            momentum_space: true,
            ..Default::default()
        };

        // from Kaapo: m=1 muv=5 4.37688e-03 m=2 muv=5 	2.48100e-03	 m=3 muv=5 1.07231e-03
        cli.run_command("set model mass_scalar_1=1.0")?;

        fn string_with_prefactor(rs: &[Complex<f64>]) -> String {
            let mut out = String::new();
            let prefactor = -(2. * std::f64::consts::PI).powi(6);
            for r in rs {
                let re = r.re * prefactor;
                let im = r.im * prefactor;
                writeln!(&mut out, "{re:.5e}+i{im:.5e}").unwrap();
            }
            out
        }

        let (jac, rall_1) = ins.run(&mut cli.state)?;
        ins.point = ins.point.iter().map(|a| a * 10.).collect_vec();
        let (jac, rall_10) = ins.run(&mut cli.state)?;
        ins.point = ins.point.iter().map(|a| a * 10.).collect_vec();
        let (jac, rall_100) = ins.run(&mut cli.state)?;
        ins.point = point.clone();
        cli.run_command("set process -p 0 -i default kv general.additional_param_values=[1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]")?;
        let (jac, r1_1) = ins.run(&mut cli.state)?;
        ins.point = ins.point.iter().map(|a| a * 10.).collect_vec();
        let (jac, r1_10) = ins.run(&mut cli.state)?;
        ins.point = ins.point.iter().map(|a| a * 10.).collect_vec();
        let (jac, r1_100) = ins.run(&mut cli.state)?;
        ins.point = point.clone();
        cli.run_command("set process -p 0 -i default kv general.additional_param_values=[0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0]")?;
        let (jac, r2_1) = ins.run(&mut cli.state)?;
        ins.point = ins.point.iter().map(|a| a * 10.).collect_vec();
        let (jac, r2_10) = ins.run(&mut cli.state)?;
        ins.point = ins.point.iter().map(|a| a * 10.).collect_vec();
        let (jac, r2_100) = ins.run(&mut cli.state)?;
        ins.point = point.clone();
        cli.run_command("set process -p 0 -i default kv general.additional_param_values=[0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0]")?;
        let (jac, r3_1) = ins.run(&mut cli.state)?;
        ins.point = ins.point.iter().map(|a| a * 10.).collect_vec();
        let (jac, r3_10) = ins.run(&mut cli.state)?;
        ins.point = ins.point.iter().map(|a| a * 10.).collect_vec();
        let (jac, r3_100) = ins.run(&mut cli.state)?;
        ins.point = point.clone();
        cli.run_command("set process -p 0 -i default kv general.additional_param_values=[0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0]")?;
        let (jac, r4_1) = ins.run(&mut cli.state)?;
        ins.point = ins.point.iter().map(|a| a * 10.).collect_vec();
        let (jac, r4_10) = ins.run(&mut cli.state)?;
        ins.point = ins.point.iter().map(|a| a * 10.).collect_vec();
        let (jac, r4_100) = ins.run(&mut cli.state)?;
        ins.point = point.clone();
        cli.run_command("set process -p 0 -i default kv general.additional_param_values=[0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0]")?;
        let (jac, r5_1) = ins.run(&mut cli.state)?;
        ins.point = ins.point.iter().map(|a| a * 10.).collect_vec();
        let (jac, r5_10) = ins.run(&mut cli.state)?;
        ins.point = ins.point.iter().map(|a| a * 10.).collect_vec();
        let (jac, r5_100) = ins.run(&mut cli.state)?;
        ins.point = point.clone();
        cli.run_command("set process -p 0 -i default kv general.additional_param_values=[0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0]")?;
        let (jac, r6_1) = ins.run(&mut cli.state)?;
        ins.point = ins.point.iter().map(|a| a * 10.).collect_vec();
        let (jac, r6_10) = ins.run(&mut cli.state)?;
        ins.point = ins.point.iter().map(|a| a * 10.).collect_vec();
        let (jac, r6_100) = ins.run(&mut cli.state)?;

        ins.point = point.clone();
        cli.run_command("set process -p 0 -i default kv general.additional_param_values=[0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0]")?;
        let (jac, r7_1) = ins.run(&mut cli.state)?;
        ins.point = ins.point.iter().map(|a| a * 10.).collect_vec();
        let (jac, r7_10) = ins.run(&mut cli.state)?;
        ins.point = ins.point.iter().map(|a| a * 10.).collect_vec();
        let (jac, r7_100) = ins.run(&mut cli.state)?;
        ins.point = point.clone();
        cli.run_command("set process -p 0 -i default kv general.additional_param_values=[0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0]")?;
        let (jac, r8_1) = ins.run(&mut cli.state)?;
        ins.point = ins.point.iter().map(|a| a * 10.).collect_vec();
        let (jac, r8_10) = ins.run(&mut cli.state)?;
        ins.point = ins.point.iter().map(|a| a * 10.).collect_vec();
        let (jac, r8_100) = ins.run(&mut cli.state)?;

        insta::assert_snapshot!(string_with_prefactor(&[r1_1,r2_1,r3_1,r4_1,r5_1,r6_1,r7_1,r8_1,rall_1]),@"
        2.18603e-4+i-0.00000e0
        -4.41097e-5+i-0.00000e0
        -1.54032e-4+i-0.00000e0
        -1.57503e-4+i-0.00000e0
        -7.17631e-5+i-0.00000e0
        4.21936e-5+i-0.00000e0
        1.40322e-4+i-0.00000e0
        8.50437e-5+i-0.00000e0
        5.87544e-5+i-0.00000e0
        ");
        insta::assert_snapshot!(string_with_prefactor(&[r1_10,r2_10,r3_10,r4_10,r5_10,r6_10,r7_10,r8_10,rall_10]),@"
        2.66555e-8+i-0.00000e0
        -1.11736e-8+i-0.00000e0
        -3.96106e-7+i-0.00000e0
        -4.55447e-8+i-0.00000e0
        -2.65802e-8+i-0.00000e0
        1.11735e-8+i-0.00000e0
        3.96096e-7+i-0.00000e0
        4.54492e-8+i-0.00000e0
        -3.03929e-11+i-0.00000e0
        ");
        insta::assert_snapshot!(string_with_prefactor(&[r1_100,r2_100,r3_100,r4_100,r5_100,r6_100,r7_100,r8_100,rall_100]),@"
        2.67150e-12+i-0.00000e0
        -1.13180e-12+i-0.00000e0
        -4.46155e-11+i-0.00000e0
        -4.62050e-12+i-0.00000e0
        -2.67150e-12+i-0.00000e0
        1.13180e-12+i-0.00000e0
        4.46155e-11+i-0.00000e0
        4.62050e-12+i-0.00000e0
        -3.63173e-19+i-0.00000e0
        ");
        // clean_test(&cli.cli_settings.state.folder);

        Ok(())
    }

    #[test]
    fn test_epem_tth_inspect_nlo_gl18() -> Result<()> {
        let mut cli = get_test_cli(
            Some("test_epem_tth_inspect_nlo_gl18.toml".into()),
            get_tests_workspace_path().join("test_epem_tth_inspect_nlo_gl18"),
            None,
            true,
        )
        .unwrap();

        let (_, inspect) = Inspect {
            process: None,
            integrand_name: None,
            point: vec![0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.823, 0.214],
            momentum_space: false,
            ..Default::default()
        }
        .run(&mut cli)
        .unwrap();

        let target = Complex::new(-9.487984855932107e-6, 3.610476200052732e-5);
        assert_eq!(inspect, target);
        Ok(())
    }

    #[test]
    fn test_qqx_aaa_ir_tree_unprocessed_inspect() -> Result<()> {
        let mut cli = get_test_cli(
            Some("generate_qqx_aaa_tree_unprocessed.toml".into()),
            get_tests_workspace_path().join("qqx_aaa_tree_unprocessed"),
            None,
            true,
        )
        .unwrap();

        let (_, inspect) = Inspect {
            process: None,
            integrand_name: None,
            point: vec![],
            momentum_space: false,
            ..Default::default()
        }
        .run(&mut cli)
        .unwrap();

        let target = Complex::new(0.00014727604164105595, -0.001150313936913021);
        assert_eq!(inspect, target);
        Ok(())
    }

    #[test]
    fn test_qqx_aaa_ir_tree_user_numerator_unprocessed_with_momtrop_table_inspect() -> Result<()> {
        let mut cli = get_test_cli(
            Some("generate_qqx_aaa_tree_user_numerator_unprocessed_with_momtrop_table.toml".into()),
            get_tests_workspace_path()
                .join("qqx_aaa_tree_user_numerator_unprocessed_with_momtrop_table"),
            None,
            true,
        )
        .unwrap();

        let (_, inspect) = Inspect {
            process: None,
            integrand_name: None,
            point: vec![],
            momentum_space: false,
            ..Default::default()
        }
        .run(&mut cli)
        .unwrap();

        let target = Complex::new(1.47276041641056e-4, -1.1503139369130214e-3);
        assert_eq!(inspect, target);
        Ok(())
    }

    #[test]
    fn test_qqx_aaa_ir_tree_user_numerator_inspect() -> Result<()> {
        let mut cli = get_test_cli(
            Some("generate_qqx_aaa_tree_user_numerator.toml".into()),
            get_tests_workspace_path().join("qqx_aaa_tree_user_numerator"),
            None,
            true,
        )
        .unwrap();

        let (_, inspect) = Inspect {
            process: None,
            integrand_name: None,
            point: vec![],
            momentum_space: false,
            ..Default::default()
        }
        .run(&mut cli)
        .unwrap();

        let target = Complex::new(1.47276041641056e-4, -1.1503139369130214e-3);
        assert_eq!(inspect, target);
        Ok(())
    }
}

mod slow {
    use super::*;

    #[test]
    fn test_qqx_aaa_ir_subtracted_inspect() -> Result<()> {
        let mut cli = get_test_cli(
            Some("generate_qqx_aaa_ir_subtracted_physical.toml".into()),
            get_tests_workspace_path().join("qqx_aaa_ir_subtracted_physical"),
            None,
            true,
        )
        .unwrap();

        let (_, inspect) = Inspect {
            process: None,
            integrand_name: None,
            point: vec![0.0, 0.10001, 0.10001],
            momentum_space: true,
            ..Default::default()
        }
        .run(&mut cli)
        .unwrap();

        let target = Complex::new(2.3159767780905335e-1, -1.8547720156633686e-4);
        assert_eq!(inspect, target);
        Ok(())
    }
}
