use std::fs;

use super::utils::*;
use super::*;

fn setup_explicit_sum_scalar_topologies_cli(
    test_name: &str,
) -> Result<gammaloop_integration_tests::CLIState> {
    setup_explicit_sum_scalar_topologies_cli_with_representation(test_name, "CFF", true, true)
}

fn setup_explicit_sum_scalar_topologies_cli_with_representation(
    test_name: &str,
    representation: &str,
    subtract_uv: bool,
    enable_thresholds: bool,
) -> Result<gammaloop_integration_tests::CLIState> {
    let mut cli = get_test_cli(
        None,
        get_tests_workspace_path().join(test_name),
        Some(test_name.to_string()),
        true,
    )?;

    run_commands(
        &mut cli,
        &[
            "import model scalars-default.json",
            "remove processes",
            &format!(
                "set global kv global.3d_representation={representation} global.generation.explicit_orientation_sum_only=true global.generation.evaluator.compile=false global.generation.uv.subtract_uv={subtract_uv} global.generation.threshold_subtraction.enable_thresholds={enable_thresholds}"
            ),
            "set default-runtime kv general.evaluator_method=Summed",
            "generate amp scalar_1 > scalar_0 scalar_0 [{1}] --allowed-vertex-interactions V_3_SCALAR_022 V_3_SCALAR_122 -p triangle -i scalar_tri",
            "generate amp scalar_0 scalar_0 > scalar_0 scalar_0 [{1}] --allowed-vertex-interactions V_3_SCALAR_022 -p box -i scalar_box --select-graphs GL0",
            "generate",
            "set model mass_scalar_2=2.0",
            "set model mass_scalar_1=1.0",
            SCALAR_TRIANGLE_EXTERNALS,
            SCALAR_BOX_ABOVE_EXTERNALS,
        ],
    )?;

    Ok(cli)
}

struct ImportedExternalTreeCase {
    graph_file: &'static str,
    momenta: &'static str,
    helicities: &'static str,
    x_space_point: Vec<f64>,
}

fn setup_imported_external_tree_cli(
    test_name: &str,
    case: &ImportedExternalTreeCase,
    explicit_orientation_sum_only: bool,
) -> Result<gammaloop_integration_tests::CLIState> {
    setup_imported_external_tree_cli_with_representation(
        test_name,
        case,
        explicit_orientation_sum_only,
        "CFF",
        true,
    )
}

fn setup_imported_external_tree_cli_with_representation(
    test_name: &str,
    case: &ImportedExternalTreeCase,
    explicit_orientation_sum_only: bool,
    representation: &str,
    subtract_uv: bool,
) -> Result<gammaloop_integration_tests::CLIState> {
    let mut cli = get_test_cli(
        None,
        get_tests_workspace_path().join(test_name),
        Some(test_name.to_string()),
        true,
    )?;

    run_commands(
        &mut cli,
        &[
            "import model scalars-default.json",
            "remove processes",
            &format!(
                "set global kv global.3d_representation={representation} global.generation.evaluator.compile=false global.generation.evaluator.summed=true global.generation.uv.subtract_uv={subtract_uv} global.generation.threshold_subtraction.enable_thresholds=false"
            ),
            r#"set default-runtime string '
[general]
evaluator_method = "Summed"

[subtraction]
disable_threshold_subtraction = true

[sampling]
graphs = "summed"
orientations = "summed"
lmb_multichanneling = false
lmb_channels = "summed"
coordinate_system = "spherical"
mapping = "linear"
'"#,
        ],
    )?;

    if explicit_orientation_sum_only {
        cli.run_command("set global kv global.generation.explicit_orientation_sum_only=true")?;
    }

    cli.run_command(&format!(
        r#"set default-runtime string '
[kinematics.externals]
type = "constant"

[kinematics.externals.data]
momenta = [
{}
]
helicities = {}
'"#,
        case.momenta, case.helicities
    ))?;
    cli.run_command(&format!(
        "import graphs ./tests/resources/graphs/{}",
        case.graph_file
    ))?;
    cli.run_command("generate")?;

    Ok(cli)
}

fn inspect_single_imported_process(
    cli: &mut gammaloop_integration_tests::CLIState,
    point: &[f64],
) -> Result<Complex<f64>> {
    let (_, inspect) = Inspect {
        process: None,
        integrand_name: None,
        point: point.to_vec(),
        momentum_space: false,
        ..Default::default()
    }
    .run(cli)?;
    Ok(inspect)
}

fn threedreps_dot_path(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("resources/graphs/threedreps")
        .join(name)
}

fn write_threedreps_dot_with_edge_numerator(
    test_name: &str,
    dot_name: &str,
    numerator: &str,
) -> Result<PathBuf> {
    let test_root = get_tests_workspace_path().join(test_name);
    clean_test(&test_root);
    fs::create_dir_all(&test_root)?;
    let source = fs::read_to_string(threedreps_dot_path(dot_name))?;
    let escaped = numerator.replace('\\', "\\\\").replace('"', "\\\"");
    let updated = source.replacen(
        "particle=\"scalar_1\"",
        &format!("num=\"{escaped}\" particle=\"scalar_1\""),
        1,
    );
    let dot_path = test_root.join(format!(
        "{}_num.dot",
        dot_name.strip_suffix(".dot").unwrap_or(dot_name)
    ));
    fs::write(&dot_path, updated)?;
    Ok(dot_path)
}

fn imported_q_dot_q(lhs_edge_id: usize, rhs_edge_id: usize) -> String {
    format!(
        "(Q({lhs_edge_id},spenso::cind(0))*Q({rhs_edge_id},spenso::cind(0))-Q({lhs_edge_id},spenso::cind(1))*Q({rhs_edge_id},spenso::cind(1))-Q({lhs_edge_id},spenso::cind(2))*Q({rhs_edge_id},spenso::cind(2))-Q({lhs_edge_id},spenso::cind(3))*Q({rhs_edge_id},spenso::cind(3)))"
    )
}

struct HighPowerImportedScalarCase {
    dot_name: &'static str,
    process_name: &'static str,
    numerator: String,
    momenta: &'static str,
    helicities: &'static str,
}

fn setup_high_power_imported_scalar_cli_with_representation(
    test_name: &str,
    case: &HighPowerImportedScalarCase,
    representation: &str,
) -> Result<gammaloop_integration_tests::CLIState> {
    let dot_path =
        write_threedreps_dot_with_edge_numerator(test_name, case.dot_name, &case.numerator)?;
    let state_path = get_tests_workspace_path().join(test_name).join("state");
    let mut cli = get_test_cli(None, &state_path, Some(test_name.to_string()), true)?;

    run_commands(
        &mut cli,
        &[
            "import model scalars-default.json",
            "remove processes",
            &format!(
                "set global kv global.3d_representation={representation} global.generation.explicit_orientation_sum_only=true global.generation.evaluator.compile=false global.generation.uv.subtract_uv=false global.generation.threshold_subtraction.enable_thresholds=false"
            ),
            &format!(
                r#"set default-runtime string '
[general]
evaluator_method = "Summed"
enable_cache = false
debug_cache = false
generate_events = true
store_additional_weights_in_event = true
integral_unit = "picobarn"

[subtraction]
disable_threshold_subtraction = true

[sampling]
graphs = "summed"
orientations = "summed"
lmb_multichanneling = false
lmb_channels = "summed"
coordinate_system = "spherical"
mapping = "linear"

[kinematics.externals]
type = "constant"

[kinematics.externals.data]
momenta = [
{}
]
helicities = {}
'"#,
                case.momenta, case.helicities
            ),
            &format!(
                "import graphs {} -p {} -i default -o",
                dot_path.display(),
                case.process_name
            ),
            "generate",
        ],
    )?;

    Ok(cli)
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
#[serial]
fn explicit_orientation_sum_inspect_matches_standard_cff() -> Result<()> {
    let mut standard = setup_scalar_topologies_cli("scalar_explicit_orientation_standard")?;
    let mut explicit = setup_explicit_sum_scalar_topologies_cli("scalar_explicit_orientation_sum")?;

    for (process, integrand, point) in [
        ("triangle", "scalar_tri", vec![0.1, 0.2, 0.3]),
        ("box", "scalar_box", vec![0.1, 0.2, 0.3]),
    ] {
        let standard_value = inspect_xspace_process(&mut standard, process, integrand, &point)?;
        let explicit_value = inspect_xspace_process(&mut explicit, process, integrand, &point)?;
        assert_complex_approx_eq(
            explicit_value,
            standard_value,
            &format!("explicit orientation sum inspect for {process}/{integrand}"),
        );
    }

    clean_test(&standard.cli_settings.state.folder);
    clean_test(&explicit.cli_settings.state.folder);
    Ok(())
}

#[test]
#[serial]
fn ltd_bare_explicit_orientation_sum_inspect_matches_cff() -> Result<()> {
    let mut cff = setup_explicit_sum_scalar_topologies_cli_with_representation(
        "scalar_ltd_bare_cff_reference",
        "CFF",
        false,
        false,
    )?;
    let mut ltd = setup_explicit_sum_scalar_topologies_cli_with_representation(
        "scalar_ltd_bare_ltd",
        "LTD",
        false,
        false,
    )?;

    for (process, integrand, point) in [
        ("triangle", "scalar_tri", vec![0.1, 0.2, 0.3]),
        ("box", "scalar_box", vec![0.1, 0.2, 0.3]),
    ] {
        let cff_value = inspect_xspace_process(&mut cff, process, integrand, &point)?;
        let ltd_value = inspect_xspace_process(&mut ltd, process, integrand, &point)?;
        assert_complex_approx_eq(
            ltd_value,
            cff_value,
            &format!("bare LTD explicit orientation sum inspect for {process}/{integrand}"),
        );
    }

    clean_test(&cff.cli_settings.state.folder);
    clean_test(&ltd.cli_settings.state.folder);
    Ok(())
}

#[test]
#[serial]
fn ltd_threshold_explicit_orientation_sum_inspect_matches_cff() -> Result<()> {
    let mut cff = setup_explicit_sum_scalar_topologies_cli_with_representation(
        "scalar_ltd_threshold_cff_reference",
        "CFF",
        false,
        true,
    )?;
    let mut ltd = setup_explicit_sum_scalar_topologies_cli_with_representation(
        "scalar_ltd_threshold_ltd",
        "LTD",
        false,
        true,
    )?;

    for (process, integrand, point) in [
        ("triangle", "scalar_tri", vec![0.1, 0.2, 0.3]),
        ("box", "scalar_box", vec![0.1, 0.2, 0.3]),
    ] {
        let cff_value = inspect_xspace_process(&mut cff, process, integrand, &point)?;
        let ltd_value = inspect_xspace_process(&mut ltd, process, integrand, &point)?;
        assert_complex_approx_eq(
            ltd_value,
            cff_value,
            &format!("threshold-subtracted LTD inspect for {process}/{integrand}"),
        );
    }

    clean_test(&cff.cli_settings.state.folder);
    clean_test(&ltd.cli_settings.state.folder);
    Ok(())
}

#[test]
#[serial]
fn ltd_with_uv_subtraction_errors_cleanly() -> Result<()> {
    let test_name = "scalar_ltd_uv_subtraction_not_supported";
    let mut cli = get_test_cli(
        None,
        get_tests_workspace_path().join(test_name),
        Some(test_name.to_string()),
        true,
    )?;

    run_commands(
        &mut cli,
        &[
            "import model scalars-default.json",
            "remove processes",
            "set global kv global.3d_representation=LTD global.generation.explicit_orientation_sum_only=true global.generation.evaluator.compile=false global.generation.threshold_subtraction.enable_thresholds=false",
        ],
    )?;

    let error = cli
        .run_command("generate amp scalar_1 > scalar_0 scalar_0 [{1}] --allowed-vertex-interactions V_3_SCALAR_022 V_3_SCALAR_122 -p triangle -i scalar_tri")
        .expect_err("LTD should reject local UV subtraction from 3D expansions");
    let rendered = format!("{error:#}");
    assert!(
        rendered.contains("with local UV counterterms from 3D expansions is not supported"),
        "{rendered}"
    );

    clean_test(&cli.cli_settings.state.folder);
    Ok(())
}

#[test]
#[serial]
fn explicit_orientation_sum_rejects_individual_orientation_runtime_modes() -> Result<()> {
    let mut cli =
        setup_explicit_sum_scalar_topologies_cli("scalar_explicit_orientation_runtime_rejections")?;

    cli.run_command(
        "set process -p box -i scalar_box kv general.evaluator_method=SingleParametric",
    )?;
    let error = Inspect {
        process: Some(ProcessRef::Unqualified("box".to_string())),
        integrand_name: Some("scalar_box".to_string()),
        point: vec![0.1, 0.2, 0.3],
        momentum_space: false,
        ..Default::default()
    }
    .run(&mut cli)
    .expect_err("explicit orientation sum should reject non-Summed runtime evaluators");
    let rendered = format!("{error:#}");
    assert!(
        rendered.contains("requires runtime `general.evaluator_method = Summed`"),
        "{rendered}"
    );

    cli.run_command("set process -p box -i scalar_box kv general.evaluator_method=Summed")?;
    cli.run_command(
        r#"set process -p box -i scalar_box string '
[sampling]
graphs = "monte_carlo"
orientations = "monte_carlo"
lmb_multichanneling = false
lmb_channels = "summed"
coordinate_system = "spherical"
mapping = "linear"
'"#,
    )?;
    let error = Inspect {
        process: Some(ProcessRef::Unqualified("box".to_string())),
        integrand_name: Some("scalar_box".to_string()),
        point: vec![0.1, 0.2, 0.3],
        momentum_space: false,
        discrete_dim: vec![0, 0],
        ..Default::default()
    }
    .run(&mut cli)
    .expect_err("explicit orientation sum should reject orientation Monte Carlo");
    let rendered = format!("{error:#}");
    assert!(
        rendered.contains("does not support runtime individual-orientation Monte Carlo sampling"),
        "{rendered}"
    );

    cli.run_command(
        r#"set process -p box -i scalar_box string '
[general]
evaluator_method = "Summed"

[sampling]
graphs = "summed"
orientations = "summed"
lmb_multichanneling = false
lmb_channels = "summed"
coordinate_system = "spherical"
mapping = "linear"
'"#,
    )?;
    let error = Inspect {
        process: Some(ProcessRef::Unqualified("box".to_string())),
        integrand_name: Some("scalar_box".to_string()),
        point: vec![0.1, 0.2, 0.3],
        momentum_space: true,
        graph_id: Some(0),
        orientation_id: Some(0),
        ..Default::default()
    }
    .run(&mut cli)
    .expect_err("explicit orientation sum should reject explicit orientation selection");
    let rendered = format!("{error:#}");
    assert!(
        rendered.contains("explicit orientation selection is not supported"),
        "{rendered}"
    );

    clean_test(&cli.cli_settings.state.folder);
    Ok(())
}

#[test]
#[serial]
fn imported_external_tree_inspects_match_explicit_orientation_sum() -> Result<()> {
    let cases = [
        ImportedExternalTreeCase {
            graph_file: "scalar_triangle_external_tree.dot",
            momenta: r#"    [3.0, 0.0, 0.0, 3.0],
    [3.0, 0.0, 3.0, 0.0],
    [3.0, 0.0, 0.0, -3.0],
    "dependent""#,
            helicities: "[0, 0, 0, 0]",
            x_space_point: vec![0.23, 0.41, 0.67],
        },
        ImportedExternalTreeCase {
            graph_file: "scalar_box_two_external_trees.dot",
            momenta: r#"    [4.0, 0.0, 0.0, 4.0],
    [2.0, 2.0, 0.0, 0.0],
    [4.0, 0.0, 0.0, -4.0],
    [2.0, -2.0, 0.0, 0.0],
    [2.0, 0.0, 2.0, 0.0],
    "dependent""#,
            helicities: "[0, 0, 0, 0, 0, 0]",
            x_space_point: vec![0.19, 0.37, 0.61],
        },
        ImportedExternalTreeCase {
            graph_file: "scalar_pure_tree_four_point.dot",
            momenta: r#"    [3.0, 0.0, 0.0, 3.0],
    [3.0, 0.0, 0.0, -3.0],
    [3.0, 0.0, 3.0, 0.0],
    "dependent""#,
            helicities: "[0, 0, 0, 0]",
            x_space_point: Vec::new(),
        },
    ];

    for case in cases {
        let stem = case.graph_file.trim_end_matches(".dot");
        let mut standard = setup_imported_external_tree_cli(
            &format!("{stem}_standard_external_tree_inspect"),
            &case,
            false,
        )?;
        let mut explicit = setup_imported_external_tree_cli(
            &format!("{stem}_explicit_external_tree_inspect"),
            &case,
            true,
        )?;

        let standard_value = inspect_single_imported_process(&mut standard, &case.x_space_point)?;
        let explicit_value = inspect_single_imported_process(&mut explicit, &case.x_space_point)?;
        assert_complex_approx_eq(
            explicit_value,
            standard_value,
            &format!("imported external-tree inspect for {}", case.graph_file),
        );

        clean_test(&standard.cli_settings.state.folder);
        clean_test(&explicit.cli_settings.state.folder);
    }

    Ok(())
}

#[test]
#[serial]
fn ltd_bare_external_tree_inspects_match_cff() -> Result<()> {
    let cases = [
        ImportedExternalTreeCase {
            graph_file: "scalar_triangle_external_tree.dot",
            momenta: r#"    [3.0, 0.0, 0.0, 3.0],
    [3.0, 0.0, 3.0, 0.0],
    [3.0, 0.0, 0.0, -3.0],
    "dependent""#,
            helicities: "[0, 0, 0, 0]",
            x_space_point: vec![0.23, 0.41, 0.67],
        },
        ImportedExternalTreeCase {
            graph_file: "scalar_box_two_external_trees.dot",
            momenta: r#"    [4.0, 0.0, 0.0, 4.0],
    [2.0, 2.0, 0.0, 0.0],
    [4.0, 0.0, 0.0, -4.0],
    [2.0, -2.0, 0.0, 0.0],
    [2.0, 0.0, 2.0, 0.0],
    "dependent""#,
            helicities: "[0, 0, 0, 0, 0, 0]",
            x_space_point: vec![0.19, 0.37, 0.61],
        },
        ImportedExternalTreeCase {
            graph_file: "scalar_pure_tree_four_point.dot",
            momenta: r#"    [3.0, 0.0, 0.0, 3.0],
    [3.0, 0.0, 0.0, -3.0],
    [3.0, 0.0, 3.0, 0.0],
    "dependent""#,
            helicities: "[0, 0, 0, 0]",
            x_space_point: Vec::new(),
        },
    ];

    for case in cases {
        let stem = case.graph_file.trim_end_matches(".dot");
        let mut cff = setup_imported_external_tree_cli_with_representation(
            &format!("{stem}_cff_bare_external_tree_inspect"),
            &case,
            true,
            "CFF",
            false,
        )?;
        let mut ltd = setup_imported_external_tree_cli_with_representation(
            &format!("{stem}_ltd_bare_external_tree_inspect"),
            &case,
            true,
            "LTD",
            false,
        )?;

        let cff_value = inspect_single_imported_process(&mut cff, &case.x_space_point)?;
        let ltd_value = inspect_single_imported_process(&mut ltd, &case.x_space_point)?;
        assert_complex_approx_eq(
            ltd_value,
            cff_value,
            &format!("bare LTD external-tree inspect for {}", case.graph_file),
        );

        clean_test(&cff.cli_settings.state.folder);
        clean_test(&ltd.cli_settings.state.folder);
    }

    Ok(())
}

#[test]
#[serial]
fn ltd_bare_quartic_numerator_imported_scalar_inspects_match_cff() -> Result<()> {
    let cases = [
        HighPowerImportedScalarCase {
            dot_name: "box.dot",
            process_name: "high_power_box",
            numerator: format!("{}*{}", imported_q_dot_q(4, 5), imported_q_dot_q(6, 7)),
            momenta: r#"    [3.0, 0.0, 0.0, 3.0],
    [3.0, 0.0, 0.0, -3.0],
    [3.0, 0.0, 3.0, 0.0],
    "dependent""#,
            helicities: "[0, 0, 0, 0]",
        },
        HighPowerImportedScalarCase {
            dot_name: "box_pow3.dot",
            process_name: "high_power_box_pow3",
            numerator: format!("{}*{}", imported_q_dot_q(4, 6), imported_q_dot_q(7, 9)),
            momenta: r#"    [3.0, 0.0, 0.0, 3.0],
    [3.0, 0.0, 0.0, -3.0],
    [3.0, 0.0, 3.0, 0.0],
    "dependent""#,
            helicities: "[0, 0, 0, 0]",
        },
    ];

    for case in cases {
        let stem = case.dot_name.trim_end_matches(".dot");
        let cff_name = format!("{stem}_quartic_ltd_bare_cff_reference");
        let ltd_name = format!("{stem}_quartic_ltd_bare_ltd");
        let mut cff =
            setup_high_power_imported_scalar_cli_with_representation(&cff_name, &case, "CFF")?;
        let mut ltd =
            setup_high_power_imported_scalar_cli_with_representation(&ltd_name, &case, "LTD")?;

        let point = default_xspace_point_for(&cff, case.process_name, "default")?;
        let cff_value = inspect_xspace_process(&mut cff, case.process_name, "default", &point)?;
        let ltd_value = inspect_xspace_process(&mut ltd, case.process_name, "default", &point)?;
        assert_complex_approx_eq(
            ltd_value,
            cff_value,
            &format!(
                "bare LTD imported {} inspect with quartic numerator",
                case.dot_name
            ),
        );

        clean_test(get_tests_workspace_path().join(cff_name));
        clean_test(get_tests_workspace_path().join(ltd_name));
    }

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
        .run(&mut cli)?;

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
        .run(&mut cli)?;

        let magnitude = (inspect.re * inspect.re + inspect.im * inspect.im).sqrt();
        inspect_magnitudes.push(magnitude);
    }

    assert!(
        inspect_magnitudes.windows(2).all(|pair| pair[1] <= pair[0]),
        "Inspect magnitude is not monotonically decreasing as mass_scalar_2 approaches 1: {inspect_magnitudes:?}"
    );
    assert!(
        inspect_magnitudes.last().copied().unwrap_or(f64::INFINITY) < 4.0e-10,
        "Inspect did not approach zero closely enough near mass_scalar_2=1: {inspect_magnitudes:?}"
    );

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
