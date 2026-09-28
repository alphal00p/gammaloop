use super::utils::single_slot_integral;
use super::*;
use gammaloop_api::commands::integrate::RendererOption;
use gammaloop_integration_tests::CLIState;
use gammalooprs::{graph::FeynmanGraph, processes::ProcessCollection};
use three_dimensional_reps::RepresentationMode;

fn check_polygon_residues_and_density(
    cli: &mut CLIState,
    name: &str,
    propagators: usize,
    expected_density: f64,
) -> Result<()> {
    assert_eq!(
        cli.cli_settings
            .global
            .generation
            .three_dimensional_representations,
        [RepresentationMode::Ltd],
    );
    let ProcessCollection::Amplitudes(amplitudes) = &cli.state.process_list.processes[0].collection
    else {
        panic!("expected a scalar amplitude")
    };
    let [graph] = amplitudes["scalar"].graphs.as_slice() else {
        panic!("the polygon benchmark must contain exactly one graph")
    };
    assert_eq!(graph.graph.get_loop_number(), 1);
    assert_eq!(graph.derived_data.representations.len(), 1);
    let expression = &graph.derived_data.representations[&RepresentationMode::Ltd].expression;
    assert_eq!(expression.representation, RepresentationMode::Ltd);
    // Inspect the native source rows, before the complete sum becomes a single
    // runtime orientation. No CFF expression is generated even for this check.
    assert_eq!(expression.expression.orientations.len(), propagators);

    // Independent one-loop residue theorem, Eq. (3) of arXiv:1906.06138:
    // rho(k) = -sum_i [1/(2 E_i) product_{j!=i}
    //     1/((E_i-s_i^0+s_j^0)^2-E_j^2)]/(2 pi)^3.
    // Here s_i=sum_{j<=i}p_j, s_n=0, and E_i=sqrt((k+s_i)^2+m_i^2).
    // The constants below were evaluated at 100 decimal digits directly from
    // that formula and the source data documented in each card, independently
    // of either 3D generator.
    // This also catches external-leg permutations and propagator mass shifts.
    let (_, value) = Inspect {
        process: Some(ProcessRef::Unqualified(name.to_owned())),
        integrand_name: Some("scalar".to_owned()),
        point: vec![0.17, -0.23, 0.41],
        momentum_space: true,
        graph_id: Some(0),
        ..Default::default()
    }
    .run(&mut cli.state)?;
    assert!(value.re.abs() <= 1e-12 * expected_density);
    assert!(
        (value.im - expected_density).abs() <= 1e-9 * expected_density,
        "{name}: {value:e}, residue-theorem density {expected_density:e}i",
    );
    Ok(())
}

#[test]
#[serial]
fn ltd_decagon_reproduces_published_scalar_integral() -> Result<()> {
    // Table I, topology a)* of https://arxiv.org/pdf/1906.06138:
    // I = +i 4.31638e-7 with d^4k/(2pi)^4, numerator one, m_i=0.1 i.
    let root = get_tests_workspace_path().join("ltd_decagon_paper");
    let mut cli = get_test_cli(
        Some(PathBuf::from("decagon_ltd.toml")),
        root.clone(),
        Some("ltd_decagon_paper".to_owned()),
        true,
    )?;
    run_commands(&mut cli, &["run generate"])?;
    check_polygon_residues_and_density(&mut cli, "decagon", 10, 1.588_242_516_512_272e-12)?;

    let output = Integrate {
        process: vec![ProcessRef::Unqualified("decagon".to_owned())],
        integrand_name: vec!["scalar".to_owned()],
        workspace_path: Some(root.join("integration_workspace")),
        n_cores: Some(1),
        restart: true,
        renderer: RendererOption::Tabled,
        show_max_weight_info: false,
        no_stream_iterations: true,
        no_stream_updates: true,
        ..Default::default()
    }
    .run(&mut cli.state, &cli.cli_settings)?;
    let integral = single_slot_integral(&output);
    let reference = 4.31638e-7;
    assert!(integral.result.re.0.abs() < 1e-12 * reference);
    assert!(integral.error.re.0.abs() < 1e-12 * reference);
    let value = integral.result.im.0;
    let error = integral.error.im.0;
    assert!(value.is_finite() && error.is_finite() && error > 0.0);
    assert!(error < 0.01 * reference, "{integral:?}");
    // Keep the statistical precision requirement separate from agreement with
    // the independent published integral, whose printed last digit is rounded.
    assert!(
        (value - reference).abs() < 5.0 * error + 5e-13,
        "decagon: {value:e} +/- {error:e}; published {reference:e}",
    );
    clean_test(&root);
    Ok(())
}

#[test]
#[serial]
fn ltd_triacontagon_example_has_thirty_residues() -> Result<()> {
    let root = get_tests_workspace_path().join("ltd_triacontagon_paper");
    let mut cli = get_example_cli(
        "scalar_topologies/triacontagon_ltd.toml",
        &["generate"],
        Some(root.clone()),
        Some("ltd_triacontagon_paper".to_owned()),
        true,
    )?;
    check_polygon_residues_and_density(&mut cli, "triacontagon", 30, 9.898566154759974e-26)?;
    let diagnostic_path = root.join("diagnostic-ltd.json");
    cli.run_command(&format!(
        "3drep build -p triacontagon -i scalar -g 0 --representation ltd --workspace-path {} --json-out {} --no-pretty",
        root.join("diagnostic").display(), diagnostic_path.display(),
    ))?;
    let diagnostic: JsonValue = serde_json::from_str(&std::fs::read_to_string(diagnostic_path)?)?;
    assert_eq!(diagnostic["family"], "ltd");
    assert!(diagnostic["energy_degree_bounds"].is_null());
    let expression: three_dimensional_reps::ThreeDExpression<
        three_dimensional_reps::OrientationID,
    > = serde_json::from_value(diagnostic["expression"].clone())?;
    assert_eq!(expression.orientations.len(), 30);
    clean_test(&root);
    Ok(())
}
