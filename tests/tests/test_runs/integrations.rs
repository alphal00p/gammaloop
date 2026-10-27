use super::utils::*;
use super::*;
use gammaloop_api::commands::integrate::{IntegrationOutput, RendererOption};
use gammalooprs::integrate::workspace_result_archive_path;

/// Integrate the run card's 100k samples in 5k-sample iterations and load the
/// estimate archived after its first 20k samples.
fn integrate_with_20k_sample_archive(
    cli: &mut gammaloop_integration_tests::CLIState,
    name: &str,
) -> Result<(
    IntegrationOutput,
    RuntimeIntegrationResult,
    std::time::Duration,
)> {
    let started = std::time::Instant::now();
    let output = Integrate {
        renderer: RendererOption::Tabled,
        show_max_weight_info: false,
        no_stream_iterations: true,
        no_stream_updates: true,
        write_results_for_each_iteration: true,
        ..default_integrate_for(name)
    }
    .run(&mut cli.state, &cli.cli_settings)?;
    let elapsed = started.elapsed();
    let partial =
        load_integration_result(&workspace_result_archive_path(&output.workspace_path, 4))?;
    Ok((output, partial, elapsed))
}

#[test]
fn fermi_surface_one_loop_integration_converges_to_analytic_target() -> Result<()> {
    let name = "fermi_surface_one_loop_integration";
    let mut cli = get_test_cli(
        Some("fermi_surface_1l_integration.toml".into()),
        get_tests_workspace_path().join(name),
        Some(name.to_string()),
        true,
    )?;
    // For D=q0²-k²-m²+i0 and |mu|>m, the vacuum-subtracted double pole integrates to
    // I2=-i acosh(|mu|/m)/(8pi²). Since D^-3 = (1/2) d(D^-2)/d(m²),
    // I3=(1/2) dI2/d(m²)=i |mu|/(32pi² m² sqrt(mu²-m²)).
    // At m=3, mu=muB/3=5, pF=4 this fixes the absolute normalization and
    // phase independently of generated coefficients and h(t). Both delta
    // and delta' contribute: dropping delta' alone shifts the result by 10%.
    let expected = 5.0 / (32.0 * std::f64::consts::PI.powi(2) * 3.0_f64.powi(2) * 4.0);
    for (profile, sigma) in [("exponential", 1.0), ("poly_exponential", 1.3)] {
        cli.run_command(&format!(
            "set process kv h_function.function=\"{profile}\" h_function.sigma={sigma}"
        ))?;
        // The run card integrates 100k samples in 5k-sample iterations. The
        // estimate archived after four iterations covers the first 20k, so
        // the convergence check repeats no samples. Check decreasing
        // uncertainty, not monotonic central values: statistical fluctuations
        // need not approach the target at every iteration. Early stopping is
        // disabled in the run card.
        let (output, partial, elapsed) = integrate_with_20k_sample_archive(&mut cli, name)?;
        let mut errors = Vec::new();
        for (samples, result) in [(20_000, &partial), (100_000, &output.result)] {
            let estimate = single_slot_integral(result);
            let (value, error) = (estimate.result.im.0, estimate.error.im.0);
            info!(
                profile,
                sigma,
                samples = estimate.neval,
                value,
                error,
                expected,
                ?elapsed,
                "One-loop Fermi-surface integration convergence"
            );
            assert_eq!(estimate.neval, samples, "the sample budget must be fixed");
            assert!(
                [value, error, estimate.result.re.0, estimate.error.re.0]
                    .into_iter()
                    .all(f64::is_finite)
                    && value > 0.0
                    && error > 0.0
                    && estimate.error.re.0 >= 0.0,
                "{profile}, samples={samples}: {estimate}"
            );
            assert!(
                estimate.result.re.0.abs() < 1.0e-12 * expected.abs()
                    && estimate.error.re.0 < 1.0e-12 * expected.abs(),
                "the real part must vanish: {estimate}"
            );
            assert!(
                (value - expected).abs() <= 5.0 * error,
                "{profile}, samples={samples}: {value:e} +/- {error:e}, target={expected:e}"
            );
            errors.push(error);
        }
        assert!(
            errors[1] < 0.005 * expected.abs() && errors[1] < 0.7 * errors[0],
            "{profile}: increased statistics must reduce uncertainty below 0.5%: {errors:?}"
        );
    }
    clean_test(&cli.cli_settings.state.folder);
    Ok(())
}

#[test]
fn fermi_surface_two_loop_integration_converges_in_nonaligned_basis() -> Result<()> {
    use gammalooprs::momentum::SignOrZero;
    use gammalooprs::processes::ProcessCollection;

    let name = "fermi_surface_two_loop_integration";
    let mut cli = get_test_cli(
        Some("fermi_surface_2l_integration.toml".into()),
        get_tests_workspace_path().join(name),
        Some(name.to_string()),
        true,
    )?;
    let ProcessCollection::Amplitudes(amplitudes) = &cli.state.process_list.processes[0].collection
    else {
        panic!("expected the dotted-sunset amplitude")
    };
    let [sunset] = amplitudes["default"].graphs.as_slice() else {
        panic!("expected a single genuine two-loop graph")
    };
    let basis = &sunset.graph.loop_momentum_basis;
    assert_eq!(
        basis
            .loop_edges
            .iter()
            .map(|edge| edge.0)
            .collect::<Vec<_>>(),
        [2, 3],
        "the input basis must use the undotted fermion and the scalar chain"
    );
    for (edge, signature) in &basis.edge_signatures {
        if edge.0 < 2 {
            assert!(
                signature
                    .internal
                    .iter()
                    .all(|sign| *sign != SignOrZero::Zero),
                "the raised fermion route must involve both input loop momenta"
            );
        }
    }
    // Independent finite-density cut integrals, differentiated with respect
    // to the common fermion mass squared. The derivation and quadrature
    // convergence are recorded in fermi_surface_2l_integration.toml. Parameters are
    // m_b=3, mu_b=5, M_H=2; all UV subgraphs converge. Two loops and five
    // propagator factors give Wick phase i²(-1)^5=+1, so this target is real.
    let expected = -4.880_698_568_451e-6_f64;
    for (profile, sigma) in [("exponential", 1.0), ("poly_exponential", 1.3)] {
        cli.run_command(&format!(
            "set process kv h_function.function=\"{profile}\" h_function.sigma={sigma}"
        ))?;
        // As in the one-loop test, the 20k-sample estimate is archived by
        // the 100k-sample run.
        let (output, partial, elapsed) = integrate_with_20k_sample_archive(&mut cli, name)?;
        let mut errors = Vec::new();
        for (samples, result) in [(20_000, &partial), (100_000, &output.result)] {
            let estimate = single_slot_integral(result);
            let (value, error) = (estimate.result.re.0, estimate.error.re.0);
            info!(
                profile,
                sigma,
                samples = estimate.neval,
                value,
                error,
                expected,
                ?elapsed,
                "Non-aligned two-loop Fermi-surface integration convergence"
            );
            assert_eq!(estimate.neval, samples, "the sample budget must be fixed");
            assert!(
                [value, error, estimate.result.im.0, estimate.error.im.0]
                    .into_iter()
                    .all(f64::is_finite)
                    && value < 0.0
                    && error > 0.0
                    && estimate.error.im.0 >= 0.0,
                "{profile}, samples={samples}: {estimate}"
            );
            assert!(
                estimate.result.im.0.abs() < 1.0e-12 * expected.abs()
                    && estimate.error.im.0 < 1.0e-12 * expected.abs(),
                "the imaginary part must vanish: {estimate}"
            );
            assert!(
                (value - expected).abs() <= 5.0 * error,
                "{profile}, samples={samples}: {value:e} +/- {error:e}, target={expected:e}"
            );
            errors.push(error);
        }
        assert!(
            errors[1] < 0.01 * expected.abs() && errors[1] < 0.7 * errors[0],
            "{profile}: increased statistics must reduce uncertainty below 1%: {errors:?}"
        );
    }
    clean_test(&cli.cli_settings.state.folder);
    Ok(())
}

#[test]
fn v_diag() -> Result<()> {
    let cli = get_test_cli(
        Some("v_diag.toml".into()),
        get_tests_workspace_path().join("v_diag"),
        Some("v_diag".to_string()),
        true,
    )?;
    Ok(())
}

#[test]
fn scalar_bubble() -> Result<()> {
    let mut cli = get_test_cli(
        Some("scalar_bubble.toml".into()),
        get_tests_workspace_path().join("scalar_bubble"),
        Some("scalar_bubble".to_string()),
        false,
    )?;

    cli.run_command("run import_graph")?;

    cli.run_command("run no_integrated")?;

    cli.run_command("generate")?;

    let integrate_command = Integrate {
        ..default_integrate_for("scalar_bubble")
    };

    let profile_cmd = Profile::UltraViolet(UltraVioletProfile {
        ..Default::default()
    });

    // Kaapo reference magnitudes: m=1 muv=5 2.03838e-02 m=2 muv=5 	1.16050e-02	 m=3 muv=5 6.46968e-03

    cli.run_command("set model mass_scalar_1=1.0")?;
    let res = profile_cmd
        .run(&mut cli.state, &cli.cli_settings)?
        .unwrap_uv();
    assert_eq!(res.pass_fail(-0.9).failed, 0);
    let integral_no_cache = integrate_command.run(&mut cli.state, &cli.cli_settings)?;
    assert!(integral_no_cache.is_compatible_with_target(Complex::new(F(0.0), F(-2.03838e-02)), 2));

    cli.run_command("set model mass_scalar_1=2.0")?;
    let integral_no_cache = integrate_command.run(&mut cli.state, &cli.cli_settings)?;
    assert!(integral_no_cache.is_compatible_with_target(Complex::new(F(0.0), F(-1.16050e-02)), 2));

    cli.run_command("set model mass_scalar_1=3.0")?;
    let integral_no_cache = integrate_command.run(&mut cli.state, &cli.cli_settings)?;
    assert!(integral_no_cache.is_compatible_with_target(Complex::new(F(0.0), F(-6.46968e-03)), 2));
    let renorm_command = Renormalize::default();

    let res = renorm_command.run(&mut cli.state, &cli.cli_settings)?;

    println!("{}", res[0]);

    clean_test(&cli.cli_settings.state.folder);

    Ok(())
}

/// Breit-Wheeler pair production a a > t t~ at sqrt(s) = 1 TeV with mt = 173 GeV
/// and alpha = 1/132.507, generated with the extra `generate xs` options
/// `generation_options`. Stability checks stay in double precision, so a basis
/// that covers the cut phase space badly fails within seconds instead of
/// escalating its unstable points to arbitrary precision for many minutes.
fn aa_ttx_born_cli(
    test_name: &str,
    generation_options: &str,
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
            "import model sm-default",
            "set global kv global.generation.evaluator.compile=false global.generation.evaluator.summed=false global.generation.evaluator.iterative_orientation_optimization=false global.generation.evaluator.store_atom=false global.generation.uv.subtract_uv=false global.generation.uv.generate_integrated=false global.generation.threshold_subtraction.enable_thresholds=false global.generation.tropical_subgraph_table.disable_tropical_generation=true",
            "set model MT=173.0",
            "set model WT=0.0",
            "set model ymt=173.0",
            "set model aEWM1=132.507",
            r#"set default-runtime string '
[general]
evaluator_method = "SingleParametric"
enable_cache = false
generate_events = false
disable_flux_factor = true

[subtraction]
disable_threshold_subtraction = true

[[stability.levels]]
precision = "Double"
required_precision_for_re = 1e-5
required_precision_for_im = 1e-5
escalate_for_large_weight_threshold = 0.9

[sampling]
graphs = "summed"
orientations = "summed"
lmb_multichanneling = false
coordinate_system = "spherical"
mapping = "linear"
b = 1.0

[kinematics]
e_cm = 1000.0

[kinematics.externals]
type = "constant"

[kinematics.externals.data]
momenta = [[500.0, 0.0, 0.0, 500.0], [500.0, 0.0, 0.0, -500.0]]
helicities = ["summed_averaged", "summed_averaged"]

[integrator]
integrated_phase = "real"
n_start = 50000
n_increase = 50000
n_max = 50000
min_samples_for_update = 1000
seed = 7331
'"#,
            &format!(
                "generate xs a a > t t~ | a t t~ g ghG ghG~ QED^2==4 [{{{{1}}}} QCD=0] --numerator-grouping no_grouping --symmetrize-initial-states false --symmetrize-final-states false {generation_options} -p aa_ttx -i LO --filter-selfenergies true --filter-snails true --filter-tadpoles true --n-cut-blobs 1 1 --n-cut-spectators 0 0"
            ),
        ],
    )?;
    Ok(cli)
}

/// The Breit-Wheeler cross section is sigma = 0.459677 pb. Without the flux
/// factor the integrand returns sigma 2s in natural units.
#[test]
fn aa_ttx_born_cross_section_integration_matches_breit_wheeler() -> Result<()> {
    const GEV_SQUARED_TO_PICOBARN: f64 = 0.389379304e9;
    const BREIT_WHEELER_PICOBARN: f64 = 0.459677;
    const E_CM: f64 = 1000.0;
    let test_name = "aa_ttx_born_breit_wheeler";
    let mut cli = aa_ttx_born_cli(test_name, "")?;

    let integral = Integrate {
        integrand_name: vec!["LO".to_string()],
        ..default_integrate_for(test_name)
    }
    .run(&mut cli.state, &cli.cli_settings)?;
    let integral = single_slot_integral(&integral);
    let target = BREIT_WHEELER_PICOBARN / GEV_SQUARED_TO_PICOBARN * 2.0 * E_CM * E_CM;
    // A badly covered cut phase space inflates the variance first.
    assert!(
        integral.error.re.0 < 0.03 * target,
        "{integral:#} has a large error for target {target:e}"
    );
    assert!(
        integral.is_compatible_with_target(Complex::new(F(target), F(0.0)), 3),
        "{integral:#} is not compatible with target {target:e}"
    );

    clean_test(&cli.cli_settings.state.folder);

    Ok(())
}

/// The loop momenta are rescaled about the origin of the generated basis. With
/// the cut edges at K + P1 and K - P0 the cut energies there add up to
/// 2 sqrt(500^2 + 173^2) GeV, 58.2 GeV above sqrt(s), so the origin lies outside
/// every cut surface. Integrating with that basis used to give 6.0(3.6)e-3
/// instead of 2.36e-3, with 17.5% unstable points; warm-up now rejects it.
#[test]
fn aa_ttx_born_rejects_a_basis_whose_origin_lies_outside_a_cut_surface() -> Result<()> {
    let test_name = "aa_ttx_born_origin_outside_cut";
    let mut cli = aa_ttx_born_cli(test_name, "--loop-momentum-bases GL0=2 GL1=3 GL2=2 GL3=3")?;

    let Err(error) = (Integrate {
        integrand_name: vec!["LO".to_string()],
        ..default_integrate_for(test_name)
    })
    .run(&mut cli.state, &cli.cli_settings) else {
        panic!("integrating with the rescaling centre outside the cut surfaces succeeded");
    };
    assert!(
        error.chain().any(|cause| {
            let cause = cause.to_string();
            cause.contains("lies on or outside the energy surface of cut group")
                && cause.contains("--loop-momentum-bases GL")
        }),
        "unexpected error: {error:?}"
    );

    clean_test(&cli.cli_settings.state.folder);

    Ok(())
}

/// The rejected group remains in the generated integrand, but is never warmed or hashed
/// when the saved runtime selection excludes it.
#[test]
fn selected_group_integrates_and_resumes_with_an_excluded_invalid_basis() -> Result<()> {
    use gammalooprs::integrands::process::ProcessIntegrand;
    use gammalooprs::integrate::{IntegrationState, workspace_state_path};

    let name = "selected_group_excludes_invalid_basis";
    let bad_name = "selected_group_invalid_basis_fixture";
    let mut cli = aa_ttx_born_cli(name, "")?;
    let mut bad_cli = aa_ttx_born_cli(bad_name, "--loop-momentum-bases GL0=2 GL1=3 GL2=2 GL3=3")?;
    for state in [&mut cli, &mut bad_cli] {
        state.run_command(
            "set process -p aa_ttx -i LO string '[sampling]\ngraphs = \"monte_carlo\"\ngraph_names = []\n'",
        )?;
    }
    let (process_id, integrand_name) =
        cli.state.find_integrand_ref(None, Some(&"LO".to_owned()))?;
    let (bad_id, bad_integrand) = bad_cli
        .state
        .find_integrand_ref(None, Some(&"LO".to_owned()))?;
    let model = cli
        .state
        .resolve_model_for_integrand(process_id, &integrand_name)?;
    let bad_source = bad_cli
        .state
        .process_list
        .get_integrand(bad_id, &bad_integrand)?
        .require_generated()?
        .clone();
    let mut source = cli
        .state
        .process_list
        .get_integrand(process_id, &integrand_name)?
        .require_generated()?
        .clone();
    assert_eq!(
        source.graph_group_count(),
        source.graph_count(),
        "fixture explicitly disables grouping"
    );

    // Classify by the actual basis-origin check, then resolve names from those groups.
    let mut rejected_name = None;
    for master in bad_source.graph_group_master_names() {
        let mut group =
            bad_source.clone_with_selected_graph_groups(&model, &[master.to_owned()])?;
        if let Err(error) = group.warm_up(&model) {
            assert!(
                format!("{error:?}").contains("lies on or outside the energy surface of cut group")
            );
            rejected_name = Some(master.to_owned());
            break;
        }
    }
    let rejected_name = rejected_name.expect("fixture must contain an invalid basis");
    let (ProcessIntegrand::CrossSection(bad), ProcessIntegrand::CrossSection(good)) =
        (&bad_source, &mut source)
    else {
        unreachable!()
    };
    let replacement = bad
        .data
        .graph_terms
        .iter()
        .find(|term| term.graph.name == rejected_name)
        .unwrap();
    let destination = good
        .data
        .graph_terms
        .iter_mut()
        .find(|term| term.graph.name == rejected_name)
        .unwrap();
    let group_id = destination.graph.group_id;
    *destination = replacement.clone();
    destination.graph.group_id = group_id;
    destination.multi_channeling_setup.graph.group_id = group_id;
    let mut accepted_name = None;
    for master in source.graph_group_master_names() {
        if master == rejected_name {
            continue;
        }
        let mut group = source.clone_with_selected_graph_groups(&model, &[master.to_owned()])?;
        group.warm_up(&model)?;
        accepted_name = Some(master.to_owned());
        break;
    }
    let accepted_name = accepted_name.expect("fixture must retain a valid group");
    *cli.state
        .process_list
        .get_integrand_mut(process_id, &integrand_name)? = source;
    cli.run_command(&format!("set process -p aa_ttx -i LO string '[sampling]\ngraphs = \"monte_carlo\"\ngraph_names = [\"{accepted_name}\"]\n[integrator]\nn_bins = 4\nmin_samples_for_update = 1\nn_start = 8\nn_increase = 8\nn_max = 8\n'"))?;
    let before = cli
        .state
        .process_list
        .get_integrand(process_id, &integrand_name)?
        .require_generated()?
        .resume_fingerprint()?;
    let mut integrate = Integrate {
        integrand_name: vec![integrand_name.clone()],
        ..default_integrate_for(name)
    };
    let first = integrate.run(&mut cli.state, &cli.cli_settings)?;
    let state_path = workspace_state_path(&first.workspace_path);
    let (saved, _): (IntegrationState, _) =
        bincode::decode_from_slice(&std::fs::read(&state_path)?, bincode::config::standard())?;
    assert!(saved.iter > 0 && saved.num_points > 0);
    assert_eq!(
        cli.state
            .process_list
            .get_integrand(process_id, &integrand_name)?
            .require_generated()?
            .resume_fingerprint()?,
        before
    );
    // Current empty selection is ignored on resume: the saved nonempty selection is effective.
    cli.run_command(
        "set process -p aa_ttx -i LO string '[sampling]\ngraphs = \"monte_carlo\"\ngraph_names = []\n'",
    )?;
    integrate.restart = false;
    integrate.run(&mut cli.state, &cli.cli_settings)?;
    let (resumed, _): (IntegrationState, _) =
        bincode::decode_from_slice(&std::fs::read(&state_path)?, bincode::config::standard())?;
    assert_eq!(
        (resumed.iter, resumed.num_points),
        (saved.iter, saved.num_points)
    );
    assert_eq!(
        cli.state
            .process_list
            .get_integrand(process_id, &integrand_name)?
            .require_generated()?
            .resume_fingerprint()?,
        before
    );

    // A restart chooses current settings; selecting the rejected group or every group
    // retains the existing physical error and does not destroy the accepted workspace.
    let protected = std::fs::read(&state_path)?;
    for selected in [format!("\"{rejected_name}\""), String::new()] {
        cli.run_command(&format!("set process -p aa_ttx -i LO string '[sampling]\ngraphs = \"monte_carlo\"\ngraph_names = [{selected}]\n'"))?;
        integrate.restart = true;
        let error = integrate
            .run(&mut cli.state, &cli.cli_settings)
            .expect_err("invalid selected basis");
        assert!(
            format!("{error:?}").contains("lies on or outside the energy surface of cut group")
        );
        assert_eq!(std::fs::read(&state_path)?, protected);
    }
    clean_test(&cli.cli_settings.state.folder);
    clean_test(&bad_cli.cli_settings.state.folder);
    Ok(())
}

#[test]
fn vakint_compare_scalar_bubble() -> Result<()> {
    let name = "vakint_compare_scalar_bubble";
    let mut cli = get_test_cli(
        Some("scalar_bubble.toml".into()),
        get_tests_workspace_path().join(name),
        Some(name.to_string()),
        false,
    )?;

    cli.run_command("run import_graph")?;

    cli.run_command("run integrated")?;
    // At p²=0 the scalar bubble is -i[1/epsilon + log(mu_r²/m²)]/(16*pi²).
    // Local subtraction removes the same expression with m=m_uv; its finite
    // integrated addback leaves the MS-bar value -i*log(mu_r²/m²)/(16*pi²).
    // The run card fixes mu_r=1000 independently of its local UV mass m_uv=5.

    cli.run_command("generate")?;

    let integrate_command = Integrate {
        ..default_integrate_for(name)
    };

    let profile_cmd = Profile::UltraViolet(UltraVioletProfile {
        ..Default::default()
    });

    cli.run_command("set model mass_scalar_1=1.0")?;
    let res = profile_cmd
        .run(&mut cli.state, &cli.cli_settings)?
        .unwrap_uv();
    assert_eq!(res.pass_fail(-0.9).failed, 0);
    let integral_no_cache = integrate_command.run(&mut cli.state, &cli.cli_settings)?;
    assert!(
        integral_no_cache
            .is_compatible_with_target(Complex::new(F(0.0), F(-0.08748774264725967)), 2)
    );

    cli.run_command("set model mass_scalar_1=2.0")?;
    let integral_no_cache = integrate_command.run(&mut cli.state, &cli.cli_settings)?;
    assert!(
        integral_no_cache
            .is_compatible_with_target(Complex::new(F(0.0), F(-0.07870893105067431)), 2)
    );

    cli.run_command("set model mass_scalar_1=3.0")?;
    let integral_no_cache = integrate_command.run(&mut cli.state, &cli.cli_settings)?;
    assert!(
        integral_no_cache
            .is_compatible_with_target(Complex::new(F(0.0), F(-0.07357365546577585)), 2)
    );
    let renorm_command = Renormalize::default();

    let res = renorm_command.run(&mut cli.state, &cli.cli_settings)?;

    println!("{}", res[0]);

    clean_test(&cli.cli_settings.state.folder);

    Ok(())
}

#[test]
fn scalar_mercedes() -> Result<()> {
    let mut cli = get_test_cli(
        Some("scalar_mercedes.toml".into()),
        get_tests_workspace_path().join("scalar_mercedes"),
        Some("scalar_mercedes".to_string()),
        false,
    )?;

    let integrate_command = Integrate {
        ..default_integrate_for("scalar_mercedes")
    };

    let profile_cmd = Profile::UltraViolet(UltraVioletProfile {
        max_scale_exponent: 4.0,
        use_f128: false,
        ..Default::default()
    });

    //5.89551e-06	3.35645e-06	1.87120e-06

    // Kaapo magnitudes: m=1 muv=5 5.89551e-06 m=2 muv=5 	3.35645e-06	 m=3 muv=5 1.87120e-06
    // This primitive K4 vacuum graph has period 6*zeta(3), so local subtraction
    // gives i*6*zeta(3)*log(25/m^2)/(16*pi^2)^3. Vertex numerators are one;
    // six propagator i factors and three Wick rotations give +i.
    // Independent period: https://arxiv.org/abs/hep-ph/9504352.
    cli.run_command("set model mass_scalar_1=1.0")?;
    let res = profile_cmd
        .run(&mut cli.state, &cli.cli_settings)?
        .unwrap_uv();
    assert_eq!(res.pass_fail(-0.9).failed, 0);
    let integral_no_cache = integrate_command.run(&mut cli.state, &cli.cli_settings)?;
    assert!(
        integral_no_cache
            .is_compatible_with_target(Complex::new(F(0.0), F(5.895509129899469e-06)), 1),
        "Not compatible: {integral_no_cache}",
    );

    cli.run_command("set model mass_scalar_1=2.0")?;
    let integral_no_cache = integrate_command.run(&mut cli.state, &cli.cli_settings)?;
    assert!(
        integral_no_cache
            .is_compatible_with_target(Complex::new(F(0.0), F(3.3564515497441014e-06)), 3),
        "Not compatible: {integral_no_cache}",
    );

    cli.run_command("set model mass_scalar_1=3.0")?;
    let integral_no_cache = integrate_command.run(&mut cli.state, &cli.cli_settings)?;
    assert!(
        integral_no_cache
            .is_compatible_with_target(Complex::new(F(0.0), F(1.8711980781814097e-06)), 3),
        "Not compatible: {integral_no_cache}",
    );

    clean_test(&cli.cli_settings.state.folder);

    Ok(())
}

#[test]
fn scalar_sunrise() -> Result<()> {
    let mut cli = get_test_cli(
        Some("scalar_sunrise.toml".into()),
        get_tests_workspace_path().join("scalar_sunrise"),
        Some("scalar_sunrise".to_string()),
        false,
    )?;

    let integrate_command = Integrate {
        ..default_integrate_for("scalar_sunrise")
    };
    let profile_cmd = Profile::UltraViolet(UltraVioletProfile {
        use_f128: false,
        max_scale_exponent: 6.0,
        min_scale_exponent: 1.0,
        ..Default::default()
    });
    // Kaapo reference magnitudes: m=1 muv=5 4.37688e-03 m=2 muv=5 2.48100e-03 m=3 muv=5 1.07231e-03
    // At p²=0 the complete local forest is (1-T₂)[S(m²)-3*B(m_uv²)*T(m²)].
    // Writing x=m², X=m_uv² and L=log(x/X), its finite integral is
    // [6*(X-x+x*L)-3*x*L²/2]/(16*pi²)². The sunset and tadpole poles fix
    // this expression; their unknown finite master parts cancel.
    // The literal numerator is one: two loops and three propagators give
    // a positive Wick phase. The original budget and one-sigma gate remain.
    cli.run_command("set model mass_scalar_1=1.0")?;
    let res = profile_cmd
        .run(&mut cli.state, &cli.cli_settings)?
        .unwrap_uv();
    assert_eq!(res.pass_fail(-0.9).failed, 0);
    let integral_no_cache = integrate_command.run(&mut cli.state, &cli.cli_settings)?;
    assert!(
        integral_no_cache.is_compatible_with_target(Complex::new_re(F(4.37688e-03)), 1) // .result
                                                                                        // .approx_eq(&Complex::new_re(F(4.37688e-03)), &F(0.01))
    );

    // assert_snapshot!(format!("{integral_no_cache:.3}"),@"4.359e-3");

    cli.run_command("set model mass_scalar_1=2.0")?;
    let res = profile_cmd
        .run(&mut cli.state, &cli.cli_settings)?
        .unwrap_uv();
    assert_eq!(res.pass_fail(-0.9).failed, 0);
    let integral_no_cache = integrate_command.run(&mut cli.state, &cli.cli_settings)?;
    // assert_snapshot!(format!("{integral_no_cache:.3}"),@"2.474e-3");
    assert!(integral_no_cache.is_compatible_with_target(Complex::new_re(F(2.48100e-03)), 1));

    // cli.run_command("set model mass_scalar_1=3.0")?;
    // let integral_no_cache = integrate_command.run(&mut cli.state, &cli.cli_settings)?;
    // assert_snapshot!(format!("{:.8e}",integral_no_cache.result),@"(-4.5184321377520566e-4+0e0i)");

    let renorm_command = Renormalize::default();

    let res = renorm_command.run(&mut cli.state, &cli.cli_settings)?;

    println!("{}", res[0]);
    clean_test(&cli.cli_settings.state.folder);

    Ok(())
}

#[test]
fn scalar_box() -> Result<()> {
    let mut cli = get_test_cli(
        Some("scalar_box.toml".into()),
        get_tests_workspace_path().join("scalar_box"),
        Some("scalar_box".to_string()),
        false,
    )?;

    let (_, a) = Inspect {
        process: None,
        integrand_name: Some("default".to_string()),
        point: vec![0.1, 0.2, 0.3],
        momentum_space: true,
        ..Default::default()
    }
    .run(&mut cli)?;

    // In the default chart K=q5, P0=q3, P1=q0, P2=q1, the four lower
    // energy-pole residues give +i*0.00032031590000891889816516.
    // This independently fixes the point before comparing cached integrals.
    assert_snapshot!(format!("{a:.8e}"),@"(0.00000000e0+3.20315900e-4i)");

    cli.run_command("set process -p 0 -i default kv general.enable_cache=false")?;
    let integral_no_cache = Integrate {
        ..default_integrate_for("scalar_box")
    }
    .run(&mut cli.state, &cli.cli_settings)?;

    cli.run_command("set process -p 0 -i default kv general.enable_cache=true")?;
    let integral_with_cache = Integrate {
        ..default_integrate_for("scalar_box")
    }
    .run(&mut cli.state, &cli.cli_settings)?;

    info!("Integral result without caching: {:#?}", integral_no_cache);
    info!(
        "Integral result with caching.  : {:#?}",
        integral_with_cache
    );

    let integral_no_cache = single_slot_integral(&integral_no_cache);
    let integral_with_cache = single_slot_integral(&integral_with_cache);

    assert_approx_eq(
        &integral_no_cache.result.re,
        &integral_with_cache.result.re,
        &F(1e-4),
    );
    assert_approx_eq(
        &integral_no_cache.result.im,
        &integral_with_cache.result.im,
        &F(1e-4),
    );

    clean_test(&cli.cli_settings.state.folder);

    Ok(())
}

mod failing {
    use super::*;

    #[test]
    fn test_grouped_subtraction() -> Result<()> {
        let mut cli = get_test_cli(
            Some("test_grouped_subtraction.toml".into()),
            get_tests_workspace_path().join("test_grouped_subtraction"),
            None,
            false,
        )?;

        let int1 = Integrate {
            process: vec![ProcessRef::Id(0)],
            integrand_name: vec!["default".to_string()],
            workspace_path: Some(
                get_tests_workspace_path()
                    .join("test_grouped_subtraction/integration_workspace_no_group/"),
            ),
            target: vec![],
            n_cores: Some(1),
            restart: true,
            ..Default::default()
        };
        let int2 = Integrate {
            process: vec![ProcessRef::Id(1)],
            integrand_name: vec!["default".to_string()],
            workspace_path: Some(
                get_tests_workspace_path()
                    .join("test_grouped_subtraction/integration_workspace_group/"),
            ),
            target: vec![],
            n_cores: Some(1),
            restart: true,
            ..Default::default()
        };

        let integration_results_no_group = int1.run(&mut cli.state, &cli.cli_settings)?;
        let integration_results_group = int2.run(&mut cli.state, &cli.cli_settings)?;

        info!("No group result: {:#?}", integration_results_no_group);
        info!("Group result: {:#?}", integration_results_group);

        let integration_results_no_group = single_slot_integral(&integration_results_no_group);
        let integration_results_group = single_slot_integral(&integration_results_group);

        assert_approx_eq(
            &integration_results_group.result.re,
            &integration_results_no_group.result.re,
            &F(1e-1),
        );
        assert_approx_eq(
            &integration_results_group.result.im,
            &integration_results_no_group.result.im,
            &F(1e-1),
        );

        clean_test(&cli.cli_settings.state.folder);

        Ok(())
    }

    #[test]
    fn scalar_basketball() -> Result<()> {
        let mut cli = get_test_cli(
            Some("scalar_basketball.toml".into()),
            get_tests_workspace_path().join("scalar_basketball"),
            Some("scalar_basketball".to_string()),
            false,
        )?;

        let integrate_command = Integrate {
            ..default_integrate_for("scalar_basketball")
        };

        // The seed-42 rays are preasymptotic over 10²..10⁴; their slopes approach
        // -1 over higher ranges, with no nondecaying term on exact representative
        // rays. Retain 20 probes, the -0.9 gate and the million-sample MC budget.
        let profile_cmd = Profile::UltraViolet(UltraVioletProfile {
            use_f128: false,
            max_scale_exponent: 5.0,
            min_scale_exponent: 3.0,
            n_points: 20,
            analyse_analytically: false,
            ..Default::default()
        });

        // Kaapo's mUV=5 vacuum (p^2=0) values were 1.47240e-03, 7.15184e-04,
        // and 2.27485e-04 for m=1,2,3. At the card's p^2=1, the complete
        // 46-term local forest adds the subtracted momentum Taylor terms and
        // finite Schwinger remainder, independently evaluated from vacuum masters.
        // Four propagator i factors and three loops convert the historical
        // Euclidean normalization to the physical -i targets below.
        cli.run_command("set model mass_scalar_1=1.0")?;
        let res = profile_cmd
            .run(&mut cli.state, &cli.cli_settings)?
            .unwrap_uv();
        assert_eq!(res.pass_fail(-0.9).failed, 0);
        let integral_no_cache = integrate_command.run(&mut cli.state, &cli.cli_settings)?;
        assert!(
            integral_no_cache.is_compatible_with_target(Complex::new(F(0.0), F(-1.49689e-03)), 1),
            "Not compatible: {integral_no_cache}",
        );

        // cli.run_command("set model mass_scalar_1=2.0")?;
        // let res = profile_cmd.run(&mut cli.state, &cli.cli_settings)?;
        // assert_eq!(res.pass_fail(-0.9).failed, 0);
        // let integral_no_cache = integrate_command.run(&mut cli.state, &cli.cli_settings)?;
        // assert!(
        //     integral_no_cache.is_compatible_with_target(Complex::new(F(0.0), F(-7.29603e-04)), 3),
        //     "Not compatible: {integral_no_cache}",
        // );

        // cli.run_command("set model mass_scalar_1=3.0")?;
        // let res = profile_cmd.run(&mut cli.state, &cli.cli_settings)?;
        // assert_eq!(res.pass_fail(-0.9).failed, 0);
        // let integral_no_cache = integrate_command.run(&mut cli.state, &cli.cli_settings)?;
        // assert!(
        //     integral_no_cache.is_compatible_with_target(Complex::new(F(0.0), F(-2.33978e-04)), 3),
        //     "Not compatible: {integral_no_cache}",
        // );

        // clean_test(&cli.cli_settings.state.folder);

        Ok(())
    }

    #[test]
    fn test_integrate_bubble_dot_1t1b() -> Result<()> {
        let target = Complex::new(F(-0.00041875868403258484), F(0.0));

        let mut cli = get_test_cli(
            Some("bubble_dot_1t1b_generate.toml".into()),
            get_tests_workspace_path().join("bubble_dot_1t1b"),
            None,
            true,
        )?;

        let integrate_command = Integrate {
            workspace_path: Some(
                get_tests_workspace_path().join("bubble_dot_1t1b/integration_workspace"),
            ),
            n_cores: Some(1),
            ..Default::default()
        };

        let integration_result = integrate_command.run(&mut cli.state, &cli.cli_settings)?;
        assert!(
            integration_result.is_compatible_with_target(target, 2),
            "Integration result {integration_result} not compatible with target {target}"
        );

        Ok(())
    }

    #[test]
    fn test_integrate_bubble_dot_2t0b() -> Result<()> {
        let target = Complex::new(F(-2.9911664985748502e-05), F(0.0));

        let mut cli = get_test_cli(
            Some("bubble_dot_2t0b_generate.toml".into()),
            get_tests_workspace_path().join("bubble_dot_2t0b"),
            None,
            true,
        )?;

        let integrate_command = Integrate {
            workspace_path: Some(
                get_tests_workspace_path().join("bubble_dot_2t0b/integration_workspace"),
            ),
            n_cores: Some(1),

            restart: true,
            ..Default::default()
        };

        let integration_result = integrate_command.run(&mut cli.state, &cli.cli_settings)?;
        assert!(
            integration_result.is_compatible_with_target(target, 2),
            "Integration result {integration_result} not compatible with target {target}"
        );
        Ok(())
    }
}

mod slow {
    use super::*;

    #[test]
    fn scalar_mercedes_with_extra_loop() -> Result<()> {
        let mut cli = get_test_cli(
            Some("scalar_mercedes_with_extra_loop.toml".into()),
            get_tests_workspace_path().join("scalar_mercedes_with_extra_loop"),
            Some("scalar_mercedes_with_extra_loop".to_string()),
            false,
        )?;

        let integrate_command = Integrate {
            ..default_integrate_for("scalar_mercedes_with_extra_loop")
        };

        let profile_cmd = Profile::UltraViolet(UltraVioletProfile {
            max_scale_exponent: 7.,
            min_scale_exponent: 4.0,
            n_points: 15,
            use_f128: false,
            output_file: Some("uv_profile_extra_loop".into()),
            ..Default::default()
        });
        // At m=m_uv=1 and Q=-p^2=12, every proper logarithmic subgraph
        // contracts both external vertices. After regulated integration its
        // vacuum term cancels against the parent, leaving S(Q)-S(0)-Q*S'(0).
        // Independent Schwinger parameters give +i/(16*pi^2)^4 times the
        // simplex integral of U^-2*((1+Q*V/U)*ln(1+Q*V/U)-Q*V/U).
        // For m=2,3 add the vacuum forest (Laporta V5 and lower-loop masters)
        // and Q*ln(m^2/m_uv^2)*P, where P=(16*pi^2)^-4 integral(V/U^3).
        // This independent forest also reproduces the old m_uv=5,p=0 targets.
        cli.run_command("set model mass_scalar_1=1.0")?;
        let res = profile_cmd
            .run(&mut cli.state, &cli.cli_settings)?
            .unwrap_uv();
        assert_eq!(res.pass_fail(-0.9).failed, 0);

        let integral_no_cache = integrate_command.run(&mut cli.state, &cli.cli_settings)?;
        assert!(
            integral_no_cache.is_compatible_with_target(Complex::new(F(0.0), F(6.88106e-09)), 1),
            "Not compatible: {integral_no_cache}",
        );

        cli.run_command("set model mass_scalar_1=2.0")?;
        let integral_no_cache = integrate_command.run(&mut cli.state, &cli.cli_settings)?;
        assert!(
            integral_no_cache.is_compatible_with_target(Complex::new(F(0.0), F(2.38360e-07)), 3),
            "Not compatible: {integral_no_cache}",
        );

        cli.run_command("set model mass_scalar_1=3.0")?;
        let integral_no_cache = integrate_command.run(&mut cli.state, &cli.cli_settings)?;
        assert!(
            integral_no_cache.is_compatible_with_target(Complex::new(F(0.0), F(6.46503e-07)), 3),
            "Not compatible: {integral_no_cache}",
        );

        clean_test(&cli.cli_settings.state.folder);

        Ok(())
    }

    #[test]
    fn test_integrate_dotted_bubble() -> Result<()> {
        // The two-point insertion differentiates one mass in the physical
        // two-body phase space: dPhi_2/dm_1^2 = -1/(8*pi*s*beta).
        // Here s=16 and both masses are one; the graph's complete numerator
        // is one, as are its three (-i) vertices times three i propagators.
        let target = Complex::new(F(-0.002871504663657095), F(0.0));
        let mut cli = get_test_cli(
            Some("dotted_bubble_generate.toml".into()),
            get_tests_workspace_path().join("dotted_bubble"),
            None,
            true,
        )?;

        let integrate_command = Integrate {
            workspace_path: Some(
                get_tests_workspace_path().join("dotted_bubble/integration_workspace"),
            ),
            n_cores: Some(1),

            restart: true,
            ..Default::default()
        };

        let integration_result = integrate_command.run(&mut cli.state, &cli.cli_settings)?;
        assert!(
            integration_result.is_compatible_with_target(target, 2),
            "Integration result {integration_result} not compatible with target {target}"
        );

        Ok(())
    }
}
