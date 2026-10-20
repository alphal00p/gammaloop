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
    use gammalooprs::{momentum::SignOrZero, processes::ProcessCollection};

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
fn box_4e_multichannel_subtraction_reproduces_published_integral() -> Result<()> {
    let _serial_guard = gammaloop_integration_tests::SERIAL_TEST_LOCK.lock();
    use gammalooprs::{
        integrands::process::ProcessIntegrand, settings::global::RepresentationMode,
    };

    let test_name = "box_4e_multichannel_subtraction";
    let root = get_tests_workspace_path().join(test_name);
    let mut cli = get_test_cli(None, &root, Some(test_name.to_owned()), true)?;
    cli.cli_settings
        .global
        .generation
        .three_dimensional_representations = vec![RepresentationMode::Cff, RepresentationMode::Ltd];
    let graph_path = root.join("box4e.dot");
    // Eq. (3.10), Section 3.1 of https://arxiv.org/pdf/1912.09291.
    // Fix every numerator to one: the paper's Eq. (2.1) integrates 1/prod(D)
    // with d^4k/(2pi)^4, without vertex/propagator i factors or couplings.
    // GammaLoop's signed dq0/(2pi*i) contour and i/(2pi)^3 normalization
    // already give that convention; no additional phase or normalization applies.
    std::fs::write(
        &graph_path,
        r#"digraph box4e {
            num=1; overall_factor=1;
            node [num=1]; edge [num=1];
            p1 [style=invis]; p2 [style=invis];
            p3 [style=invis]; p4 [style=invis];
            p1 -> A:0 [id=0, particle=scalar_0];
            p2 -> B:1 [id=1, particle=scalar_0];
            p3 -> C:2 [id=2, particle=scalar_0];
            D:3 -> p4 [id=3, particle=scalar_0];
            A -> B [id=4, particle=scalar_0];
            B -> C [id=5, particle=scalar_0];
            C -> D [id=6, particle=scalar_0];
            D -> A [id=7, particle=scalar_0, lmb_id=0];
        }"#,
    )?;
    run_commands(
        &mut cli,
        &[
            "import model scalars-default.json",
            "set global kv global.generation.evaluator.compile=false global.generation.evaluator.iterative_orientation_optimization=false global.generation.evaluator.store_atom=false global.generation.evaluator.summed=false global.generation.evaluator.summed_function_map=true global.generation.uv.local_uv_cts_from_expanded_4d_integrands=true global.generation.uv.subtract_uv=false global.generation.uv.generate_integrated=false global.generation.threshold_subtraction.enable_thresholds=true global.generation.threshold_subtraction.disable_integrated_ct=false global.generation.threshold_subtraction.check_esurface_at_generation=true global.generation.tropical_subgraph_table.disable_tropical_generation=true",
            r#"set default-runtime string '
[general]
evaluator_method = "SummedFunctionMap"
integral_unit = "none"
disable_flux_factor = true
[kinematics]
e_cm = 100.0
[kinematics.externals]
type = "constant"
[kinematics.externals.data]
momenta = [[14.0,-6.6,-40.0,0.0],[-43.0,15.2,33.0,0.0],[-17.9,-50.0,11.8,0.0],"dependent"]
helicities = [0,0,0,0]
[subtraction]
disable_threshold_subtraction = false
[stability]
rotation_axis = [{ type = "euler_angles", alpha = 0.1, beta = 0.2, gamma = 0.3 }]
levels = [
 { precision = "Double", required_precision_for_re = 1e-7, required_precision_for_im = 1e-7, escalate_for_large_weight_threshold = -1.0 },
 { precision = "Quad", required_precision_for_re = 1e-7, required_precision_for_im = 1e-7, escalate_for_large_weight_threshold = -1.0 },
]
[sampling]
graphs = "summed"
orientations = "summed"
lmb_multichanneling = false
lmb_channels = "summed"
[integrator]
integrated_phase = "imag"
n_start = 2500
n_increase = 2500
n_max = 25000
n_bins = 16
min_samples_for_update = 250
seed = 20260929
'"#,
            &format!("import graphs {} -p cff -i scalar", graph_path.display()),
            "generate existing -p cff -i scalar",
            "duplicate integrand -p cff -i scalar --output_process_name ltd --output_integrand_name scalar",
        ],
    )?;
    // The scalar box is UV finite, but both the local and integrated threshold
    // counterterms are essential here. These two slots differ only by 3D mode.
    for (process, representation) in [RepresentationMode::Cff, RepresentationMode::Ltd]
        .into_iter()
        .enumerate()
    {
        for level in &mut cli
            .state
            .process_list
            .get_integrand_mut(process, "scalar")?
            .get_mut_settings()
            .stability
            .levels
        {
            level.three_dimensional_representation = Some(representation);
        }
    }
    // Table 1 (printed page 61): both phases share the 10^-8 exponent.
    let reference = Complex::new(F(-6.57830e-8), F(-7.43707e-8));
    let mut command = scalar_topology_integrate_command(
        test_name,
        "correlated",
        &[("cff", "scalar"), ("ltd", "scalar")],
        &[("cff@scalar", reference), ("ltd@scalar", reference)],
    );
    command.n_cores = Some(3);
    command.no_stream_iterations = true;
    command.no_stream_updates = true;
    let result = command.run(&mut cli.state, &cli.cli_settings)?;
    assert_eq!(result.slots.len(), 2);
    for (process, name) in ["cff", "ltd"].into_iter().enumerate() {
        let ProcessIntegrand::Amplitude(amplitude) = cli
            .state
            .process_list
            .get_integrand_mut(process, "scalar")?
        else {
            panic!("expected the scalar Box4E amplitude")
        };
        let [graph] = amplitude.data.graph_terms.as_slice() else {
            panic!("expected one box graph")
        };
        // These are the warmed owners used by Integrate. Four maximal two-surface
        // groups, each with a nonempty complement, exclude a common center.
        // Their prefactor evaluators implement E-surface multichanneling in the
        // counterterm group loop, independently of the disabled LMB sampling above.
        let overlap = &graph.threshold_counterterm.overlap;
        assert_eq!(overlap.existing_esurfaces.len(), 4);
        assert_eq!(overlap.overlap_groups.len(), 4);
        let mut groups = std::collections::BTreeSet::new();
        let mut centers = Vec::new();
        for group in &overlap.overlap_groups {
            assert_eq!(group.existing_esurfaces.len(), 2);
            assert_eq!(group.complement.len(), 2);
            assert!(
                group
                    .prefactor_evaluator
                    .as_ref()
                    .is_some_and(|e| !e.is_empty())
            );
            let mut ids = group
                .existing_esurfaces
                .iter()
                .map(|id| usize::from(*id))
                .collect::<Vec<_>>();
            ids.sort_unstable();
            assert!(groups.insert(ids), "overlap channels must be distinct");
            let [center] = group.center.0.as_slice() else {
                panic!("Box4E has one loop momentum")
            };
            let coordinates = [center.px.0, center.py.0, center.pz.0];
            assert!(coordinates.iter().all(|coordinate| coordinate.is_finite()));
            centers.push(coordinates);
        }
        for (index, left) in centers.iter().enumerate() {
            for right in &centers[..index] {
                let separation = left
                    .iter()
                    .zip(right)
                    .map(|(a, b)| (a - b).powi(2))
                    .sum::<f64>()
                    .sqrt();
                assert!(
                    separation > 1.0e-6 * 100.0,
                    "Box4E centers must differ on the E_cm=100 scale: {centers:?}",
                );
            }
        }

        let slot = result.slot(&format!("{name}@scalar")).unwrap();
        assert_eq!(slot.integration_statistics.nan_or_unstable_percentage, 0.0);
        let estimate = &slot.integral;
        assert_eq!(estimate.neval, 25_000);
        // Five decimal places in the tabulated mantissa give 5e-14 rounding.
        for (value, error, reference) in [
            (estimate.result.re.0, estimate.error.re.0, reference.re.0),
            (estimate.result.im.0, estimate.error.im.0, reference.im.0),
        ] {
            assert!(value.is_finite() && error.is_finite() && error > 0.0);
            assert!(error < 0.10 * reference.abs(), "{name}: {estimate:?}");
            assert!(
                (value - reference).abs() < 5.0 * error + 5.0e-14,
                "{name}: {value:e} +/- {error:e}, published {reference:e}",
            );
        }
    }
    let cff = &result.slot("cff@scalar").unwrap().integral;
    let ltd = &result.slot("ltd@scalar").unwrap().integral;
    for (left, right, left_error, right_error) in [
        (
            cff.result.re.0,
            ltd.result.re.0,
            cff.error.re.0,
            ltd.error.re.0,
        ),
        (
            cff.result.im.0,
            ltd.result.im.0,
            cff.error.im.0,
            ltd.error.im.0,
        ),
    ] {
        assert!(
            (left - right).abs() < 5.0 * left_error.hypot(right_error),
            "Box4E CFF/LTD integral disagreement: {cff:?}, {ltd:?}",
        );
    }
    clean_test(&root);
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
