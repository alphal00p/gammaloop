use super::utils::*;
use super::*;

#[test]
#[serial]
fn box_4e_multichannel_subtraction_reproduces_published_integral() -> Result<()> {
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
[subtraction.overlap_settings]
enable_heuristics = false
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
#[serial]
fn overlap_center_objectives_preserve_the_integrated_amplitude() -> Result<()> {
    use gammalooprs::{
        integrands::process::ProcessIntegrand, settings::runtime::OverlapCenterObjective,
    };

    let test_name = "overlap_center_objectives";
    let root = get_tests_workspace_path().join(test_name);
    let mut cli = get_test_cli(None, &root, Some(test_name.to_owned()), true)?;
    let graph_path = root.join("triangle_box.dot");
    // The UV-finite triangle-box topology has an asymmetric mass assignment.
    // Its off-shell external legs also avoid massless soft/collinear divergences.
    // Unequal routing norms distinguish energy depth from geometric clearance;
    // massless focal segments leave a genuine optimal face for the signed sum.
    std::fs::write(
        &graph_path,
        r#"digraph triangle_box {
            num=1; overall_factor=1;
            node [num=1]; edge [num=1];
            p1 [style=invis]; p2 [style=invis]; p3 [style=invis];
            p1 -> A:0 [id=0, particle=scalar_0];
            C:1 -> p2 [id=1, particle=scalar_0];
            D:2 -> p3 [id=2, particle=scalar_0];
            A -> B [id=3, particle=scalar_0];
            B -> C [id=4, particle=scalar_1];
            C -> D [id=5, particle=scalar_0];
            D -> E [id=6, particle=scalar_2, lmb_id=0];
            B -> E [id=7, particle=scalar_0, lmb_id=1];
            E -> A [id=8, particle=scalar_1];
        }"#,
    )?;
    run_commands(
        &mut cli,
        &[
            "import model scalars-default.json",
            "set model mass_scalar_1=0.5 mass_scalar_2=1.0",
            "set global kv global.generation.evaluator.compile=false global.generation.evaluator.iterative_orientation_optimization=false global.generation.evaluator.store_atom=false global.generation.evaluator.summed=false global.generation.evaluator.summed_function_map=true global.generation.uv.subtract_uv=false global.generation.uv.generate_integrated=false global.generation.threshold_subtraction.enable_thresholds=true global.generation.threshold_subtraction.check_esurface_at_generation=true global.generation.tropical_subgraph_table.disable_tropical_generation=true",
            r#"set default-runtime string '
[general]
evaluator_method = "SummedFunctionMap"
integral_unit = "none"
disable_flux_factor = true
[kinematics]
e_cm = 2.5
[kinematics.externals]
type = "constant"
[kinematics.externals.data]
momenta = [[2.5,1.0,0.0,0.0],[3.0,2.0,0.0,0.0],"dependent"]
helicities = [0,0,0]
[subtraction]
disable_threshold_subtraction = false
[subtraction.overlap_settings]
enable_heuristics = false
[stability]
rotation_axis = []
[sampling]
graphs = "summed"
orientations = "summed"
lmb_multichanneling = false
lmb_channels = "summed"
[integrator]
integrated_phase = "imag"
n_start = 2500
n_increase = 2500
n_max = 75000
n_bins = 16
min_samples_for_update = 250
seed = 20260927
'"#,
            &format!(
                "import graphs {} -p triangle_box -i scalar",
                graph_path.display()
            ),
            "generate existing -p triangle_box -i scalar",
            "duplicate integrand -p triangle_box -i scalar --output_process_name chebyshev --output_integrand_name scalar",
            "duplicate integrand -p triangle_box -i scalar --output_process_name min_sum --output_integrand_name scalar",
        ],
    )?;

    let objectives = [
        OverlapCenterObjective::MaxMinDepth,
        OverlapCenterObjective::RelaxedChebyshev,
        OverlapCenterObjective::MinSum,
    ];
    for (process, objective) in objectives.into_iter().enumerate() {
        cli.state
            .process_list
            .get_integrand_mut(process, "scalar")?
            .get_mut_settings()
            .subtraction
            .overlap_settings
            .objective = objective;
    }
    // The existing correlated integration owner adapts a single grid and evaluates
    // all three amplitudes on exactly the same samples, with the same seed/settings.
    let mut command = scalar_topology_integrate_command(
        test_name,
        "correlated",
        &[
            ("triangle_box", "scalar"),
            ("chebyshev", "scalar"),
            ("min_sum", "scalar"),
        ],
        &[],
    );
    command.n_cores = Some(3);
    command.no_stream_iterations = true;
    command.no_stream_updates = true;
    let result = command.run(&mut cli.state, &cli.cli_settings)?;
    assert_eq!(result.slots.len(), 3);
    let mut centers = Vec::new();
    let mut estimates = Vec::new();
    for (process, name) in ["triangle_box", "chebyshev", "min_sum"]
        .into_iter()
        .enumerate()
    {
        // Integrate warmed these exact owners before cloning its workers. Reading
        // them after the run checks centers actually used, without a second solve.
        let integrand = cli
            .state
            .process_list
            .get_integrand_mut(process, "scalar")?;
        let ProcessIntegrand::Amplitude(amplitude) = integrand else {
            panic!("expected the scalar triangle-box amplitude")
        };
        let [graph] = amplitude.data.graph_terms.as_slice() else {
            panic!("expected one graph")
        };
        let overlap = &graph.threshold_counterterm.overlap;
        assert_eq!(overlap.existing_esurfaces.len(), 6);
        let [group] = overlap.overlap_groups.as_slice() else {
            panic!("all six thresholds must have one common overlap")
        };
        assert_eq!(group.existing_esurfaces.len(), 6);
        assert!(group.complement.is_empty());
        centers.push(
            group
                .center
                .0
                .iter()
                .flat_map(|momentum| [momentum.px.0, momentum.py.0, momentum.pz.0])
                .collect::<Vec<_>>(),
        );

        let estimate = result
            .slot(&format!("{name}@scalar"))
            .expect("each objective must have an integration slot")
            .integral
            .clone();
        assert_eq!(estimate.neval, 75_000);
        assert!(estimate.result.re.0.is_finite() && estimate.result.im.0.is_finite());
        estimates.push(estimate);
    }

    // This short regression resolves the dominant imaginary part to 15%. The
    // smaller real part has a separate absolute error budget on the same physical
    // amplitude scale; neither error is allowed to hide a discrepancy in the other.
    let amplitude_scale = estimates
        .iter()
        .map(|estimate| estimate.result.im.0.abs())
        .sum::<f64>()
        / 3.0;
    assert!(amplitude_scale > 1.0e-6);
    for estimate in &estimates {
        assert!(
            estimate.error.im.0.is_finite()
                && estimate.error.im.0 > 0.0
                && estimate.error.im.0 < 0.15 * estimate.result.im.0.abs(),
            "imaginary precision: {estimate:?}"
        );
        assert!(
            estimate.error.re.0.is_finite()
                && estimate.error.re.0 > 0.0
                && estimate.error.re.0 < 0.20 * amplitude_scale,
            "real precision: {estimate:?}"
        );
    }

    // Independent geometry: energy depth fixes k_x=-2/3. The optimal relaxed
    // radius is 1/(2+sqrt(2)), requiring k_x<=-3/4; its flat face admits the
    // signed-sum optimum k_x=-3/4. Do not compare solver-dependent free coordinates.
    assert!((centers[0][0] + 2.0 / 3.0).abs() < 1.0e-4);
    assert!(centers[1][0] <= -0.75 + 1.0e-5);
    assert!((centers[2][0] + 0.75).abs() < 1.0e-4);
    let optimal_radius = 1.0 / (2.0 + 2.0_f64.sqrt());
    for (index, center) in centers.iter().enumerate() {
        let k = center[0];
        let l = center[3];
        // P=|k+l+1|+sqrt(1+k²)+|l|-2.5, with L_P=2+sqrt(2).
        let radius = (2.5 - (k + l + 1.0).abs() - (1.0 + k * k).sqrt() - l.abs()) * optimal_radius;
        if index == 0 {
            assert!(radius < optimal_radius - 0.005);
        } else {
            assert!(radius >= optimal_radius - 1.0e-6);
        }
    }
    for left in 0..3 {
        for right in 0..left {
            let separation = centers[left]
                .iter()
                .zip(&centers[right])
                .map(|(a, b)| (a - b).powi(2))
                .sum::<f64>()
                .sqrt();
            assert!(
                separation > 0.0025,
                "objectives must use meaningfully different centers: {centers:?}"
            );
            for (a, b, a_error, b_error) in [
                (
                    estimates[left].result.re.0,
                    estimates[right].result.re.0,
                    estimates[left].error.re.0,
                    estimates[right].error.re.0,
                ),
                (
                    estimates[left].result.im.0,
                    estimates[right].result.im.0,
                    estimates[left].error.im.0,
                    estimates[right].error.im.0,
                ),
            ] {
                assert!(
                    (a - b).abs() <= 5.0 * a_error.hypot(b_error),
                    "center objectives changed the integral: {estimates:?}"
                );
            }
        }
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
    fn vakint_compare_scalar_bubble() -> Result<()> {
        let mut cli = get_test_cli(
            Some("scalar_bubble.toml".into()),
            get_tests_workspace_path().join("scalar_bubble"),
            Some("scalar_bubble".to_string()),
            false,
        )?;

        cli.run_command("run import_graph")?;

        cli.run_command("run integrated")?;
        // These historical targets describe the local-only subtraction at m_uv=5.
        // Enabling the integrated addback instead gives the MS-bar result at mu_r;
        // its reference scale must be reconciled before this failing case is restored.

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
        assert!(
            integral_no_cache.is_compatible_with_target(Complex::new(F(0.0), F(-2.03838e-02)), 2)
        );

        cli.run_command("set model mass_scalar_1=2.0")?;
        let integral_no_cache = integrate_command.run(&mut cli.state, &cli.cli_settings)?;
        assert!(
            integral_no_cache.is_compatible_with_target(Complex::new(F(0.0), F(-1.16050e-02)), 2)
        );

        cli.run_command("set model mass_scalar_1=3.0")?;
        let integral_no_cache = integrate_command.run(&mut cli.state, &cli.cli_settings)?;
        assert!(
            integral_no_cache.is_compatible_with_target(Complex::new(F(0.0), F(-6.46968e-03)), 2)
        );
        let renorm_command = Renormalize::default();

        let res = renorm_command.run(&mut cli.state, &cli.cli_settings)?;

        println!("{}", res[0]);

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
        // from Kaapo: m=1 muv=5 4.37688e-03 m=2 muv=5 	2.48100e-03	 m=3 muv=5 1.07231e-03
        cli.run_command("set model mass_scalar_1=1.0")?;
        let res = profile_cmd
            .run(&mut cli.state, &cli.cli_settings)?
            .unwrap_uv();
        assert_eq!(res.pass_fail(-0.9).failed, 0);
        let integral_no_cache = integrate_command.run(&mut cli.state, &cli.cli_settings)?;
        assert!(
            integral_no_cache.is_compatible_with_target(Complex::new_re(F(-4.37688e-03)), 1) // .result
                                                                                             // .approx_eq(&Complex::new_re(F(-4.37688e-03)), &F(0.01))
        );

        // assert_snapshot!(format!("{integral_no_cache:.3}"),@"-4.359e-3");

        cli.run_command("set model mass_scalar_1=2.0")?;
        let res = profile_cmd
            .run(&mut cli.state, &cli.cli_settings)?
            .unwrap_uv();
        assert_eq!(res.pass_fail(-0.9).failed, 0);
        let integral_no_cache = integrate_command.run(&mut cli.state, &cli.cli_settings)?;
        // assert_snapshot!(format!("{integral_no_cache:.3}"),@"-2.474e-3");
        assert!(integral_no_cache.is_compatible_with_target(Complex::new_re(F(-2.48100e-03)), 1));

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

        let profile_cmd = Profile::UltraViolet(UltraVioletProfile {
            use_f128: false,
            max_scale_exponent: 4.,
            min_scale_exponent: 2.0,
            analyse_analytically: false,
            ..Default::default()
        });

        //1.47240e-03	7.15184e-04	2.27485e-04

        // from Kaapo: m=1 muv=5 1.47240e-03 m=2 muv=5 	7.15184e-04	 m=3 muv=5 2.27485e-04
        // Four propagator i factors and three loops convert the historical
        // Euclidean targets below to the physical -i axis. Their subdivergence
        // prescription and p^2=1 magnitudes remain unverified; this phase migration
        // does not certify the old baselines for the current subtraction settings.
        cli.run_command("set model mass_scalar_1=1.0")?;
        let res = profile_cmd
            .run(&mut cli.state, &cli.cli_settings)?
            .unwrap_uv();
        assert_eq!(res.pass_fail(-0.9).failed, 0);
        let integral_no_cache = integrate_command.run(&mut cli.state, &cli.cli_settings)?;
        assert!(
            integral_no_cache.is_compatible_with_target(Complex::new(F(0.0), F(-1.47240e-03)), 1),
            "Not compatible: {integral_no_cache}",
        );

        // cli.run_command("set model mass_scalar_1=2.0")?;
        // let res = profile_cmd.run(&mut cli.state, &cli.cli_settings)?;
        // assert_eq!(res.pass_fail(-0.9).failed, 0);
        // let integral_no_cache = integrate_command.run(&mut cli.state, &cli.cli_settings)?;
        // assert!(
        //     integral_no_cache.is_compatible_with_target(Complex::new(F(0.0), F(-7.15184e-04)), 3),
        //     "Not compatible: {integral_no_cache}",
        // );

        // cli.run_command("set model mass_scalar_1=3.0")?;
        // let res = profile_cmd.run(&mut cli.state, &cli.cli_settings)?;
        // assert_eq!(res.pass_fail(-0.9).failed, 0);
        // let integral_no_cache = integrate_command.run(&mut cli.state, &cli.cli_settings)?;
        // assert!(
        //     integral_no_cache.is_compatible_with_target(Complex::new(F(0.0), F(-2.27485e-04)), 3),
        //     "Not compatible: {integral_no_cache}",
        // );

        // clean_test(&cli.cli_settings.state.folder);

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

        // At these subthreshold inputs the four-vertex/four-propagator scalar
        // box has positive imaginary phase. Preserve the historical magnitude
        // on that axis; it still needs an independent pointwise certificate.
        assert_snapshot!(format!("{a:.8e}"),@"(0e0+1.3485885914334373e-3i)");

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
        //2.90078e-06	1.59168e-06	6.86001e-07

        // from Kaapo: m=1 muv=5 2.90078e-06 m=2 muv=5 	1.59168e-06	 m=3 muv=5 6.86001e-07
        // These historical targets have an unresolved input mismatch: the card
        // uses m_uv=1 and p^2=-12. Seven propagator i factors and four loops
        // convert the historical Euclidean targets to the physical +i axis;
        // this does not certify their signs or magnitudes for the current inputs.
        cli.run_command("set model mass_scalar_1=1.0")?;
        let res = profile_cmd
            .run(&mut cli.state, &cli.cli_settings)?
            .unwrap_uv();
        assert_eq!(res.pass_fail(-0.9).failed, 0);

        let integral_no_cache = integrate_command.run(&mut cli.state, &cli.cli_settings)?;
        assert!(
            integral_no_cache.is_compatible_with_target(Complex::new(F(0.0), F(2.90078e-06)), 1),
            "Not compatible: {integral_no_cache}",
        );

        cli.run_command("set model mass_scalar_1=2.0")?;
        let integral_no_cache = integrate_command.run(&mut cli.state, &cli.cli_settings)?;
        assert!(
            integral_no_cache.is_compatible_with_target(Complex::new(F(0.0), F(1.59168e-06)), 3),
            "Not compatible: {integral_no_cache}",
        );

        cli.run_command("set model mass_scalar_1=3.0")?;
        let integral_no_cache = integrate_command.run(&mut cli.state, &cli.cli_settings)?;
        assert!(
            integral_no_cache.is_compatible_with_target(Complex::new(F(0.0), F(6.86001e-07)), 3),
            "Not compatible: {integral_no_cache}",
        );

        clean_test(&cli.cli_settings.state.folder);

        Ok(())
    }

    #[test]
    fn test_integrate_dotted_bubble() -> Result<()> {
        let target = Complex::new(F(0.0), F(0.002871504663657095));
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
