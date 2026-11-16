use color_eyre::Result;
use gammaloop_api::commands::Profile;
use gammaloop_api::commands::profile::{InfraRedProfile, UltraVioletProfile};
use gammaloop_api::state::{ProcessRef, SyncSettings};
use gammaloop_integration_tests as integration_tests;
use gammaloop_integration_tests::CLIState;
use gammalooprs::processes::ProcessCollection;
use gammalooprs::uv::ApproximationType;
use gammalooprs::uv::profile::{UV_PROFILE_MAX_DOD, UVLimitSelection};
use gammalooprs::uv::uv_graph::UVE;
use linnet::half_edge::involution::EdgeIndex;

// A logarithmic endpoint can fit a tiny positive exponent on a finite ray.
// These limits have integer power counting, so require a margin from zero.
const MIN_SOFT_SCALING: f64 = 0.5;

#[test]
fn thermal_quark_bubble_soft_ct_preserves_uv_finiteness() -> Result<()> {
    let _serial = integration_tests::SERIAL_TEST_LOCK.lock();
    let name = "thermal_quark_bubble";
    let state_path = integration_tests::get_tests_workspace_path().join(name);
    let mut cli = integration_tests::get_test_cli(
        Some(format!("{name}_soft_ir.toml").into()),
        &state_path,
        None,
        true,
    )?;
    assert!(!cli.cli_settings.global.generation.uv.generate_integrated);
    assert_eq!(
        cli.cli_settings
            .global
            .generation
            .uv
            .renormalization_prescription
            .massless_power_divergent,
        ApproximationType::IR,
    );
    cli.run_command("run generate")?;

    for mu_b in [0.0, 3.0] {
        cli.run_command(&format!("set model muB={mu_b}"))?;
        // This two-loop graph isolates the UV regression caused by extending
        // a thermal soft subtraction away from its soft external momenta.
        // All three proper divergent cycle unions and the all-hard limit
        // must remain finite, including the alternative physical LMB.
        assert_thermal_uv_finiteness(&mut cli, name, mu_b, 4)?;

        let report = Profile::InfraRed(InfraRedProfile {
            process: Some(ProcessRef::Unqualified(name.to_owned())),
            integrand_name: Some(name.to_owned()),
            n_points: 15,
            min_scale_exponent: -2.0,
            max_scale_exponent: -5.0,
            seed: Some(1337),
            select: Some(format!("{name} S(e0)")),
            ..Default::default()
        })
        .run(&mut cli.state, &cli.cli_settings)?
        .unwrap_ir();
        assert_eq!(report.results_per_graph.len(), 1, "{report}");
        let graph = &report.results_per_graph[0];
        assert_eq!(graph.graph_name, name);
        assert_eq!(graph.single_limit_reports.len(), 1, "{graph}");
        let soft = &graph.single_limit_reports[0];
        assert_eq!(soft.limit_name, "S(e0)");
        assert!(soft.orientation_label.is_none());
        assert!(
            soft.passed
                && soft.scaling.is_finite()
                && soft.scaling > MIN_SOFT_SCALING
                && soft.power_law_fit.r_squared() >= 0.98,
            "the bubble's soft remainder is not integrable at muB={mu_b}: {soft}"
        );
    }

    integration_tests::clean_test(state_path);
    Ok(())
}

#[test]
fn thermal_quark_ring_soft_ct_cancels_bose_enhancement() -> Result<()> {
    assert_thermal_ring_soft_cancellation("thermal_quark_ring", 1)
}

#[test]
fn thermal_gluon_quark_ring_soft_ct_cancels_bose_enhancement() -> Result<()> {
    assert_thermal_ring_soft_cancellation("thermal_gluon_quark_ring", 21)
}

fn assert_thermal_ring_soft_cancellation(name: &str, second_bubble_pdg: isize) -> Result<()> {
    let _serial = integration_tests::SERIAL_TEST_LOCK.lock();
    let mut profiles = Vec::new();

    for prescription in [ApproximationType::MUV, ApproximationType::IR] {
        let state_path =
            integration_tests::get_tests_workspace_path().join(format!("{name}_{prescription:?}"));
        let mut cli = integration_tests::get_test_cli(
            Some(format!("{name}_soft_ir.toml").into()),
            &state_path,
            None,
            true,
        )?;
        assert!(!cli.cli_settings.global.generation.uv.generate_integrated);
        cli.cli_settings
            .global
            .generation
            .uv
            .renormalization_prescription
            .massless_power_divergent = prescription;
        cli.cli_settings.sync_settings()?;
        cli.run_command("run generate")?;

        let process_id = cli
            .state
            .resolve_process_ref(Some(&ProcessRef::Unqualified(name.to_owned())))?;
        let ProcessCollection::Amplitudes(amplitudes) =
            &cli.state.process_list.processes[process_id].collection
        else {
            panic!("the thermal ring must be a vacuum amplitude");
        };
        let amplitude = amplitudes.get(name).expect("the ring was generated");
        assert_eq!(amplitude.graphs.len(), 1);
        let graph = &amplitude.graphs[0].graph;
        assert_eq!(
            graph
                .loop_momentum_basis
                .loop_edges
                .iter()
                .copied()
                .map(usize::from)
                .collect::<Vec<_>>(),
            [0, 2, 4],
            "the ring must retain [q, p, k] so S(e0) scales both gluon carriers"
        );
        assert_eq!(
            graph.loop_momentum_basis.edge_signatures[EdgeIndex::from(0)],
            graph.loop_momentum_basis.edge_signatures[EdgeIndex::from(1)],
            "the two ring gluons must carry the same q"
        );
        for edge in [4, 5] {
            assert_eq!(
                graph[EdgeIndex::from(edge)]
                    .particle_pdg_code()
                    .expect("ring edges retain their model particles"),
                second_bubble_pdg,
                "{name} must contain the requested second polarization loop"
            );
        }

        let profile = |limit| {
            Profile::InfraRed(InfraRedProfile {
                process: Some(ProcessRef::Unqualified(name.to_owned())),
                integrand_name: Some(name.to_owned()),
                n_points: 8,
                min_scale_exponent: -2.0,
                max_scale_exponent: -5.0,
                seed: Some(1337),
                select: Some(format!("{name} {limit}")),
                ..Default::default()
            })
        };
        let mut chemical_potential_profiles = Vec::new();
        for mu_b in [0.0, 3.0] {
            cli.run_command(&format!("set model muB={mu_b}"))?;
            let mut report = profile("S(e0)")
                .run(&mut cli.state, &cli.cli_settings)?
                .unwrap_ir();
            assert_eq!(report.results_per_graph.len(), 1, "{report}");
            let mut graph = report.results_per_graph.pop().expect("one ring report");
            assert_eq!(graph.graph_name, name);
            assert_eq!(graph.single_limit_reports.len(), 1, "{graph}");
            let soft = graph
                .single_limit_reports
                .pop()
                .expect("one selected soft limit");
            assert_eq!(soft.limit_name, "S(e0)");
            assert!(soft.orientation_label.is_none());
            assert!(
                soft.scaling.is_finite() && soft.power_law_fit.r_squared() >= 0.98,
                "unresolved {prescription:?} {name} soft fit at muB={mu_b}: {soft}"
            );
            chemical_potential_profiles.push((mu_b, soft));

            if second_bubble_pdg == 21 && prescription == ApproximationType::IR {
                // The quark loop stays hard while both gluon loop momenta
                // vanish. This checks uniformity beyond the single-q limit.
                let joint = profile("S(e0)S(e4)")
                    .run(&mut cli.state, &cli.cli_settings)?
                    .unwrap_ir();
                assert_eq!(joint.results_per_graph.len(), 1, "{joint}");
                let graph = &joint.results_per_graph[0];
                assert_eq!(graph.single_limit_reports.len(), 1, "{graph}");
                let soft = &graph.single_limit_reports[0];
                assert_eq!(soft.limit_name, "S(e0)S(e4)");
                assert!(soft.orientation_label.is_none());
                assert!(
                    soft.passed
                        && soft.scaling.is_finite()
                        && soft.scaling > MIN_SOFT_SCALING
                        && soft.power_law_fit.r_squared() >= 0.98,
                    "the mixed ring's joint soft remainder is not integrable at muB={mu_b}: {soft}"
                );
            }

            if prescription == ApproximationType::IR {
                assert_thermal_uv_finiteness(&mut cli, name, mu_b, 8)?;
            }
        }
        profiles.push(chemical_potential_profiles);
        integration_tests::clean_test(state_path);
    }

    for ((muv_mu, muv), (ir_mu, ir)) in profiles[0].iter().zip(&profiles[1]) {
        assert_eq!(muv_mu, ir_mu);
        assert!(muv.ray_fingerprint.is_some());
        assert_eq!(
            muv.ray_fingerprint, ir.ray_fingerprint,
            "MUV and IR must evaluate the identical {name} ray at muB={ir_mu}"
        );
        // The profiler includes the three-dimensional soft volume factor.
        // The control must retain a divergence; the local thermal subtraction
        // must cross the integrability boundary without an integrated addback.
        assert!(
            muv.scaling < -0.5,
            "the MUV control lost the {name} soft divergence at muB={muv_mu}: {muv}"
        );
        assert!(
            ir.scaling > MIN_SOFT_SCALING && ir.passed,
            "the {name} soft remainder is not integrable at muB={ir_mu}: {ir}"
        );
        assert!(
            ir.scaling - muv.scaling > 1.0,
            "thermal IR subtraction did not improve {name} at muB={ir_mu}: MUV={muv}; IR={ir}"
        );
    }
    Ok(())
}

fn assert_thermal_uv_finiteness(
    cli: &mut CLIState,
    name: &str,
    mu_b: f64,
    expected_limits: usize,
) -> Result<()> {
    let report = Profile::UltraViolet(UltraVioletProfile {
        process: Some(ProcessRef::Unqualified(name.to_owned())),
        integrand_name: Some(name.to_owned()),
        selected_limits: UVLimitSelection::OnlyDivergent,
        n_points: 8,
        min_scale_exponent: 3.0,
        max_scale_exponent: 6.0,
        seed: Some(1337),
        use_f128: true,
        // The vacuum-subtracted all-hard limit decays exponentially.
        allow_vanishing_missing_fits: true,
        ..Default::default()
    })
    .run(&mut cli.state, &cli.cli_settings)?
    .unwrap_uv();
    assert!(!report.stopped_early);
    assert_eq!(report.graphs.len(), 1);
    assert_eq!(report.graphs[0].graph_name, name);
    let verdict = report.pass_fail(UV_PROFILE_MAX_DOD);
    assert_eq!(verdict.total, expected_limits, "{verdict}");
    assert_eq!(
        verdict.failed, 0,
        "the thermal soft CT introduced a UV divergence in {name} at muB={mu_b}: {verdict}"
    );
    Ok(())
}
