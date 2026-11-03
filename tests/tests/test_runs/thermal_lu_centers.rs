use super::*;
use ndarray::Array2;
use std::collections::BTreeMap;

/// Probe the spacelike double-triangle overlap of the LL and HH two-particle absorption
/// cuts. At `P=(0.5,0,0,1.5)`, `K0=s0 (sin 0.3, 0, cos 0.3)` lies on
/// `E_L(K0)-E_L(P-K0)=E` and `K1` on `E_H(P-K1)-E_H(K1)=E`, so both cuts vanish at `t=1`;
/// off the probe axis no ray is degenerate. Origin flow crosses the two sheets with opposite
/// slopes and leaves the overlap pole uncancelled; the searched center orients both and keeps
/// the complete cut sum finite. Scaling `K0` moves its root away from the other one.
#[test]
fn thermal_lu_center_search_cancels_the_spacelike_double_triangle_overlap() -> Result<()> {
    let _serial_guard = gammaloop_integration_tests::SERIAL_TEST_LOCK.lock();
    let state_path =
        get_tests_workspace_path().join(format!("thermal_lu_centers_{}", std::process::id()));
    let mut cli = get_example_cli(
        "thermal_rates/lu_centers/run.toml",
        &["generate"],
        Some(state_path),
        None,
        true,
    )?;
    let (s0, angle) = (1.0870138095413315, 0.3_f64);
    let k1 = [0.0, 0.22539501322019227, 0.41258280423477773];
    let displacements = [1e-3, -1e-3, 1e-6, -1e-6];
    let points = Array2::from_shape_fn((displacements.len(), 6), |(row, column)| {
        let scale = s0 * (1.0 + displacements[row]);
        match column {
            0 => scale * angle.sin(),
            2 => scale * angle.cos(),
            3..=5 => k1[column - 3],
            _ => 0.0,
        }
    });
    let discrete = Array2::<usize>::zeros((displacements.len(), 0));
    let mut values = BTreeMap::new();
    for centers in ["origin", "search"] {
        cli.run_command(&format!(
            "set process -p dt -i NLO kv kinematics.externals.data.momenta=[[0.5,0.0,0.0,1.5]] kinematics.e_cm={} lu_flow.centers={centers}",
            2.0_f64.sqrt()
        ))?;
        let (process_id, integrand_name) = cli.state.find_integrand_ref(
            Some(&ProcessRef::Unqualified("dt".into())),
            Some(&"NLO".into()),
        )?;
        let result = EvaluateSamples {
            process_id: Some(process_id),
            integrand_name: Some(integrand_name),
            use_arb_prec: true,
            minimal_output: true,
            return_generated_events: Some(false),
            momentum_space: true,
            points: points.view(),
            integrator_weights: None,
            discrete_dims: Some(discrete.view()),
            graph_names: None,
            orientations: None,
        }
        .run(&mut cli.state)?;
        let sums = result
            .samples
            .iter()
            .map(|sample| sample.evaluation.integrand_result.re.0)
            .collect::<Vec<_>>();
        assert!(
            sums.iter().all(|value| value.is_finite()),
            "{centers}: {sums:?}"
        );
        gammalooprs::debug_tags!(#integration, #thermal, #benchmark, #summary;
            stage = "thermal_lu_center_overlap_probe", centers, sums = ?sums,
            "Complete double-triangle cut sums approaching the spacelike overlap"
        );
        values.insert(centers, sums);
    }
    // The uncancelled origin pole grows like 1/delta and changes sign across the overlap.
    let origin = &values["origin"];
    assert!(
        origin[2].abs() > 300.0 * origin[0].abs(),
        "origin {origin:?}"
    );
    assert!(origin[2] * origin[3] < 0.0, "origin {origin:?}");
    // With the searched center both sides approach the same finite limit.
    let search = &values["search"];
    let scale = search[0].abs().max(search[1].abs());
    assert!(
        (search[2] - search[3]).abs() < 1e-2 * scale,
        "search {search:?}"
    );
    assert!(search[2].abs() < 10.0 * scale, "search {search:?}");
    clean_test(cli.cli_settings.state.folder);
    Ok(())
}
