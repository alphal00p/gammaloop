use super::*;
use gammaloop_api::state::RunHistory;
use gammaloop_integration_tests::workspace_root;
use std::{fs, process::Command};

fn is_generated_state_dir(path: &Path) -> bool {
    path.file_name()
        .and_then(|name| name.to_str())
        .is_some_and(|name| name.contains("state"))
}

fn is_ignored_old_card(path: &Path) -> bool {
    path.file_name()
        .and_then(|name| name.to_str())
        .is_some_and(|name| name.starts_with("OLD_"))
}

fn discover_toml_cards(root: &Path, cards: &mut Vec<PathBuf>) -> Result<()> {
    for entry in fs::read_dir(root)? {
        let entry = entry?;
        let path = entry.path();
        if path.is_dir() {
            if !is_generated_state_dir(&path) {
                discover_toml_cards(&path, cards)?;
            }
        } else if !is_ignored_old_card(&path)
            && path.extension().and_then(|extension| extension.to_str()) == Some("toml")
        {
            cards.push(path);
        }
    }
    Ok(())
}

fn estimate_delta_sigma(lhs: &IntegralEstimate, rhs: &IntegralEstimate) -> f64 {
    let delta_re = lhs.result.re.0 - rhs.result.re.0;
    let delta_im = lhs.result.im.0 - rhs.result.im.0;
    let delta_norm = delta_re.hypot(delta_im);
    let error_norm = lhs
        .error
        .re
        .0
        .hypot(rhs.error.re.0)
        .hypot(lhs.error.im.0.hypot(rhs.error.im.0));

    if error_norm == 0.0 {
        if delta_norm == 0.0 {
            0.0
        } else {
            f64::INFINITY
        }
    } else {
        delta_norm / error_norm
    }
}

fn integration_slot<'a>(
    result: &'a RuntimeIntegrationResult,
    key: &str,
) -> Result<&'a IntegralEstimate> {
    Ok(&result
        .slot(key)
        .ok_or_else(|| eyre::eyre!("Missing integration-result slot '{key}'"))?
        .integral)
}

fn uv_profile_json_failed_count(profile: &JsonValue, max_dod: f64) -> Result<usize> {
    let slope_fails = |analysis: Option<&JsonValue>| {
        analysis
            .and_then(|analysis| analysis.pointer("/result/slope"))
            .and_then(JsonValue::as_f64)
            .is_none_or(|slope| slope.is_nan() || slope > max_dod)
    };

    let mut failed = 0;
    let graphs = profile
        .get("graphs")
        .and_then(JsonValue::as_array)
        .ok_or_else(|| eyre::eyre!("UV profile JSON is missing graph entries"))?;
    for graph in graphs {
        let lmbs = graph
            .get("lmbs")
            .and_then(JsonValue::as_array)
            .ok_or_else(|| eyre::eyre!("UV profile JSON is missing LMB entries"))?;
        for lmb in lmbs {
            let subsets = lmb
                .get("subsets")
                .and_then(JsonValue::as_array)
                .ok_or_else(|| eyre::eyre!("UV profile JSON is missing subset entries"))?;
            for subset in subsets {
                let per_orientation_entries = subset
                    .get("per_orientation_inspect_entries")
                    .and_then(JsonValue::as_array);
                let per_orientation_failed = per_orientation_entries
                    .map(|entries| {
                        entries
                            .iter()
                            .filter(|entry| slope_fails(entry.get("analysis")))
                            .count()
                    })
                    .unwrap_or(0);
                let all_orientations_pass =
                    per_orientation_entries.is_some() && per_orientation_failed == 0;

                if slope_fails(subset.pointer("/analysis/inspect_level")) && !all_orientations_pass
                {
                    failed += 1;
                }
                failed += per_orientation_failed;
            }
        }
    }
    Ok(failed)
}

#[test]
fn all_example_toml_cards_are_loadable() -> Result<()> {
    let workspace_root = workspace_root();
    let mut cards = Vec::new();
    for examples_dir in ["examples/api", "examples/cli"] {
        discover_toml_cards(&workspace_root.join(examples_dir), &mut cards)?;
    }
    cards.sort();
    if cards.is_empty() {
        return Err(eyre::eyre!("No example TOML cards were discovered"));
    }

    let failures = cards
        .iter()
        .filter_map(|card| {
            RunHistory::load(card).err().map(|err| {
                let display_path = card.strip_prefix(&workspace_root).unwrap_or(card);
                format!("{}:\n{err:?}", display_path.display())
            })
        })
        .collect_vec();

    if !failures.is_empty() {
        return Err(eyre::eyre!(
            "Failed to load {} example TOML card(s):\n\n{}",
            failures.len(),
            failures.join("\n\n")
        ));
    }
    Ok(())
}

#[test]
fn test_scalar_bubble_example_cli() -> Result<()> {
    let state_path = get_tests_workspace_path().join("scalar_bubble_example");
    let cli = get_example_cli(
        "scalar_topologies/bubble.toml",
        &["generate"],
        Some(state_path.clone()),
        None,
        true,
    )?;

    assert!(
        !cli.state.process_list.processes.is_empty(),
        "No processes were generated"
    );
    assert_eq!(
        cli.state.model.name, "scalars",
        "Expected scalars model to be loaded"
    );
    assert!(
        state_path.join("state_manifest.toml").exists(),
        "Expected saved state manifest in {}",
        state_path.display()
    );

    clean_test(&state_path);
    Ok(())
}

#[test]
#[serial]
fn test_sunrise_2loop_dod2_example_uv() -> Result<()> {
    let state_path = get_tests_workspace_path().join("sunrise_2loop_dod2_example");
    clean_test(&state_path);
    let workspace_root = workspace_root();
    let make_gammaloop_command = || {
        if workspace_root.join("target/debug/gammaloop").is_file()
            || workspace_root.join("target/dev-optim/gammaloop").is_file()
            || workspace_root.join("target/release/gammaloop").is_file()
        {
            Command::new("./gammaloop")
        } else {
            let mut command = Command::new("cargo");
            command.args(["run", "-p", "gammaloop-api", "--bin", "gammaloop", "--"]);
            command
        }
    };

    let run_example_block = |clean_state: bool, block_name: &str| -> Result<()> {
        let mut command = make_gammaloop_command();
        if clean_state {
            command.arg("--clean-state");
        } else {
            command.arg("-n");
        }
        let output = command
            .current_dir(&workspace_root)
            .args(["-s"])
            .arg(&state_path)
            .args([
                "-l",
                "warn",
                "examples/cli/scalar_topologies/multi_scalar_integrands.toml",
                "run",
                block_name,
            ])
            .env("RUST_MIN_STACK", "67108864")
            .output()
            .map_err(|err| eyre::eyre!("Failed to run gammaloop block '{block_name}': {err}"))?;

        if !output.status.success() {
            return Err(eyre::eyre!(
                "gammaloop block '{block_name}' failed with status {}\nstdout:\n{}\nstderr:\n{}",
                output.status,
                String::from_utf8_lossy(&output.stdout),
                String::from_utf8_lossy(&output.stderr)
            ));
        }
        Ok(())
    };
    let run_uv_profile = |process: &str, output_dir: &Path| -> Result<()> {
        let _ = fs::remove_dir_all(output_dir);
        let output = make_gammaloop_command()
            .current_dir(&workspace_root)
            .arg("-n")
            .args(["-s"])
            .arg(&state_path)
            .args([
                "-l",
                "warn",
                "examples/cli/scalar_topologies/multi_scalar_integrands.toml",
                "profile",
                "ultra-violet",
                "-p",
                process,
                "-i",
                "scalar_sunrise_below_thres",
                "--n-points",
                "100",
                "--min-scaling",
                "4",
                "--max-scaling",
                "5",
                "--per-orientation",
                "-o",
            ])
            .arg(output_dir)
            .env("RUST_MIN_STACK", "67108864")
            .output()
            .map_err(|err| eyre::eyre!("Failed to run UV profile for '{process}': {err}"))?;

        if !output.status.success() {
            return Err(eyre::eyre!(
                "UV profile for '{process}' failed with status {}\nstdout:\n{}\nstderr:\n{}",
                output.status,
                String::from_utf8_lossy(&output.stdout),
                String::from_utf8_lossy(&output.stderr)
            ));
        }
        Ok(())
    };

    run_example_block(true, "generate")?;

    for process in ["sunrise_2loop_dod2", "sunrise_2loop_dod2_no_integrated_UV"] {
        let output_dir = state_path.join(format!("uv_profile_{process}"));
        run_uv_profile(process, &output_dir)?;
        let uv: JsonValue =
            serde_json::from_reader(fs::File::open(output_dir.join("uv_profile.json"))?)?;
        assert_eq!(
            uv_profile_json_failed_count(&uv, -0.9)?,
            0,
            "{process} UV profile should be locally convergent"
        );
    }

    run_example_block(false, "integrate_sunrise_2loop_dod2_muv_independence")?;

    let result = super::utils::load_integration_result(
        &state_path.join("integration_workspace/integration_result.json"),
    )?;
    let integrated = integration_slot(&result, "sunrise_2loop_dod2@scalar_sunrise_below_thres")?;
    let integrated_muv3 = integration_slot(
        &result,
        "sunrise_2loop_dod2@scalar_sunrise_below_thres_muv3",
    )?;
    let no_integrated = integration_slot(
        &result,
        "sunrise_2loop_dod2_no_integrated_UV@scalar_sunrise_below_thres",
    )?;
    let no_integrated_muv3 = integration_slot(
        &result,
        "sunrise_2loop_dod2_no_integrated_UV@scalar_sunrise_below_thres_muv3",
    )?;

    let integrated_sigma = estimate_delta_sigma(integrated, integrated_muv3);
    assert!(
        integrated_sigma <= 2.0,
        "integrated result should be m_uv independent within 2 sigma, got {integrated_sigma:.2} sigma: {integrated} vs {integrated_muv3}"
    );

    let no_integrated_sigma = estimate_delta_sigma(no_integrated, no_integrated_muv3);
    assert!(
        no_integrated_sigma >= 5.0,
        "no-integrated result should show visible m_uv dependence, got {no_integrated_sigma:.2} sigma: {no_integrated} vs {no_integrated_muv3}"
    );

    clean_test(&state_path);
    Ok(())
}

#[test]
fn test_aa_aa_example_card() -> Result<()> {
    let run_history = example_run_card("aa_aa/1L/aa_aa.toml")?;

    assert_eq!(
        run_history.command_blocks[0].commands[0]
            .raw_string
            .as_deref(),
        Some("import model sm-default.json")
    );
    assert_eq!(
        run_history.cli_settings.state.folder,
        PathBuf::from("./examples/cli/aa_aa/1L/gammaloop_state")
    );
    assert_eq!(run_history.default_runtime_settings.kinematics.e_cm, 300.0);
    Ok(())
}

#[test]
fn test_epem_a_tth_nlo_example_cli() -> Result<()> {
    let state_path = get_tests_workspace_path().join("epem_a_tth_nlo_example");
    let cli = get_example_cli(
        "epem_a_ttxh/NLO/epem_a_tth_NLO.toml",
        &["generate_diagrams"],
        Some(state_path.clone()),
        None,
        true,
    )?;

    assert_eq!(cli.state.model.name, "sm", "Expected SM model to be loaded");
    assert!(
        !cli.state.process_list.processes.is_empty(),
        "No processes were generated"
    );
    assert_eq!(
        cli.state.process_list.processes[0].definition.folder_name,
        "epem_a_tth"
    );
    assert!(
        state_path.join("justfile").exists(),
        "Expected dot export helper files to be created in {}",
        state_path.display()
    );

    clean_test(&state_path);
    Ok(())
}
