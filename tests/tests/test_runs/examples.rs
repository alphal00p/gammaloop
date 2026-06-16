use super::*;
use gammaloop_api::{StateLoadOption, state::RunHistory};
use gammaloop_integration_tests::workspace_root;
use std::fs;

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
    let mut cli = get_example_cli(
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

    let inspections = ["bubble", "bubble_no_integrated_UV"].map(|process| Inspect {
        process: Some(ProcessRef::Unqualified(process.to_string())),
        integrand_name: Some("scalar_bubble_below_thres".to_string()),
        point: vec![0.11, 0.23, 0.37],
        discrete_dim: vec![0, 0],
        ..Default::default()
    });
    let expected = inspections
        .iter()
        .map(|inspect| inspect.run(&mut cli.state))
        .collect::<Result<Vec<_>>>()?;
    let mut reloaded = StateLoadOption::read_only(&state_path).load()?;
    for (inspect, expected) in inspections.iter().zip(expected) {
        assert!(
            expected.1.re.is_finite()
                && expected.1.im.is_finite()
                && expected.1.re.hypot(expected.1.im) > 0.0,
            "state roundtrip requires a finite, nonzero reference"
        );
        assert_eq!(
            inspect.run(&mut reloaded.state)?,
            expected,
            "saved amplitude changed after reload"
        );
    }

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
