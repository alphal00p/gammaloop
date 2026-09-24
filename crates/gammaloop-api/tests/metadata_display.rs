#![cfg(feature = "cli")]

use std::{collections::BTreeMap, fs, path::PathBuf, process::Command};

use gammaloop_api::state::{CommandHistory, RunHistory};
use gammalooprs::{
    graph::threshold_counterterms::ThresholdCountertermSpec,
    settings::runtime::SamplingSettingsParser,
};
use linnet::parser::DotGraph;
use tempfile::tempdir;
use walkdir::WalkDir;

#[test]
fn generated_metadata_stdout_roundtrips_without_state_writes() {
    let temp = tempdir().unwrap();
    let state = temp.path().join("state");
    // Nextest remaps these paths when running an archived build in the Nix sandbox.
    let binary = std::env::var_os("NEXTEST_BIN_EXE_gammaloop")
        .unwrap_or_else(|| env!("CARGO_BIN_EXE_gammaloop").into());
    let fixture = PathBuf::from(std::env::var_os("CARGO_MANIFEST_DIR").unwrap())
        .join("../../tests/resources/graphs/scalar_box.dot");
    let mut history = RunHistory::default();
    history.cli_settings.global = toml::from_str(
        r#"
[n_cores]
generate = 1
[generation.uv]
subtract_uv = false
generate_integrated = false
[generation.evaluator]
compile = false
iterative_orientation_optimization = false
[generation.threshold_subtraction]
enable_thresholds = true
[generation.tropical_subgraph_table]
disable_tropical_generation = true
"#,
    )
    .unwrap();
    history.default_runtime_settings = toml::from_str(
        r#"
[kinematics]
e_cm = 0.2
[kinematics.externals]
type = "constant"
[kinematics.externals.data]
momenta = [[0.1, 0.0, 0.0, 0.1], [0.1, 0.0, 0.0, -0.1], [0.1, 0.1, 0.0, 0.0], "dependent"]
helicities = [0, 0, 0, 0]
[sampling]
sampling_multichanneling = true
default_channel_selection = ["auto:lmb"]
"#,
    )
    .unwrap();
    history.commands = [
        "import model scalars".to_owned(),
        format!("import graphs '{}' -p metadata -i LO", fixture.display()),
        "generate".to_owned(),
    ]
    .iter()
    .map(|command| CommandHistory::from_raw_string(command).unwrap())
    .collect();
    let card = temp.path().join("generate.toml");
    fs::write(&card, toml::to_string_pretty(&history).unwrap()).unwrap();
    let mut generate = Command::new(&binary);
    generate
        .current_dir(temp.path())
        .env("GL_DISPLAY_FILTER", "off")
        .env("GL_LOGFILE_FILTER", "off")
        .arg("-s")
        .arg(&state)
        .arg(&card)
        .args(["quit", "-o"]);
    let generated = generate.output().unwrap();
    assert!(
        generated.status.success(),
        "generation failed:\nstdout:\n{}\nstderr:\n{}",
        String::from_utf8_lossy(&generated.stdout),
        String::from_utf8_lossy(&generated.stderr)
    );
    let snapshot = || {
        WalkDir::new(&state)
            .into_iter()
            .map(Result::unwrap)
            .filter(|entry| entry.file_type().is_file())
            .map(|entry| (entry.path().to_owned(), fs::read(entry.path()).unwrap()))
            .collect::<BTreeMap<_, _>>()
    };
    let before = snapshot();
    let display = |arguments: &[&str]| {
        let output = Command::new(&binary)
            .current_dir(temp.path())
            .env_remove("GL_NO_HARD_WARNINGS")
            .env("GL_DISPLAY_FILTER", "info")
            .env("GL_LOGFILE_FILTER", "off")
            .env("CLICOLOR_FORCE", "1")
            .env("SYMBOLICA_COLOR", "1")
            .arg("--read-only-state")
            .arg("-s")
            .arg(&state)
            .args(arguments)
            .output()
            .unwrap();
        assert!(
            output.status.success(),
            "display failed:\nstdout:\n{}\nstderr:\n{}",
            String::from_utf8_lossy(&output.stdout),
            String::from_utf8_lossy(&output.stderr)
        );
        let stdout = String::from_utf8(output.stdout).unwrap();
        assert!(!stdout.contains('\u{1b}'), "{stdout}");
        for stream in [&*stdout, &*String::from_utf8_lossy(&output.stderr)] {
            assert!(!stream.contains("compile enabled"), "{stream}");
            assert!(!stream.contains("Threshold esurfaces"), "{stream}");
        }
        stdout
    };
    let sampling = display(&[
        "display",
        "integrands",
        "-p",
        "metadata",
        "-i",
        "LO",
        "--show_sampling",
        "toml",
        "--category",
        "generation",
        "--hide-non-existing-thresholds",
        "--show-threshold-functions",
    ]);
    let sampling_document: BTreeMap<String, SamplingSettingsParser> =
        toml::from_str(&sampling).unwrap();
    let selection = &sampling_document["sampling"].channel_selection["scalar_box"];
    assert!(!selection.is_empty());
    assert!(selection
        .iter()
        .all(|selector| selector.starts_with("lmb(")));
    let exported = temp.path().join("sampling.toml");
    fs::write(&exported, &sampling).unwrap();
    let merged = display(&[
        "run", "-c", &format!(
            "set process -p metadata -i LO file '{}'; display integrands -p metadata -i LO --show_sampling toml; quit",
            exported.display()
        ),
    ]);
    assert_eq!(sampling, merged);
    let thresholds = display(&[
        "display",
        "integrands",
        "-p",
        "metadata",
        "-i",
        "LO",
        "--graph",
        "scalar_box",
        "--show_threshold_subtraction",
        "toml",
    ]);
    assert!(thresholds.starts_with("threshold_counterterms = "));
    let graph: DotGraph =
        DotGraph::from_string(format!("digraph exported {{\n{thresholds}\n}}")).unwrap();
    let _: ThresholdCountertermSpec =
        toml::from_str(&graph.global_data.statements["threshold_counterterms"]).unwrap();
    assert_eq!(
        snapshot(),
        before,
        "metadata display changed the saved state"
    );
}
