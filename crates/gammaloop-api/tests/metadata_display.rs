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

#[test]
fn saved_gl15_amplitude_preserves_model_parameters_in_a_fresh_process() {
    let temp = tempdir().unwrap();
    let state = temp.path().join("state");
    let before = temp.path().join("before.json");
    let after = temp.path().join("after.json");
    let binary = std::env::var_os("NEXTEST_BIN_EXE_gammaloop")
        .unwrap_or_else(|| env!("CARGO_BIN_EXE_gammaloop").into());
    let mut history = RunHistory::default();
    history.cli_settings.global = toml::from_str(
        r#"
[n_cores]
feyngen = 1
generate = 1
compile = 1
integrate = 1
[generation.evaluator]
compile = false
store_atom = false
iterative_orientation_optimization = false
summed = false
summed_function_map = true
[generation.threshold_subtraction]
enable_thresholds = true
check_esurface_at_generation = true
"#,
    )
    .unwrap();
    history.default_runtime_settings = toml::from_str(
        r#"
[general]
evaluator_method = "SingleParametric"
enable_cache = false
generate_events = true
store_additional_weights_in_event = true
[kinematics.externals]
type = "constant"
[kinematics.externals.data]
momenta = [
    [500.0, 0.0, 0.0, 500.0],
    [500.0, 0.0, 0.0, -500.0],
    [438.5555662246945, 155.3322001835378, 348.0160396513587, -177.3773615718412],
    [356.3696374921922, -16.802389008511, -318.7291102436005, 97.48719163688098],
    "dependent",
]
helicities = [1, 1, 0, 0, 0]
[sampling]
graphs = "monte_carlo"
orientations = "monte_carlo"
sampling_multichanneling = true
sampling_channels = "monte_carlo"
[stability]
rotation_axis = []
"#,
    )
    .unwrap();
    history.commands = [
        "import model sm-default".to_owned(),
        r#"generate amp g g > h h h / u d c s b QED==3 [{1}]
            --only-diagrams --numerator-grouping only_detect_zeroes
            --select-graphs GL15 --loop-momentum-bases GL15=8
            --global-prefactor-projector 'gammalooprs::ϵ(0,spenso::mink(4,gammalooprs::hedge(0)))
                * gammalooprs::ϵ(1,spenso::mink(4,gammalooprs::hedge(1)))
                * (1/8)*spenso::g(spenso::coad(8,gammalooprs::hedge(0)),spenso::coad(8,gammalooprs::hedge(1)))'
            -p gg_hhh -i 1L"#
            .to_owned(),
        "generate".to_owned(),
        "set model MT=173.0".to_owned(),
        "set model WT=0.0".to_owned(),
        "set model ymt=173.0".to_owned(),
        format!(
            "inspect -p gg_hhh -i 1L -x 0.23 0.41 0.67 -d 0 0 0 --json-output '{}'",
            before.display()
        ),
    ]
    .iter()
    .map(|command| CommandHistory::from_raw_string(command).unwrap())
    .collect();
    let card = temp.path().join("generate.toml");
    fs::write(&card, toml::to_string_pretty(&history).unwrap()).unwrap();

    let generated = Command::new(&binary)
        .current_dir(temp.path())
        .env("GL_DISPLAY_FILTER", "off")
        .env("GL_LOGFILE_FILTER", "off")
        .env("RAYON_NUM_THREADS", "1")
        .arg("-s")
        .arg(&state)
        .arg(&card)
        .args(["quit", "-o"])
        .output()
        .unwrap();
    assert!(
        generated.status.success(),
        "generation and in-memory inspection failed:\nstdout:\n{}\nstderr:\n{}",
        String::from_utf8_lossy(&generated.stdout),
        String::from_utf8_lossy(&generated.stderr)
    );
    // A new process has a new Symbolica registry. Model values must follow
    // their saved parameter names rather than the new registry's map order.
    let loaded = Command::new(&binary)
        .current_dir(temp.path())
        .env("GL_DISPLAY_FILTER", "off")
        .env("GL_LOGFILE_FILTER", "off")
        .env("RAYON_NUM_THREADS", "1")
        .arg("--read-only-state")
        .arg("-s")
        .arg(&state)
        .args([
            "inspect",
            "-p",
            "gg_hhh",
            "-i",
            "1L",
            "-x",
            "0.23",
            "0.41",
            "0.67",
            "-d",
            "0",
            "0",
            "0",
            "--json-output",
        ])
        .arg(&after)
        .output()
        .unwrap();
    assert!(
        loaded.status.success(),
        "fresh-process inspection failed:\nstdout:\n{}\nstderr:\n{}",
        String::from_utf8_lossy(&loaded.stdout),
        String::from_utf8_lossy(&loaded.stderr)
    );
    let [reference, reloaded] = [before, after].map(|path| {
        let output: serde_json::Value = serde_json::from_slice(&fs::read(path).unwrap()).unwrap();
        let evaluation = &output["evaluation"];
        (
            evaluation["integrand_result"]["re"].as_f64().unwrap(),
            evaluation["integrand_result"]["im"].as_f64().unwrap(),
            evaluation["parameterization_jacobian"].as_f64().unwrap(),
            evaluation["integrator_weight"].as_f64().unwrap(),
        )
    });
    for value in [reference, reloaded] {
        assert!(
            [value.0, value.1, value.2, value.3]
                .iter()
                .all(|component| component.is_finite())
                && value.0.hypot(value.1) > 0.0,
            "saved-state comparison requires finite, nonzero amplitudes: {value:?}"
        );
    }
    assert_eq!(
        reloaded, reference,
        "model parameters changed the GL15 amplitude after a fresh-process reload"
    );
}
