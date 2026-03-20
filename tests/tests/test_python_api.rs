use std::{
    env,
    ffi::OsString,
    path::Path,
    process::{Command, Output},
};

use color_eyre::{Result, eyre::eyre};
use gammaloop_integration_tests::new_test_artifact_dir;
use serde_json::Value as JsonValue;
use serial_test::serial;

const JSON_SENTINEL: &str = "__GL_JSON__";

fn python_interpreter() -> Result<OsString> {
    if let Some(python) = env::var_os("PYTHON") {
        return Ok(python);
    }

    if let Some(virtual_env) = env::var_os("VIRTUAL_ENV") {
        let candidate = Path::new(&virtual_env).join("bin/python");
        if candidate.exists() {
            return Ok(candidate.into_os_string());
        }
    }

    for candidate in ["python", "python3"] {
        if which::which(candidate).is_ok() {
            return Ok(candidate.into());
        }
    }

    Err(eyre!(
        "Could not locate a Python interpreter. Set $PYTHON or activate a virtual environment before running test_python_api"
    ))
}

fn python_command() -> Result<Command> {
    let command = Command::new(python_interpreter()?);
    // The Python API integration tests are meant to exercise the public
    // `import gammaloop` flow from a user-prepared Python environment.
    // They intentionally do not force a local extension artifact into
    // `PYTHONPATH`, because that can pick up a stale build with a mismatched
    // compile-time Symbolica/OEM configuration.
    // command.env("NO_SYMBOLICA_OEM_LICENSE", "1");
    Ok(command)
}

fn extract_json_from_output(output: &Output) -> Result<JsonValue> {
    let stdout = String::from_utf8_lossy(&output.stdout);
    let payload = stdout
        .lines()
        .rev()
        .find_map(|line| line.strip_prefix(JSON_SENTINEL))
        .ok_or_else(|| {
            eyre!(
                "Python API test did not emit a JSON payload marker.\nstdout:\n{stdout}\nstderr:\n{}",
                String::from_utf8_lossy(&output.stderr)
            )
        })?;
    Ok(serde_json::from_str(payload)?)
}

fn run_python_case(state_name: &str, commands: &[String], body: &str) -> Result<JsonValue> {
    let state_dir = new_test_artifact_dir(state_name)?;
    let commands_json = serde_json::to_string(commands)?;
    let script = format!(
        r#"
import json
import os

try:
    import gammaloop
except Exception as exc:
    raise SystemExit(
        "Failed to import gammaloop for GammaLoop Python API tests. "
        "Run `maturin develop` or `just build-api` first in the active Python environment. "
        "If you build the extension yourself for these tests, make sure it was compiled with NO_SYMBOLICA_OEM_LICENSE=1 so it does not try the gammalooprs OEM path from gammaloop-api. "
        f"Original import error: {{exc}}"
    )

STATE_DIR = os.environ["GL_STATE_DIR"]

def build_api():
    return gammaloop.GammaLoopAPI(state_folder=STATE_DIR)

def run_commands(api, commands):
    for command in commands:
        api.run(command)

def summarize_result(result):
    return {{
        "formatted": str(result),
        "event_groups_len": len(result.event_groups),
        "group_sizes": [len(group.events) for group in result.event_groups],
        "first_additional_weights_len": (
            len(result.event_groups[0].events[0].additional_weights)
            if result.event_groups and result.event_groups[0].events
            else 0
        ),
        "observable_keys": sorted(result.observables.keys()),
        "histogram_bin_counts": {{
            name: len(hist.bins) for name, hist in result.observables.items()
        }},
        "generated_event_count": result.generated_event_count,
        "accepted_event_count": result.accepted_event_count,
        "event_processing_time_seconds": result.event_processing_time_seconds,
        "parameterization_jacobian": result.parameterization_jacobian,
        "stability_results": [
            {{
                "precision": entry.precision,
                "estimated_relative_accuracy": entry.estimated_relative_accuracy,
                "accepted_as_stable": entry.accepted_as_stable,
                "total_time_seconds": entry.total_time_seconds,
            }}
            for entry in (result.stability_results or [])
        ],
    }}

def summarize_batch_result(result):
    return {{
        "formatted": str(result),
        "samples": [summarize_result(sample) for sample in result.samples],
        "observable_keys": sorted(result.observables.keys()),
        "histogram_bin_counts": {{
            name: len(hist.bins) for name, hist in result.observables.items()
        }},
        "histogram_sample_counts": {{
            name: hist.sample_count for name, hist in result.observables.items()
        }},
    }}

api = build_api()
run_commands(api, {commands_json})

{body}

print("{JSON_SENTINEL}" + json.dumps(payload, sort_keys=True))
"#
    );

    let output = python_command()?
        .env("GL_STATE_DIR", &state_dir)
        .arg("-c")
        .arg(script)
        .output()?;

    if !output.status.success() {
        return Err(eyre!(
            "Python API test subprocess failed.\nstdout:\n{}\nstderr:\n{}",
            String::from_utf8_lossy(&output.stdout),
            String::from_utf8_lossy(&output.stderr)
        ));
    }

    extract_json_from_output(&output)
}

fn base_setup_commands() -> Vec<String> {
    vec![
        "import model sm-default".to_string(),
        r#"set default-runtime kv kinematics.externals='{"type":"constant","data":{"momenta":[[32.0,0.0,0.0,32.0],[32.0,0.0,0.0,-32.0]],"helicities":[1,1]}}'"#
            .to_string(),
        "set default-runtime kv subtraction.disable_threshold_subtraction=true".to_string(),
        "generate xs e+ e- > d d~ g | e- a d g QED^2==4 [{{2}} QCD=0] --numerator-grouping group_identical_graphs_up_to_sign --clear-existing-processes --only-diagrams".to_string(),
        "generate".to_string(),
    ]
}

fn default_xspace_point() -> &'static str {
    "[0.17, 0.31, 0.53, 0.23, 0.41, 0.67]"
}

fn default_momentum_space_point() -> &'static str {
    "[0.11, -0.07, 0.19, -0.13, 0.05, 0.29]"
}

#[test]
#[serial]
fn python_evaluate_sample_preserves_graph_grouping_and_incoming_pdgs() -> Result<()> {
    let mut commands = base_setup_commands();
    commands.push(
        "set process kv general.generate_events=true general.store_additional_weights_in_event=true"
            .to_string(),
    );

    let payload = run_python_case(
        "python_api_event_grouping",
        &commands,
        &format!(
            r#"
point = {}
result = api.evaluate_sample(point)
payload = {{
    "group_sizes": [len(group.events) for group in result.event_groups],
    "groups": [
        [
            {{
                "graph_id": event.cut_info.graph_id,
                "cut_id": event.cut_info.cut_id,
                "incoming_pdgs": list(event.cut_info.incoming_pdgs),
            }}
            for event in group.events
        ]
        for group in result.event_groups
    ],
    "generated_event_count": result.generated_event_count,
    "accepted_event_count": result.accepted_event_count,
    "formatted": str(result),
}}
"#,
            default_xspace_point()
        ),
    )?;

    assert_eq!(
        payload["group_sizes"].as_array().unwrap(),
        &vec![JsonValue::from(1_u64), JsonValue::from(2_u64)]
    );
    assert_eq!(payload["generated_event_count"].as_u64(), Some(3));
    assert_eq!(payload["accepted_event_count"].as_u64(), Some(3));

    let groups = payload["groups"].as_array().unwrap();
    assert_eq!(groups[0][0]["graph_id"].as_u64(), Some(0));
    assert_eq!(groups[0][0]["cut_id"].as_u64(), Some(0));
    assert_eq!(
        groups[0][0]["incoming_pdgs"]
            .as_array()
            .unwrap()
            .iter()
            .filter_map(|value| value.as_i64())
            .collect::<Vec<_>>(),
        vec![-11, 11]
    );
    assert_eq!(groups[1][0]["graph_id"].as_u64(), Some(1));
    assert_eq!(groups[1][0]["cut_id"].as_u64(), Some(0));
    assert_eq!(groups[1][1]["graph_id"].as_u64(), Some(1));
    assert_eq!(groups[1][1]["cut_id"].as_u64(), Some(1));
    let formatted = payload["formatted"].as_str().unwrap();
    assert!(formatted.contains("Kinematics"));
    assert!(formatted.contains("PDG"));
    assert!(formatted.contains("state"));
    assert!(!formatted.contains("incoming PDGs"));

    Ok(())
}

#[test]
#[serial]
fn python_evaluate_sample_honors_generate_events_and_observable_snapshots() -> Result<()> {
    let mut commands = base_setup_commands();
    commands.push(
        r#"set process string '
[quantities.leading_jet_pt]
type = "jet_pt"
dR = 0.4

[quantities.jet_count]
type = "jet_count"
dR = 0.4

[selectors.leading_jet_pt_cut]
quantity = "leading_jet_pt"
selector = "value_range"
entry_selection = "leading_only"
min = 0.0

[observables.leading_jet_pt_hist]
quantity = "leading_jet_pt"
entry_selection = "leading_only"
x_min = 0.0
x_max = 1000.0
n_bins = 8

[observables.jet_count_hist]
quantity = "jet_count"
entry_selection = "all"
x_min = 0.0
x_max = 6.0
n_bins = 6
'"#
        .to_string(),
    );
    commands.push("set process kv general.generate_events=false".to_string());

    let payload = run_python_case(
        "python_api_evaluate_sample",
        &commands,
        &format!(
            r#"
point = {}
without_events = summarize_result(api.evaluate_sample(point))
api.run("set process kv general.generate_events=true general.store_additional_weights_in_event=true")
with_events = summarize_result(api.evaluate_sample(point))
minimal = summarize_result(api.evaluate_sample(point, minimal_output=True))
payload = {{
    "without_events": without_events,
    "with_events": with_events,
    "minimal": minimal,
}}
"#,
            default_xspace_point()
        ),
    )?;

    let without_events = &payload["without_events"];
    assert_eq!(without_events["event_groups_len"].as_u64(), Some(0));
    assert_eq!(
        without_events["observable_keys"].as_array().map(Vec::len),
        Some(2)
    );
    assert_eq!(
        without_events["observable_keys"][0].as_str(),
        Some("jet_count_hist")
    );
    assert_eq!(
        without_events["observable_keys"][1].as_str(),
        Some("leading_jet_pt_hist")
    );
    assert!(without_events["generated_event_count"].as_u64().unwrap() > 0);
    assert!(without_events["accepted_event_count"].as_u64().unwrap() > 0);
    assert!(
        without_events["accepted_event_count"].as_u64().unwrap()
            <= without_events["generated_event_count"].as_u64().unwrap()
    );
    assert!(
        without_events["event_processing_time_seconds"]
            .as_f64()
            .unwrap()
            >= 0.0
    );
    assert!(without_events["parameterization_jacobian"].is_number());
    assert!(
        without_events["stability_results"]
            .as_array()
            .map(Vec::len)
            .unwrap_or(0)
            > 0
    );

    let with_events = &payload["with_events"];
    assert!(with_events["event_groups_len"].as_u64().unwrap() > 0);
    assert!(
        with_events["first_additional_weights_len"]
            .as_u64()
            .unwrap()
            > 0
    );
    assert!(
        with_events["formatted"]
            .as_str()
            .unwrap()
            .contains("Stability results")
    );

    let minimal = &payload["minimal"];
    assert!(minimal["generated_event_count"].is_null());
    assert!(minimal["accepted_event_count"].is_null());
    assert!(minimal["event_processing_time_seconds"].is_null());
    assert!(minimal["event_groups_len"].as_u64().unwrap() > 0);
    assert_eq!(minimal["observable_keys"].as_array().map(Vec::len), Some(2));
    assert_eq!(
        minimal["observable_keys"][0].as_str(),
        Some("jet_count_hist")
    );
    assert_eq!(
        minimal["observable_keys"][1].as_str(),
        Some("leading_jet_pt_hist")
    );
    assert!(
        !minimal["formatted"]
            .as_str()
            .unwrap()
            .contains("Stability results")
    );

    Ok(())
}

#[test]
#[serial]
fn python_evaluate_samples_batch_and_momentum_space_have_expected_shape() -> Result<()> {
    let mut commands = base_setup_commands();
    commands.push(
        r#"set process string '
[quantities.leading_jet_pt]
type = "jet_pt"
dR = 0.4

[quantities.jet_count]
type = "jet_count"
dR = 0.4

[observables.leading_jet_pt_hist]
quantity = "leading_jet_pt"
entry_selection = "leading_only"
x_min = 0.0
x_max = 1000.0
n_bins = 8

[observables.jet_count_hist]
quantity = "jet_count"
entry_selection = "all"
x_min = 0.0
x_max = 6.0
n_bins = 6
'"#
        .to_string(),
    );
    commands.push("set process kv general.generate_events=true".to_string());

    let payload = run_python_case(
        "python_api_evaluate_samples_batch",
        &commands,
        &format!(
            r#"
try:
    import numpy as np
except Exception as exc:
    raise SystemExit(
        f"Failed to import numpy for GammaLoop Python batch API tests: {{exc}}. "
        "Install numpy in the active Python environment used for the API tests."
    )
points = np.array([{}, {}], dtype=float)
batch = summarize_batch_result(api.evaluate_samples(points))
momentum = summarize_result(api.evaluate_sample({}, momentum_space=True))
payload = {{
    "batch": batch,
    "momentum": momentum,
}}
"#,
            default_xspace_point(),
            default_xspace_point(),
            default_momentum_space_point()
        ),
    )?;

    let batch = &payload["batch"];
    let samples = batch["samples"]
        .as_array()
        .expect("batch samples payload must be a list");
    assert_eq!(samples.len(), 2);
    for sample in samples {
        assert!(sample["event_groups_len"].as_u64().unwrap() > 0);
    }
    assert_eq!(batch["observable_keys"].as_array().map(Vec::len), Some(2));
    assert_eq!(batch["observable_keys"][0].as_str(), Some("jet_count_hist"));
    assert_eq!(
        batch["observable_keys"][1].as_str(),
        Some("leading_jet_pt_hist")
    );
    assert_eq!(
        batch["histogram_bin_counts"]["leading_jet_pt_hist"].as_u64(),
        Some(8)
    );
    assert_eq!(
        batch["histogram_bin_counts"]["jet_count_hist"].as_u64(),
        Some(6)
    );
    assert_eq!(
        batch["histogram_sample_counts"]["leading_jet_pt_hist"].as_u64(),
        Some(2)
    );
    assert_eq!(
        batch["histogram_sample_counts"]["jet_count_hist"].as_u64(),
        Some(2)
    );

    let momentum = &payload["momentum"];
    assert!(momentum["parameterization_jacobian"].is_null());
    assert!(momentum["event_processing_time_seconds"].as_f64().unwrap() >= 0.0);
    let formatted = momentum["formatted"].as_str().unwrap();
    assert!(formatted.contains("parameterization jacobian"));
    assert!(formatted.contains("None"));

    Ok(())
}
