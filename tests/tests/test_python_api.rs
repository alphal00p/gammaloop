use std::{
    env,
    ffi::OsString,
    path::Path,
    process::{Command, Output},
};

use color_eyre::{Result, eyre::eyre};
use gammaloop_integration_tests::{new_test_artifact_dir, workspace_root};
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
    // Symbolica version. Library licensing uses the user's runtime environment;
    // NO_SYMBOLICA_OEM_LICENSE only affects compilation of the CLI binary.
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
        "Set SYMBOLICA_LICENSE to your user license when running these tests; the extension does not activate the CLI's OEM license. "
        f"Original import error: {{exc}}"
    )

STATE_DIR = os.environ["GL_STATE_DIR"]

def build_api():
    return gammaloop.GammaLoopAPI(state_folder=STATE_DIR)

def run_commands(api, commands):
    for command in commands:
        api.run(command)

def summarize_sample_result(result):
    return {{
        "formatted": str(result),
        "event_groups_len": len(result.event_groups),
        "group_sizes": [len(group.events) for group in result.event_groups],
        "first_additional_weights_len": (
            len(result.event_groups[0].events[0].additional_weights)
            if result.event_groups and result.event_groups[0].events
            else 0
        ),
        "generated_event_count": result.generated_event_count,
        "accepted_event_count": result.accepted_event_count,
        "event_processing_time_seconds": result.event_processing_time_seconds,
        "parameterization_jacobian": result.parameterization_jacobian,
        "stability_results": [
            {{
                "precision": entry.precision,
                "estimated_relative_accuracy": entry.estimated_relative_accuracy,
                "status": entry.status,
                "sample_count": entry.sample_count,
                "total_time_seconds": entry.total_time_seconds,
            }}
            for entry in (result.stability_results or [])
        ],
    }}

def summarize_result(result):
    summary = summarize_sample_result(result)
    summary.update({{
        "observable_keys": sorted(result.observables.keys()),
        "histogram_bin_counts": {{
            name: len(hist.bins) for name, hist in result.observables.items()
        }},
    }})
    return summary

def summarize_batch_result(result):
    return {{
        "formatted": str(result),
        "samples": [summarize_sample_result(sample) for sample in result.samples],
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
fn python_settings_wrappers_expose_nested_runtime_and_global_settings() -> Result<()> {
    let mut commands = base_setup_commands();
    commands.push(
        "set global kv global.display_directive=warn global.generation.threshold_subtraction.enable_thresholds=false"
            .to_string(),
    );

    let full_payload = run_python_case(
        "python_api_settings_wrapper",
        &commands,
        r#"
integrand_settings = api.get_integrand_settings()
import tomllib
global_settings = tomllib.loads(api.get_global_settings())
default_runtime_settings = api.get_default_runtime_settings()

payload = {
    "integrand_disable_threshold_subtraction": integrand_settings.subtraction.disable_threshold_subtraction,
    "integrand_external_type": integrand_settings.get("kinematics.externals.type"),
    "integrand_to_dict_disable_threshold_subtraction": (
        integrand_settings.to_dict()["subtraction"]["disable_threshold_subtraction"]
    ),
    "integrand_keys": integrand_settings.keys(),
    "global_display_directive": global_settings["global"]["display_directive"],
    "global_thresholds_enabled": global_settings["global"]["generation"]["threshold_subtraction"]["enable_thresholds"],
    "default_runtime_disable_threshold_subtraction": (
        default_runtime_settings.subtraction.disable_threshold_subtraction
    ),
    "integrand_repr": repr(integrand_settings),
}
"#,
    );

    match full_payload {
        Ok(payload) => {
            assert_eq!(
                payload["integrand_disable_threshold_subtraction"].as_bool(),
                Some(true)
            );
            assert_eq!(
                payload["integrand_external_type"].as_str(),
                Some("constant")
            );
            assert_eq!(
                payload["integrand_to_dict_disable_threshold_subtraction"].as_bool(),
                Some(true)
            );
            assert!(
                payload["integrand_keys"]
                    .as_array()
                    .unwrap()
                    .iter()
                    .any(|value| value.as_str() == Some("subtraction"))
            );
            assert_eq!(payload["global_display_directive"].as_str(), Some("warn"));
            assert_eq!(payload["global_thresholds_enabled"].as_bool(), Some(false));
            assert_eq!(
                payload["default_runtime_disable_threshold_subtraction"].as_bool(),
                Some(true)
            );
            assert!(
                payload["integrand_repr"]
                    .as_str()
                    .unwrap()
                    .contains("SettingsValue(")
            );
        }
        Err(err)
            if err
                .to_string()
                .contains("Cannot start new unlicensed Symbolica instance") =>
        {
            let fallback_commands = vec![
                "set global kv global.display_directive=warn global.generation.threshold_subtraction.enable_thresholds=false"
                    .to_string(),
                "set default-runtime kv subtraction.disable_threshold_subtraction=true"
                    .to_string(),
            ];
            let payload = run_python_case(
                "python_api_settings_wrapper_fallback",
                &fallback_commands,
                r#"
import tomllib
global_settings = tomllib.loads(api.get_global_settings())
default_runtime_settings = api.get_default_runtime_settings()

payload = {
    "global_display_directive": global_settings["global"]["display_directive"],
    "global_thresholds_enabled": global_settings["global"]["generation"]["threshold_subtraction"]["enable_thresholds"],
    "default_runtime_disable_threshold_subtraction": (
        default_runtime_settings.subtraction.disable_threshold_subtraction
    ),
    "default_runtime_keys": default_runtime_settings.keys(),
    "default_runtime_to_dict_disable_threshold_subtraction": (
        default_runtime_settings.to_dict()["subtraction"]["disable_threshold_subtraction"]
    ),
}
"#,
            )?;

            assert_eq!(payload["global_display_directive"].as_str(), Some("warn"));
            assert_eq!(payload["global_thresholds_enabled"].as_bool(), Some(false));
            assert_eq!(
                payload["default_runtime_disable_threshold_subtraction"].as_bool(),
                Some(true)
            );
            assert_eq!(
                payload["default_runtime_to_dict_disable_threshold_subtraction"].as_bool(),
                Some(true)
            );
            assert!(
                payload["default_runtime_keys"]
                    .as_array()
                    .unwrap()
                    .iter()
                    .any(|value| value.as_str() == Some("subtraction"))
            );
        }
        Err(err) => return Err(err),
    }

    Ok(())
}

#[test]
#[serial]
fn python_settings_expose_ecm_relative_stability_tolerances() -> Result<()> {
    let payload = run_python_case(
        "python_api_ecm_relative_stability_tolerances",
        &[r#"set default-runtime string '
[stability]
integrated_energy_dimension = -2
[[stability.levels]]
precision = "Arb"
required_precision_for_re = 1e-12
required_precision_for_im = 1e-12
ecm_relative_tolerance_for_re = 1e-100
escalate_for_large_weight_threshold = -1.0
'"#
        .to_string()],
        r#"
settings = api.get_default_runtime_settings()
level = settings.stability.levels[0]
payload = {
    "dimension": settings.stability.integrated_energy_dimension,
    "serialized_dimension": settings.to_dict()["stability"]["integrated_energy_dimension"],
    "real": level.ecm_relative_tolerance_for_re,
    "imaginary": level.ecm_relative_tolerance_for_im,
    "serialized": settings.to_dict()["stability"]["levels"][0],
}
"#,
    )?;
    assert_eq!(payload["dimension"].as_i64(), Some(-2));
    assert_eq!(payload["serialized_dimension"], payload["dimension"]);
    assert_eq!(payload["real"].as_f64(), Some(1e-100));
    assert_eq!(payload["imaginary"].as_f64(), Some(0.0));
    assert_eq!(
        payload["serialized"]["ecm_relative_tolerance_for_re"],
        payload["real"]
    );
    assert_eq!(
        payload["serialized"]["ecm_relative_tolerance_for_im"],
        payload["imaginary"]
    );
    Ok(())
}

#[test]
#[serial]
fn python_settings_preserve_ose_and_named_channel_weight() -> Result<()> {
    let payload = run_python_case(
        "python_api_ose_channel_weights",
        &[r#"set default-runtime string '
[sampling]
sampling_multichanneling = true
sampling_channel_weight = "ose"
alpha = 2.5
default_channel_selection = ["auto:optimized_lmb"]
'"#
        .to_string()],
        r#"
legacy = api.get_default_runtime_settings().to_dict()["sampling"]
run_commands(api, ["""set default-runtime string '
[sampling]
sampling_multichanneling = true
sampling_channel_weight = "map_density"
alpha = 3.0
default_channel_selection = ["soft"]
[sampling.channel_definitions.G.soft]
around = "lmb(1,2)"
parent_lmb = [1,2]
channel_weight = "ose"
'"""])
mixed = api.get_default_runtime_settings().to_dict()["sampling"]
payload = {"legacy": legacy, "mixed": mixed}
"#,
    )?;
    assert_eq!(payload["legacy"]["sampling_channel_weight"], "ose");
    assert_eq!(payload["legacy"]["alpha"], 2.5);
    assert_eq!(payload["mixed"]["sampling_channel_weight"], "map_density");
    assert_eq!(payload["mixed"]["alpha"], 3.0);
    assert_eq!(
        payload["mixed"]["channel_definitions"]["G"]["soft"]["channel_weight"],
        "ose"
    );
    Ok(())
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
fn python_session_getters_return_live_toml_and_command_blocks() -> Result<()> {
    let commands = vec!["set global kv global.display_directive=warn".to_string()];

    let payload = run_python_case(
        "python_session_getters_return_live_toml_and_command_blocks",
        &commands,
        r#"
import tomllib

api.run("start_commands_block demo")
api.run("set global kv global.logfile_directive=error")
api.run("finish_commands_block")
api.run("run demo -c 'set default-runtime kv general.mu_r=3.0'")

run_history = tomllib.loads(api.get_run_history())
global_settings = tomllib.loads(api.get_global_settings())
command_blocks = api.get_active_command_blocks()

payload = {
    "history_commands": run_history["commands"],
    "global_display_directive": global_settings["global"]["display_directive"],
    "global_logfile_directive": global_settings["global"]["logfile_directive"],
    "command_blocks": command_blocks,
}
"#,
    )?;

    assert_eq!(
        payload["history_commands"],
        serde_json::json!([
            "set global kv global.display_directive=warn",
            "run demo -c 'set default-runtime kv general.mu_r=3.0'",
        ])
    );
    assert_eq!(payload["global_display_directive"].as_str(), Some("warn"));
    assert_eq!(payload["global_logfile_directive"].as_str(), Some("error"));
    assert_eq!(
        payload["command_blocks"],
        serde_json::json!({
            "demo": ["set global kv global.logfile_directive=error"]
        })
    );

    Ok(())
}

#[test]
#[serial]
fn python_evaluate_sample_honors_generate_events_and_observable_snapshots() -> Result<()> {
    let mut commands = base_setup_commands();
    commands.push(
        r#"set process string '
[quantities.leading_jet_pt]
type = "jet"
quantity = "PT"
dR = 0.4

[quantities.jet_count]
type = "jet"
computation = "count"
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
kind = "discrete"
domain = { type = "explicit_range", min = 0, max = 5 }
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
type = "jet"
quantity = "PT"
dR = 0.4

[quantities.jet_count]
type = "jet"
computation = "count"
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
kind = "discrete"
domain = { type = "explicit_range", min = 0, max = 5 }
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
try:
    api.evaluate_sample([0.1, 0.2, 0.3, 0.4], momentum_space=True)
except Exception as exc:
    momentum_shape_error = str(exc)
else:
    raise AssertionError("momentum-space evaluation accepted an incomplete triplet")
payload = {{
    "batch": batch,
    "momentum": momentum,
    "momentum_shape_error": momentum_shape_error,
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
    assert!(payload["momentum_shape_error"]
        .as_str()
        .unwrap()
        .contains(
            "Momentum-space evaluation expects flattened (px, py, pz) triplets, so the coordinate count must be a multiple of 3; got 4."
        ));

    Ok(())
}

#[test]
#[serial]
fn python_residue_map_preserves_native_affine_maps_and_saved_selection() -> Result<()> {
    let graph = workspace_root().join("tests/resources/graphs/scalar_box.dot");
    let commands = vec![
        "import model scalars-default.json".to_string(),
        "set global kv global.display_directive=warn global.generation.three_dimensional_representations='[\"ltd\",\"cff\"]' global.generation.uv.local_uv_cts_from_expanded_4d_integrands=true global.generation.uv.subtract_uv=false global.generation.uv.generate_integrated=false global.generation.threshold_subtraction.enable_thresholds=false global.generation.evaluator.compile=false global.generation.evaluator.iterative_orientation_optimization=false global.generation.tropical_subgraph_table.disable_tropical_generation=true".to_string(),
        r#"set default-runtime kv kinematics.externals='{"type":"constant","data":{"momenta":[[5.0,0.0,0.0,5.0],[5.0,0.0,0.0,-5.0],[5.0,3.0,0.0,4.0],"dependent"],"helicities":[0,0,0,0]}}' subtraction.disable_threshold_subtraction=true"#.to_string(),
        format!("import graphs '{}' -p residue_box -i scalar", graph.display()),
    ];
    let payload = run_python_case(
        "python_residue_map_affine_snapshot",
        &commands,
        r#"
from fractions import Fraction

assert not hasattr(api, "get_orientations")
try:
    api.get_residue_map("scalar_box")
except ValueError:
    pass
else:
    raise AssertionError("an imported graph has no generated residue map")
api.run("generate existing -p residue_box -i scalar")
info = api.get_integrand_info()
assert info.kind == "amplitude"
assert info.generated_representations == ["ltd", "cff"]
graphs = [graph for group in info.graph_groups for graph in group.graphs]
assert len(graphs) == 1 and graphs[0].name == "scalar_box"
native_counts = dict(graphs[0].native_residue_counts)

def contents(snapshot):
    # Keep exact Python objects in this comparison; JSON is only the final
    # subprocess summary, never the residue-map transport or persistence oracle.
    entries = {}
    for key, residues in snapshot.entries.items():
        assert type(key).__name__ == "ResidueMapKey"
        assert isinstance(key.directions, tuple)
        assert all(type(direction) is int and direction in (-1, 0, 1)
                   for direction in key.directions)
        assert isinstance(key.loop_energy_map, tuple)
        assert isinstance(key.edge_energy_map, tuple)
        assert key.canonical_string
        for energy in key.loop_energy_map + key.edge_energy_map:
            assert type(energy).__name__ == "LinearEnergyExpression"
            assert isinstance(energy.internal_terms, tuple)
            assert isinstance(energy.external_terms, tuple)
            assert isinstance(energy.uniform_scale_coeff, Fraction)
            assert isinstance(energy.constant, Fraction)
            for edge, coefficient in energy.internal_terms + energy.external_terms:
                assert type(edge) is int
                assert isinstance(coefficient, Fraction)
            assert energy.canonical_string
        assert isinstance(residues, list) and residues
        records = []
        for residue in residues:
            assert type(residue).__name__ == "Residue"
            assert type(residue.native_id) is int
            assert residue.label is None or isinstance(residue.label, str)
            assert residue.numerator_map_index is None or type(residue.numerator_map_index) is int
            variants = []
            for variant in residue.variants:
                assert type(variant).__name__ == "ResidueVariant"
                assert isinstance(variant.prefactor, Fraction)
                tree = variant.denominator
                assert tree["root"] == 0
                nodes = {node["node_id"]: node for node in tree["nodes"]}
                assert len(nodes) == len(tree["nodes"]) and nodes
                assert nodes[0]["parent"] is None
                for node in nodes.values():
                    surface = node["surface"]
                    assert isinstance(surface, tuple)
                    assert surface[0] in ("esurface", "hsurface", "linear", "unit", "infinite")
                    if surface[0] in ("unit", "infinite"):
                        assert surface[1] is None
                    else:
                        assert surface in snapshot.surfaces
                    for child in node["children"]:
                        assert nodes[child]["parent"] == node["node_id"]
                variants.append((variant.origin, variant.prefactor, variant.half_edges,
                                 variant.denominator_edges, variant.denominator_surface_signs,
                                 variant.denominator_edge_support_signs, variant.uniform_scale_power,
                                 variant.numerator_surfaces, tree))
            assert variants
            records.append((residue.native_id, residue.label, residue.numerator_map_index, variants))
        entries[key] = records
    assert entries
    surfaces = snapshot.surfaces
    for reference, surface in surfaces.items():
        if reference[0] in ("esurface", "hsurface"):
            assert surface["kind"] == reference[0]
            assert isinstance(surface["external_shift"], tuple)
            assert all(type(edge) is int and isinstance(coefficient, Fraction)
                       for edge, coefficient in surface["external_shift"])
            energy_fields = ("energies",) if reference[0] == "esurface" else ("positive_energies", "negative_energies")
            assert all(type(edge) is int for field in energy_fields for edge in surface[field])
            assert "vertex_set" in surface
        elif reference[0] == "linear":
            assert surface["kind"] in ("esurface", "hsurface")
            assert surface["origin"] in ("physical", "helper")
            assert type(surface["numerator_only"]) is bool
            assert type(surface["expression"]).__name__ == "LinearEnergyExpression"
    assert snapshot.energy_factor_ownership in ("global_source_product", "variant_local")
    for component in snapshot.energy_factor_components:
        assert component["ownership"] in ("global_source_product", "variant_local")
        assert all(type(edge) is int for edge in component["internal_edge_ids"])
        assert component["denominator_only_global_prefactor_sign"] in (-1, 1)
        assert component["core_global_prefactor_sign"] in (-1, 1)
    assert snapshot.denominator_only_global_prefactor_sign in (-1, 1)
    assert snapshot.core_global_prefactor_sign in (-1, 1)
    assert snapshot.residual_denominators == []
    return (entries, surfaces, snapshot.residual_denominators,
            snapshot.energy_factor_ownership, snapshot.energy_factor_components,
            snapshot.denominator_only_global_prefactor_sign, snapshot.core_global_prefactor_sign)

snapshots = {
    mode: api.get_residue_map("scalar_box", process_id=info.process_id,
                              integrand_name="scalar", three_dimensional_representation=mode)
    for mode in ("ltd", "cff")
}
before = {mode: contents(snapshot) for mode, snapshot in snapshots.items()}
for mode, snapshot in snapshots.items():
    assert snapshot.graph_name == "scalar_box" and snapshot.representation == mode
    assert sum(len(rows) for rows in snapshot.entries.values()) == native_counts[mode]
    assert isinstance(snapshot.entries, dict) and isinstance(snapshot.surfaces, dict)
assert native_counts["ltd"] == 4
assert contents(api.get_residue_map("scalar_box")) == before["ltd"]
assert any(energy.internal_terms and energy.external_terms
           for key in snapshots["ltd"].entries
           for energy in key.edge_energy_map)

key = next(iter(snapshots["ltd"].entries))
original_hash = hash(key)
for obj, field, value in ((key, "directions", ()),
                          (key.edge_energy_map[0], "constant", Fraction(99))):
    try:
        setattr(obj, field, value)
    except AttributeError:
        pass
    else:
        raise AssertionError(f"{type(obj).__name__}.{field} must be immutable")
assert hash(key) == original_hash
assert {key: "native"}[next(k for k in api.get_residue_map("scalar_box").entries if k == key)] == "native"

for kwargs in ({"graph_name": "missing_graph"},
               {"graph_name": "scalar_box", "three_dimensional_representation": "unknown"}):
    try:
        api.get_residue_map(**kwargs)
    except ValueError as exc:
        assert str(exc)
    else:
        raise AssertionError(f"invalid residue-map selection was accepted: {kwargs}")

# Future generation settings must not select a representation in a saved payload.
api.run("set global kv global.generation.three_dimensional_representations='[\"cff\"]' global.generation.explicit_orientation_sum_only=true")
assert contents(api.get_residue_map("scalar_box")) == before["ltd"]
api.run("save state -o")
loaded = gammaloop.GammaLoopAPI(state_folder=STATE_DIR, read_only_state=True)
assert loaded.get_integrand_info().generated_representations == ["ltd", "cff"]
for mode in ("ltd", "cff"):
    restored = loaded.get_residue_map("scalar_box", three_dimensional_representation=mode)
    assert contents(restored) == before[mode]
    assert set(restored.entries) == set(snapshots[mode].entries)
assert contents(loaded.get_residue_map("scalar_box")) == before["ltd"]

# Python containers are detached too: mutating them cannot erase native rows.
detached = loaded.get_residue_map("scalar_box")
detached_key = next(iter(detached.entries))
detached.entries[detached_key].clear()
detached.entries.clear()
assert contents(loaded.get_residue_map("scalar_box")) == before["ltd"]
del loaded

# Rebuild only CFF to exercise unavailable-mode errors and snapshot ownership.
api.run("generate existing -p residue_box -i scalar")
assert api.get_residue_map("scalar_box").representation == "cff"
try:
    api.get_residue_map("scalar_box", three_dimensional_representation="ltd")
except ValueError as exc:
    assert "ltd" in str(exc).lower()
else:
    raise AssertionError("an unavailable representation was silently substituted")
assert contents(snapshots["ltd"]) == before["ltd"]
assert contents(snapshots["cff"]) == before["cff"]
payload = {"native_counts": native_counts, "saved_order": info.generated_representations}
"#,
    )?;
    assert_eq!(payload["native_counts"]["ltd"], 4);
    assert!(payload["native_counts"]["cff"].as_u64().unwrap() > 0);
    assert_eq!(payload["saved_order"], serde_json::json!(["ltd", "cff"]));
    Ok(())
}

#[test]
#[serial]
fn python_residue_map_keeps_raised_cross_section_denominator_structure() -> Result<()> {
    let root = workspace_root();
    let commands = vec![
        format!(
            "import model '{}'",
            root.join("assets/models/json/scalars/scalars_2p_3p.json")
                .display()
        ),
        "set global kv global.display_directive=warn global.generation.three_dimensional_representations='[\"cff\",\"ltd\"]' global.generation.uv.local_uv_cts_from_expanded_4d_integrands=true global.generation.uv.subtract_uv=false global.generation.uv.generate_integrated=false global.generation.threshold_subtraction.enable_thresholds=false global.generation.evaluator.compile=false global.generation.evaluator.iterative_orientation_optimization=false global.generation.tropical_subgraph_table.disable_tropical_generation=true".to_string(),
        r#"set default-runtime kv kinematics.externals='{"type":"constant","data":{"momenta":[[4.0,0.0,0.0,0.0]],"helicities":[0]}}' subtraction.disable_threshold_subtraction=true"#.to_string(),
        format!(
            "import graphs '{}' -p residue_dotted -i scalar",
            root.join("tests/resources/graphs/dotted_bubble.dot").display()
        ),
        "generate existing -p residue_dotted -i scalar".to_string(),
    ];
    let payload = run_python_case(
        "python_residue_map_raised_cross_section",
        &commands,
        r#"
from fractions import Fraction

info = api.get_integrand_info()
assert info.kind == "cross section"
assert info.generated_representations == ["cff", "ltd"]
assert api.get_residue_map("dotted_bubble").representation == "cff"
graphs = [graph for group in info.graph_groups for graph in group.graphs]
assert len(graphs) == 1 and graphs[0].name == "dotted_bubble"
counts = dict(graphs[0].native_residue_counts)
maps = {mode: api.get_residue_map("dotted_bubble", three_dimensional_representation=mode)
        for mode in ("cff", "ltd")}
max_poles = {}
repeated_energy_factors = False
for mode, snapshot in maps.items():
    assert snapshot.graph_name == "dotted_bubble" and snapshot.representation == mode
    assert isinstance(snapshot.entries, dict) and snapshot.entries
    assert sum(len(rows) for rows in snapshot.entries.values()) == counts[mode]
    native_ids = [row.native_id for rows in snapshot.entries.values() for row in rows]
    assert len(set(native_ids)) == len(native_ids)
    maximum = 0
    for key, residues in snapshot.entries.items():
        assert len(key.directions) == len(key.edge_energy_map)
        if mode == "ltd":
            assert key.loop_energy_map
        for residue in residues:
            for variant in residue.variants:
                assert isinstance(variant.prefactor, Fraction)
                assert variant.denominator_edges
                assert all(type(edge) is int for edge in variant.denominator_edges)
                assert all(type(sign) is int and sign in (-1, 1)
                           for sign in variant.denominator_surface_signs.values())
                repeated_energy_factors |= len(set(variant.half_edges)) < len(variant.half_edges)
                tree = variant.denominator
                nodes = {node["node_id"]: node for node in tree["nodes"]}
                pending = [(tree["root"], 0)]
                while pending:
                    node_id, poles = pending.pop()
                    node = nodes[node_id]
                    surface = node["surface"]
                    if surface[0] not in ("unit", "infinite"):
                        assert surface in snapshot.surfaces
                        poles += 1
                    if not node["children"]:
                        maximum = max(maximum, poles)
                    pending.extend((child, poles) for child in node["children"])
    max_poles[mode] = maximum

# The serial pair in this one-loop graph produces a raised pole. A flattened
# set of surfaces or half edges would erase its multiplicity.
assert max_poles["cff"] >= 2
assert repeated_energy_factors
assert any(energy.internal_terms for key in maps["ltd"].entries for energy in key.loop_energy_map)
assert any(energy.external_terms for key in maps["ltd"].entries for energy in key.edge_energy_map)
payload = {"native_counts": counts, "max_denominator_poles": max_poles,
           "repeated_energy_factors": repeated_energy_factors}
"#,
    )?;
    assert!(payload["native_counts"]["cff"].as_u64().unwrap() > 0);
    assert!(payload["native_counts"]["ltd"].as_u64().unwrap() > 0);
    assert!(payload["max_denominator_poles"]["cff"].as_u64().unwrap() >= 2);
    assert_eq!(payload["repeated_energy_factors"].as_bool(), Some(true));
    Ok(())
}
