#![cfg(all(feature = "feynkit-ingress", unix))]
//! Public monitor calls preserve the ordinary mathematical and error results.

#[allow(dead_code)]
mod five_loop_closure;
#[path = "support/rustred.rs"]
#[allow(dead_code)]
mod rustred_support;

use std::{fs, path::Path, time::Duration};

use feynkit_kinematics::Kinematics;
use serde_json::Value;
use symbolica::atom::Atom;
use vakint::{
    LoopNormalizationFactor, MasterSymbol, RustRedMasterMode, RustRedMonitorOptions,
    RustRedProgressSink, Vakint, VakintExpression, VakintSettings,
    diagram_integral::DiagramIntegralOptions, vakint_parse,
};

type Result<T> = std::result::Result<T, Box<dyn std::error::Error>>;

fn settings(mode: RustRedMasterMode) -> VakintSettings {
    VakintSettings {
        epsilon_symbol: "ep".into(),
        integral_normalization_factor: LoopNormalizationFactor::FMFTandMATAD,
        use_dot_product_notation: true,
        ..rustred_support::rustred_only(mode, "/monitor-test-must-not-use-form")
    }
}

fn snapshot_options(path: &Path) -> RustRedMonitorOptions {
    RustRedMonitorOptions {
        // A destination alone opts in. Terminal output remains disabled.
        progress_snapshot: Some(path.to_owned()),
        interval: Duration::ZERO,
        no_color: true,
        ..Default::default()
    }
}

fn latest(path: &Path, status: &str) -> Result<Value> {
    // Parsing the whole file also rejects a JSONL stream or a partial append.
    let value: Value = serde_json::from_slice(&fs::read(path)?)?;
    assert_eq!(value["schema"], "rustred.saved-application-progress.v1");
    assert_eq!(value["event"], "application_progress");
    assert_eq!(value["observation_scope"], "caller");
    assert_eq!(value["activity"], "caller");
    assert_eq!(value["request_finished"], true);
    assert_eq!(value["status"], status);
    assert_eq!(value.get("closure_claim"), Some(&Value::Bool(false)));
    assert_eq!(
        value.get("global_remaining_work_known"),
        Some(&Value::Bool(false))
    );
    assert_eq!(value.get("eta_seconds"), Some(&Value::Null));
    assert_eq!(value.get("percent_complete"), Some(&Value::Null));
    assert_eq!(
        value.get("arithmetic_internals_known"),
        Some(&Value::Bool(false))
    );
    Ok(value)
}

#[test]
fn ordinary_one_loop_monitor_defaults_are_off_and_path_alone_is_observational() -> Result<()> {
    let directory = tempfile::tempdir()?;
    let engine = Vakint::new()?;
    let settings = settings(RustRedMasterMode::Retain);
    let input = vakint_parse!("topo(I1L(1,2))-(1-ep)*topo(I1L(1,1))")?;
    let reference = engine.evaluate(&settings, input.as_view())?;
    assert!(reference.is_zero(), "unit-mass tadpole recurrence changed");

    let default = RustRedMonitorOptions::default();
    assert!(!default.enabled());
    assert_eq!(
        engine.evaluate_with_monitor(&settings, input.as_view(), &default)?,
        reference
    );
    assert_eq!(fs::read_dir(directory.path())?.count(), 0);

    let path = directory.path().join("latest.json");
    let options = snapshot_options(&path);
    assert!(options.enabled());
    assert_eq!(
        engine.evaluate_with_monitor(&settings, input.as_view(), &options)?,
        reference
    );
    let first = latest(&path, "succeeded")?;
    // The direct Route/Apply observer does not claim to measure the 1L backend.
    assert_eq!(first["core"], Value::Null);
    fs::write(&path, b"{\"stale\":true}\n")?;
    assert_eq!(
        engine.evaluate_with_monitor(&settings, input.as_view(), &options)?,
        reference
    );
    let latest = latest(&path, "succeeded")?;
    assert!(latest.get("stale").is_none());
    let size = &latest["monitor"]["caller_expression"];
    assert_eq!(size["phase"], "evaluation_output");
    assert_eq!(size["encoded_bytes"], reference.as_view().get_byte_size());
    assert_eq!(size["stored_summands"], reference.nterms());
    assert!(latest["monitor"]["stored_symbolic_state"].is_null());
    assert_eq!(
        fs::read_dir(directory.path())?
            .map(|entry| entry.map(|entry| entry.path()))
            .collect::<std::io::Result<Vec<_>>>()?,
        vec![path]
    );
    Ok(())
}

#[test]
fn progress_sink_alone_enables_monitor_and_ends_with_final_model() -> Result<()> {
    let directory = tempfile::tempdir()?;
    let engine = Vakint::new()?;
    let settings = settings(RustRedMasterMode::Retain);
    let input = vakint_parse!("topo(I1L(1,2))-(1-ep)*topo(I1L(1,1))")?;
    let reference = engine.evaluate(&settings, input.as_view())?;
    let sink = RustRedProgressSink::default();
    let options = RustRedMonitorOptions {
        // No terminal and no file: only the in-memory model.
        progress_sink: Some(sink.clone()),
        interval: Duration::ZERO,
        ..Default::default()
    };
    assert!(options.enabled());
    options.validate()?;
    let cwd = std::env::current_dir()?;
    let before = fs::read_dir(&cwd)?.count();
    assert_eq!(
        engine.evaluate_with_monitor(&settings, input.as_view(), &options)?,
        reference
    );
    let (sequence, model) = sink.latest_since(0).expect("at least one model");
    assert!(sequence >= 1);
    // The final drain stores the finished request's model.
    assert_eq!(model["status"], "Succeeded");
    assert_eq!(model["api_status"], "Succeeded");
    assert!(model["elapsed_seconds"].as_f64().is_some_and(|s| s >= 0.0));
    assert!(sink.latest_since(sequence).is_none());
    // A second call on the same sink continues its sequence.
    assert_eq!(
        engine.evaluate_integral_with_monitor(&settings, input.as_view(), &options)?,
        engine.evaluate_integral(&settings, input.as_view())?
    );
    let (next, model) = sink.latest_since(sequence).expect("a new model");
    assert!(next > sequence);
    assert_eq!(model["status"], "Succeeded");
    // Nothing was written: neither in the test directory nor the cwd.
    assert_eq!(fs::read_dir(directory.path())?.count(), 0);
    assert_eq!(fs::read_dir(&cwd)?.count(), before);
    Ok(())
}

#[test]
fn snapshot_publication_fault_does_not_replace_the_native_result() -> Result<()> {
    let directory = tempfile::tempdir()?;
    let occupied = directory.path().join("occupied");
    fs::create_dir(&occupied)?;
    let engine = Vakint::new()?;
    let settings = settings(RustRedMasterMode::Retain);
    let input = vakint_parse!("topo(I1L(1,2))")?;
    let reference = engine.evaluate(&settings, input.as_view())?;
    assert!(!reference.is_zero());
    // Configuration is valid, but publishing a regular file over this directory
    // fails in the real sink. That presentation fault must not become a CAS error.
    let options = snapshot_options(&occupied);
    options.validate()?;
    assert_eq!(
        engine.evaluate_with_monitor(&settings, input.as_view(), &options)?,
        reference
    );
    assert!(occupied.is_dir());
    assert_eq!(
        fs::read_dir(directory.path())?
            .map(|entry| entry.map(|entry| entry.path()))
            .collect::<std::io::Result<Vec<_>>>()?,
        vec![occupied]
    );
    Ok(())
}

#[test]
fn five_loop_monitor_preserves_full_results_and_final_substitution_errors() -> Result<()> {
    // One authenticated load serves every 5L assertion in this test process.
    let mut engine = Vakint::new()?;
    let report = engine.load_official_external_data(5, &rustred_support::pinned_release())?;
    assert!(report.loops.iter().any(|pack| pack.loop_count == 5));
    let settings = settings(RustRedMasterMode::Retain);
    let fixture = five_loop_closure::factorized::fixtures(&engine)?
        .into_iter()
        .find(|fixture| fixture.label == "single_positive_line")
        .ok_or("missing existing positive-line fixture")?;
    let directory = tempfile::tempdir()?;
    let path = directory.path().join("latest.json");
    let options = snapshot_options(&path);

    let ordinary = engine.evaluate(&settings, fixture.cancelled.as_view())?;
    assert!(!ordinary.is_zero());
    assert_eq!(
        engine.evaluate_with_monitor(&settings, fixture.cancelled.as_view(), &options)?,
        ordinary
    );
    let request = latest(&path, "succeeded")?;
    assert_eq!(request["core"]["status"], "succeeded");
    // The ordered scheduler's final native state: every phase completed, an
    // empty frontier and the collected terminal keys.
    let ordered = &request["monitor"]["ordered_reduction"];
    assert_eq!(ordered["completed_phases"], ordered["total_phases"]);
    assert_eq!(ordered["frontier_keys"], 0);
    assert!(
        ordered["terminal_keys"]
            .as_u64()
            .is_some_and(|keys| keys > 0)
    );
    assert_eq!(
        request["monitor"]["caller_expression"]["encoded_bytes"],
        ordinary.as_view().get_byte_size()
    );

    let reference = engine.evaluate_rustred_symbolic(&settings, fixture.cancelled.as_view())?;
    let observed = engine.evaluate_rustred_symbolic_with_monitor(
        &settings,
        fixture.cancelled.as_view(),
        &options,
    )?;
    assert_eq!(observed.value, reference.value);
    assert_eq!(observed.nonzero_conditions, reference.nonzero_conditions);
    assert_eq!(observed.families, reference.families);
    assert_eq!(observed.value, ordinary);
    assert!(!reference.nonzero_conditions.is_empty());
    assert!(!reference.families.is_empty());
    assert_eq!(latest(&path, "succeeded")?["core"]["status"], "succeeded");

    // Inspect an already declared, positive-power scalar master: the last one
    // the authenticated 5L manifest declares (P131 of the v0.2.1 pack, P640 of
    // the 640-terminal v0.2.0 campaign). Its stored graph/routing/powers supply
    // the public ingress; no new reduction search or replacement master
    // fixture is introduced to provoke an error. A terminal without
    // coefficients must fail its final substitution as a declared terminal; a
    // computed one must fail closed when the request needs one Laurent order
    // beyond its last published coefficient.
    let manifest: Value = serde_json::from_slice(&fs::read(report.path.join("5L/manifest.json"))?)?;
    let metadata = manifest["masters"]
        .as_array()
        .ok_or("5L manifest has no masters")?
        .iter()
        .rfind(|master| {
            master["powers"].as_array().is_some_and(|powers| {
                powers
                    .iter()
                    .all(|power| power.as_i64().is_some_and(|power| power >= 0))
            })
        })
        .ok_or("authenticated campaign declares no positive-power scalar master")?;
    let name = metadata["name"].as_str().ok_or("5L master has no name")?;
    let identity = metadata["master_id"]
        .as_str()
        .ok_or("5L master has no native identity")?;
    let computed = match metadata["status"].as_str() {
        Some("not_computed") => false,
        Some("accepted") => true,
        status => return Err(format!("{name}: unexpected master status {status:?}").into()),
    };
    eprintln!("selected 5L master {name} ({identity}), computed: {computed}");
    let symbol: MasterSymbol = identity.parse()?;
    let master = engine.get_master(&symbol, 16)?;
    assert_eq!(master.loop_count(), 5);
    assert_eq!(master.is_computed(), computed, "{name}");
    assert_eq!(master.coefficients().is_empty(), !computed, "{name}");
    assert!(master.graph().numerator().expression().is_one());
    // Depth whose last requested order is one beyond the last coefficient:
    // the target order is number_of_terms_in_epsilon_expansion - loops - 1.
    let depth = master
        .coefficients()
        .keys()
        .next_back()
        .map(|last| last + 1 + 5 + 1);
    let dimension = Atom::var(feynkit_graph::symbols::dimension());
    let kinematics = Kinematics::in_dimension(&dimension)?;
    let family = master.graph().integral_family(&kinematics, &[])?;
    let input = Atom::from(VakintExpression::from_diagram(
        master.graph(),
        &family,
        master.graph().numerator().expression(),
        &DiagramIntegralOptions::default(),
    )?);
    // A request made of package masters only is its own reduction: it is
    // built without the native request (master transport), so that snapshot
    // must claim no native core. To keep a request whose native reduction
    // succeeds before its final substitution fails, the same master is also
    // requested next to the positive-line fixture above. That request as a
    // whole runs a native reduction: the fixture's terms still need the
    // native walk (their ordered phases were checked above), while the
    // master itself may still be reduced directly when it falls in a
    // (mass, spectator, exponent) channel of its own.
    let beside_native = input.clone() + &fixture.cancelled;
    for (mode, phase) in [
        (RustRedMasterMode::Symbolic, "symbolic_master_substitution"),
        (
            RustRedMasterMode::Numerical,
            "numerical_master_substitution",
        ),
    ] {
        let mut settings = self::settings(mode);
        if let Some(depth) = depth {
            settings.number_of_terms_in_epsilon_expansion = depth;
        }
        for (request_input, native) in [(&input, false), (&beside_native, true)] {
            let reference = engine
                .evaluate(&settings, request_input.as_view())
                .unwrap_err();
            let text = reference.to_string();
            if computed {
                assert!(
                    matches!(reference, vakint::VakintError::MissingMasterOrders(_)),
                    "{name}: {text}"
                );
                assert!(
                    text.contains("beyond the available expansion depth"),
                    "{name}: {text}"
                );
                assert!(text.contains(identity), "{name}: {text}");
            } else {
                assert!(
                    text.contains("declared terminal has no computed coefficients"),
                    "{name}: {text}"
                );
            }
            assert!(
                !text.contains("undeclared residual master"),
                "{name}: {text}"
            );
            let observed = engine
                .evaluate_with_monitor(&settings, request_input.as_view(), &options)
                .unwrap_err();
            assert_eq!(
                std::mem::discriminant(&observed),
                std::mem::discriminant(&reference)
            );
            assert_eq!(observed.to_string(), reference.to_string());
            let request = latest(&path, "failed")?;
            assert_eq!(request["caller_phase"], phase);
            if native {
                // The request's native reduction (at least the fixture's)
                // succeeded; the whole request failed later. Core success
                // must never overwrite this final caller-owned error status.
                assert_eq!(request["core"]["status"], "succeeded", "{name}");
            } else {
                // No native request ran: the failed status is the caller's own.
                // This deliberately pins the direct path for a request made of
                // the package's own masters only (built without the native
                // request, so no core snapshot is recorded). If that path ever
                // gains a native request, this assertion fails by design: that
                // is a behaviour change of master transport, not a monitor
                // regression.
                assert_eq!(request["core"], Value::Null, "{name}");
            }
        }
    }
    Ok(())
}
