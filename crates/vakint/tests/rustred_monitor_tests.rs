#![cfg(all(feature = "feynkit-ingress", unix))]
//! Public monitor calls preserve the ordinary mathematical and error results.

#[path = "support/external_data.rs"]
mod external_data_fixture;
#[allow(dead_code)]
mod five_loop_closure;

use std::{fs, path::Path, sync::atomic::AtomicBool, time::Duration};

use feynkit_kinematics::Kinematics;
use rustred::family::IntegralKey;
use rustred_app::CachePreparationStatus;
use serde_json::Value;
use symbolica::atom::Atom;
use vakint::{
    EvaluationMethod, EvaluationOrder, LoopNormalizationFactor, MasterSymbol,
    OfficialExternalDataOptions, RustRedCachePreparationOptions, RustRedCacheRoot,
    RustRedCacheWorkload, RustRedEvaluationOptions, RustRedMasterMode, RustRedMonitorOptions,
    Vakint, VakintExpression, VakintSettings, diagram_integral::DiagramIntegralOptions,
    vakint_parse,
};

type Result<T> = std::result::Result<T, Box<dyn std::error::Error>>;

fn settings(mode: RustRedMasterMode) -> VakintSettings {
    VakintSettings {
        epsilon_symbol: "ep".into(),
        evaluation_order: EvaluationOrder(vec![EvaluationMethod::RustRed(
            RustRedEvaluationOptions { master_mode: mode },
        )]),
        integral_normalization_factor: LoopNormalizationFactor::FMFTandMATAD,
        form_exe_path: "/monitor-test-must-not-use-form".into(),
        use_dot_product_notation: true,
        ..Default::default()
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
    assert!(latest(&path, "succeeded")?.get("stale").is_none());
    assert_eq!(
        fs::read_dir(directory.path())?
            .map(|entry| entry.map(|entry| entry.path()))
            .collect::<std::io::Result<Vec<_>>>()?,
        vec![path]
    );
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
fn one_loop_cache_monitor_keeps_complete_and_cancelled_report_status() -> Result<()> {
    let directory = tempfile::tempdir()?;
    let engine = Vakint::new()?;
    let path = directory.path().join("latest.json");
    let monitor = snapshot_options(&path);
    let options = RustRedCachePreparationOptions {
        workload: RustRedCacheWorkload::Explicit(vec![RustRedCacheRoot {
            namespace: 0,
            integral: IntegralKey::try_new([2])?,
        }]),
        ..Default::default()
    };
    for cancelled in [false, true] {
        let cancellation = AtomicBool::new(cancelled);
        let reference = engine.prepare_cache(
            1,
            directory.path().join(format!("ordinary-{cancelled}")),
            &options,
            &cancellation,
            |_| Ok(()),
        )?;
        let observed = engine.prepare_cache_with_monitor(
            1,
            directory.path().join(format!("observed-{cancelled}")),
            &options,
            &cancellation,
            |_| Ok(()),
            &monitor,
        )?;
        assert_eq!(observed.status, reference.status);
        assert_eq!(observed.total_roots, reference.total_roots);
        assert_eq!(observed.completed_roots, reference.completed_roots);
        assert_eq!(observed.total_roots, 1);
        let document = latest(&path, "succeeded")?;
        assert_eq!(document["core"], Value::Null);
        if cancelled {
            assert_eq!(observed.status, CachePreparationStatus::Cancelled);
            assert_eq!(observed.completed_roots, 0);
            // The call returned its real partial report. Successful return is
            // not completion of the selected workload; the phase stays explicit.
            assert_eq!(document["caller_phase"], "cache_preparation_cancelled");
            assert_eq!(document["monitor"]["view"]["status"], "Partial");
        } else {
            assert_eq!(observed.status, CachePreparationStatus::Complete);
            assert_eq!(observed.completed_roots, 1);
            assert_eq!(document["caller_phase"], "cache_preparation_complete");
            assert_eq!(document["monitor"]["view"]["status"], "Succeeded");
        }
        assert!(!engine.cache_statistics(1)?.attached);
    }
    Ok(())
}

#[test]
fn five_loop_monitor_preserves_full_results_and_final_substitution_errors() -> Result<()> {
    // One authenticated load serves every 5L assertion in this test process.
    let mut engine = Vakint::new()?;
    let report = engine.load_official_external_data(
        5,
        &OfficialExternalDataOptions {
            version: external_data_fixture::required_version(),
            local_artifacts_path: std::env::var_os("VAKINT_TEST_EXTERNAL_DATA").map(Into::into),
            load_prepared_caches: false,
            ..Default::default()
        },
    )?;
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

    // Inspect an already declared, positive-power scalar master. Its stored
    // graph/routing/powers supply the public ingress; no new reduction search or
    // replacement master fixture is introduced to provoke an error.
    let manifest: Value = serde_json::from_slice(&fs::read(report.path.join("5L/manifest.json"))?)?;
    let metadata = manifest["masters"]
        .as_array()
        .ok_or("5L manifest has no masters")?
        .iter()
        .find(|master| master["name"] == "P640")
        .ok_or("authenticated campaign lacks P640")?;
    let powers: Vec<i64> = serde_json::from_value(metadata["powers"].clone())?;
    assert_eq!(powers, [4, 2, 1, 0, 1, 1, 1, 0, 0, 2, 1, 1, 1, 1, 1]);
    assert_eq!(metadata["status"], "not_computed");
    let symbol: MasterSymbol = metadata["master_id"]
        .as_str()
        .ok_or("P640 has no native identity")?
        .parse()?;
    let master = engine.get_master(&symbol, 16)?;
    assert_eq!(master.loop_count(), 5);
    assert!(!master.is_computed());
    assert!(master.coefficients().is_empty());
    assert!(master.graph().numerator().expression().is_one());
    let dimension = Atom::var(feynkit_graph::symbols::dimension());
    let kinematics = Kinematics::in_dimension(&dimension)?;
    let family = master.graph().integral_family(&kinematics, &[])?;
    let input = Atom::from(VakintExpression::from_diagram(
        master.graph(),
        &family,
        master.graph().numerator().expression(),
        &DiagramIntegralOptions::default(),
    )?);
    for (mode, phase) in [
        (RustRedMasterMode::Symbolic, "symbolic_master_substitution"),
        (
            RustRedMasterMode::Numerical,
            "numerical_master_substitution",
        ),
    ] {
        let settings = self::settings(mode);
        let reference = engine.evaluate(&settings, input.as_view()).unwrap_err();
        assert!(
            reference
                .to_string()
                .contains("declared terminal has no computed coefficients")
        );
        assert!(!reference.to_string().contains("undeclared residual master"));
        let observed = engine
            .evaluate_with_monitor(&settings, input.as_view(), &options)
            .unwrap_err();
        assert_eq!(
            std::mem::discriminant(&observed),
            std::mem::discriminant(&reference)
        );
        assert_eq!(observed.to_string(), reference.to_string());
        let request = latest(&path, "failed")?;
        assert_eq!(request["caller_phase"], phase);
        // Native reduction succeeded; the whole request failed later. Core
        // success must never overwrite this final caller-owned error status.
        assert_eq!(request["core"]["status"], "succeeded");
    }
    Ok(())
}
