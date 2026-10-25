//! RR-11 vakint acceptance V1-V4 for a staged (scratch) five-loop pack.
//!
//! Driven by environment variables (all paths are scratch copies; nothing is
//! written into any data checkout):
//! - `VAKINT_LAYER_ROOT`: staged external-data root (`manifest.json`, 1L-5L);
//! - `VAKINT_LAYER_BASE_ROOT`: optional unmodified root for V4 (old vs new);
//! - `VAKINT_LAYER_EXPECTED_MASTERS`: N' (V2);
//! - `VAKINT_LAYER_SAMPLES`: raw terminal witnesses reduced in V3/V4 (default 12);
//! - `VAKINT_LAYER_SYMBOLIC`: vakint-level symbolic evaluations compared (default 2);
//! - `VAKINT_LAYER_REPORT`: JSON receipt path.
use super::*;
use rustred_app::SavedApplicationIntegral;
use serde_json::{Value, json};
use std::time::Instant;

fn memory() -> Value {
    let status = std::fs::read_to_string("/proc/self/status").unwrap_or_default();
    let field = |name: &str| {
        status
            .lines()
            .find(|line| line.starts_with(name))
            .and_then(|line| line.split_whitespace().nth(1))
            .and_then(|kb| kb.parse::<u64>().ok())
            .map(|kb| kb * 1024)
    };
    json!({"vm_hwm_bytes": field("VmHWM:"), "vm_rss_bytes": field("VmRSS:")})
}

fn settings() -> VakintSettings {
    VakintSettings {
        form_exe_path: "/five-loop-layer-acceptance-must-not-use-form".into(),
        evaluation_order: crate::EvaluationOrder(vec![crate::EvaluationMethod::RustRed(
            crate::RustRedEvaluationOptions {
                master_mode: crate::RustRedMasterMode::Retain,
            },
        )]),
        ..Default::default()
    }
}

struct Loaded {
    vakint: Vakint,
    settings: VakintSettings,
    application: Arc<SavedApplication>,
    data: Arc<super::super::data::RustRedData>,
    receipt: Value,
}

fn load(root: &Path) -> Loaded {
    let started = Instant::now();
    let mut vakint = Vakint::new().unwrap();
    vakint.load_external_data(root, 5).unwrap();
    let data_seconds = started.elapsed().as_secs_f64();
    let settings = settings();
    let native_started = Instant::now();
    let application = vakint.saved_application(&settings).unwrap();
    let native_seconds = native_started.elapsed().as_secs_f64();
    let data = vakint.rustred_data().unwrap();
    let receipt = json!({
        "root": root,
        "load_external_data_seconds": data_seconds,
        "native_application_seconds": native_seconds,
        "memory_after_load": memory(),
        "metadata_remaining_terminals":
            application.prepared.lock().unwrap().metadata()["remaining_terminals"],
        "refinement_layer_blake3":
            application.prepared.lock().unwrap().metadata()["refinement_layer"]["layer_blake3"],
        "products": format!("{:?}", application.products),
        "product_charts": format!("{:?}", application.product_charts),
    });
    Loaded {
        vakint,
        settings,
        application,
        data,
        receipt,
    }
}

type Terms = BTreeMap<SavedApplicationIntegral, Coefficient>;

fn reduce(application: &SavedApplication, key: &IntegralKey) -> Result<(Terms, usize), String> {
    let mut prepared = application.prepared.lock().unwrap();
    let one = prepared.family().coefficient_context().one();
    let reduced = prepared
        .reduce_linear_combination([(key.clone(), one)])
        .map_err(|e| e.to_string())?;
    Ok((reduced.terms().clone(), reduced.unresolved().len()))
}

/// Admitted raw terminal witnesses, stratified over the sorted raw keys,
/// plus the recorded LRU and banana controls when admitted.
fn witnesses(application: &SavedApplication, count: usize) -> Vec<IntegralKey> {
    let prepared = application.prepared.lock().unwrap();
    let mut keys: Vec<IntegralKey> = [
        [0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 5],
        [0, 0, 0, 0, 1, 1, 0, 0, 1, 0, 0, 1, 0, 1, 4],
        [0, 0, 0, 0, 4, 1, 0, 0, 1, 0, 0, 1, 0, 1, 1],
    ]
    .into_iter()
    .map(|powers| IntegralKey::try_new(powers).unwrap())
    .filter(|key| prepared.routed_entry_is_admitted(key).unwrap_or(false))
    .collect();
    let raw: Vec<_> = prepared
        .terminal_sources()
        .map(|key| key.integral().clone())
        .collect();
    let mut index = 0;
    let stride = (raw.len() / count.max(1)).max(1);
    while keys.len() < count + 3 && index < raw.len() {
        let key = &raw[index];
        if !keys.contains(key) && prepared.routed_entry_is_admitted(key).unwrap_or(false) {
            keys.push(key.clone());
        }
        index += stride;
    }
    keys
}

#[test]
#[ignore = "RR-11 acceptance: requires VAKINT_LAYER_ROOT pointing to a staged scratch pack"]
fn staged_layer_pack_loads_reduces_onto_its_basis_and_matches_the_base() {
    let started = Instant::now();
    let root = PathBuf::from(std::env::var_os("VAKINT_LAYER_ROOT").expect("VAKINT_LAYER_ROOT"));
    let base_root = std::env::var_os("VAKINT_LAYER_BASE_ROOT").map(PathBuf::from);
    let expected: Option<usize> = std::env::var("VAKINT_LAYER_EXPECTED_MASTERS")
        .ok()
        .map(|v| v.parse().unwrap());
    let samples: usize = std::env::var("VAKINT_LAYER_SAMPLES")
        .ok()
        .map_or(12, |v| v.parse().unwrap());
    let symbolic: usize = std::env::var("VAKINT_LAYER_SYMBOLIC")
        .ok()
        .map_or(2, |v| v.parse().unwrap());
    let report_path = std::env::var_os("VAKINT_LAYER_REPORT").map(PathBuf::from);
    let mut report = json!({"schema": "vakint.rr11-layer-acceptance.v1", "root": root,
        "base_root": base_root, "checks": {}});

    // Base pass first (sequential, to bound memory): reference reductions.
    let mut base_results = Vec::new();
    let mut base_symbolic = Vec::new();
    let mut base_remaining = BTreeSet::new();
    let mut sample_keys = Vec::new();
    if let Some(base_root) = &base_root {
        let base = load(base_root);
        report["base_load"] = base.receipt.clone();
        sample_keys = witnesses(&base.application, samples);
        base_remaining = base
            .application
            .prepared
            .lock()
            .unwrap()
            .remaining_terminals()
            .clone();
        for key in &sample_keys {
            base_results.push(reduce(&base.application, key));
        }
        for key in sample_keys.iter().take(symbolic) {
            let input = super::admission_census_tests::input_for_key(key)
                .unwrap()
                .expression();
            let t = Instant::now();
            let value = base
                .vakint
                .evaluate_rustred_symbolic(&base.settings, input.as_view())
                .map(|e| e.value.to_canonical_string());
            base_symbolic.push((value, t.elapsed().as_secs_f64()));
        }
        drop(base);
    }

    // V1: load.
    let staged = load(&root);
    report["checks"]["V1"] = json!({"passed": true, "load": staged.receipt});
    println!("{}", json!({"event": "V1", "receipt": staged.receipt}));

    // V2: basis and graph inventory.
    let prepared = staged.application.prepared.lock().unwrap();
    let remaining = prepared.remaining_terminals().clone();
    let topology_count = staged
        .vakint
        .topologies
        .five_loop_catalogue_statistics()
        .map(|s| s.topology_count);
    let pack = &staged.data.loops[&5];
    let families: BTreeMap<String, String> = pack
        .manifest
        .families
        .iter()
        .map(|f| {
            (
                f["sha256"].as_str().unwrap().to_owned(),
                f["fingerprint"].as_str().unwrap().to_owned(),
            )
        })
        .collect();
    let mut manifest_keys = BTreeSet::new();
    let mut missing_graphs = Vec::new();
    for master in &pack.manifest.masters {
        let powers: Vec<i64> = serde_json::from_value(master["powers"].clone()).unwrap();
        let fingerprint = &families[master["family_sha256"].as_str().unwrap()];
        manifest_keys.insert((fingerprint.clone(), powers));
        let graph = master["files"]["graph.dot"]["path"].as_str().unwrap_or("");
        if !pack
            .path(graph)
            .ok()
            .flatten()
            .is_some_and(|p| p.is_file())
        {
            missing_graphs.push(master["name"].clone());
        }
    }
    let remaining_keys: BTreeSet<_> = remaining
        .iter()
        .map(|k| (k.family_fingerprint().to_owned(), k.integral().powers().to_vec()))
        .collect();
    let v2 = expected.is_none_or(|n| n == remaining.len())
        && topology_count == Some(67)
        && manifest_keys == remaining_keys
        && missing_graphs.is_empty();
    report["checks"]["V2"] = json!({"passed": v2, "remaining": remaining.len(),
        "expected": expected, "topology_count": topology_count,
        "manifest_masters": pack.manifest.masters.len(),
        "manifest_equals_remaining": manifest_keys == remaining_keys,
        "missing_graphs": missing_graphs});
    drop(prepared);

    // V3 (only basis keys, nothing unresolved) and V4 (old vs new exact).
    if sample_keys.is_empty() {
        sample_keys = witnesses(&staged.application, samples);
    }
    let mut v3_cases = Vec::new();
    let mut v4_cases = Vec::new();
    let mut v3 = true;
    let mut v4 = base_root.is_some();
    let mut source_rows = BTreeMap::<SavedApplicationIntegral, Result<Terms, String>>::new();
    for (index, key) in sample_keys.iter().enumerate() {
        let t = Instant::now();
        let new = reduce(&staged.application, key);
        let seconds = t.elapsed().as_secs_f64();
        let ok3 = matches!(&new, Ok((terms, 0)) if terms.keys().all(|k| remaining.contains(k)));
        v3 &= ok3;
        v3_cases.push(json!({"key": key.powers(), "passed": ok3, "seconds": seconds,
            "terms": new.as_ref().map(|(t, _)| t.len()).ok(),
            "error": new.as_ref().err()}));
        if let Some(Ok((old, 0))) = base_results.get(index) {
            // Substitute every base terminal that left the basis by its
            // reduction in the staged pack (the layer row), exactly.
            let prepared = staged.application.prepared.lock().unwrap();
            let context = prepared.family().coefficient_context().clone();
            drop(prepared);
            let limits = Default::default();
            let mut expected_terms = Terms::new();
            let mut skipped = None;
            for (terminal, coefficient) in old {
                let row = if remaining.contains(terminal) {
                    Ok(Terms::from([(terminal.clone(), context.one())]))
                } else {
                    source_rows
                        .entry(terminal.clone())
                        .or_insert_with(|| {
                            reduce(&staged.application, terminal.integral()).and_then(
                                |(terms, unresolved)| {
                                    if unresolved == 0 {
                                        Ok(terms)
                                    } else {
                                        Err("unresolved".into())
                                    }
                                },
                            )
                        })
                        .clone()
                };
                match row {
                    Err(error) => {
                        skipped = Some(format!(
                            "source {:?}: {error}",
                            terminal.integral().powers()
                        ));
                        break;
                    }
                    Ok(row) => {
                        for (master, b) in row {
                            let value = context.try_mul(coefficient, &b, limits).unwrap();
                            let sum = match expected_terms.remove(&master) {
                                Some(previous) => context.try_add(&previous, &value, limits).unwrap(),
                                None => value,
                            };
                            if !sum.is_zero() {
                                expected_terms.insert(master, sum);
                            }
                        }
                    }
                }
            }
            let equal = skipped.is_none()
                && matches!(&new, Ok((terms, 0)) if *terms == expected_terms);
            v4 &= equal;
            v4_cases.push(json!({"key": key.powers(), "passed": equal, "skipped": skipped,
                "base_terms": old.len(), "base_terms_outside_new_basis":
                    old.keys().filter(|k| !remaining.contains(*k)).count()}));
        } else if base_root.is_some() {
            v4 = false;
            v4_cases.push(json!({"key": key.powers(), "passed": false,
                "base_error": format!("{:?}", base_results.get(index).map(|r| r.as_ref().err()))}));
        }
    }
    // Vakint-level symbolic evaluation (Retain mode) on the first inputs.
    let mut symbolic_cases = Vec::new();
    for (index, key) in sample_keys.iter().take(symbolic).enumerate() {
        let input = super::admission_census_tests::input_for_key(key)
            .unwrap()
            .expression();
        let t = Instant::now();
        let value = staged
            .vakint
            .evaluate_rustred_symbolic(&staged.settings, input.as_view())
            .map(|e| e.value.to_canonical_string());
        let seconds = t.elapsed().as_secs_f64();
        let identical = base_symbolic
            .get(index)
            .map(|(old, _)| old.as_ref().ok() == value.as_ref().ok() && value.is_ok());
        symbolic_cases.push(json!({"key": key.powers(), "ok": value.is_ok(),
            "seconds": seconds, "base_seconds": base_symbolic.get(index).map(|(_, s)| s),
            "identical_to_base": identical, "error": value.as_ref().err().map(|e| e.to_string())}));
        v3 &= value.is_ok();
    }
    // Admitted, non-factorized basis terminals that reduce to themselves:
    // candidates for the exact_basis_expression alias prototype.
    let candidates: Vec<_> = {
        let prepared = staged.application.prepared.lock().unwrap();
        remaining
            .iter()
            .filter(|t| {
                prepared.terminal_product(t).is_none()
                    && prepared.terminal_product_chart(t).is_none()
                    && prepared.routed_entry_is_admitted(t.integral()).unwrap_or(false)
            })
            .take(24)
            .cloned()
            .collect()
    };
    let mut self_reducing = Vec::new();
    for terminal in candidates {
        if self_reducing.len() >= 6 {
            break;
        }
        if let Ok((terms, 0)) = reduce(&staged.application, terminal.integral()) {
            if terms.len() == 1 && terms.contains_key(&terminal) {
                self_reducing.push(json!(terminal.integral().powers()));
            }
        }
    }
    report["self_reducing_terminals"] = json!(self_reducing);
    let basis_changed = !base_remaining.is_empty() && base_remaining != remaining;
    report["checks"]["V3"] = json!({"passed": v3, "cases": v3_cases,
        "symbolic": symbolic_cases});
    report["checks"]["V4"] = json!({"passed": v4, "evaluated": base_root.is_some(),
        "basis_changed": basis_changed, "cases": v4_cases,
        "method": "base reduction, each base terminal outside the new basis substituted by its staged reduction, exact equality with the staged reduction"});
    report["memory_final"] = memory();
    report["seconds"] = json!(started.elapsed().as_secs_f64());
    let text = serde_json::to_string_pretty(&report).unwrap();
    println!("{text}");
    if let Some(path) = report_path {
        std::fs::write(path, &text).unwrap();
    }
    assert!(v2, "V2 failed");
    assert!(v3, "V3 failed");
    assert!(base_root.is_none() || v4, "V4 failed");
}

/// Prototype of the `exact_basis_expression` alias route at 5L (scratch pack
/// prepared by `alias_toy_pack.py`): the basis master carries toy numerical
/// coefficients and the primitive symbol `Q5TOY1`; the aliased master carries
/// no coefficients and `RATIO*Q5TOY1`. Numerical evaluation of the aliased
/// master must then equal RATIO times the basis evaluation; Retain mode must
/// keep the aliased symbol.
/// Env: `VAKINT_ALIAS_ROOT`, `VAKINT_ALIAS_ALIASED` and `VAKINT_ALIAS_BASIS`
/// (JSON power lists), `VAKINT_ALIAS_RATIO` (a rational), `VAKINT_ALIAS_REPORT`.
#[test]
#[ignore = "RR-11 alias prototype: requires VAKINT_ALIAS_ROOT pointing to a scratch toy pack"]
fn exact_basis_expression_alias_route_prototype_on_scratch_pack() {
    let root = PathBuf::from(std::env::var_os("VAKINT_ALIAS_ROOT").expect("VAKINT_ALIAS_ROOT"));
    let powers = |name: &str| -> IntegralKey {
        let values: Vec<i64> = serde_json::from_str(&std::env::var(name).unwrap()).unwrap();
        IntegralKey::try_new(values).unwrap()
    };
    let aliased = powers("VAKINT_ALIAS_ALIASED");
    let basis = powers("VAKINT_ALIAS_BASIS");
    let ratio = crate::utils::vakint_macros::vk_parse!(&std::env::var("VAKINT_ALIAS_RATIO").unwrap())
        .unwrap();
    let mut vakint = Vakint::new().unwrap();
    let started = Instant::now();
    vakint.load_external_data(&root, 5).unwrap();
    let load_seconds = started.elapsed().as_secs_f64();
    let mut report = json!({"root": root, "load_seconds": load_seconds});
    let input_a = super::admission_census_tests::input_for_key(&aliased)
        .unwrap()
        .expression();
    let input_b = super::admission_census_tests::input_for_key(&basis)
        .unwrap()
        .expression();
    let mode_settings = |mode| {
        let mut settings = settings();
        settings.evaluation_order = crate::EvaluationOrder(vec![crate::EvaluationMethod::RustRed(
            crate::RustRedEvaluationOptions { master_mode: mode },
        )]);
        settings
    };
    let retain = mode_settings(crate::RustRedMasterMode::Retain);
    let symbolic = mode_settings(crate::RustRedMasterMode::Symbolic);
    let numerical = mode_settings(crate::RustRedMasterMode::Numerical);
    let show = |value: &Result<Atom, VakintError>| -> Value {
        json!(value.as_ref().map(|a| a.to_canonical_string()).map_err(|e| e.to_string()))
    };
    let mut timed = |label: &str, settings: &VakintSettings, input: &Atom| {
        let t = Instant::now();
        let value = vakint.evaluate_integral(settings, input.as_view());
        report[label] = json!({"seconds": t.elapsed().as_secs_f64(), "value": show(&value)});
        value
    };
    // Retain keeps the 640 symbols (aliases are never applied there).
    let _ = timed("retain_aliased", &retain, &input_a);
    let _ = timed("retain_basis", &retain, &input_b);
    // Symbolic: the aliased master has no coefficients; finalize applies the
    // alias and evaluates the basis master's exact toy coefficients.
    let sym_a = timed("symbolic_aliased", &symbolic, &input_a);
    let sym_b = timed("symbolic_basis", &symbolic, &input_b);
    // Numerical: shown for the constraint (unresolved-master check runs
    // before any alias is tried).
    let _ = timed("numerical_aliased", &numerical, &input_a);
    let mut agreement = None;
    if let (Ok(a), Ok(b)) = (&sym_a, &sym_b) {
        let difference = (a.clone() - ratio.clone() * b.clone()).expand();
        agreement = Some(json!({"exact_zero": difference.is_zero(),
            "difference": difference.to_canonical_string()}));
    }
    report["symbolic_agreement"] = json!(agreement);
    let text = serde_json::to_string_pretty(&report).unwrap();
    println!("{text}");
    if let Some(path) = std::env::var_os("VAKINT_ALIAS_REPORT") {
        std::fs::write(path, &text).unwrap();
    }
    assert!(
        agreement.is_some_and(|a| a["exact_zero"] == true),
        "symbolic alias evaluation must equal RATIO times the basis evaluation"
    );
}
