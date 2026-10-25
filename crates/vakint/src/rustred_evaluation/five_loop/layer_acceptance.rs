//! RR-11 vakint acceptance V1-V4 for a staged (scratch) five-loop pack.
//!
//! Driven by environment variables (all paths are scratch copies; nothing is
//! written into any data checkout):
//! - `VAKINT_LAYER_ROOT`: staged external-data root (`manifest.json`, 1L-5L);
//! - `VAKINT_LAYER_BASE_ROOT`: optional unmodified root for V4 (old vs new);
//! - `VAKINT_LAYER_EXPECTED_MASTERS`: N' (V2);
//! - `VAKINT_LAYER_SAMPLES`: raw terminal witnesses reduced in V3/V4 (default 12);
//! - `VAKINT_LAYER_SYMBOLIC`: vakint-level symbolic evaluations compared (default 2);
//! - `VAKINT_LAYER_REPORT`: JSON receipt path;
//! - `VAKINT_LAYER_FILE`: the layer's `layer.json` with TEXT coefficients (V4);
//! - `VAKINT_LAYER_CHECK_VALUES`: compose-layer `--check-values-out` file whose
//!   `touched_rows` supply the V4 witnesses (raw rows the layer touches);
//! - `VAKINT_LAYER_V4_WITNESSES`: number of touched witnesses (default 16).
//!
//! V4 (old vs new, exact) is independent of the staged pack's own rows: the
//! expected reduction of a witness is its BASE reduction with every layer
//! source replaced by the layer row parsed here from the layer file's text
//! (Symbolica parser, not the rustred-app layer reader), compared exactly with
//! the STAGED reduction. Witnesses are raw rows the layer touches that are
//! not sources and have >= 2 terms (plus the V3 samples).
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

/// Text layer rows: source powers -> [(master powers, coefficient text)].
type TextRows = BTreeMap<Vec<i64>, Vec<(Vec<i64>, String)>>;

fn read_text_layer(path: &Path) -> Result<TextRows, String> {
    let layer: Value = serde_json::from_slice(&std::fs::read(path).map_err(|e| e.to_string())?)
        .map_err(|e| e.to_string())?;
    let powers = |v: &Value| -> Result<Vec<i64>, String> {
        serde_json::from_value(v.clone()).map_err(|e| e.to_string())
    };
    let mut rows = TextRows::new();
    for row in layer["rows"].as_array().ok_or("layer has no rows")? {
        let mut terms = Vec::new();
        for term in row["terms"].as_array().ok_or("row without terms")? {
            let text = term["coefficient"]
                .as_str()
                .ok_or("V4 needs text coefficients (a native-table layer cannot be re-parsed)")?;
            terms.push((powers(&term["master"])?, text.to_owned()));
        }
        if rows.insert(powers(&row["source"])?, terms).is_some() {
            return Err("duplicate layer source".into());
        }
    }
    Ok(rows)
}

/// Parse a coefficient text into the family's coefficient context with
/// Symbolica directly (independent of the rustred-app layer reader).
fn parse_text_coefficient(
    context: &rustred::algebra::CoefficientContext,
    text: &str,
) -> Result<Coefficient, String> {
    use symbolica::prelude::{Q, Z};
    let atom = symbolica::try_parse!(text, default_namespace = "rustred")
        .map_err(|e| format!("{text}: {e}"))?;
    let coefficient = atom
        .as_view()
        .try_to_rational_polynomial(&Q, &Z, Some(context.one().get_variables().clone()))
        .map_err(|e| format!("{text}: {e}"))?;
    if !context.contains(&coefficient) {
        return Err(format!("{text}: symbols outside the family parameters"));
    }
    Ok(coefficient)
}

/// Touched raw rows (not sources, >= 2 terms) from the check-values file:
/// the rows with the most source terms and the most terms first, then a
/// stratified stride over the rest (sorted keys).
fn touched_witnesses(path: &Path, count: usize) -> Result<(Vec<IntegralKey>, Value), String> {
    let values: Value = serde_json::from_slice(&std::fs::read(path).map_err(|e| e.to_string())?)
        .map_err(|e| e.to_string())?;
    let mut eligible: Vec<(Vec<i64>, u64, u64)> = Vec::new();
    let all = values["touched_rows"].as_array().ok_or("no touched_rows")?;
    for row in all {
        let key: Vec<i64> =
            serde_json::from_value(row["key"].clone()).map_err(|e| e.to_string())?;
        let terms = row["terms"].as_u64().unwrap_or(0);
        let sources = row["source_terms"].as_u64().unwrap_or(0);
        if row["is_source"] == false && terms >= 2 && sources >= 1 {
            eligible.push((key, terms, sources));
        }
    }
    eligible.sort();
    let mut chosen: Vec<Vec<i64>> = Vec::new();
    if let Some(row) = eligible.iter().max_by_key(|r| (r.2, r.1)) {
        chosen.push(row.0.clone());
    }
    if let Some(row) = eligible.iter().max_by_key(|r| (r.1, r.2)) {
        if !chosen.contains(&row.0) {
            chosen.push(row.0.clone());
        }
    }
    let stride = (eligible.len() / count.max(1)).max(1);
    let mut index = 0;
    while chosen.len() < count && index < eligible.len() {
        if !chosen.contains(&eligible[index].0) {
            chosen.push(eligible[index].0.clone());
        }
        index += stride;
    }
    let summary = json!({"touched_rows": all.len(), "eligible": eligible.len(),
        "chosen": chosen.len(), "check_values_layer_blake3": values["layer_blake3"]});
    let keys = chosen
        .into_iter()
        .map(|powers| IntegralKey::try_new(powers).map_err(|e| e.to_string()))
        .collect::<Result<_, _>>()?;
    Ok((keys, summary))
}

/// Terms keyed by (family fingerprint, powers): new layer masters cannot be
/// built as `SavedApplicationIntegral` outside rustred-app.
type PlainTerms = BTreeMap<(String, Vec<i64>), Coefficient>;

fn plain(terms: &Terms) -> PlainTerms {
    terms
        .iter()
        .map(|(k, c)| {
            (
                (
                    k.family_fingerprint().to_owned(),
                    k.integral().powers().to_vec(),
                ),
                c.clone(),
            )
        })
        .collect()
}

/// Expected staged reduction: the base reduction with every layer source
/// replaced by its (independently parsed) text row; base terminals that are
/// neither sources nor staged basis elements are an error.
fn expected_from_text_layer(
    old: &Terms,
    layer: &TextRows,
    remaining: &BTreeSet<SavedApplicationIntegral>,
    context: &rustred::algebra::CoefficientContext,
    parsed: &mut BTreeMap<(Vec<i64>, usize), Coefficient>,
) -> Result<(PlainTerms, usize), String> {
    let limits = Default::default();
    let mut expected = PlainTerms::new();
    let mut source_terms = 0;
    let add = |expected: &mut PlainTerms, key: (String, Vec<i64>), value: Coefficient| {
        let sum = match expected.remove(&key) {
            Some(previous) => context.try_add(&previous, &value, limits).unwrap(),
            None => value,
        };
        if !sum.is_zero() {
            expected.insert(key, sum);
        }
    };
    for (terminal, coefficient) in old {
        let powers = terminal.integral().powers().to_vec();
        if let Some(row) = layer.get(&powers) {
            source_terms += 1;
            for (index, (master, text)) in row.iter().enumerate() {
                let b = match parsed.get(&(powers.clone(), index)) {
                    Some(b) => b.clone(),
                    None => {
                        let b = parse_text_coefficient(context, text)?;
                        parsed.insert((powers.clone(), index), b.clone());
                        b
                    }
                };
                let key = (terminal.family_fingerprint().to_owned(), master.clone());
                add(
                    &mut expected,
                    key,
                    context.try_mul(coefficient, &b, limits).unwrap(),
                );
            }
        } else if remaining.contains(terminal) {
            add(
                &mut expected,
                (terminal.family_fingerprint().to_owned(), powers.clone()),
                coefficient.clone(),
            );
        } else {
            return Err(format!(
                "base terminal {powers:?} is neither a layer source nor a staged basis element"
            ));
        }
    }
    Ok((expected, source_terms))
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
    let layer_file = std::env::var_os("VAKINT_LAYER_FILE").map(PathBuf::from);
    let check_values = std::env::var_os("VAKINT_LAYER_CHECK_VALUES").map(PathBuf::from);
    let touched_count: usize = std::env::var("VAKINT_LAYER_V4_WITNESSES")
        .ok()
        .map_or(16, |v| v.parse().unwrap());
    let mut report = json!({"schema": "vakint.rr11-layer-acceptance.v2", "root": root,
        "base_root": base_root, "layer_file": layer_file, "check_values": check_values,
        "checks": {}});
    let text_layer = layer_file.as_deref().map(|p| read_text_layer(p).unwrap());
    let (touched_keys, touched_summary) = match &check_values {
        Some(path) => {
            let (keys, summary) = touched_witnesses(path, touched_count).unwrap();
            (keys, summary)
        }
        None => (Vec::new(), Value::Null),
    };

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
        // V3 samples first, then the touched V4 witnesses (deduplicated).
        for key in &touched_keys {
            if !sample_keys.contains(key) {
                sample_keys.push(key.clone());
            }
        }
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
        if !pack.path(graph).ok().flatten().is_some_and(|p| p.is_file()) {
            missing_graphs.push(master["name"].clone());
        }
    }
    let remaining_keys: BTreeSet<_> = remaining
        .iter()
        .map(|k| {
            (
                k.family_fingerprint().to_owned(),
                k.integral().powers().to_vec(),
            )
        })
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
    let mut v4 = base_root.is_some() && text_layer.is_some();
    let mut parsed = BTreeMap::new();
    let mut touched_checked = 0usize;
    for (index, key) in sample_keys.iter().enumerate() {
        let t = Instant::now();
        let new = reduce(&staged.application, key);
        let seconds = t.elapsed().as_secs_f64();
        let ok3 = matches!(&new, Ok((terms, 0)) if terms.keys().all(|k| remaining.contains(k)));
        v3 &= ok3;
        v3_cases.push(
            json!({"key": key.powers(), "passed": ok3, "seconds": seconds,
            "terms": new.as_ref().map(|(t, _)| t.len()).ok(),
            "error": new.as_ref().err()}),
        );
        let touched = touched_keys.contains(key);
        match (base_results.get(index), &text_layer) {
            (Some(Ok((old, 0))), Some(layer)) => {
                let context = staged
                    .application
                    .prepared
                    .lock()
                    .unwrap()
                    .family()
                    .coefficient_context()
                    .clone();
                let outcome =
                    expected_from_text_layer(old, layer, &remaining, &context, &mut parsed);
                let (equal, source_terms, error) = match &outcome {
                    Ok((expected, sources)) => (
                        matches!(&new, Ok((terms, 0)) if plain(terms) == *expected),
                        Some(*sources),
                        None,
                    ),
                    Err(e) => (false, None, Some(e.clone())),
                };
                // A touched witness must really touch the layer in the base.
                let witness_ok = !touched || source_terms.is_some_and(|n| n >= 1);
                v4 &= equal && witness_ok;
                if touched && source_terms.is_some_and(|n| n >= 1) {
                    touched_checked += 1;
                }
                v4_cases.push(json!({"key": key.powers(), "touched_witness": touched,
                    "passed": equal && witness_ok, "base_terms": old.len(),
                    "base_layer_source_terms": source_terms,
                    "is_layer_source": layer.contains_key(key.powers()),
                    "staged_terms": new.as_ref().map(|(t, _)| t.len()).ok(), "error": error}));
            }
            (Some(Ok((_, 0))), None) => {}
            _ if base_root.is_some() => {
                v4 = false;
                v4_cases.push(json!({"key": key.powers(), "touched_witness": touched,
                    "passed": false, "base_error":
                        format!("{:?}", base_results.get(index).map(|r| r.as_ref().err()))}));
            }
            _ => {}
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
                    && prepared
                        .routed_entry_is_admitted(t.integral())
                        .unwrap_or(false)
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
    // Require the touched witnesses that were available (at least one when
    // the layer touches any non-source row).
    v4 &= touched_keys.is_empty() || touched_checked == touched_keys.len();
    v4 &= check_values.is_some();
    report["checks"]["V4"] = json!({"passed": v4,
        "evaluated": base_root.is_some() && text_layer.is_some(),
        "independent": true, "basis_changed": basis_changed,
        "touched_witnesses": touched_summary, "touched_checked": touched_checked,
        "text_rows_parsed": parsed.len(), "cases": v4_cases,
        "method": "expected = BASE reduction with every layer source replaced by its row parsed \
            from the layer file text (Symbolica, independent of the staged pack and of the \
            rustred-app layer reader); exact equality with the STAGED reduction; witnesses = \
            raw rows the layer touches (not sources, >= 2 terms) from compose-layer \
            check-values, plus the V3 samples"});
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

/// V5 survey (non-panicking variant of the generation-54 factorized test):
/// every remaining terminal with a product chart or a scalar product is
/// reduced through the lower-loop product route; per-terminal outcome.
/// Env: `VAKINT_LAYER_ROOT`, `VAKINT_LAYER_REPORT`.
#[test]
#[ignore = "RR-11 acceptance: requires VAKINT_LAYER_ROOT pointing to a staged scratch pack"]
fn staged_layer_pack_factorized_terminals_product_route_survey() {
    use rustred::reduction::terminal_normalization::VacuumProductProjection;
    let started = Instant::now();
    let root = PathBuf::from(std::env::var_os("VAKINT_LAYER_ROOT").expect("VAKINT_LAYER_ROOT"));
    let report_path = std::env::var_os("VAKINT_LAYER_REPORT").map(PathBuf::from);
    let loaded = load(&root);
    let application = loaded.application.clone();
    let prepared = application.prepared.lock().unwrap();
    let mut reducer = application.product_reducer.lock().unwrap();
    let mut cases = Vec::new();
    for terminal in prepared.remaining_terminals() {
        let t = Instant::now();
        let mut crosses = false;
        let (route, projection) = if let Some(chart) = prepared.terminal_product_chart(terminal) {
            let projection = chart
                .project(
                    S.dot,
                    (1..=chart.combined_family().loop_count())
                        .map(|axis| S.k.call(axis))
                        .collect(),
                    Default::default(),
                    |input| {
                        crosses = true;
                        factor_projection::project(input)
                    },
                )
                .map_err(|e| format!("{e:?}"));
            ("chart", projection)
        } else if let Some(product) = prepared.terminal_product(terminal) {
            let factors: Result<Vec<_>, String> = product
                .factors()
                .iter()
                .map(|factor| {
                    Ok(VacuumFactor {
                        family: prepared
                            .output_family_shared(factor.family_fingerprint())
                            .ok_or_else(|| "factor family absent".to_string())?,
                        integral: factor.integral().clone(),
                    })
                })
                .collect();
            (
                "product",
                factors.map(|factors| VacuumProductProjection {
                    terms: vec![(factors, product.coefficient().clone())],
                    nonzero_conditions: Vec::new(),
                }),
            )
        } else {
            continue;
        };
        let outcome = projection.and_then(|projection| {
            reducer
                .reduce_products(projection.terms, projection.nonzero_conditions)
                .map_err(|e| format!("{e:?}"))
        });
        let error = outcome
            .as_ref()
            .err()
            .map(|e| e.chars().take(600).collect::<String>());
        cases.push(json!({"key": terminal.integral().powers(), "route": route,
            "cross_projected": crosses, "ok": outcome.is_ok(),
            "product_terms": outcome.as_ref().map(|r| r.terms().len()).ok(),
            "nonzero_conditions": outcome.as_ref().map(|r| r.nonzero_conditions().len()).ok(),
            "error": error, "seconds": t.elapsed().as_secs_f64()}));
        let _ = reducer.clear_cache();
    }
    let failed = cases.iter().filter(|c| c["ok"] == false).count();
    let report = json!({"schema": "vakint.rr11-v5-survey.v1", "root": root,
        "load": loaded.receipt, "remaining": prepared.remaining_terminals().len(),
        "factorized": cases.len(), "failed": failed, "cases": cases,
        "seconds": started.elapsed().as_secs_f64()});
    let text = serde_json::to_string_pretty(&report).unwrap();
    println!("{text}");
    if let Some(path) = report_path {
        std::fs::write(path, &text).unwrap();
    }
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
    let ratio =
        crate::utils::vakint_macros::vk_parse!(&std::env::var("VAKINT_ALIAS_RATIO").unwrap())
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
        json!(
            value
                .as_ref()
                .map(|a| a.to_canonical_string())
                .map_err(|e| e.to_string())
        )
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
