//! Read-only validation of a local five-loop external-data replacement.
//!
//! cargo run --release -p vakint --example external_data_acceptance -- DATA REPORT

#[path = "external_data_acceptance/dot.rs"]
mod dot;

use std::collections::{BTreeMap, BTreeSet};
use std::fs;
use std::path::PathBuf;
use std::time::Instant;

use rustred::family::IntegralKey;
use rustred_app::{SavedApplicationOptions, load_saved_application};
use serde_json::{Value, json};
use sha2::{Digest, Sha256};
use symbolica::atom::{Atom, AtomCore, AtomView};
use vakint::{
    EvaluationMethod, EvaluationOrder, LoopNormalizationFactor, RustRedEvaluationOptions,
    RustRedMasterMode, Vakint, VakintSettings, vakint_parse, vakint_symbol,
};

type Result<T> = std::result::Result<T, Box<dyn std::error::Error>>;

fn main() -> Result<()> {
    let mut args = std::env::args_os().skip(1);
    let root = PathBuf::from(args.next().ok_or("missing external-data directory")?);
    let report = PathBuf::from(args.next().ok_or("missing output report")?);
    let graphs_only = match args.next() {
        Some(value) if value == "--graphs-only" => true,
        None => false,
        Some(_) => return Err("unknown option; expected --graphs-only".into()),
    };
    if args.next().is_some() || report.exists() {
        return Err("usage: external_data_acceptance DATA NEW_REPORT [--graphs-only]".into());
    }
    Vakint::initialize_vakint_symbols();
    let started = Instant::now();
    let manifest_bytes = fs::read(root.join("5L/manifest.json"))?;
    let manifest: Value = serde_json::from_slice(&manifest_bytes)?;
    let masters = manifest["masters"].as_array().ok_or("missing masters")?;
    if masters.len() != 640 {
        return Err("expected 640 five-loop masters".into());
    }
    for master in masters {
        if master["status"] != "not_computed"
            || master["coefficients"]
                .as_array()
                .is_none_or(|values| !values.is_empty())
        {
            return Err("five-loop master has unexpected numerical coefficients".into());
        }
    }
    eprintln!("validating all 640 DOT graphs through native FeynKit and Vakint ingress");
    let graphs_started = Instant::now();
    let graph_report = dot::validate(&root.join("5L"), &manifest)?;
    let mut output = json!({
        "schema": "vakint.external-data.acceptance.v1",
        "measurement_kind": "validation run; timings are not performance benchmarks",
        "data": root,
        "loop_manifest_sha256": format!("{:x}", Sha256::digest(&manifest_bytes)),
        "graphs": graph_report,
        "remaining_masters": masters.len(),
        "cache_available": manifest["cache"]["available"],
        "all_five_loop_masters_unevaluated": true,
    });
    output["graphs"]["elapsed_seconds"] = json!(graphs_started.elapsed().as_secs_f64());
    fs::write(
        report.with_extension("graphs.json"),
        serde_json::to_vec_pretty(&output)?,
    )?;
    if graphs_only {
        output["elapsed_seconds"] = json!(started.elapsed().as_secs_f64());
        fs::write(report, serde_json::to_vec_pretty(&output)?)?;
        return Ok(());
    }

    eprintln!("loading portable native application with strict starting-domain admission");
    let native_started = Instant::now();
    let native = load_saved_application::<15>(
        &root.join("5L/ibps"),
        SavedApplicationOptions {
            preparation_workers: 16,
            ..Default::default()
        },
    )?;
    let native_load_seconds = native_started.elapsed().as_secs_f64();
    eprintln!(
        "portable native load completed in {native_load_seconds:.3}s with 16 preparation workers"
    );
    let families = manifest["families"]
        .as_array()
        .ok_or("missing families")?
        .iter()
        .map(|family| {
            Ok((
                family["sha256"].as_str().ok_or("missing family SHA256")?,
                family["fingerprint"]
                    .as_str()
                    .ok_or("missing family fingerprint")?,
            ))
        })
        .collect::<Result<BTreeMap<_, _>>>()?;
    let expected = masters
        .iter()
        .map(|master| {
            Ok((
                families[master["family_sha256"]
                    .as_str()
                    .ok_or("missing master family")?]
                .to_owned(),
                serde_json::from_value::<Vec<i64>>(master["powers"].clone())?,
            ))
        })
        .collect::<Result<BTreeSet<_>>>()?;
    let actual = native
        .remaining_terminals()
        .iter()
        .map(|terminal| {
            (
                terminal.family_fingerprint().to_owned(),
                terminal.integral().powers().to_vec(),
            )
        })
        .collect::<BTreeSet<_>>();
    if actual != expected || native.terminal_sources().len() != 939 {
        return Err("portable native terminal roster disagrees with the manifest".into());
    }
    let canary = native
        .remaining_terminals()
        .iter()
        .filter(|terminal| {
            terminal.integral().powers().iter().all(|power| *power >= 0)
                && native.validate_entry(terminal.integral()).is_ok()
        })
        .min_by_key(|terminal| terminal.integral().powers().iter().sum::<i64>())
        .ok_or("no admitted scalar master canary")?;
    let canary = masters
        .iter()
        .find(|master| {
            families[master["family_sha256"].as_str().unwrap()] == canary.family_fingerprint()
                && serde_json::from_value::<Vec<i64>>(master["powers"].clone()).unwrap()
                    == canary.integral().powers()
        })
        .ok_or("native canary absent from master definitions")?;
    let native_canary = dot::input(&root.join("5L"), canary)?;
    let queries: Value = serde_json::from_str(native.queries_json())?;
    let scalar_d_values = queries["queries"]
        .as_array()
        .ok_or("missing recorded queries")?
        .iter()
        .filter_map(|query| query["power_bounds"]["min_power_difference"].as_i64())
        .collect::<BTreeSet<_>>();
    let mut banana = None;
    for master in masters {
        if !dot::is_scalar_banana(&root.join("5L"), master)? {
            continue;
        }
        let base: Vec<i64> = serde_json::from_value(master["powers"].clone())?;
        let sum = base.iter().sum::<i64>();
        let slot = base
            .iter()
            .position(|power| *power > 0)
            .ok_or("banana has no positive line")?;
        for d in &scalar_d_values {
            if *d < sum {
                continue;
            }
            let mut powers = base.clone();
            powers[slot] += d - sum;
            let key = IntegralKey::try_new(powers.clone())?;
            if native.validate_entry(&key).is_ok() {
                banana = Some((
                    master["name"].clone(),
                    dot::input_powers(&root.join("5L"), master, Some(&powers))?,
                    powers,
                ));
                break;
            }
        }
        if banana.is_some() {
            break;
        }
    }
    let banana = banana.ok_or("no admitted dotted scalar two-vertex banana control")?;
    let unsupported = IntegralKey::try_new(vec![100; native.family().denominators().len()])?;
    let outside_error = native
        .validate_entry(&unsupported)
        .expect_err("out-of-scope key must be refused");
    output["native"] = json!({
        "load_seconds": native_load_seconds,
        "preparation_workers": 16,
        "source_maps": native.terminal_sources().len(),
        "remaining_terminals": actual.len(),
        "all_output_keys_match_manifest": true,
        "outside_scope_error": outside_error.to_string(),
        "dotted_banana_input_powers": banana.2,
        "metadata": native.metadata(),
    });
    fs::write(
        report.with_extension("native.json"),
        serde_json::to_vec_pretty(&output)?,
    )?;
    drop(native);

    eprintln!("loading all 1L through 5L packs through Vakint::load_external_data");
    let load_started = Instant::now();
    let mut engine = Vakint::new()?;
    engine.load_external_data(&root, 5)?;
    let load_seconds = load_started.elapsed().as_secs_f64();
    eprintln!(
        "Vakint cumulative load completed in {:.3}s (public loader defaults)",
        load_seconds
    );
    let inspection_started = Instant::now();
    for metadata in masters {
        let symbol = metadata["master_id"]
            .as_str()
            .ok_or("missing master identity")?
            .parse()?;
        let master = engine.get_master(&symbol, 16)?;
        if master.loop_count() != 5 || master.is_computed() || !master.coefficients().is_empty() {
            return Err("public inspection misclassified an uncomputed five-loop master".into());
        }
        let dot = master.graph().to_dot()?;
        let restored = feynkit_graph::FeynmanDiagram::from_dot(master.graph().model_arc(), &dot)?;
        if restored.to_dot()? != dot {
            return Err("public master graph changed during canonical DOT round-trip".into());
        }
    }
    output["public_inspection"] = json!({
        "masters_checked": masters.len(),
        "uncomputed_coefficient_maps_empty": true,
        "canonical_dot_roundtrip": true,
        "seconds": inspection_started.elapsed().as_secs_f64(),
    });
    let mut settings = VakintSettings {
        evaluation_order: EvaluationOrder(vec![EvaluationMethod::RustRed(
            RustRedEvaluationOptions {
                master_mode: RustRedMasterMode::Retain,
            },
        )]),
        integral_normalization_factor: LoopNormalizationFactor::FMFTandMATAD,
        form_exe_path: "/640-master-acceptance-must-not-use-FORM".into(),
        ..Default::default()
    };
    let first_master = &masters[0];
    let powers: Vec<i64> = serde_json::from_value(first_master["powers"].clone())?;
    let master_symbol = vakint_parse!(format!(
        "RustRedMaster(rustred_family_{},{})",
        first_master["family_sha256"].as_str().unwrap(),
        powers
            .iter()
            .map(ToString::to_string)
            .collect::<Vec<_>>()
            .join(","),
    ))?;
    let unavailable = engine
        .numerical_evaluation(
            &settings,
            master_symbol.as_view(),
            &Default::default(),
            &Default::default(),
            None,
        )
        .expect_err("five-loop numerical evaluation must report absent coefficients")
        .to_string();
    if !unavailable.contains("declared terminal has no computed coefficients")
        || !unavailable.contains(first_master["master_id"].as_str().unwrap())
    {
        return Err(format!("unexpected numerical-unavailable diagnostic: {unavailable}").into());
    }
    let mut declared_symbols = BTreeSet::new();
    for loop_count in 1..=5 {
        let pack: Value = serde_json::from_slice(&fs::read(
            root.join(format!("{loop_count}L/manifest.json")),
        )?)?;
        for master in pack["masters"]
            .as_array()
            .ok_or("missing lower-loop masters")?
        {
            let powers: Vec<i64> = serde_json::from_value(master["powers"].clone())?;
            declared_symbols.insert(vakint_parse!(format!(
                "RustRedMaster(rustred_family_{},{})",
                master["family_sha256"]
                    .as_str()
                    .ok_or("missing lower-loop family")?,
                powers
                    .iter()
                    .map(ToString::to_string)
                    .collect::<Vec<_>>()
                    .join(",")
            ))?);
        }
    }
    let mut canary_reports = Vec::new();
    for (name, input, kind) in [
        (
            canary["name"].clone(),
            native_canary,
            "admitted minimum-power scalar",
        ),
        (
            banana.0,
            banana.1,
            "admitted dotted scalar two-vertex banana, one cycle-matroid block",
        ),
    ] {
        eprintln!("running five-loop {kind} control {name} and exact warm repeat");
        let first_started = Instant::now();
        let first = engine.evaluate_integral(&settings, input.as_view())?;
        let first_seconds = first_started.elapsed().as_secs_f64();
        let warm_started = Instant::now();
        let warm = engine.evaluate_integral(&settings, input.as_view())?;
        let warm_seconds = warm_started.elapsed().as_secs_f64();
        if !(&first - &warm).expand().together().is_zero() {
            return Err("five-loop repeated symbolic reduction differs exactly".into());
        }
        let mut output_masters = BTreeSet::new();
        first.visitor(&mut |part| {
            if let AtomView::Fun(function) = part
                && function.get_symbol() == vakint_symbol!("RustRedMaster")
            {
                output_masters.insert(part.to_owned());
            }
            true
        });
        if output_masters.is_empty() || !output_masters.is_subset(&declared_symbols) {
            return Err("symbolic canary contains undeclared or no master integrals".into());
        }
        canary_reports.push(json!({
            "master": name, "kind": kind,
            "input": input.to_canonical_string(), "output": first.to_canonical_string(),
            "all_output_masters_declared": true,
            "output_master_symbols": output_masters.iter().map(Atom::to_canonical_string).collect::<Vec<_>>(),
            "first_seconds": first_seconds, "warm_seconds": warm_seconds,
            "warm_exact_agreement": true,
        }));
    }
    output["vakint"] = json!({
        "load_seconds": load_seconds,
        "five_loop_catalogue": format!("{:?}", engine.topologies.five_loop_catalogue_statistics()),
        "numerical_unavailable": unavailable,
        "canaries": canary_reports,
        "recorded_scope_policy": "native strict admission checked separately; Vakint retains its existing Preserve policy",
        "cache_statistics": format!("{:?}", engine.rustred_application_cache_statistics()?),
        "cache_work_statistics": format!("{:?}", engine.rustred_application_cache_work_statistics()?),
    });
    let lower_inputs = [
        "topo(I1L(1,1))",
        "topo(I2L(1,1,1,1))",
        "topo(I3L(1,1,1,1,1,1,1))",
    ];
    let mut lower = Vec::new();
    for source in lower_inputs {
        let input = vakint_parse!(source)?;
        let now = Instant::now();
        let value = engine.evaluate_integral(&settings, input.as_view())?;
        if value.is_zero() {
            return Err(format!("unexpected zero for lower-loop control {source}").into());
        }
        lower.push(json!({"input": source, "seconds": now.elapsed().as_secs_f64(), "output": value.to_canonical_string()}));
    }
    let cases: Value =
        serde_json::from_str(include_str!("../tests/fixtures/benchmark_reduction.json"))?;
    for case in cases.as_array().unwrap().iter().take(4) {
        let input = Atom::parse(
            case["input"].as_str().unwrap(),
            "vakint",
            Default::default(),
        )?;
        let now = Instant::now();
        let value = engine.evaluate_integral(&settings, input.as_view())?;
        if value.is_zero() {
            return Err("unexpected zero for four-loop control".into());
        }
        lower.push(json!({"case": case["case"], "family": case["family"], "seconds": now.elapsed().as_secs_f64(), "output": value.to_canonical_string()}));
    }
    settings.allow_unknown_integrals = false;
    let malformed = vakint_parse!("topo(I99L(1,1))")?;
    output["unsupported_input_error"] = json!(
        engine
            .evaluate_integral(&settings, malformed.as_view())
            .expect_err("unsupported topology must fail")
            .to_string()
    );
    output["lower_loop_controls"] = json!(lower);
    output["elapsed_seconds"] = json!(started.elapsed().as_secs_f64());
    fs::write(report, serde_json::to_vec_pretty(&output)?)?;
    println!("all 640 graph definitions and the native/Vakint package acceptance checks passed");
    Ok(())
}
