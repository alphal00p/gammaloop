//! V7 pinch consistency (track C6) for a staged five-loop pack.
//!
//! A case pairs LHS = (cancelling numerator D_i...) * I(a), written in another
//! routing of the same graph (a parent automorphism k -> M k), with RHS = the
//! pinched/lowered integral I(a - e_i...) in the original routing. Both sides
//! are lowered by RustRed's scalar-numerator service (which performs the
//! numerator/denominator cancellation algebraically) and reduced through the
//! loaded saved application (the call vakint's `evaluate_detailed` makes).
//! Their coefficient vectors over the remaining terminals are compared
//! exactly (exact rational functions of d; `try_sub(..).is_zero()`).
//! A subset also runs vakint's full frontend (`evaluate_rustred_symbolic`).
//!
//! Environment: `VAKINT_V7_ROOT` (staged root), `VAKINT_V7_CASES` (case list
//! from v7-pinch/scripts/gen_cases.py), `VAKINT_V7_REPORT` (JSONL, appended;
//! cases already present are skipped, so a rerun resumes), optional
//! `VAKINT_V7_BUDGET_SECONDS` (default 1500), `VAKINT_V7_VK` (0 disables the
//! frontend path), `VAKINT_V7_ONLY` (comma-separated case ids), `VAKINT_V7_CASE_TIMEOUT` (seconds,
//! default 0 = none): a watchdog records a case that exceeds it as outcome "timeout" and EXITS the process
//! with code 3 (a reduction cannot be interrupted safely); the caller reruns to resume after that case.
use super::*;
use rustred_app::SavedApplicationIntegral;
use serde_json::{Value, json};
use std::io::Write as _;
use std::time::Instant;

use crate::topologies::five_loop::data as census;

type Terms = BTreeMap<SavedApplicationIntegral, Coefficient>;

fn settings() -> VakintSettings {
    VakintSettings {
        form_exe_path: "/five-loop-v7-must-not-use-form".into(),
        evaluation_order: crate::EvaluationOrder(vec![crate::EvaluationMethod::RustRed(
            crate::RustRedEvaluationOptions {
                master_mode: crate::RustRedMasterMode::Retain,
            },
        )]),
        ..Default::default()
    }
}

fn momentum(row: &[i64]) -> Atom {
    row.iter()
        .enumerate()
        .filter(|(_, c)| **c != 0)
        .map(|(axis, c)| Atom::num(*c) * S.k.call(axis + 1))
        .sum::<Atom>()
}

struct Side {
    target: IntegralKey,
    numerator: Atom,
}

fn side(spec: &Value) -> Result<Side, String> {
    let target: Vec<i64> =
        serde_json::from_value(spec["target"].clone()).map_err(|e| e.to_string())?;
    let mut numerator = Atom::one();
    for factor in spec["numerator"].as_array().ok_or("numerator list")? {
        let row: Vec<i64> =
            serde_json::from_value(factor[0].clone()).map_err(|e| e.to_string())?;
        let exponent = factor[1].as_u64().ok_or("exponent")? as u32;
        let q = momentum(&row);
        numerator *= (S.dot.call((&q, &q)) - Atom::one()).pow(exponent);
    }
    Ok(Side {
        target: IntegralKey::try_new(target).map_err(|e| e.to_string())?,
        numerator,
    })
}

/// vakint input topo(...) * numerator for a target, on the first compiled
/// parent that contains the positive support (zero-power rows retained), as
/// `admission_census_tests::input_for_key` builds it.
fn vk_input(side: &Side) -> Result<Atom, String> {
    let powers = side.target.powers();
    let parent = census::PARENTS
        .iter()
        .find(|parent| {
            powers
                .iter()
                .enumerate()
                .all(|(slot, power)| *power <= 0 || parent.physical_slots.contains(&slot))
        })
        .ok_or("no compiled parent contains the support")?;
    let rows = common_momentum_rows();
    let mut integral = Atom::one();
    for (physical, slot) in parent.physical_slots.iter().enumerate() {
        let signed = (*slot + 1) as i64;
        let source = parent
            .vertices
            .iter()
            .position(|v| v.contains(&signed))
            .ok_or("source incidence")?;
        let target = parent
            .vertices
            .iter()
            .position(|v| v.contains(&-signed))
            .ok_or("target incidence")?;
        integral *= S.prop.call((
            physical + 1,
            S.edge.call((source + 1, target + 1)),
            &momentum(&rows[*slot]),
            Atom::one(),
            Atom::num(powers[*slot]),
        ));
    }
    Ok(&side.numerator * &S.topo.call(integral))
}

struct Lowered {
    terms: Vec<(IntegralKey, Coefficient)>,
    not_admitted: Option<Vec<i64>>,
}

fn lower(application: &SavedApplication, side: &Side) -> Result<Lowered, String> {
    let prepared = application.prepared.lock().unwrap();
    let scalar = FamilyScalarNumeratorService::try_new(
        prepared.family(),
        S.dot,
        (1..=5).map(|axis| S.k.call(axis)).collect(),
        Default::default(),
    )
    .map_err(|e| e.to_string())?;
    let lowered = scalar
        .lower_weighted(&side.numerator, &side.target, &Atom::one())
        .map_err(|e| e.to_string())?;
    let mut terms = Vec::new();
    let mut not_admitted = None;
    for term in lowered.terms() {
        if term.scalar_spectator() != &Atom::one() {
            return Err(format!("non-unit scalar spectator {}", term.scalar_spectator()));
        }
        if not_admitted.is_none()
            && !prepared
                .routed_entry_is_admitted(term.integral())
                .map_err(|e| e.to_string())?
        {
            not_admitted = Some(term.integral().powers().to_vec());
        }
        terms.push((term.integral().clone(), term.coefficient().clone()));
    }
    Ok(Lowered {
        terms,
        not_admitted,
    })
}

fn reduce(
    application: &SavedApplication,
    terms: Vec<(IntegralKey, Coefficient)>,
) -> Result<Terms, String> {
    let mut prepared = application.prepared.lock().unwrap();
    let reduced = prepared
        .reduce_routed_linear_combination(terms)
        .map_err(|e| e.to_string())?;
    if !reduced.unresolved().is_empty() {
        return Err(format!("{} unresolved terminals", reduced.unresolved().len()));
    }
    Ok(reduced.terms().clone())
}

fn keys(terms: &Terms) -> Vec<Vec<i64>> {
    terms.keys().map(|k| k.integral().powers().to_vec()).collect()
}

/// Exact difference LHS - RHS over the remaining terminals (nonzero entries).
fn difference(
    application: &SavedApplication,
    left: &Terms,
    right: &Terms,
) -> Result<Vec<Value>, String> {
    let context = application
        .prepared
        .lock()
        .unwrap()
        .family()
        .coefficient_context()
        .clone();
    let limits = Default::default();
    let mut out = Vec::new();
    let all: BTreeSet<&SavedApplicationIntegral> = left.keys().chain(right.keys()).collect();
    for key in all {
        let d = match (left.get(key), right.get(key)) {
            (Some(a), Some(b)) => context.try_sub(a, b, limits).map_err(|e| e.to_string())?,
            (Some(a), None) => a.clone(),
            (None, Some(b)) => context.try_neg(b, limits).map_err(|e| e.to_string())?,
            (None, None) => unreachable!(),
        };
        if !d.is_zero() {
            let text = d.to_expression().to_canonical_string();
            // full coefficient (no truncation): the C6 gate applies later-staged rows to the difference exactly
            out.push(json!({"terminal": key.integral().powers(), "difference": text}));
        }
    }
    Ok(out)
}

fn rr_case(application: &SavedApplication, case: &Value) -> Value {
    let started = Instant::now();
    let run = || -> Result<Value, String> {
        let lhs = side(&case["lhs"])?;
        let rhs = side(&case["rhs"])?;
        let l = lower(application, &lhs)?;
        let r = lower(application, &rhs)?;
        let lhs_keys: BTreeSet<Vec<i64>> = l.terms.iter().map(|(k, _)| k.powers().to_vec()).collect();
        let rhs_keys: BTreeSet<Vec<i64>> = r.terms.iter().map(|(k, _)| k.powers().to_vec()).collect();
        let nontrivial = lhs_keys != rhs_keys;
        if l.not_admitted.is_some() || r.not_admitted.is_some() {
            return Ok(json!({"outcome": "not_admitted", "lhs_not_admitted": l.not_admitted,
                "rhs_not_admitted": r.not_admitted, "lhs_lowered": lhs_keys.len(),
                "rhs_lowered": rhs_keys.len(), "nontrivial": nontrivial}));
        }
        let lt = reduce(application, l.terms)?;
        let rt = reduce(application, r.terms)?;
        let diff = difference(application, &lt, &rt)?;
        Ok(json!({"outcome": if diff.is_empty() {"pass"} else {"fail"},
            "nontrivial": nontrivial, "lhs_lowered": lhs_keys.len(), "rhs_lowered": rhs_keys.len(),
            "lhs_lowered_keys": lhs_keys.iter().take(64).collect::<Vec<_>>(),
            "lhs_terminals": keys(&lt), "rhs_terminals": keys(&rt), "difference": diff,
            "difference_full": true}))
    };
    let mut value = run().unwrap_or_else(|e| json!({"outcome": "error", "error": e}));
    value["seconds"] = json!(started.elapsed().as_secs_f64());
    value
}

fn vk_case(vakint: &Vakint, settings: &VakintSettings, lhs: &Value, rhs: &Value) -> Value {
    let started = Instant::now();
    let run = || -> Result<Value, String> {
        let li = vk_input(&side(lhs)?)?;
        let ri = vk_input(&side(rhs)?)?;
        let lv = vakint
            .evaluate_rustred_symbolic(settings, li.as_view())
            .map_err(|e| format!("lhs: {e:?}"))?
            .value;
        let rv = vakint
            .evaluate_rustred_symbolic(settings, ri.as_view())
            .map_err(|e| format!("rhs: {e:?}"))?
            .value;
        let zero = (lv.clone() - &rv).expand().together().is_zero();
        Ok(json!({"outcome": if zero {"pass"} else {"fail"},
            "identical_strings": lv.to_canonical_string() == rv.to_canonical_string(),
            "lhs": lv.to_canonical_string().chars().take(300).collect::<String>(),
            "rhs": rv.to_canonical_string().chars().take(300).collect::<String>()}))
    };
    let mut value = run().unwrap_or_else(|e| json!({"outcome": "error", "error": e}));
    value["seconds"] = json!(started.elapsed().as_secs_f64());
    value
}

#[test]
#[ignore = "V7 pinch consistency: requires VAKINT_V7_ROOT, VAKINT_V7_CASES, VAKINT_V7_REPORT"]
fn staged_pack_pinch_consistency_v7() {
    let started = Instant::now();
    let root = PathBuf::from(std::env::var_os("VAKINT_V7_ROOT").expect("VAKINT_V7_ROOT"));
    let cases_path = PathBuf::from(std::env::var_os("VAKINT_V7_CASES").expect("VAKINT_V7_CASES"));
    let report = PathBuf::from(std::env::var_os("VAKINT_V7_REPORT").expect("VAKINT_V7_REPORT"));
    let budget: f64 = std::env::var("VAKINT_V7_BUDGET_SECONDS")
        .ok()
        .map_or(1500.0, |v| v.parse().unwrap());
    let vk_enabled = std::env::var("VAKINT_V7_VK").map_or(true, |v| v != "0");
    let only: Option<BTreeSet<String>> = std::env::var("VAKINT_V7_ONLY")
        .ok()
        .map(|v| v.split(',').map(str::to_owned).collect());
    let cases: Value = serde_json::from_slice(&std::fs::read(&cases_path).unwrap()).unwrap();
    let done: BTreeSet<String> = std::fs::read_to_string(&report)
        .unwrap_or_default()
        .lines()
        .filter_map(|l| serde_json::from_str::<Value>(l).ok())
        .filter_map(|v| v["id"].as_str().map(str::to_owned))
        .collect();

    let mut vakint = Vakint::new().unwrap();
    vakint.load_external_data(&root, 5).unwrap();
    let settings = settings();
    let application = vakint.saved_application(&settings).unwrap();
    let remaining = application.prepared.lock().unwrap().remaining_terminals().len();
    let load_seconds = started.elapsed().as_secs_f64();
    println!("{}", json!({"event": "loaded", "root": root, "remaining": remaining, "seconds": load_seconds}));

    let mut out = std::fs::OpenOptions::new()
        .create(true)
        .append(true)
        .open(&report)
        .unwrap();
    let case_timeout: u64 = std::env::var("VAKINT_V7_CASE_TIMEOUT")
        .ok()
        .map_or(0, |v| v.parse().unwrap());
    let current: Arc<Mutex<Option<(Value, Instant)>>> = Arc::new(Mutex::new(None));
    if case_timeout > 0 {
        let current = current.clone();
        let report = report.clone();
        let root = root.clone();
        std::thread::spawn(move || {
            loop {
                std::thread::sleep(std::time::Duration::from_secs(5));
                let guard = current.lock().unwrap();
                if let Some((case, since)) = guard.as_ref()
                    && since.elapsed().as_secs() > case_timeout
                {
                    let line = json!({"id": case["id"], "group": case["group"], "class": case["class"],
                        "variant": case["variant"], "seed": case["seed"], "root": root,
                        "rr": {"outcome": "timeout", "seconds": since.elapsed().as_secs_f64()}});
                    let mut f = std::fs::OpenOptions::new().append(true).open(&report).unwrap();
                    writeln!(f, "{line}").unwrap();
                    f.flush().unwrap();
                    println!("{}", json!({"event": "case_timeout", "id": case["id"]}));
                    std::process::exit(3);
                }
            }
        });
    }
    let (mut ran, mut skipped_budget) = (0usize, 0usize);
    for case in cases["cases"].as_array().unwrap() {
        let id = case["id"].as_str().unwrap().to_owned();
        if done.contains(&id) || only.as_ref().is_some_and(|o| !o.contains(&id)) {
            continue;
        }
        if started.elapsed().as_secs_f64() > budget {
            skipped_budget += 1;
            continue;
        }
        *current.lock().unwrap() = Some((case.clone(), Instant::now()));
        let rr = rr_case(&application, case);
        let mut line = json!({"id": id, "group": case["group"], "class": case["class"],
            "variant": case["variant"], "seed": case["seed"], "root": root,
            "remaining": remaining, "rr": rr});
        let admitted = line["rr"]["outcome"] != "not_admitted";
        if vk_enabled && admitted && case["vk"].as_bool().unwrap_or(false) {
            line["vk_sigma"] = vk_case(&vakint, &settings, &case["lhs"], &case["rhs"]);
            line["vk_identity"] = vk_case(&vakint, &settings, &case["lhs_identity"], &case["rhs"]);
        }
        let mut guard = current.lock().unwrap();
        *guard = None;
        writeln!(out, "{line}").unwrap();
        out.flush().unwrap();
        drop(guard);
        ran += 1;
    }
    println!("{}", json!({"event": "done", "ran": ran, "skipped_budget": skipped_budget,
        "seconds": started.elapsed().as_secs_f64()}));
}
