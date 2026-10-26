//! Bind completed public lowering receipts to native verified owner transports.
use eyre::{Result, ensure};
use rustred::{family::IntegralKey, solver::UnresolvedPolicy};
use rustred_app::{
    SavedApplicationEntryCoveragePolicy, SavedApplicationOptions, load_saved_application,
};
use serde_json::{Value, json};
use sha2::{Digest, Sha256};
use std::{
    collections::{BTreeMap, BTreeSet},
    fs::{self, File},
    io::{Read, Write},
    path::Path,
};
use symbolica::atom::AtomCore;

#[path = "build_cache/owner_proof.rs"]
#[cfg_attr(
    not(test),
    expect(
        dead_code,
        reason = "The mapper produces frontier certificates; consumers own legacy proof loading"
    )
)]
mod owner_proof;
#[path = "build_cache/workload.rs"]
#[cfg_attr(
    not(test),
    expect(
        dead_code,
        reason = "The mapper shares the existing workload types; selection validation is owned by the proof reader"
    )
)]
mod workload;
use hash as digest;
use owner_proof::{
    FrontierSelection, FrontierWork, validate_frontier_certificate, validate_hit_diagnostic,
};

fn check(ok: bool, why: &str) -> Result<()> {
    ensure!(ok, "{why}");
    Ok(())
}
fn hash(path: &Path) -> Result<String> {
    let mut file = File::open(path)?;
    let mut h = Sha256::new();
    let mut block = [0; 65536];
    loop {
        let n = file.read(&mut block)?;
        if n == 0 {
            break;
        }
        h.update(&block[..n]);
    }
    Ok(format!("{:x}", h.finalize()))
}
fn required<'a>(v: &'a Value, key: &str) -> Result<&'a str> {
    v[key]
        .as_str()
        .ok_or_else(|| eyre::eyre!("missing string {key}"))
}
fn required_hash<'a>(value: &'a Value, field: &str) -> Result<&'a str> {
    let hash = required(value, field)?;
    check(
        hash.len() == 64 && hash.bytes().all(|byte| byte.is_ascii_hexdigit()),
        "missing or invalid actual SHA256 binding",
    )?;
    Ok(hash)
}
fn compiled_binding(name: &str, value: Option<&'static str>) -> Result<&'static str> {
    value
        .filter(|value| !value.trim().is_empty())
        .ok_or_else(|| eyre::eyre!("mapping requires compile-time {name} provenance"))
}
fn available_executable(path: Option<&Path>, expected: &str) -> Result<()> {
    if let Some(path) = path.filter(|path| path.exists()) {
        check(
            hash(path)? == expected,
            "available measurement executable changed",
        )?;
    }
    Ok(())
}
fn compiled_cohort(input: &Value) -> Result<Value> {
    compiled_binding(
        "VAKINT_CACHE_BUILD_RUSTC",
        option_env!("VAKINT_CACHE_BUILD_RUSTC"),
    )?;
    compiled_binding(
        "VAKINT_CACHE_BUILD_PROFILE",
        option_env!("VAKINT_CACHE_BUILD_PROFILE"),
    )?;
    let path = Path::new(required(input, "cohort_receipt_path")?);
    check(
        hash(path)? == required_hash(input, "cohort_receipt_sha256")?,
        "compiled-cohort receipt changed",
    )?;
    let receipt = read(path)?;
    for (field, actual) in [
        (
            "main_revision",
            compiled_binding(
                "VAKINT_CACHE_BUILD_GIT_REVISION",
                option_env!("VAKINT_CACHE_BUILD_GIT_REVISION"),
            )?,
        ),
        (
            "native_revision",
            compiled_binding(
                "VAKINT_CACHE_BUILD_RUSTRED_REVISION",
                option_env!("VAKINT_CACHE_BUILD_RUSTRED_REVISION"),
            )?,
        ),
        (
            "kernel_revision",
            compiled_binding(
                "VAKINT_CACHE_BUILD_SYMBOLICA_REVISION",
                option_env!("VAKINT_CACHE_BUILD_SYMBOLICA_REVISION"),
            )?,
        ),
    ] {
        let revision = required(&receipt, field)?;
        check(
            revision.len() == 40
                && revision.bytes().all(|byte| byte.is_ascii_hexdigit())
                && revision == actual,
            "compiled revision is missing or differs from final receipt",
        )?;
    }
    for (field, actual) in [
        (
            "main_source_snapshot_sha256",
            compiled_binding(
                "VAKINT_CACHE_BUILD_MAIN_SOURCE_SNAPSHOT_SHA256",
                option_env!("VAKINT_CACHE_BUILD_MAIN_SOURCE_SNAPSHOT_SHA256"),
            )?,
        ),
        (
            "native_source_snapshot_sha256",
            compiled_binding(
                "VAKINT_CACHE_BUILD_NATIVE_SOURCE_SNAPSHOT_SHA256",
                option_env!("VAKINT_CACHE_BUILD_NATIVE_SOURCE_SNAPSHOT_SHA256"),
            )?,
        ),
    ] {
        check(
            required_hash(&receipt, field)? == actual,
            "actual compiled source snapshot differs from final receipt",
        )?;
    }
    for field in [
        "measurement_executable_sha256",
        "data_manifest_sha256",
        "loop_manifest_sha256",
        "application_sha256",
    ] {
        check(
            required_hash(&receipt, field)? == required_hash(input, field)?,
            "compiled-cohort data/binary binding differs",
        )?;
    }
    check(
        required_hash(&receipt, "mapper_source_sha256")?
            == format!(
                "{:x}",
                Sha256::digest(include_bytes!("authenticate_owner_roots.rs"))
            )
            && required_hash(&receipt, "mapper_executable_sha256")?
                == hash(&std::env::current_exe()?)?,
        "actual mapper source/binary differs from compiled-cohort receipt",
    )?;
    check(
        receipt["semantic_binding"].is_object(),
        "missing native program binding",
    )?;
    Ok(receipt)
}
fn read(path: &Path) -> Result<Value> {
    Ok(serde_json::from_slice(&fs::read(path)?)?)
}
fn validate_data(input: &Value) -> Result<()> {
    let root = Path::new(required(input, "data_path")?);
    for (path, field) in [
        ("manifest.json", "data_manifest_sha256"),
        ("5L/manifest.json", "loop_manifest_sha256"),
        ("5L/ibps/application.json", "application_sha256"),
    ] {
        check(
            hash(&root.join(path))? == required(input, field)?,
            "data changed or differs from owner-map input",
        )?;
    }
    Ok(())
}

fn certify_dispatch(
    prepared: &rustred_app::PreparedSavedApplication<15>,
    selection: &FrontierSelection,
    dispatch: &Value,
) -> Result<Value> {
    fn visit(
        source: FrontierWork,
        prepared: &rustred_app::PreparedSavedApplication<15>,
        selection: &FrontierSelection,
        active: &mut BTreeSet<FrontierWork>,
        steps: &mut BTreeMap<FrontierWork, Value>,
        frontier: &mut BTreeSet<FrontierWork>,
    ) -> Result<()> {
        frontier.insert(source.clone());
        if steps.contains_key(&source) {
            return Ok(());
        }
        check(active.insert(source.clone()), "active routing cycle")?;
        let mapped =
            prepared.route_entry_to_owner(&IntegralKey::try_new(source.powers.clone())?)?;
        let mut edges = Vec::new();
        for (endpoint, coefficient) in mapped.terms() {
            if coefficient.is_zero() {
                continue;
            }
            let target = selection.endpoint(&source, endpoint.powers().to_vec())?;
            edges.push(json!({"target":target.document(),"coefficient":coefficient.to_expression().to_canonical_string()}));
            frontier.insert(target.clone());
            if target.stage == "route" {
                visit(target, prepared, selection, active, steps, frontier)?;
            }
        }
        let guards = mapped
            .nonzero_conditions()
            .iter()
            .map(|guard| guard.to_expression().to_canonical_string())
            .collect();
        let step = selection.step(&source, edges, guards)?;
        active.remove(&source);
        steps.insert(source, step);
        Ok(())
    }
    let mut steps = BTreeMap::new();
    let mut frontier = BTreeSet::new();
    for seed in dispatch["seeds"]
        .as_array()
        .ok_or_else(|| eyre::eyre!("missing collected seeds"))?
    {
        visit(
            FrontierWork::parse(seed)?,
            prepared,
            selection,
            &mut BTreeSet::new(),
            &mut steps,
            &mut frontier,
        )?;
    }
    Ok(
        json!({"ordinal":dispatch["ordinal"],"family_fingerprint":dispatch["family_fingerprint"],"seeds":dispatch["seeds"],
        "seed_sha256":format!("{:x}",Sha256::digest(serde_json::to_vec(&dispatch["seeds"])?)),
        "steps":steps.into_values().collect::<Vec<_>>(),"frontier":frontier.iter().map(FrontierWork::document).collect::<Vec<_>>()}),
    )
}

/// Join one already completed cached diagnostic, never repeat its reduction.
/// The legacy uncached receipt and its conservative owner sets remain intact.
fn initial_frontier(
    case: &Value,
    reference: &Value,
    input: &Value,
    binding: &Value,
    compiled: &Value,
    prepared: &rustred_app::PreparedSavedApplication<15>,
    bound_files: &mut Vec<(std::path::PathBuf, String)>,
) -> Result<Value> {
    let Some(path) = case["cache_hit_receipt_path"].as_str() else {
        check(
            case["cache_hit_receipt_sha256"].is_null(),
            "cache-hit receipt hash has no path",
        )?;
        return Ok(Value::Null);
    };
    let path = Path::new(path);
    let receipt_hash = required_hash(case, "cache_hit_receipt_sha256")?;
    check(
        hash(path)? == receipt_hash,
        "saved diagnostic receipt changed",
    )?;
    let receipt = read(path)?;
    owner_proof::validate_frontier_outputs(reference, &receipt)?;
    check(
        receipt["success"] == true
            && receipt["worker_success"] == true
            && receipt["fresh_process"] == true
            && receipt["supervision"]["success"] == true
            && receipt["supervision"]["exit_code"] == 0
            && receipt["supervision"]["deadline_exhausted"] == false
            && receipt["case_id"] == reference["case_id"]
            && receipt["input_sha256"] == reference["input_sha256"]
            && receipt["identity"]["workload_sha256"] == reference["identity"]["workload_sha256"]
            && receipt["identity"]["source_census_sha256"]
                == reference["identity"]["source_census_sha256"]
            && receipt["identity"]["data_manifest_sha256"] == input["data_manifest_sha256"]
            && receipt["identity"]["executable_sha256"] == input["measurement_executable_sha256"]
            && receipt["cache_statistics"]["binding"] == *binding,
        "saved diagnostic differs from complete case/cohort/application",
    )?;
    for (file, field) in [
        ("value.symbolica", "value_sha256"),
        ("guards.json", "guards_sha256"),
    ] {
        let output = path.with_file_name(file);
        check(
            hash(&output)? == required_hash(&receipt, field)? && receipt[field] == reference[field],
            "saved diagnostic output or canonical guard-set file differs",
        )?;
        bound_files.push((output, required_hash(&receipt, field)?.to_owned()));
    }
    let trace_ref = &receipt["cache_hit_diagnostic"];
    let trace_path = Path::new(required(trace_ref, "path")?);
    let trace_path = if trace_path.is_absolute() {
        trace_path.to_owned()
    } else {
        path.parent().unwrap_or(Path::new(".")).join(trace_path)
    };
    let trace_hash = required_hash(trace_ref, "sha256")?;
    check(hash(&trace_path)? == trace_hash, "saved hit stream changed")?;
    bound_files.push((path.to_owned(), receipt_hash.to_owned()));
    bound_files.push((trace_path.clone(), trace_hash.to_owned()));
    let mut certificate = json!({"schema":"vakint.initial-routing-frontier-hit-proof.v1","complete":false,"reason":null,
        "trace_path":trace_path,"trace_sha256":trace_hash,"cache_hit_receipt_path":path,"cache_hit_receipt_sha256":receipt_hash,
        "semantic_binding":binding,"compiled_cohort":compiled,"input_sha256":receipt["input_sha256"],
        "selection_sha256":format!("{:x}",Sha256::digest(prepared.selection_json().as_bytes())),
        "selection":null,"dispatches":[],"derived":null});
    let build = (|| -> Result<()> {
        let normalized = validate_hit_diagnostic(&read(&trace_path)?)?;
        for field in [
            "lookups",
            "route_hits",
            "apply_hits",
            "requested_root_hits",
            "intermediate_hits",
            "misses",
        ] {
            check(
                normalized["persistent"][field] == receipt["cache_statistics"]["persistent"][field],
                "diagnostic dispatch counters differ from final raw receipt",
            )?;
        }
        let selection =
            FrontierSelection::parse(&serde_json::from_str(prepared.selection_json())?, binding)?;
        certificate["selection"] = selection.document();
        let dispatches = normalized["dispatches"]
            .as_array()
            .ok_or_else(|| eyre::eyre!("missing normalized dispatches"))?;
        certificate["dispatches"] = json!(
            dispatches
                .iter()
                .map(|dispatch| certify_dispatch(prepared, &selection, dispatch))
                .collect::<Result<Vec<_>>>()?
        );
        certificate["complete"] = json!(true);
        certificate["derived"] = validate_frontier_certificate(
            &certificate,
            &normalized,
            binding,
            compiled,
            trace_hash,
        )?;
        Ok(())
    })();
    if let Err(error) = build {
        certificate["complete"] = json!(false);
        certificate["reason"] = json!(error.to_string());
        certificate["derived"] = json!({"complete":false,"root_frontier_hits":null,"proper_descendant_hits":null,"scope":"unclassified"});
    }
    Ok(certificate)
}
fn main() -> Result<()> {
    let args = std::env::args().collect::<Vec<_>>();
    check(
        args.len() == 3,
        "usage: authenticate_owner_roots INPUT.json OUTPUT.json",
    )?;
    let input_path = Path::new(&args[1]);
    let output = Path::new(&args[2]);
    check(!output.exists(), "preserve the existing proof output")?;
    let input_hash = hash(input_path)?;
    let input = read(input_path)?;
    check(
        input["schema"] == "vakint.public-lowering-owner-map-input.v1",
        "wrong mapping-input schema",
    )?;
    let cohort = compiled_cohort(&input)?;
    let cohort_path = Path::new(required(&input, "cohort_receipt_path")?);
    let measurement_executable = input["measurement_executable_path"]
        .as_str()
        .filter(|path| !path.trim().is_empty())
        .map(Path::new);
    available_executable(
        measurement_executable,
        required_hash(&input, "measurement_executable_sha256")?,
    )?;
    validate_data(&input)?;
    for (path, expected) in [
        ("frozen_suite_path", "frozen_suite_sha256"),
        ("reserve_path", "reserve_sha256"),
        ("protocol_path", "protocol_sha256"),
    ] {
        check(
            hash(Path::new(required(&input, path)?))? == required(&input, expected)?,
            "frozen identity or selection protocol changed",
        )?;
    }
    let frozen_suite = read(Path::new(required(&input, "frozen_suite_path")?))?;
    let frozen_reserves = read(Path::new(required(&input, "reserve_path")?))?;

    let options = SavedApplicationOptions {
        unresolved_policy: UnresolvedPolicy::Preserve,
        entry_coverage_policy: SavedApplicationEntryCoveragePolicy::Advisory,
        preparation_workers: 1,
        ..Default::default()
    };
    let prepared = load_saved_application::<15>(
        &Path::new(required(&input, "data_path")?).join("5L/ibps"),
        options,
    )?;
    let mut rows = Vec::new();
    let mut ids = BTreeSet::new();
    let mut observed_binding: Option<Value> = None;
    let mut bound_files = Vec::new();
    for case in input["cases"]
        .as_array()
        .ok_or_else(|| eyre::eyre!("missing case index"))?
    {
        let id = required(case, "id")?;
        check(ids.insert(id.to_owned()), "duplicate indexed case")?;
        let path = Path::new(required(case, "receipt_path")?);
        let receipt_hash = hash(path)?;
        check(
            receipt_hash == required(case, "receipt_sha256")?,
            "measurement receipt changed",
        )?;
        let receipt = read(path)?;
        required_hash(&receipt, "input_sha256")?;
        bound_files.push((path.to_owned(), receipt_hash.clone()));
        check(
            receipt["case_id"] == id && receipt["variant"] == "uncached" && receipt["sample"] == 0,
            "wrong case/variant/sample",
        )?;
        check(
            receipt["fresh_process"] == true
                && receipt["identity"]["cache_sha256"].is_null()
                && receipt["success"] == true
                && receipt["worker_success"] == true
                && receipt["supervision"]["success"] == true
                && receipt["supervision"]["exit_code"] == 0
                && receipt["supervision"]["deadline_exhausted"] == false,
            "only complete successful uncached fresh-process reports authenticate roots",
        )?;
        check(
            required_hash(&receipt["identity"], "data_manifest_sha256")?
                == required_hash(&input, "data_manifest_sha256")?
                && required_hash(&receipt["identity"], "executable_sha256")?
                    == required_hash(&input, "measurement_executable_sha256")?,
            "wrong data or public-evaluator binary",
        )?;
        let workload_path = Path::new(required(case, "workload_path")?);
        check(
            hash(workload_path)? == required(case, "workload_sha256")?
                && receipt["identity"]["workload_sha256"] == case["workload_sha256"],
            "workload provenance mismatch",
        )?;
        let workload = read(workload_path)?;
        bound_files.push((
            workload_path.to_owned(),
            required(case, "workload_sha256")?.to_owned(),
        ));
        let declarations = workload["tuning"]
            .as_array()
            .ok_or_else(|| eyre::eyre!("missing tuning"))?
            .iter()
            .chain(
                workload["held_out"]
                    .as_array()
                    .ok_or_else(|| eyre::eyre!("missing held_out"))?,
            )
            .filter(|r| r["id"] == id)
            .collect::<Vec<_>>();
        check(
            declarations.len() == 1,
            "case declaration is ambiguous or absent",
        )?;
        let declared = declarations[0];
        let frozen_cases = match required(case, "logical_cohort")? {
            "tuning" => frozen_suite["tuning"]
                .as_array()
                .ok_or_else(|| eyre::eyre!("missing frozen tuning cases"))?
                .iter()
                .collect::<Vec<_>>(),
            "held-out" => frozen_suite["held_out"]
                .as_array()
                .ok_or_else(|| eyre::eyre!("missing frozen held-out cases"))?
                .iter()
                .chain(
                    frozen_reserves["held_out_reserves"]
                        .as_array()
                        .ok_or_else(|| eyre::eyre!("missing frozen reserve cases"))?,
                )
                .collect(),
            _ => eyre::bail!("unknown frozen case role"),
        };
        let frozen = frozen_cases
            .into_iter()
            .filter(|frozen| frozen["id"] == id)
            .collect::<Vec<_>>();
        check(
            frozen.len() == 1
                && frozen[0]["public_input"] == declared["public_input"]
                && frozen[0]["powers"] == declared["powers"],
            "public input differs from its frozen identity or role",
        )?;
        check(
            format!(
                "{:x}",
                Sha256::digest(required(declared, "public_input")?.as_bytes())
            ) == required(case, "public_input_source_sha256")?,
            "frozen public input changed",
        )?;
        let source_census = Path::new(required(case, "source_census_path")?);
        check(
            hash(source_census)? == required(&workload, "source_census_sha256")?
                && receipt["identity"]["source_census_sha256"] == workload["source_census_sha256"],
            "source census changed or belongs to another attempt",
        )?;
        bound_files.push((
            source_census.to_owned(),
            required(&workload, "source_census_sha256")?.to_owned(),
        ));
        check(
            hash(&path.with_file_name("value.symbolica"))? == required(&receipt, "value_sha256")?
                && hash(&path.with_file_name("guards.json"))?
                    == required(&receipt, "guards_sha256")?,
            "successful output or guard receipt changed",
        )?;
        bound_files.push((
            path.with_file_name("value.symbolica"),
            required(&receipt, "value_sha256")?.to_owned(),
        ));
        bound_files.push((
            path.with_file_name("guards.json"),
            required(&receipt, "guards_sha256")?.to_owned(),
        ));
        let binding = receipt["cache_statistics"]["binding"].clone();
        let digest = serde_json::from_value::<Vec<u8>>(binding["artifact_digest"].clone())?;
        let digest: [u8; 32] = digest
            .try_into()
            .map_err(|_| eyre::eyre!("invalid artifact digest"))?;
        let native = prepared.cache_binding(digest);
        let native_binding = json!({"artifact_digest":native.artifact_digest,"family_fingerprint":native.family_fingerprint,"physical_arity":native.physical_arity,"layer":native.layer});
        check(
            binding == native_binding && binding == cohort["semantic_binding"],
            "native mapper differs from the public semantic binding",
        )?;
        if let Some(previous) = &observed_binding {
            check(
                previous == &binding,
                "public cases used different semantic applications",
            )?;
        } else {
            observed_binding = Some(binding.clone());
        }
        let keys =
            serde_json::from_value::<Vec<Vec<i64>>>(receipt["actual_lowered_roots"].clone())?;
        check(
            !keys.is_empty() && keys.iter().all(|key| key.len() == 15),
            "missing, empty or malformed actual five-loop roots",
        )?;
        let mut protected = BTreeSet::new();
        let mut transports = Vec::new();
        let mut mapping_error = None;
        for powers in &keys {
            let key = IntegralKey::try_new(powers.clone())?;
            match prepared.route_entry_to_owner(&key) {
                Ok(mapped) => {
                    let mut endpoints = Vec::new();
                    for (owner, coefficient) in mapped.terms() {
                        if coefficient.is_zero() {
                            continue;
                        }
                        protected.insert(owner.powers().to_vec());
                        endpoints.push(json!({"powers":owner.powers(),"coefficient":coefficient.to_expression().to_canonical_string(),"recorded_coverage_admitted":prepared.entry_is_admitted(owner)?}));
                    }
                    transports.push(json!({"powers":powers,"owner_endpoints":endpoints,"nonzero_conditions":mapped.nonzero_conditions().iter().map(|guard|guard.to_expression().to_canonical_string()).collect::<Vec<_>>() }));
                }
                Err(error) => {
                    mapping_error = Some(error.to_string());
                    break;
                }
            }
        }
        let protected = protected.into_iter().collect::<Vec<_>>();
        let signature = format!("{:x}", Sha256::digest(serde_json::to_vec(&protected)?));
        let mut frontier = match initial_frontier(
            case,
            &receipt,
            &input,
            &binding,
            &cohort,
            &prepared,
            &mut bound_files,
        ) {
            Ok(frontier) => frontier,
            Err(error) => json!({"schema":"vakint.initial-routing-frontier-hit-proof.v1",
                "complete":false,"reason":error.to_string(),"semantic_binding":binding,
                "compiled_cohort":cohort,"derived":{"complete":false,"scope":"unclassified",
                    "root_frontier_hits":null,"proper_descendant_hits":null}}),
        };
        if frontier["complete"] == true {
            match owner_proof::declared_training_matches(&frontier, &workload) {
                Ok(matches) => frontier["declared_training_matches"] = matches,
                Err(error) => {
                    frontier["complete"] = json!(false);
                    frontier["reason"] = json!(error.to_string());
                    frontier["derived"] = json!({"complete":false,"scope":"unclassified",
                        "root_frontier_hits":null,"proper_descendant_hits":null});
                }
            }
        }
        let mut row = json!({"id":id,"logical_cohort":case["logical_cohort"],"receipt_path":path,"receipt_sha256":receipt_hash,"workload_path":workload_path,"workload_sha256":case["workload_sha256"],"source_census_path":source_census,"source_census_sha256":workload["source_census_sha256"],"application_binding":binding,"input_sha256":receipt["input_sha256"],"public_input_source_sha256":case["public_input_source_sha256"],"actual_lowered_roots":keys,"protected_owner_roots":protected,"protected_root_set_sha256":signature,"complete_owner_transport":mapping_error.is_none(),"mapping_error":mapping_error,"native_transports":transports,"reduction_succeeded":true,"measurement_claim_allowed":false});
        if !frontier.is_null() {
            row["initial_routing_frontier"] = frontier;
        }
        rows.push(row);
    }
    validate_data(&input)?;

    for (path, expected) in [
        ("frozen_suite_path", "frozen_suite_sha256"),
        ("reserve_path", "reserve_sha256"),
        ("protocol_path", "protocol_sha256"),
    ] {
        check(
            hash(Path::new(required(&input, path)?))? == required(&input, expected)?,
            "frozen identity or selection protocol changed during mapping",
        )?;
    }
    for (path, expected_sha256) in &bound_files {
        check(
            hash(path)? == *expected_sha256,
            "owner-mapping evidence changed during execution",
        )?;
    }
    check(
        hash(input_path)? == input_hash,
        "owner-map index changed during execution",
    )?;
    available_executable(
        measurement_executable,
        required_hash(&input, "measurement_executable_sha256")?,
    )?;
    check(
        hash(cohort_path)? == required_hash(&input, "cohort_receipt_sha256")?,
        "compiled-cohort receipt changed during mapping",
    )?;
    check(
        observed_binding.as_ref() == Some(&cohort["semantic_binding"]),
        "no complete current-program case was authenticated",
    )?;
    let value = json!({"schema":"vakint.public-lowering-owner-map-proof.v1","authority":"Native verified owner transport of every observed key in complete successful public evaluation reports; no additional IBP or campaign-closure claim. Protected union is conservative before possible owner-combination cancellation.","input_path":input_path,"input_sha256":input_hash,"measurement_executable_path":measurement_executable,"measurement_executable_sha256":input["measurement_executable_sha256"],"mapper_executable_path":std::env::current_exe()?,"mapper_executable_sha256":hash(&std::env::current_exe()?)?,"mapper_source_sha256":format!("{:x}",Sha256::digest(include_bytes!("authenticate_owner_roots.rs"))),"mapper_rustc":compiled_binding("VAKINT_CACHE_BUILD_RUSTC", option_env!("VAKINT_CACHE_BUILD_RUSTC"))?,"mapper_profile":compiled_binding("VAKINT_CACHE_BUILD_PROFILE", option_env!("VAKINT_CACHE_BUILD_PROFILE"))?,"mapper_git_revision":cohort["main_revision"],"kernel_revision":cohort["kernel_revision"],"rustred_revision":cohort["native_revision"],"compiled_cohort":cohort,"cohort_receipt_path":cohort_path,"cohort_receipt_sha256":input["cohort_receipt_sha256"],"native_options":format!("{options:?}"),"semantic_binding":observed_binding,"data_manifest_sha256":input["data_manifest_sha256"],"loop_manifest_sha256":input["loop_manifest_sha256"],"application_sha256":input["application_sha256"],"frozen_suite_path":input["frozen_suite_path"],"frozen_suite_sha256":input["frozen_suite_sha256"],"reserve_path":input["reserve_path"],"reserve_sha256":input["reserve_sha256"],"protocol_path":input["protocol_path"],"protocol_sha256":input["protocol_sha256"],"rows":rows});
    let mut file = fs::OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(output)?;
    file.write_all(&serde_json::to_vec_pretty(&value)?)?;
    file.sync_all()?;
    File::open(output.parent().unwrap_or(Path::new(".")))?.sync_all()?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::compiled_binding;

    #[test]
    fn mapping_requires_nonempty_actual_compile_time_bindings() {
        for value in [None, Some(""), Some(" \n\t")] {
            let error = compiled_binding("TEST_REVISION", value).unwrap_err();
            assert!(error.to_string().contains("compile-time TEST_REVISION"));
        }
        assert_eq!(
            compiled_binding("TEST_REVISION", Some("abc")).unwrap(),
            "abc"
        );
    }
}
