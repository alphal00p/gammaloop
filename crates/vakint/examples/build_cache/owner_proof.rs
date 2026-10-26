//! Pinned native owner maps from complete successful public evaluations only.
use super::{
    digest,
    workload::{Root, Workload},
};
use eyre::{Result, ensure};
use serde_json::{Value, json};
use sha2::{Digest, Sha256};
use std::{
    collections::{BTreeMap, BTreeSet},
    fs,
    path::{Path, PathBuf},
};

pub(super) struct OwnerProof {
    document: Value,
    identity: Value,
    correlations: Value,
}

fn string<'a>(value: &'a Value, field: &str) -> Result<&'a str> {
    value[field]
        .as_str()
        .ok_or_else(|| eyre::eyre!("missing owner-proof {field}"))
}
fn sha256<'a>(value: &'a Value, field: &str) -> Result<&'a str> {
    let hash = string(value, field)?;
    ensure!(
        hash.len() == 64 && hash.bytes().all(|byte| byte.is_ascii_hexdigit()),
        "invalid or missing owner-proof SHA256 {field}"
    );
    Ok(hash)
}
fn available_executable(document: &Value, base: &Path, path: &str, hash: &str) -> Result<()> {
    let expected = sha256(document, hash)?;
    if let Some(path) = document[path]
        .as_str()
        .filter(|path| !path.trim().is_empty())
    {
        let path = resolved(base, path);
        if path.exists() {
            ensure!(
                digest(&path)? == expected,
                "available executable differs from provenance"
            );
        }
    }
    Ok(())
}
fn compiled_bindings() -> [(&'static str, Option<&'static str>); 5] {
    [
        (
            "main_revision",
            option_env!("VAKINT_CACHE_BUILD_GIT_REVISION"),
        ),
        (
            "native_revision",
            option_env!("VAKINT_CACHE_BUILD_RUSTRED_REVISION"),
        ),
        (
            "kernel_revision",
            option_env!("VAKINT_CACHE_BUILD_SYMBOLICA_REVISION"),
        ),
        (
            "main_source_snapshot_sha256",
            option_env!("VAKINT_CACHE_BUILD_MAIN_SOURCE_SNAPSHOT_SHA256"),
        ),
        (
            "native_source_snapshot_sha256",
            option_env!("VAKINT_CACHE_BUILD_NATIVE_SOURCE_SNAPSHOT_SHA256"),
        ),
    ]
}
fn cohort(document: &Value, binding: &Value, expected: &[(&str, Option<&str>); 5]) -> Result<()> {
    let compiled = &document["compiled_cohort"];
    ensure!(
        compiled.is_object(),
        "missing final compiled-cohort receipt"
    );
    for field in ["main_revision", "native_revision", "kernel_revision"] {
        let revision = string(compiled, field)?;
        ensure!(
            revision.len() == 40 && revision.bytes().all(|byte| byte.is_ascii_hexdigit()),
            "missing actual compiled revision {field}"
        );
    }
    for field in [
        "main_source_snapshot_sha256",
        "native_source_snapshot_sha256",
    ] {
        sha256(compiled, field)?;
    }
    for field in [
        "mapper_source_sha256",
        "mapper_executable_sha256",
        "measurement_executable_sha256",
        "data_manifest_sha256",
        "loop_manifest_sha256",
        "application_sha256",
    ] {
        ensure!(
            sha256(compiled, field)? == sha256(document, field)?,
            "compiled-cohort binding differs: {field}"
        );
    }
    sha256(document, "cohort_receipt_sha256")?;
    ensure!(
        binding.is_object()
            && binding["artifact_digest"]
                .as_array()
                .is_some_and(|bytes| bytes.len() == 32
                    && bytes
                        .iter()
                        .all(|byte| byte.as_u64().is_some_and(|byte| byte <= 255)))
            && !string(binding, "family_fingerprint")?.trim().is_empty()
            && compiled["semantic_binding"] == *binding,
        "missing or mismatched native program/family binding"
    );
    for &(field, actual) in expected {
        ensure!(
            actual.is_some_and(
                |actual| !actual.trim().is_empty() && compiled[field].as_str() == Some(actual)
            ),
            "owner proof differs from actual compiled consumer: {field}"
        );
    }
    Ok(())
}
pub(super) fn validate_completed_outputs(
    receipt: &Value,
    directory: &Path,
    root: &Root,
) -> Result<()> {
    if receipt["success"] != true {
        // Failed attempts remain excluded without retry; they are never adopted.
        return Ok(());
    }
    ensure!(
        receipt["worker_success"] == true
            && receipt["fresh_process"] == true
            && receipt["supervision"]["success"] == true
            && receipt["supervision"]["exit_code"] == 0
            && receipt["supervision"]["deadline_exhausted"] == false,
        "incomplete successful worker receipt"
    );
    for (name, field) in [
        ("value.symbolica", "value_sha256"),
        ("guards.json", "guards_sha256"),
    ] {
        let expected = receipt[field]
            .as_str()
            .filter(|hash| hash.len() == 64 && hash.bytes().all(|byte| byte.is_ascii_hexdigit()))
            .ok_or_else(|| eyre::eyre!("missing completed-output SHA256 {field}"))?;
        ensure!(
            digest(&directory.join(name))? == expected,
            "completed output changed: {name}"
        );
    }
    sha256(receipt, "input_sha256")?;
    if receipt["identity"]["owner_proof"].is_object() {
        let pinned = &receipt["identity"]["owner_proof"];
        ensure!(
            receipt["native_keys_authenticated"] == true
                && receipt["actual_lowered_roots"] == json!(root.lowered_roots)
                && receipt["owner_proof"] == *pinned
                && receipt["cache_statistics"]["binding"]
                    == pinned["compiled_cohort"]["semantic_binding"],
            "resumed public key/program evidence differs from the pinned case"
        );
    }
    Ok(())
}

fn array<'a>(value: &'a Value, field: &str) -> Result<&'a Vec<Value>> {
    value[field]
        .as_array()
        .ok_or_else(|| eyre::eyre!("missing owner-proof {field}"))
}
fn roots(value: &Value) -> Result<BTreeSet<Vec<i64>>> {
    let keys: Vec<Vec<i64>> = serde_json::from_value(value.clone())?;
    let set = keys.iter().cloned().collect::<BTreeSet<_>>();
    ensure!(
        keys.iter().all(|key| key.len() == 15) && keys == set.iter().cloned().collect::<Vec<_>>(),
        "owner roots must be sorted unique fifteen-axis keys"
    );
    Ok(set)
}
fn resolved(base: &Path, path: &str) -> PathBuf {
    let path = Path::new(path);
    if path.is_absolute() {
        path.into()
    } else {
        base.join(path)
    }
}
fn evidence(document: &Value, base: &Path, path: &str, hash: &str) -> Result<PathBuf> {
    let path = resolved(base, string(document, path)?);
    ensure!(
        digest(&path)? == sha256(document, hash)?,
        "owner-proof evidence changed: {}",
        path.display()
    );
    Ok(path)
}

/// Full native work identity. A Route and an Apply of the same powers are
/// different cache records; an Apply also carries its literal owner support.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub(super) struct FrontierWork {
    pub(super) stage: String,
    pub(super) powers: Vec<i64>,
    owner: Option<Vec<bool>>,
}

impl FrontierWork {
    pub(super) fn route(powers: Vec<i64>) -> Result<Self> {
        Self::parse(&json!({"stage":"route","powers":powers,"owner_sector":null}))
    }

    pub(super) fn parse(value: &Value) -> Result<Self> {
        let stage = string(value, "stage")?.to_owned();
        let powers: Vec<i64> = serde_json::from_value(value["powers"].clone())?;
        ensure!(powers.len() == 15, "work identity must have fifteen axes");
        ensure!(value.get("owner_sector").is_some(), "omitted owner support");
        let owner = match stage.as_str() {
            "route" => {
                ensure!(
                    value["owner_sector"].is_null(),
                    "Route carries an Apply owner"
                );
                None
            }
            "apply" => {
                let owner: Vec<bool> = serde_json::from_value(value["owner_sector"].clone())?;
                ensure!(
                    owner.len() == 15 && owner == Self::support(&powers),
                    "Apply owner differs from literal support"
                );
                Some(owner)
            }
            _ => eyre::bail!("unknown cache work stage"),
        };
        Ok(Self {
            stage,
            powers,
            owner,
        })
    }

    fn support(powers: &[i64]) -> Vec<bool> {
        powers.iter().map(|power| *power > 0).collect()
    }

    pub(super) fn document(&self) -> Value {
        json!({"stage":self.stage,"powers":self.powers,"owner_sector":self.owner})
    }
}

fn count(value: &Value, field: &str) -> Result<u64> {
    value[field]
        .as_u64()
        .ok_or_else(|| eyre::eyre!("missing nonnegative counter {field}"))
}

/// A successful retained result can still contain unknown integrals. This
/// optional capability requires actual closed outputs on both sides of its join.
pub(super) fn validate_frontier_outputs(original: &Value, cached: &Value) -> Result<()> {
    for receipt in [original, cached] {
        ensure!(
            receipt["retained_result_available"] == true
                && receipt["catalogue_closure"]["closed"] == true
                && receipt["catalogue_closure"]["unknown_count"] == 0
                && receipt["catalogue_closure"]["unknown_masters"]
                    .as_array()
                    .is_some_and(Vec::is_empty),
            "frontier capability requires complete catalogue-closed retained outputs"
        );
    }
    Ok(())
}

/// Reconcile a lossless diagnostic with native per-dispatch counters. A latest
/// monitor snapshot or a root union is deliberately insufficient evidence.
pub(super) fn validate_hit_diagnostic(document: &Value) -> Result<Value> {
    ensure!(
        document["schema"] == "vakint.cache-hit-diagnostic.v1"
            && document["interval_ns"] == 0
            && document["success"] == true
            && document.get("error").is_some_and(Value::is_null),
        "incomplete cache-hit diagnostic"
    );
    let expected = count(document, "input_dispatch_count")?;
    let mut dispatches = Vec::new();
    let mut current: Option<Value> = None;
    let mut roots = BTreeSet::new();
    let mut hits = Vec::new();
    let mut last_sequence = None;
    let mut finished = false;
    for event in array(document, "events")? {
        ensure!(
            !finished && event["schema"] == "vakint.cache-hit-events.v1",
            "event after observer finish or wrong event schema"
        );
        let kind = string(event, "event")?;
        if kind == "observer_finish" {
            ensure!(
                current.is_none()
                    && event["observer_status"] == "Active"
                    && event["operation_succeeded"] == true
                    && event["interval_ns"] == 0
                    && event.get("sequence").is_some_and(Value::is_null),
                "observer detached, faulted or operation incomplete"
            );
            finished = true;
            continue;
        }
        let sequence = count(event, "sequence")?;
        ensure!(
            last_sequence.is_none_or(|previous| sequence > previous),
            "native event sequence is not increasing"
        );
        last_sequence = Some(sequence);
        let family = string(event, "family_fingerprint")?;
        ensure!(!family.trim().is_empty(), "omitted native family");
        for field in [
            "work",
            "rootwork",
            "requested_roots",
            "completed_roots",
            "persistent",
            "status",
        ] {
            ensure!(
                event.get(field).is_some(),
                "omitted native event field {field}"
            );
        }
        if !event["work"].is_null() {
            FrontierWork::parse(&event["work"])?;
        }
        ensure!(
            event["status"]
                == if kind == "dispatch_finish" {
                    "succeeded"
                } else {
                    "running"
                },
            "invalid native event status"
        );
        ensure!(
            kind == "root_collected" || event["rootwork"].is_null(),
            "unexpected root identity"
        );
        if kind == "dispatch_start" {
            ensure!(current.is_none(), "nested or missing dispatch finish");
            current = Some(json!({"ordinal":dispatches.len(),"family_fingerprint":family}));
            roots.clear();
            hits.clear();
            continue;
        }
        let active = current
            .as_ref()
            .ok_or_else(|| eyre::eyre!("event without a dispatch start"))?;
        ensure!(
            active["family_fingerprint"] == family,
            "dispatch family changed"
        );
        match kind {
            "backing_hit" => hits.push(FrontierWork::parse(&event["work"])?),
            "root_collected" => {
                let root = FrontierWork::parse(&event["rootwork"])?;
                let current = FrontierWork::parse(&event["work"])?;
                ensure!(
                    root.stage == "route" && current == root && roots.insert(root),
                    "omitted, duplicate or non-Route collected root"
                );
            }
            "dispatch_finish" => {
                ensure!(event["status"] == "succeeded", "native dispatch failed");
                let requested = count(event, "requested_roots")?;
                ensure!(
                    requested == count(event, "completed_roots")?
                        && requested == u64::try_from(roots.len())?,
                    "collected/requested/completed roots differ"
                );
                let persistent = &event["persistent"];
                count(persistent, "lookups")?;
                count(persistent, "misses")?;
                let route = count(persistent, "route_hits")?;
                let apply = count(persistent, "apply_hits")?;
                ensure!(
                    route == u64::try_from(hits.iter().filter(|hit| hit.stage == "route").count())?
                        && apply
                            == u64::try_from(
                                hits.iter().filter(|hit| hit.stage == "apply").count()
                            )?
                        && route.checked_add(apply).is_some_and(|total| count(
                            persistent,
                            "requested_root_hits"
                        )
                        .ok()
                        .and_then(
                            |roots| roots.checked_add(count(persistent, "intermediate_hits").ok()?)
                        ) == Some(total)),
                    "lost cache-hit event or inconsistent native counters"
                );
                let mut normalized = current.take().unwrap();
                normalized["seeds"] =
                    json!(roots.iter().map(FrontierWork::document).collect::<Vec<_>>());
                normalized["hits"] =
                    json!(hits.iter().map(FrontierWork::document).collect::<Vec<_>>());
                normalized["persistent"] = persistent.clone();
                dispatches.push(normalized);
            }
            _ => eyre::bail!("unknown cache-hit event"),
        }
    }
    ensure!(
        finished && current.is_none() && u64::try_from(dispatches.len())? == expected,
        "missing dispatch or observer completion"
    );
    let mut totals = serde_json::Map::new();
    for field in [
        "lookups",
        "route_hits",
        "apply_hits",
        "requested_root_hits",
        "intermediate_hits",
        "misses",
    ] {
        let total = dispatches.iter().try_fold(0_u64, |total, dispatch| {
            total
                .checked_add(count(&dispatch["persistent"], field)?)
                .ok_or_else(|| eyre::eyre!("persistent counter sum overflow"))
        })?;
        totals.insert(field.into(), json!(total));
    }
    Ok(
        json!({"interval_ns":0,"input_dispatch_count":expected,"dispatches":dispatches,"persistent":totals}),
    )
}

/// Only the labels needed to classify a native route result. Transport matrices
/// remain owned and verified by Native, never interpreted by this proof reader.
pub(super) struct FrontierSelection {
    family: String,
    owners: BTreeSet<Vec<bool>>,
    routes: BTreeMap<Vec<bool>, Vec<bool>>,
}

impl FrontierSelection {
    pub(super) fn parse(document: &Value, binding: &Value) -> Result<Self> {
        let family = string(document, "family_fingerprint")?.to_owned();
        ensure!(
            binding["physical_arity"] == 15 && binding["family_fingerprint"] == family,
            "selection family/arity differs from loaded application"
        );
        let mut owners = BTreeSet::new();
        for row in array(document, "owners")? {
            ensure!(
                owners.insert(Self::mask(string(row, "mask")?)?),
                "duplicate owner label"
            );
        }
        ensure!(!owners.is_empty(), "missing native owners");
        let mut routes = BTreeMap::new();
        for row in array(document, "initial_frontier_routes")? {
            let source = Self::mask(string(row, "source_mask")?)?;
            let owner = Self::mask(string(row, "owner_mask")?)?;
            ensure!(
                owners.contains(&owner) && routes.insert(source, owner).is_none(),
                "unknown route owner or duplicate route label"
            );
        }
        Ok(Self {
            family,
            owners,
            routes,
        })
    }

    fn mask(value: &str) -> Result<Vec<bool>> {
        ensure!(
            value.len() == 15 && value.bytes().all(|byte| byte == b'0' || byte == b'1'),
            "invalid fifteen-axis native mask"
        );
        Ok(value.bytes().map(|byte| byte == b'1').collect())
    }

    fn label(mask: &[bool]) -> String {
        mask.iter()
            .map(|active| if *active { '1' } else { '0' })
            .collect()
    }

    pub(super) fn owner(&self, source: &FrontierWork) -> Result<(Vec<bool>, bool)> {
        let support = FrontierWork::support(&source.powers);
        if self.owners.contains(&support) {
            return Ok((support, true));
        }
        Ok((
            self.routes
                .get(&support)
                .cloned()
                .ok_or_else(|| eyre::eyre!("missing native route label"))?,
            false,
        ))
    }

    pub(super) fn endpoint(&self, source: &FrontierWork, powers: Vec<i64>) -> Result<FrontierWork> {
        let (owner, literal) = self.owner(source)?;
        let support = FrontierWork::support(&powers);
        ensure!(powers.len() == 15, "native endpoint has wrong arity");
        if support == owner {
            ensure!(
                !literal || powers == source.powers,
                "literal owner changed its key"
            );
            return Ok(FrontierWork {
                stage: "apply".into(),
                powers,
                owner: Some(owner),
            });
        }
        ensure!(
            !literal
                && support
                    .iter()
                    .zip(&owner)
                    .all(|(child, owner)| !child || *owner)
                && support.iter().filter(|active| **active).count()
                    < FrontierWork::support(&source.powers)
                        .iter()
                        .filter(|active| **active)
                        .count(),
            "native endpoint is neither its owner nor a strict pinch"
        );
        FrontierWork::route(powers)
    }

    pub(super) fn step(
        &self,
        source: &FrontierWork,
        edges: Vec<Value>,
        guards: Vec<String>,
    ) -> Result<Value> {
        let zero = edges.is_empty();
        // Loaded Native owners cannot overlap installed zero sectors. Unknown
        // nonowner sources may still have a guarded exact zero before routing.
        ensure!(
            !zero || !self.owners.contains(&FrontierWork::support(&source.powers)),
            "literal installed owner unexpectedly has a zero outcome"
        );
        let labels = if zero {
            None
        } else {
            Some(self.owner(source)?)
        };
        Ok(
            json!({"source":source.document(),"zero":zero,"edges":edges,"nonzero_conditions":guards,
            "owner_label":labels.as_ref().map(|(owner,_)|Self::label(owner)),
            "route_label":labels.as_ref().and_then(|(_,literal)|(!literal).then(|| Self::label(&FrontierWork::support(&source.powers))))}),
        )
    }

    pub(super) fn document(&self) -> Value {
        json!({"family_fingerprint":self.family,
            "owners":self.owners.iter().map(|mask|json!({"mask":Self::label(mask)})).collect::<Vec<_>>(),
            "initial_frontier_routes":self.routes.iter().map(|(source,owner)|json!({"source_mask":Self::label(source),"owner_mask":Self::label(owner)})).collect::<Vec<_>>()})
    }
}

/// Independently check graph closure/reachability and join saved hits only after
/// every root is known. Never infer a zero root count from incomplete evidence.
pub(super) fn validate_frontier_certificate(
    certificate: &Value,
    normalized: &Value,
    binding: &Value,
    compiled: &Value,
    trace_sha256: &str,
) -> Result<Value> {
    ensure!(
        certificate["schema"] == "vakint.initial-routing-frontier-hit-proof.v1"
            && certificate["complete"] == true
            && certificate["reason"].is_null()
            && certificate["semantic_binding"] == *binding
            && certificate["compiled_cohort"] == *compiled
            && sha256(certificate, "trace_sha256")? == trace_sha256,
        "incomplete or mismatched frontier certificate"
    );
    // The outer OwnerProof pins the whole certificate and actual executable/
    // source cohort. Its native outcomes are authenticated mapper output, not
    // mathematical authority inferred from a structurally plausible graph.
    ensure!(
        compiled["mapper_source_sha256"]
            == format!(
                "{:x}",
                Sha256::digest(include_bytes!("../authenticate_owner_roots.rs"))
            ),
        "recursive frontier requires this actual compiled mapper source"
    );
    // Mapper-recorded Native-selection provenance. Authority is the pinned
    // whole proof plus exact compiled mapper/data/cohort chain, not standalone
    // replay of transport matrices from the projected labels below.
    sha256(certificate, "selection_sha256")?;
    let selection = FrontierSelection::parse(&certificate["selection"], binding)?;
    let observed = array(normalized, "dispatches")?;
    let dispatches = array(certificate, "dispatches")?;
    ensure!(
        dispatches.len() == observed.len(),
        "frontier omits a dispatch"
    );
    let mut derived = Vec::new();
    for (proof, actual) in dispatches.iter().zip(observed) {
        ensure!(
            proof["ordinal"] == actual["ordinal"]
                && proof["family_fingerprint"] == actual["family_fingerprint"]
                && proof["family_fingerprint"] == selection.family
                && proof["seeds"] == actual["seeds"],
            "frontier dispatch or postcombine seed signature differs"
        );
        let seeds = array(actual, "seeds")?
            .iter()
            .map(FrontierWork::parse)
            .collect::<Result<BTreeSet<_>>>()?;
        ensure!(
            proof["seed_sha256"]
                == format!(
                    "{:x}",
                    Sha256::digest(serde_json::to_vec(&actual["seeds"])?)
                ),
            "frontier seed hash differs"
        );
        let mut steps = BTreeMap::new();
        for step in array(proof, "steps")? {
            let source = FrontierWork::parse(&step["source"])?;
            ensure!(
                source.stage == "route" && steps.insert(source, step).is_none(),
                "duplicate step or illegal Apply expansion"
            );
        }
        let mut visited = BTreeSet::new();
        let mut frontier = BTreeSet::new();
        let mut active = BTreeSet::new();
        fn visit(
            source: &FrontierWork,
            selection: &FrontierSelection,
            steps: &BTreeMap<FrontierWork, &Value>,
            visited: &mut BTreeSet<FrontierWork>,
            active: &mut BTreeSet<FrontierWork>,
            frontier: &mut BTreeSet<FrontierWork>,
        ) -> Result<()> {
            frontier.insert(source.clone());
            if visited.contains(source) {
                return Ok(());
            }
            ensure!(
                active.insert(source.clone()),
                "cycle in initial routing frontier"
            );
            let step = steps
                .get(source)
                .ok_or_else(|| eyre::eyre!("missing reachable Route step"))?;
            ensure!(
                array(step, "nonzero_conditions")?.iter().all(|guard| guard
                    .as_str()
                    .is_some_and(|guard| !guard.trim().is_empty() && guard.trim() != "0")),
                "invalid native ordered guard"
            );
            let edges = array(step, "edges")?;
            ensure!(
                step["zero"] == edges.is_empty(),
                "zero outcome conflicts with native endpoints"
            );
            if edges.is_empty() {
                ensure!(
                    !selection
                        .owners
                        .contains(&FrontierWork::support(&source.powers))
                        && step["owner_label"].is_null()
                        && step["route_label"].is_null(),
                    "zero outcome invents a routing label"
                );
            } else {
                let (owner, literal) = selection.owner(source)?;
                ensure!(
                    step["owner_label"] == FrontierSelection::label(&owner)
                        && step["route_label"]
                            == if literal {
                                Value::Null
                            } else {
                                json!(FrontierSelection::label(&FrontierWork::support(
                                    &source.powers
                                )))
                            },
                    "routing labels differ"
                );
                ensure!(
                    !literal || edges.len() == 1,
                    "literal owner is not identity"
                );
                let mut local = BTreeSet::new();
                for edge in edges {
                    let target = FrontierWork::parse(&edge["target"])?;
                    ensure!(
                        target == selection.endpoint(source, target.powers.clone())?
                            && local.insert(target.clone())
                            && !string(edge, "coefficient")?.trim().is_empty()
                            && string(edge, "coefficient")?.trim() != "0"
                            && (!literal || edge["coefficient"] == "1"),
                        "invalid nonzero native routing edge"
                    );
                    frontier.insert(target.clone());
                    if target.stage == "route" {
                        visit(&target, selection, steps, visited, active, frontier)?;
                    }
                }
            }
            active.remove(source);
            visited.insert(source.clone());
            Ok(())
        }
        for seed in &seeds {
            visit(
                seed,
                &selection,
                &steps,
                &mut visited,
                &mut active,
                &mut frontier,
            )?;
        }
        ensure!(
            visited.len() == steps.len(),
            "extraneous unreachable routing step"
        );
        let claimed = array(proof, "frontier")?
            .iter()
            .map(FrontierWork::parse)
            .collect::<Result<Vec<_>>>()?;
        ensure!(
            claimed == frontier.iter().cloned().collect::<Vec<_>>(),
            "frontier union is missing, duplicate or extraneous"
        );
        let hits = array(actual, "hits")?
            .iter()
            .map(FrontierWork::parse)
            .collect::<Result<Vec<_>>>()?;
        let root = hits.iter().filter(|hit| frontier.contains(*hit)).count();
        derived.push(json!({"ordinal":actual["ordinal"],"family_fingerprint":actual["family_fingerprint"],
            "root_frontier_hits":root,"proper_descendant_hits":hits.len()-root,"observed_hits":hits.len(),"native_persistent":actual["persistent"]}));
    }
    Ok(
        json!({"complete":true,"interpretation":"benchmark initial-routing-frontier v1; Native counters unchanged",
        "dispatches":derived,"trace_sha256":trace_sha256}),
    )
}

pub(super) fn declared_training_matches(certificate: &Value, workload: &Value) -> Result<Value> {
    let mut matches = Vec::new();
    for dispatch in array(certificate, "dispatches")? {
        for work in array(dispatch, "frontier")? {
            let work = FrontierWork::parse(work)?;
            for training in array(workload, "training")? {
                let powers: Vec<i64> = serde_json::from_value(training["powers"].clone())?;
                ensure!(powers.len() == 15, "declared training key has wrong arity");
                if work.powers == powers {
                    matches.push(json!({"ordinal":dispatch["ordinal"],
                        "training_id":string(training,"id")?,"work":work.document()}));
                }
            }
        }
    }
    Ok(
        json!({"scope":"exact-declared-training-powers-only","positive_matches":matches,
        "absence_establishes_training_alias_disjointness":false,
        "changes_workload_roles_or_exclusions":false}),
    )
}

impl OwnerProof {
    fn read(workload: &Workload, workload_path: &Path) -> Result<Option<(Value, PathBuf, Value)>> {
        let Some(reference) = &workload.owner_proof else {
            return Ok(None);
        };
        let path = resolved(
            workload_path.parent().unwrap_or(Path::new(".")),
            &reference.path.to_string_lossy(),
        );
        let bytes = fs::read(&path)?;
        ensure!(
            format!("{:x}", Sha256::digest(&bytes)) == reference.sha256,
            "pinned owner-proof hash differs"
        );
        let document: Value = serde_json::from_slice(&bytes)?;
        ensure!(
            document["schema"] == "vakint.public-lowering-owner-map-proof.v1",
            "unknown completed owner-proof schema"
        );
        let identity = json!({"path":path,"sha256":reference.sha256,"mapper_executable_sha256":document["mapper_executable_sha256"],"mapper_source_sha256":document["mapper_source_sha256"],"cohort_receipt_sha256":document["cohort_receipt_sha256"],"compiled_cohort":document["compiled_cohort"]});
        Ok(Some((document, path, identity)))
    }

    pub(super) fn pinned_identity(
        workload: &Workload,
        workload_path: &Path,
        artifacts: &Path,
    ) -> Result<Value> {
        Self::pinned_identity_with_bindings(
            workload,
            workload_path,
            artifacts,
            &compiled_bindings(),
        )
    }

    fn pinned_identity_with_bindings(
        workload: &Workload,
        workload_path: &Path,
        artifacts: &Path,
        expected: &[(&str, Option<&str>); 5],
    ) -> Result<Value> {
        let Some((document, path, identity)) = Self::read(workload, workload_path)? else {
            return Ok(Value::Null);
        };
        Self::validate_data(&document, artifacts)?;
        let binding = document["semantic_binding"].clone();
        let proof = Self::validate(document, workload, &binding, identity, expected)?;
        proof.validate_evidence(path.parent().unwrap_or(Path::new(".")), workload)?;
        Ok(proof.identity)
    }

    pub(super) fn load(
        workload: &Workload,
        workload_path: &Path,
        artifacts: &Path,
        binding: &Value,
    ) -> Result<Option<Self>> {
        let Some((document, path, identity)) = Self::read(workload, workload_path)? else {
            return Ok(None);
        };
        let base = path.parent().unwrap_or(Path::new("."));
        Self::validate_data(&document, artifacts)?;
        let proof = Self::validate(document, workload, binding, identity, &compiled_bindings())?;
        proof.validate_evidence(base, workload)?;
        ensure!(
            digest(&path)? == workload.owner_proof.as_ref().unwrap().sha256,
            "owner proof changed during validation"
        );
        Ok(Some(proof))
    }

    fn validate_data(document: &Value, artifacts: &Path) -> Result<()> {
        for (relative, field) in [
            ("manifest.json", "data_manifest_sha256"),
            ("5L/manifest.json", "loop_manifest_sha256"),
            ("5L/ibps/application.json", "application_sha256"),
        ] {
            ensure!(
                digest(&artifacts.join(relative))? == sha256(document, field)?,
                "owner-proof data or application differs"
            );
        }
        Ok(())
    }

    fn validate_evidence(&self, base: &Path, workload: &Workload) -> Result<()> {
        let mut frozen = Vec::new();
        for (path, hash) in [
            ("frozen_suite_path", "frozen_suite_sha256"),
            ("reserve_path", "reserve_sha256"),
            ("protocol_path", "protocol_sha256"),
        ] {
            let path = evidence(&self.document, base, path, hash)?;
            frozen.push(serde_json::from_slice::<Value>(&fs::read(path)?)?);
        }
        let cohort_path = evidence(
            &self.document,
            base,
            "cohort_receipt_path",
            "cohort_receipt_sha256",
        )?;
        let receipt: Value = serde_json::from_slice(&fs::read(cohort_path)?)?;
        ensure!(
            receipt == self.document["compiled_cohort"],
            "compiled-cohort receipt differs"
        );
        // Executable hashes are mandatory provenance, but archived ELFs are not.
        // If a binary is available, still authenticate its bytes.
        available_executable(
            &self.document,
            base,
            "mapper_executable_path",
            "mapper_executable_sha256",
        )?;
        available_executable(
            &self.document,
            base,
            "measurement_executable_path",
            "measurement_executable_sha256",
        )?;
        for root in workload.tuning.iter().chain(&workload.held_out) {
            let row = self.row(&root.id)?;
            let workload_path = evidence(row, base, "workload_path", "workload_sha256")?;
            evidence(row, base, "source_census_path", "source_census_sha256")?;
            let source_workload: Value = serde_json::from_slice(&fs::read(workload_path)?)?;
            let role = string(row, "logical_cohort")?;
            let source_cases = array(
                &source_workload,
                if role == "tuning" {
                    "tuning"
                } else {
                    "held_out"
                },
            )?;
            let declarations = source_cases
                .iter()
                .filter(|case| case["id"] == root.id)
                .collect::<Vec<_>>();
            let fixed_cases = if role == "tuning" {
                array(&frozen[0], "tuning")?.iter().collect::<Vec<_>>()
            } else {
                array(&frozen[0], "held_out")?
                    .iter()
                    .chain(array(&frozen[1], "held_out_reserves")?)
                    .collect()
            };
            let fixed = fixed_cases
                .into_iter()
                .filter(|case| case["id"] == root.id)
                .collect::<Vec<_>>();
            ensure!(
                declarations.len() == 1
                    && fixed.len() == 1
                    && declarations[0]["public_input"] == json!(root.public_input)
                    && declarations[0]["powers"] == json!(root.powers)
                    && fixed[0]["public_input"] == declarations[0]["public_input"]
                    && fixed[0]["powers"] == declarations[0]["powers"],
                "historical source identity or role differs from the pinned case"
            );
            let receipt_path = evidence(row, base, "receipt_path", "receipt_sha256")?;
            let receipt: Value = serde_json::from_slice(&fs::read(&receipt_path)?)?;
            ensure!(
                receipt["success"] == true
                    && receipt["worker_success"] == true
                    && receipt["fresh_process"] == true
                    && receipt["variant"] == "uncached"
                    && receipt["sample"] == 0
                    && receipt["supervision"]["success"] == true
                    && receipt["supervision"]["exit_code"] == 0
                    && receipt["supervision"]["deadline_exhausted"] == false,
                "owner proof requires a completed successful public report"
            );
            ensure!(
                receipt["case_id"] == root.id
                    && sha256(&receipt, "input_sha256")? == sha256(row, "input_sha256")?
                    && receipt["actual_lowered_roots"] == row["actual_lowered_roots"]
                    && receipt["cache_statistics"]["binding"] == self.document["semantic_binding"]
                    && receipt["identity"]["data_manifest_sha256"] == workload.data_manifest_sha256
                    && sha256(&receipt["identity"], "workload_sha256")?
                        == sha256(row, "workload_sha256")?
                    && sha256(&receipt["identity"], "source_census_sha256")?
                        == sha256(row, "source_census_sha256")?
                    && sha256(&receipt["identity"], "executable_sha256")?
                        == sha256(&self.document, "measurement_executable_sha256")?
                    && receipt["identity"]["cache_sha256"].is_null(),
                "owner-proof public report differs"
            );
            ensure!(
                digest(&receipt_path.with_file_name("value.symbolica"))?
                    == string(&receipt, "value_sha256")?
                    && digest(&receipt_path.with_file_name("guards.json"))?
                        == string(&receipt, "guards_sha256")?,
                "completed public result or guards changed"
            );
        }
        Ok(())
    }

    fn validate(
        document: Value,
        workload: &Workload,
        binding: &Value,
        identity: Value,
        expected: &[(&str, Option<&str>); 5],
    ) -> Result<Self> {
        workload.validate_final_selection()?;
        cohort(&document, binding, expected)?;
        ensure!(
            document["semantic_binding"] == *binding
                && binding.is_object()
                && binding["physical_arity"] == 15
                && binding["layer"] == "saved-application-route-apply-v1"
                && document["data_manifest_sha256"] == workload.data_manifest_sha256,
            "owner proof differs from the loaded application"
        );
        let mut indexed = BTreeMap::new();
        for row in array(&document, "rows")? {
            ensure!(
                indexed.insert(string(row, "id")?, row).is_none(),
                "duplicate owner-proof case"
            );
        }
        let mut signatures = BTreeSet::new();
        let mut selected = Vec::new();
        for (cohort, cases) in [
            ("tuning", &workload.tuning),
            ("held-out", &workload.held_out),
        ] {
            for root in cases {
                let row = indexed
                    .get(root.id.as_str())
                    .ok_or_else(|| eyre::eyre!("owner proof omits {}", root.id))?;
                let source = format!(
                    "{:x}",
                    Sha256::digest(
                        root.public_input
                            .as_ref()
                            .ok_or_else(|| eyre::eyre!("missing public input"))?
                            .as_bytes()
                    )
                );
                ensure!(
                    row["logical_cohort"] == cohort
                        && row["public_input_source_sha256"] == source
                        && row["reduction_succeeded"] == true
                        && row["complete_owner_transport"] == true
                        && row["mapping_error"].is_null()
                        && row["application_binding"] == *binding,
                    "incomplete or mismatched owner case {}",
                    root.id
                );
                let raw = roots(&row["actual_lowered_roots"])?;
                let protected = roots(&row["protected_owner_roots"])?;
                ensure!(
                    !raw.is_empty()
                        && !protected.is_empty()
                        && raw == roots(&json!(root.lowered_roots))?
                        && protected == roots(&json!(root.protected_roots))?,
                    "frozen owner sets differ"
                );
                let mut covered = BTreeSet::new();
                let mut endpoints = BTreeSet::new();
                for transport in array(row, "native_transports")? {
                    let key: Vec<i64> = serde_json::from_value(transport["powers"].clone())?;
                    ensure!(covered.insert(key), "duplicate native source transport");
                    let mut local = BTreeSet::new();
                    for endpoint in array(transport, "owner_endpoints")? {
                        let key: Vec<i64> = serde_json::from_value(endpoint["powers"].clone())?;
                        ensure!(
                            key.len() == 15
                                && local.insert(key.clone())
                                && !string(endpoint, "coefficient")?.trim().is_empty()
                                && string(endpoint, "coefficient")?.trim() != "0"
                                && endpoint["recorded_coverage_admitted"].is_boolean(),
                            "invalid native endpoint"
                        );
                        endpoints.insert(key);
                    }
                    ensure!(
                        array(transport, "nonzero_conditions")?
                            .iter()
                            .all(|guard| guard
                                .as_str()
                                .is_some_and(
                                    |guard| !guard.trim().is_empty() && guard.trim() != "0"
                                )),
                        "invalid transport guard"
                    );
                }
                ensure!(
                    covered == raw && endpoints == protected,
                    "native transports do not cover the exact frozen sets"
                );
                let signature = format!(
                    "{:x}",
                    Sha256::digest(serde_json::to_vec(&row["protected_owner_roots"])?)
                );
                ensure!(
                    row["protected_root_set_sha256"] == signature,
                    "owner-set signature differs"
                );
                if cohort == "held-out" {
                    signatures.insert(signature);
                }
                ensure!(
                    workload
                        .training
                        .iter()
                        .all(|root| !protected.contains(&root.powers)),
                    "training contains a protected owner root"
                );
                selected.push((root.id.clone(), cohort, protected));
            }
        }
        ensure!(
            signatures.len() >= 12,
            "require twelve distinct nonempty held-out owner-root sets"
        );
        let mut correlations = Vec::new();
        for (index, (a, cohort_a, left)) in selected.iter().enumerate() {
            for (b, cohort_b, right) in selected.iter().skip(index + 1) {
                let shared = left.intersection(right).collect::<Vec<_>>();
                if !shared.is_empty() {
                    correlations.push(json!({"case_a":a,"cohort_a":cohort_a,"case_b":b,"cohort_b":cohort_b,"shared_owner_roots":shared,"identical_owner_sets":left==right}));
                }
            }
        }
        Ok(Self {
            document,
            identity,
            correlations: json!(correlations),
        })
    }
    fn row(&self, id: &str) -> Result<&Value> {
        array(&self.document, "rows")?
            .iter()
            .find(|row| row["id"] == id)
            .ok_or_else(|| eyre::eyre!("missing owner case"))
    }
    pub(super) fn observe(
        &self,
        id: &str,
        raw: &BTreeSet<Vec<i64>>,
        protected: &[Vec<i64>],
    ) -> Result<()> {
        let row = self.row(id)?;
        ensure!(
            raw == &roots(&row["actual_lowered_roots"])?
                && roots(&json!(protected))? == roots(&row["protected_owner_roots"])?,
            "observed public roots differ from the pinned owner proof"
        );
        Ok(())
    }
    pub(super) fn identity(&self) -> &Value {
        &self.identity
    }
    pub(super) fn correlations(&self) -> &Value {
        &self.correlations
    }

    /// Optional v1 evidence is checked independently from the conservative
    /// legacy owner sets. Its absence never upgrades legacy counter semantics.
    fn validate_initial_frontier(&self, row: &Value, base: &Path) -> Result<Value> {
        let proof = &row["initial_routing_frontier"];
        let receipt_path = evidence(
            proof,
            base,
            "cache_hit_receipt_path",
            "cache_hit_receipt_sha256",
        )?;
        let receipt: Value = serde_json::from_slice(&fs::read(&receipt_path)?)?;
        ensure!(
            receipt["success"] == true
                && receipt["worker_success"] == true
                && receipt["fresh_process"] == true
                && receipt["supervision"]["success"] == true
                && receipt["supervision"]["exit_code"] == 0
                && receipt["supervision"]["deadline_exhausted"] == false
                && receipt["case_id"] == row["id"]
                && receipt["input_sha256"] == row["input_sha256"]
                && proof["input_sha256"] == row["input_sha256"]
                && receipt["cache_statistics"]["binding"] == self.document["semantic_binding"]
                && receipt["identity"]["data_manifest_sha256"]
                    == self.document["data_manifest_sha256"]
                && receipt["identity"]["executable_sha256"]
                    == self.document["measurement_executable_sha256"]
                && receipt["identity"]["workload_sha256"] == row["workload_sha256"]
                && receipt["identity"]["source_census_sha256"] == row["source_census_sha256"],
            "saved frontier diagnostic differs from original case/cohort"
        );
        let original_path = evidence(row, base, "receipt_path", "receipt_sha256")?;
        let original: Value = serde_json::from_slice(&fs::read(original_path)?)?;
        validate_frontier_outputs(&original, &receipt)?;
        for (name, field) in [
            ("value.symbolica", "value_sha256"),
            ("guards.json", "guards_sha256"),
        ] {
            ensure!(
                sha256(&receipt, field)? == sha256(&original, field)?
                    && digest(&receipt_path.with_file_name(name))? == sha256(&receipt, field)?,
                "diagnostic output/canonical guard set differs"
            );
        }
        let trace = evidence(proof, base, "trace_path", "trace_sha256")?;
        let from_receipt = resolved(
            receipt_path.parent().unwrap_or(Path::new(".")),
            string(&receipt["cache_hit_diagnostic"], "path")?,
        );
        ensure!(
            trace == from_receipt
                && sha256(&receipt["cache_hit_diagnostic"], "sha256")?
                    == sha256(proof, "trace_sha256")?,
            "frontier trace differs from successful receipt binding"
        );
        let normalized =
            validate_hit_diagnostic(&serde_json::from_slice::<Value>(&fs::read(trace)?)?)?;
        for field in [
            "lookups",
            "route_hits",
            "apply_hits",
            "requested_root_hits",
            "intermediate_hits",
            "misses",
        ] {
            ensure!(
                normalized["persistent"][field] == receipt["cache_statistics"]["persistent"][field],
                "saved dispatch counters differ from final raw receipt"
            );
        }
        let mut derived = validate_frontier_certificate(
            proof,
            &normalized,
            &self.document["semantic_binding"],
            &self.document["compiled_cohort"],
            sha256(proof, "trace_sha256")?,
        )?;
        ensure!(
            proof["derived"] == derived,
            "claimed derived hit counts differ from independent join"
        );
        if proof.get("declared_training_matches").is_some() {
            let workload_path = evidence(row, base, "workload_path", "workload_sha256")?;
            let workload: Value = serde_json::from_slice(&fs::read(workload_path)?)?;
            let matches = declared_training_matches(proof, &workload)?;
            ensure!(
                proof["declared_training_matches"] == matches,
                "declared training overlap differs"
            );
            derived["declared_training_matches"] = matches;
        }
        Ok(derived)
    }

    pub(super) fn initial_routing_frontier(&self, id: &str) -> Result<Value> {
        let row = self.row(id)?;
        if row.get("initial_routing_frontier").is_none() {
            return Ok(Value::Null);
        }
        let base = Path::new(string(&self.identity, "path")?)
            .parent()
            .unwrap_or(Path::new("."));
        Ok(match self.validate_initial_frontier(row, base) {
            Ok(derived) => derived,
            Err(error) => {
                json!({"complete":false,"reason":error.to_string(),"scope":"unclassified",
                "root_frontier_hits":null,"proper_descendant_hits":null})
            }
        })
    }
}

#[cfg(test)]
pub(super) fn hit_diagnostic_fixture() -> Value {
    let route = json!({"stage":"route","powers":vec![1_i64;15],"owner_sector":null});
    let apply = json!({"stage":"apply","powers":vec![1_i64;15],"owner_sector":vec![true;15]});
    let event = |name: &str, sequence: u64| {
        json!({"schema":"vakint.cache-hit-events.v1","event":name,
        "sequence":sequence,"family_fingerprint":"fixture","work":null,"rootwork":null,
        "requested_roots":null,"completed_roots":null,"persistent":null,"status":"running"})
    };
    let start = event("dispatch_start", 0);
    let mut hit = event("backing_hit", 2);
    hit["work"] = apply;
    let mut root = event("root_collected", 3);
    root["work"] = route.clone();
    root["rootwork"] = route;
    let mut finish = event("dispatch_finish", 4);
    finish["status"] = json!("succeeded");
    finish["requested_roots"] = json!(1);
    finish["completed_roots"] = json!(1);
    finish["persistent"] = json!({"lookups":1,"route_hits":0,"apply_hits":1,"requested_root_hits":0,"intermediate_hits":1,"misses":0});
    json!({"schema":"vakint.cache-hit-diagnostic.v1","interval_ns":0,"input_dispatch_count":1,
        "success":true,"error":null,"events":[start,hit,root,finish,
            {"schema":"vakint.cache-hit-events.v1","event":"observer_finish","sequence":null,
                "observer_status":"Active","operation_succeeded":true,"interval_ns":0}]})
}

#[cfg(test)]
pub(super) mod tests {
    use super::*;

    fn frontier_fixture(normalized: &Value) -> (Value, Value, Value) {
        let binding = json!({"physical_arity":15,"family_fingerprint":"fixture"});
        let compiled = json!({"native_revision":"4190f2c6f3158e3fc3cbcb1937e12e5971da9175",
            "mapper_source_sha256":format!("{:x}",Sha256::digest(include_bytes!("../authenticate_owner_roots.rs")))});
        let selection = FrontierSelection::parse(
            &json!({"family_fingerprint":"fixture",
            "owners":[{"mask":"111111111111111"}],"initial_frontier_routes":[]}),
            &binding,
        )
        .unwrap();
        let dispatches=normalized["dispatches"].as_array().unwrap().iter().map(|dispatch| {
            let mut frontier=BTreeSet::new();
            let steps=dispatch["seeds"].as_array().unwrap().iter().map(|seed|{
                let source=FrontierWork::parse(seed).unwrap();
                let target=selection.endpoint(&source,source.powers.clone()).unwrap();
                frontier.insert(source.clone());frontier.insert(target.clone());
                selection.step(&source,vec![json!({"target":target.document(),"coefficient":"1"})],vec!["d-2".into()]).unwrap()
            }).collect::<Vec<_>>();
            json!({"ordinal":dispatch["ordinal"],"family_fingerprint":"fixture","seeds":dispatch["seeds"],
                "seed_sha256":format!("{:x}",Sha256::digest(serde_json::to_vec(&dispatch["seeds"]).unwrap())),
                "steps":steps,"frontier":frontier.iter().map(FrontierWork::document).collect::<Vec<_>>()})
        }).collect::<Vec<_>>();
        (
            json!({"schema":"vakint.initial-routing-frontier-hit-proof.v1","complete":true,"reason":null,
            "selection":selection.document(),"selection_sha256":"a".repeat(64),"trace_sha256":"b".repeat(64),
            "semantic_binding":binding,"compiled_cohort":compiled,"dispatches":dispatches}),
            binding,
            compiled,
        )
    }

    #[test]
    fn optional_frontier_refuses_missing_or_unclosed_outputs_and_event_fields() {
        let closed = json!({"retained_result_available":true,"catalogue_closure":{
            "closed":true,"unknown_count":0,"unknown_masters":[]}});
        validate_frontier_outputs(&closed, &closed).unwrap();
        for side in 0..2 {
            for field in ["retained_result_available", "catalogue_closure"] {
                let mut pair = [closed.clone(), closed.clone()];
                pair[side].as_object_mut().unwrap().remove(field);
                assert!(validate_frontier_outputs(&pair[0], &pair[1]).is_err());
            }
            for (field, value) in [
                ("closed", json!(false)),
                ("unknown_count", json!(1)),
                ("unknown_masters", json!(["unknown"])),
            ] {
                let mut pair = [closed.clone(), closed.clone()];
                pair[side]["catalogue_closure"][field] = value;
                assert!(validate_frontier_outputs(&pair[0], &pair[1]).is_err());
            }
        }
        for field in [
            "work",
            "rootwork",
            "requested_roots",
            "completed_roots",
            "persistent",
            "status",
        ] {
            let mut missing = hit_diagnostic_fixture();
            missing["events"][0].as_object_mut().unwrap().remove(field);
            assert!(validate_hit_diagnostic(&missing).is_err(), "{field}");
        }
    }

    #[test]
    fn collected_root_requires_its_current_clone_and_actual_mapper_authority() {
        for current in [
            Value::Null,
            json!({"stage":"route","powers":vec![2_i64;15],"owner_sector":null}),
        ] {
            let mut missing = hit_diagnostic_fixture();
            missing["events"][2]["work"] = current;
            assert!(validate_hit_diagnostic(&missing).is_err());
        }
        let normalized = validate_hit_diagnostic(&hit_diagnostic_fixture()).unwrap();
        let (mut proof, binding, mut compiled) = frontier_fixture(&normalized);
        compiled["mapper_source_sha256"] = json!("f".repeat(64));
        proof["compiled_cohort"] = compiled.clone();
        assert!(
            validate_frontier_certificate(
                &proof,
                &normalized,
                &binding,
                &compiled,
                &"b".repeat(64)
            )
            .is_err()
        );
    }

    #[test]
    fn declared_training_matches_report_only_exact_positive_stage_identities() {
        let normalized = validate_hit_diagnostic(&hit_diagnostic_fixture()).unwrap();
        let (proof, _, _) = frontier_fixture(&normalized);
        let workload = json!({"training":[{"id":"selected","powers":vec![1_i64;15]},
            {"id":"different","powers":vec![2_i64;15]}]});
        let matches = declared_training_matches(&proof, &workload).unwrap();
        assert_eq!(matches["scope"], "exact-declared-training-powers-only");
        assert_eq!(matches["positive_matches"].as_array().unwrap().len(), 2);
        assert!(
            matches["positive_matches"]
                .as_array()
                .unwrap()
                .iter()
                .all(|entry| entry["training_id"] == "selected")
        );
        assert_eq!(
            matches["absence_establishes_training_alias_disjointness"],
            false
        );
        assert_eq!(matches["changes_workload_roles_or_exclusions"], false);
    }

    #[test]
    fn hit_stream_requires_every_dispatch_identity_counter_and_acknowledgement() {
        let original = hit_diagnostic_fixture();
        let normalized = validate_hit_diagnostic(&original).unwrap();
        assert_eq!(
            normalized["dispatches"][0]["seeds"]
                .as_array()
                .unwrap()
                .len(),
            1
        );
        for (path, value) in [
            (vec!["events", "1", "work"], Value::Null),
            (
                vec!["events", "1", "work", "owner_sector"],
                json!(vec![false; 15]),
            ),
            (vec!["events", "3", "persistent", "apply_hits"], json!(0)),
            (vec!["events", "3", "completed_roots"], json!(0)),
            (vec!["events", "4", "observer_status"], json!("Detached")),
            (vec!["input_dispatch_count"], json!(2)),
        ] {
            let mut changed = original.clone();
            let mut cursor = &mut changed;
            for key in path {
                cursor = if let Ok(index) = key.parse::<usize>() {
                    &mut cursor[index]
                } else {
                    &mut cursor[key]
                };
            }
            *cursor = value;
            assert!(validate_hit_diagnostic(&changed).is_err());
        }
        let mut missing = original.clone();
        missing["events"].as_array_mut().unwrap().remove(0);
        assert!(validate_hit_diagnostic(&missing).is_err());
        let mut duplicate = original.clone();
        let root = duplicate["events"][2].clone();
        duplicate["events"].as_array_mut().unwrap().insert(3, root);
        assert!(validate_hit_diagnostic(&duplicate).is_err());
        let mut dropped = original;
        dropped["events"] = json!([dropped["events"][4]]);
        assert!(validate_hit_diagnostic(&dropped).is_err());
    }

    #[test]
    fn native_intermediate_apply_hit_is_frontier_root_without_rewriting_native_counters() {
        let normalized = validate_hit_diagnostic(&hit_diagnostic_fixture()).unwrap();
        let (proof, binding, compiled) = frontier_fixture(&normalized);
        let joined = validate_frontier_certificate(
            &proof,
            &normalized,
            &binding,
            &compiled,
            &"b".repeat(64),
        )
        .unwrap();
        assert_eq!(joined["dispatches"][0]["root_frontier_hits"], 1);
        assert_eq!(joined["dispatches"][0]["proper_descendant_hits"], 0);
        assert_eq!(
            joined["dispatches"][0]["native_persistent"]["requested_root_hits"],
            0
        );
        assert_eq!(
            joined["dispatches"][0]["native_persistent"]["intermediate_hits"],
            1
        );
        for field in ["compiled_cohort", "semantic_binding", "trace_sha256"] {
            let mut changed = proof.clone();
            changed[field] = Value::Null;
            assert!(
                validate_frontier_certificate(
                    &changed,
                    &normalized,
                    &binding,
                    &compiled,
                    &"b".repeat(64)
                )
                .is_err()
            );
        }
        let mut changed = proof;
        changed["dispatches"][0]["steps"][0]["edges"][0]["coefficient"] = json!("2");
        assert!(
            validate_frontier_certificate(
                &changed,
                &normalized,
                &binding,
                &compiled,
                &"b".repeat(64)
            )
            .is_err()
        );
    }

    #[test]
    fn join_uses_global_collected_roots_not_hit_or_caller_order() {
        let mut trace = hit_diagnostic_fixture();
        let mut a = vec![1_i64; 15];
        a[0] = 2;
        trace["events"][2]["rootwork"]["powers"] = json!(a.clone());
        trace["events"][2]["work"]["powers"] = json!(a.clone());
        let only_a = validate_hit_diagnostic(&trace).unwrap();
        let (proof, binding, compiled) = frontier_fixture(&only_a);
        assert_eq!(
            validate_frontier_certificate(&proof, &only_a, &binding, &compiled, &"b".repeat(64))
                .unwrap()["dispatches"][0]["proper_descendant_hits"],
            1
        );
        let mut b = trace["events"][2].clone();
        b["sequence"] = json!(4);
        b["rootwork"]["powers"] = json!(vec![1_i64; 15]);
        b["work"]["powers"] = json!(vec![1_i64; 15]);
        trace["events"].as_array_mut().unwrap().insert(3, b);
        trace["events"][4]["sequence"] = json!(5);
        trace["events"][4]["requested_roots"] = json!(2);
        trace["events"][4]["completed_roots"] = json!(2);
        let batch = validate_hit_diagnostic(&trace).unwrap();
        let (proof, binding, compiled) = frontier_fixture(&batch);
        assert_eq!(
            validate_frontier_certificate(&proof, &batch, &binding, &compiled, &"b".repeat(64))
                .unwrap()["dispatches"][0]["root_frontier_hits"],
            1
        );
        trace["events"][2]["rootwork"]["powers"] = json!(vec![1_i64; 15]);
        trace["events"][2]["work"]["powers"] = json!(vec![1_i64; 15]);
        trace["events"][3]["rootwork"]["powers"] = json!(a.clone());
        trace["events"][3]["work"]["powers"] = json!(a);
        assert_eq!(validate_hit_diagnostic(&trace).unwrap(), batch);
    }

    #[test]
    fn zero_route_and_strict_permuted_pinch_preserve_stage_and_guards() {
        let binding = json!({"physical_arity":15,"family_fingerprint":"fixture"});
        let owner = "011000000000000";
        let selection=FrontierSelection::parse(&json!({"family_fingerprint":"fixture","owners":[{"mask":owner},{"mask":"010000000000000"}],
            "initial_frontier_routes":[{"source_mask":"110000000000000","owner_mask":owner}]}),&binding).unwrap();
        let mut powers = vec![0; 15];
        powers[0] = 1;
        powers[1] = 1;
        let source = FrontierWork::route(powers).unwrap();
        let mut pinched = vec![0; 15];
        pinched[1] = 1;
        let endpoint = selection.endpoint(&source, pinched.clone()).unwrap();
        assert_eq!(endpoint.stage, "route");
        assert_eq!(
            selection.endpoint(&endpoint, pinched).unwrap().stage,
            "apply"
        );
        assert!(selection.step(&endpoint, vec![], vec![]).is_err());
        let zero = selection
            .step(
                &FrontierWork::route(vec![0; 15]).unwrap(),
                vec![],
                vec!["d-3".into(), "d-2".into()],
            )
            .unwrap();
        assert_eq!(zero["zero"], true);
        assert!(zero["owner_label"].is_null());
        assert_eq!(zero["nonzero_conditions"], json!(["d-3", "d-2"]));
        let mut invalid = vec![0; 15];
        invalid[0] = 1;
        invalid[2] = 1;
        assert!(selection.endpoint(&source, invalid).is_err());
        assert!(FrontierSelection::parse(&json!({"family_fingerprint":"fixture","owners":[{"mask":owner},{"mask":owner}],"initial_frontier_routes":[]}),&binding).is_err());
    }

    #[test]
    fn certificate_rejects_missing_extraneous_or_apply_expanded_steps() {
        let normalized = validate_hit_diagnostic(&hit_diagnostic_fixture()).unwrap();
        let (proof, binding, compiled) = frontier_fixture(&normalized);
        for field in ["steps", "frontier", "seeds"] {
            let mut changed = proof.clone();
            changed["dispatches"][0][field] = json!([]);
            assert!(
                validate_frontier_certificate(
                    &changed,
                    &normalized,
                    &binding,
                    &compiled,
                    &"b".repeat(64)
                )
                .is_err()
            );
        }
        let mut changed = proof.clone();
        let extra = changed["dispatches"][0]["steps"][0].clone();
        changed["dispatches"][0]["steps"]
            .as_array_mut()
            .unwrap()
            .push(extra);
        assert!(
            validate_frontier_certificate(
                &changed,
                &normalized,
                &binding,
                &compiled,
                &"b".repeat(64)
            )
            .is_err()
        );
        let mut changed = proof;
        changed["dispatches"][0]["steps"][0]["source"] =
            changed["dispatches"][0]["steps"][0]["edges"][0]["target"].clone();
        assert!(
            validate_frontier_certificate(
                &changed,
                &normalized,
                &binding,
                &compiled,
                &"b".repeat(64)
            )
            .is_err()
        );
    }
    fn expected_bindings() -> [(&'static str, Option<&'static str>); 5] {
        [
            (
                "main_revision",
                Some("1111111111111111111111111111111111111111"),
            ),
            (
                "native_revision",
                Some("2222222222222222222222222222222222222222"),
            ),
            (
                "kernel_revision",
                Some("3333333333333333333333333333333333333333"),
            ),
            (
                "main_source_snapshot_sha256",
                Some("aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa"),
            ),
            (
                "native_source_snapshot_sha256",
                Some("aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa"),
            ),
        ]
    }
    fn cohort(document: &Value, binding: &Value) -> Result<()> {
        super::cohort(document, binding, &expected_bindings())
    }
    fn validate(
        document: Value,
        workload: &Workload,
        binding: &Value,
        identity: Value,
    ) -> Result<OwnerProof> {
        OwnerProof::validate(document, workload, binding, identity, &expected_bindings())
    }
    fn pinned_identity(workload: &Workload, path: &Path, artifacts: &Path) -> Result<Value> {
        OwnerProof::pinned_identity_with_bindings(workload, path, artifacts, &expected_bindings())
    }
    #[test]
    fn absent_or_empty_current_bindings_cannot_accept_a_valid_proof() {
        let workload = super::super::workload::tests::selected();
        let (document, binding) = document(&workload);
        for index in 0..5 {
            for absent in [None, Some(""), Some(" \n\t")] {
                let mut expected = expected_bindings();
                expected[index].1 = absent;
                assert!(super::cohort(&document, &binding, &expected).is_err());
            }
        }
    }

    fn document(workload: &Workload) -> (Value, Value) {
        let binding = json!({"artifact_digest":vec![1;32],"family_fingerprint":"fixture","physical_arity":15,"layer":"saved-application-route-apply-v1"});
        let rows=[("tuning",&workload.tuning),("held-out",&workload.held_out)].into_iter().flat_map(|(cohort,cases)|cases.iter().map(move |root|(cohort,root))).map(|(cohort,root)|{
            json!({"id":root.id,"logical_cohort":cohort,"public_input_source_sha256":format!("{:x}",Sha256::digest(root.public_input.as_ref().unwrap().as_bytes())),
                "reduction_succeeded":true,"complete_owner_transport":true,"mapping_error":null,"application_binding":binding,
                "actual_lowered_roots":root.lowered_roots,"protected_owner_roots":root.protected_roots,
                "protected_root_set_sha256":format!("{:x}",Sha256::digest(serde_json::to_vec(&root.protected_roots).unwrap())),
                "native_transports":root.lowered_roots.iter().map(|key|json!({"powers":key,"owner_endpoints":root.protected_roots.iter().map(|key|json!({"powers":key,"coefficient":"1","recorded_coverage_admitted":true})).collect::<Vec<_>>(),"nonzero_conditions":[]})).collect::<Vec<_>>()})
        }).collect::<Vec<_>>();
        (
            {
                let hash = "a".repeat(64);
                let compiled = json!({"main_revision":"1".repeat(40),"native_revision":"2".repeat(40),"kernel_revision":"3".repeat(40),
                    "main_source_snapshot_sha256":hash,"native_source_snapshot_sha256":hash,
                    "mapper_source_sha256":hash,"mapper_executable_sha256":hash,"measurement_executable_sha256":hash,
                    "data_manifest_sha256":workload.data_manifest_sha256,"loop_manifest_sha256":hash,"application_sha256":hash,"semantic_binding":binding});
                json!({"schema":"vakint.public-lowering-owner-map-proof.v1","semantic_binding":binding,
                    "data_manifest_sha256":workload.data_manifest_sha256,"loop_manifest_sha256":hash,"application_sha256":hash,
                    "mapper_source_sha256":hash,"mapper_executable_sha256":hash,"measurement_executable_sha256":hash,
                    "cohort_receipt_path":"cohort.json","cohort_receipt_sha256":hash,"compiled_cohort":compiled,"rows":rows})
            },
            binding,
        )
    }
    pub(crate) fn fixture() -> OwnerProof {
        let workload = super::super::workload::tests::selected();
        let (document, binding) = document(&workload);
        validate(
            document,
            &workload,
            &binding,
            json!({"sha256":"synthetic fixture"}),
        )
        .unwrap()
    }
    #[test]
    fn missing_null_or_mismatched_cohort_bindings_never_authenticate() {
        let workload = super::super::workload::tests::selected();
        let (original, binding) = document(&workload);
        cohort(&original, &binding).unwrap();
        for field in [
            "mapper_source_sha256",
            "mapper_executable_sha256",
            "measurement_executable_sha256",
            "data_manifest_sha256",
            "loop_manifest_sha256",
            "application_sha256",
        ] {
            for value in [Value::Null, json!(""), json!("wrong")] {
                let mut changed = original.clone();
                changed[field] = value.clone();
                changed["compiled_cohort"][field] = value;
                assert!(cohort(&changed, &binding).is_err(), "{field}");
            }
        }
        for field in [
            "main_revision",
            "native_revision",
            "kernel_revision",
            "main_source_snapshot_sha256",
            "native_source_snapshot_sha256",
        ] {
            let mut changed = original.clone();
            changed["compiled_cohort"][field] = Value::Null;
            assert!(cohort(&changed, &binding).is_err(), "{field}");
        }
        let mut changed = original;
        changed["compiled_cohort"]["semantic_binding"] = Value::Null;
        assert!(cohort(&changed, &binding).is_err());
    }

    #[test]
    fn self_consistent_other_compiled_cohort_is_rejected() {
        let workload = super::super::workload::tests::selected();
        let (original, binding) = document(&workload);
        for field in ["main_revision", "native_revision", "kernel_revision"] {
            let mut changed = original.clone();
            changed["compiled_cohort"][field] = json!("f".repeat(40));
            assert!(cohort(&changed, &binding).is_err());
        }
        for field in [
            "main_source_snapshot_sha256",
            "native_source_snapshot_sha256",
        ] {
            let mut changed = original.clone();
            changed["compiled_cohort"][field] = json!("f".repeat(64));
            assert!(cohort(&changed, &binding).is_err());
        }
    }

    #[test]
    fn source_only_provenance_does_not_require_archived_elf_files() {
        let document = json!({"binary_path":"certainly-absent-owner-proof-test-binary",
            "binary_sha256":"a".repeat(64)});
        available_executable(&document, Path::new("."), "binary_path", "binary_sha256").unwrap();
        let mut missing = document;
        missing["binary_sha256"] = Value::Null;
        assert!(
            available_executable(&missing, Path::new("."), "binary_path", "binary_sha256").is_err()
        );
    }

    #[test]
    fn resumed_success_cannot_omit_or_change_output_evidence() {
        let workload = super::super::workload::tests::selected();
        let root = &workload.held_out[0];
        let missing = json!({"success":true,"worker_success":true,"fresh_process":true,
            "supervision":{"success":true,"exit_code":0,"deadline_exhausted":false}});
        assert!(validate_completed_outputs(&missing, Path::new("."), root).is_err());
        let failed = json!({"success":false,"worker_success":false});
        validate_completed_outputs(&failed, Path::new("."), root).unwrap();
    }

    #[test]
    fn completed_output_sidecars_are_rechecked_before_resume() {
        let workload = super::super::workload::tests::selected();
        let root = &workload.held_out[0];
        let stamp = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        let path =
            std::env::temp_dir().join(format!("owner-proof-output-{}-{stamp}", std::process::id()));
        fs::create_dir(&path).unwrap();
        fs::write(path.join("value.symbolica"), b"complete").unwrap();
        fs::write(path.join("guards.json"), b"[]").unwrap();
        let receipt = json!({"success":true,"worker_success":true,"fresh_process":true,
            "input_sha256":"a".repeat(64),"supervision":{"success":true,"exit_code":0,"deadline_exhausted":false},
            "value_sha256":digest(&path.join("value.symbolica")).unwrap(),"guards_sha256":digest(&path.join("guards.json")).unwrap()});
        validate_completed_outputs(&receipt, &path, root).unwrap();
        fs::write(path.join("value.symbolica"), b"changed").unwrap();
        assert!(validate_completed_outputs(&receipt, &path, root).is_err());
        fs::remove_dir_all(path).unwrap();
    }

    #[test]
    fn pinned_resume_rechecks_receipt_chain_without_archived_binaries() {
        let mut workload = super::super::workload::tests::selected();
        let stamp = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        let base =
            std::env::temp_dir().join(format!("owner-proof-chain-{}-{stamp}", std::process::id()));
        let artifacts = base.join("data");
        fs::create_dir_all(artifacts.join("5L/ibps")).unwrap();
        for path in [
            "manifest.json",
            "5L/manifest.json",
            "5L/ibps/application.json",
        ] {
            fs::write(artifacts.join(path), path.as_bytes()).unwrap();
        }
        workload.data_manifest_sha256 = digest(&artifacts.join("manifest.json")).unwrap();
        let (mut document, binding) = document(&workload);
        for (path, field) in [
            ("5L/manifest.json", "loop_manifest_sha256"),
            ("5L/ibps/application.json", "application_sha256"),
        ] {
            let hash = digest(&artifacts.join(path)).unwrap();
            document[field] = json!(hash);
            document["compiled_cohort"][field] = json!(hash);
        }
        fs::write(
            base.join("cohort.json"),
            serde_json::to_vec(&document["compiled_cohort"]).unwrap(),
        )
        .unwrap();
        document["cohort_receipt_sha256"] = json!(digest(&base.join("cohort.json")).unwrap());
        let declared = |root: &Root| json!({"id":root.id,"powers":root.powers,"public_input":root.public_input});
        let suite = json!({"tuning":workload.tuning.iter().map(declared).collect::<Vec<_>>(),
            "held_out":workload.held_out.iter().map(declared).collect::<Vec<_>>()});
        for (name, path_field, hash_field, value) in [
            (
                "suite.json",
                "frozen_suite_path",
                "frozen_suite_sha256",
                suite,
            ),
            (
                "reserves.json",
                "reserve_path",
                "reserve_sha256",
                json!({"held_out_reserves":[]}),
            ),
            (
                "protocol.json",
                "protocol_path",
                "protocol_sha256",
                json!({"fixture":true}),
            ),
        ] {
            fs::write(base.join(name), serde_json::to_vec(&value).unwrap()).unwrap();
            document[path_field] = json!(name);
            document[hash_field] = json!(digest(&base.join(name)).unwrap());
        }
        fs::write(base.join("source-census.json"), b"source-only fixture").unwrap();
        document["mapper_executable_path"] = json!("absent-mapper-ELF");
        document["measurement_executable_path"] = json!("absent-measurement-ELF");
        for (index, row) in document["rows"]
            .as_array_mut()
            .unwrap()
            .iter_mut()
            .enumerate()
        {
            let directory = base.join(format!("case-{index}"));
            fs::create_dir(&directory).unwrap();
            fs::write(directory.join("value.symbolica"), b"complete").unwrap();
            fs::write(directory.join("guards.json"), b"[]").unwrap();
            let root = workload
                .tuning
                .iter()
                .chain(&workload.held_out)
                .find(|root| row["id"] == root.id)
                .unwrap();
            let source_workload = json!({"tuning":if row["logical_cohort"]=="tuning" {vec![declared(root)]}else{vec![]},
                "held_out":if row["logical_cohort"]=="held-out" {vec![declared(root)]}else{vec![]}});
            fs::write(
                directory.join("workload.json"),
                serde_json::to_vec(&source_workload).unwrap(),
            )
            .unwrap();
            row["workload_path"] = json!(format!("case-{index}/workload.json"));
            row["workload_sha256"] = json!(digest(&directory.join("workload.json")).unwrap());
            row["source_census_path"] = json!("source-census.json");
            row["source_census_sha256"] = json!(digest(&base.join("source-census.json")).unwrap());
            let receipt = json!({"case_id":row["id"],"variant":"uncached","sample":0,
                "success":true,"worker_success":true,"fresh_process":true,"input_sha256":"b".repeat(64),
                "actual_lowered_roots":row["actual_lowered_roots"],"cache_statistics":{"binding":binding},
                "identity":{"data_manifest_sha256":workload.data_manifest_sha256,"executable_sha256":"a".repeat(64),"cache_sha256":null,
                    "workload_sha256":row["workload_sha256"],"source_census_sha256":row["source_census_sha256"]},
                "supervision":{"success":true,"exit_code":0,"deadline_exhausted":false},
                "value_sha256":digest(&directory.join("value.symbolica")).unwrap(),"guards_sha256":digest(&directory.join("guards.json")).unwrap()});
            fs::write(
                directory.join("measurement.json"),
                serde_json::to_vec(&receipt).unwrap(),
            )
            .unwrap();
            row["input_sha256"] = receipt["input_sha256"].clone();
            row["receipt_path"] = json!(format!("case-{index}/measurement.json"));
            row["receipt_sha256"] = json!(digest(&directory.join("measurement.json")).unwrap());
        }
        let proof_path = base.join("proof.json");
        fs::write(&proof_path, serde_json::to_vec(&document).unwrap()).unwrap();
        let reference = workload.owner_proof.as_mut().unwrap();
        reference.path = proof_path.clone();
        reference.sha256 = digest(&proof_path).unwrap();
        pinned_identity(&workload, &base.join("workload.json"), &artifacts).unwrap();
        fs::write(
            base.join("case-0/value.symbolica"),
            b"changed after successful receipt",
        )
        .unwrap();
        assert!(pinned_identity(&workload, &base.join("workload.json"), &artifacts).is_err());
        fs::write(base.join("case-0/value.symbolica"), b"complete").unwrap();
        fs::write(base.join("source-census.json"), b"changed source evidence").unwrap();
        assert!(pinned_identity(&workload, &base.join("workload.json"), &artifacts).is_err());
        fs::write(base.join("source-census.json"), b"source-only fixture").unwrap();
        let source_path = base.join("case-0/workload.json");
        let mut wrong_role: Value =
            serde_json::from_slice(&fs::read(&source_path).unwrap()).unwrap();
        wrong_role["held_out"] = wrong_role["tuning"].clone();
        wrong_role["tuning"] = json!([]);
        fs::write(&source_path, serde_json::to_vec(&wrong_role).unwrap()).unwrap();
        let source_hash = digest(&source_path).unwrap();
        document["rows"][0]["workload_sha256"] = json!(source_hash);
        let receipt_path = base.join("case-0/measurement.json");
        let mut wrong_receipt: Value =
            serde_json::from_slice(&fs::read(&receipt_path).unwrap()).unwrap();
        wrong_receipt["identity"]["workload_sha256"] = json!(source_hash);
        fs::write(&receipt_path, serde_json::to_vec(&wrong_receipt).unwrap()).unwrap();
        document["rows"][0]["receipt_sha256"] = json!(digest(&receipt_path).unwrap());
        fs::write(&proof_path, serde_json::to_vec(&document).unwrap()).unwrap();
        workload.owner_proof.as_mut().unwrap().sha256 = digest(&proof_path).unwrap();
        assert!(pinned_identity(&workload, &base.join("workload.json"), &artifacts).is_err());
        fs::remove_dir_all(base).unwrap();
    }

    #[test]
    fn requires_complete_matching_sources_and_loaded_binding() {
        let workload = super::super::workload::tests::selected();
        let (original, binding) = document(&workload);
        assert!(validate(original.clone(), &workload, &binding, Value::Null).is_ok());
        for (field, value) in [
            ("reduction_succeeded", json!(false)),
            ("complete_owner_transport", json!(false)),
            ("public_input_source_sha256", json!("wrong")),
        ] {
            let mut changed = original.clone();
            changed["rows"][0][field] = value;
            assert!(validate(changed, &workload, &binding, Value::Null).is_err());
        }
        let mut wrong = binding.clone();
        wrong["artifact_digest"][0] = json!(2);
        assert!(validate(original, &workload, &wrong, Value::Null).is_err());
    }
    #[test]
    fn rejects_missing_duplicate_or_inconsistent_native_transports() {
        let workload = super::super::workload::tests::selected();
        let (original, binding) = document(&workload);
        let mut missing = original.clone();
        missing["rows"][0]["native_transports"] = json!([]);
        assert!(validate(missing, &workload, &binding, Value::Null).is_err());
        let mut duplicate = original.clone();
        let row = duplicate["rows"][0].clone();
        duplicate["rows"].as_array_mut().unwrap().push(row);
        assert!(validate(duplicate, &workload, &binding, Value::Null).is_err());
        let mut changed = original;
        changed["rows"][0]["native_transports"][0]["owner_endpoints"][0]["powers"][0] = json!(100);
        assert!(validate(changed, &workload, &binding, Value::Null).is_err());
    }
    #[test]
    fn aliases_do_not_count_as_twelve_and_correlations_are_reported() {
        let mut workload = super::super::workload::tests::selected();
        workload.held_out[11].protected_roots = workload.held_out[0].protected_roots.clone();
        let (aliased, binding) = document(&workload);
        assert!(validate(aliased, &workload, &binding, Value::Null).is_err());
        let mut workload = super::super::workload::tests::selected();
        workload.tuning[0].protected_roots = workload.held_out[0].protected_roots.clone();
        let (correlated, binding) = document(&workload);
        let proof = validate(correlated, &workload, &binding, Value::Null).unwrap();
        assert_eq!(proof.correlations().as_array().unwrap().len(), 1);
    }
}
