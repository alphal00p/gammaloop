//! Lossless native event reconciliation, independent of workload authority.
use eyre::{Result, ensure};
use serde_json::{Value, json};
use std::collections::BTreeSet;

fn string<'a>(value: &'a Value, field: &str) -> Result<&'a str> {
    value[field]
        .as_str()
        .ok_or_else(|| eyre::eyre!("missing owner-proof {field}"))
}

fn array<'a>(value: &'a Value, field: &str) -> Result<&'a Vec<Value>> {
    value[field]
        .as_array()
        .ok_or_else(|| eyre::eyre!("missing owner-proof {field}"))
}

fn count(value: &Value, field: &str) -> Result<u64> {
    value[field]
        .as_u64()
        .ok_or_else(|| eyre::eyre!("missing nonnegative counter {field}"))
}

/// Full native work identity. A Route and an Apply of the same powers are
/// different cache records; an Apply also carries its literal owner support.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub(crate) struct FrontierWork {
    pub(crate) stage: String,
    pub(crate) powers: Vec<i64>,
    pub(crate) owner: Option<Vec<bool>>,
}

impl FrontierWork {
    pub(crate) fn parse(value: &Value) -> Result<Self> {
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

    pub(crate) fn support(powers: &[i64]) -> Vec<bool> {
        powers.iter().map(|power| *power > 0).collect()
    }

    pub(crate) fn document(&self) -> Value {
        json!({"stage":self.stage,"powers":self.powers,"owner_sector":self.owner})
    }
}

/// Reconcile a lossless diagnostic with native per-dispatch counters. A latest
/// monitor snapshot or a root union is deliberately insufficient evidence.
pub(crate) fn validate_hit_diagnostic(document: &Value) -> Result<Value> {
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
