//! Scalar observations of the actual traversal, with disjoint self accounting.
//! No coefficient, guard, cache record, application handle or event JSON is retained.
use rustred::solver::{
    CandidateRoutedActivity as Activity, CandidateRoutedApplicationEnvelope as Envelope,
    CandidateRoutedApplicationEvent as Event, CandidateRoutedApplicationPhase as Phase,
    CandidateRoutedApplicationStatus as Status, CandidateRoutedCompletionSource as Source,
    CandidateRoutedObserverControl as Control, CandidateRoutedObserverStatus as ObserverStatus,
    CandidateRoutedWorkRef as Work,
};
use serde::Serialize;
use serde_json::{Value, json};
use std::collections::HashMap;

const COUNTER_NAMES: [&str; 19] = [
    "rule_applications",
    "coalescing_additions",
    "cache_hits",
    "transport_calls",
    "transport_operations",
    "transport_endpoints",
    "persistent_lookups",
    "persistent_route_hits",
    "persistent_apply_hits",
    "persistent_requested_root_hits",
    "persistent_intermediate_hits",
    "persistent_misses",
    "captured_records",
    "skipped_records",
    "recomputed_nodes",
    "evicted_entries",
    "evicted_coefficient_terms",
    "evicted_coefficient_bytes",
    "observed_descent_checks",
];
type DiagnosticResult<T> = Result<T, &'static str>;

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Serialize)]
struct Counters([u64; COUNTER_NAMES.len()]);

impl Counters {
    fn read(event: &Envelope<'_>) -> DiagnosticResult<Self> {
        let snapshot = event.snapshot.ok_or("native event omitted snapshot")?;
        let s = snapshot.statistics;
        let p = s.persistent_cache;
        let values = [
            s.work.rule_applications(),
            s.work.coalescing_additions(),
            s.work.cache_hits(),
            s.transport_calls,
            s.transport_operations,
            s.transport_endpoints,
            p.lookups,
            p.route_hits,
            p.apply_hits,
            p.requested_root_hits,
            p.intermediate_hits,
            p.misses,
            p.captured_records,
            p.skipped_records,
            s.cache.recomputed_nodes,
            s.cache.evicted_entries,
            s.cache.evicted_coefficient_terms,
            s.cache.evicted_coefficient_bytes,
            event.metrics.observed_descent_checks,
        ];
        let mut result = Self::default();
        for (target, value) in result.0.iter_mut().zip(values) {
            // Native cumulative counters saturate. Once at MAX their differences
            // cannot authenticate all work, even if no later sample changes.
            if value == usize::MAX {
                return Err("native counter saturated");
            }
            *target = u64::try_from(value).map_err(|_| "native counter exceeds u64")?;
        }
        Ok(result)
    }

    fn difference(self, before: Self) -> DiagnosticResult<Self> {
        let mut result = Self::default();
        for ((target, after), before) in result.0.iter_mut().zip(self.0).zip(before.0) {
            *target = after
                .checked_sub(before)
                .ok_or("native counter decreased")?;
        }
        Ok(result)
    }

    fn add(&mut self, other: Self) -> DiagnosticResult<()> {
        for (target, value) in self.0.iter_mut().zip(other.0) {
            *target = target
                .checked_add(value)
                .ok_or("diagnostic counter overflow")?;
        }
        Ok(())
    }
}

#[derive(Debug, PartialEq, Eq, Hash, Serialize)]
struct Identity {
    family_id: usize,
    stage: &'static str,
    powers: Vec<i64>,
    owner_sector: Option<Vec<bool>>,
}

impl Identity {
    fn new(family_id: usize, work: Work<'_>) -> DiagnosticResult<Self> {
        let mut powers = Vec::new();
        powers
            .try_reserve(work.target().powers().len())
            .map_err(|_| "key allocation failed")?;
        powers.extend_from_slice(work.target().powers());
        let (stage, owner_sector) = match work {
            Work::Route(_) => ("route", None),
            Work::Apply { owner_sector, .. } => {
                let mut owner = Vec::new();
                owner
                    .try_reserve(owner_sector.len())
                    .map_err(|_| "owner allocation failed")?;
                owner.extend_from_slice(owner_sector);
                ("apply", Some(owner))
            }
        };
        Ok(Self {
            family_id,
            stage,
            powers,
            owner_sector,
        })
    }
}

#[derive(Default, Serialize)]
struct KeySummary {
    visits: usize,
    completions_by_source: [usize; 6],
    #[serde(skip)]
    last_completed_dispatch: Option<usize>,
}

#[derive(Serialize)]
struct Attempt {
    id: usize,
    dispatch_id: usize,
    key_id: usize,
    parent_id: Option<usize>,
    enter_sequence: u64,
    completion_sequence: Option<u64>,
    enter_elapsed_ns: u64,
    completion_elapsed_ns: Option<u64>,
    source: Option<&'static str>,
    start_counters: Counters,
    end_counters: Option<Counters>,
    inclusive_counters: Option<Counters>,
    self_counters: Counters,
    direct_child_inclusive_counters: Counters,
    inclusive_wall_ns: Option<u64>,
    self_wall_ns: u64,
    direct_child_inclusive_wall_ns: u64,
    reconciled: bool,
}

#[derive(Serialize)]
struct Dispatch {
    id: usize,
    family_id: usize,
    start_sequence: u64,
    end_sequence: Option<u64>,
    start_elapsed_ns: u64,
    end_elapsed_ns: Option<u64>,
    start_counters: Counters,
    total_counters: Option<Counters>,
    attributed_self_counters: Counters,
    completed_attempt_self_sum: Counters,
    unscoped_counters: Counters,
    total_wall_ns: Option<u64>,
    attributed_self_wall_ns: u64,
    completed_attempt_self_wall_sum: u64,
    unscoped_wall_ns: u64,
    entered_attempts: usize,
    completed_attempts: usize,
    completed_unique_keys: usize,
    native_completion_attempts: Option<usize>,
    native_completed_unique_nodes: Option<usize>,
    succeeded: bool,
    reconciled: bool,
    #[serde(skip)]
    last_counters: Counters,
    #[serde(skip)]
    last_elapsed_ns: u64,
    #[serde(skip)]
    last_native_sequence: u64,
    #[serde(skip)]
    family_pointer: usize,
    #[serde(skip)]
    family_length: usize,
}

pub(super) struct CostObservation {
    binding: Value,
    families: HashMap<String, usize>,
    identities: HashMap<Identity, usize>,
    keys: Vec<KeySummary>,
    attempts: Vec<Attempt>,
    dispatches: Vec<Dispatch>,
    active_dispatch: Option<usize>,
    stack: Vec<usize>,
    last_sequence: Option<u64>,
    fault: Option<&'static str>,
    native_events: usize,
    max_open_attempts: usize,
}

impl CostObservation {
    pub(super) fn new(binding: Value) -> Self {
        Self {
            binding,
            families: HashMap::new(),
            identities: HashMap::new(),
            keys: Vec::new(),
            attempts: Vec::new(),
            dispatches: Vec::new(),
            active_dispatch: None,
            stack: Vec::new(),
            last_sequence: None,
            fault: None,
            native_events: 0,
            max_open_attempts: 0,
        }
    }

    pub(super) fn observe(&mut self, event: &Envelope<'_>) -> Result<Control, ()> {
        if self.fault.is_some() {
            return Err(());
        }
        if let Err(error) = self.observe_native(event) {
            self.fault = Some(error);
            return Err(());
        }
        Ok(Control::Continue)
    }

    fn observe_native(&mut self, event: &Envelope<'_>) -> DiagnosticResult<()> {
        // Saved/caller envelopes may forward an old native snapshot. They are
        // outside core cost accounting even when their family matches.
        if event.activity == Activity::Caller || event.caller_phase.is_some() {
            return Ok(());
        }
        let snapshot = event.snapshot.ok_or("native event omitted snapshot")?;
        if event.sequence == 0
            || self
                .last_sequence
                .is_some_and(|last| event.sequence <= last)
        {
            return Err("observer sequence did not increase");
        }
        self.last_sequence = Some(event.sequence);
        self.native_events = self
            .native_events
            .checked_add(1)
            .ok_or("event count overflow")?;
        let now = u64::try_from(event.elapsed.as_nanos())
            .map_err(|_| "wall clock exceeds u64 nanoseconds")?;
        let counters = Counters::read(event)?;
        let initial = event.event == Event::Snapshot
            && snapshot.sequence == 1
            && snapshot.phase == Phase::Admission;
        if initial {
            return self.start_dispatch(event, now, counters);
        }
        let dispatch_id = self
            .active_dispatch
            .ok_or("native event without dispatch start")?;
        let family_id = self.dispatches[dispatch_id].family_id;
        let dispatch = &mut self.dispatches[dispatch_id];
        let family_pointer = event.family_fingerprint.as_ptr() as usize;
        if family_pointer != dispatch.family_pointer
            || event.family_fingerprint.len() != dispatch.family_length
        {
            // Native normally borrows the same immutable family string for the
            // whole request. Only an address change needs a full-content lookup.
            if self.families.get(event.family_fingerprint) != Some(&family_id) {
                return Err("family changed inside dispatch");
            }
            dispatch.family_pointer = family_pointer;
            dispatch.family_length = event.family_fingerprint.len();
        }
        if self.dispatches[dispatch_id]
            .last_native_sequence
            .checked_add(1)
            != Some(snapshot.sequence)
        {
            return Err("native snapshot sequence was not contiguous");
        }
        self.dispatches[dispatch_id].last_native_sequence = snapshot.sequence;
        self.charge(dispatch_id, counters, now)?;
        match event.event {
            Event::Snapshot if snapshot.phase == Phase::Expansion => {
                let work = event.current.ok_or("enter omitted current work")?;
                self.enter(dispatch_id, family_id, work, event.sequence, now, counters)?;
            }
            Event::NodeCompleted(source) => {
                let work = event.current.ok_or("completion omitted current work")?;
                self.complete(dispatch_id, work, source, event.sequence, now, counters)?;
            }
            Event::CoreFinished => self.end_dispatch(dispatch_id, event, now, counters)?,
            _ => {}
        }
        Ok(())
    }

    fn start_dispatch(
        &mut self,
        event: &Envelope<'_>,
        now: u64,
        counters: Counters,
    ) -> DiagnosticResult<()> {
        if self.active_dispatch.is_some() || !self.stack.is_empty() {
            return Err("new dispatch before previous core finished");
        }
        if event.family_fingerprint.is_empty() || event.status != Status::Running {
            return Err("invalid dispatch start");
        }
        let family_id = if let Some(id) = self.families.get(event.family_fingerprint) {
            *id
        } else {
            self.families
                .try_reserve(1)
                .map_err(|_| "family table allocation failed")?;
            let mut family = String::new();
            family
                .try_reserve(event.family_fingerprint.len())
                .map_err(|_| "family allocation failed")?;
            family.push_str(event.family_fingerprint);
            let id = self.families.len();
            self.families.insert(family, id);
            id
        };
        self.dispatches
            .try_reserve(1)
            .map_err(|_| "dispatch allocation failed")?;
        let id = self.dispatches.len();
        self.dispatches.push(Dispatch {
            id,
            family_id,
            start_sequence: event.sequence,
            end_sequence: None,
            start_elapsed_ns: now,
            end_elapsed_ns: None,
            start_counters: counters,
            total_counters: None,
            attributed_self_counters: Counters::default(),
            completed_attempt_self_sum: Counters::default(),
            unscoped_counters: Counters::default(),
            total_wall_ns: None,
            attributed_self_wall_ns: 0,
            completed_attempt_self_wall_sum: 0,
            unscoped_wall_ns: 0,
            entered_attempts: 0,
            completed_attempts: 0,
            completed_unique_keys: 0,
            native_completion_attempts: None,
            native_completed_unique_nodes: None,
            succeeded: false,
            reconciled: false,
            last_counters: counters,
            last_elapsed_ns: now,
            last_native_sequence: 1,
            family_pointer: event.family_fingerprint.as_ptr() as usize,
            family_length: event.family_fingerprint.len(),
        });
        self.active_dispatch = Some(id);
        Ok(())
    }

    fn charge(&mut self, id: usize, counters: Counters, now: u64) -> DiagnosticResult<()> {
        let dispatch = &mut self.dispatches[id];
        let delta = counters.difference(dispatch.last_counters)?;
        let wall = now
            .checked_sub(dispatch.last_elapsed_ns)
            .ok_or("wall clock decreased")?;
        if let Some(attempt) = self.stack.last().copied() {
            self.attempts[attempt].self_counters.add(delta)?;
            self.attempts[attempt].self_wall_ns = self.attempts[attempt]
                .self_wall_ns
                .checked_add(wall)
                .ok_or("attempt wall overflow")?;
            dispatch.attributed_self_counters.add(delta)?;
            dispatch.attributed_self_wall_ns = dispatch
                .attributed_self_wall_ns
                .checked_add(wall)
                .ok_or("dispatch wall overflow")?;
        } else {
            dispatch.unscoped_counters.add(delta)?;
            dispatch.unscoped_wall_ns = dispatch
                .unscoped_wall_ns
                .checked_add(wall)
                .ok_or("unscoped wall overflow")?;
        }
        dispatch.last_counters = counters;
        dispatch.last_elapsed_ns = now;
        Ok(())
    }

    fn enter(
        &mut self,
        dispatch_id: usize,
        family_id: usize,
        work: Work<'_>,
        sequence: u64,
        now: u64,
        counters: Counters,
    ) -> DiagnosticResult<()> {
        let identity = Identity::new(family_id, work)?;
        let key_id = if let Some(id) = self.identities.get(&identity) {
            *id
        } else {
            self.identities
                .try_reserve(1)
                .map_err(|_| "key table allocation failed")?;
            self.keys
                .try_reserve(1)
                .map_err(|_| "key summary allocation failed")?;
            let id = self.keys.len();
            self.identities.insert(identity, id);
            self.keys.push(KeySummary::default());
            id
        };
        self.attempts
            .try_reserve(1)
            .map_err(|_| "attempt allocation failed")?;
        self.stack
            .try_reserve(1)
            .map_err(|_| "stack allocation failed")?;
        let id = self.attempts.len();
        self.attempts.push(Attempt {
            id,
            dispatch_id,
            key_id,
            parent_id: self.stack.last().copied(),
            enter_sequence: sequence,
            completion_sequence: None,
            enter_elapsed_ns: now,
            completion_elapsed_ns: None,
            source: None,
            start_counters: counters,
            end_counters: None,
            inclusive_counters: None,
            self_counters: Counters::default(),
            direct_child_inclusive_counters: Counters::default(),
            inclusive_wall_ns: None,
            self_wall_ns: 0,
            direct_child_inclusive_wall_ns: 0,
            reconciled: false,
        });
        self.stack.push(id);
        self.max_open_attempts = self.max_open_attempts.max(self.stack.len());
        self.keys[key_id].visits = self.keys[key_id]
            .visits
            .checked_add(1)
            .ok_or("visit count overflow")?;
        self.dispatches[dispatch_id].entered_attempts += 1;
        Ok(())
    }

    fn complete(
        &mut self,
        dispatch_id: usize,
        work: Work<'_>,
        source: Source,
        sequence: u64,
        now: u64,
        counters: Counters,
    ) -> DiagnosticResult<()> {
        let family_id = self.dispatches[dispatch_id].family_id;
        let id = *self
            .stack
            .last()
            .ok_or("completion without an open enter")?;
        let identity = Identity::new(family_id, work)?;
        let key_id = *self
            .identities
            .get(&identity)
            .ok_or("completion key was never entered")?;
        let attempt = &mut self.attempts[id];
        if attempt.key_id != key_id || attempt.dispatch_id != dispatch_id {
            return Err("completion differs from active attempt");
        }
        attempt.completion_sequence = Some(sequence);
        attempt.completion_elapsed_ns = Some(now);
        attempt.end_counters = Some(counters);
        attempt.inclusive_counters = Some(counters.difference(attempt.start_counters)?);
        attempt.inclusive_wall_ns = Some(
            now.checked_sub(attempt.enter_elapsed_ns)
                .ok_or("attempt wall decreased")?,
        );
        let (source_index, source_name) = match source {
            Source::Computed => (0, "computed"),
            Source::RamHit => (1, "ram_hit"),
            Source::BackingHit => (2, "backing_hit"),
            Source::Zero => (3, "zero"),
            Source::Terminal => (4, "terminal"),
            Source::UnresolvedPreserved => (5, "unresolved_preserved"),
        };
        attempt.source = Some(source_name);
        let inclusive = attempt.inclusive_counters.unwrap();
        let inclusive_wall = attempt.inclusive_wall_ns.unwrap();
        let mut sum = attempt.self_counters;
        sum.add(attempt.direct_child_inclusive_counters)?;
        attempt.reconciled = sum == inclusive
            && attempt
                .self_wall_ns
                .checked_add(attempt.direct_child_inclusive_wall_ns)
                == Some(inclusive_wall);
        if !attempt.reconciled {
            return Err("attempt self and child accounting did not reconcile");
        }
        self.stack.pop();
        if let Some(parent_id) = self.stack.last().copied() {
            let parent = &mut self.attempts[parent_id];
            parent.direct_child_inclusive_counters.add(inclusive)?;
            parent.direct_child_inclusive_wall_ns = parent
                .direct_child_inclusive_wall_ns
                .checked_add(inclusive_wall)
                .ok_or("child wall overflow")?;
        }
        let summary = &mut self.keys[key_id];
        summary.completions_by_source[source_index] = summary.completions_by_source[source_index]
            .checked_add(1)
            .ok_or("source count overflow")?;
        let dispatch = &mut self.dispatches[dispatch_id];
        dispatch.completed_attempts += 1;
        if summary.last_completed_dispatch != Some(dispatch_id) {
            summary.last_completed_dispatch = Some(dispatch_id);
            dispatch.completed_unique_keys += 1;
        }
        Ok(())
    }

    fn end_dispatch(
        &mut self,
        id: usize,
        event: &Envelope<'_>,
        now: u64,
        counters: Counters,
    ) -> DiagnosticResult<()> {
        let dispatch = &mut self.dispatches[id];
        dispatch.end_sequence = Some(event.sequence);
        dispatch.end_elapsed_ns = Some(now);
        dispatch.total_counters = Some(counters.difference(dispatch.start_counters)?);
        dispatch.total_wall_ns = Some(
            now.checked_sub(dispatch.start_elapsed_ns)
                .ok_or("dispatch wall decreased")?,
        );
        dispatch.native_completion_attempts = Some(event.metrics.completion_attempts);
        dispatch.native_completed_unique_nodes =
            Some(event.metrics.completed_unique_operational_nodes);
        dispatch.succeeded = event.status == Status::Succeeded;
        let mut reconciled = dispatch.attributed_self_counters;
        reconciled.add(dispatch.unscoped_counters)?;
        dispatch.reconciled = Some(reconciled) == dispatch.total_counters
            && dispatch
                .attributed_self_wall_ns
                .checked_add(dispatch.unscoped_wall_ns)
                == dispatch.total_wall_ns;
        if !self.stack.is_empty() {
            return Err("core finished with uncompleted attempts");
        }
        if dispatch.entered_attempts != dispatch.completed_attempts
            || dispatch.completed_attempts != event.metrics.completion_attempts
            || dispatch.completed_unique_keys != event.metrics.completed_unique_operational_nodes
        {
            return Err("native completion census differs from observed attempts");
        }
        if !dispatch.reconciled {
            return Err("self and unscoped accounting did not reconcile");
        }
        self.active_dispatch = None;
        Ok(())
    }

    fn reconcile_completed_attempts(&mut self) -> DiagnosticResult<()> {
        // Re-sum retained completed attempts independently of the callback's
        // dispatch accumulators. Inclusive spans deliberately are not summed.
        for attempt in self
            .attempts
            .iter()
            .filter(|attempt| attempt.completion_sequence.is_some())
        {
            let dispatch = &mut self.dispatches[attempt.dispatch_id];
            dispatch
                .completed_attempt_self_sum
                .add(attempt.self_counters)?;
            dispatch.completed_attempt_self_wall_sum = dispatch
                .completed_attempt_self_wall_sum
                .checked_add(attempt.self_wall_ns)
                .ok_or("completed attempt wall sum overflow")?;
        }
        for dispatch in &mut self.dispatches {
            dispatch.reconciled &= dispatch.completed_attempt_self_sum
                == dispatch.attributed_self_counters
                && dispatch.completed_attempt_self_wall_sum == dispatch.attributed_self_wall_ns;
            if !dispatch.reconciled {
                return Err("completed attempt self sum did not reconcile");
            }
        }
        Ok(())
    }

    pub(super) fn finish(
        mut self,
        operation_succeeded: bool,
        observer_status: ObserverStatus,
    ) -> Value {
        if self.fault.is_none() && (self.active_dispatch.is_some() || !self.stack.is_empty()) {
            self.fault = Some("observation ended before core completion");
        }
        if let Err(error) = self.reconcile_completed_attempts() {
            self.fault.get_or_insert(error);
        }
        let complete = operation_succeeded
            && observer_status == ObserverStatus::Active
            && self.fault.is_none()
            && !self.dispatches.is_empty()
            && self
                .dispatches
                .iter()
                .all(|d| d.succeeded && d.reconciled && d.end_sequence.is_some());
        let family_utf8_bytes: usize = self.families.keys().map(String::len).sum();
        let key_power_slots: usize = self.identities.keys().map(|key| key.powers.len()).sum();
        let key_owner_slots: usize = self
            .identities
            .keys()
            .filter_map(|key| key.owner_sector.as_ref())
            .map(Vec::len)
            .sum();
        let family_table_capacity = self.families.capacity();
        let key_table_capacity = self.identities.capacity();
        let mut families: Vec<_> = self.families.into_iter().collect();
        families.sort_unstable_by_key(|(_, id)| *id);
        let family_count = families.len();
        let mut keys: Vec<_> = self.identities.into_iter().collect();
        keys.sort_unstable_by_key(|(_, id)| *id);
        let keys: Vec<_> = keys
            .into_iter()
            .map(|(work, id)| {
                json!({"id": id, "work": work, "visits": self.keys[id].visits,
                "completions_by_source": self.keys[id].completions_by_source})
            })
            .collect();
        json!({
            "schema": "vakint.native-cost-observation.v1", "binding": self.binding,
            "complete": complete, "operation_succeeded": operation_succeeded,
            "observer_status": format!("{observer_status:?}"), "error": self.fault,
            "interval_ns": 0, "max_frontier_frames": 0, "max_path_nodes": 0,
            "counter_names": COUNTER_NAMES,
            "completion_sources": ["computed", "ram_hit", "backing_hit", "zero", "terminal", "unresolved_preserved"],
            "families": families.into_iter().map(|(fingerprint, id)| json!({"id": id, "fingerprint": fingerprint})).collect::<Vec<_>>(),
            "metadata": {"native_events": self.native_events, "families": family_count,
                "keys": keys.len(), "attempts": self.attempts.len(), "dispatches": self.dispatches.len(),
                "max_open_attempts": self.max_open_attempts, "uncompleted_attempts": self.stack.len(),
                "attempt_capacity": self.attempts.capacity(), "stack_capacity": self.stack.capacity(),
                "family_table_capacity": family_table_capacity, "key_table_capacity": key_table_capacity,
                "family_utf8_bytes": family_utf8_bytes, "key_power_slots": key_power_slots,
                "key_owner_slots": key_owner_slots},
            "keys": keys, "attempts": self.attempts, "dispatches": self.dispatches,
            "scope": "actual observed core traversal; additive counters and disjoint self intervals; wall spans include observation overhead; no counterfactual saved-work claim",
            "eligibility_changed": false, "performance_claim": false,
        })
    }
}

#[cfg(test)]
#[path = "cost_observation_tests.rs"]
mod tests;
