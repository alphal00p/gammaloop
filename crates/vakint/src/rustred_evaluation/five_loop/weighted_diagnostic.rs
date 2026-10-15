//! Ignored same-build physical witness gate. All hooks are passive libtest-only.
use super::*;
use rustred::scalar_numerator::{
    LoweredScalarNumeratorTerm, WeightedScalarNumeratorLowering,
    lowering_diagnostics::{self, Counts},
};
use std::{cell::RefCell, sync::Arc, thread::ThreadId, time::Instant};

const MAX_RECORDS: usize = 100_000;
const MAX_PAYLOAD_BYTES: usize = 2 * 1024 * 1024 * 1024;
type Key = Vec<(IntegralKey, Atom)>;

#[derive(Clone)]
enum Record {
    Begin {
        integral: Atom,
        numerator: Atom,
    },
    Candidate {
        search: bool,
        witness: ReplacementRules,
        numerator: Atom,
        target: IntegralKey,
        mass: Atom,
        family: String,
        terms: Vec<LoweredScalarNumeratorTerm>,
        conditions: Vec<Atom>,
        rejected: Option<usize>,
        key: Option<Key>,
    },
    End {
        selected: Option<ReplacementRules>,
    },
}
fn same_witness(a: &ReplacementRules, b: &ReplacementRules) -> bool {
    let left = a.canonical_topology.get_integral();
    let right = b.canonical_topology.get_integral();
    a.edge_ids_canonical_to_input_map == b.edge_ids_canonical_to_input_map
        && a.canonical_expression_substitutions == b.canonical_expression_substitutions
        && a.numerator_substitutions == b.numerator_substitutions
        && left.canonical_expression == right.canonical_expression
        && left.parent_routing == right.parent_routing
        && left.short_expression == right.short_expression
}
fn same_optional_witness(a: Option<&ReplacementRules>, b: Option<&ReplacementRules>) -> bool {
    match (a, b) {
        (Some(a), Some(b)) => same_witness(a, b),
        (None, None) => true,
        _ => false,
    }
}
impl Record {
    fn agrees(&self, other: &Self) -> bool {
        match (self, other) {
            (
                Self::Begin {
                    integral: a,
                    numerator: b,
                },
                Self::Begin {
                    integral: c,
                    numerator: d,
                },
            ) => a == c && b == d,
            (Self::End { selected: a }, Self::End { selected: b }) => {
                same_optional_witness(a.as_ref(), b.as_ref())
            }
            (
                Self::Candidate {
                    search: a,
                    witness: aw,
                    numerator: an,
                    target: at,
                    mass: am,
                    family: af,
                    terms: av,
                    conditions: ac,
                    rejected: ar,
                    key: ak,
                },
                Self::Candidate {
                    search: b,
                    witness: bw,
                    numerator: bn,
                    target: bt,
                    mass: bm,
                    family: bf,
                    terms: bv,
                    conditions: bc,
                    rejected: br,
                    key: bk,
                },
            ) => {
                a == b
                    && same_witness(aw, bw)
                    && an == bn
                    && at == bt
                    && am == bm
                    && af == bf
                    && av == bv
                    && ac == bc
                    && ar == br
                    && ak == bk
            }
            _ => false,
        }
    }
}
struct State {
    reference: Option<Arc<Vec<Record>>>,
    records: Vec<Record>,
    cursor: usize,
    bytes: usize,
    owner: ThreadId,
    in_search: bool,
    best: Option<(Key, ReplacementRules)>,
    searches: usize,
    candidates: usize,
    rechecks: usize,
    observer_seconds: f64,
}
thread_local! {static STATE:RefCell<Option<State>>=const {RefCell::new(None)};}

fn error(message: impl Into<String>) -> VakintError {
    failure(message.into()).into()
}
fn with_state<T>(
    operation: impl FnOnce(&mut State) -> Result<T, VakintError>,
) -> Result<Option<T>, VakintError> {
    STATE.with(|state| {
        let mut state = state.borrow_mut();
        let Some(state) = state.as_mut() else {
            return Ok(None);
        };
        if state.owner != std::thread::current().id() || lowering_diagnostics::current().is_none() {
            return Err(error(
                "weighted witness diagnostic lost its synchronous scope",
            ));
        }
        let start = Instant::now();
        let result = operation(state);
        state.observer_seconds += start.elapsed().as_secs_f64();
        result.map(Some)
    })
}
impl State {
    fn complete(&self) -> Result<(), VakintError> {
        if self.in_search
            || self.searches == 0
            || self
                .reference
                .as_ref()
                .is_some_and(|rows| rows.len() != self.cursor)
        {
            return Err(error("weighted witness recording incomplete"));
        }
        Ok(())
    }
    fn consider(&mut self, key: &Key, witness: &ReplacementRules) {
        if self.best.as_ref().is_none_or(|(old, _)| key < old) {
            self.best = Some((key.clone(), witness.clone()));
        }
    }
    fn charge(&mut self, bytes: usize) -> Result<(), VakintError> {
        if self.cursor >= MAX_RECORDS {
            return Err(error("weighted witness recording count exhausted"));
        }
        let bytes = bytes
            .checked_add(std::mem::size_of::<Record>())
            .ok_or_else(|| error("weighted witness record metadata overflow"))?;
        self.bytes = self
            .bytes
            .checked_add(bytes)
            .ok_or_else(|| error("weighted witness recording byte overflow"))?;
        if self.bytes > MAX_PAYLOAD_BYTES {
            return Err(error("weighted witness recording payload exhausted"));
        }
        Ok(())
    }
    fn record(&mut self, value: Record) -> Result<(), VakintError> {
        if let Some(reference) = &self.reference {
            if !reference
                .get(self.cursor)
                .is_some_and(|expected| expected.agrees(&value))
            {
                return Err(error(format!(
                    "weighted witness exact record mismatch at{}",
                    self.cursor
                )));
            }
        } else {
            self.records.push(value);
        }
        self.cursor += 1;
        Ok(())
    }
}
pub(super) fn begin_search(
    integral: AtomView,
    numerator: AtomView,
) -> Result<Option<ThreadId>, VakintError> {
    with_state(|state| {
        if state.in_search {
            return Err(error("nested weighted witness searches"));
        }
        state.charge(
            integral
                .get_byte_size()
                .checked_add(numerator.get_byte_size())
                .ok_or_else(|| error("input byte overflow"))?,
        )?;
        state.record(Record::Begin {
            integral: integral.to_owned(),
            numerator: numerator.to_owned(),
        })?;
        state.in_search = true;
        state.best = None;
        state.searches += 1;
        Ok(state.owner)
    })
}
pub(super) fn assert_thread(expected: Option<ThreadId>) -> Result<(), VakintError> {
    if let Some(expected) = expected {
        if expected != std::thread::current().id()
            || lowering_diagnostics::current().is_none()
            || !STATE.with(|state| state.borrow().is_some())
        {
            return Err(error(
                "weighted witness candidate crossed diagnostic threads",
            ));
        }
    }
    Ok(())
}
pub(super) fn before_lowering() -> Result<Option<Counts>, VakintError> {
    with_state(|_| {
        lowering_diagnostics::current().ok_or_else(|| error("missing scalar diagnostic scope"))
    })
}
pub(super) fn candidate(
    witness: &ReplacementRules,
    input: &MatchedScalarInput,
    lowered: &WeightedScalarNumeratorLowering,
    rejected: Option<usize>,
    key: Option<&Key>,
    before: Option<Counts>,
) -> Result<(), VakintError> {
    let Some(before) = before else {
        return Ok(());
    };
    with_state(|state| {
        let after = lowering_diagnostics::current()
            .ok_or_else(|| error("missing scalar diagnostic counts"))?;
        let rational = after
            .rational
            .checked_sub(before.rational)
            .ok_or_else(|| error("diagnostic rational counter regression"))?;
        let generic = after
            .generic
            .checked_sub(before.generic)
            .ok_or_else(|| error("diagnostic generic counter regression"))?;
        if rational + generic != 1 {
            return Err(error(
                "candidate lowering did not execute exactly one conversion on the observed thread",
            ));
        }
        let mut bytes = super::witness_cache::witness_atom_bytes(witness, 0)
            .ok_or_else(|| error("witness byte overflow"))?;
        let coordinates = witness
            .canonical_topology
            .get_integral()
            .parent_routing
            .as_ref()
            .map_or(0, |(_, coordinates)| coordinates.len());
        for count in [lowered.nonzero_conditions().len(), coordinates] {
            bytes = bytes
                .checked_add(
                    count
                        .checked_mul(std::mem::size_of::<Atom>())
                        .ok_or_else(|| error("condition/coordinate metadata overflow"))?,
                )
                .ok_or_else(|| error("condition/coordinate payload overflow"))?;
        }
        for count in [
            witness.canonical_expression_substitutions.len(),
            witness.numerator_substitutions.len(),
        ] {
            bytes = bytes
                .checked_add(
                    count
                        .checked_mul(std::mem::size_of::<(Atom, Atom)>())
                        .ok_or_else(|| error("witness map metadata overflow"))?,
                )
                .ok_or_else(|| error("witness map payload overflow"))?;
        }
        bytes = bytes
            .checked_add(
                witness
                    .edge_ids_canonical_to_input_map
                    .len()
                    .checked_mul(std::mem::size_of::<(usize, usize)>())
                    .ok_or_else(|| error("witness edge map overflow"))?,
            )
            .and_then(|n| n.checked_add(lowered.family_fingerprint().len()))
            .and_then(|n| n.checked_add(8 * input.target.powers().len()))
            .ok_or_else(|| error("candidate metadata overflow"))?;
        for atom in [&input.numerator, &input.mass_squared]
            .into_iter()
            .chain(lowered.nonzero_conditions())
            .chain(lowered.terms().iter().map(|term| term.scalar_spectator()))
        {
            bytes = bytes
                .checked_add(atom.as_view().get_byte_size())
                .ok_or_else(|| error("lowered atom byte overflow"))?;
        }
        bytes = bytes
            .checked_add(lowering_diagnostics::coefficient_payload_bytes(lowered).map_err(error)?)
            .ok_or_else(|| error("coefficient byte overflow"))?;
        for term in lowered.terms() {
            bytes = bytes
                .checked_add(
                    std::mem::size_of::<LoweredScalarNumeratorTerm>()
                        + 8 * term.integral().powers().len(),
                )
                .ok_or_else(|| error("term metadata overflow"))?;
        }
        if let Some(key) = key {
            bytes = bytes
                .checked_add(
                    key.len()
                        .checked_mul(std::mem::size_of::<(IntegralKey, Atom)>())
                        .ok_or_else(|| error("accepted key tuple metadata overflow"))?,
                )
                .ok_or_else(|| error("accepted key metadata payload overflow"))?;
            for (integral, weight) in key {
                bytes = bytes
                    .checked_add(8 * integral.powers().len() + weight.as_view().get_byte_size())
                    .ok_or_else(|| error("accepted key byte overflow"))?;
            }
        }
        // Also cover the transient cloned best-key/witness and compare record;
        // allocator bucket overhead is still covered only by the outer RSS cap.
        state.charge(
            bytes
                .checked_mul(2)
                .ok_or_else(|| error("candidate retained payload overflow"))?,
        )?;
        if state.in_search {
            state.candidates += 1;
            if let Some(key) = key {
                state.consider(key, witness);
            }
        } else {
            state.rechecks += 1;
        }
        state.record(Record::Candidate {
            search: state.in_search,
            witness: witness.clone(),
            numerator: input.numerator.clone(),
            target: input.target.clone(),
            mass: input.mass_squared.clone(),
            family: lowered.family_fingerprint().to_owned(),
            terms: lowered.terms().to_vec(),
            conditions: lowered.nonzero_conditions().to_vec(),
            rejected,
            key: key.cloned(),
        })
    })?
    .ok_or_else(|| error("candidate record lost its diagnostic scope"))
}
pub(super) fn end_search(selected: Option<&ReplacementRules>) -> Result<(), VakintError> {
    with_state(|state| {
        if !state.in_search {
            return Err(error("unmatched weighted search end"));
        }
        let expected = state.best.as_ref().map(|(_, witness)| witness);
        if expected.is_none() || !same_optional_witness(expected, selected) {
            return Err(error(
                "weighted diagnostic winner differs from exact strict minimum accepted key",
            ));
        }
        state.charge(
            selected
                .map_or(Some(0), |witness| {
                    super::witness_cache::witness_atom_bytes(witness, 0)
                })
                .ok_or_else(|| error("selected witness byte overflow"))?,
        )?;
        state.record(Record::End {
            selected: selected.cloned(),
        })?;
        state.in_search = false;
        state.best = None;
        Ok(())
    })?;
    Ok(())
}

fn measured<T>(
    generic: bool,
    reference: Option<Arc<Vec<Record>>>,
    run: impl FnOnce() -> Result<T, VakintError>,
) -> Result<(T, Arc<Vec<Record>>, Counts, usize, usize, usize, f64, usize), VakintError> {
    struct Clear;
    impl Drop for Clear {
        fn drop(&mut self) {
            STATE.with(|state| {
                state.borrow_mut().take();
            });
        }
    }
    let (result, counts) = lowering_diagnostics::scoped(generic, || {
        STATE.with(|slot| {
            assert!(slot.borrow().is_none());
            *slot.borrow_mut() = Some(State {
                reference,
                records: Vec::new(),
                cursor: 0,
                bytes: 0,
                owner: std::thread::current().id(),
                in_search: false,
                best: None,
                searches: 0,
                candidates: 0,
                rechecks: 0,
                observer_seconds: 0.0,
            });
        });
        let _clear = Clear;
        let value = run()?;
        let state = STATE.with(|state| state.borrow_mut().take()).unwrap();
        state.complete()?;
        Ok::<_, VakintError>((
            value,
            Arc::new(state.records),
            state.searches,
            state.candidates,
            state.rechecks,
            state.observer_seconds,
            state.bytes,
        ))
    });
    let (value, records, searches, candidates, rechecks, seconds, bytes) = result?;
    if generic && counts.rational != 0 || counts.rational + counts.generic < candidates + rechecks {
        return Err(error(
            "weighted conversion counters do not cover the complete witness search",
        ));
    }
    Ok((
        value, records, counts, searches, candidates, rechecks, seconds, bytes,
    ))
}

fn check_lane(generic: bool, counts: Counts) {
    assert!(
        counts.rational + counts.generic > 0,
        "phase never exercised weighted conversion"
    );
    if generic {
        assert_eq!(counts.rational, 0, "forced generic phase used Q");
    } else {
        assert!(
            counts.rational > 0,
            "Q phase never exercised rational conversion"
        );
    }
}

#[test]
#[ignore = "actual external R2 package; full weighted witness parity under external30min/RSS guard"]
fn actual_mixed_witness_generic_and_rational_pair() -> Result<(), Box<dyn std::error::Error>> {
    let mut reference = None;
    let mut previous: Option<RustRedSymbolicEvaluation> = None;
    for generic in [true, false] {
        let initialized = Instant::now();
        let (engine, settings) = super::r2_acceptance::engine()?;
        let fixture = super::r2_acceptance::cancellation_fixtures::fixtures(&engine)?
            .into_iter()
            .find(|fixture| fixture.label == "mixed_dotted_lines")
            .ok_or("missing unchanged mixed fixture")?;
        assert_eq!(fixture.positive_cancellations, vec![(1, 2), (6, 3)]);
        assert_eq!(fixture.mass_squared, Atom::one());
        assert!(
            engine
                .saved_application(&settings)?
                .prepared
                .lock()
                .map_err(|_| "prepared lock poisoned")?
                .routed_entry_is_admitted(&fixture.expected_key)?
        );
        println!(
            "weighted_vk phase=initialized lane={} seconds={}",
            if generic { "generic" } else { "rational" },
            initialized.elapsed().as_secs_f64()
        );
        let started = Instant::now();
        let (first, records, counts, searches, candidates, rechecks, observation, bytes) =
            measured(generic, reference.clone(), || {
                engine.evaluate_rustred_symbolic(&settings, fixture.cancelled.as_view())
            })?;
        println!(
            "weighted_vk phase=first_complete lane={} seconds={} observer_seconds={observation} searches={searches} candidates={candidates} rechecks={rechecks} rational={} generic={} charged_payload={bytes}",
            if generic { "generic" } else { "rational" },
            started.elapsed().as_secs_f64(),
            counts.rational,
            counts.generic
        );
        assert!(!first.value.is_zero());
        check_lane(generic, counts);
        if generic {
            reference = Some(records);
        } else {
            assert!(counts.rational > 0);
        }
        if let Some(previous) = &previous {
            assert_eq!(previous.value, first.value);
            assert_eq!(previous.nonzero_conditions, first.nonzero_conditions);
            assert_eq!(previous.families, first.families);
        }
        // Separately timed controls: no witness record is allowed to masquerade
        // as a first-application speedup or independent master evaluation.
        let started = Instant::now();
        let (pinched, pinched_counts) = lowering_diagnostics::scoped(generic, || {
            engine.evaluate_rustred_symbolic(&settings, fixture.pinched.as_view())
        });
        let pinched = pinched?;
        check_lane(generic, pinched_counts);
        assert!(
            (&first.value - &pinched.value)
                .expand()
                .together()
                .is_zero()
        );
        assert_eq!(first.nonzero_conditions, pinched.nonzero_conditions);
        assert_eq!(first.families, pinched.families);
        println!(
            "weighted_vk phase=pinched_complete lane={} seconds={} counts={pinched_counts:?}",
            if generic { "generic" } else { "rational" },
            started.elapsed().as_secs_f64()
        );
        let before_warm = engine.rustred_application_witness_cache_statistics()?;
        let started = Instant::now();
        let (warm, warm_counts) = lowering_diagnostics::scoped(generic, || {
            engine.evaluate_rustred_symbolic(&settings, fixture.cancelled.as_view())
        });
        let warm = warm?;
        check_lane(generic, warm_counts);
        assert_eq!(warm.value, first.value);
        assert_eq!(warm.nonzero_conditions, first.nonzero_conditions);
        assert_eq!(warm.families, first.families);
        let after_warm = engine.rustred_application_witness_cache_statistics()?;
        assert!(
            after_warm.hits > before_warm.hits,
            "warm witness memo was not reused"
        );
        assert_eq!(
            after_warm.misses, before_warm.misses,
            "warm call started a new witness search"
        );
        assert_eq!(
            after_warm.invalidations, before_warm.invalidations,
            "warm admission invalidated the witness"
        );
        println!(
            "weighted_vk phase=warm_complete lane={} seconds={} counts={warm_counts:?} witness_cache={:?}",
            if generic { "generic" } else { "rational" },
            started.elapsed().as_secs_f64(),
            after_warm
        );
        previous = Some(first);
    }
    println!("weighted_vk phase=complete exact_full_candidate_and_minimum_witness_parity=true");
    Ok(())
}

#[cfg(test)]
mod protocol_tests {
    use super::*;
    fn state(reference: Option<Arc<Vec<Record>>>) -> State {
        State {
            reference,
            records: Vec::new(),
            cursor: 0,
            bytes: 0,
            owner: std::thread::current().id(),
            in_search: false,
            best: None,
            searches: 1,
            candidates: 0,
            rechecks: 0,
            observer_seconds: 0.0,
        }
    }
    #[test]
    fn diagnostic_scope_restores_after_panic_and_refuses_thread_handoff() {
        let panic = std::panic::catch_unwind(|| {
            let _: Result<((), _, _, _, _, _, _, _), _> =
                measured(true, None, || panic!("intentional scoped probe panic"));
        });
        assert!(panic.is_err());
        assert!(lowering_diagnostics::current().is_none());
        assert!(STATE.with(|state| state.borrow().is_none()));
        let (_, counts) = lowering_diagnostics::scoped(true, || {
            assert!(lowering_diagnostics::current().is_some());
            let nested = std::panic::catch_unwind(|| lowering_diagnostics::scoped(false, || ()));
            assert!(nested.is_err());
            assert!(lowering_diagnostics::current().is_some());
            let owner = std::thread::current().id();
            assert!(
                std::thread::spawn(move || assert_thread(Some(owner)).is_err())
                    .join()
                    .unwrap()
            );
        });
        assert_eq!(
            counts,
            Counts {
                rational: 0,
                generic: 0
            }
        );
        assert!(lowering_diagnostics::current().is_none());
    }
    #[test]
    fn diagnostic_records_refuse_limits_missing_extra_and_changed_values() {
        let row = Record::Begin {
            integral: Atom::one(),
            numerator: Atom::num(2),
        };
        let mut measured = state(Some(Arc::new(vec![row.clone()])));
        assert!(measured.complete().is_err());
        measured.record(row.clone()).unwrap();
        assert!(measured.complete().is_ok());
        assert!(measured.record(row.clone()).is_err());
        let mut changed = state(Some(Arc::new(vec![row.clone()])));
        assert!(
            changed
                .record(Record::Begin {
                    integral: Atom::one(),
                    numerator: Atom::num(3)
                })
                .is_err()
        );
        let mut limited = state(None);
        limited.cursor = MAX_RECORDS;
        assert!(limited.charge(0).is_err());
        limited.cursor = 0;
        limited.bytes = MAX_PAYLOAD_BYTES;
        assert!(limited.charge(1).is_err());
        let mut unfinished = state(None);
        unfinished.in_search = true;
        assert!(unfinished.complete().is_err());
    }
    #[test]
    fn diagnostic_minimum_keeps_first_equal_key_witness() {
        let key = vec![(IntegralKey::try_new([1]).unwrap(), Atom::one())];
        let first = ReplacementRules::default();
        let mut second = first.clone();
        second.edge_ids_canonical_to_input_map.insert(1, 2);
        let mut observed = state(None);
        observed.consider(&key, &first);
        observed.consider(&key, &second);
        assert!(same_witness(&observed.best.as_ref().unwrap().1, &first));
        let lower = vec![(IntegralKey::try_new([0]).unwrap(), Atom::one())];
        observed.consider(&lower, &second);
        assert!(same_witness(&observed.best.as_ref().unwrap().1, &second));
    }
}
