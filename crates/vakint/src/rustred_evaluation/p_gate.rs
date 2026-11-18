//! Gate hooks of the reference path R2 (helper VH, numerator assignment, lane
//! `refpath`): what RustRed receives from vakint's ordinary p-integral
//! evaluation, as a comparison view that another engine can be gated against
//! key by key (`N/notes/DESIGN-0-ANATOMY.md` section 2.4).
//!
//! * [`observe`] runs a closure with an observer on this thread: every batch
//!   of folded channels that [`super::p_integral::reduce_requests_folded`]
//!   hands to RustRed is cloned into it ([`observe_folded`]), with the
//!   statistics of the reduction ([`observe_stats`]) and the finalizer's
//!   seconds ([`observe_finalize`]). Without an observer every hook is one
//!   thread-local check: values and outputs never depend on it.
//! * [`gate_view`] turns folded channels into the view RustRed receives: per
//!   key `(loops, scale_power, presentation, powers)` the coefficient
//!   `Σ_{channel, basis} channel.coefficient · basis(ε → (4 - D)/2) · value`,
//!   normalised by `together()`. It does not depend on how the coefficient is
//!   split between the channel's spectator monomial and the ε basis, nor on
//!   how channels are folded.
//! * [`drop_zero_sectors`] removes the keys whose sector the library's census
//!   marks `Zero` (RustRed returns 0 for them without work, so the gate drops
//!   them on both sides); a support missing from the census is kept.
//! * [`write_gate_view`] / [`read_gate_view`] write and read a canonical text
//!   form (`vakint.p-gate-view.v1`), [`compare_gate_views`] compares two
//!   views exactly (string equality, else an exact rational check).

use std::cell::RefCell;
use std::collections::BTreeMap;
use std::fmt::Write as _;
use std::time::Duration;

use rustred::propagator::PresentationId;
use rustred::propagator::census::PSectorFate;
use symbolica::atom::{Atom, AtomCore, Symbol};

use super::p_data::PLibraryHandle;
use super::p_integral::{Channel, FoldedInputs, PReduceStats};

/// One key of the gate view: what RustRed reduces, without the external head
/// (a two-point function has one).
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct GateKey {
    pub loops: usize,
    pub scale_power: i64,
    pub presentation: u32,
    pub powers: Vec<i64>,
}

/// The gate view: key → exact coefficient (a rational function of `D` times
/// spectator symbols), normalised by `together()`.
pub type GateView = BTreeMap<GateKey, Atom>;

/// The symbol `D` of the gate view (the space-time dimension, `ε = (4 - D)/2`).
pub fn gate_dimension() -> Symbol {
    symbolica::symbol!("vakint_gate::D")
}

/// What one observed evaluation handed to RustRed and how long the stages
/// took (sums over the reduction batches of the evaluation).
#[derive(Clone, Debug, Default)]
pub struct GateObservation {
    /// Every observed batch of folded channels, with the census verdicts of
    /// its library.
    pub(crate) batches: Vec<ObservedBatch>,
    /// Seconds spent cloning the folded inputs for the observer (to be
    /// subtracted from a wall time that contains them).
    pub observer_seconds: f64,
    /// [`PReduceStats`] summed over the observed reductions.
    pub split_seconds: f64,
    pub lower_seconds: f64,
    pub prepare_seconds: f64,
    pub library_seconds: f64,
    pub render_seconds: f64,
    pub finalize_seconds: f64,
    pub reductions: usize,
    pub pairs: usize,
    pub channels: usize,
    pub folded: usize,
    pub memo_hits: usize,
    /// Plain requests the direct lowering (`p_lowering::lower_request`) took,
    /// and those it declined by reason (the ordinary split and lowering ran);
    /// requests served by the request memo are in neither.
    pub direct_lowered: usize,
    pub direct_declined: BTreeMap<&'static str, usize>,
}

/// One observed batch: its folded channels and the zero sectors of its library.
#[derive(Clone, Debug)]
pub(crate) struct ObservedBatch {
    pub(crate) epsilon: Symbol,
    pub(crate) folded: BTreeMap<Channel, FoldedInputs>,
    /// Per presentation: (forced axes, support → zero?) of every support of
    /// the census. Only presentations the batch uses.
    pub(crate) zero: BTreeMap<PresentationId, (u64, BTreeMap<u64, bool>)>,
}

impl GateObservation {
    /// The gate view of everything observed, every key kept.
    pub fn view(&self) -> GateView {
        let mut view = GateView::new();
        for batch in &self.batches {
            add_view(&mut view, &batch.folded, batch.epsilon);
        }
        normalise(view)
    }

    /// [`Self::view`] with the zero-sector keys dropped, and the number of
    /// keys dropped and of keys whose support the census does not list.
    pub fn view_without_zero_sectors(&self) -> (GateView, ZeroSectorCount) {
        let view = self.view();
        let mut zero = BTreeMap::<PresentationId, (u64, BTreeMap<u64, bool>)>::new();
        for batch in &self.batches {
            for (presentation, entry) in &batch.zero {
                zero.entry(*presentation).or_insert_with(|| entry.clone());
            }
        }
        drop_zero_sectors_with(view, |key| {
            let (forced, supports) = zero.get(&PresentationId(key.presentation))?;
            supports.get(&support(*forced, &key.powers)).copied()
        })
    }

    /// Distinct keys of the folded inputs (before the zero filter): the keys
    /// RustRed receives, over all channels.
    pub fn distinct_keys(&self) -> usize {
        self.view().len()
    }
}

/// Keys dropped by [`drop_zero_sectors`] and keys whose support has no census
/// entry (kept).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct ZeroSectorCount {
    pub total: usize,
    pub zero: usize,
    pub unknown: usize,
}

thread_local! {
    static OBSERVER: RefCell<Option<GateObservation>> = const { RefCell::new(None) };
}

fn observing() -> bool {
    OBSERVER.with(|observer| observer.borrow().is_some())
}

/// Run `work` with a gate observer on this thread and return its result and
/// the observation. Nested calls are not supported (the inner one replaces
/// the outer observer for its duration and restores it afterwards).
pub fn observe<T>(work: impl FnOnce() -> T) -> (T, GateObservation) {
    let outer = OBSERVER.with(|observer| observer.replace(Some(GateObservation::default())));
    let value = work();
    let observation = OBSERVER
        .with(|observer| observer.replace(outer))
        .unwrap_or_default();
    (value, observation)
}

/// The census verdicts of the presentations of `folded`: for each, its forced
/// axes and every listed support's zero verdict.
fn zero_sectors(
    handle: &PLibraryHandle,
    folded: &BTreeMap<Channel, FoldedInputs>,
) -> BTreeMap<PresentationId, (u64, BTreeMap<u64, bool>)> {
    let library = &handle.library;
    let census = library.census();
    let mut presentations = std::collections::BTreeSet::new();
    for inputs in folded.values() {
        for keys in inputs.by_basis.values() {
            for (presentation, _) in keys.keys() {
                presentations.insert(*presentation);
            }
        }
    }
    presentations
        .into_iter()
        .map(|presentation| {
            let forced = library.presentation(presentation).forced_axes();
            let supports = census
                .sectors
                .range((presentation, 0)..=(presentation, u64::MAX))
                .map(|((_, support), fate)| (*support, matches!(fate, PSectorFate::Zero)))
                .collect();
            (presentation, (forced, supports))
        })
        .collect()
}

/// Hook of [`super::p_integral::reduce_requests_folded`]: clone the folded
/// channels into this thread's observer, if any.
pub(crate) fn observe_folded(
    handle: &PLibraryHandle,
    epsilon: Symbol,
    folded: &BTreeMap<Channel, FoldedInputs>,
) {
    if !observing() {
        return;
    }
    let started = std::time::Instant::now();
    let batch = ObservedBatch {
        epsilon,
        folded: folded.clone(),
        zero: zero_sectors(handle, folded),
    };
    let seconds = started.elapsed().as_secs_f64();
    OBSERVER.with(|observer| {
        if let Some(observation) = observer.borrow_mut().as_mut() {
            observation.batches.push(batch);
            observation.observer_seconds += seconds;
        }
    });
}

/// Hook of [`super::p_integral::reduce_requests`]: add the statistics of the
/// reduction that just ended to this thread's observer, if any.
pub(crate) fn observe_stats(stats: &PReduceStats) {
    OBSERVER.with(|observer| {
        if let Some(observation) = observer.borrow_mut().as_mut() {
            observation.split_seconds += stats.split_seconds;
            observation.lower_seconds += stats.lower_seconds;
            observation.prepare_seconds += stats.prepare_seconds;
            observation.library_seconds += stats.library_seconds;
            observation.render_seconds += stats.render_seconds;
            observation.reductions += stats.reductions;
            observation.pairs += stats.pairs;
            observation.channels += stats.channels;
            observation.folded += stats.folded;
            observation.memo_hits += stats.memo_hits;
        }
    });
}

/// Hook of `fold_request`: the direct lowering took a plain request (`None`)
/// or declined it (`Some(reason)`).
pub(crate) fn observe_direct(declined: Option<&'static str>) {
    OBSERVER.with(|observer| {
        if let Some(observation) = observer.borrow_mut().as_mut() {
            match declined {
                None => observation.direct_lowered += 1,
                Some(reason) => *observation.direct_declined.entry(reason).or_default() += 1,
            }
        }
    });
}

/// Hook of [`super::p_integral::evaluate_ordered_batch`]: the finalizer's time.
pub(crate) fn observe_finalize(elapsed: Duration) {
    OBSERVER.with(|observer| {
        if let Some(observation) = observer.borrow_mut().as_mut() {
            observation.finalize_seconds += elapsed.as_secs_f64();
        }
    });
}

/// Whether a gate observer is active on this thread (lets callers skip
/// timing work that only the observer reads).
pub(crate) fn active() -> bool {
    observing()
}

/// `raw | forced` of `powers` (RustRed's sector support: slots with a positive
/// power, and the presentation's forced ε axes).
fn support(forced: u64, powers: &[i64]) -> u64 {
    powers
        .iter()
        .enumerate()
        .filter(|&(_, &power)| power > 0)
        .fold(forced, |mask, (slot, _)| mask | 1 << slot)
}

fn add_view(view: &mut GateView, folded: &BTreeMap<Channel, FoldedInputs>, epsilon: Symbol) {
    let dimension = Atom::var(gate_dimension());
    let epsilon_at_d = (Atom::num(4) - &dimension) / Atom::num(2);
    for (channel, inputs) in folded {
        for (basis, keys) in &inputs.by_basis {
            let weight = &channel.coefficient
                * basis
                    .replace(Atom::var(epsilon).to_pattern())
                    .with(epsilon_at_d.to_pattern());
            for ((presentation, powers), value) in keys {
                if value.is_zero() {
                    continue;
                }
                let key = GateKey {
                    loops: channel.loops,
                    scale_power: channel.scale_power,
                    presentation: presentation.0,
                    powers: powers.clone(),
                };
                *view.entry(key).or_default() += &weight * Atom::num(value.clone());
            }
        }
    }
}

/// The canonical value of a gate coefficient: `together()` of the expanded
/// sum (a rational function in `D` and the spectators, numerator and
/// denominator in Symbolica's canonical form).
pub fn normalise_value(value: &Atom) -> Atom {
    value.expand().together()
}

fn normalise(view: GateView) -> GateView {
    view.into_iter()
        .map(|(key, value)| (key, normalise_value(&value)))
        .filter(|(_, value)| !value.is_zero())
        .collect()
}

/// The gate view of folded channels (every key kept; `epsilon` is the
/// settings' ε symbol of the basis functions).
#[cfg_attr(not(test), allow(dead_code))]
pub(crate) fn gate_view(folded: &BTreeMap<Channel, FoldedInputs>, epsilon: Symbol) -> GateView {
    let mut view = GateView::new();
    add_view(&mut view, folded, epsilon);
    normalise(view)
}

/// Drop the keys whose sector `zero(key)` says is zero (`Some(true)`); keep
/// `Some(false)` and count `None` (no census entry) as unknown.
pub fn drop_zero_sectors_with(
    view: GateView,
    mut zero: impl FnMut(&GateKey) -> Option<bool>,
) -> (GateView, ZeroSectorCount) {
    let mut count = ZeroSectorCount {
        total: view.len(),
        ..Default::default()
    };
    let kept = view
        .into_iter()
        .filter(|(key, _)| match zero(key) {
            Some(true) => {
                count.zero += 1;
                false
            }
            Some(false) => true,
            None => {
                count.unknown += 1;
                true
            }
        })
        .collect();
    (kept, count)
}

/// [`drop_zero_sectors_with`] by the census of `handle`'s library.
#[cfg_attr(not(test), allow(dead_code))]
pub(crate) fn drop_zero_sectors(
    handle: &PLibraryHandle,
    view: GateView,
) -> (GateView, ZeroSectorCount) {
    let library = &handle.library;
    let census = library.census();
    drop_zero_sectors_with(view, |key| {
        let presentation = PresentationId(key.presentation);
        let forced = library.presentation(presentation).forced_axes();
        census
            .sectors
            .get(&(presentation, support(forced, &key.powers)))
            .map(|fate| matches!(fate, PSectorFate::Zero))
    })
}

/// The schema line of the text form.
pub const GATE_VIEW_SCHEMA: &str = "vakint.p-gate-view.v1";

/// The canonical text of a gate view: the schema line, `# key value` header
/// lines (sorted by key), then one line per key, sorted:
/// `loops<TAB>scale_power<TAB>presentation<TAB>p1,p2,..<TAB>value` with the
/// value's `to_canonical_string()`.
pub fn write_gate_view(view: &GateView, header: &BTreeMap<String, String>) -> String {
    let mut text = String::new();
    let _ = writeln!(text, "{GATE_VIEW_SCHEMA}");
    for (key, value) in header {
        let _ = writeln!(text, "# {key} {}", value.replace('\n', " "));
    }
    let _ = writeln!(text, "# keys {}", view.len());
    let _ = writeln!(text, "# dimension {}", Atom::var(gate_dimension()).to_canonical_string());
    for (key, value) in view {
        let powers = key
            .powers
            .iter()
            .map(i64::to_string)
            .collect::<Vec<_>>()
            .join(",");
        let _ = writeln!(
            text,
            "{}\t{}\t{}\t{}\t{}",
            key.loops,
            key.scale_power,
            key.presentation,
            powers,
            value.to_canonical_string()
        );
    }
    text
}

/// Read the text of [`write_gate_view`]: the view and the header lines.
pub fn read_gate_view(text: &str) -> Result<(GateView, BTreeMap<String, String>), String> {
    let mut lines = text.lines();
    match lines.next() {
        Some(GATE_VIEW_SCHEMA) => {}
        other => return Err(format!("not a {GATE_VIEW_SCHEMA} file (first line {other:?})")),
    }
    let mut view = GateView::new();
    let mut header = BTreeMap::new();
    for (number, line) in lines.enumerate() {
        if let Some(rest) = line.strip_prefix("# ") {
            let (key, value) = rest.split_once(' ').unwrap_or((rest, ""));
            header.insert(key.to_string(), value.to_string());
            continue;
        }
        if line.trim().is_empty() {
            continue;
        }
        let fields = line.split('\t').collect::<Vec<_>>();
        let [loops, scale_power, presentation, powers, value] = fields[..] else {
            return Err(format!("line {}: {} fields", number + 2, fields.len()));
        };
        let parse = |field: &str| -> Result<i64, String> {
            field
                .parse::<i64>()
                .map_err(|error| format!("line {}: {field}: {error}", number + 2))
        };
        let key = GateKey {
            loops: parse(loops)? as usize,
            scale_power: parse(scale_power)?,
            presentation: parse(presentation)? as u32,
            powers: if powers.is_empty() {
                Vec::new()
            } else {
                powers.split(',').map(parse).collect::<Result<_, _>>()?
            },
        };
        let value = symbolica::try_parse!(value)
            .map_err(|error| format!("line {}: {error}", number + 2))?;
        if view.insert(key, value).is_some() {
            return Err(format!("line {}: repeated key", number + 2));
        }
    }
    if let Some(keys) = header.get("keys")
        && keys.parse::<usize>().ok() != Some(view.len())
    {
        return Err(format!("header says {keys} keys, read {}", view.len()));
    }
    Ok((view, header))
}

/// The differences between two gate views.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct GateDifference {
    /// Keys only in the left (reference) view.
    pub only_left: Vec<GateKey>,
    /// Keys only in the right view.
    pub only_right: Vec<GateKey>,
    /// Keys in both whose coefficients differ exactly: (key, left − right).
    pub differ: Vec<(GateKey, Atom)>,
    /// Keys equal as text.
    pub identical: usize,
    /// Keys equal only after the exact rational check.
    pub equal: usize,
}

impl GateDifference {
    pub fn is_empty(&self) -> bool {
        self.only_left.is_empty() && self.only_right.is_empty() && self.differ.is_empty()
    }
}

/// Compare two gate views key by key, exactly: a coefficient pair is equal
/// when its canonical strings agree or `(left - right).together()` is zero.
pub fn compare_gate_views(left: &GateView, right: &GateView) -> GateDifference {
    let mut difference = GateDifference::default();
    for (key, value) in left {
        match right.get(key) {
            None => difference.only_left.push(key.clone()),
            Some(other) => {
                if value.to_canonical_string() == other.to_canonical_string() {
                    difference.identical += 1;
                } else {
                    let delta = (value - other).expand().together();
                    if delta.is_zero() {
                        difference.equal += 1;
                    } else {
                        difference.differ.push((key.clone(), delta));
                    }
                }
            }
        }
    }
    difference.only_right = right
        .keys()
        .filter(|key| !left.contains_key(*key))
        .cloned()
        .collect();
    difference
}

#[cfg(test)]
mod tests {
    use super::*;

    fn key(loops: usize, powers: &[i64]) -> GateKey {
        GateKey {
            loops,
            scale_power: 1,
            presentation: 3,
            powers: powers.to_vec(),
        }
    }

    /// The text form round-trips, and the comparison tells equal values in
    /// another form from different values and missing keys.
    #[test]
    fn gate_views_round_trip_and_compare_exactly() {
        crate::Vakint::initialize_vakint_symbols();
        let d = Atom::var(gate_dimension());
        let ca = Atom::var(symbolica::symbol!("vakint_gate_test::CA"));
        let mut view = GateView::new();
        view.insert(key(2, &[1, 1, 0, -1]), normalise_value(&(&ca * (&d - Atom::num(4)) / (&d - Atom::num(1)))));
        view.insert(key(2, &[1, 1, 1, 0]), normalise_value(&(Atom::num(7) / Atom::num(3))));
        let header = BTreeMap::from([("graph".to_string(), "test".to_string())]);
        let text = write_gate_view(&view, &header);
        let (read, read_header) = read_gate_view(&text).unwrap();
        assert_eq!(read, view);
        assert_eq!(read_header["graph"], "test");
        assert_eq!(write_gate_view(&read, &header), text);
        assert!(compare_gate_views(&view, &read).is_empty());

        // The same coefficient written otherwise is equal, not identical.
        let mut other = view.clone();
        other.insert(
            key(2, &[1, 1, 0, -1]),
            &ca * (&d - Atom::num(4)) * (&d - Atom::num(1)).pow(Atom::num(-1)) + Atom::Zero,
        );
        let difference = compare_gate_views(&view, &other);
        assert!(difference.is_empty(), "{difference:?}");
        // A changed value and a missing key are reported.
        other.insert(key(2, &[1, 1, 1, 0]), Atom::num(8) / Atom::num(3));
        other.insert(key(2, &[0, 1, 1, 0]), Atom::num(1));
        let difference = compare_gate_views(&view, &other);
        assert_eq!(difference.differ.len(), 1);
        assert_eq!(difference.only_right, vec![key(2, &[0, 1, 1, 0])]);
        assert!(difference.only_left.is_empty());
    }

    /// Zero sectors are dropped by the verdict, unknown supports kept and
    /// counted; the support is `raw | forced`.
    #[test]
    fn zero_sectors_are_dropped_by_the_support_verdict() {
        assert_eq!(support(0, &[1, 0, -2, 3]), 0b1001);
        assert_eq!(support(0b10, &[1, 0, -2, 3]), 0b1011);
        let mut view = GateView::new();
        view.insert(key(1, &[1, 1]), Atom::num(1));
        view.insert(key(1, &[1, 0]), Atom::num(2));
        view.insert(key(1, &[0, 1]), Atom::num(3));
        let (kept, count) = drop_zero_sectors_with(view, |key| match support(0, &key.powers) {
            0b11 => Some(false),
            0b01 => Some(true),
            _ => None,
        });
        assert_eq!(kept.len(), 2);
        assert_eq!(
            count,
            ZeroSectorCount {
                total: 3,
                zero: 1,
                unknown: 1
            }
        );
    }
}
