//! IRH-DEPTH: the Laurent depth of an IR hard part on the child route
//! (`GAMMALOOP_IR_HARD=child`, `Evaluation::ir_hard_child`).
//!
//! The child route integrates the hard part `H(γ) = Σ_j h_j E_j` of one IR
//! correction with `L + extra` Laurent terms (`L` the remainder's loops):
//! every `h_j` is exact through `ε^(extra-1)`, whatever its own poles, since
//! the single-mass integration truncates the product with the insertion
//! coefficients there (`integrate_single_mass`). It then multiplies each row
//! by its quotient factors `z_j = Σ_F Z̃_F(G/γ; E_j) · closure`. When `z_j`
//! starts at `ε^-p`, the poles of `h_j z_j` need `h_j` through `ε^(p-1)`, so
//! the depth is right once `extra >= p` for every row. The planner starts at
//! `extra` = the quotient's loops and measures `p` on the rows the truncated
//! `H` shows, retrying with a larger `extra` when a row needs it.
//!
//! A hard group (one remainder loop count) whose truncated value vanishes
//! shows no row at all, and a row whose coefficient vanishes is skipped
//! before its quotient factors are built (or is dropped by the moment
//! extraction): their `p` is never measured. With a true valuation `v` of
//! the hard coefficient and `extra <= v < p` they still contribute a pole.
//!
//! [`quotient_pole_bound`] bounds every quotient factor's pole order before
//! the hard part is truncated, independently of which rows survive:
//!
//! * `Z̃_F` is the c-operation over the forest's steps; each step is the UV
//!   pole part (a `K`, through R* on the decorated quotient or the massive R
//!   operation) of a region of the quotient, its coefficients kept outside
//!   the next step. The steps' regions partition the quotient `G/γ`, so their
//!   loops add up to the quotient's loops. A `K` of an `l`-loop region whose
//!   factors are regular at `ε = 0` has poles of order at most `l` (every
//!   term of its R-bar is a product of subgraph `K`s and a remaining
//!   integral, each bounded by its loops); a factor `f` of order `ord_ε f`
//!   raises that bound by at most `max(0, -ord_ε f)`, in whichever term the
//!   factor lands. Every quotient vertex and propagator numerator is such a
//!   factor (the decorated quotient keeps them; the contracted `γ` becomes
//!   one vertex carrying `E_j`). A propagator-insertion surrogate's
//!   coefficient vertex is one of them, regular by construction (an N-loop
//!   insertion over an N-loop scalar chain).
//! * `E_j` holds only the covariant tensor structure (scalars, `ε` among
//!   them, go to `h_j`, [`super::moments::hard_moments`]); a Dirac closure
//!   adds its gamma structure (`ir_dirac::reassemble` puts its scalar
//!   normalisation into the coefficient). Its own order is counted per row.
//!   A row the moment extraction DROPPED has no insertion to read: it is
//!   bounded by `B` only through this invariant (`E_j` free of `ε` outside
//!   representation dimensions). The planner checks it on every row the
//!   extraction returns ([`dropped_row_bound`]); one row with an insertion
//!   excess means a dropped row's could be larger, so the dropped rows are
//!   unbounded and the correction is refused (`unbounded`).
//! * The root closure multiplies `z_j` (a metric or Dirac closure, regular
//!   at `d = 4`; counted anyway), the outer Taylor step takes polynomial
//!   orders in the root momenta, and the fixed-mass compensation
//!   `Σ_k (-1)^k binom(a-εL, k)` is a polynomial in `ε`.
//! * Dimension factors are read at `d = 4-2ε` (as the valuation of `z_j`
//!   is), so a `1/(d-4)` would count as a pole. Vacuum tensor reduction and
//!   the UV-child closures only divide by `d + 2k` or `d - 1`.
//!
//! Hence `p_j <= B = quotient loops + excess(quotient numerator) +
//! excess(closure) + excess(E_j)`, with `excess` of [`excess_poles`]. When
//! `B <= extra`, no row can hide a pole.
//!
//! The planner measures the quotient factors of every row it is given,
//! also of a row whose coefficient vanishes at the current depth, and takes
//! `B` for the valuation of what it cannot measure (a vanished group, a row
//! the moment extraction dropped): below `B` it raises `extra` to `B` through
//! the existing retry. At depth `B` every group and row is either shown, and
//! measured, or of valuation at least `B`, out of reach of every pole.
//!
//! The planner starts at `max(quotient loops, B)` (no attempt below `B` is
//! ever useful: it would be raised to `B` after a whole integration of `H`);
//! on production `B` is the quotient's loops, the old first depth.
//!
//! The guard stays as a backstop: it refuses the child route for a
//! correction whose accepted depth stays below `B`, for every correction
//! whose `B` is unknown (a row of a nonzero group that the moment extraction
//! dropped is counted nowhere, so without `B` nothing certifies that no pole
//! hides; refused before `H` is integrated when the quotient's own bound is
//! unknown, and before the quotient factors are built when a row's or the
//! dropped rows' bound is unknown), and for a measured row above its own
//! bound (`bound_exceeded`, compared row by row: the proof of `B` failed;
//! refused at once, an alarm of its own, fatal under
//! `GAMMALOOP_IR_HARD_DEPTH_STRICT=1`, also under the opt-in closed-path
//! guard `GAMMALOOP_IR_HARD_CHILD_GUARD=1`, and at any nesting depth: a
//! nested evaluation's alarm is a typed error, [`DepthAlarm`], that no guard,
//! vanished-row valuation or counterterm cache turns into a decline). The
//! correction is then
//! evaluated on today's path (the closed or projected moments with the
//! generic forest step, which measures every quotient factor before its
//! hard coefficient is integrated), so a refusal costs time, never a value.
//! Every refusal is visible three times: its `IR_HARD_DEPTH_REFUSED` line,
//! the correction's `IR_HARD_CHILD … outcome=declined:depth_refused:<reason>`
//! line (the line every log filter keeps), and the per-reason counts of
//! `IR_HARD_DEPTH_SUMMARY`.
//!
//! # Domain of `B`
//!
//! `B` reads the `ε` dependence of the quotient where it is explicit: the
//! vertex and propagator NUMERATORS outside the hard part
//! (`Graph::numerator`, the same source the R* integrand uses), the root
//! closure and each row's insertion. It assumes that nothing else of the
//! quotient depends on `ε`:
//! * propagator powers and denominators: on this route every propagator is
//!   `1/p²` with power one (feynkit's `propagator_power` is an `i32` and the
//!   runtime graph rejects any power but 1, `graph/parse`; R* requires a
//!   massless graph, `Evaluation::new`), so an `ε`-dependent power or mass
//!   cannot reach the quotient;
//! * opaque coefficient symbols: a symbol whose value, substituted later,
//!   has poles in `ε` is invisible to `B` (it is read as regular). No
//!   production input has one (no `ε` and no such symbol in any 5L/6L
//!   diagram numerator, census-6l.md). Production numerators DO carry `ε`,
//!   though: with `--insertions extract` the propagator-insertion
//!   surrogates are built at run time, and their coefficient vertices carry
//!   the normalised insertion explicitly, a polynomial in `ε` (with rational
//!   and `ζ` coefficients and rational functions of `d` = `4-2ε` that are
//!   regular at `d = 4`), regular by construction. So `B` reads it in full
//!   (the series per factor) and finds no excess
//!   (`tests::a_surrogate_coefficient_has_no_excess`); B = quotient loops
//!   holds on production because those series are regular, not because `ε`
//!   is absent. A future input with such a symbol is outside the domain;
//! * a function of `ε` (`f(ε)`, any function whose arguments carry `ε` after
//!   `d = 4-2ε`): its order is not read from a series; detectable, so `B` is
//!   unknown and the correction is refused ([`excess_poles`]). A surrogate
//!   coefficient that kept an unexpanded function of `ε` would therefore
//!   refuse every correction of its diagram into today's path: right value,
//!   today's cost, counted as `unbounded=` in `IR_HARD_DEPTH_SUMMARY`.
//!
//! Round 5 ([`super::regularity`]): `B` is read on the NORMALISED factors
//! of a closed list (closed traces contracted, `d = 4-2ε`), the graph
//! prefactor and projector included; a factor that is not provably regular
//! refuses the whole correction before the hard part is integrated.

use symbolica::{
    atom::{Atom, AtomCore, AtomView},
    poly::series::SeriesDepth,
};

use crate::{utils::GS, uv::approx::integrated::Integrated};

/// `Σ max(0, -ord_ε f)` over the product factors `f` of `expression`, read
/// at `d = 4-2ε`: how many poles beyond its loops a `K` of a region carrying
/// these factors can gain. `Some(0)` without `ε` (or `d`); `None` when a
/// factor's order is not an integer or cannot be determined, and when a
/// function's argument carries `ε` (its order would be guessed from a
/// series of an opaque function: the module doc's domain).
pub(super) fn excess_poles(expression: &Atom) -> Option<usize> {
    let epsilon = GS.dim_epsilon;
    if !expression.contains_symbol(epsilon) && !expression.contains_symbol(GS.dim) {
        return Some(0);
    }
    let regulated = Integrated::dimensionally_regularized(expression);
    if !regulated.contains_symbol(epsilon) {
        return Some(0);
    }
    // A representation's dimension (`mink(4-2ε, i)`) is not a function of
    // `ε` in this sense: drop the representations before looking.
    let representation = &spenso::network::tags::SPENSO_TAG.representation;
    let scalars = regulated.replace_map(|view, _, out| {
        if let AtomView::Fun(function) = view
            && function.get_symbol().has_tag(representation)
        {
            **out = Atom::num(1);
        }
    });
    let mut opaque = false;
    scalars.visitor(&mut |view| {
        if let AtomView::Fun(function) = view
            && function
                .iter()
                .any(|argument| argument.contains_symbol(epsilon))
        {
            opaque = true;
        }
        !opaque
    });
    if opaque {
        return None;
    }
    let factors = match regulated.as_view() {
        AtomView::Mul(product) => product.iter().collect::<Vec<_>>(),
        view => vec![view],
    };
    let mut excess = 0_usize;
    for factor in factors {
        if !factor.contains_symbol(epsilon) {
            continue;
        }
        let series = factor.series(epsilon, 0, SeriesDepth::relative(1)).ok()?;
        if series.is_zero() {
            continue;
        }
        let exponent = series.get_trailing_exponent();
        if !exponent.is_integer() {
            return None;
        }
        let order = exponent.numerator().to_i64()?;
        excess += usize::try_from(-order.min(0)).ok()?;
    }
    Some(excess)
}

/// The pole bound `B` of every quotient factor of one IR correction
/// (module doc): the quotient's loops plus the excess of its numerator (the
/// vertex and propagator factors outside the hard part) and of the root
/// closure. `None` when an excess is unknown.
pub(super) fn quotient_pole_bound(
    quotient_loops: usize,
    quotient_numerator: &Atom,
    closure: Option<&Atom>,
) -> Option<usize> {
    let numerator = excess_poles(quotient_numerator)?;
    let closure = match closure {
        Some(closure) => excess_poles(closure)?,
        None => 0,
    };
    Some(quotient_loops + numerator + closure)
}

/// The bound of one row: `bound` raised by the excess of its insertion
/// (tests; the planner reads the excesses once, [`dropped_row_bound`]).
#[cfg(test)]
pub(super) fn row_bound(bound: Option<usize>, insertion: &Atom) -> Option<usize> {
    Some(bound? + excess_poles(insertion)?)
}

/// The bound of the rows the moment extraction dropped, from the insertion
/// excesses of the rows it returned (module doc): `bound` when every one is
/// `Some(0)` (the invariant: `E_j` is free of `ε` outside representation
/// dimensions), else `None` (a dropped row's excess could be larger).
pub(super) fn dropped_row_bound(
    bound: Option<usize>,
    excesses: impl IntoIterator<Item = Option<usize>>,
) -> Option<usize> {
    let bound = bound?;
    excesses
        .into_iter()
        .all(|excess| excess == Some(0))
        .then_some(bound)
}

/// The bound of the whole correction: the dropped rows' and every shown
/// row's (`None` when one is unknown).
pub(super) fn correction_bound(dropped: Option<usize>, rows: &[Option<usize>]) -> Option<usize> {
    rows.iter()
        .try_fold(dropped?, |largest, row| Some(largest.max((*row)?)))
}

/// The first row whose measured need exceeds its own bound (`needs[i]` and
/// `bounds[i]` per row), as `(need, bound)`.
pub(super) fn exceeded_row(needs: &[usize], bounds: &[Option<usize>]) -> Option<(usize, usize)> {
    needs.iter().zip(bounds).find_map(|(need, bound)| {
        bound
            .filter(|bound| need > bound)
            .map(|bound| (*need, bound))
    })
}

/// What the truncated hard part of one correction shows at depth `extra`.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub(super) struct Truncation {
    /// Hard groups (remainder loop counts) in the hard part.
    pub groups: usize,
    /// Groups whose truncated value vanishes (no row to measure).
    pub vanished_groups: usize,
    /// Rows of the surviving groups.
    pub rows: usize,
    /// Rows whose truncated coefficient vanishes.
    pub zero_rows: usize,
}

impl Truncation {
    pub fn vanished(&self) -> usize {
        self.vanished_groups + self.zero_rows
    }
}

/// The guard (module doc), at the depth `extra` the planner would accept
/// with `needed` measured on the shown rows and `bound` = `B` (the largest
/// row bound included): `Some(reason)` refuses the child route.
/// * `unbounded`: `B` is unknown. Refused whatever the truncation shows: a
///   row of a nonzero group that the moment extraction dropped is counted
///   nowhere (not in [`Truncation::vanished`]), so nothing certifies that no
///   pole hides.
/// * `bound_exceeded`: a measured row needs more than `B`; the bound does not
///   hold, so nothing it certifies can be trusted.
/// * `vanished`: a vanished group or row may carry up to `B > extra` poles.
/// * `unverified`: `B > extra` with every group and row shown; a row of a
///   nonzero group whose coefficient vanished at the truncation is dropped
///   by the moment extraction and cannot be told apart from an absent one.
pub(super) fn refusal(
    bound: Option<usize>,
    extra: usize,
    needed: usize,
    truncation: &Truncation,
) -> Option<&'static str> {
    match bound {
        None => Some("unbounded"),
        Some(bound) if needed > bound => Some("bound_exceeded"),
        Some(bound) if bound <= extra => None,
        Some(_) if truncation.vanished() > 0 => Some("vanished"),
        Some(_) => Some("unverified"),
    }
}

/// `bound` as a log field.
pub(super) fn bound_label(bound: Option<usize>) -> String {
    bound.map_or_else(|| "unknown".to_owned(), |bound| bound.to_string())
}

/// `GAMMALOOP_IR_HARD_DEPTH_STRICT=1` (default off): a measured row above
/// the quotient bound (`bound_exceeded`, the proof of `B` failed) is an
/// error instead of a refusal into today's path. Off, the refusal is safe
/// (today's path measures every quotient factor before integrating) and
/// loud: its own `IR_HARD_DEPTH_REFUSED reason=bound_exceeded` line and the
/// `exceeded=` field of `IR_HARD_DEPTH_SUMMARY`.
pub(super) fn strict() -> bool {
    #[cfg(test)]
    if let Some(strict) = TEST.with(|slot| slot.get()).strict {
        return strict;
    }
    static STRICT: std::sync::LazyLock<bool> = std::sync::LazyLock::new(|| {
        matches!(
            std::env::var("GAMMALOOP_IR_HARD_DEPTH_STRICT")
                .as_deref()
                .map(str::trim),
            Ok("1")
        )
    });
    *STRICT
}

/// The fatal [`strict`] alarm, a typed error: every boundary that turns an
/// error into a decline, a missing measurement or a stored message (the
/// closed-path guard, a vanished row's valuation, the counterterm cache) lets
/// it through ([`is_alarm`]), so STRICT is fatal at any nesting depth, with
/// and without `GAMMALOOP_IR_HARD_CHILD_GUARD=1`.
#[derive(Debug)]
pub(super) struct DepthAlarm(String);

impl std::fmt::Display for DepthAlarm {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&self.0)
    }
}

impl std::error::Error for DepthAlarm {}

/// The [`DepthAlarm`] with `message`.
pub(super) fn alarm(message: String) -> eyre::Report {
    eyre::Report::new(DepthAlarm(message))
}

/// A check-mode mismatch that must fail the diagram (MATRIX-R45 review,
/// nested-and-cache finding 2): `IR_DIRAC_CHECK` (the Dirac-closed IR hard
/// part differs from the open one). Fatal through the same boundaries as
/// the [`DepthAlarm`] ([`is_alarm`]): a guarded child attempt never
/// declines it to today's terms.
#[derive(Debug)]
pub(super) struct CheckFailure(String);

impl std::fmt::Display for CheckFailure {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&self.0)
    }
}

impl std::error::Error for CheckFailure {}

/// The [`CheckFailure`] with `message`. On the G-G stack `IR_DIRAC_CHECK`
/// raises its mismatch as a diagnostic failure
/// (`super::ir_dirac::diagnostic_failure`, returned by every guard as well),
/// so this constructor is used by the guard's unit test only.
#[cfg(test)]
pub(super) fn check_failure(message: String) -> eyre::Report {
    eyre::Report::new(CheckFailure(message))
}

/// Whether `error` carries a fatal error anywhere in its chain, also one
/// replayed by the counterterm cache: the [`DepthAlarm`], a
/// [`CheckFailure`], or (U7b) a failure of the uncontracted switch that
/// fails the diagram by design ([`super::uncontracted::Failure`]: a K3 CHECK
/// mismatch, an undesigned refusal, vakint's unknown masters).
pub(super) fn is_alarm(error: &eyre::Report) -> bool {
    error.chain().any(|cause| {
        cause.is::<DepthAlarm>()
            || cause.is::<CheckFailure>()
            || cause.is::<super::uncontracted::Failure>()
            // IRH-DEPTH round 8b (C1): a scalar d at an ε expansion that
            // could not be replaced by 4 − 2ε refuses the diagram.
            || cause.is::<crate::uv::approx::integrated::DimensionAtTruncation>()
            || cause
                .downcast_ref::<crate::uv::counterterm_cache::CachedError>()
                .is_some_and(|cached| is_alarm(cached.cause()))
    })
}

/// IRH-DEPTH counts of one evaluation's child-route IR corrections
/// (`IR_HARD_DEPTH_SUMMARY`).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub(super) struct DepthCounts {
    /// Corrections integrated on the child route (not skipped by selection).
    pub corrections: usize,
    /// Corrections whose quotient bound exceeds the quotient's loops (or is
    /// unknown).
    pub bounded: usize,
    /// Corrections integrated more than once (a depth retry).
    pub retried: usize,
    /// Retries whose depth the bound raised beyond the measured need.
    pub raised: usize,
    /// Corrections of an evaluation on a proper support (a UV child's or a
    /// decorated quotient's R*): nested child calls.
    pub nested: usize,
    /// Nested corrections whose quotient bound exceeds the quotient's loops.
    pub nested_bounded: usize,
    /// Corrections with a measured row above the bound (also refused).
    pub exceeded: usize,
    /// Corrections refused by the guard (`IR_HARD_DEPTH_REFUSED`).
    pub refused: usize,
    /// Refusals by reason (`bound_exceeded` is [`Self::exceeded`]): an
    /// unknown bound, a vanished group or row below the bound, a depth below
    /// the bound with everything shown.
    pub unbounded: usize,
    pub vanished: usize,
    pub unverified: usize,
    /// Refusals by the regularity certificate (IRH-DEPTH round 5,
    /// [`super::regularity`]): a factor of the closed list is not provably
    /// regular, refused before the hard part is integrated.
    pub unresolved_tensor: usize,
    pub opaque_function: usize,
    pub symbol_in_denominator: usize,
    pub outside_list: usize,
    /// A fixed factor that must be regular has a pole (MATRIX-R45 review:
    /// the settings normalisations, the IR Dirac closure normalisation).
    pub pole_in_fixed_factor: usize,
    /// Corrections whose child route itself gave the value
    /// (`outcome=child`), and those of them in a nested evaluation: a
    /// correction counted above can still decline to today's terms later
    /// (`moments:`, `dirac_`, `zero_row:`, `depth:`), so a test that calls a
    /// cell RIGHT needs these (MATRIX-R45 review).
    pub taken: usize,
    pub nested_taken: usize,
    /// Corrections with no shown row whose root projector callback the unit
    /// insertion probe could not evaluate (admitted, certified by contract
    /// only: `IR_HARD_DEPTH_CALLBACK_UNPROBED`).
    pub callback_unprobed: usize,
    /// Tests: the last refusal's reason.
    pub reason: Option<&'static str>,
}

impl DepthCounts {
    /// Tests: fold an inner [`collect`] scope into the enclosing one.
    #[cfg(test)]
    fn add(&mut self, other: &Self) {
        self.corrections += other.corrections;
        self.bounded += other.bounded;
        self.retried += other.retried;
        self.raised += other.raised;
        self.nested += other.nested;
        self.nested_bounded += other.nested_bounded;
        self.exceeded += other.exceeded;
        self.refused += other.refused;
        self.unbounded += other.unbounded;
        self.vanished += other.vanished;
        self.unverified += other.unverified;
        self.unresolved_tensor += other.unresolved_tensor;
        self.opaque_function += other.opaque_function;
        self.symbol_in_denominator += other.symbol_in_denominator;
        self.outside_list += other.outside_list;
        self.pole_in_fixed_factor += other.pole_in_fixed_factor;
        self.taken += other.taken;
        self.nested_taken += other.nested_taken;
        self.callback_unprobed += other.callback_unprobed;
        self.reason = other.reason.or(self.reason);
    }
}

/// Count one event in `counts` (the evaluation's) and, in tests, in the
/// enclosing [`collect`] scope of this thread.
pub(super) fn note(counts: &mut DepthCounts, field: fn(&mut DepthCounts) -> &mut usize) {
    *field(counts) += 1;
    #[cfg(test)]
    COLLECTED.with(|slot| {
        if let Some(collected) = slot.borrow_mut().as_mut() {
            *field(collected) += 1;
        }
    });
}

/// Record a refusal's reason (tests: in the enclosing [`collect`] scope).
pub(super) fn refused(counts: &mut DepthCounts, reason: &'static str) {
    note(counts, |counts| &mut counts.refused);
    match reason {
        "unbounded" => note(counts, |counts| &mut counts.unbounded),
        "vanished" => note(counts, |counts| &mut counts.vanished),
        "unverified" => note(counts, |counts| &mut counts.unverified),
        "unresolved_tensor" => note(counts, |counts| &mut counts.unresolved_tensor),
        "opaque_function" => note(counts, |counts| &mut counts.opaque_function),
        "symbol_in_denominator" => note(counts, |counts| &mut counts.symbol_in_denominator),
        "outside_list" => note(counts, |counts| &mut counts.outside_list),
        "pole_in_fixed_factor" => note(counts, |counts| &mut counts.pole_in_fixed_factor),
        _ => {}
    }
    counts.reason = Some(reason);
    #[cfg(test)]
    COLLECTED.with(|slot| {
        if let Some(collected) = slot.borrow_mut().as_mut() {
            collected.reason = Some(reason);
        }
    });
}

/// Tests: what the planner does (`Fixed` is production).
#[cfg(test)]
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub(super) enum Planner {
    /// The planner of this module (production).
    #[default]
    Fixed,
    /// The fixed planner without the raise to the bound: it starts at the
    /// quotient's loops and retries only to the measured need, so the guard
    /// alone stands between a hidden pole and the value.
    NoRaise,
    /// The planner before IRH-DEPTH (072f9952cb): no raise, no measurement
    /// of vanished rows, no guard (reproduces the missed pole).
    Legacy,
}

/// Tests: the quotient bound the planner sees.
#[cfg(test)]
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub(super) enum BoundOverride {
    #[default]
    Computed,
    /// As if `B` could not be computed.
    Unknown,
    /// `B` replaced by this value (below the truth: `bound_exceeded`).
    Fixed(usize),
    /// `B` replaced by `bound` in nested evaluations only (a decorated
    /// quotient's or a UV child's R*), and only inside `within` of the
    /// top-level evaluation (the graph's root).
    Nested { bound: usize, within: Within },
}

/// Tests: where [`BoundOverride::Nested`] applies. Both scopes are marked by
/// the top-level evaluation only, so a nested evaluation's own attempt (or
/// vanished row) inside the root's fallback to today's terms is outside them.
#[cfg(test)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) enum Within {
    /// A child-route attempt of the top-level evaluation
    /// ([`child_attempt_scope`]; not its fallback to today's terms).
    ChildAttempt,
    /// The valuation of a vanished row of the top-level evaluation
    /// ([`zero_row_scope`]).
    ZeroRow,
}

/// Tests: the planner settings of [`with_depth_test`].
#[cfg(test)]
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub(super) struct DepthTest {
    pub planner: Planner,
    pub bound: BoundOverride,
    /// [`strict`] read as this value.
    pub strict: Option<bool>,
    /// `GAMMALOOP_IR_HARD_CHILD_GUARD` read as this value ([`child_guard`]).
    pub guard: Option<bool>,
}

/// Tests: where a child-route attempt calls the [`RowHook`].
#[cfg(test)]
#[derive(Clone, Debug)]
pub(super) struct HookSite {
    pub hard_loops: usize,
    pub quotient_loops: usize,
    pub extra: usize,
    /// The evaluation is on a proper support (a nested child call).
    pub nested: bool,
    /// The hard part's edges (as the `edges=` log field).
    pub edges: String,
}

/// Tests: a rewrite of the rows `(h_j, E_j, L)` of every child-route attempt,
/// applied before the quotient factors are built (simulates a hard row whose
/// coefficient vanishes at the first depth).
#[cfg(test)]
pub(super) type RowHook = Box<dyn Fn(&HookSite, &mut Vec<(Atom, Atom, usize)>)>;

#[cfg(test)]
thread_local! {
    static ROW_HOOK: std::cell::RefCell<Option<RowHook>> = const { std::cell::RefCell::new(None) };
    static START_EXTRA: std::cell::Cell<Option<usize>> = const { std::cell::Cell::new(None) };
    static TEST: std::cell::Cell<DepthTest> = std::cell::Cell::new(DepthTest::default());
    static COLLECTED: std::cell::RefCell<Option<DepthCounts>> = const { std::cell::RefCell::new(None) };
    /// Child-route attempts running on this thread ([`child_attempt_scope`]).
    static CHILD_ATTEMPT: std::cell::Cell<usize> = const { std::cell::Cell::new(0) };
    /// Vanished rows being valued on this thread ([`zero_row_scope`]).
    static ZERO_ROW: std::cell::Cell<usize> = const { std::cell::Cell::new(0) };
}

/// Tests: run `f` and return the IRH-DEPTH counts of every child-route
/// correction it evaluates on this thread, nested evaluations (UV children
/// computed by R*) included. Per thread and scoped, so tests running in
/// parallel (`cargo test`) do not see each other's counts; an enclosing
/// scope also receives them. The child route itself is selected per thread
/// (`with_ir_hard`), so no child-route correction of `f` runs elsewhere.
#[cfg(test)]
pub(super) fn collect<T>(f: impl FnOnce() -> T) -> (T, DepthCounts) {
    // Restores the enclosing scope (also on a panic), handing it this
    // scope's counts.
    struct Scope(Option<DepthCounts>);
    impl Drop for Scope {
        fn drop(&mut self) {
            let mut outer = self.0.take();
            COLLECTED.with(|slot| {
                let inner = slot.borrow_mut().take().unwrap_or_default();
                if let Some(outer) = outer.as_mut() {
                    outer.add(&inner);
                }
                *slot.borrow_mut() = outer;
            });
        }
    }
    let scope = Scope(COLLECTED.with(|slot| slot.borrow_mut().replace(DepthCounts::default())));
    let value = f();
    let counts = COLLECTED.with(|slot| slot.borrow().unwrap_or_default());
    drop(scope);
    (value, counts)
}

/// Tests: run `f` with the planner settings `test` on this thread.
#[cfg(test)]
pub(super) fn with_depth_test<T>(test: DepthTest, f: impl FnOnce() -> T) -> T {
    struct Restore(DepthTest);
    impl Drop for Restore {
        fn drop(&mut self) {
            let previous = self.0;
            TEST.with(|slot| slot.set(previous));
        }
    }
    let _restore = Restore(TEST.with(|slot| slot.replace(test)));
    f()
}

/// Tests: a counted scope on this thread (restored on unwinding too).
#[cfg(test)]
fn counted_scope(slot: &'static std::thread::LocalKey<std::cell::Cell<usize>>) -> impl Drop {
    struct Scope(&'static std::thread::LocalKey<std::cell::Cell<usize>>);
    impl Drop for Scope {
        fn drop(&mut self) {
            self.0.with(|slot| slot.set(slot.get() - 1));
        }
    }
    slot.with(|slot| slot.set(slot.get() + 1));
    Scope(slot)
}

/// Tests: mark a child-route attempt (not its fallback) on this thread while
/// the returned scope lives ([`Within::ChildAttempt`]).
#[cfg(test)]
pub(super) fn child_attempt_scope() -> impl Drop {
    counted_scope(&CHILD_ATTEMPT)
}

/// Tests: mark the valuation of a vanished row (the quotient factors of a
/// row whose hard coefficient vanished) on this thread while the returned
/// scope lives ([`Within::ZeroRow`]).
#[cfg(test)]
pub(super) fn zero_row_scope() -> impl Drop {
    counted_scope(&ZERO_ROW)
}

/// Tests: run `f` with `hook` rewriting the child route's rows on this thread.
#[cfg(test)]
pub(super) fn with_row_hook<T>(hook: RowHook, f: impl FnOnce() -> T) -> T {
    struct Restore;
    impl Drop for Restore {
        fn drop(&mut self) {
            ROW_HOOK.with(|slot| slot.borrow_mut().take());
        }
    }
    ROW_HOOK.with(|slot| *slot.borrow_mut() = Some(hook));
    let _restore = Restore;
    f()
}

#[cfg(test)]
pub(super) fn apply_row_hook(site: &HookSite, rows: &mut Vec<(Atom, Atom, usize)>) {
    ROW_HOOK.with(|slot| {
        if let Some(hook) = slot.borrow().as_ref() {
            hook(site, rows);
        }
    });
}

/// Tests: run `f` with the child route's first depth at least `extra`.
#[cfg(test)]
pub(super) fn with_start_extra<T>(extra: usize, f: impl FnOnce() -> T) -> T {
    struct Restore(Option<usize>);
    impl Drop for Restore {
        fn drop(&mut self) {
            let previous = self.0;
            START_EXTRA.with(|slot| slot.set(previous));
        }
    }
    let _restore = Restore(START_EXTRA.with(|slot| slot.replace(Some(extra))));
    f()
}

/// The quotient bound the planner uses: `bound` (tests: overridden by
/// [`with_depth_test`]; `nested`: the evaluation is on a proper support).
pub(super) fn planner_bound(bound: Option<usize>, nested: bool) -> Option<usize> {
    #[cfg(test)]
    match TEST.with(|slot| slot.get()).bound {
        BoundOverride::Computed => {}
        BoundOverride::Unknown => return None,
        BoundOverride::Fixed(fixed) => return Some(fixed),
        BoundOverride::Nested {
            bound: fixed,
            within,
        } => {
            let scope = match within {
                Within::ChildAttempt => &CHILD_ATTEMPT,
                Within::ZeroRow => &ZERO_ROW,
            };
            if nested && scope.with(|slot| slot.get()) > 0 {
                return Some(fixed);
            }
        }
    }
    #[cfg(not(test))]
    let _ = nested;
    bound
}

/// Whether the child route runs under the closed path's guard: `env`, the
/// value of `GAMMALOOP_IR_HARD_CHILD_GUARD=1` (tests: overridden by
/// [`with_depth_test`]).
pub(super) fn child_guard(env: bool) -> bool {
    #[cfg(test)]
    if let Some(guard) = TEST.with(|slot| slot.get()).guard {
        return guard;
    }
    env
}

/// Whether the planner raises the depth to the bound (tests: not for
/// [`Planner::NoRaise`] and [`Planner::Legacy`]).
pub(super) fn raises() -> bool {
    #[cfg(test)]
    if TEST.with(|slot| slot.get()).planner != Planner::Fixed {
        return false;
    }
    true
}

/// Whether the planner measures the quotient factors of a vanished row and
/// applies the guard (tests: not for [`Planner::Legacy`]).
pub(super) fn guards() -> bool {
    #[cfg(test)]
    if TEST.with(|slot| slot.get()).planner == Planner::Legacy {
        return false;
    }
    true
}

/// The child route's first depth: the quotient's loops raised to the
/// quotient bound when known (module doc; tests: at least the
/// [`with_start_extra`] value, and without the raise for the planners
/// that do not raise).
pub(super) fn start_extra(quotient_loops: usize, bound: Option<usize>) -> usize {
    let start = match bound {
        Some(bound) if raises() => quotient_loops.max(bound),
        _ => quotient_loops,
    };
    #[cfg(test)]
    if let Some(extra) = START_EXTRA.with(|slot| slot.get()) {
        return start.max(extra);
    }
    start
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::initialisation::test_initialise;
    use spenso::structure::representation::RepName;
    use symbolica::parse;

    fn eps() -> Atom {
        Atom::var(GS.dim_epsilon)
    }

    #[test]
    fn excess_counts_the_negative_orders_of_each_factor() -> eyre::Result<()> {
        test_initialise()?;
        let x = parse!("x");
        assert_eq!(excess_poles(&x), Some(0));
        assert_eq!(excess_poles(&(&x * eps().pow(3))), Some(0));
        assert_eq!(excess_poles(&(&x / eps())), Some(1));
        // Factorwise: the positive order of one factor does not offset the
        // pole of another (the two can land in different K's).
        let mixed = (&x + Atom::one() / eps()) * eps().pow(2) / (Atom::one() + eps());
        assert_eq!(excess_poles(&mixed), Some(1));
        assert_eq!(excess_poles(&(&x / eps().pow(2) / eps())), Some(3));
        // d = 4-2ε: 1/(d-4) is a pole, 1/(d-1) is not.
        let d = Atom::var(GS.dim);
        assert_eq!(excess_poles(&(&x / (&d - Atom::num(1)))), Some(0));
        assert_eq!(excess_poles(&(&x / (&d - Atom::num(4)))), Some(1));
        // Outside the domain, detectable: a non-integer order, a function of
        // ε (symbolica would expand an opaque f(ε) as regular).
        assert_eq!(
            excess_poles(&(&x * eps().pow(Atom::num(1) / Atom::num(2)))),
            None
        );
        let opaque = symbolica::function!(symbolica::symbol!("irh_depth_test_f"), eps());
        assert_eq!(excess_poles(&(&x * &opaque)), None);
        assert_eq!(excess_poles(&(&x * opaque / eps())), None);
        // A representation's dimension is no function of ε.
        let mink = spenso::structure::representation::Minkowski {}.new_rep(GS.dim);
        let g = spenso::network::library::symbolic::ETS.metric(
            mink.to_symbolic([Atom::num(1)]),
            mink.to_symbolic([Atom::num(2)]),
        );
        assert_eq!(excess_poles(&(&g * &d / eps())), Some(1));
        Ok(())
    }

    #[test]
    fn scoped_counts_see_only_their_own_scope() {
        let mut counts = DepthCounts::default();
        let ((), outer) = collect(|| {
            note(&mut counts, |counts| &mut counts.corrections);
            let ((), inner) = collect(|| {
                note(&mut counts, |counts| &mut counts.retried);
                refused(&mut counts, "vanished");
            });
            assert_eq!(
                (
                    inner.corrections,
                    inner.retried,
                    inner.refused,
                    inner.reason
                ),
                (0, 1, 1, Some("vanished"))
            );
        });
        assert_eq!(
            (
                outer.corrections,
                outer.retried,
                outer.refused,
                outer.reason
            ),
            (1, 1, 1, Some("vanished"))
        );
        // Outside every scope only the evaluation's own counts change.
        note(&mut counts, |counts| &mut counts.raised);
        assert_eq!(
            (counts.corrections, counts.raised, counts.refused),
            (1, 1, 1)
        );
        let ((), empty) = collect(|| ());
        assert_eq!(empty, DepthCounts::default());
    }

    #[test]
    fn the_bound_adds_the_quotient_loops_and_every_excess() -> eyre::Result<()> {
        test_initialise()?;
        let x = parse!("x");
        assert_eq!(quotient_pole_bound(1, &x, None), Some(1));
        assert_eq!(quotient_pole_bound(1, &(&x / eps()), None), Some(2));
        let closure = Atom::one() / (Atom::var(GS.dim) - Atom::num(1));
        assert_eq!(quotient_pole_bound(2, &x, Some(&closure)), Some(2));
        assert_eq!(row_bound(Some(2), &(&x / eps())), Some(3));
        assert_eq!(row_bound(None, &x), None);
        // The dropped rows: bounded by B only while every returned insertion
        // is free of ε (correctness review, finding 1).
        assert_eq!(dropped_row_bound(Some(2), [Some(0), Some(0)]), Some(2));
        assert_eq!(dropped_row_bound(Some(2), []), Some(2));
        assert_eq!(dropped_row_bound(Some(2), [Some(0), Some(1)]), None);
        assert_eq!(dropped_row_bound(Some(2), [Some(0), None]), None);
        assert_eq!(dropped_row_bound(None, [Some(0)]), None);
        assert_eq!(correction_bound(Some(1), &[Some(1), Some(3)]), Some(3));
        assert_eq!(correction_bound(Some(1), &[]), Some(1));
        assert_eq!(correction_bound(None, &[Some(1)]), None);
        assert_eq!(correction_bound(Some(1), &[Some(1), None]), None);
        // Row by row (correctness review, finding 2): a row at its bound B
        // measuring B + 1 is flagged although another row's bound is higher.
        assert_eq!(exceeded_row(&[2, 3], &[Some(1), Some(3)]), Some((2, 1)));
        assert_eq!(exceeded_row(&[1, 3], &[Some(1), Some(3)]), None);
        assert_eq!(exceeded_row(&[5], &[None]), None);
        Ok(())
    }

    /// A propagator-insertion surrogate's coefficient vertex, as
    /// `insertions.rs` builds it at run time (`with_coefficient` of the
    /// normalised insertion): an ε polynomial with rational and ζ
    /// coefficients, rational functions of d regular at d = 4, and the
    /// insertion's tensor structure. It is read in full (production review,
    /// finding 2) and has no excess: B stays the quotient's loops.
    #[test]
    fn a_surrogate_coefficient_has_no_excess() -> eyre::Result<()> {
        test_initialise()?;
        let d = Atom::var(GS.dim);
        let z3 = parse!("zeta(3)");
        let mink = spenso::structure::representation::Minkowski {}.new_rep(GS.dim);
        let g = spenso::network::library::symbolic::ETS.metric(
            mink.to_symbolic([Atom::num(1)]),
            mink.to_symbolic([Atom::num(2)]),
        );
        let p2 = parse!("p2");
        let series = q(-5, 3) + q(31, 9) * eps() + (q(7, 2) - Atom::num(6) * &z3) * eps().pow(2);
        let coefficient = (&d - Atom::num(2)) / (&d - Atom::num(1))
            * (Atom::num(3) * &d - Atom::num(2))
            / (&d + Atom::num(2))
            * series
            * (&g * &p2);
        assert_eq!(excess_poles(&coefficient), Some(0));
        assert_eq!(quotient_pole_bound(2, &coefficient, None), Some(2));
        // The same coefficient with a pole left in it is counted.
        assert_eq!(excess_poles(&(&coefficient / eps())), Some(1));
        Ok(())
    }

    fn q(n: i64, d: i64) -> Atom {
        Atom::num(n) / Atom::num(d)
    }

    #[test]
    fn the_guard_refuses_only_a_depth_below_the_bound() {
        let shown = Truncation {
            groups: 3,
            vanished_groups: 0,
            rows: 18,
            zero_rows: 0,
        };
        let vanished = Truncation {
            vanished_groups: 3,
            rows: 0,
            ..shown
        };
        let zero_row = Truncation {
            zero_rows: 1,
            ..shown
        };
        // Production: B = quotient loops = extra.
        assert_eq!(refusal(Some(1), 1, 1, &shown), None);
        assert_eq!(refusal(Some(1), 1, 0, &vanished), None);
        // L3_006 with v0·ε³, v1·ε⁻¹: every group vanishes at extra 1, B = 2.
        assert_eq!(refusal(Some(2), 1, 0, &vanished), Some("vanished"));
        assert_eq!(refusal(Some(2), 1, 1, &zero_row), Some("vanished"));
        assert_eq!(refusal(Some(2), 1, 1, &shown), Some("unverified"));
        // After the retry to the bound.
        assert_eq!(refusal(Some(2), 2, 2, &vanished), None);
        assert_eq!(refusal(Some(2), 3, 3, &shown), Some("bound_exceeded"));
        // An unknown bound refuses whatever the truncation shows: a row of a
        // nonzero group dropped by the moment extraction is in no count.
        assert_eq!(refusal(None, 1, 1, &shown), Some("unbounded"));
        assert_eq!(refusal(None, 1, 1, &zero_row), Some("unbounded"));
        assert_eq!(refusal(None, 5, 0, &shown), Some("unbounded"));
    }
}
