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
//!   adds its gamma structure. Its own order is counted per row.
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
//! The guard stays as a backstop: it refuses the child route for a
//! correction whose accepted depth stays below `B` (only possible when `B`
//! is unknown, or when a measured row exceeds `B`); the correction is then
//! evaluated on today's path (the closed or projected moments with the
//! generic forest step, which measures every quotient factor before its
//! hard coefficient is integrated).

use std::sync::atomic::{AtomicUsize, Ordering};

use symbolica::{
    atom::{Atom, AtomCore, AtomView},
    poly::series::SeriesDepth,
};

use crate::{utils::GS, uv::approx::integrated::Integrated};

/// `Σ max(0, -ord_ε f)` over the product factors `f` of `expression`, read
/// at `d = 4-2ε`: how many poles beyond its loops a `K` of a region carrying
/// these factors can gain. `Some(0)` without `ε` (or `d`); `None` when a
/// factor's order is not an integer or cannot be determined.
pub(super) fn excess_poles(expression: &Atom) -> Option<usize> {
    let epsilon = GS.dim_epsilon;
    if !expression.contains_symbol(epsilon) && !expression.contains_symbol(GS.dim) {
        return Some(0);
    }
    let regulated = Integrated::dimensionally_regularized(expression);
    if !regulated.contains_symbol(epsilon) {
        return Some(0);
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

/// The bound of one row: `bound` raised by the excess of its insertion.
pub(super) fn row_bound(bound: Option<usize>, insertion: &Atom) -> Option<usize> {
    Some(bound? + excess_poles(insertion)?)
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
/// * `bound_exceeded`: a measured row needs more than `B`; the bound does not
///   hold, so nothing it certifies can be trusted.
/// * `vanished`: a vanished group or row may carry up to `B > extra` poles.
/// * `unverified`: `B > extra` with every group and row shown; a row of a
///   nonzero group whose coefficient vanished at the truncation is dropped
///   by the moment extraction and cannot be told apart from an absent one.
/// * `unbounded`: `B` is unknown and something vanished.
pub(super) fn refusal(
    bound: Option<usize>,
    extra: usize,
    needed: usize,
    truncation: &Truncation,
) -> Option<&'static str> {
    match bound {
        None if truncation.vanished() > 0 => Some("unbounded"),
        None => None,
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

/// IRH-DEPTH counts of one evaluation's child-route IR corrections
/// (`IR_HARD_CHILD_DEPTH_SUMMARY`).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub(super) struct DepthCounts {
    /// Corrections integrated on the child route (not skipped by selection).
    pub corrections: usize,
    /// Corrections whose quotient bound exceeds the quotient's loops.
    pub bounded: usize,
    /// Corrections integrated more than once (a depth retry).
    pub retried: usize,
    /// Retries whose depth the bound raised beyond the measured need.
    pub raised: usize,
    /// Corrections refused by the guard (`IR_HARD_CHILD_DEPTH_REFUSED`).
    pub refused: usize,
}

/// Process-wide totals (tests, and every refusal also has its own line).
pub(super) static CORRECTIONS: AtomicUsize = AtomicUsize::new(0);
pub(super) static RETRIED: AtomicUsize = AtomicUsize::new(0);
pub(super) static RAISED: AtomicUsize = AtomicUsize::new(0);
pub(super) static REFUSED: AtomicUsize = AtomicUsize::new(0);

pub(super) fn count(counter: &AtomicUsize) {
    counter.fetch_add(1, Ordering::Relaxed);
}

/// Tests: the process-wide `(corrections, retried, raised, refused)`.
#[cfg(test)]
pub(super) fn totals() -> [usize; 4] {
    [&CORRECTIONS, &RETRIED, &RAISED, &REFUSED].map(|counter| counter.load(Ordering::Relaxed))
}

/// Tests: where a child-route attempt calls the [`RowHook`].
#[cfg(test)]
#[derive(Clone, Copy, Debug)]
pub(super) struct HookSite {
    pub hard_loops: usize,
    pub quotient_loops: usize,
    pub extra: usize,
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

/// The child route's first depth: the quotient's loops (tests: at least the
/// [`with_start_extra`] value).
pub(super) fn start_extra(quotient_loops: usize) -> usize {
    #[cfg(test)]
    if let Some(extra) = START_EXTRA.with(|slot| slot.get()) {
        return quotient_loops.max(extra);
    }
    quotient_loops
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::initialisation::test_initialise;
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
        Ok(())
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
        Ok(())
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
        assert_eq!(refusal(None, 1, 1, &shown), None);
        assert_eq!(refusal(None, 1, 1, &zero_row), Some("unbounded"));
    }
}
