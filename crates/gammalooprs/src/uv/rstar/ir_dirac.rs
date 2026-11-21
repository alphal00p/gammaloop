//! `GAMMALOOP_IR_DIRAC_CLOSURE=on|check` (unset: `on`): the Dirac closure of IR
//! hard parts with one open Dirac line, as in FORM (projector -> trace).
//!
//! The hard part `H` of an IR correction with one open Dirac line (row `r`,
//! column `c`) whose vacuum average has a unique Dirac structure is computed
//! as a closed trace instead of the open matrix:
//!
//! * hard Taylor order `omega = 0` with one open Lorentz slot `mu` (`vertex`):
//!   a constant matrix with one Lorentz index is `c gamma^mu` in `d`
//!   dimensions, so `<H^mu> = Tr[gamma_mu <H^mu>] / (Tr[1] d) gamma^mu`;
//! * `omega = 1` with no Lorentz slot (`slash`): linear in the boundary
//!   momenta, `<H> = sum_e c_e slash(Q_e) = Tr[gamma_nu <H>] / Tr[1] gamma^nu`
//!   with a fresh index `nu`, provided the line is odd in the gamma matrices
//!   (an even line would add a unit matrix at order 0; such hard parts decline).
//!
//! The closure `gamma(c, r, mu)` is momentum free: it is multiplied into every
//! routed hard remainder before the soft Taylor expansion (it commutes with the
//! expansion, the contraction, the vacuum projection and the integration), so
//! the contraction produces a closed trace and the closed scalar paths take it.
//! Each hard moment `(h_j, E_j)` of the closed hard part is then reassembled as
//! `(h_j N, E_j gamma(r, c, mu))` (vertex, `N = 1/(Tr[1] d)`) or
//! `(h_j N, gamma(r, c, Q_e))` for `E_j = Q_e^nu` (slash, `N = 1/Tr[1]`). `d`
//! is regular and stays in the coefficient, outside the quotient's `K`, as the
//! closed path's Gram inverse does; the insertion is a plain gamma matrix, so a
//! Dirac-closed quotient (`GAMMALOOP_UV_CHILD_CLOSURE=dirac`) can close too.
//!
//! Any other hard part with spinor slots (two Dirac segments, Lorentz slots at
//! `omega >= 1`, colour slots, parity-odd content) declines to the open path
//! with a logged reason. `off` (or empty, `0`): the open path, unchanged.

use std::sync::LazyLock;

use eyre::{Result, eyre};
use idenso::{
    CookMode, CookSettings,
    color::ColorSimplifySettings,
    dirac::GammaSimplifySettings,
    tensor::{AlgebraSettings, SymbolicTensor},
};
use linnet::half_edge::subgraph::SuBitGraph;
use spenso::structure::partial::PartialStructureExt;
#[cfg(test)]
use spenso::structure::representation::{Minkowski, RepName};
use symbolica::{
    atom::{Atom, AtomCore, AtomView, Symbol},
    function, symbol,
};

use super::expression::{DiracBoundary, DiracKind, RStarIntegrand};
use crate::{graph::Graph, utils::GS, uv::approx::integrated::Integrated};

/// `GAMMALOOP_IR_DIRAC_CLOSURE`.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum IrDiracMode {
    /// Empty, `0` or `off`: the open hard part (the path before the
    /// production default).
    Off,
    /// `on`, `1` or unset: close the shapes with a unique Dirac structure.
    On,
    /// `check`: as `on`, and also compute the open hard moments, which must
    /// equal the closed, reassembled ones (`IR_DIRAC_CHECK`, an error
    /// otherwise; on the child route through the moments' justified order,
    /// [`agree`]). The closed moments are used.
    Check,
}

impl IrDiracMode {
    pub fn parse(value: Option<&str>) -> Result<Self> {
        match value.map(str::trim) {
            Some("" | "0" | "off") => Ok(Self::Off),
            None | Some("1" | "on") => Ok(Self::On),
            Some("check") => Ok(Self::Check),
            Some(other) => Err(eyre!(
                "GAMMALOOP_IR_DIRAC_CLOSURE must be unset, `off`, `on` or `check`, not `{other}`"
            )),
        }
    }

    pub fn label(self) -> &'static str {
        match self {
            Self::Off => "off",
            Self::On => "on",
            Self::Check => "check",
        }
    }
}

static MODE: LazyLock<std::result::Result<IrDiracMode, String>> =
    LazyLock::new(|| match std::env::var("GAMMALOOP_IR_DIRAC_CLOSURE") {
        Ok(value) => IrDiracMode::parse(Some(&value)).map_err(|error| format!("{error}")),
        Err(std::env::VarError::NotPresent) => Ok(IrDiracMode::On),
        Err(std::env::VarError::NotUnicode(value)) => Err(format!(
            "GAMMALOOP_IR_DIRAC_CLOSURE={value:?} is not unicode"
        )),
    });

#[cfg(test)]
thread_local! {
    static TEST_MODE: std::cell::Cell<Option<IrDiracMode>> = const { std::cell::Cell::new(None) };
}

/// The mode of `GAMMALOOP_IR_DIRAC_CLOSURE`, read once; an invalid value is
/// an error at first use.
pub fn ir_dirac_mode() -> Result<IrDiracMode> {
    #[cfg(test)]
    if let Some(mode) = TEST_MODE.with(|slot| slot.get()) {
        return Ok(mode);
    }
    MODE.clone().map_err(|error| eyre!(error))
}

/// Tests: run `f` with the mode set on this thread, never through the
/// environment.
#[cfg(test)]
pub(crate) fn with_ir_dirac_mode<T>(mode: IrDiracMode, f: impl FnOnce() -> T) -> T {
    struct Restore(Option<IrDiracMode>);
    impl Drop for Restore {
        fn drop(&mut self) {
            let previous = self.0;
            TEST_MODE.with(|slot| slot.set(previous));
        }
    }
    let _restore = Restore(TEST_MODE.with(|slot| slot.replace(Some(mode))));
    f()
}

#[cfg(test)]
thread_local! {
    static TEST_CHILD_GUARD: std::cell::Cell<Option<bool>> = const { std::cell::Cell::new(None) };
    static TEST_MUTATION: std::cell::Cell<Option<(usize, i64)>> = const { std::cell::Cell::new(None) };
}

/// Tests: run `f` with `GAMMALOOP_IR_HARD_CHILD_GUARD` (the closed path's
/// scoped budget around an IR hard part taking the child route) set to
/// `guard` on this thread, never through the environment.
#[cfg(test)]
pub(crate) fn with_child_guard<T>(guard: bool, f: impl FnOnce() -> T) -> T {
    struct Restore(Option<bool>);
    impl Drop for Restore {
        fn drop(&mut self) {
            let previous = self.0;
            TEST_CHILD_GUARD.with(|slot| slot.set(previous));
        }
    }
    let _restore = Restore(TEST_CHILD_GUARD.with(|slot| slot.replace(Some(guard))));
    f()
}

/// Tests: the guard set by [`with_child_guard`] on this thread.
#[cfg(test)]
pub(super) fn test_child_guard() -> Option<bool> {
    TEST_CHILD_GUARD.with(std::cell::Cell::get)
}

/// Tests: run `f` with a controlled change of the Dirac-closed IR hard part
/// `H(γ)` of the child route (the closed arm of `IR_DIRAC_CHECK`, used
/// downstream too) on this thread: the group of `loops` remainder loops
/// becomes `H (1 + ε^(order - t) / 7)`, `t` its leading Laurent order, so
/// its `ε^order` coefficient changes (by a seventh of the leading one) and
/// no coefficient below `ε^order`.
#[cfg(test)]
pub(crate) fn with_closed_mutation<T>(loops: usize, order: i64, f: impl FnOnce() -> T) -> T {
    struct Restore(Option<(usize, i64)>);
    impl Drop for Restore {
        fn drop(&mut self) {
            let previous = self.0;
            TEST_MUTATION.with(|slot| slot.set(previous));
        }
    }
    let _restore = Restore(TEST_MUTATION.with(|slot| slot.replace(Some((loops, order)))));
    f()
}

/// Tests: the closed `H(γ)` group `h` of `loops` remainder loops after the
/// change of [`with_closed_mutation`] (unchanged without one, or for another
/// loop count). Counted as `mutated`.
#[cfg(test)]
pub(super) fn mutate_closed(loops: usize, h: &Atom) -> Result<Atom> {
    let Some((target, order)) = TEST_MUTATION.with(std::cell::Cell::get) else {
        return Ok(h.clone());
    };
    if target != loops || h.is_zero() {
        return Ok(h.clone());
    }
    let series = h.series(
        GS.dim_epsilon,
        0,
        symbolica::poly::series::SeriesDepth::relative(1),
    )?;
    let exponent = series.get_trailing_exponent();
    let leading = exponent
        .is_integer()
        .then(|| exponent.numerator().to_i64())
        .flatten()
        .ok_or_else(|| eyre!("mutate_closed: a non-integer leading order {exponent}"))?;
    let shift = order - leading;
    eyre::ensure!(
        shift >= 0,
        "mutate_closed: the order {order} is below the leading order {leading}"
    );
    count("mutated".into());
    let factor = Atom::var(GS.dim_epsilon).pow(Atom::num(shift)) / Atom::num(7);
    Ok((h + &(h * &factor)).expand())
}

/// The Δ-cache key of an R* child under `mode`: unchanged when off, else
/// with `#ir-dirac-closure=<mode>` (the child's IR hard parts are closed).
pub(super) fn child_key(mut key: String, mode: IrDiracMode) -> String {
    if mode != IrDiracMode::Off {
        key.push_str("#ir-dirac-closure=");
        key.push_str(mode.label());
    }
    key
}

/// The fresh index `nu` of the slash closure of an IR hard part: distinct from
/// the UV child's ([`super::expression::DIRAC_CLOSURE_INDEX`]), and never left
/// in an insertion ([`reassemble`] writes `gamma(r, c, Q_e)`).
pub(super) static IR_DIRAC_INDEX: LazyLock<Symbol> =
    LazyLock::new(|| symbol!("gammalooprs::ir_dirac_closure_nu"));

/// The Dirac closure of the IR hard part `hard` at Taylor order `omega`, or
/// the decline reason.
pub(super) fn boundary(
    graph: &Graph,
    hard: &SuBitGraph,
    omega: i32,
) -> Result<std::result::Result<DiracBoundary, String>> {
    RStarIntegrand::dirac_boundary_with(graph, hard, omega, *IR_DIRAC_INDEX, true)
}

/// `vertex`, `slash` or `none` (the `kind=` of an `IR_DIRAC` line).
pub(super) fn kind_label(boundary: Option<&DiracBoundary>) -> &'static str {
    match boundary.map(|boundary| boundary.kind) {
        Some(DiracKind::Vertex) => "vertex",
        Some(DiracKind::SelfEnergy) => "slash",
        None => "none",
    }
}

fn cooking() -> CookSettings {
    CookSettings::indices()
        .with_mode(CookMode::ReversibleEncoding)
        .with_representation_payloads(true, true)
}

/// The number of open slots of `expression`.
fn open_slot_count(expression: &Atom) -> std::result::Result<usize, String> {
    let cooked = cooking()
        .try_cook(expression.as_view())
        .map_err(|error| format!("cook:{error}"))?;
    let tensor = SymbolicTensor::infer(cooked).map_err(|error| format!("infer:{error}"))?;
    Ok(tensor.structure().logical_slots().len())
}

/// `Q(e, mink(D, nu))` -> `gamma(r, c, Q(e, mink(D)))` for the closure index
/// `nu` of `boundary`; `None` when `nu` remains anywhere else.
fn slash_insertion(insertion: &Atom, boundary: &DiracBoundary) -> Option<Atom> {
    let AtomView::Fun(index) = boundary.index.as_view() else {
        return None;
    };
    let representation = index.get_symbol();
    let label = index.iter().last()?.to_owned();
    let mut replaced = 0_usize;
    let result = insertion.replace_map(|view, _, out| {
        let AtomView::Fun(momentum) = view else {
            return;
        };
        if momentum.get_symbol() != GS.emr_mom || momentum.get_nargs() != 2 {
            return;
        }
        let AtomView::Fun(slot) = momentum.get(1) else {
            return;
        };
        if slot.get_symbol() != representation || slot.iter().last() != Some(label.as_view()) {
            return;
        }
        // The bare vector: the slot's representation without its index.
        let arguments = slot.iter().collect::<Vec<_>>();
        let mut bare = symbolica::atom::FunctionBuilder::new(representation);
        for argument in &arguments[..arguments.len() - 1] {
            bare = bare.add_arg(*argument);
        }
        let vector = function!(GS.emr_mom, momentum.get(0), bare.finish());
        replaced += 1;
        **out = function!(
            idenso::dirac::AGS.gamma,
            &boundary.row,
            &boundary.column,
            vector
        );
    });
    let left = {
        let mut found = false;
        result.visitor(&mut |view| {
            found |= view == label.as_view();
            !found
        });
        found
    };
    (replaced > 0 && !left).then_some(result)
}

/// The closed hard moments `(h_j, E_j)` reassembled with the open Dirac
/// structure: `(h_j N, E_j gamma(r, c, mu))` for a vertex (each `E_j` free of
/// open slots), `(h_j N, gamma(r, c, Q_e))` for `E_j = Q_e^nu` (slash), with
/// `d` regulated in the coefficient. A row of any other form declines.
pub(super) fn reassemble(
    rows: &[(Atom, Atom)],
    boundary: &DiracBoundary,
) -> std::result::Result<Vec<(Atom, Atom)>, String> {
    rows.iter()
        .map(|(coefficient, insertion)| {
            let insertion = match boundary.kind {
                DiracKind::Vertex => match open_slot_count(insertion)? {
                    0 => insertion * &boundary.structure,
                    count => return Err(format!("reassembly:open_slots:{count}")),
                },
                DiracKind::SelfEnergy => match open_slot_count(insertion)? {
                    1 => slash_insertion(insertion, boundary)
                        .ok_or_else(|| "reassembly:slash_form".to_owned())?,
                    count => return Err(format!("reassembly:open_slots:{count}")),
                },
            };
            let coefficient =
                Integrated::dimensionally_regularized(&(coefficient * &boundary.normalization));
            Ok((coefficient, insertion))
        })
        .collect()
}

/// Whether the tensor expression `difference` vanishes: by expansion, else
/// after gamma and colour simplification (chains and traces undone) and
/// tensor canonization, as `RStarPlan::dirac_closure_agrees` compares.
fn tensor_zero(difference: &Atom) -> Result<bool> {
    let difference = Integrated::dimensionally_regularized(difference).expand();
    if difference.is_zero() {
        return Ok(true);
    }
    let cooking = cooking();
    let tensor = SymbolicTensor::infer_indexed(cooking.try_cook(difference.as_view())?)?
        .undo_trace()?
        .undo_chain()?
        .simplify_algebra(&AlgebraSettings {
            gamma: Some(GammaSimplifySettings::default()),
            color: Some(ColorSimplifySettings::default()),
            ..Default::default()
        })?
        .to_dots()?;
    let reduced =
        Integrated::dimensionally_regularized(&cooking.uncook(tensor.expression().as_view()))
            .expand();
    if reduced.is_zero() {
        return Ok(true);
    }
    let canonical = SymbolicTensor::infer(cooking.try_cook(reduced.as_view())?)?.canonize()?;
    let canonical = cooking.uncook(canonical.expression().as_view()).expand();
    Ok(canonical.together().is_zero())
}

/// Whether `difference` vanishes: at the check's random rational point, else
/// by [`tensor_zero`].
fn vanishes(difference: &Atom) -> Result<bool> {
    if difference.is_zero() || super::closed::random_point(difference).is_zero() {
        return Ok(true);
    }
    tensor_zero(difference)
}

/// The check of the closed moments against the open ones of one hard part.
#[derive(Debug, Default, Clone, Copy, PartialEq, Eq)]
pub(super) struct Agreement {
    /// The reassembled rows are the open rows, atom by atom.
    pub identical: bool,
    /// Every open insertion is a multiple of the closed structure (the unique
    /// structure the closure assumes), or its row vanishes through
    /// [`Self::through`].
    pub structure: bool,
    /// The closed coefficients equal the open ones projected on the structure
    /// (through [`Self::through`]).
    pub coefficients: bool,
    /// The Laurent order `ε^k` through which structure and coefficients are
    /// compared: the hard part's justified depth (the child route, `ε^(extra
    /// - 1)` as the K3 records' `valid_through`); `None`: in full.
    pub through: Option<i64>,
    /// The whole projected difference vanishes (the comparison in full;
    /// informational when [`Self::through`] is set).
    pub full: bool,
    /// The lowest order at which the projected coefficients differ, scanned
    /// up to [`FIRST_DIFFERENCE_SCAN`] orders above [`Self::through`] when it
    /// is set; `None`: no difference found (or not scanned).
    pub first_difference: Option<i64>,
}

impl Agreement {
    pub fn equal(&self) -> bool {
        self.identical || (self.structure && self.coefficients)
    }

    /// `eps^k`, or `all` for a comparison in full.
    pub fn through_label(&self) -> String {
        self.through
            .map_or_else(|| "all".to_owned(), |order| format!("eps^{order}"))
    }

    /// `eps^k`, or `none`.
    pub fn first_difference_label(&self) -> String {
        self.first_difference
            .map_or_else(|| "none".to_owned(), |order| format!("eps^{order}"))
    }
}

/// The orders scanned above the justified one for the first difference.
const FIRST_DIFFERENCE_SCAN: i64 = 8;

/// The Laurent part of `value` (`d` regulated) through `ε^order`.
fn laurent_through(value: &Atom, order: i64) -> Result<Atom> {
    Ok(Integrated::dimensionally_regularized(value)
        .series(GS.dim_epsilon, 0, order)?
        .to_atom()
        .expand())
}

/// Whether the projected difference `difference` vanishes through
/// `ε^through` (in full for `None`, as before the justified order), whether
/// it vanishes in full, and its lowest differing order (scanned only for a
/// set `through` whose full difference does not vanish).
fn coefficients_agree(
    difference: &Atom,
    through: Option<i64>,
) -> Result<(bool, bool, Option<i64>)> {
    let full = vanishes(difference)?;
    let Some(order) = through else {
        return Ok((full, full, None));
    };
    if full {
        return Ok((true, true, None));
    }
    let regulated = Integrated::dimensionally_regularized(difference).expand();
    let series = regulated.series(
        GS.dim_epsilon,
        0,
        symbolica::poly::series::SeriesDepth::relative(1),
    )?;
    let exponent = series.get_trailing_exponent();
    let start = exponent
        .is_integer()
        .then(|| exponent.numerator().to_i64())
        .flatten()
        .map_or(order - FIRST_DIFFERENCE_SCAN, |leading| {
            leading.min(order + 1)
        });
    for k in start..=order + FIRST_DIFFERENCE_SCAN {
        if !vanishes(&laurent_through(&regulated, k)?)? {
            return Ok((k > order, false, Some(k)));
        }
    }
    Ok((true, false, None))
}

/// Compare the closed moments `closed` (before [`reassemble`]) and their
/// reassembly `reassembled` with the open moments `open` of the same hard part
/// (`sum_j h_j E_j`). Each open insertion `E_o` is projected on the structure,
/// `t = N Tr[closure E_o]` (a scalar, or a vector in `nu`), and must equal
/// `t * structure` (unique structure); then `sum h_o t` must equal
/// `sum h_c N E_c` (the closed coefficients). With `through = Some(k)` both
/// through `ε^k`, the justified depth of the moments (an open row off the
/// structure passes when `h_o (E_o - t structure)` vanishes through `ε^k`;
/// the projected coefficients' difference is truncated at `ε^k`): the orders
/// above it are incomplete in either arm (the K3 records stop at their
/// `valid_through`, today's integration carries more), and never reach the
/// counterterm. `None`: both in full (today's path, unchanged).
pub(super) fn agree(
    closed: &[(Atom, Atom)],
    reassembled: &[(Atom, Atom)],
    open: &[(Atom, Atom)],
    boundary: &DiracBoundary,
    through: Option<i64>,
) -> Result<Agreement> {
    let identical = {
        let key = |rows: &[(Atom, Atom)]| {
            let mut keys = rows
                .iter()
                .map(|(h, e)| (h.to_canonical_string(), e.to_canonical_string()))
                .collect::<Vec<_>>();
            keys.sort();
            keys
        };
        key(reassembled) == key(open)
    };
    let normalization = Integrated::dimensionally_regularized(&boundary.normalization);
    let mut structure = true;
    let mut projected_open = Atom::Zero;
    for (coefficient, insertion) in open {
        let trace = RStarIntegrand::simplify(&(insertion * &boundary.closure))?;
        let projection = Integrated::dimensionally_regularized(&(trace * &normalization));
        if structure {
            let residual = insertion - &projection * &boundary.structure;
            structure = tensor_zero(&residual)?
                || match through {
                    Some(order) => {
                        tensor_zero(&laurent_through(&(coefficient * &residual), order)?)?
                    }
                    None => false,
                };
        }
        projected_open += coefficient * projection;
    }
    let projected_closed = closed
        .iter()
        .fold(Atom::Zero, |sum, (coefficient, insertion)| {
            sum + Integrated::dimensionally_regularized(&(coefficient * &normalization)) * insertion
        });
    let (coefficients, full, first_difference) =
        coefficients_agree(&(projected_closed - projected_open), through)?;
    Ok(Agreement {
        identical,
        structure,
        coefficients,
        through,
        full,
        first_difference,
    })
}

/// A failure no fallback may decline: a failed identity comparison
/// (`IR_DIRAC_CHECK`, `UNCONTRACTED_CHECK`, `DIRAC_CLOSURE_CHECK`,
/// `INSERTION_CLOSURE_CHECK`, `BOUNDARY_CHECK`, `GAMMALOOP_IR_OUTER=check`) or
/// a record refusal that is not designed in a strict or check run. Its text
/// is the failure's. The Δ-cache keeps it one on every reuse of a stored
/// failure (`uv::counterterm_cache`): a decorated quotient's R* runs
/// inside the child guard.
#[derive(Debug)]
pub(super) struct DiagnosticFailure(String);

impl std::fmt::Display for DiagnosticFailure {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&self.0)
    }
}

impl std::error::Error for DiagnosticFailure {}

/// `report` (its whole text) as a [`DiagnosticFailure`].
pub(crate) fn diagnostic_failure(report: eyre::Report) -> eyre::Report {
    eyre::Report::new(DiagnosticFailure(format!("{report:#}")))
}

/// Whether `report` carries a [`DiagnosticFailure`].
pub(crate) fn is_diagnostic_failure(report: &eyre::Report) -> bool {
    report.chain().any(|cause| cause.is::<DiagnosticFailure>())
}

/// [`super::closed::guarded`] (`GAMMALOOP_IR_HARD_CHILD_GUARD=1`, the scoped
/// budget around an IR hard part on the child route): an exhausted budget, a
/// panic, an ordinary error or a decline of `work` declines (`Ok(Err)`), as
/// before; a [`DiagnosticFailure`] is returned as the error, so a failed
/// comparison fails the diagram whether the guard is on or off.
pub(super) fn guarded_unless_diagnostic<T>(
    limit: Option<u64>,
    release: Option<fn()>,
    work: impl FnOnce() -> Result<std::result::Result<T, String>>,
) -> Result<std::result::Result<T, String>> {
    let mut failure = None;
    let outcome = super::closed::guarded(limit, release, || match work() {
        Err(error) if is_diagnostic_failure(&error) => {
            failure = Some(error);
            Ok(Err("diagnostic".to_owned()))
        }
        other => other,
    });
    match failure {
        Some(error) => Err(error),
        None => Ok(outcome),
    }
}

/// The bare Minkowski vector `Q(e, mink(d))` of soft edge `edge`.
#[cfg(test)]
pub(super) fn bare_vector(edge: usize) -> Atom {
    function!(
        GS.emr_mom,
        edge,
        Minkowski {}.new_rep(GS.dim).to_symbolic([])
    )
}

#[cfg(test)]
thread_local! {
    static COUNTS: std::cell::RefCell<std::collections::BTreeMap<String, usize>> =
        const { std::cell::RefCell::new(std::collections::BTreeMap::new()) };
}

/// Tests: count one outcome on this thread (`<path>:<kind>:<outcome>`, with
/// the decline reason cut at its first colon-separated field).
#[cfg(test)]
pub(super) fn count(key: String) {
    COUNTS.with(|counts| *counts.borrow_mut().entry(key).or_default() += 1);
}

/// Tests: count one `IR_DIRAC_CHECK` (`check`, `check_identical`,
/// `check_equal`).
#[cfg(test)]
pub(super) fn count_check(agreement: Agreement) {
    count("check".into());
    if agreement.identical {
        count("check_identical".into());
    }
    if agreement.equal() {
        count("check_equal".into());
    }
    count(format!("check_through:{}", agreement.through_label()));
    count(format!(
        "check_first:{}",
        agreement.first_difference_label()
    ));
}

/// Tests: the outcomes counted on this thread since the last call.
#[cfg(test)]
pub(crate) fn take_counts() -> std::collections::BTreeMap<String, usize> {
    COUNTS.with(|counts| std::mem::take(&mut *counts.borrow_mut()))
}
