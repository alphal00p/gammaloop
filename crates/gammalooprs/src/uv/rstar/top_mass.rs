//! Opt-in switches of the top-level R* root (a plan with `root_boundary`
//! `None`); children and their counterterms are never affected. All three are
//! inert when unset.
//!
//! * `GAMMALOOP_RSTAR_TERMS=all|bare|nonbare` (default `all`): `bare` evaluates
//!   only the empty-spinney ordinary term at the root (no proper UV child, no
//!   IR correction); `nonbare` evaluates every other R* term. A partial result
//!   keeps `rstar_mass_squared` symbolic (the mass terms cancel only in the
//!   sum of both).
//! * `GAMMALOOP_RSTAR_TOP_MASS_TERMS=keep|drop|drop-odd` (default `keep`): at
//!   the root, the outer Taylor operation `T_comp` (root edges scaled by `λ`,
//!   the reference mass compensated `M² -> (1-λ²) M²`) has two pieces whose
//!   `K` cancels in the sum of all R* terms: the `λ^0` order `F(0, M²)` and the
//!   `λ²` compensation `-M² ∂_{M²} F(0, M²)` (`K R̄* F(0, M²)` is local and
//!   homogeneous of degree one in `M²`). `drop` generates neither: no
//!   compensation substitution, and the `λ^0` coefficient is discarded right
//!   after the series, before any simplification, projection, conversion or
//!   integration; only the `p`-dependent orders at fixed reference mass remain.
//!   `drop-odd` also discards the odd orders (they integrate to zero by
//!   odd-rank vacuum symmetry); it acts only on a projected root (beta's
//!   background-field self-energy, a scalar channel after the root closure)
//!   and is `drop` elsewhere.
//! * `GAMMALOOP_TOP_MASS_SPLIT=1` (instrumentation, no value changes): at the
//!   root, each outer Taylor polynomial is split by origin (`l0`, `l1`,
//!   `l2_mom`, `comp_den`, `comp_ins`, and `l3` at degree three) and each kept
//!   piece is simplified, converted and integrated on its own, with one
//!   `TOP_MASS_PIECE` line per piece; the remainder-boundary path logs the
//!   same split of its projected coefficients (`TOP_MASS_BOUNDARY`).

use std::{cell::RefCell, sync::LazyLock};

use eyre::{Result, eyre};
use symbolica::atom::{Atom, AtomCore, AtomView, Symbol};

use crate::utils::GS;

use super::ScopedSlot;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum RStarTerms {
    All,
    Bare,
    Nonbare,
}

impl RStarTerms {
    pub fn as_str(self) -> &'static str {
        match self {
            RStarTerms::All => "all",
            RStarTerms::Bare => "bare",
            RStarTerms::Nonbare => "nonbare",
        }
    }

    /// Whether a root term with (`true`) or without proper UV children or IR
    /// counterterms is evaluated.
    pub(super) fn keeps(self, nonbare_term: bool) -> bool {
        match self {
            RStarTerms::All => true,
            RStarTerms::Bare => !nonbare_term,
            RStarTerms::Nonbare => nonbare_term,
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum TopMassTerms {
    Keep,
    Drop,
    DropOdd,
}

impl TopMassTerms {
    pub fn as_str(self) -> &'static str {
        match self {
            TopMassTerms::Keep => "keep",
            TopMassTerms::Drop => "drop",
            TopMassTerms::DropOdd => "drop-odd",
        }
    }

    /// Whether the outer Taylor order `order` (from the edge scaling alone, at
    /// fixed reference mass) is generated.
    pub(super) fn keeps_order(self, order: usize) -> bool {
        match self {
            TopMassTerms::Keep => true,
            TopMassTerms::Drop => order >= 1,
            TopMassTerms::DropOdd => order >= 1 && order.is_multiple_of(2),
        }
    }

    /// Whether the mass compensation is applied.
    pub(super) fn compensates(self) -> bool {
        self == TopMassTerms::Keep
    }
}

fn parse<T>(name: &str, values: &[(&str, T)], default: T) -> std::result::Result<T, String>
where
    T: Copy,
{
    match std::env::var(name) {
        Err(std::env::VarError::NotPresent) => Ok(default),
        Err(std::env::VarError::NotUnicode(value)) => {
            Err(format!("{name}={value:?} is not unicode"))
        }
        Ok(value) => {
            let value = value.trim();
            if value.is_empty() {
                return Ok(default);
            }
            values
                .iter()
                .find(|(name, _)| *name == value)
                .map(|(_, mode)| *mode)
                .ok_or_else(|| {
                    format!(
                        "{name}={value:?}: expected unset, {}",
                        values
                            .iter()
                            .map(|(name, _)| *name)
                            .collect::<Vec<_>>()
                            .join(", ")
                    )
                })
        }
    }
}

static TERMS: LazyLock<std::result::Result<RStarTerms, String>> = LazyLock::new(|| {
    parse(
        "GAMMALOOP_RSTAR_TERMS",
        &[
            ("all", RStarTerms::All),
            ("bare", RStarTerms::Bare),
            ("nonbare", RStarTerms::Nonbare),
        ],
        RStarTerms::All,
    )
});

static MASS_TERMS: LazyLock<std::result::Result<TopMassTerms, String>> = LazyLock::new(|| {
    parse(
        "GAMMALOOP_RSTAR_TOP_MASS_TERMS",
        &[
            ("keep", TopMassTerms::Keep),
            ("drop", TopMassTerms::Drop),
            ("drop-odd", TopMassTerms::DropOdd),
        ],
        TopMassTerms::Keep,
    )
});

static SPLIT: LazyLock<bool> = LazyLock::new(|| {
    matches!(
        std::env::var("GAMMALOOP_TOP_MASS_SPLIT").as_deref(),
        Ok("1" | "on")
    )
});

/// `GAMMALOOP_RSTAR_TERMS`, read once; a bad value is an error at first use.
pub fn rstar_terms() -> Result<RStarTerms> {
    TERMS.clone().map_err(|error| eyre!(error))
}

/// `GAMMALOOP_RSTAR_TOP_MASS_TERMS`, read once; a bad value is an error at
/// first use.
pub fn top_mass_terms() -> Result<TopMassTerms> {
    MASS_TERMS.clone().map_err(|error| eyre!(error))
}

/// `GAMMALOOP_TOP_MASS_SPLIT=1`.
pub fn split_enabled() -> bool {
    *SPLIT
}

/// Both switches parse (call at start-up for an early error).
pub fn validate_top_level_switches() -> Result<()> {
    rstar_terms()?;
    top_mass_terms()?;
    Ok(())
}

/// The switches as one start-up/provenance fragment.
pub fn top_level_switches_line() -> String {
    let show = |value: std::result::Result<&'static str, String>| match value {
        Ok(value) => value.to_owned(),
        Err(error) => format!("error({error})"),
    };
    format!(
        "rstar_terms={} rstar_top_mass_terms={} top_mass_split={}",
        show(TERMS.clone().map(RStarTerms::as_str)),
        show(MASS_TERMS.clone().map(TopMassTerms::as_str)),
        u8::from(split_enabled())
    )
}

/// The switches as a provenance value (`null` when both are at their default
/// and the split is off, so that default records do not change).
pub fn top_level_switches_provenance() -> serde_json::Value {
    let terms = rstar_terms().unwrap_or(RStarTerms::All);
    let mass = top_mass_terms().unwrap_or(TopMassTerms::Keep);
    if terms == RStarTerms::All && mass == TopMassTerms::Keep && !split_enabled() {
        return serde_json::Value::Null;
    }
    serde_json::json!({
        "rstar_terms": terms.as_str(),
        "top_mass_terms": mass.as_str(),
        "split": split_enabled(),
    })
}

/// The auxiliary mass squared of the R* outer Taylor operation
/// (`rstar_mass_squared`), as the integrated results carry it.
pub fn auxiliary_mass_squared() -> Symbol {
    symbolica::symbol!("rstar_mass_squared"; Scalar)
}

/// A partial R* result split by its power of the auxiliary mass squared,
/// `[M²^0, M²^1]` (each `together`); `log(M²)` is kept inside the
/// coefficients (a partial sum is not local: its logarithms of `M²` and
/// `μ²` and its `p²` dependence cancel only in the complete sum). An error
/// when a higher power or another dependence on `M²` remains. For
/// `GAMMALOOP_RSTAR_TERMS=bare|nonbare`.
pub fn split_auxiliary_mass(expression: &Atom) -> Result<[Atom; 2]> {
    use symbolica::{function, id::Replacement};
    let mass = auxiliary_mass_squared();
    let log = function!(Symbol::LOG, Atom::var(mass));
    let placeholder = Atom::var(symbolica::symbol!("top_mass_log_rstar_mass_squared"));
    let hidden =
        expression.replace_multiple(&[Replacement::new(log.to_pattern(), placeholder.clone())]);
    let restore = |atom: Atom| {
        atom.replace_multiple(&[Replacement::new(placeholder.to_pattern(), log.clone())])
            .together()
    };
    let constant = hidden.replace(mass).with(Atom::Zero).together();
    let linear = hidden
        .derivative(mass)
        .replace(mass)
        .with(Atom::Zero)
        .together();
    let rest = (&hidden - &constant - Atom::var(mass) * &linear).together();
    if !rest.is_zero() {
        return Err(eyre!(
            "the result is not linear in rstar_mass_squared: rest {rest}"
        ));
    }
    Ok([restore(constant), restore(linear)])
}

thread_local! {
    /// The mode of the root outer integration running on this thread; `None`
    /// outside a top-level root (children, quotients, tests).
    static ROOT: RefCell<Option<TopMassTerms>> = const { RefCell::new(None) };
}

/// Run `f` as an outer integration of the top-level root with mode `mode`.
pub(super) fn in_root_scope<T>(mode: Option<TopMassTerms>, f: impl FnOnce() -> T) -> T {
    match mode {
        None => f(),
        Some(mode) => {
            let scope = ScopedSlot::enter(&ROOT, mode);
            let result = f();
            scope.exit();
            result
        }
    }
}

/// The root mode of the enclosing [`in_root_scope`], if any.
pub(super) fn root_mode() -> Option<TopMassTerms> {
    ROOT.with(|slot| *slot.borrow())
}

/// The origin of one piece of the split outer Taylor polynomial.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) enum Origin {
    /// `[λ^n]` of the edge scaling at fixed reference mass.
    Momentum(usize),
    /// The compensation of the reference mass in `den` (the marked line).
    CompensationDen,
    /// The compensation of the reference mass outside `den` (child
    /// insertions at the reference mass).
    CompensationInsertion,
}

impl Origin {
    pub(super) fn label(self) -> String {
        match self {
            Origin::Momentum(2) => "l2_mom".into(),
            Origin::Momentum(order) => format!("l{order}"),
            Origin::CompensationDen => "comp_den".into(),
            Origin::CompensationInsertion => "comp_ins".into(),
        }
    }

    pub(super) fn kept(self, mode: TopMassTerms) -> bool {
        match self {
            Origin::Momentum(order) => mode.keeps_order(order),
            Origin::CompensationDen | Origin::CompensationInsertion => mode.compensates(),
        }
    }
}

/// `expression` with the reference mass squared replaced by `in_den` inside
/// every `den` and by `outside` everywhere else.
pub(super) fn compensate_where(
    expression: &Atom,
    reference_mass_squared: Symbol,
    in_den: &Atom,
    outside: &Atom,
) -> Atom {
    expression.replace_map(|view, _, output| match view {
        AtomView::Fun(function) if function.get_symbol() == GS.den => {
            **output = view
                .to_owned()
                .replace(reference_mass_squared)
                .with(in_den.clone());
        }
        AtomView::Var(variable) if variable.get_symbol() == reference_mass_squared => {
            **output = outside.clone();
        }
        _ => {}
    })
}

/// `GAMMALOOP_IR_QUOTIENT=rstar|massive|off|auto|check` (unset or empty:
/// `rstar`, production since r5k-nf1). The decorated IR quotient
/// counterterm `Z̃_F(G/γ; E)`: by R* on the same decorated graph as the
/// massive route (the recursive-R* children's route: massless p-integrals),
/// or by the massive R operation (mUV tadpoles, at most four loops):
/// `massive`, and `off` (the massive route without the per-quotient
/// `IR_QUOTIENT` log line, the behaviour before the production default);
/// `auto` = R* only above four loops; `check` computes both, compares them
/// and uses the massive value.
///
/// The R* route writes a quotient's external momentum through one leg
/// (`Q(1)·Q(1)` where the massive route has `Q(0)·Q(1)`); on two-point
/// kinematics (`Q(0) = Q(1)`, imposed by the root closure and the beta
/// projection) every quotient checked equals the massive one
/// (results/ir-quotient-rstar). The uv::rstar test
/// `four_loop_nested_ghost_hard_taylor_is_hard_line_independent` therefore
/// compares on those kinematics: before them, its marked line 2 carries
/// `log(rstar_mass_squared)` terms, all proportional to `(Q(0) - Q(1))·Q(1)`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum IrQuotient {
    /// `off`: the massive route, unlogged (the unset behaviour before the
    /// production default).
    Off,
    Massive,
    /// `rstar`, unset or empty.
    Rstar,
    Auto,
    Check,
}

static IR_QUOTIENT: LazyLock<std::result::Result<IrQuotient, String>> = LazyLock::new(|| {
    parse(
        "GAMMALOOP_IR_QUOTIENT",
        &[
            ("massive", IrQuotient::Massive),
            ("off", IrQuotient::Off),
            ("rstar", IrQuotient::Rstar),
            ("auto", IrQuotient::Auto),
            ("check", IrQuotient::Check),
        ],
        IrQuotient::Rstar,
    )
});

/// `GAMMALOOP_IR_QUOTIENT`, read once; a bad value is an error at first use.
pub fn ir_quotient() -> Result<IrQuotient> {
    IR_QUOTIENT.clone().map_err(|error| eyre!(error))
}

/// `GAMMALOOP_IR_OUTER=generic|analytic|check` (unset: generic, no log).
/// `analytic`: an IR correction term (of any R* root) integrates each
/// distinct closed hard coefficient once, as a fixed-mass vacuum integral
/// `h(ε) (M²)^{a-εL}`, multiplies it by the hard-loop-free quotient
/// counterterm times the root closure, and applies the compensated outer
/// Taylor analytically: the `p`-orders `b ≤ D` of that factor, times
/// `Σ_{k ≤ (D-b)/2} (-1)^k binom(a-εL, k)` (results/ir-hard-design §2.4,
/// §2.5). Declines fall back to the generic path per term. `check`: both,
/// with the pole parts compared (an error when they differ).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum IrOuter {
    Unset,
    Generic,
    Analytic,
    Check,
}

static IR_OUTER: LazyLock<std::result::Result<IrOuter, String>> = LazyLock::new(|| {
    parse(
        "GAMMALOOP_IR_OUTER",
        &[
            ("generic", IrOuter::Generic),
            ("analytic", IrOuter::Analytic),
            ("check", IrOuter::Check),
        ],
        IrOuter::Unset,
    )
});

/// `GAMMALOOP_IR_OUTER`, read once; a bad value is an error at first use.
pub fn ir_outer() -> Result<IrOuter> {
    IR_OUTER.clone().map_err(|error| eyre!(error))
}

/// The power `a` of `mass` in a single-mass vacuum integral `h (M²)^{a-εL}`
/// as integrated (its `ε` expansion carries `log(M²)`): the unique power of
/// `mass` over the summands once `log(mass)` is hidden; `None` when `mass`
/// occurs in any other way or with two powers.
pub(super) fn mass_power(integral: &Atom, mass: Symbol) -> Option<i64> {
    use symbolica::{function, id::Replacement};
    let log = function!(Symbol::LOG, Atom::var(mass));
    let placeholder = Atom::var(symbolica::symbol!("top_mass_log_mass_power"));
    let hidden = integral
        .replace_multiple(&[Replacement::new(log.to_pattern(), placeholder)])
        .expand();
    let mut powers = std::collections::BTreeSet::new();
    let summands = match hidden.as_view() {
        AtomView::Add(sum) => sum.iter().collect::<Vec<_>>(),
        view if view.is_zero() => Vec::new(),
        view => vec![view],
    };
    for summand in summands {
        let factors = match summand {
            AtomView::Mul(product) => product.iter().collect::<Vec<_>>(),
            view => vec![view],
        };
        let mut power = 0_i64;
        for factor in factors {
            match factor {
                AtomView::Var(variable) if variable.get_symbol() == mass => power += 1,
                AtomView::Pow(pow) => {
                    let (base, exponent) = pow.get_base_exp();
                    match base {
                        AtomView::Var(variable) if variable.get_symbol() == mass => {
                            power += i64::try_from(exponent).ok()?;
                        }
                        _ if factor.contains_symbol(mass) => return None,
                        _ => {}
                    }
                }
                _ if factor.contains_symbol(mass) => return None,
                _ => {}
            }
        }
        powers.insert(power);
    }
    match powers.len() {
        1 => powers.into_iter().next(),
        _ => None,
    }
}

/// The analytic compensated outer Taylor factor of a `p`-order `b` of an
/// integral `∝ (M²)^{a-εL}` at total degree `degree` under `mode`:
/// `Σ_k (-1)^k binom(a-εL, k)` over the `k ≤ (degree-b)/2` the mode keeps
/// (`keep`: all; `drop`/`drop-odd`: `k = 0` when the order `b` is kept).
pub(super) fn compensation_factor(
    mode: TopMassTerms,
    a: i64,
    loops: usize,
    order: usize,
    degree: usize,
    epsilon: Symbol,
) -> Atom {
    let x = Atom::num(a) - Atom::num(loops as i64) * Atom::var(epsilon);
    let mut sum = Atom::Zero;
    let mut binomial = Atom::one();
    let mut k = 0;
    while order + 2 * k <= degree {
        if mode.keeps_order(order + 2 * k) && (k == 0 || mode.compensates()) {
            let sign = if k % 2 == 0 { 1 } else { -1 };
            sum += Atom::num(sign) * &binomial;
        }
        // binom(x, k+1) = binom(x, k) (x-k)/(k+1)
        binomial = binomial * (&x - Atom::num(k as i64)) / Atom::num(k as i64 + 1);
        k += 1;
    }
    sum
}

#[cfg(test)]
mod tests {
    use super::*;
    use symbolica::{function, parse, symbol};

    #[test]
    fn compensation_splits_into_den_and_insertions() {
        let m = symbol!("top_mass_test_m");
        let (a, b) = (parse!("top_mass_test_a"), parse!("top_mass_test_b"));
        let den = function!(GS.den, Atom::num(1), parse!("q"), parse!("q"), Atom::var(m));
        let expression = (Atom::var(m) * &den + Atom::var(m).pow(2)) * den.pow(2);
        let den_a = function!(GS.den, Atom::num(1), parse!("q"), parse!("q"), &a);
        let den_b = function!(GS.den, Atom::num(1), parse!("q"), parse!("q"), &b);
        assert_eq!(
            compensate_where(&expression, m, &a, &b),
            (&b * &den_a + b.pow(2)) * den_a.pow(2)
        );
        assert_eq!(
            compensate_where(&expression, m, &b, &a),
            (&a * &den_b + a.pow(2)) * den_b.pow(2)
        );
    }

    #[test]
    fn modes_keep_the_documented_orders() {
        assert!((0..4).all(|order| TopMassTerms::Keep.keeps_order(order)));
        assert_eq!(
            (0..4)
                .map(|order| TopMassTerms::Drop.keeps_order(order))
                .collect::<Vec<_>>(),
            [false, true, true, true]
        );
        assert_eq!(
            (0..4)
                .map(|order| TopMassTerms::DropOdd.keeps_order(order))
                .collect::<Vec<_>>(),
            [false, false, true, false]
        );
        assert!(RStarTerms::Nonbare.keeps(true) && !RStarTerms::Nonbare.keeps(false));
        assert!(!RStarTerms::Bare.keeps(true) && RStarTerms::Bare.keeps(false));
    }
}

#[cfg(test)]
mod ir_outer_tests {
    use super::*;
    use symbolica::{parse, symbol};

    #[test]
    fn compensation_factor_is_the_truncated_binomial_series() {
        let eps = symbol!("top_mass_test_eps");
        // keep, D = 2: b = 0 -> 1 - (a - εL); b = 1, 2 -> 1.
        let x = Atom::num(1) - Atom::num(3) * Atom::var(eps);
        assert_eq!(
            compensation_factor(TopMassTerms::Keep, 1, 3, 0, 2, eps),
            Atom::one() - &x
        );
        assert_eq!(
            compensation_factor(TopMassTerms::Keep, 1, 3, 1, 2, eps),
            Atom::one()
        );
        assert_eq!(
            compensation_factor(TopMassTerms::Keep, 1, 3, 2, 2, eps),
            Atom::one()
        );
        // drop: no b = 0 term and no compensation; drop-odd: no odd b.
        assert!(compensation_factor(TopMassTerms::Drop, 1, 3, 0, 2, eps).is_zero());
        assert_eq!(
            compensation_factor(TopMassTerms::Drop, 1, 3, 1, 2, eps),
            Atom::one()
        );
        assert!(compensation_factor(TopMassTerms::DropOdd, 1, 3, 1, 2, eps).is_zero());
        // D = 4, b = 0: 1 - x + x(x-1)/2.
        let expected = (Atom::one() - &x + &x * (&x - Atom::one()) / Atom::num(2)).expand();
        assert_eq!(
            compensation_factor(TopMassTerms::Keep, 1, 3, 0, 4, eps).expand(),
            expected
        );
    }

    #[test]
    fn mass_power_reads_the_unique_power_beside_logarithms() {
        let m = symbol!("top_mass_test_m2");
        assert_eq!(
            mass_power(&parse!("top_mass_test_m2^2*(3+log(top_mass_test_m2)*e)"), m),
            Some(2)
        );
        assert_eq!(mass_power(&parse!("5*e"), m), Some(0));
        assert_eq!(mass_power(&parse!("top_mass_test_m2+e"), m), None);
        assert_eq!(mass_power(&parse!("sqrt(top_mass_test_m2)"), m), None);
    }
}
