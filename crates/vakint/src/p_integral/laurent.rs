//! Exact Laurent expansion of finalized p-integrals.
//!
//! Every term of a p-integral is `rational(ε) × Π Γ(n + cε)^e × exp(...) ×
//! PMaster(...)`. Symbolica's generic series multiplies the exponential
//! series of each Γ function with unexpanded transcendental coefficients,
//! which explodes at four loops (minutes for one ε^0 coefficient). Here every
//! term is split into a rational part, expanded by Symbolica, and a single
//! exponential whose exponent collects `log Γ(1+cε) = -γcε + Σ_k (-1)^k ζ(k)(cε)^k/k`
//! and the logarithms of the scale factors. The exponential is expanded by
//! the recurrence `f_n = Σ_k k e_k f_(n-k) / n` with expanded coefficients,
//! each factor exactly to the depth its partners' valuations require.

use symbolica::atom::{Atom, AtomCore, AtomView, Symbol};
use symbolica::coefficient::CoefficientView;
use symbolica::domains::rational::Rational;
use symbolica::function;
use symbolica::poly::series::SeriesDepth;
use symbolica::transcendental::{euler_gamma, gamma, zeta};

use super::{masters, p_gamma_symbol, p_master_symbol, pg_symbol};
use crate::{VakintError, get_integer_from_atom};

#[path = "rational_series.rs"]
mod rational_series;
#[cfg(test)]
#[path = "rational_series_tests.rs"]
mod rational_series_tests;

fn failure(detail: impl std::fmt::Display) -> VakintError {
    VakintError::SymbolicaError(format!("p-integral expansion: {detail}"))
}

/// `Σ_(k=valuation)^(order) coefficients[k-valuation] ε^k`, exact through
/// `order = valuation + coefficients.len() - 1`.
#[derive(Clone, Debug)]
pub(crate) struct Laurent {
    valuation: i64,
    coefficients: Vec<Atom>,
}

impl Laurent {
    fn order(&self) -> i64 {
        self.valuation + self.coefficients.len() as i64 - 1
    }

    /// The product, exact through `order` or as far as the factors allow.
    fn mul(&self, other: &Self, order: i64) -> Self {
        let valuation = self.valuation + other.valuation;
        let order = order
            .min(self.order() + other.valuation)
            .min(other.order() + self.valuation);
        let length = (order - valuation + 1).max(0) as usize;
        let coefficients = (0..length)
            .map(|n| {
                let mut total = Atom::Zero;
                for i in 0..=n {
                    if let (Some(a), Some(b)) =
                        (self.coefficients.get(i), other.coefficients.get(n - i))
                        && !a.is_zero()
                        && !b.is_zero()
                    {
                        total += a * b;
                    }
                }
                total.expand()
            })
            .collect();
        Laurent {
            valuation,
            coefficients,
        }
    }

    fn to_atom(&self, epsilon: Symbol) -> Atom {
        let e = Atom::var(epsilon);
        self.coefficients
            .iter()
            .enumerate()
            .fold(Atom::Zero, |sum, (index, coefficient)| {
                sum + coefficient * e.pow(Atom::num(self.valuation + index as i64))
            })
    }

    /// An expanded polynomial in `ε` and `1/ε` (a tabulated master).
    fn from_polynomial(expression: &Atom, epsilon: Symbol) -> Result<Self, VakintError> {
        let mut terms = std::collections::BTreeMap::<i64, Atom>::new();
        let expanded = expression.expand();
        let parts = match expanded.as_view() {
            AtomView::Add(sum) => sum.iter().map(|term| term.to_owned()).collect(),
            _ => vec![expanded.clone()],
        };
        for term in parts {
            let power = epsilon_power(term.as_view(), epsilon)
                .ok_or_else(|| failure(format!("not a Laurent polynomial term: {term}")))?;
            let coefficient = term / Atom::var(epsilon).pow(Atom::num(power));
            *terms.entry(power).or_default() += coefficient;
        }
        let Some((&valuation, _)) = terms.first_key_value() else {
            return Ok(Laurent {
                valuation: 0,
                coefficients: Vec::new(),
            });
        };
        let order = *terms.last_key_value().unwrap().0;
        Ok(Laurent {
            valuation,
            coefficients: (valuation..=order)
                .map(|power| terms.remove(&power).unwrap_or_default().expand())
                .collect(),
        })
    }

    /// Symbolica's Laurent series of `expression` through `order`.
    fn series(expression: &Atom, epsilon: Symbol, order: i64) -> Result<Self, VakintError> {
        let series = expression
            .series(epsilon, Atom::Zero.as_view(), Rational::from(order))
            .map_err(failure)?;
        let mut terms = std::collections::BTreeMap::<i64, Atom>::new();
        for (exponent, coefficient) in series.terms() {
            if !exponent.is_integer() {
                return Err(failure(format!("fractional ε power in {expression}")));
            }
            let power = exponent
                .numerator()
                .to_i64()
                .ok_or_else(|| failure("ε power"))?;
            if power <= order {
                terms.insert(power, coefficient.expand());
            }
        }
        let valuation = terms
            .first_key_value()
            .map_or(order + 1, |(power, _)| *power);
        Ok(Laurent {
            valuation,
            coefficients: (valuation..=order)
                .map(|power| terms.remove(&power).unwrap_or_default())
                .collect(),
        })
    }

    /// `exp(Σ_(k>=1) exponent[k-1] ε^k)` through `ε^order`.
    fn exp(exponent: &[Atom], order: i64) -> Self {
        let mut coefficients: Vec<Atom> = vec![Atom::num(1)];
        for n in 1..=order.max(0) as usize {
            let mut total = Atom::Zero;
            for k in 1..=n.min(exponent.len()) {
                if !exponent[k - 1].is_zero() && !coefficients[n - k].is_zero() {
                    total += Atom::num(k as i64) * &exponent[k - 1] * &coefficients[n - k];
                }
            }
            coefficients.push((total / Atom::num(n as i64)).expand());
        }
        coefficients.truncate(order.max(-1) as usize + 1);
        Laurent {
            valuation: 0,
            coefficients,
        }
    }
}

/// A tabulated master `PMaster(name)`: its known Laurent coefficients and,
/// when the table stops, its first unknown order `k` (the table's tail
/// `ε^k·Oep(k, PMaster(name))`).
///
/// The tail stands for every order from `k` on, not for one coefficient:
/// [`TabulatedMaster::through`] writes each unknown order `j >= k` that a
/// product needs as its own marker `ε^j·Oep(j, PMaster(name))`. Distinct
/// unknown coefficients are distinct symbols, so a marker cancels between
/// terms only where the unknown coefficient itself cancels, and every other
/// one survives truncation for the guard
/// ([`crate::truncation_guard::ensure_known_orders`]).
#[derive(Clone, Debug)]
struct TabulatedMaster {
    /// The known coefficients (orders below `k`).
    known: Laurent,
    /// `(k, PMaster(name))` when the table has an unknown tail.
    tail: Option<(i64, Atom)>,
    /// The lowest order of the table (known or the tail): a lower bound of
    /// the master's valuation.
    valuation: i64,
}

impl TabulatedMaster {
    /// Split a table entry (an expanded Laurent polynomial in `ε`, at most
    /// one term `ε^k·Oep(k, master)`) into known coefficients and tail.
    fn new(expression: &Atom, epsilon: Symbol) -> Result<Self, VakintError> {
        let expanded = expression.expand();
        let parts: Vec<Atom> = match expanded.as_view() {
            AtomView::Add(sum) => sum.iter().map(|term| term.to_owned()).collect(),
            _ => vec![expanded.clone()],
        };
        let marker = crate::truncation_guard::unknown_tail_symbol();
        let mut known = Atom::Zero;
        let mut tail = None;
        for term in parts {
            if !term.contains_symbol(marker) {
                known += term;
                continue;
            }
            let found = Self::tail_of(&term, epsilon)
                .ok_or_else(|| failure(format!("malformed unknown master tail {term}")))?;
            if tail.replace(found).is_some() {
                return Err(failure(format!(
                    "more than one unknown tail in the master table entry {expression}"
                )));
            }
        }
        let known = Laurent::from_polynomial(&known, epsilon)?;
        let known_valuation = (!known.coefficients.is_empty()).then_some(known.valuation);
        if let (Some((k, _)), false) = (&tail, known.coefficients.is_empty())
            && known.order() >= *k
        {
            return Err(failure(format!(
                "the master table entry {expression} has known coefficients at or beyond its unknown tail"
            )));
        }
        let valuation = match (known_valuation, &tail) {
            (Some(v), _) => v,
            (None, Some((k, _))) => *k,
            (None, None) => 0,
        };
        Ok(TabulatedMaster {
            known,
            tail,
            valuation,
        })
    }

    /// `(k, master)` of a term exactly `ε^k·Oep(k, master)`.
    fn tail_of(term: &Atom, epsilon: Symbol) -> Option<(i64, Atom)> {
        let mut marker = None;
        let mut markers = 0;
        term.visitor(&mut |part| {
            if let AtomView::Fun(function) = part
                && function.get_symbol() == crate::truncation_guard::unknown_tail_symbol()
            {
                markers += 1;
                marker = Some(part.to_owned());
                return false;
            }
            true
        });
        let marker = marker.filter(|_| markers == 1)?;
        let AtomView::Fun(function) = marker.as_view() else {
            return None;
        };
        let arguments: Vec<_> = function.iter().collect();
        let [k, master] = arguments.as_slice() else {
            return None;
        };
        let k = get_integer_from_atom(*k)?;
        let rest = (term / &marker - Atom::var(epsilon).pow(Atom::num(k))).expand();
        (rest.is_zero() && !master.contains_symbol(epsilon)).then(|| (k, master.to_owned()))
    }

    /// The master through `ε^through`: known coefficients, zero between the
    /// last known order and the tail, and `Oep(j, master)` for every unknown
    /// order `j` up to `through`. Empty (exactly zero through `through`) when
    /// `through` is below the valuation.
    fn through(&self, through: i64) -> Laurent {
        let length = (through - self.valuation + 1).max(0) as usize;
        let coefficients = (0..length)
            .map(|index| {
                let power = self.valuation + index as i64;
                match &self.tail {
                    Some((k, master)) if power >= *k => {
                        crate::truncation_guard::unknown_tail_marker(power, master)
                    }
                    _ => self
                        .known
                        .coefficients
                        .get((power - self.known.valuation) as usize)
                        .filter(|_| power >= self.known.valuation)
                        .cloned()
                        .unwrap_or_default(),
                }
            })
            .collect();
        Laurent {
            valuation: self.valuation,
            coefficients,
        }
    }
}

/// The power of `ε` in a monomial, if it is an integer.
fn epsilon_power(term: AtomView, epsilon: Symbol) -> Option<i64> {
    match term {
        AtomView::Var(variable) if variable.get_symbol() == epsilon => Some(1),
        AtomView::Pow(power) => {
            let (base, exponent) = power.get_base_exp();
            if base == Atom::var(epsilon).as_view() {
                get_integer_from_atom(exponent)
            } else if term.contains_symbol(epsilon) {
                None
            } else {
                Some(0)
            }
        }
        AtomView::Mul(product) => product
            .iter()
            .try_fold(0, |sum, factor| Some(sum + epsilon_power(factor, epsilon)?)),
        _ if !term.contains_symbol(epsilon) => Some(0),
        _ => None,
    }
}

/// The coefficients `x_0, x_1, ...` of a polynomial in `ε`, or `None`.
fn polynomial_coefficients(expression: AtomView, epsilon: Symbol) -> Option<Vec<Atom>> {
    let laurent = Laurent::from_polynomial(&expression.to_owned(), epsilon).ok()?;
    if laurent.coefficients.is_empty() {
        return Some(vec![Atom::Zero]);
    }
    if laurent.valuation < 0 {
        return None;
    }
    let mut coefficients = vec![Atom::Zero; laurent.valuation as usize];
    coefficients.extend(laurent.coefficients);
    Some(coefficients)
}

/// A source of exponent coefficients.
enum ExponentPart {
    /// `Σ_(k>=1) coefficients[k] ε^k` (index 0 unused).
    Polynomial(Vec<Atom>),
    /// `multiplicity × log Γ(1 + c·ε)`.
    LogGamma { c: Atom, multiplicity: Atom },
}

/// A product split into `rational × exp(exponent)`.
struct Exponential {
    rational: Atom,
    exponent: Vec<ExponentPart>,
}

impl Exponential {
    /// The exponent coefficients of `ε^1..=ε^depth`.
    fn exponent(&self, depth: i64) -> Vec<Atom> {
        let mut coefficients = vec![Atom::Zero; depth.max(0) as usize];
        for part in &self.exponent {
            match part {
                ExponentPart::Polynomial(polynomial) => {
                    for (k, value) in polynomial.iter().enumerate().skip(1) {
                        if let Some(slot) = coefficients.get_mut(k - 1) {
                            *slot += value;
                        }
                    }
                }
                ExponentPart::LogGamma { c, multiplicity } => {
                    for k in 1..=depth {
                        let value = if k == 1 {
                            -Atom::var(euler_gamma()) * c
                        } else {
                            let sign = if k % 2 == 0 { 1 } else { -1 };
                            Atom::num(sign) * zeta().call(Atom::num(k)) * c.pow(Atom::num(k))
                                / Atom::num(k)
                        };
                        coefficients[k as usize - 1] += multiplicity * value;
                    }
                }
            }
        }
        coefficients.into_iter().map(|c| c.expand()).collect()
    }
}

/// `n + c·ε` of a linear argument, with integer `n` and numeric `c`.
fn linear_argument(argument: AtomView, epsilon: Symbol) -> Option<(i64, Atom)> {
    let coefficients = polynomial_coefficients(argument, epsilon)?;
    if coefficients.len() > 2 {
        return None;
    }
    let n = get_integer_from_atom(coefficients[0].as_view())?;
    let c = coefficients.get(1).cloned().unwrap_or_default();
    matches!(c.as_view(), AtomView::Num(_)).then_some((n, c))
}

/// The logarithm of a positive ε-free base, exact for `𝑒^y` and powers.
fn logarithm(base: AtomView) -> Atom {
    match base {
        AtomView::Var(variable) if variable.get_symbol() == Symbol::E => Atom::num(1),
        AtomView::Pow(power) => {
            let (base, exponent) = power.get_base_exp();
            exponent.to_owned() * logarithm(base)
        }
        AtomView::Fun(function) if function.get_symbol() == Symbol::EXP => {
            function.iter().next().unwrap().to_owned()
        }
        _ => function!(Symbol::LOG, base.to_owned()),
    }
}

/// Split a product into rational factors and one exponential.
fn split(expression: &Atom, epsilon: Symbol) -> Exponential {
    let e = Atom::var(epsilon);
    let mut rational = Atom::num(1);
    let mut exponent = Vec::new();
    let factors = match expression.as_view() {
        AtomView::Mul(product) => product.iter().collect::<Vec<_>>(),
        _ => vec![expression.as_view()],
    };
    for factor in factors {
        let (base, power) = match factor {
            AtomView::Pow(pow) => {
                let (base, exponent) = pow.get_base_exp();
                (base, exponent.to_owned())
            }
            _ => (factor, Atom::num(1)),
        };
        if !factor.contains_symbol(epsilon) {
            rational *= factor.to_owned();
            continue;
        }
        // Γ(n + cε)^power with an integer power.
        if let AtomView::Fun(function) = base
            && function.get_symbol() == gamma()
            && let Some(p) = get_integer_from_atom(power.as_view())
            && let Some((n, c)) = linear_argument(function.iter().next().unwrap(), epsilon)
        {
            let z = &c * &e;
            let mut shift = Atom::num(1);
            for m in 1..n {
                shift *= Atom::num(m) + &z;
            }
            for m in n..=0 {
                shift /= Atom::num(m) + &z;
            }
            rational *= shift.pow(Atom::num(p));
            if !c.is_zero() {
                exponent.push(ExponentPart::LogGamma {
                    c,
                    multiplicity: Atom::num(p),
                });
            }
            continue;
        }
        // base^power with an ε-free base and a polynomial exponent in ε
        // (Symbolica writes exp(y) as 𝑒^y).
        if base.contains_symbol(epsilon) {
            rational *= factor.to_owned();
            continue;
        }
        let (log_base, total_power) = (logarithm(base), power.clone());
        let Some(coefficients) =
            polynomial_coefficients((&log_base * &total_power).as_view(), epsilon)
        else {
            rational *= factor.to_owned();
            continue;
        };
        let constant_power = polynomial_coefficients(power.as_view(), epsilon)
            .map(|c| c[0].clone())
            .unwrap_or_default();
        rational *= base.to_owned().pow(constant_power);
        exponent.push(ExponentPart::Polynomial(coefficients));
    }
    Exponential { rational, exponent }
}

/// Capture a generic input before expansion, so a stalled call remains replayable.
fn trace_generic<T>(
    expression: &Atom,
    epsilon: Symbol,
    monomial: &Atom,
    order: i64,
    phase: &'static str,
    exact_plan_available: bool,
    evaluate: impl FnOnce() -> Result<T, VakintError>,
) -> Result<T, VakintError> {
    // diagnostics: the declined input to a file named by the caller
    if let Some(path) = std::env::var_os("VAKINT_SERIES_DECLINE_LOG") {
        use std::io::Write;
        if let Ok(mut file) = std::fs::OpenOptions::new()
            .create(true)
            .append(true)
            .open(path)
        {
            let _ = writeln!(
                file,
                "DECLINE phase={phase} order={order} exact_plan={exact_plan_available}\nmonomial={}\nexpression={}",
                monomial.to_canonical_string(),
                expression.to_canonical_string()
            );
        }
    }
    if !tracing::enabled!(target: "vakint::p_integral::series_input", tracing::Level::INFO) {
        return evaluate();
    }
    static NEXT_CALL: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);
    let call_id = NEXT_CALL.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
    let thread_id = std::thread::current().id();
    tracing::info!(target: "vakint::p_integral::series_input",
        stage = "generic_series_start",
        call_id,
        thread_id = ?thread_id,
        phase,
        target_order = order,
        exact_plan_available,
        input_bytes = expression.as_view().get_byte_size(),
        file.expr = %expression.to_canonical_string(),
        file.epsilon = %Atom::var(epsilon).to_canonical_string(),
        file.master_monomial = %monomial.to_canonical_string(),
        "Starting generic epsilon expansion"
    );
    let started = std::time::Instant::now();
    let result = evaluate();
    let seconds = started.elapsed().as_secs_f64();
    tracing::info!(target: "vakint::p_integral::series_input",
        stage = "generic_series_end",
        call_id,
        thread_id = ?thread_id,
        phase,
        target_order = order,
        success = result.is_ok(),
        seconds,
        "Finished generic epsilon expansion"
    );
    result
}

/// The valuation of `expression` in `ε`, `None` for an exact zero.
fn valuation(expression: &Atom, epsilon: Symbol) -> Result<Option<i64>, VakintError> {
    if expression.is_zero() {
        return Ok(None);
    }
    let series = expression
        .series(epsilon, Atom::Zero.as_view(), SeriesDepth::relative(1))
        .map_err(failure)?;
    if series.is_zero() {
        return Ok(None);
    }
    let exponent = series.get_trailing_exponent();
    if !exponent.is_integer() {
        return Err(failure(format!("fractional ε power in {expression}")));
    }
    Ok(Some(
        exponent
            .numerator()
            .to_i64()
            .ok_or_else(|| failure("ε power"))?,
    ))
}

/// The first power in `expression` of a function in `symbols` whose exponent
/// is not an `i16` integer. The collection in `i16` exponents (the range of
/// RustRed's `GammaTag`) panics or leaves such a power in a coefficient.
pub(crate) fn unsupported_power(expression: &Atom, symbols: &[Symbol]) -> Option<Atom> {
    let mut unsupported = None;
    expression.visitor(&mut |part| {
        if let AtomView::Pow(power) = part {
            let (base, exponent) = power.get_base_exp();
            let i16_integer = get_integer_from_atom(exponent)
                .is_some_and(|exponent| i16::try_from(exponent).is_ok());
            if !i16_integer
                && matches!(base, AtomView::Fun(function) if symbols.contains(&function.get_symbol()))
            {
                unsupported = Some(part.to_owned());
            }
        }
        unsupported.is_none()
    });
    unsupported
}

/// The functions `PG`, `PMaster` and `PGamma` of the master basis.
fn master_symbols() -> [Symbol; 3] {
    [pg_symbol(), p_master_symbol(), p_gamma_symbol()]
}

/// The master functions of `reduced`, sorted and deduplicated, and
/// `reduced` expanded; an unsupported power of a master function fails.
fn masters_and_expansion(
    reduced: &Atom,
    symbols: &[Symbol; 3],
) -> Result<(Vec<Atom>, Atom), VakintError> {
    let mut masters = Vec::new();
    reduced.visitor(&mut |part| {
        if let AtomView::Fun(function) = part
            && symbols.contains(&function.get_symbol())
        {
            masters.push(part.to_owned());
            return false;
        }
        true
    });
    masters.sort();
    masters.dedup();
    let expanded = reduced.expand();
    if let Some(power) = unsupported_power(&expanded, symbols) {
        return Err(failure(format!("{power} is not an i16 integer power")));
    }
    Ok((masters, expanded))
}

/// How the rational factors of one expansion were expanded.
#[derive(Default)]
struct SeriesCounts {
    exact: usize,
    generic: usize,
}

/// Master monomial terms (process total) in which a tabulated master's
/// unknown tail lies inside the orders the term needs: before lane
/// laurent-depth these were the terms that the old expansion clipped below the
/// requested order (or that the guard refused). A diagnostic for exposure
/// scans: a result with none of them is byte-identical to the old expansion.
static TAIL_IN_RANGE_TERMS: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);

/// See [`TAIL_IN_RANGE_TERMS`].
#[cfg_attr(not(test), allow(dead_code))]
pub(crate) fn tail_in_range_terms() -> u64 {
    TAIL_IN_RANGE_TERMS.load(std::sync::atomic::Ordering::Relaxed)
}

#[cfg(test)]
thread_local! {
    /// Tests only: every rational factor of [`monomial_term`] takes the
    /// generic series, the reference of the exact route.
    static FORCE_GENERIC_SERIES: std::cell::Cell<bool> = const { std::cell::Cell::new(false) };
}

/// The exact plan of a rational factor, `None` for the generic series.
fn prepare_rational(rational: &Atom, epsilon: Symbol) -> Option<rational_series::RationalSeries> {
    #[cfg(test)]
    if FORCE_GENERIC_SERIES.with(std::cell::Cell::get) {
        return None;
    }
    rational_series::RationalSeries::new(rational, epsilon)
}

/// The Laurent series through `ε^order` of `prefactor × coefficient ×
/// monomial` for one master monomial of [`expand`], `None` when it is an
/// exact zero or starts above `order`.
fn monomial_term(
    prefactor: &Atom,
    monomial: &Atom,
    coefficient: &Atom,
    symbols: &[Symbol; 3],
    epsilon: Symbol,
    order: i64,
    counts: &mut SeriesCounts,
    master: &dyn Fn(&str, Symbol) -> Result<Atom, VakintError>,
) -> Result<Option<Laurent>, VakintError> {
    if symbols.iter().any(|s| coefficient.contains_symbol(*s)) {
        return Err(failure(format!(
            "{coefficient} is not free of PG, PMaster and PGamma"
        )));
    }
    // The Γ form of the PG factors joins the prefactor; PMaster factors
    // are tabulated Laurent polynomials.
    let mut product = prefactor * coefficient;
    // Each occurrence of a tabulated master, expanded only once every
    // factor's valuation is known (see `TabulatedMaster`).
    let mut tabulated: Vec<TabulatedMaster> = Vec::new();
    let factors = match monomial.as_view() {
        AtomView::Mul(product) => product.iter().map(|f| f.to_owned()).collect(),
        _ => vec![monomial.clone()],
    };
    for factor in factors {
        let (base, power) = match factor.as_view() {
            AtomView::Pow(pow) => {
                let (base, exponent) = pow.get_base_exp();
                (base.to_owned(), get_integer_from_atom(exponent))
            }
            _ => (factor.clone(), Some(1)),
        };
        if let AtomView::Fun(function) = base.as_view()
            && function.get_symbol() == symbols[2]
        {
            // A Gamma tag: any nonzero integer power of Γ(1+mε), m ≠ 0.
            let m = match function.iter().collect::<Vec<_>>().as_slice() {
                [m] => get_integer_from_atom(*m).filter(|m| *m != 0),
                _ => None,
            };
            let (Some(m), Some(power)) = (m, power) else {
                return Err(failure(format!("invalid Gamma tag {factor}")));
            };
            product *=
                super::finalize::gamma_one_plus(m, &Atom::var(epsilon)).pow(Atom::num(power));
            continue;
        }
        let power = power
            .filter(|p| *p > 0)
            .ok_or_else(|| failure(format!("invalid master power {factor}")))?;
        let AtomView::Fun(function) = base.as_view() else {
            if base == Atom::num(1) {
                continue;
            }
            return Err(failure(format!("unexpected master factor {factor}")));
        };
        if function.get_symbol() == symbols[0] {
            let arguments = function
                .iter()
                .map(get_integer_from_atom)
                .collect::<Option<Vec<_>>>();
            let Some([x1, x2]) = arguments.as_deref() else {
                return Err(failure(format!("invalid {base}")));
            };
            product *= super::finalize::pg_value(*x1, *x2, epsilon).pow(Atom::num(power));
        } else {
            let name = masters::p_master_name(base.as_view())
                .ok_or_else(|| failure(format!("invalid {base}")))?;
            let value = TabulatedMaster::new(&master(&name, epsilon)?, epsilon)?;
            for _ in 0..power {
                tabulated.push(value.clone());
            }
        }
    }
    let split = split(
        &crate::Vakint::normalize_gamma_functions(product.as_view()),
        epsilon,
    );
    let prepared_rational = {
        let _stage = crate::stage_timing::StageSpan::new("vakint.p_integral.rational_prepare");
        prepare_rational(&split.rational, epsilon)
    };
    let rational_valuation = {
        let _stage = crate::stage_timing::StageSpan::new("vakint.p_integral.rational_valuation");
        match &prepared_rational {
            Some(prepared) => prepared.valuation(),
            None => trace_generic(
                &split.rational,
                epsilon,
                monomial,
                order - tabulated.iter().map(|t| t.valuation).sum::<i64>(),
                "relative_valuation",
                false,
                || valuation(&split.rational, epsilon),
            )?,
        }
    };
    let Some(rational_valuation) = rational_valuation else {
        return Ok(None);
    };
    let tabulated_valuation = tabulated.iter().map(|t| t.valuation).sum::<i64>();
    let depth = order - rational_valuation - tabulated_valuation;
    if depth < 0 {
        return Ok(None);
    }
    let rational = {
        let _stage = crate::stage_timing::StageSpan::new("vakint.p_integral.rational_series");
        match prepared_rational
            .as_ref()
            .and_then(|prepared| prepared.series(order - tabulated_valuation))
        {
            Some(series) => {
                counts.exact += 1;
                series
            }
            None => {
                counts.generic += 1;
                trace_generic(
                    &split.rational,
                    epsilon,
                    monomial,
                    order - tabulated_valuation,
                    "absolute_series",
                    prepared_rational.is_some(),
                    || Laurent::series(&split.rational, epsilon, order - tabulated_valuation),
                )?
            }
        }
    };
    let exponential = Laurent::exp(&split.exponent(depth), depth);
    let mut term = rational.mul(&exponential, order - tabulated_valuation);
    // Every occurrence of a master is needed through its own valuation plus
    // `depth` = order minus the valuations of ALL other factors (rational
    // part, Γ and PG factors, the prefactor and the other masters): each is
    // expanded exactly that far, an unknown order as its own marker
    // `Oep(j, master)`, so that the product is exact through `ε^order`.
    if tabulated.iter().any(|master| {
        master
            .tail
            .as_ref()
            .is_some_and(|(k, _)| *k <= master.valuation + depth)
    }) {
        TAIL_IN_RANGE_TERMS.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
    }
    // The masters' product first, then the term. No cap on a partial
    // product: with factors through v_i + depth, a product of any prefix is
    // exact through (its valuations) + depth, which Laurent::mul keeps; a cap
    // at the final valuation would cut a prefix whose remaining factors have
    // negative valuations (reviewer, three factors of PMaster(no1)).
    let mut masters: Option<Laurent> = None;
    for master in &tabulated {
        let factor = master.through(master.valuation + depth);
        masters = Some(match masters {
            None => factor,
            Some(product) => product.mul(&factor, i64::MAX / 4),
        });
    }
    if let Some(masters) = masters {
        term = term.mul(&masters, order);
    }
    if term.order() < order {
        return Err(failure(format!(
            "the expansion of {monomial} ends at ε^{} below the requested ε^{order}",
            term.order()
        )));
    }
    Ok(Some(term))
}

/// The Laurent expansion through `ε^order` of `prefactor × reduced`, where
/// `reduced` is a polynomial in the masters `PG(x1,x2)` (replaced by
/// `pg_value`) and `PMaster(name)` (replaced by its tabulated expansion and
/// unknown tail) and a Laurent polynomial in the Gamma tags `PGamma(m)`, every
/// `PGamma(m)^e = Γ(1+mε)^e` with an `i16` exponent `e` expanded exactly like
/// the Γ functions of the prefactor, and `prefactor` is a product of rational
/// factors, Γ functions and exponentials. Any other dependence on these
/// symbols is an error, never a result that still contains them.
/// (Production callers use [`expand_with`] with their master source.)
#[cfg(test)]
pub(crate) fn expand(
    prefactor: &Atom,
    reduced: &Atom,
    epsilon: Symbol,
    order: i64,
) -> Result<Atom, VakintError> {
    expand_with(
        prefactor,
        reduced,
        epsilon,
        order,
        &masters::p_master_expansion,
    )
}

/// [`expand`], with `PMaster(name)` replaced by `master(name, epsilon)` (its
/// known coefficients and unknown tail): the Forcer table for [`expand`], a
/// p master catalogue for `master_catalogue::PMasterSource::expansion`.
pub(crate) fn expand_with(
    prefactor: &Atom,
    reduced: &Atom,
    epsilon: Symbol,
    order: i64,
    master: &dyn Fn(&str, Symbol) -> Result<Atom, VakintError>,
) -> Result<Atom, VakintError> {
    let symbols = master_symbols();
    let (masters, expanded) = masters_and_expansion(reduced, &symbols)?;
    let mut total = Atom::Zero;
    let mut counts = SeriesCounts::default();
    for (monomial, coefficient) in expanded.coefficient_list::<i16>(&masters) {
        if let Some(term) = monomial_term(
            prefactor,
            &monomial,
            &coefficient,
            &symbols,
            epsilon,
            order,
            &mut counts,
            master,
        )? {
            total += term.to_atom(epsilon);
        }
    }
    tracing::info!(
        uv = true,
        rstar = true,
        timing = true,
        exact_rational_series = counts.exact,
        generic_rational_series = counts.generic,
        target_order = order,
        master_variables = masters.len(),
        "Expanded p-integral Laurent coefficients"
    );
    Ok(total.expand())
}

/// A spectator of [`expand_factored`]: an exact rational function of `ε`
/// (spectator symbols allowed) without Γ functions or exponentials, prepared
/// once for its exact valuation and series.
pub(crate) struct Spectator(rational_series::RationalSeries);

impl Spectator {
    /// `None` when `spectator` is not such a rational function, or a
    /// denominator depends on another symbol (the caller expands the product).
    pub(crate) fn new(spectator: &Atom, epsilon: Symbol) -> Option<Self> {
        if !Self::epsilon_denominators(spectator.as_view(), epsilon) {
            return None;
        }
        let split = split(
            &crate::Vakint::normalize_gamma_functions(spectator.as_view()),
            epsilon,
        );
        if !split.exponent.is_empty() {
            return None;
        }
        let _stage = crate::stage_timing::StageSpan::new("vakint.p_integral.spectator_prepare");
        rational_series::RationalSeries::new(&split.rational, epsilon).map(Self)
    }

    /// A conservative guard for atom identity, not merely algebraic equality.
    /// RationalSeries keeps the ε-free denominator content as one atom:
    /// expanding a product can produce `(3*nf-6)^-1`, while separately
    /// expanding its factors produces `1/3*(nf-2)^-1`. Admit only rational
    /// denominators in ε, whose ε-free content is numeric. This also declines
    /// harmless symbolic monomial denominators rather than assuming a more
    /// general normalization identity.
    fn epsilon_denominators(expression: AtomView, epsilon: Symbol) -> bool {
        Self::scalar_denominators(expression, epsilon, false)
    }

    /// Inner factors also admit symbolic monomial denominator content. An
    /// integer inverse of `c*x^a*y^b` normalizes to `c^-1*x^-a*y^-b`, so its
    /// numerical content combines identically before or after factorization.
    /// A symbolic sum such as `3*nf-6` has no such normalization rule.
    fn scalar_denominators(
        expression: AtomView,
        epsilon: Symbol,
        symbolic_monomials: bool,
    ) -> bool {
        fn epsilon_only(base: AtomView, epsilon: Symbol) -> bool {
            let mut only = true;
            base.visitor(&mut |factor| {
                only &= match factor {
                    AtomView::Var(variable) => variable.get_symbol() == epsilon,
                    AtomView::Fun(_) => false,
                    _ => true,
                };
                only
            });
            only
        }
        fn monomial_factor(factor: AtomView, epsilon: Symbol) -> bool {
            match factor {
                AtomView::Num(_) | AtomView::Var(_) => true,
                AtomView::Pow(power) => {
                    let (base, exponent) = power.get_base_exp();
                    (matches!(base, AtomView::Var(_)) && get_integer_from_atom(exponent).is_some())
                        || epsilon_only(factor, epsilon)
                }
                _ => epsilon_only(factor, epsilon),
            }
        }
        let mut safe = true;
        expression.visitor(&mut |part| {
            if !safe {
                return false;
            }
            // Polynomial coefficients are numeric atoms whose variables the
            // visitor cannot see. As in RationalSeries, accept only ordinary
            // exact numbers here; generic expansion handles all other domains.
            if let AtomView::Num(number) = part {
                safe = matches!(
                    number.get_coeff_view(),
                    CoefficientView::Natural(..) | CoefficientView::Large(..)
                );
            }
            if let AtomView::Pow(power) = part {
                let (base, exponent) = power.get_base_exp();
                let integer = get_integer_from_atom(exponent);
                if !integer.is_some_and(|n| n >= 0) {
                    safe = epsilon_only(base, epsilon)
                        || (symbolic_monomials
                            && integer.is_some()
                            && match base {
                                AtomView::Mul(product) => product
                                    .iter()
                                    .all(|factor| monomial_factor(factor, epsilon)),
                                _ => monomial_factor(base, epsilon),
                            });
                }
            }
            safe
        });
        safe
    }
}

/// [`expand`]`(prefactor, reduced · spectator, ε, order)` with the factors
/// expanded apart: `prefactor × reduced` through `ε^(order - v_S)` and the
/// spectator through `ε^(order - v_R)`, where `v_S` and `v_R` are the
/// valuations of the spectator and of that inner series, and the truncated
/// series multiplied through `ε^order`. Both valuations are lower bounds
/// (a valuation moves up only on a proof that a coefficient vanishes, and
/// a sum takes the lower valuation of its terms), so each factor is exact
/// through every order the product needs and the result equals the exact
/// expansion of the product.
///
/// `None` (the caller then uses [`expand`] on the product) when a
/// rational denominator of `prefactor` or a master coefficient has
/// non-monomial content in symbols other than ε, when a term of the inner series ends
/// below `ε^(order - v_S)`, which only a tabulated master's unknown tail
/// does, or when the spectator series exceeds its bounds.
pub(crate) fn expand_factored(
    prefactor: &Atom,
    reduced: &Atom,
    spectator: &Spectator,
    epsilon: Symbol,
    order: i64,
) -> Result<Option<Atom>, VakintError> {
    let Some(spectator_valuation) = spectator.0.valuation() else {
        return Ok(Some(Atom::Zero));
    };
    // Gamma functions and scale exponentials are handled by the existing
    // split, not mistaken for symbolic rational denominators.
    if !Spectator::scalar_denominators(split(prefactor, epsilon).rational.as_view(), epsilon, true)
    {
        return Ok(None);
    }
    let inner_order = order - spectator_valuation;
    let symbols = master_symbols();
    let (masters, expanded) = masters_and_expansion(reduced, &symbols)?;
    let mut counts = SeriesCounts::default();
    // Coefficients of the inner series by power of ε.
    let mut inner = std::collections::BTreeMap::<i64, Vec<Atom>>::new();
    for (monomial, coefficient) in expanded.coefficient_list::<i16>(&masters) {
        if !Spectator::scalar_denominators(coefficient.as_view(), epsilon, true) {
            return Ok(None);
        }
        let Some(term) = monomial_term(
            prefactor,
            &monomial,
            &coefficient,
            &symbols,
            epsilon,
            inner_order,
            &mut counts,
            &masters::p_master_expansion,
        )?
        else {
            continue;
        };
        if term.order() < inner_order {
            return Ok(None);
        }
        for (index, coefficient) in term.coefficients.into_iter().enumerate() {
            if !coefficient.is_zero() {
                inner
                    .entry(term.valuation + index as i64)
                    .or_default()
                    .push(coefficient);
            }
        }
    }
    tracing::info!(
        uv = true,
        rstar = true,
        timing = true,
        exact_rational_series = counts.exact,
        generic_rational_series = counts.generic,
        target_order = inner_order,
        master_variables = masters.len(),
        spectator_valuation,
        "Expanded p-integral Laurent coefficients"
    );
    let Some(inner_valuation) = inner.first_key_value().map(|(power, _)| *power) else {
        return Ok(Some(Atom::Zero));
    };
    let inner = Laurent {
        valuation: inner_valuation,
        coefficients: (inner_valuation..=inner_order)
            .map(|power| {
                inner
                    .remove(&power)
                    .map_or(Atom::Zero, |parts| Atom::add_many(parts.iter()).expand())
            })
            .collect(),
    };
    let Some(outer) = ({
        let _stage = crate::stage_timing::StageSpan::new("vakint.p_integral.spectator_series");
        spectator.0.series(order - inner_valuation)
    }) else {
        return Ok(None);
    };
    Ok(Some(inner.mul(&outer, order).to_atom(epsilon).expand()))
}

#[cfg(test)]
mod tests {
    use symbolica::atom::{Atom, AtomCore};
    use symbolica::function;

    use super::super::finalize::{gamma_one_plus, gamma_series, pg_value};
    use super::super::{masters, p_gamma_symbol, p_master_symbol, pg_symbol};
    use super::{expand, valuation};
    use crate::Vakint;
    use crate::utils::vakint_macros::{vk_parse, vk_symbol};

    fn assert_same(lhs: Atom, rhs: Atom) {
        // Symbolica's own series keeps log(𝑒^x).
        let x = Atom::var(crate::symbols::S.x_);
        let rhs = rhs
            .replace(
                function!(
                    symbolica::atom::Symbol::LOG,
                    Atom::var(symbolica::atom::Symbol::E).pow(&x)
                )
                .to_pattern(),
            )
            .with(x.to_pattern())
            .expand();
        let difference = (&lhs - &rhs).expand();
        assert!(
            difference.is_zero(),
            "{lhs}\n!=\n{rhs}\ndifference {difference}"
        );
    }

    #[test]
    fn tabulated_and_g_function_masters_expand_exactly() {
        Vakint::initialize_vakint_symbols();
        let e = vk_symbol!("ε");
        let one = Atom::num(1);
        // A bare tabulated master is its table through the requested order.
        let no = function!(p_master_symbol(), Atom::var(vk_symbol!("no")));
        let table = masters::p_master_expansion("no", e).unwrap();
        assert_same(
            expand(&one, &no, e, 2).unwrap(),
            table
                .series(
                    e,
                    Atom::Zero.as_view(),
                    symbolica::domains::rational::Rational::from(2),
                )
                .unwrap()
                .to_atom(),
        );
        let no1 = function!(p_master_symbol(), Atom::var(vk_symbol!("no1")));
        assert!(!expand(&one, &no1, e, -1).unwrap().is_zero());
        // G-function masters, products and poles against the direct Γ series.
        let pg = |x1: i64, x2: i64| function!(pg_symbol(), Atom::num(x1), Atom::num(x2));
        let ep = Atom::var(e);
        let reduced = (Atom::num(1) - Atom::num(2) * &ep) / &ep * pg(0, 0).pow(Atom::num(2))
            + (Atom::num(3) * &ep - Atom::num(1)) / &ep * pg(0, 0) * pg(1, 0)
            + pg(2, 1) * vk_parse!("mass^2").unwrap();
        let reference = (Atom::num(1) - Atom::num(2) * &ep) / &ep
            * pg_value(0, 0, e).pow(Atom::num(2))
            + (Atom::num(3) * &ep - Atom::num(1)) / &ep * pg_value(0, 0, e) * pg_value(1, 0, e)
            + pg_value(2, 1, e) * vk_parse!("mass^2").unwrap();
        let prefactor = vk_parse!(
            "(𝑖*𝜋^(2-ε)*Gamma(1+ε)*Gamma(1-ε)^2/Gamma(2-2*ε))^2*exp(-logQ2mu-log_mu_sq)^(2*ε)\
             *(2*𝜋)^(-8)*exp(log_mu_sq)^(2*ε)*𝜋^(2*ε)*exp(EulerGamma)^(2*ε)"
        )
        .unwrap();
        let prefactor = Vakint::normalize_gamma_functions(prefactor.as_view());
        for order in [-3, -1, 0] {
            assert_same(
                expand(&prefactor, &reduced, e, order).unwrap(),
                gamma_series(&(&prefactor * &reference), e, order).unwrap(),
            );
        }
        assert_eq!(
            valuation(&(Atom::num(3) / ep.pow(Atom::num(2))), e).unwrap(),
            Some(-2)
        );
    }

    /// The factored expansion is the same atom as the expansion of the
    /// product: masters with spurious poles, tabulated masters and Gamma
    /// tags against spectators with poles, zeros, spectator symbols and a
    /// projector denominator, at several depths.
    #[test]
    fn factored_expansion_is_the_expansion_of_the_product() {
        Vakint::initialize_vakint_symbols();
        let e = vk_symbol!("ε");
        let ep = Atom::var(e);
        let pg = |x1: i64, x2: i64| function!(pg_symbol(), Atom::num(x1), Atom::num(x2));
        let no = function!(p_master_symbol(), Atom::var(vk_symbol!("no")));
        let tag = |m: i64| function!(p_gamma_symbol(), Atom::num(m));
        let masters = [
            (Atom::num(1) - Atom::num(2) * &ep) / ep.pow(Atom::num(2)) * pg(0, 0).pow(Atom::num(2))
                + (Atom::num(3) * &ep - Atom::num(1)) / (Atom::num(2) - Atom::num(3) * &ep)
                    * pg(0, 0)
                    * pg(1, 0)
                + Atom::num(5) / (Atom::num(4) * &ep) * tag(2).pow(Atom::num(-1)) * pg(2, 1),
            Atom::num(7) / ((Atom::num(1) - &ep) * &ep) * &no
                + (Atom::num(1) + &ep) / ep.pow(Atom::num(3)) * tag(-1) * pg(1, 1),
            // ε-denominators with a non-unit constant term.
            vk_parse!("(nf-ε)/(3-2*ε)").unwrap() * pg(1, 0)
                + vk_parse!("ca/((6-4*ε)*ε)").unwrap() * pg(0, 0).pow(Atom::num(2)),
            Atom::Zero,
        ];
        let spectators = [
            Atom::num(1),
            Atom::num(-3) / ep.pow(Atom::num(3)) + Atom::num(2) / &ep,
            vk_parse!("(nf*ε^2-2*ca)/(3-2*ε)").unwrap(),
            vk_parse!("ε*(1+ε)/(1-2*ε)^2+mass^2").unwrap(),
            Atom::num(-1) / Atom::num(3),
            vk_parse!("(nf+3*ca)*(2+𝑖)").unwrap(),
            vk_parse!("mass^2*(nf+ε)/(ε*(3-2*ε))").unwrap(),
        ];
        let prefactor = Vakint::normalize_gamma_functions(
            vk_parse!("(Gamma(1+ε)*Gamma(1-ε)^2/Gamma(2-2*ε))^3*exp(-logMassiveClosure)^(4*ε)*Gamma(2-4*ε)*Gamma(1+4*ε)")
                .unwrap()
                .as_view(),
        );
        for masters in &masters {
            for spectator in &spectators {
                let prepared = super::Spectator::new(spectator, e).unwrap();
                for order in [-3, -1, 0, 1] {
                    let product = (masters * spectator).expand();
                    let old = expand(&prefactor, &product, e, order).unwrap();
                    match super::expand_factored(&prefactor, masters, &prepared, e, order).unwrap()
                    {
                        Some(new) => assert_eq!(
                            old.to_canonical_string(),
                            new.to_canonical_string(),
                            "{masters} · {spectator} at order {order}"
                        ),
                        // Only an unknown master tail inside the range declines.
                        None => assert!(crate::truncation_guard::has_unknown_tail(old.as_view())),
                    }
                }
            }
        }
        // Γ functions in the spectator are not admitted.
        assert!(super::Spectator::new(&vk_parse!("Gamma(1+ε)/ε").unwrap(), e).is_none());
    }

    #[test]
    fn factored_expansion_declines_symbolic_denominators() {
        Vakint::initialize_vakint_symbols();
        let e = vk_symbol!("ε");
        let masters = vk_parse!("PG(0,0)^2/(3-2*ε)+PG(1,0)/(ε*(6-4*ε))").unwrap();
        let prefactor = Atom::num(1);
        // Without the guard the separately normalized denominator content
        // gives algebraically equal values but different canonical atoms.
        let rational_spectators = [
            vk_parse!("(nf+1)/(nf-2)").unwrap(),
            vk_parse!("(nf+1)/(1-nf)").unwrap(),
        ];
        for spectator in &rational_spectators {
            assert!(super::Spectator::new(spectator, e).is_none(), "{spectator}");
            let unguarded = super::Spectator(
                super::rational_series::RationalSeries::new(spectator, e).unwrap(),
            );
            let old = expand(&prefactor, &(&masters * spectator).expand(), e, -1).unwrap();
            let new = super::expand_factored(&prefactor, &masters, &unguarded, e, -1)
                .unwrap()
                .unwrap();
            assert_ne!(old.to_canonical_string(), new.to_canonical_string());
            assert!((&old - &new).together().is_zero());
        }
        for expression in [
            "ε/(3-2*ε)+ca/(nf-2*ca)",
            "(nf+3*ca)*nf^-2*(2+𝑖)",
            "(nf+ε)/(mass^2*ε*(3-2*ε))",
            "1/f(nf)",
        ] {
            let expression = vk_parse!(expression).unwrap();
            assert!(super::Spectator::new(&expression, e).is_none());
        }
        for expression in [
            "ε/(3-2*ε)+ca/(nf-2*ca)",
            "1/f(nf)",
            "1/(nf*(nf-2))",
            "1/(nf+ε)",
            "nf^(-1/2)",
            "1/(nf^-2+ca^-1)",
        ] {
            let expression = vk_parse!(expression).unwrap();
            let spectator = super::Spectator::new(&vk_parse!("1/(3-2*ε)").unwrap(), e).unwrap();
            assert!(
                super::expand_factored(&prefactor, &(&masters * &expression), &spectator, e, -1)
                    .unwrap()
                    .is_none()
            );
            assert!(
                super::expand_factored(&expression, &masters, &spectator, e, -1)
                    .unwrap()
                    .is_none()
            );
        }
        let packed = vk_parse!("(nf+1)/(nf-2)")
            .unwrap()
            .set_coefficient_ring(vk_symbol!("nf"));
        assert!(matches!(
            packed.as_view(),
            symbolica::atom::AtomView::Num(_)
        ));
        assert!(!super::Spectator::epsilon_denominators(packed.as_view(), e));
        assert!(super::Spectator::new(&packed, e).is_none());
        let spectator = super::Spectator::new(&Atom::num(1), e).unwrap();
        assert!(
            super::expand_factored(&packed, &masters, &spectator, e, -1)
                .unwrap()
                .is_none()
        );
        assert!(
            super::expand_factored(&prefactor, &(&masters * packed), &spectator, e, -1)
                .unwrap()
                .is_none()
        );
    }

    #[test]
    fn factored_expansion_preserves_inner_symbolic_monomial_denominators() {
        Vakint::initialize_vakint_symbols();
        let e = vk_symbol!("ε");
        let masters = vk_parse!("PG(0,0)^2/(3-2*ε)+PG(1,0)/(ε*(6-4*ε))").unwrap();
        let one = Atom::num(1);
        let spectators = [
            vk_parse!("(nf+1)/(ε*(3-2*ε))").unwrap(),
            vk_parse!("mass^2*nf^2/(ε*(3-2*ε))").unwrap(),
        ];
        // Include the pi denominator of the production MSbar normalization,
        // symbolic masses, numeric content and the review's complex monomial.
        for factor in [
            "1/(1024*𝜋^10)",
            "1/(3*mass^2*nf^2)",
            "(2+𝑖)*nf^-2",
            "1/(mass^2*(3-2*ε))",
            "1/mass^2+1/nf^2",
        ] {
            let factor = vk_parse!(factor).unwrap();
            assert!(super::Spectator::new(&factor, e).is_none());
            assert!(super::Spectator::scalar_denominators(
                factor.as_view(),
                e,
                true
            ));
            for spectator in &spectators {
                let prepared = super::Spectator::new(spectator, e).unwrap();
                for (prefactor, reduced) in [(&factor, masters.clone()), (&one, &factor * &masters)]
                {
                    for order in [-1, 0, 1] {
                        let old =
                            expand(prefactor, &(&reduced * spectator).expand(), e, order).unwrap();
                        let new = super::expand_factored(prefactor, &reduced, &prepared, e, order)
                            .unwrap()
                            .unwrap();
                        assert_eq!(
                            old.to_canonical_string(),
                            new.to_canonical_string(),
                            "{factor} times {spectator} at order {order}"
                        );
                    }
                }
            }
        }
    }

    #[test]
    fn gamma_tags_expand_exactly_to_any_depth() {
        Vakint::initialize_vakint_symbols();
        let e = vk_symbol!("ε");
        let ep = Atom::var(e);
        let one = Atom::num(1);
        let tag = |m: i64| p_gamma_symbol().call((m,));
        let pg = |x1: i64, x2: i64| function!(pg_symbol(), Atom::num(x1), Atom::num(x2));
        let gamma = |m: i64| gamma_one_plus(m, &ep);
        // Products of tag powers (negative ones included), with poles, a G-scheme
        // constant and a coefficient that is not a number.
        let reduced = tag(1).pow(Atom::num(2)) * tag(-2).pow(Atom::num(-1)) / ep.pow(Atom::num(2))
            + (&one - Atom::num(2) * &ep) / &ep * tag(3) * tag(-1).pow(Atom::num(-3)) * pg(1, 0)
            + tag(2) * tag(-3).pow(Atom::num(2)) * vk_parse!("mass^2").unwrap()
            + tag(4).pow(Atom::num(-2)) * Atom::num(7);
        let reference = gamma(1).pow(Atom::num(2)) * gamma(-2).pow(Atom::num(-1))
            / ep.pow(Atom::num(2))
            + (&one - Atom::num(2) * &ep) / &ep
                * gamma(3)
                * gamma(-1).pow(Atom::num(-3))
                * pg_value(1, 0, e)
            + gamma(2) * gamma(-3).pow(Atom::num(2)) * vk_parse!("mass^2").unwrap()
            + gamma(4).pow(Atom::num(-2)) * Atom::num(7);
        let prefactor = Vakint::normalize_gamma_functions(
            vk_parse!("(𝑖*𝜋^(2-ε)*Gamma(1+ε)*Gamma(1-ε)^2/Gamma(2-2*ε))^3")
                .unwrap()
                .as_view(),
        );
        for order in [-2, 0, 3, 6] {
            assert_same(
                expand(&prefactor, &reduced, e, order).unwrap(),
                gamma_series(&(&prefactor * &reference), e, order).unwrap(),
            );
        }
        // A tag times a tabulated master: the tag series times the table.
        let t105 = function!(p_master_symbol(), Atom::var(vk_symbol!("t105")));
        let table = masters::p_master_expansion("t105", e).unwrap();
        assert_same(
            expand(&one, &(tag(-2).pow(Atom::num(3)) * &t105), e, 6).unwrap(),
            gamma_series(&(gamma(-2).pow(Atom::num(3)) * table), e, 6).unwrap(),
        );
        // Malformed tags fail closed.
        for bad in [p_gamma_symbol().call((0,)), p_gamma_symbol().call((1, 2))] {
            assert!(expand(&one, &bad, e, 0).is_err(), "{bad}");
        }
    }

    #[test]
    fn gamma_tag_powers_expand_over_the_whole_i16_range() {
        Vakint::initialize_vakint_symbols();
        let e = vk_symbol!("ε");
        let ep = Atom::var(e);
        let one = Atom::num(1);
        let tag = p_gamma_symbol().call((1,));
        for power in [127, 128, 200, -128, -129, -200] {
            assert_same(
                expand(&one, &tag.pow(Atom::num(power)), e, 2).unwrap(),
                gamma_series(&gamma_one_plus(1, &ep).pow(Atom::num(power)), e, 2).unwrap(),
            );
        }
        // The ends of the range, by the exponential series Γ(1+ε) = exp(-γε + O(ε²)).
        for power in [i64::from(i16::MAX), i64::from(i16::MIN)] {
            let value = expand(&one, &tag.pow(Atom::num(power)), e, 1).unwrap();
            let gamma_e = Atom::var(symbolica::transcendental::euler_gamma());
            assert_same(value, &one - Atom::num(power) * gamma_e * &ep);
        }
    }

    #[test]
    fn anything_but_a_laurent_polynomial_in_the_tags_fails_closed() {
        Vakint::initialize_vakint_symbols();
        let e = vk_symbol!("ε");
        let ep = Atom::var(e);
        let one = Atom::num(1);
        let tag = p_gamma_symbol().call((1,));
        let no = function!(p_master_symbol(), Atom::var(vk_symbol!("no")));
        let foreign_no = function!(
            p_master_symbol(),
            Atom::var(symbolica::symbol!("other::no"))
        );
        for bad in [
            tag.pow(Atom::num(32768)),
            tag.pow(Atom::num(-32769)),
            tag.pow(Atom::num(40000)),
            (tag.pow(Atom::num(30000)) + &one).pow(Atom::num(2)),
            tag.pow(Atom::num(1) / Atom::num(2)),
            tag.pow(ep.clone()),
            &one / (&one + &tag),
            function!(vk_symbol!("foo"), &tag),
            no.pow(Atom::num(-1)),
            &one / (&one + &no),
            pg_symbol().call((1, 0)).pow(ep.clone()),
            foreign_no,
            function!(p_master_symbol(), Atom::var(vk_symbol!("no")), 2),
        ] {
            let error = expand(&one, &bad, e, 2).expect_err(&bad.to_string());
            assert!(
                matches!(error, crate::VakintError::SymbolicaError(_)),
                "{bad}: {error}"
            );
        }
    }

    #[test]
    fn msbar_cancels_the_scale_logarithm_in_the_exponent() {
        Vakint::initialize_vakint_symbols();
        let e = vk_symbol!("ε");
        let settings = crate::VakintSettings {
            integral_normalization_factor: crate::LoopNormalizationFactor::MSbar,
            ..Default::default()
        };
        let normalization = settings
            .get_integral_normalization_factor_atom()
            .unwrap()
            .replace(crate::symbols::S.n_loops.to_pattern())
            .with(Atom::num(1).to_pattern());
        let measure =
            vk_parse!("(𝑖*𝜋^(2-ε)*Gamma(1+ε)*Gamma(1-ε)^2/Gamma(2-2*ε))*exp(-logQ2mu-log_mu_sq)^ε")
                .unwrap();
        let prefactor = Vakint::normalize_gamma_functions((normalization * measure).as_view());
        let split = super::split(&prefactor, e);
        let exponent = split.exponent(2);
        assert!(
            !exponent[0].contains_symbol(vk_symbol!("log_mu_sq")),
            "{exponent:?}"
        );
    }

    /// Lane laurent-depth (reviewer O2 "orders", 2026-11-18): a table's tail
    /// `ε^k·Oep(k, M)` stands for every order from `k` on. Two terms of the
    /// same master whose leading unknown orders cancel must not lose the
    /// orders after it: `ε^-4·M·(ε·PG(0,0) - 1 + ε^9) = ε^5·M` exactly, and
    /// `ε^-4·M·(ε·PG(0,0) - 1 + ε) = ε^-3·M` is refused like `ε^-3·M`.
    #[test]
    fn cancelling_leading_unknown_orders_keep_the_later_ones() {
        Vakint::initialize_vakint_symbols();
        let e = vk_symbol!("ε");
        let one = Atom::num(1);
        let no = function!(p_master_symbol(), Atom::var(vk_symbol!("no")));
        let pg = function!(pg_symbol(), Atom::num(0), Atom::num(0));
        let ep = Atom::var(e);
        let pow = |n: i64| ep.pow(Atom::num(n));
        // The orders follow from the table: valuation v, first unknown k.
        let table =
            super::TabulatedMaster::new(&masters::p_master_expansion("no", e).unwrap(), e).unwrap();
        let (v, k) = (table.valuation, table.tail.as_ref().unwrap().0);
        let n = 4;
        // ε^(n-v)·no needs only the leading (known) order of no.
        let control = expand(&one, &(pow(n - v) * &no), e, n).unwrap();
        assert!(!control.is_zero());
        assert!(!crate::truncation_guard::has_unknown_tail(
            control.as_view()
        ));
        // ε^(n-k-1)·no·(ε·PG(0,0) - 1 + ε^(k+1-v)) = ε^(n-v)·no: each term
        // needs no through k+1, the unknown orders cancel exactly.
        let cancel =
            (pow(n - k - 1) * &no * (&ep * &pg - (Atom::num(1) - pow(k + 1 - v)))).expand();
        let value = expand(&one, &cancel, e, n).unwrap();
        assert_eq!(value.to_canonical_string(), control.to_canonical_string());
        let guard = |value: &Atom| {
            crate::truncation_guard::ensure_known_orders(
                value,
                e,
                Some(n),
                "test",
                crate::truncation_guard::plain_name,
            )
        };
        // ε^(n-k)·no needs no's first unknown order k at ε^n.
        let alarm = guard(&expand(&one, &(pow(n - k) * &no), e, n).unwrap()).unwrap_err();
        let cancel_missing = (pow(n - k - 1) * &no * (&ep * &pg - (Atom::num(1) - &ep))).expand();
        let missing = guard(&expand(&one, &cancel_missing, e, n).unwrap()).unwrap_err();
        assert_eq!(missing, alarm);
        assert_eq!(alarm.masters.len(), 1);
        assert_eq!(
            (
                alarm.masters[0].first_missing_order,
                alarm.masters[0].last_needed_order,
                alarm.masters[0].first_affected_result_order
            ),
            (k, Some(k), n)
        );
        // Two orders deeper: one entry, k..k+2.
        let deep = guard(&expand(&one, &(pow(n - k - 2) * &no), e, n).unwrap()).unwrap_err();
        assert_eq!(deep.masters.len(), 1);
        assert_eq!(
            (
                deep.masters[0].first_missing_order,
                deep.masters[0].last_needed_order,
                deep.masters[0].first_affected_result_order
            ),
            (k, Some(k + 2), n - 2)
        );
        // A master squared needs each factor through k + 1 here.
        let square = (pow(n - k - v - 1) * no.pow(Atom::num(2))).expand();
        let squared = guard(&expand(&one, &square, e, n).unwrap()).unwrap_err();
        assert_eq!(squared.masters[0].first_missing_order, k);
        assert_eq!(squared.masters[0].last_needed_order, Some(k + 1));
    }

    /// Property: for random small sums of products of rational functions of
    /// ε, PG factors and tabulated masters with unknown tails (with terms
    /// built to cancel), the expansion through `ε^N` is the exact product
    /// truncated at `N`, where the reference writes every unknown master
    /// coefficient as its own symbol `Oep(j, M)` far beyond any needed
    /// order; the guard fires exactly when an unknown order survives.
    #[test]
    fn products_with_unknown_tails_expand_exactly() {
        Vakint::initialize_vakint_symbols();
        let e = vk_symbol!("ε");
        let ep = Atom::var(e);
        let pow = |n: i64| ep.pow(Atom::num(n));
        let mut state: u64 = 0x5eed_1a0e_0de9_7a11;
        let mut next = move |n: u64| {
            state = state
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            ((state >> 33) % n) as i64
        };
        // Three masters: (valuation, first unknown order k), known
        // coefficients random integers (some zero), or fully exact (no tail).
        let shapes = [("ma", -2, 0), ("mb", 0, 2), ("mc", -1, 2), ("md", -1, -1)];
        let mut tables = std::collections::BTreeMap::<String, (Atom, Atom)>::new();
        for (name, valuation, k) in shapes {
            let master = function!(p_master_symbol(), Atom::var(vk_symbol!(name)));
            let last = if k < 0 { 2 } else { k - 1 };
            let mut known = Atom::Zero;
            for j in valuation..=last {
                known += Atom::num(next(7) - 2) * pow(j);
            }
            known += pow(valuation) * Atom::num(1 + next(3));
            let (table, reference) = if k < 0 {
                (known.clone(), known)
            } else {
                let mut reference = known.clone();
                // Beyond any order a trial needs: ε^order (<= 1) minus the
                // lowest valuation of the other factors (>= -9).
                for j in k..=10 {
                    reference += pow(j) * crate::truncation_guard::unknown_tail_marker(j, &master);
                }
                (
                    known + crate::truncation_guard::unknown_tail(e, k, &master),
                    reference,
                )
            };
            tables.insert(name.to_string(), (table, reference));
        }
        let table = |name: &str, epsilon: symbolica::atom::Symbol| {
            assert_eq!(epsilon, e);
            Ok(tables[name].0.clone())
        };
        let master = |name: &str| function!(p_master_symbol(), Atom::var(vk_symbol!(name)));
        let names = ["ma", "mb", "mc", "md"];
        let pg = |x1: i64, x2: i64| function!(pg_symbol(), Atom::num(x1), Atom::num(x2));
        let mut fired = 0;
        let mut clean = 0;
        for trial in 0..64 {
            let mut reduced = Atom::Zero;
            for _ in 0..1 + next(3) {
                // coefficient ε^p·(a + bε + cε²)/(1 - dε)
                let coefficient = pow(next(5) - 2)
                    * (Atom::num(1 + next(3))
                        + Atom::num(next(5) - 2) * &ep
                        + Atom::num(next(5) - 2) * pow(2))
                    / (Atom::num(1) - Atom::num(next(3)) * &ep);
                // Products of 1 to 4 master factors (repeats = powers).
                let mut monomial = master(names[next(4) as usize]);
                for _ in 0..next(4) {
                    monomial *= master(names[next(4) as usize]);
                }
                match next(4) {
                    0 => monomial *= pg(0, 0),
                    1 => monomial *= pg(1, 0),
                    _ => {}
                }
                let term = &coefficient * &monomial;
                reduced += &term;
                // A partner cancelling the leading orders of the term
                // (ε·PG(0,0) = 1 exactly), but not all of them.
                if next(2) == 1 {
                    let shift = next(4);
                    reduced -=
                        &coefficient * &monomial * &ep * pg(0, 0) * (Atom::num(1) - pow(1 + shift));
                }
            }
            let reduced = reduced.expand();
            let order = next(3) - 1;
            let prefactor = if trial % 3 == 0 {
                Vakint::normalize_gamma_functions(
                    vk_parse!("Gamma(1+ε)*Gamma(1-ε)^2/Gamma(2-2*ε)*exp(-logQ2mu)^ε")
                        .unwrap()
                        .as_view(),
                )
            } else {
                Atom::num(1)
            };
            let started = std::time::Instant::now();
            let value = super::expand_with(&prefactor, &reduced, e, order, &table).unwrap();
            let expanded = started.elapsed().as_secs_f64();
            let mut exact = reduced.clone();
            for name in names {
                exact = exact
                    .replace(master(name).to_pattern())
                    .with(tables[name].1.to_pattern());
            }
            for (x1, x2) in [(0, 0), (1, 0)] {
                exact = exact
                    .replace(pg(x1, x2).to_pattern())
                    .with(pg_value(x1, x2, e).to_pattern());
            }
            let reference = gamma_series(&(&prefactor * &exact), e, order).unwrap();
            eprintln!(
                "trial {trial}: expansion {expanded:.2} s, reference {:.2} s",
                started.elapsed().as_secs_f64() - expanded
            );
            assert_same(value.clone(), reference.clone());
            let guard = crate::truncation_guard::ensure_known_orders(
                &value,
                e,
                Some(order),
                "test",
                crate::truncation_guard::plain_name,
            );
            let unknown = crate::truncation_guard::has_unknown_tail(reference.expand().as_view());
            assert_eq!(
                guard.is_err(),
                unknown,
                "trial {trial}: {reduced} through ε^{order}"
            );
            if unknown {
                fired += 1;
            } else {
                clean += 1;
            }
        }
        // Both outcomes are exercised.
        assert!(fired > 0 && clean > 0, "fired {fired}, clean {clean}");
    }
}
