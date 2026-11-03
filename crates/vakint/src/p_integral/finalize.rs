//! The shared p-integral finalizer: normalization, Laurent expansion and
//! master substitution of [`ReducedPropagatorValue`]s, for every backend.
//!
//! With `Q² = -dot(p,p)` and the per-loop normalization `F` of the settings,
//!
//! ```text
//! I = Σ_w [F·i·π^(d/2)·εG(1,1)]^L · (Q²)^w · (Q²/μ²)^(-Lε) · (μ²)^(-Lε) · reduced_w
//! ```
//!
//! (`docs/p_integrals.md`, section 4). The logarithms of `(Q²/μ²)^(-Lε)` and
//! `(μ²)^(-Lε)` stay separate placeholders until the end, mirroring the vacuum
//! `logmUVmu`/`log_mu_sq` split, so that normalizations carrying
//! `exp(log_mu_sq)^(εL)` (MSbar, GScheme) cancel `(μ²)^(-Lε)` symbolically.

use symbolica::atom::{Atom, AtomCore, AtomView, Symbol};
use symbolica::domains::rational::Rational;
use symbolica::function;
use symbolica::id::Replacement;
use symbolica::transcendental::{euler_gamma, gamma, zeta};

use super::{ReducedPropagatorValue, external_head, masters, p_master_symbol, pg_symbol};
use crate::rustred_evaluation::RustRedMasterMode;
use crate::symbols::S;
use crate::utils::simplify_real;
use crate::utils::vakint_macros::{vk_parse, vk_symbol};
use crate::{Vakint, VakintError, VakintSettings, get_integer_from_atom};

/// Placeholder for `log(Q²/μ²) = log(-dot(p,p)/mursq)`.
fn log_q_mu() -> Symbol {
    vk_symbol!("logQ2mu")
}

/// Placeholder for `log(μ²)`, as in the normalization factors.
fn log_mu() -> Symbol {
    vk_symbol!("log_mu_sq")
}

/// `Γ(1 + c·ε)` with a rational `c` (`Γ(1)` stays exact as `1`).
fn gamma_one_plus(c: i64, epsilon: &Atom) -> Atom {
    gamma().call(Atom::num(1) + Atom::num(c) * epsilon)
}

/// `PG(x1, x2) = (1+x1+x2)·G(1+x1ε, 1+x2ε)/(ε·G(1,1))`, Forcer's
/// `GschemeConstants(x1,x2)`, written with `Γ(1 + cε)` only:
///
/// ```text
/// Γ(1+(1+s)ε) Γ(1-(1+x1)ε) Γ(1-(1+x2)ε) (1-2ε) Γ(1-2ε)
/// / (ε Γ(1+x1ε) Γ(1+x2ε) (1-(1+s)ε) Γ(1-(2+s)ε) Γ(1+ε) Γ(1-ε)²),   s = x1+x2.
/// ```
pub(crate) fn pg_value(x1: i64, x2: i64, epsilon: Symbol) -> Atom {
    let e = Atom::var(epsilon);
    let s = x1 + x2;
    let numerator = gamma_one_plus(1 + s, &e)
        * gamma_one_plus(-(1 + x1), &e)
        * gamma_one_plus(-(1 + x2), &e)
        * (Atom::num(1) - Atom::num(2) * &e)
        * gamma_one_plus(-2, &e);
    let denominator = &e
        * gamma_one_plus(x1, &e)
        * gamma_one_plus(x2, &e)
        * (Atom::num(1) - Atom::num(2 + s) * &e)
        * gamma_one_plus(-(2 + s), &e)
        * gamma_one_plus(1, &e)
        * gamma_one_plus(-1, &e).pow(Atom::num(2));
    numerator / denominator
}

/// Marks the first omitted order of a truncated `log Γ(1+cε)`.
fn gamma_tail() -> Symbol {
    vk_symbol!("GammaSeriesTail")
}

/// `n + c·ε` of a linear Γ argument, with integer `n` and rational `c`.
fn linear_argument(argument: AtomView, epsilon: Symbol) -> Option<(i64, Atom)> {
    let argument = argument.to_owned().expand();
    let e = Atom::var(epsilon);
    let slope = argument.coefficient(e.as_view()).expand();
    let constant = argument
        .replace(e.to_pattern())
        .with(Atom::Zero.to_pattern())
        .expand();
    let n = get_integer_from_atom(constant.as_view())?;
    (matches!(slope.as_view(), AtomView::Num(_))
        && (Atom::num(n) + &slope * &e - &argument).expand().is_zero())
    .then_some((n, slope))
}

/// Every `Γ(n + cε)` with integer `n` and rational `c` written as a rational
/// function of `ε` times `Γ(1+cε) = exp(-γcε + Σ_(k=2..=depth) (-1)^k ζ(k)(cε)^k/k)`,
/// the exponent ending in `ε^(depth+1)·GammaSeriesTail`. Symbolica then
/// expands the exponentials; other Γ functions are left unchanged.
pub(crate) fn expand_gamma_functions(expression: &Atom, epsilon: Symbol, depth: i64) -> Atom {
    let e = Atom::var(epsilon);
    expression
        .replace(function!(gamma(), Atom::var(S.x_)).to_pattern())
        .with_map(move |matched| {
            let argument = matched.get(S.x_).unwrap().to_atom();
            let Some((n, c)) = linear_argument(argument.as_view(), epsilon) else {
                return function!(gamma(), argument);
            };
            let z = &c * &e;
            let mut rational = Atom::num(1);
            for m in 1..n {
                rational *= Atom::num(m) + &z;
            }
            for m in n..=0 {
                rational /= Atom::num(m) + &z;
            }
            if c.is_zero() {
                return rational;
            }
            let mut exponent = -Atom::var(euler_gamma()) * &z;
            for k in 2..=depth {
                let sign = if k % 2 == 0 { 1 } else { -1 };
                exponent += Atom::num(sign) * zeta().call(Atom::num(k)) * z.pow(Atom::num(k))
                    / Atom::num(k);
            }
            exponent += e.pow(Atom::num(depth + 1)) * Atom::var(gamma_tail());
            rational * function!(Symbol::EXP, exponent)
        })
}

/// The Laurent series of `expression` through `ε^order`, with every linear
/// Γ function expanded deep enough: the depth grows until no truncated
/// `log Γ` tail reaches the requested order.
pub(crate) fn gamma_series(
    expression: &Atom,
    epsilon: Symbol,
    order: i64,
) -> Result<Atom, VakintError> {
    let expression = Vakint::normalize_gamma_functions(expression.as_view());
    let mut depth = order.max(0) + 8;
    loop {
        let series = expand_gamma_functions(&expression, epsilon, depth)
            .series(epsilon, Atom::Zero.as_view(), Rational::from(order))
            .map_err(|error| VakintError::SymbolicaError(error.to_string()))?
            .to_atom()
            .expand();
        if !series.contains_symbol(gamma_tail()) {
            return Ok(series);
        }
        if depth > 256 {
            return Err(VakintError::EvaluationError(
                "Γ expansion beyond depth 256 in a p-integral".into(),
            ));
        }
        depth *= 2;
    }
}

/// Replace every `PG(x1,x2)` by its exact Γ form and every `PMaster(name)`
/// by its tabulated expansion with an unknown tail.
fn substitute_masters(expression: &Atom, epsilon: Symbol) -> Result<Atom, VakintError> {
    let mut replacements = Vec::new();
    let mut error = None;
    expression.visitor(&mut |part| {
        let AtomView::Fun(function) = part else {
            return true;
        };
        if function.get_symbol() == pg_symbol() {
            let arguments = function
                .iter()
                .map(get_integer_from_atom)
                .collect::<Option<Vec<_>>>();
            match arguments.as_deref() {
                Some([x1, x2]) => {
                    replacements.push((part.to_owned(), pg_value(*x1, *x2, epsilon)));
                }
                _ => error = Some(VakintError::EvaluationError(format!("invalid {part}"))),
            }
            return false;
        }
        if function.get_symbol() == p_master_symbol() {
            let name = match function.iter().next() {
                Some(AtomView::Var(name)) if function.get_nargs() == 1 => {
                    name.get_symbol().get_stripped_name().to_string()
                }
                _ => {
                    error = Some(VakintError::EvaluationError(format!("invalid {part}")));
                    return false;
                }
            };
            match masters::p_master_expansion(&name, epsilon) {
                Ok(value) => replacements.push((part.to_owned(), value)),
                Err(e) => error = Some(e),
            }
            return false;
        }
        true
    });
    if let Some(error) = error {
        return Err(error);
    }
    replacements.sort();
    replacements.dedup();
    Ok(expression.replace_multiple(
        &replacements
            .into_iter()
            .map(|(source, target)| Replacement::new(source.to_pattern(), target))
            .collect::<Vec<_>>(),
    ))
}

/// Expand the arguments of all exponentials so that merged factors such as
/// `exp(εL·log_mu_sq)·exp(-εL·(logQ2mu+log_mu_sq))` cancel exactly.
fn expand_exponents(expression: &Atom) -> Atom {
    expression
        .replace(function!(Symbol::EXP, Atom::var(S.x_)).to_pattern())
        .with_map(|matched| {
            function!(Symbol::EXP, matched.get(S.x_).unwrap().to_atom().expand())
        })
}

/// The finalized p-integral `Σ_w ...` of the module documentation, in
/// `mode`: `Retain` keeps the exact `ε` dependence and the masters `PG` and
/// `PMaster`; `Symbolic` and `Numerical` expand through
/// `ε^(number_of_terms_in_epsilon_expansion - L - 1)` with `PG` expanded
/// exactly and `PMaster` substituted from the Forcer table (constants stay
/// symbolic; the numerical boundary of the caller evaluates them). A surviving
/// unknown master tail is a [`crate::truncation_guard::MissingMasterOrders`]
/// error naming `backend`.
pub(crate) fn finalize(
    values: &[ReducedPropagatorValue],
    loops: usize,
    external: &Atom,
    settings: &VakintSettings,
    mode: RustRedMasterMode,
    backend: &str,
) -> Result<Atom, VakintError> {
    let epsilon = vk_symbol!(&settings.epsilon_symbol);
    let (head, _) = external_head(external)?;
    let square = function!(S.dot, &head, &head);
    let mut reduced = Atom::Zero;
    for value in values {
        let sign = if value.scale_power % 2 == 0 { 1 } else { -1 };
        reduced += Atom::num(sign) * square.pow(Atom::num(value.scale_power)) * &value.reduced;
    }
    if reduced.is_zero() {
        return Ok(Atom::Zero);
    }
    let l = loops as i64;
    let normalization = settings
        .get_staged_integral_normalization_factor_atom()?
        .replace(S.n_loops.to_pattern())
        .with(Atom::num(l).to_pattern());
    let measure = vk_parse!(&format!(
        "(𝑖*𝜋^(2-({eps}))*Gamma(1+{eps})*Gamma(1-{eps})^2/Gamma(2-2*{eps}))^{l}\
         *exp(-logQ2mu-log_mu_sq)^({eps}*{l})",
        eps = settings.epsilon_symbol
    ))
    .map_err(|error| VakintError::SymbolicaError(error.to_string()))?;
    let total = Vakint::normalize_gamma_functions((normalization * measure * reduced).as_view());

    let mut result = if mode == RustRedMasterMode::Retain {
        expand_exponents(&simplify_real(total.as_view()))
    } else {
        let target = settings.number_of_terms_in_epsilon_expansion - l - 1;
        let expanded = gamma_series(&substitute_masters(&total, epsilon)?, epsilon, target)?;
        crate::truncation_guard::ensure_known_orders(
            &expanded,
            epsilon,
            Some(target),
            backend,
            crate::truncation_guard::plain_name,
        )?;
        simplify_real(expanded.as_view())
    };
    if !settings.use_dot_product_notation {
        result = Vakint::convert_from_dot_notation(result.as_view());
    }
    let mursq = Atom::var(vk_symbol!(&settings.mu_r_sq_symbol));
    let log_q = function!(Symbol::LOG, -square / &mursq);
    Ok(result
        .replace(Atom::var(log_q_mu()).to_pattern())
        .with(log_q.to_pattern())
        .replace(Atom::var(log_mu()).to_pattern())
        .with(function!(Symbol::LOG, mursq).to_pattern()))
}

#[cfg(test)]
pub(crate) mod tests;
