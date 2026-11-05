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
use symbolica::domains::rational::Rational;
use symbolica::function;
use symbolica::poly::series::SeriesDepth;
use symbolica::transcendental::{euler_gamma, gamma, zeta};

use super::{masters, p_gamma_symbol, p_master_symbol, pg_symbol};
use crate::{VakintError, get_integer_from_atom};

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

/// The Laurent expansion through `ε^order` of `prefactor × reduced`, where
/// `reduced` is a polynomial in the masters `PG(x1,x2)` (replaced by
/// `pg_value`) and `PMaster(name)` (replaced by its tabulated expansion and
/// unknown tail) and a Laurent polynomial in the Gamma tags `PGamma(m)`, every
/// `PGamma(m)^e = Γ(1+mε)^e` with an `i16` exponent `e` expanded exactly like
/// the Γ functions of the prefactor, and `prefactor` is a product of rational
/// factors, Γ functions and exponentials. Any other dependence on these
/// symbols is an error, never a result that still contains them.
pub(crate) fn expand(
    prefactor: &Atom,
    reduced: &Atom,
    epsilon: Symbol,
    order: i64,
) -> Result<Atom, VakintError> {
    let mut masters = Vec::new();
    reduced.visitor(&mut |part| {
        if let AtomView::Fun(function) = part
            && [pg_symbol(), p_master_symbol(), p_gamma_symbol()].contains(&function.get_symbol())
        {
            masters.push(part.to_owned());
            return false;
        }
        true
    });
    masters.sort();
    masters.dedup();
    let symbols = [pg_symbol(), p_master_symbol(), p_gamma_symbol()];
    let expanded = reduced.expand();
    if let Some(power) = unsupported_power(&expanded, &symbols) {
        return Err(failure(format!("{power} is not an i16 integer power")));
    }
    let mut total = Atom::Zero;
    for (monomial, coefficient) in expanded.coefficient_list::<i16>(&masters) {
        if symbols.iter().any(|s| coefficient.contains_symbol(*s)) {
            return Err(failure(format!(
                "{coefficient} is not free of PG, PMaster and PGamma"
            )));
        }
        // The Γ form of the PG factors joins the prefactor; PMaster factors
        // are tabulated Laurent polynomials.
        let mut product = prefactor * &coefficient;
        let mut tabulated: Option<Laurent> = None;
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
                && function.get_symbol() == p_gamma_symbol()
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
            if function.get_symbol() == pg_symbol() {
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
                let value = Laurent::from_polynomial(
                    &masters::p_master_expansion(&name, epsilon)?,
                    epsilon,
                )?;
                for _ in 0..power {
                    tabulated = Some(match tabulated {
                        None => value.clone(),
                        Some(tabulated) => tabulated.mul(&value, i64::MAX / 4),
                    });
                }
            }
        }
        let split = split(
            &crate::Vakint::normalize_gamma_functions(product.as_view()),
            epsilon,
        );
        let Some(rational_valuation) = valuation(&split.rational, epsilon)? else {
            continue;
        };
        let tabulated_valuation = tabulated.as_ref().map_or(0, |t| t.valuation);
        let depth = order - rational_valuation - tabulated_valuation;
        if depth < 0 {
            continue;
        }
        let rational = Laurent::series(&split.rational, epsilon, order - tabulated_valuation)?;
        let exponential = Laurent::exp(&split.exponent(depth), depth);
        let mut term = rational.mul(&exponential, order - tabulated_valuation);
        if let Some(mut tabulated) = tabulated {
            // Through ε^(order - rational_valuation); a tabulated unknown
            // tail inside that range survives as its Oep marker.
            tabulated.coefficients.truncate(depth as usize + 1);
            term = term.mul(&tabulated, order);
        }
        total += term.to_atom(epsilon);
    }
    Ok(total.expand())
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
}
