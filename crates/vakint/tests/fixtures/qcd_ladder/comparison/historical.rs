//! Complete scalar identities against the preserved q0 result, without
//! rationalizing a polynomial whose variables include all formal masters.

use super::{
    Atom, AtomCore, AtomView, CompletedRun, Result, canonical_retained,
    convert_historical_retained, sha256, vakint_parse, vakint_symbol, write_identity,
};
use rustred::algebra::{Coefficient, CoefficientContext, ExactAlgebraLimits};
use std::collections::{BTreeMap, BTreeSet};
use std::path::Path;
use std::sync::Arc;
use std::time::Instant;
use symbolica::domains::float::{Float, RoundingDirection};
use symbolica::domains::integer::Integer;
use symbolica::domains::rational::Rational;
use symbolica::poly::{PolyVariable, PolynomialConversionError};
use symbolica::prelude::{Q, Z};

const HISTORICAL_BYTES: usize = 6_683_821;
const HISTORICAL_SHA256: &str = "da909dc0b8f199a9ce46c27368bd2d1c88b277d72e765ab70f9a83680caa1212";
const INPUT_SHA256: &str = "9814c3b49a1b486b9790bbe5458d8bc582c3fb3c299ba6a1a0362db8ba71cb16";

struct Comparison {
    converted: Atom,
    factor: Atom,
    rows: Vec<serde_json::Value>,
    master_count: usize,
    masters: Vec<serde_json::Value>,
}

/// Leading Laurent terms shown per master-product coefficient.
const LAURENT_TERMS: usize = 3;

/// Taylor coefficients in ε of `Σ_k dense[k] d^k` at `d = 4 - 2ε` (Horner's
/// scheme; index = power of ε).
fn horner_at_four_dimensions(dense: &[Integer]) -> Vec<Integer> {
    let mut series = Vec::<Integer>::new();
    for coefficient in dense.iter().rev() {
        // series <- series * (4 - 2ε) + coefficient
        let mut next = vec![Integer::zero(); series.len() + 1];
        for (power, value) in series.iter().enumerate() {
            next[power] = &next[power] + &(value * &Integer::from(4i64));
            next[power + 1] = &next[power + 1] - &(value * &Integer::from(2i64));
        }
        next[0] = &next[0] + coefficient;
        series = next;
    }
    series
}

/// Taylor coefficients in ε of an integer polynomial in `d` at `d = 4 - 2ε`
/// as Gaussian integers `(re, im)` (index = power of ε). An optional second
/// variable is the imaginary-unit stand-in of [`gaussian_numbers_as_unit`]
/// and is evaluated at `𝑖`.
fn at_four_dimensions(
    polynomial: &rustred::algebra::CoefficientPolynomial,
) -> Result<Vec<(Integer, Integer)>> {
    if polynomial.nvars() > 2 {
        return Err(
            "Laurent summary expects coefficients in the dimension and the imaginary unit alone"
                .into(),
        );
    }
    // Dense coefficients of d^k, real and imaginary: 𝑖^j cycles through 1, 𝑖, -1, -𝑖.
    let mut dense = [Vec::<Integer>::new(), Vec::<Integer>::new()];
    for term in polynomial {
        let degree = term.exponents.first().copied().unwrap_or(0) as usize;
        let unit = term.exponents.get(1).copied().unwrap_or(0) % 4;
        let part = &mut dense[usize::from(unit % 2)];
        if part.len() <= degree {
            part.resize(degree + 1, Integer::zero());
        }
        part[degree] = if unit < 2 {
            &part[degree] + term.coefficient
        } else {
            &part[degree] - term.coefficient
        };
    }
    let [real, imaginary] = dense.map(|part| horner_at_four_dimensions(&part));
    let at =
        |series: &[Integer], index: usize| series.get(index).cloned().unwrap_or_else(Integer::zero);
    Ok((0..real.len().max(imaginary.len()))
        .map(|index| (at(&real, index), at(&imaginary, index)))
        .collect())
}

/// The first `terms` Laurent coefficients of an exact coefficient at
/// `d = 4 - 2ε` (see [`at_four_dimensions`]): the power of ε of the leading
/// term and the exact Gaussian-rational values `(re, im)`.
fn laurent_at_four_dimensions(
    coefficient: &Coefficient,
    terms: usize,
) -> Result<(i64, Vec<(Rational, Rational)>)> {
    let numerator = at_four_dimensions(&coefficient.numerator)?;
    let denominator = at_four_dimensions(&coefficient.denominator)?;
    let lowest = |series: &[(Integer, Integer)]| {
        series
            .iter()
            .position(|(re, im)| !re.is_zero() || !im.is_zero())
    };
    let (Some(a), Some(b)) = (lowest(&numerator), lowest(&denominator)) else {
        return Err("zero numerator or denominator in a retained coefficient".into());
    };
    let at = |series: &[(Integer, Integer)], index: usize| -> (Rational, Rational) {
        series.get(index).map_or_else(
            || (Rational::zero(), Rational::zero()),
            |(re, im)| (Rational::from(re), Rational::from(im)),
        )
    };
    // Division by the leading denominator term u + 𝑖v: times (u - 𝑖v)/(u² + v²).
    let (u, v) = at(&denominator, b);
    let norm = &(&u * &u) + &(&v * &v);
    let mut values = Vec::<(Rational, Rational)>::with_capacity(terms);
    for order in 0..terms {
        let (mut x, mut y) = at(&numerator, a + order);
        for (index, (p, q)) in values.iter().enumerate() {
            // (x + 𝑖y) -= (p + 𝑖q)(s + 𝑖t)
            let (s, t) = at(&denominator, b + order - index);
            x = &x - &(&(p * &s) - &(q * &t));
            y = &y - &(&(p * &t) + &(q * &s));
        }
        values.push((
            &(&(&x * &u) + &(&y * &v)) / &norm,
            &(&(&y * &u) - &(&x * &v)) / &norm,
        ));
    }
    Ok((a as i64 - b as i64, values))
}

/// `+1.23456789012e2` (twelve significant digits) or `0`.
fn decimal_text(value: &Rational) -> String {
    if value.is_zero() {
        return "0".into();
    }
    let float = Float::from_rational_round(value, 128, RoundingDirection::Nearest);
    let magnitude = if float.is_negative() { -float } else { float };
    let sign = if value.is_negative() { "-" } else { "+" };
    // The precision of a multiple-precision float's `e` format counts significant digits.
    format!("{sign}{magnitude:.12e}")
}

/// `+1.23456789012e2 ε^-5  -3.4e1 ε^-4 …` for real coefficients,
/// `+1.2e0 i ε^-5 …` for purely imaginary ones and `(+1.2e0 -3.4e-1 i) ε^-5 …`
/// otherwise.
fn laurent_text(leading_power: i64, values: &[(Rational, Rational)]) -> String {
    let real = values.iter().all(|(_, im)| im.is_zero());
    let imaginary = values.iter().all(|(re, _)| re.is_zero());
    values
        .iter()
        .enumerate()
        .map(|(offset, (re, im))| {
            let power = leading_power + offset as i64;
            match (real, imaginary) {
                (true, _) => format!("{} ε^{power}", decimal_text(re)),
                (false, true) if im.is_zero() => format!("0 ε^{power}"),
                (false, true) => format!("{} i ε^{power}", decimal_text(im)),
                (false, false) => {
                    format!("({} {} i) ε^{power}", decimal_text(re), decimal_text(im))
                }
            }
        })
        .collect::<Vec<_>>()
        .join("  ")
}

/// `scalar` with every non-real number `a + b𝑖` written as `a + b unit`, also
/// inside powers: rustred factors denominators over Q(𝑖), e.g.
/// `(ε-5/3+𝑖/3)^-1 (ε-5/3-𝑖/3)^-1`. The result converts to a rational function
/// over Q in `d` and `unit`; [`at_four_dimensions`] evaluates `unit` at `𝑖`.
fn gaussian_numbers_as_unit(scalar: &Atom, unit: &Atom) -> Atom {
    scalar.replace_map(|part, _, output| {
        if let AtomView::Num(number) = part
            && let symbolica::coefficient::Coefficient::Complex(value) =
                number.get_coeff_view().to_owned()
            && !value.im.is_zero()
        {
            **output = &Atom::num(value.re) + &(&Atom::num(value.im) * unit);
        }
    })
}

/// Short names `M1, M2, …` for formal masters, ordered by loop count (the
/// number of family powers is L(L+1)/2), with a printable legend.
fn master_names(symbols: &BTreeSet<Atom>) -> (BTreeMap<Atom, String>, Vec<serde_json::Value>) {
    let describe = |master: &Atom| -> (usize, String, String) {
        let AtomView::Fun(function) = master.as_view() else {
            return (0, master.to_canonical_string(), String::new());
        };
        let arguments = function.iter().collect::<Vec<_>>();
        let family = arguments
            .first()
            .map(|family| {
                let text = family.to_canonical_string();
                let name = text.rsplit("::").next().unwrap_or(&text).to_owned();
                match name.strip_prefix("rustred_family_") {
                    Some(hash) => hash.chars().take(8).collect(),
                    None => name,
                }
            })
            .unwrap_or_default();
        let powers = arguments
            .iter()
            .skip(1)
            .map(|power| power.to_canonical_string())
            .collect::<Vec<_>>();
        let count = powers.len();
        let loops = (1..=8)
            .find(|loops| loops * (loops + 1) / 2 == count)
            .unwrap_or(0);
        (loops, family, powers.join(","))
    };
    let mut ordered = symbols
        .iter()
        .map(|master| (describe(master), master.clone()))
        .collect::<Vec<_>>();
    ordered.sort();
    let mut names = BTreeMap::new();
    let mut legend = Vec::new();
    for (index, ((loops, family, powers), master)) in ordered.into_iter().enumerate() {
        let name = format!("M{}", index + 1);
        legend.push(
            serde_json::json!({"name":name,"loops":loops,"family":family,"powers":powers,
            "master":master.to_canonical_string()}),
        );
        names.insert(master, name);
    }
    (names, legend)
}

/// `M1^2*M7`: a master product in the short names of [`master_names`].
fn product_label(product: AtomView<'_>, names: &BTreeMap<Atom, String>) -> String {
    let factor = |view: AtomView<'_>| -> String {
        if let AtomView::Pow(power) = view
            && let Some(name) = names.get(&power.get_base().to_owned())
        {
            return format!("{name}^{}", power.get_exp().to_canonical_string());
        }
        names
            .get(&view.to_owned())
            .cloned()
            .unwrap_or_else(|| view.to_canonical_string())
    };
    match product {
        AtomView::Mul(factors) => factors.iter().map(factor).collect::<Vec<_>>().join("*"),
        view if view == Atom::one().as_view() => "1".into(),
        view => factor(view),
    }
}

fn masters(value: &Atom) -> BTreeSet<Atom> {
    let mut result = BTreeSet::new();
    value.visitor(&mut |part| {
        if let AtomView::Fun(function) = part
            && function.get_symbol() == vakint_symbol!("RustRedMaster")
        {
            result.insert(part.to_owned());
            return false;
        }
        true
    });
    result
}

fn product_is_polynomial(product: AtomView<'_>, symbols: &BTreeSet<Atom>) -> bool {
    if symbols.contains(&product.to_owned()) || product == Atom::one().as_view() {
        return true;
    }
    match product {
        AtomView::Mul(factors) => factors
            .iter()
            .all(|factor| product_is_polynomial(factor, symbols)),
        AtomView::Pow(power) => {
            symbols.contains(&power.get_base().to_owned())
                && i64::try_from(power.get_exp()).is_ok_and(|exponent| exponent >= 0)
        }
        _ => false,
    }
}

fn unit_point(value: &Atom, report: &serde_json::Value) -> Result<Atom> {
    let model = super::super::Model::from_json(include_str!("../../qcd.json"))?;
    let scale = report["mu_r_sq_symbol"]
        .as_str()
        .ok_or("target report has no renormalization-scale symbol")?;
    Ok(model
        .expand_couplings(value)
        .replace(Atom::var(vakint_symbol!("mUV")))
        .with(Atom::one())
        .replace(Atom::var(vakint_symbol!(scale)))
        .with(Atom::one())
        .replace(Atom::var(symbolica::symbol!("UFO::G")))
        .with(Atom::one()))
}

/// Exact Q(d) conversion of one expanded scalar. Terms sharing their negative
/// powers are summed as one polynomial numerator before the group is divided
/// by its denominator, instead of normalizing every term as a separate rational
/// function against the same large denominator. Both compute the same rational
/// function; any other shape, refusal or variable outside `variables` takes the
/// direct conversion, so results, errors and diagnostics are unchanged.
fn rational_coefficient(
    scalar: AtomView<'_>,
    variables: &Arc<Vec<PolyVariable>>,
) -> std::result::Result<Coefficient, PolynomialConversionError> {
    let direct = || scalar.try_to_rational_polynomial(&Q, &Z, Some(variables.clone()));
    let AtomView::Add(sum) = scalar else {
        return direct();
    };
    let mut groups = BTreeMap::<Atom, Vec<Atom>>::new();
    for term in sum.iter() {
        let factors = match term {
            AtomView::Mul(product) => product.iter().collect::<Vec<_>>(),
            view => vec![view],
        };
        let (denominator, numerator): (Vec<_>, Vec<_>) = factors.into_iter().partition(|factor| {
            let AtomView::Pow(power) = factor else {
                return false;
            };
            i64::try_from(power.get_exp()).is_ok_and(|exponent| exponent < 0)
        });
        let denominator = if denominator.is_empty() {
            Atom::num(1)
        } else {
            Atom::mul_many(denominator)
        };
        let numerator = if numerator.is_empty() {
            Atom::num(1)
        } else {
            Atom::mul_many(numerator)
        };
        groups.entry(denominator).or_default().push(numerator);
    }
    let mut result: Option<Coefficient> = None;
    for (denominator, numerators) in groups {
        let convert = |part: Atom| -> Option<Coefficient> {
            part.as_view()
                .try_to_rational_polynomial(&Q, &Z, Some(variables.clone()))
                .ok()
        };
        let (Some(numerator), Some(denominator)) =
            (convert(Atom::add_many(numerators)), convert(denominator))
        else {
            return direct();
        };
        if numerator.get_variables() != variables || denominator.get_variables() != variables {
            return direct();
        }
        let part = &numerator * &denominator;
        result = Some(match result {
            Some(sum) => &sum + &part,
            None => part,
        });
    }
    result.map_or_else(direct, Ok)
}

/// Collect only in the formal masters; rational normalization is scalar-local.
/// Products are converted independently on up to `threads` workers and merged
/// in collection order, exactly as one thread would.
fn coefficients(
    value: &Atom,
    symbols: &BTreeSet<Atom>,
    common: &Atom,
    epsilon: &Atom,
    context: &CoefficientContext,
    threads: usize,
) -> Result<BTreeMap<Atom, Coefficient>> {
    let variables = context.one().get_variables().clone();
    let dimension = context.parameter("d").unwrap().to_expression();
    let physical_epsilon = (Atom::num(4) - &dimension) / Atom::num(2);
    let parameters = symbols.iter().cloned().collect::<Vec<_>>();
    let archived_euler = Atom::var(symbolica::transcendental::euler_gamma());
    let current_euler = Atom::var(vakint_symbol!("EulerGamma"));
    let rows = value.coefficient_list::<i32>(&parameters);
    let convert = |(product, scalar): &(Atom, Atom)| -> std::result::Result<Coefficient, String> {
        if !product_is_polynomial(product.as_view(), symbols) || !masters(scalar).is_empty() {
            return Err(format!(
                "retained result is not a polynomial in formal masters at {product}"
            ));
        }
        // The small universal loop factor is never included in the native
        // polynomial variable map. Expand only this one scalar coefficient,
        // allowing its common Euler/phase/pi factors to cancel exactly.
        // The archive uses Symbolica's exact builtin gamma; current MSbar
        // normalization uses Vakint's spelling of the same Euler constant.
        // Bind only this known identity, leaving all foreign symbols visible
        // to strict Q(d) admission below.
        let scalar = (scalar / common)
            .replace(archived_euler.clone())
            .with(current_euler.clone());
        let scalar = canonical_retained(&scalar)
            .map_err(|error| error.to_string())?
            .expand();
        let scalar = scalar
            .replace(epsilon.clone())
            .with(physical_epsilon.clone());
        let coefficient = rational_coefficient(scalar.as_view(), &variables)
            .map_err(|error| format!("coefficient of {product} is not exact Q(d): {error}"))?;
        context
            .validate_with_limits(&coefficient, ExactAlgebraLimits::default())
            .map_err(|error| {
                let describe = |variables: &[symbolica::poly::PolyVariable]| {
                    variables
                        .iter()
                        .map(|variable| match variable {
                            symbolica::poly::PolyVariable::Temporary(_) => format!("{variable:?}"),
                            _ => Atom::from(variable.clone())
                                .to_canonical_string()
                                .chars()
                                .take(240)
                                .collect::<String>(),
                        })
                        .collect::<Vec<_>>()
                };
                let preview = scalar.to_canonical_string().chars().take(240).collect::<String>();
                format!(
                    "coefficient of {product} failed exact Q(d) admission: {error}; expected variables {:?}; numerator variables {:?}; denominator variables {:?}; scalar preview {preview}",
                    describe(&variables),
                    describe(coefficient.numerator.variables()),
                    describe(coefficient.denominator.variables()),
                )
            })?;
        Ok(coefficient)
    };
    let converted = super::map_ordered(threads, &rows, convert);
    let mut result = BTreeMap::new();
    for ((product, _), coefficient) in rows.into_iter().zip(converted) {
        let coefficient = coefficient?;
        if !coefficient.is_zero() {
            if let Some(previous) = result.get_mut(&product) {
                *previous =
                    context.try_add(previous, &coefficient, ExactAlgebraLimits::default())?;
            } else {
                result.insert(product, coefficient);
            }
        }
    }
    result.retain(|_, coefficient| !coefficient.is_zero());
    Ok(result)
}

impl Comparison {
    fn checked(old: &Atom, new: &Atom, report: &serde_json::Value) -> Result<Self> {
        let (converted, factor) = convert_historical_retained(old, report)?;
        let loops = report["loop_count"]
            .as_i64()
            .ok_or("target report has no integer loop count")?;
        let epsilon = Atom::var(vakint_symbol!(
            report["epsilon_symbol"]
                .as_str()
                .ok_or("target report has no epsilon symbol")?
        ));
        let common = vakint_parse!("𝑖/(16*𝜋^2)")?.pow(loops);
        let left = unit_point(&converted, report)?;
        let right = unit_point(new, report)?;
        let symbols = masters(&left)
            .union(&masters(&right))
            .cloned()
            .collect::<BTreeSet<_>>();
        let context = CoefficientContext::try_new(["d"])?;
        let threads = super::fixture_threads();
        let left = coefficients(&left, &symbols, &common, &epsilon, &context, threads)?;
        let right = coefficients(&right, &symbols, &common, &epsilon, &context, threads)?;
        let products = left
            .keys()
            .chain(right.keys())
            .cloned()
            .collect::<BTreeSet<_>>();
        let product_masters = products.iter().flat_map(masters).collect::<BTreeSet<_>>();
        let master_count = product_masters.len();
        let (names, legend) = master_names(&product_masters);
        let mut rows = Vec::with_capacity(products.len());
        for product in products {
            let zero = context.zero();
            let before = left.get(&product).unwrap_or(&zero);
            let after = right.get(&product).unwrap_or(&zero);
            if before != after
                && !context
                    .try_sub(before, after, ExactAlgebraLimits::default())?
                    .is_zero()
            {
                return Err(
                    format!("complete historical coefficient differs for {product}").into(),
                );
            }
            // Display only: a failure is reported in the row, never as an error.
            let laurent = |coefficient: &Coefficient| -> serde_json::Value {
                if coefficient.is_zero() {
                    return serde_json::Value::Null;
                }
                match laurent_at_four_dimensions(coefficient, LAURENT_TERMS) {
                    Ok((power, values)) => serde_json::json!(laurent_text(power, &values)),
                    Err(error) => {
                        serde_json::json!(format!("(Laurent terms unavailable: {error})"))
                    }
                }
            };
            rows.push(serde_json::json!({
                "master_product":product.to_canonical_string(),"exactly_equal":true,
                "label":product_label(product.as_view(), &names),
                "historical_terms":{"numerator":before.numerator.nterms(),"denominator":before.denominator.nterms()},
                "ordered_terms":{"numerator":after.numerator.nterms(),"denominator":after.denominator.nterms()},
                "historical_laurent":laurent(before),"ordered_laurent":laurent(after),
            }));
        }
        Ok(Self {
            converted,
            factor,
            rows,
            master_count,
            masters: legend,
        })
    }
}

/// Readable summary of a completed retained result without a reference: its
/// formal masters (short names and legend) and, per master product, the
/// leading Laurent terms of the exact Q(𝑖)(d) coefficient at `d = 4 - 2ε`
/// (unit point, common `(i/(16π²))^L` removed). Display only; never a gate.
pub fn master_product_summary(completed: &CompletedRun) -> Result<serde_json::Value> {
    let report = &completed.report;
    let loops = report["loop_count"]
        .as_i64()
        .ok_or("target report has no integer loop count")?;
    let epsilon = Atom::var(vakint_symbol!(
        report["epsilon_symbol"]
            .as_str()
            .ok_or("target report has no epsilon symbol")?
    ));
    let common = vakint_parse!("𝑖/(16*𝜋^2)")?.pow(loops);
    let value = unit_point(&completed.value, report)?;
    let symbols = masters(&value);
    // `imaginary_unit` stands for 𝑖 (see `gaussian_numbers_as_unit`).
    let context = CoefficientContext::try_new(["d", "imaginary_unit"])?;
    let variables = context.one().get_variables().clone();
    let dimension = context.parameter("d").unwrap().to_expression();
    let unit = context.parameter("imaginary_unit").unwrap().to_expression();
    let physical_epsilon = (Atom::num(4) - &dimension) / Atom::num(2);
    let archived_euler = Atom::var(symbolica::transcendental::euler_gamma());
    let current_euler = Atom::var(vakint_symbol!("EulerGamma"));
    let mut grouped = BTreeMap::<Atom, Vec<Atom>>::new();
    for (product, scalar) in
        value.coefficient_list::<i32>(&symbols.iter().cloned().collect::<Vec<_>>())
    {
        grouped.entry(product).or_default().push(scalar);
    }
    let grouped = grouped.into_iter().collect::<Vec<_>>();
    // Same scalar preparation as the exact comparison; a coefficient may be
    // complex relative to the common loop factor, so it is converted with 𝑖 as
    // a second variable.
    let convert = |(product, scalars): &(Atom, Vec<Atom>)| -> std::result::Result<_, String> {
        let scalar = (Atom::add_many(scalars.iter().cloned()) / &common)
            .replace(archived_euler.clone())
            .with(current_euler.clone());
        let scalar = canonical_retained(&scalar)
            .map_err(|error| error.to_string())?
            .expand()
            .replace(epsilon.clone())
            .with(physical_epsilon.clone())
            .expand();
        if !masters(&scalar).is_empty() {
            return Err(format!("{product} is not a polynomial in formal masters"));
        }
        let scalar = gaussian_numbers_as_unit(&scalar, &unit);
        if scalar.is_zero() {
            return Ok(None);
        }
        rational_coefficient(scalar.as_view(), &variables)
            .map(|coefficient| (!coefficient.is_zero()).then_some(coefficient))
            .map_err(|error| format!("coefficient of {product} is not in Q(𝑖)(d): {error}"))
    };
    let converted = super::map_ordered(super::fixture_threads(), &grouped, convert);
    let mut products = Vec::new();
    for ((product, _), coefficient) in grouped.into_iter().zip(converted) {
        if let Some(coefficient) = coefficient? {
            products.push((product, coefficient));
        }
    }
    let product_masters = products
        .iter()
        .flat_map(|(product, _)| masters(product))
        .collect::<BTreeSet<_>>();
    let (names, legend) = master_names(&product_masters);
    let mut rows = Vec::with_capacity(products.len());
    for (product, coefficient) in &products {
        let laurent = match laurent_at_four_dimensions(coefficient, LAURENT_TERMS) {
            Ok((power, values)) => laurent_text(power, &values),
            Err(error) => format!("(Laurent terms unavailable: {error})"),
        };
        rows.push(serde_json::json!({
            "master_product":product.to_canonical_string(),
            "label":product_label(product.as_view(), &names),
            "laurent":laurent,
        }));
    }
    Ok(
        serde_json::json!({"master_count":product_masters.len(),"master_product_count":rows.len(),
        "masters":legend,"products":rows}),
    )
}

/// Authenticate the preserved raw result and compare all complete rational
/// coefficients at unit scales. Never reads a cache or archived guard payload.
pub fn compare_saved_zero_momentum(
    payload: &[u8],
    completed: &CompletedRun,
    directory: &Path,
) -> Result<serde_json::Value> {
    if payload.len() != HISTORICAL_BYTES || sha256(payload) != HISTORICAL_SHA256 {
        return Err("historical q0 result failed its accepted byte-length/SHA256 binding".into());
    }
    let report = &completed.report;
    if report["schema"] != "vakint.ordered-rung-comparison.v1"
        || report["loop_count"] != 5
        || report["projection"] != "zero-momentum-metric"
        || report["strategy"] != "Ordered"
        || report["input_sha256"] != INPUT_SHA256
        || report["normalization"] != "MSbar"
        || report["reduction_completed"] != true
        || report["catalogue_closure"]["closed"] != true
        || report["catalogue_closure"]["unknown_count"] != 0
    {
        return Err("historical q0 comparison requires the authenticated completed MSbar ordered zero-momentum cohort".into());
    }
    let actual = completed.value.to_canonical_string();
    if report["retained"]["bytes"].as_u64() != Some(actual.len() as u64)
        || report["retained"]["sha256"].as_str() != Some(sha256(actual.as_bytes()).as_str())
    {
        return Err("ordered result failed its retained payload binding".into());
    }
    let started = Instant::now();
    let historical = Atom::parse(std::str::from_utf8(payload)?, "vakint", Default::default())?;
    let comparison = Comparison::checked(&historical, &completed.value, report)?;
    if comparison.master_count != 55 || comparison.rows.len() != 58 {
        return Err(
            "historical q0 comparison did not preserve its 55 masters and 58 complete products"
                .into(),
        );
    }
    let converted = write_identity(
        directory,
        "historical-converted-msbar.symbolica",
        &comparison.converted.to_canonical_string(),
    )?;
    let factor = comparison.factor.to_canonical_string();
    let certificate = serde_json::json!({
        "schema":"vakint.historical-q0-coefficient-comparison.v1","status":"PASS",
        "historical":{"bytes":HISTORICAL_BYTES,"sha256":HISTORICAL_SHA256,"normalization":"pySecDec"},
        "input_sha256":INPUT_SHA256,"ordered_retained":report["retained"],
        "converted_historical":converted,
        "normalization_conversion":{"api":"VakintSettings::get_integral_normalization_factor_atom","source":"pySecDec","target":"MSbar","loop_count":5,
            "factor":factor,"sha256":sha256(factor.as_bytes()),"before_laurent_truncation":true},
        "dimension_binding":report["dimension_binding"],
        "point":{"mUV":1,"mu_r_squared":1,"UFO_G":1,"coupling_definition":"authenticated QCD model, GC_11=i*G",
            "model_sha256":sha256(include_bytes!("../../qcd.json"))},
        "exact_full_rational_coefficients":true,"laurent_truncation_used":false,"master_values_evaluated":false,
        "master_count":comparison.master_count,"master_product_count":comparison.rows.len(),"products":comparison.rows,
        "masters":comparison.masters,"laurent_scope":"leading Laurent terms of each exact master-product coefficient at d=4-2ε (unit point, common (i/(16π²))^L removed); display only",
        "guards":{"identical_guard_sets_claimed":false,"ordered":report["guards"],
            "historical_sha256":"a266f432acf77ad596282574c8b6485b84eebace226d96ac122535ed308d29c1",
            "scope":"equality of complete rational functions on the intersection of sufficient generic domains; archived guard payload not reread"},
        "comparison_seconds":started.elapsed().as_secs_f64(),"new_reductions":0,
    });
    std::fs::write(
        directory.join("historical-coefficient-comparison.json"),
        serde_json::to_vec_pretty(&certificate)?,
    )?;
    Ok(certificate)
}

#[cfg(test)]
mod tests;
