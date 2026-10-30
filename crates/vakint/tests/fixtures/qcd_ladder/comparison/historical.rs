//! Complete scalar identities against the preserved q0 result, without
//! rationalizing a polynomial whose variables include all formal masters.

use super::{
    Atom, AtomCore, AtomView, CompletedRun, Result, canonical_retained,
    convert_historical_retained, sha256, vakint_parse, vakint_symbol, write_identity,
};
use rustred::algebra::{Coefficient, CoefficientContext, ExactAlgebraLimits};
use std::collections::{BTreeMap, BTreeSet};
use std::path::Path;
use std::time::Instant;
use symbolica::prelude::{Q, Z};

const HISTORICAL_BYTES: usize = 6_683_821;
const HISTORICAL_SHA256: &str = "da909dc0b8f199a9ce46c27368bd2d1c88b277d72e765ab70f9a83680caa1212";
const INPUT_SHA256: &str = "9814c3b49a1b486b9790bbe5458d8bc582c3fb3c299ba6a1a0362db8ba71cb16";

struct Comparison {
    converted: Atom,
    factor: Atom,
    rows: Vec<serde_json::Value>,
    master_count: usize,
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

/// Collect only in the formal masters; rational normalization is scalar-local.
fn coefficients(
    value: &Atom,
    symbols: &BTreeSet<Atom>,
    common: &Atom,
    epsilon: &Atom,
    context: &CoefficientContext,
) -> Result<BTreeMap<Atom, Coefficient>> {
    let variables = context.one().get_variables().clone();
    let dimension = context.parameter("d").unwrap().to_expression();
    let physical_epsilon = (Atom::num(4) - &dimension) / Atom::num(2);
    let parameters = symbols.iter().cloned().collect::<Vec<_>>();
    let archived_euler = Atom::var(symbolica::transcendental::euler_gamma());
    let current_euler = Atom::var(vakint_symbol!("EulerGamma"));
    let mut result = BTreeMap::new();
    for (product, scalar) in value.coefficient_list::<i32>(&parameters) {
        if !product_is_polynomial(product.as_view(), symbols) || !masters(&scalar).is_empty() {
            return Err(format!(
                "retained result is not a polynomial in formal masters at {product}"
            )
            .into());
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
        let scalar = canonical_retained(&scalar)?.expand();
        let scalar = scalar
            .replace(epsilon.clone())
            .with(physical_epsilon.clone());
        let coefficient = scalar
            .as_view()
            .try_to_rational_polynomial(&Q, &Z, Some(variables.clone()))
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
        let left = coefficients(&left, &symbols, &common, &epsilon, &context)?;
        let right = coefficients(&right, &symbols, &common, &epsilon, &context)?;
        let products = left
            .keys()
            .chain(right.keys())
            .cloned()
            .collect::<BTreeSet<_>>();
        let master_count = products
            .iter()
            .flat_map(masters)
            .collect::<BTreeSet<_>>()
            .len();
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
            rows.push(serde_json::json!({
                "master_product":product.to_canonical_string(),"exactly_equal":true,
                "historical_terms":{"numerator":before.numerator.nterms(),"denominator":before.denominator.nterms()},
                "ordered_terms":{"numerator":after.numerator.nterms(),"denominator":after.denominator.nterms()},
            }));
        }
        Ok(Self {
            converted,
            factor,
            rows,
            master_count,
        })
    }
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
