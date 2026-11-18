//! Instrumentation (opt-in, no behaviour change): per-call counts of
//! [`super::integrated::to_vakint_integrand`] and of the outer Taylor product
//! `Σ_{i+j<=degree} P_i D_j` that it expands.
//!
//! The counts are collected only when the `vakint::stage_timing` INFO tracing
//! target is enabled, `PINT6_CLOSURE_DUMP` is set or the Vakint lowering census
//! (`VAKINT_LOWERING_CENSUS`) is on. Each `to_vakint_integrand` call then emits
//! one INFO event on the `vakint::stage_timing` target (message
//! `to_vakint counts`, field `counts` = one JSON object), adds the same keys to
//! the census `to_vakint` record, and hands the object to every PINT6 closure
//! record written for that call (line `to_vakint_counts <json>`).
//!
//! The outer Taylor site ([`taylor_shape`]) records `P`, the product of the
//! `den`-free factors of the compensated remainder, and `D`, the product of
//! its `den` factors, with their coefficients `P_i`, `D_j` of the Taylor
//! parameter; the next `to_vakint_integrand` call on the same thread reports
//! them under `taylor`.

use std::{
    cell::RefCell,
    collections::{HashMap, HashSet},
    sync::atomic::{AtomicUsize, Ordering},
};

use symbolica::atom::{Atom, AtomCore, AtomView, Symbol};

use crate::utils::GS;

thread_local! {
    static TAYLOR_SHAPE: RefCell<Option<serde_json::Value>> = const { RefCell::new(None) };
    static LAST_COUNTS: RefCell<Option<String>> = const { RefCell::new(None) };
}

static CALLS: AtomicUsize = AtomicUsize::new(0);

/// Whether the counts are collected (see the module documentation).
pub(crate) fn enabled() -> bool {
    tracing::enabled!(target: "vakint::stage_timing", tracing::Level::INFO)
        || crate::uv::rstar::pint6_dump::enabled()
        || vakint::census::enabled()
        || crate::uv::rstar::top_mass_split_enabled()
}

/// The sequence number of a counted `to_vakint_integrand` call (per process).
pub(crate) fn next_call() -> usize {
    CALLS.fetch_add(1, Ordering::Relaxed)
}

/// Keep the outer Taylor shape for the next counted conversion on this thread.
pub(crate) fn set_taylor_shape(shape: serde_json::Value) {
    TAYLOR_SHAPE.with(|cell| *cell.borrow_mut() = Some(shape));
}

/// The outer Taylor shape left by [`set_taylor_shape`], if any (cleared).
pub(crate) fn take_taylor_shape() -> Option<serde_json::Value> {
    TAYLOR_SHAPE.with(|cell| cell.borrow_mut().take())
}

/// Emit the counts of one conversion (tracing event) and keep them for the
/// closure records of this thread's following massive closures.
pub(crate) fn publish(counts: &serde_json::Value) {
    let line = counts.to_string();
    tracing::info!(
        target: "vakint::stage_timing",
        counts = %line,
        "to_vakint counts"
    );
    LAST_COUNTS.with(|cell| *cell.borrow_mut() = Some(line));
}

/// The counts of this thread's last counted conversion (one JSON line).
pub(crate) fn last_counts() -> Option<String> {
    LAST_COUNTS.with(|cell| cell.borrow().clone())
}

/// Forget this thread's last counted conversion (a caller that reads
/// [`last_counts`] for one conversion of its own).
pub(crate) fn clear_last_counts() {
    LAST_COUNTS.with(|cell| *cell.borrow_mut() = None);
}

fn unit(degree: usize) -> Vec<f64> {
    let mut orders = vec![0.0; degree + 1];
    orders[0] = 1.0;
    orders
}

fn convolve(left: &[f64], right: &[f64]) -> Vec<f64> {
    let mut product = vec![0.0; left.len()];
    for (i, a) in left.iter().enumerate() {
        if *a == 0.0 {
            continue;
        }
        for (j, b) in right.iter().enumerate().take(left.len() - i) {
            product[i + j] += a * b;
        }
    }
    product
}

/// Expanded summand counts of the coefficients of `parameter^0 ..
/// parameter^degree` in `view`, before like terms merge: sums add, products
/// convolve, and anything else that depends on the parameter (a function of
/// it, a negative power) is expanded as a series once.
fn orders(view: AtomView<'_>, parameter: Symbol, degree: usize) -> Vec<f64> {
    if !view.contains_symbol(parameter) {
        return unit(degree);
    }
    match view {
        AtomView::Var(_) => {
            let mut orders = vec![0.0; degree + 1];
            if degree >= 1 {
                orders[1] = 1.0;
            }
            orders
        }
        AtomView::Add(sum) => sum.iter().fold(vec![0.0; degree + 1], |mut total, term| {
            for (sum, value) in total.iter_mut().zip(orders(term, parameter, degree)) {
                *sum += value;
            }
            total
        }),
        AtomView::Mul(product) => product.iter().fold(unit(degree), |total, factor| {
            convolve(&total, &orders(factor, parameter, degree))
        }),
        AtomView::Pow(power) => {
            let (base, exponent) = power.get_base_exp();
            match i64::try_from(exponent) {
                Ok(n) if n > 0 => {
                    let base = orders(base, parameter, degree);
                    (0..n).fold(unit(degree), |total, _| convolve(&total, &base))
                }
                _ => series_orders(view, parameter, degree).0,
            }
        }
        _ => series_orders(view, parameter, degree).0,
    }
}

/// The coefficients of `parameter^0 .. parameter^degree` of `view` as a
/// series, with their expanded summand counts.
fn series_orders(view: AtomView<'_>, parameter: Symbol, degree: usize) -> (Vec<f64>, Vec<Atom>) {
    let mut counts = vec![0.0; degree + 1];
    let mut coefficients = vec![Atom::Zero; degree + 1];
    if let Ok(series) = AtomCore::series(&view, parameter, 0, degree) {
        for (power, coefficient) in series.terms() {
            let index = power
                .is_integer()
                .then(|| power.numerator_ref().to_i64())
                .flatten()
                .and_then(|power| usize::try_from(power).ok())
                .filter(|power| *power <= degree);
            if let Some(index) = index {
                counts[index] = super::streamed::expansion_estimate(coefficient.as_view()).0;
                coefficients[index] = coefficient.clone();
            }
        }
    }
    (counts, coefficients)
}

fn is_denominator(factor: AtomView<'_>) -> bool {
    match factor {
        AtomView::Fun(function) => function.get_symbol() == GS.den,
        AtomView::Pow(power) => matches!(
            power.get_base_exp().0,
            AtomView::Fun(function) if function.get_symbol() == GS.den
        ),
        _ => false,
    }
}

/// The shape of the outer Taylor product of the compensated remainder
/// `compensated` (before truncation and before the parameter is set to one).
/// Each summand (the remainder itself if it is not a sum) is split into `D`,
/// its factors that are `den` powers, and `P`, its other factors (a factor
/// that only contains a `den`, e.g. a sum of products with `den`s, stays in
/// `P` and is counted without expanding it: `mixed_factors`). Only `D`, a
/// product of a few `den` powers, is expanded:
///
/// - `p_summands`, `p_bytes`: summands and bytes of `P` (of the whole
///   remainder if it is a sum of several products);
/// - `p_orders[i]`: expanded summands of `P_i` (coefficient of `parameter^i`);
/// - `d_orders[j]`: summands of the expanded `D_j` (summed over the distinct
///   `D` of a non-factorized remainder);
/// - `den_patterns`: distinct products of `den` powers in the expanded `D_j`;
/// - `pair_by_order[n]`: `Σ_{i+j=n} |P_i| |D_j|` (the share of each Taylor order);
/// - `pair_estimate`: `Σ_{i+j<=degree} |P_i| |D_j|`, the summands of the
///   expanded Taylor polynomial before like terms merge.
pub(crate) fn taylor_shape(
    compensated: AtomView<'_>,
    parameter: Symbol,
    degree: usize,
) -> serde_json::Value {
    let summands = match compensated {
        AtomView::Add(sum) => sum.iter().collect::<Vec<_>>(),
        _ => vec![compensated],
    };
    let mut p_orders = vec![0.0; degree + 1];
    let mut d_orders = vec![0.0; degree + 1];
    let mut pair_estimate = 0.0;
    let mut pair_by_order = vec![0.0; degree + 1];
    let mut patterns = HashSet::<Vec<Vec<u8>>>::new();
    let mut distinct_d = HashMap::<Vec<u8>, Vec<f64>>::new();
    let (mut p_summands, mut p_bytes, mut d_factors) = (0_usize, 0_usize, 0_usize);
    let mut mixed_factors = 0_usize;
    for summand in &summands {
        let factors = match summand {
            AtomView::Mul(product) => product.iter().collect::<Vec<_>>(),
            other => vec![*other],
        };
        let (dens, rest): (Vec<_>, Vec<_>) = factors
            .into_iter()
            .partition(|factor| is_denominator(*factor));
        d_factors = d_factors.max(dens.len());
        mixed_factors += rest
            .iter()
            .filter(|factor| factor.contains_symbol(GS.den))
            .count();
        let p = rest.iter().fold(unit(degree), |total, factor| {
            convolve(&total, &orders(*factor, parameter, degree))
        });
        if summands.len() == 1 {
            p_summands = rest
                .iter()
                .map(|factor| vakint::census::summands(*factor))
                .max()
                .unwrap_or(1);
            p_bytes = rest.iter().map(|factor| factor.get_byte_size()).sum();
        }
        if dens.is_empty() {
            let d = unit(degree);
            for (i, value) in p.iter().enumerate() {
                p_orders[i] += value;
            }
            pair_estimate += p.iter().sum::<f64>();
            for (i, value) in p.iter().enumerate() {
                pair_by_order[i] += value;
            }
            if summands.len() == 1 {
                d_orders = d;
            }
            continue;
        }
        let d_atom = Atom::mul_many(dens.iter().copied());
        let d = distinct_d
            .entry(d_atom.as_view().get_data().to_vec())
            .or_insert_with(|| {
                let (_, coefficients) = series_orders(d_atom.as_view(), parameter, degree);
                coefficients
                    .iter()
                    .map(|coefficient| {
                        let expanded = coefficient.expand();
                        let terms = match expanded.as_view() {
                            AtomView::Add(sum) => sum.iter().collect::<Vec<_>>(),
                            view if view.is_zero() => Vec::new(),
                            view => vec![view],
                        };
                        for term in &terms {
                            let mut pattern = match term {
                                AtomView::Mul(product) => product
                                    .iter()
                                    .filter(|factor| is_denominator(*factor))
                                    .map(|factor| factor.get_data().to_vec())
                                    .collect::<Vec<_>>(),
                                factor if is_denominator(*factor) => {
                                    vec![factor.get_data().to_vec()]
                                }
                                _ => Vec::new(),
                            };
                            pattern.sort_unstable();
                            patterns.insert(pattern);
                        }
                        terms.len() as f64
                    })
                    .collect()
            })
            .clone();
        for (i, value) in p.iter().enumerate() {
            p_orders[i] += value;
            for (j, d_value) in d.iter().enumerate().take(degree + 1 - i) {
                pair_estimate += value * d_value;
                pair_by_order[i + j] += value * d_value;
            }
        }
        if summands.len() == 1 || distinct_d.len() == 1 {
            d_orders.clone_from(&d);
        }
    }
    if summands.len() != 1 {
        p_summands = summands.len();
        p_bytes = compensated.get_byte_size();
        if distinct_d.len() != 1 {
            d_orders = distinct_d
                .values()
                .fold(vec![0.0; degree + 1], |mut total, d| {
                    for (sum, value) in total.iter_mut().zip(d) {
                        *sum += value;
                    }
                    total
                });
        }
    }
    serde_json::json!({
        "degree": degree,
        "factorized": summands.len() == 1,
        "remainder_summands": summands.len(),
        "p_summands": p_summands,
        "p_bytes": p_bytes,
        "p_orders": p_orders,
        "d_factors": d_factors,
        "mixed_factors": mixed_factors,
        "distinct_d": distinct_d.len(),
        "d_orders": d_orders,
        "den_patterns": patterns.len(),
        "pair_estimate": pair_estimate,
        "pair_by_order": pair_by_order,
    })
}
