//! Joint p reduction (`VAKINT_P_JOINT`, default `on`): every channel of a
//! batch whose non-rational coefficient is `Σ_M c_M(ε)·M`, with `c_M` a
//! rational function of `ε` and `M` free of `ε` (e.g. a power of the
//! external `dot(p,p)`), is folded into ONE request per (loop count,
//! external head, scale power, `M`): `c_M` moves into the input coefficients
//! (`ε = (4-d)/2` in the library's `Q(d)`), identical integrals of different
//! channels merge on entry, and the value is `M` times the reduction. The
//! reduction is linear, so the joint request returns
//! `Σ_channel c_M · reduction(channel)` term by term.
//!
//! physics_L4 (one graph, 595 channels, 56,098 inputs over 4,208 distinct
//! keys) needs 5 requests instead of 595, and its Retain output collapses
//! from ~2 MB to ~21 KB (one value per scale power, as Forcer's): an equal
//! but differently arranged expression than one value per channel.
//!
//! `off`: one request per channel (the arrangement before the joint path).
//! `check`: both, an error on any difference of the exact `Q(d)` terms, the
//! per-channel values returned; the profiler (`VAKINT_P_PHASES=2`) prints
//! the work of both.

use std::collections::BTreeMap;
use std::time::Instant;

use rustred::algebra::{Coefficient, CoefficientContext, ExactAlgebraLimits};
use rustred::family::IntegralKey;
use rustred::propagator::PresentationId;
use rustred::propagator::library::{PLibrary, PReduction};
use symbolica::atom::{Atom, AtomCore, AtomView, Symbol};
use symbolica::coefficient::CoefficientView;

use super::{
    Channel, Channels, PRequestOptions, RustRedEvaluationError, coefficient, invalid,
    p_reduction_request,
};
use crate::p_integral::phases;
use rustred::reduction::OrderedStepCache;
use std::sync::Arc;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) enum JointMode {
    Off,
    On,
    Check,
}

pub(super) fn mode() -> JointMode {
    static MODE: std::sync::OnceLock<JointMode> = std::sync::OnceLock::new();
    *MODE.get_or_init(|| match std::env::var("VAKINT_P_JOINT").as_deref() {
        Ok("off" | "0") => JointMode::Off,
        Ok("check") => JointMode::Check,
        _ => JointMode::On,
    })
}

/// A rational function of `epsilon` alone as a coefficient of `context`
/// (`ε = (4 - d)/2`), or `None` for anything else.
fn epsilon_coefficient(
    context: &CoefficientContext,
    atom: AtomView,
    epsilon: Symbol,
    limits: ExactAlgebraLimits,
) -> Option<Coefficient> {
    match atom {
        AtomView::Num(number) => match number.get_coeff_view() {
            CoefficientView::Natural(numerator, denominator, 0, _) => {
                let numerator = context.integer(numerator);
                let denominator = context.integer(denominator);
                context.try_div(&numerator, &denominator, limits).ok()
            }
            _ => None,
        },
        AtomView::Var(variable) if variable.get_symbol() == epsilon => {
            let d = context.parameter("d")?;
            let four_minus_d = context.try_sub(&context.integer(4), &d, limits).ok()?;
            context
                .try_div(&four_minus_d, &context.integer(2), limits)
                .ok()
        }
        AtomView::Add(sum) => {
            let mut total = context.zero();
            for term in sum.iter() {
                let term = epsilon_coefficient(context, term, epsilon, limits)?;
                total = context.try_add(&total, &term, limits).ok()?;
            }
            Some(total)
        }
        AtomView::Mul(product) => {
            let mut total = context.one();
            for factor in product.iter() {
                let factor = epsilon_coefficient(context, factor, epsilon, limits)?;
                total = context.try_mul(&total, &factor, limits).ok()?;
            }
            Some(total)
        }
        AtomView::Pow(power) => {
            let (base, exponent) = power.get_base_exp();
            let AtomView::Num(exponent) = exponent else {
                return None;
            };
            let CoefficientView::Natural(exponent, 1, 0, _) = exponent.get_coeff_view() else {
                return None;
            };
            let base = epsilon_coefficient(context, base, epsilon, limits)?;
            let mut total = context.one();
            for _ in 0..exponent.unsigned_abs() {
                total = context.try_mul(&total, &base, limits).ok()?;
            }
            if exponent < 0 {
                total = context.try_div(&context.one(), &total, limits).ok()?;
            }
            Some(total)
        }
        _ => None,
    }
}

/// A channel coefficient as `Σ_M c_M(ε) · M` with `M` a product free of
/// `ε` (e.g. a power of the external `dot(p,p)`), or `None` when a factor
/// holding `ε` is not a rational function of it.
fn decompose(
    context: &CoefficientContext,
    coefficient: AtomView,
    epsilon: Symbol,
    limits: ExactAlgebraLimits,
) -> Option<Vec<(Atom, Coefficient)>> {
    let expanded = coefficient.expand();
    let terms: Vec<AtomView> = match expanded.as_view() {
        AtomView::Add(sum) => sum.iter().collect(),
        view => vec![view],
    };
    let mut parts = Vec::<(Atom, Coefficient)>::new();
    for term in terms {
        let factors: Vec<AtomView> = match term {
            AtomView::Mul(product) => product.iter().collect(),
            view => vec![view],
        };
        let mut value = context.one();
        let mut rest = Atom::num(1);
        for factor in factors {
            if matches!(factor, AtomView::Num(_)) || factor.contains_symbol(epsilon) {
                let factor = epsilon_coefficient(context, factor, epsilon, limits)?;
                value = context.try_mul(&value, &factor, limits).ok()?;
            } else {
                rest *= factor.to_owned();
            }
        }
        match parts.iter_mut().find(|(monomial, _)| *monomial == rest) {
            Some((_, entry)) => *entry = context.try_add(entry, &value, limits).ok()?,
            None => parts.push((rest, value)),
        }
    }
    parts.retain(|(_, value)| !value.is_zero());
    Some(parts)
}

/// One joint request: the summed inputs of every joinable channel of one
/// (loop count, head, scale power, ε-free monomial `M`), and the channels it
/// replaces; its value is `M` times its reduction.
pub(super) struct JointGroup {
    pub loops: usize,
    pub head: Atom,
    pub scale_power: i64,
    pub monomial: Atom,
    pub inputs: BTreeMap<(PresentationId, Vec<i64>), Coefficient>,
    /// The replaced channels with their coefficient in `Q(d)`.
    pub channels: Vec<(Channel, Coefficient)>,
}

/// The joint groups of `channels` (joinable channels only).
pub(super) fn groups(
    library: &PLibrary,
    channels: &Channels,
    epsilon: Symbol,
    limits: ExactAlgebraLimits,
) -> Result<Vec<JointGroup>, RustRedEvaluationError> {
    let context = library.coefficient_context();
    let mut groups = BTreeMap::<(usize, Atom, i64, Atom), JointGroup>::new();
    for (channel, inputs) in channels {
        let Some(parts) = decompose(context, channel.coefficient.as_view(), epsilon, limits) else {
            continue;
        };
        for (monomial, factor) in parts {
            let group = groups
                .entry((
                    channel.loops,
                    channel.head.clone(),
                    channel.scale_power,
                    monomial.clone(),
                ))
                .or_insert_with(|| JointGroup {
                    loops: channel.loops,
                    head: channel.head.clone(),
                    scale_power: channel.scale_power,
                    monomial,
                    inputs: BTreeMap::new(),
                    channels: Vec::new(),
                });
            for (key, value) in &inputs.inputs {
                if value.is_zero() {
                    continue;
                }
                let value = coefficient(context, value, limits)?;
                let value = context
                    .try_mul(&value, &factor, limits)
                    .map_err(RustRedEvaluationError::from_native)?;
                match group.inputs.get_mut(key) {
                    Some(entry) => {
                        *entry = context
                            .try_add(entry, &value, limits)
                            .map_err(RustRedEvaluationError::from_native)?;
                    }
                    None => {
                        group.inputs.insert(key.clone(), value);
                    }
                }
            }
            group.channels.push((channel.clone(), factor));
        }
    }
    for group in groups.values_mut() {
        group.inputs.retain(|_, value| !value.is_zero());
    }
    Ok(groups.into_values().collect())
}

/// The engine inputs of a joint group.
pub(super) fn group_inputs(
    group: &JointGroup,
) -> Result<Vec<(PresentationId, IntegralKey, Coefficient)>, RustRedEvaluationError> {
    group
        .inputs
        .iter()
        .map(|((presentation, powers), value)| {
            let key = IntegralKey::try_new(powers.clone())
                .map_err(|error| invalid(format!("integral key: {error}")))?;
            Ok((*presentation, key, value.clone()))
        })
        .collect()
}

/// Print the work counters of one request (profiler level 2).
pub(super) fn report(label: &str, inputs: usize, reduction: &PReduction, seconds: f64) {
    if !phases::per_request() {
        return;
    }
    let work = reduction.statistics.work;
    let progress = reduction.statistics.progress;
    eprintln!(
        "vakint-p-work {label} inputs={inputs} input_keys={} rule_applications={} coalescing_additions={} cache_hits={} peak_frontier_keys={} terms={} seconds={seconds:.6}",
        progress.input_keys,
        work.rule_applications(),
        work.coalescing_additions(),
        work.cache_hits(),
        progress.peak_frontier_keys,
        reduction.terms.len(),
    );
}

/// `Σ_channel factor · reduction(channel)` against the joint reduction, term
/// by term (basis and Gamma tag) in `Q(d)`.
pub(super) fn compare(
    context: &CoefficientContext,
    limits: ExactAlgebraLimits,
    separate: &[(Coefficient, PReduction)],
    joint: &PReduction,
) -> Result<(), RustRedEvaluationError> {
    let mut expected = Vec::<(String, Coefficient)>::new();
    let mut add = |key: String, value: Coefficient| -> Result<(), RustRedEvaluationError> {
        match expected.iter_mut().find(|(existing, _)| *existing == key) {
            Some((_, entry)) => {
                *entry = context
                    .try_add(entry, &value, limits)
                    .map_err(RustRedEvaluationError::from_native)?;
            }
            None => expected.push((key, value)),
        }
        Ok(())
    };
    for (factor, reduction) in separate {
        for term in &reduction.terms {
            let value = context
                .try_mul(factor, &term.coefficient, limits)
                .map_err(RustRedEvaluationError::from_native)?;
            add(format!("{:?}|{:?}", term.basis, term.tag), value)?;
        }
    }
    for term in &joint.terms {
        let value = context
            .try_neg(&term.coefficient, limits)
            .map_err(RustRedEvaluationError::from_native)?;
        add(format!("{:?}|{:?}", term.basis, term.tag), value)?;
    }
    let differing = expected
        .iter()
        .filter(|(_, value)| !value.is_zero())
        .count();
    if differing != 0 {
        return Err(invalid(format!(
            "VAKINT_P_JOINT=check: the joint reduction differs from the sum of the channel \
             reductions in {differing} of {} terms",
            expected.len()
        )));
    }
    Ok(())
}

/// Reduce every joint group and, in check mode, its channels one by one;
/// compare; return the joint reductions in group order.
pub(super) fn reduce_groups(
    library: &PLibrary,
    options: &PRequestOptions,
    step_cache: &Option<Arc<OrderedStepCache>>,
    groups: &[JointGroup],
    channels: &Channels,
    check: bool,
    limits: ExactAlgebraLimits,
) -> Result<Vec<PReduction>, RustRedEvaluationError> {
    let context = library.coefficient_context();
    let mut joint = Vec::with_capacity(groups.len());
    let (mut joint_seconds, mut separate_seconds) = (0.0, 0.0);
    for group in groups {
        let inputs = group_inputs(group)?;
        let count = inputs.len();
        let started = Instant::now();
        let reduction = {
            let _math = super::MathSpan::new();
            library
                .reduce(p_reduction_request(options, inputs, step_cache.clone()))
                .map_err(RustRedEvaluationError::from_native)?
        };
        let seconds = started.elapsed().as_secs_f64();
        joint_seconds += seconds;
        report(
            &format!("joint scale_power={}", group.scale_power),
            count,
            &reduction,
            seconds,
        );
        if check {
            let mut separate = Vec::with_capacity(group.channels.len());
            for (channel, factor) in &group.channels {
                let inputs =
                    super::channel_inputs(library, channels[channel].inputs.clone(), limits)?;
                if inputs.is_empty() {
                    continue;
                }
                let count = inputs.len();
                let started = Instant::now();
                let reduction = library
                    .reduce(p_reduction_request(options, inputs, step_cache.clone()))
                    .map_err(RustRedEvaluationError::from_native)?;
                let seconds = started.elapsed().as_secs_f64();
                separate_seconds += seconds;
                report("separate", count, &reduction, seconds);
                separate.push((factor.clone(), reduction));
            }
            compare(context, limits, &separate, &reduction)?;
        }
        joint.push(reduction);
    }
    if phases::enabled() {
        eprintln!(
            "vakint-phase p.joint_reduce {joint_seconds:.6} count={} channels={}",
            groups.len(),
            groups
                .iter()
                .map(|group| group.channels.len())
                .sum::<usize>()
        );
        if check {
            eprintln!("vakint-phase p.joint_check_separate_reduce {separate_seconds:.6}");
        }
    }
    Ok(joint)
}
