//! The RustRed p backend: massless propagator-type integrals reduced by the
//! engine's RustRed p library program (M2 plan, sections 3.3-3.5, M2-V3;
//! `docs/p_integrals.md`, sections 5-6).
//!
//! A batch of matched p contributions is split into the homogeneous parts of
//! the shared numerator split (`p_integral::parts`, the one Forcer renders),
//! grouped initially by (loop count, external head `p(i)`, scale power `w`,
//! non-rational coefficient), then folded by spectator monomial: rational
//! functions of ε enter the reduction as rational functions of d = 4-2ε.
//! Inputs are lowered onto the slots of the family's
//! top presentation (lines, then ISPs; a numerator scalar product is the
//! presentation's own affine form in its slots, Minkowski with `q·q = -1`) and
//! reduced with one [`PReductionRequest`] per folded channel onto the masters:
//! `PReduction` is `(εG(1,1))^-L Σ c·I_M(x)` at `Q² = 1` (convention C3), which
//! is exactly vakint's `reduced = (-1)^(Σa+r)·C_G`. The Wick sign of a part
//! (`PPart::wick_sign`) belongs to the Forcer path only and is never applied
//! here (critique-phys SF3). The terms are rendered as
//! `Σ coefficient(d = 4-2ε) · Π PGamma(m)^e · (1 | PMaster(name))` (convention
//! C4) and finalized by the shared p-finalizer with the backend name `RustRed`
//! (normalization, `(Q²)^w`, Laurent expansion and the `Oep` guard).
//!
//! Entry points (crate-visible): [`supports`] and [`has_p_data`] for the
//! dispatch, [`evaluate_ordered_batch`] for a batch of matched contributions,
//! and [`reduce_requests`], the request level below the finalizer (backend
//! requests in, reduced values per scale power out), for crate tests that
//! compare the RustRed p path with Forcer's frozen values per scale power
//! (M2-V5 Tier A).

use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::sync::Arc;

use rustred::algebra::{Coefficient, CoefficientContext, ExactAlgebraLimits};
use rustred::family::IntegralKey;
use rustred::propagator::PresentationId;
use rustred::propagator::error::PError;
use rustred::propagator::library::{PBasis, PReduction, PReductionRequest, PUncoveredPolicy};
use rustred::propagator::presentation::{AffineSlotForm, Presentation};
use rustred::reduction::{
    OrderedSchedulerOptions, OrderedStepCache, OrderedWeightArithmetic, ReductionLimits,
};
use symbolica::atom::{Atom, AtomCore, Symbol};
use symbolica::domains::integer::Integer;
use symbolica::domains::rational::Rational;
use symbolica::function;

use super::p_data::{PDataPresence, PLibraryHandle, PTop};
use super::{MatchedContribution, MathSpan, RustRedEvaluationError, RustRedMasterMode};
use crate::p_integral::parts::{self, Fraction, PDot, PMonomial};
use crate::p_integral::{
    P_FAMILIES, P_INTEGRAL_LOOP_COUNTS, PropagatorRequest, ReducedPropagatorValue, external_head,
    finalize, masters, p_gamma_symbol, p_master_symbol,
};
use crate::utils::simplify_real;
use crate::utils::vakint_macros::vk_symbol;
use crate::{EvaluationMethod, Topology, Vakint, VakintError, VakintSettings};

fn invalid(detail: impl Into<String>) -> RustRedEvaluationError {
    RustRedEvaluationError::InvalidMatchedFamily {
        detail: format!("p-integral: {}", detail.into()),
    }
}

/// Whether the RustRed backend reduces this topology structurally: a matched
/// p-integral of one to four loops. Whether a library is loaded for its loop
/// count is [`has_p_data`].
pub(crate) fn supports(topology: &Topology) -> bool {
    let Topology::Propagator(integral) = topology else {
        return false;
    };
    P_INTEGRAL_LOOP_COUNTS.contains(&integral.n_loops)
        && integral
            .applicable_evaluation_methods
            .0
            .iter()
            .any(|method| matches!(method, EvaluationMethod::RustRed(_)))
}

/// Whether the engine's p library covers `loops`. A `p/` area of the loaded
/// external data that cannot be resolved counts as present: its p-integrals
/// then fail with that error instead of silently going to another backend.
pub(crate) fn has_p_data(vakint: &Vakint, loops: usize) -> bool {
    match vakint.p_library_presence() {
        PDataPresence::None => false,
        PDataPresence::Scope(scope) => {
            u8::try_from(loops).is_ok_and(|loops| scope.contains(&loops))
        }
        PDataPresence::Broken => true,
    }
}

/// The engine options of every p request: the Ordered weight arithmetic of
/// the settings (default `modular_polynomial`, the `Z_p(d)` replay), the
/// Ordered scheduler on `n_cores` threads and the native reduction limits.
pub(super) struct PRequestOptions {
    pub arithmetic: OrderedWeightArithmetic,
    pub scheduler: OrderedSchedulerOptions,
    pub limits: ReductionLimits,
}

pub(super) fn request_options(settings: &VakintSettings) -> PRequestOptions {
    PRequestOptions {
        arithmetic: settings.rustred_ordered_weight_arithmetic(),
        scheduler: super::ordered_scheduler(settings),
        limits: settings.rustred_application_options.reduction_limits,
    }
}

/// The request [`reduce_requests`] sends for one channel: `inputs` with the
/// `options` of the settings, the library's step cache, uncovered keys
/// rejected, reduced onto the masters (module documentation).
fn p_reduction_request<'a>(
    options: &PRequestOptions,
    inputs: Vec<(PresentationId, IntegralKey, Coefficient)>,
    step_cache: Option<Arc<OrderedStepCache>>,
) -> PReductionRequest<'a> {
    PReductionRequest {
        inputs,
        arithmetic: options.arithmetic,
        scheduler: options.scheduler,
        step_cache,
        uncovered: PUncoveredPolicy::Reject,
        onto_masters: true,
        limits: options.limits,
        inspector: None,
    }
}

/// One reduction request: the inputs of one loop count, external head, scale
/// power and non-rational coefficient.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
struct Channel {
    loops: usize,
    head: Atom,
    scale_power: i64,
    coefficient: Atom,
}

/// The lowered inputs of one channel and the vakint integrals they come from
/// (`P1L(2,1)`: family and line powers), named when the library cannot reduce
/// them.
#[cfg(test)]
#[derive(Default)]
struct ChannelInputs {
    inputs: BTreeMap<(PresentationId, Vec<i64>), Rational>,
    integrals: BTreeSet<String>,
}

/// At most this many vakint integrals are named in a reduction error.
const NAMED_INTEGRALS: usize = 8;

/// A failed p reduction. A key the library does not cover (`Uncovered`,
/// `MissingOwner`) names the vakint integrals of the channel next to RustRed's
/// owner class and key; any other error maps as for the vacuum programs.
fn reduction_error(error: PError, integrals: &BTreeSet<String>) -> RustRedEvaluationError {
    if !matches!(
        error,
        PError::Uncovered { .. } | PError::MissingOwner { .. }
    ) {
        return RustRedEvaluationError::from_native(error);
    }
    let mut named = integrals
        .iter()
        .take(NAMED_INTEGRALS)
        .cloned()
        .collect::<Vec<_>>();
    if integrals.len() > NAMED_INTEGRALS {
        named.push(format!("{} more", integrals.len() - NAMED_INTEGRALS));
    }
    RustRedEvaluationError::Reduction {
        detail: format!(
            "{error}, reducing the p-integrals {} with the RustRed p library",
            named.join(", ")
        ),
    }
}

/// Exponents of the slots in a polynomial of affine slot forms.
type SlotMonomial = Vec<u32>;

/// `a·b` of a canonical scalar product as momentum rows over `(k_1..k_L, q)`.
fn dot_rows(dot: PDot, loops: usize) -> (Vec<i64>, Vec<i64>) {
    let unit = |index: usize| {
        let mut row = vec![0; loops + 1];
        row[index] = 1;
        row
    };
    match dot {
        PDot::Loops(i, j) => (unit(i - 1), unit(j - 1)),
        PDot::External(i) => (unit(i - 1), unit(loops)),
    }
}

fn fraction(value: Fraction) -> Rational {
    Rational::from((value.0, value.1))
}

/// `polynomial · form` with `form = c + Σ r_s D_s`.
fn multiply(
    polynomial: &BTreeMap<SlotMonomial, Rational>,
    form: &AffineSlotForm,
) -> BTreeMap<SlotMonomial, Rational> {
    let mut product = BTreeMap::<SlotMonomial, Rational>::new();
    let mut add = |monomial: SlotMonomial, value: Rational| {
        let entry = product.entry(monomial).or_insert_with(Rational::zero);
        *entry += &value;
    };
    for (monomial, coefficient) in polynomial {
        if !form.constant.is_zero() {
            add(monomial.clone(), coefficient * &form.constant);
        }
        for (slot, value) in &form.terms {
            let mut shifted = monomial.clone();
            shifted[usize::from(*slot)] += 1;
            add(shifted, coefficient * value);
        }
    }
    product.retain(|_, value| !value.is_zero());
    product
}

/// The inputs `Σ c·I_M(key)` of one homogeneous part on the slots of its top
/// presentation. The line slots carry the positive line powers (a negative
/// power is already a numerator factor of the part), every scalar product of
/// the numerator monomials is the presentation's affine form in its slots, and
/// a slot factor `D_s^m` lowers the slot power by `m` (ISPs at negative
/// powers). The dot forms of a presentation are cached in `forms`.
fn lower_part(
    presentation: &Presentation,
    top: PTop,
    powers: &[i64],
    polynomial: &[(Fraction, PMonomial)],
    forms: &mut BTreeMap<PDot, AffineSlotForm>,
) -> Result<BTreeMap<Vec<i64>, Rational>, RustRedEvaluationError> {
    if powers.len() != top.line_count {
        return Err(invalid(format!(
            "{} line powers for {} lines",
            powers.len(),
            top.line_count
        )));
    }
    let loops = usize::from(presentation.loops);
    let mut lowered = BTreeMap::<Vec<i64>, Rational>::new();
    for (value, monomial) in polynomial {
        let mut product = BTreeMap::from([(vec![0u32; top.arity], fraction(*value))]);
        for &(dot, exponent) in &monomial.0 {
            let form = match forms.entry(dot) {
                std::collections::btree_map::Entry::Occupied(entry) => entry.into_mut(),
                std::collections::btree_map::Entry::Vacant(entry) => {
                    let (a, b) = dot_rows(dot, loops);
                    entry.insert(
                        presentation
                            .express_scalar_product(&a, &b)
                            .map_err(|error| invalid(format!("{dot:?}: {error}")))?,
                    )
                }
            };
            for _ in 0..exponent {
                product = multiply(&product, form);
            }
        }
        for (monomial, coefficient) in product {
            let key = (0..top.arity)
                .map(|slot| {
                    let base = if slot < top.line_count {
                        powers[slot].max(0)
                    } else {
                        0
                    };
                    base - i64::from(monomial[slot])
                })
                .collect::<Vec<_>>();
            let entry = lowered.entry(key).or_insert_with(Rational::zero);
            *entry += &coefficient;
        }
    }
    lowered.retain(|_, value| !value.is_zero());
    Ok(lowered)
}

/// An exact rational as a coefficient of the library's context.
fn coefficient(
    context: &CoefficientContext,
    value: &Rational,
    limits: ExactAlgebraLimits,
) -> Result<Coefficient, RustRedEvaluationError> {
    let integer = |integer: &Integer| -> Result<Coefficient, RustRedEvaluationError> {
        let value = i128::try_from(integer.clone())
            .map_err(|_| invalid(format!("coefficient {value} exceeds 128 bits")))?;
        let magnitude = context.unsigned_integer(value.unsigned_abs());
        if value < 0 {
            context
                .try_neg(&magnitude, limits)
                .map_err(RustRedEvaluationError::from_native)
        } else {
            Ok(magnitude)
        }
    };
    context
        .try_div(
            &integer(value.numerator_ref())?,
            &integer(value.denominator_ref())?,
            limits,
        )
        .map_err(RustRedEvaluationError::from_native)
}

/// The reduction in vakint's basis (convention C4): every coefficient at
/// `d = 4 - 2ε`, the tag as `Π PGamma(m)^e`, `PBasis::Unit` as `1` and
/// `PBasis::Master(name)` as `PMaster(name)` for a tabulated Forcer master;
/// the coefficient of each tag and master is one simplified rational function
/// of `ε` (`together`). A residual terminal or an unknown master fails closed,
/// and so does a nonzero condition that vanishes identically at
/// `d = 4 - 2ε`.
fn render(
    handle: &PLibraryHandle,
    reduction: &PReduction,
    epsilon: Symbol,
) -> Result<Atom, RustRedEvaluationError> {
    Ok(collect_terms(render_terms(handle, reduction, epsilon)?))
}

/// The terms of [`render`]: per tag and master, the coefficient at
/// `d = 4 - 2ε`, not yet simplified.
fn render_terms(
    handle: &PLibraryHandle,
    reduction: &PReduction,
    epsilon: Symbol,
) -> Result<BTreeMap<Atom, Atom>, RustRedEvaluationError> {
    let dimension = Atom::num(4) - Atom::num(2) * Atom::var(epsilon);
    let at_dimension = |coefficient: &Coefficient| {
        coefficient
            .to_expression()
            .replace(handle.dimension.to_pattern())
            .with(dimension.to_pattern())
    };
    for condition in &reduction.nonzero_conditions {
        if at_dimension(condition).together().is_zero() {
            return Err(invalid(format!(
                "the reduction requires {} != 0, which vanishes at d = 4 - 2ε",
                condition.to_expression()
            )));
        }
    }
    let names =
        masters::p_master_names().map_err(|error| RustRedEvaluationError::MasterCapability {
            detail: error.to_string(),
        })?;
    let mut coefficients = BTreeMap::<Atom, Atom>::new();
    for term in &reduction.terms {
        let basis = match &term.basis {
            PBasis::Unit => Atom::num(1),
            PBasis::Master(name) if names.contains(&name.as_str()) => {
                function!(p_master_symbol(), Atom::var(vk_symbol!(name.as_str())))
            }
            PBasis::Master(name) => {
                return Err(RustRedEvaluationError::MasterCapability {
                    detail: format!(
                        "RustRed returned the p-master {name}, which the vendored Forcer \
                         Masters.prc does not tabulate"
                    ),
                });
            }
            PBasis::Residual(terminal) => {
                return Err(RustRedEvaluationError::MasterCapability {
                    detail: format!(
                        "RustRed returned the residual p terminal {} instead of Forcer masters",
                        terminal.0
                    ),
                });
            }
        };
        let mut tag = Atom::num(1);
        for &(m, e) in term.tag.exponents() {
            tag *= p_gamma_symbol()
                .call((i64::from(m),))
                .pow(Atom::num(i64::from(e)));
        }
        *coefficients.entry(tag * basis).or_default() += at_dimension(&term.coefficient);
    }
    Ok(coefficients)
}

/// The lowered inputs of a batch, by channel (the eager table of
/// [`reduce_requests_eager`]).
#[cfg(test)]
type Channels = BTreeMap<Channel, ChannelInputs>;

/// Lower one backend request (`docs/p_integrals.md`, section 5) into
/// `channels`; a scaleless request adds nothing. The eager lowering of
/// [`reduce_requests_eager`].
#[cfg(test)]
fn add_request(
    handle: &PLibraryHandle,
    request: &PropagatorRequest,
    forms: &mut BTreeMap<PresentationId, BTreeMap<PDot, AffineSlotForm>>,
    channels: &mut Channels,
) -> Result<(), VakintError> {
    let family = P_FAMILIES
        .iter()
        .position(|family| family.name == request.family.name)
        .ok_or_else(|| invalid(format!("unknown family {}", request.family.name)))?;
    let top = handle.top(family)?;
    let Some(split) = parts::split_parts(request)? else {
        return Ok(()); // scaleless: an exact zero
    };
    let (head, _) = external_head(&request.external)?;
    let presentation = handle.library.presentation(top.presentation);
    let forms = forms.entry(top.presentation).or_default();
    let integral = integral_name(request);
    for part in split {
        let lowered = lower_part(presentation, top, request.powers, &part.polynomial, forms)?;
        let channel = channels
            .entry(Channel {
                loops: request.family.loops,
                head: head.clone(),
                scale_power: part.scale_power,
                coefficient: part.coefficient,
            })
            .or_default();
        channel.integrals.insert(integral.clone());
        for (key, value) in lowered {
            let entry = channel
                .inputs
                .entry((top.presentation, key))
                .or_insert_with(Rational::zero);
            *entry += &value;
        }
    }
    Ok(())
}

/// `P1L(2,1)`: the family and line powers of a request, as a reduction
/// error names it.
fn integral_name(request: &PropagatorRequest) -> String {
    format!(
        "{}({})",
        request.family.name,
        request
            .powers
            .iter()
            .map(i64::to_string)
            .collect::<Vec<_>>()
            .join(",")
    )
}

/// A homogeneous part of a request, split but not yet lowered: the top it is
/// lowered onto, the request's line powers and the loop-dot polynomial.
struct PendingPart<'r> {
    top: PTop,
    powers: &'r [i64],
    polynomial: Vec<(Fraction, PMonomial)>,
}

/// The split parts of one channel and the vakint integrals they come from.
#[derive(Default)]
struct ChannelParts<'r> {
    parts: Vec<PendingPart<'r>>,
    integrals: BTreeSet<String>,
}

/// The split parts of a batch, by channel.
type PendingChannels<'r> = BTreeMap<Channel, ChannelParts<'r>>;

/// The checks of [`lower_part`] on one part, without lowering it: the line
/// count, and the dot form of every scalar product in the order
/// [`lower_part`] builds them, cached in `forms`. A part that passes lowers
/// without error with these forms, and the first error is the one
/// [`lower_part`] would return.
fn prepare_part(
    presentation: &Presentation,
    top: PTop,
    powers: &[i64],
    polynomial: &[(Fraction, PMonomial)],
    forms: &mut BTreeMap<PDot, AffineSlotForm>,
) -> Result<(), RustRedEvaluationError> {
    if powers.len() != top.line_count {
        return Err(invalid(format!(
            "{} line powers for {} lines",
            powers.len(),
            top.line_count
        )));
    }
    let loops = usize::from(presentation.loops);
    for (_, monomial) in polynomial {
        for &(dot, _) in &monomial.0 {
            if let std::collections::btree_map::Entry::Vacant(entry) = forms.entry(dot) {
                let (a, b) = dot_rows(dot, loops);
                entry.insert(
                    presentation
                        .express_scalar_product(&a, &b)
                        .map_err(|error| invalid(format!("{dot:?}: {error}")))?,
                );
            }
        }
    }
    Ok(())
}

/// Split one backend request (`docs/p_integrals.md`, section 5) into the
/// parts of `channels`, checked as [`lower_part`] checks them but not
/// lowered; a scaleless request adds nothing.
fn split_request<'r>(
    handle: &PLibraryHandle,
    request: &PropagatorRequest<'r>,
    forms: &mut BTreeMap<PresentationId, BTreeMap<PDot, AffineSlotForm>>,
    channels: &mut PendingChannels<'r>,
) -> Result<(), VakintError> {
    let family = P_FAMILIES
        .iter()
        .position(|family| family.name == request.family.name)
        .ok_or_else(|| invalid(format!("unknown family {}", request.family.name)))?;
    let top = handle.top(family)?;
    let Some(split) = parts::split_parts(request)? else {
        return Ok(()); // scaleless: an exact zero
    };
    let (head, _) = external_head(&request.external)?;
    let presentation = handle.library.presentation(top.presentation);
    let forms = forms.entry(top.presentation).or_default();
    let integral = integral_name(request);
    for part in split {
        prepare_part(presentation, top, request.powers, &part.polynomial, forms)?;
        let channel = channels
            .entry(Channel {
                loops: request.family.loops,
                head: head.clone(),
                scale_power: part.scale_power,
                coefficient: part.coefficient,
            })
            .or_default();
        channel.integrals.insert(integral.clone());
        channel.parts.push(PendingPart {
            top,
            powers: request.powers,
            polynomial: part.polynomial,
        });
    }
    Ok(())
}

/// The lowered inputs of one channel: its parts lowered in order and merged
/// by key, each part's polynomial dropped once it is lowered.
fn lower_channel(
    handle: &PLibraryHandle,
    parts: Vec<PendingPart>,
    forms: &mut BTreeMap<PresentationId, BTreeMap<PDot, AffineSlotForm>>,
) -> Result<BTreeMap<(PresentationId, Vec<i64>), Rational>, RustRedEvaluationError> {
    let mut inputs = BTreeMap::<(PresentationId, Vec<i64>), Rational>::new();
    for PendingPart {
        top,
        powers,
        polynomial,
    } in parts
    {
        let lowered = lower_part(
            handle.library.presentation(top.presentation),
            top,
            powers,
            &polynomial,
            forms.entry(top.presentation).or_default(),
        )?;
        drop(polynomial);
        for (key, value) in lowered {
            let entry = inputs
                .entry((top.presentation, key))
                .or_insert_with(Rational::zero);
            *entry += &value;
        }
    }
    Ok(inputs)
}

/// What [`reduce_channel`] needs of the batch.
struct ChannelReducer<'h> {
    handle: &'h PLibraryHandle,
    options: PRequestOptions,
    step_cache: Option<Arc<OrderedStepCache>>,
    epsilon: Symbol,
}

/// Reduce the lowered `inputs` of one channel onto the masters, without the
/// channel coefficient (`None` when no input is left).
fn reduce_channel(
    reducer: &ChannelReducer,
    channel: &Channel,
    inputs: BTreeMap<(PresentationId, Vec<i64>), Rational>,
    integrals: &BTreeSet<String>,
) -> Result<Option<Atom>, VakintError> {
    let library = &reducer.handle.library;
    let exact = reducer.options.limits.exact_algebra;
    let inputs = inputs
        .into_iter()
        .filter(|(_, value)| !value.is_zero())
        .map(|((presentation, powers), value)| {
            let key = IntegralKey::try_new(powers)
                .map_err(|error| invalid(format!("integral key: {error}")))?;
            Ok((
                presentation,
                key,
                coefficient(library.coefficient_context(), &value, exact)?,
            ))
        })
        .collect::<Result<Vec<_>, RustRedEvaluationError>>()?;
    if inputs.is_empty() {
        return Ok(None);
    }
    let reduction = reduce_inputs(reducer, channel, inputs, integrals)?;
    let rendering = std::time::Instant::now();
    let rendered = render(reducer.handle, &reduction, reducer.epsilon)?;
    with_stats(|stats| stats.render_seconds += rendering.elapsed().as_secs_f64());
    Ok(Some(rendered))
}

/// The reduced values of a batch of backend requests (`docs/p_integrals.md`,
/// section 5), by loop count and external head: one [`PReductionRequest`] per
/// channel (module documentation), the values of one head NOT yet finalized
/// (`reduced` per scale power, comparable with Forcer's `combine_parts` in
/// the tag basis). The request-level entry of the RustRed p path below
/// [`evaluate_ordered_batch`], for crate tests (M2-V5 Tier A).
///
/// Every request is split (and checked as its lowering checks it) before any
/// channel is reduced, but a channel is lowered only just before its
/// reduction, so only one channel's lowered inputs are alive at a time. The
/// values, their order and the errors are those of lowering every channel
/// first ([`reduce_requests_eager`]): rational sums are exact, and lowering
/// fails only where [`prepare_part`] fails, in the same order.
///
/// The channels are folded first ([`reduce_requests_folded`]) unless
/// `VAKINT_P_FOLD` is `off` (`0`, `false`), which keeps one reduction per
/// non-rational coefficient; `basis` and `coefficient` pick the fold mode
/// ([`FoldMode::DEFAULT`] when unset).
pub(crate) fn reduce_requests(
    vakint: &Vakint,
    settings: &VakintSettings,
    requests: &[PropagatorRequest],
    progress: Option<&super::monitor::CallerProgress>,
) -> Result<BTreeMap<(usize, Atom), Vec<ReducedPropagatorValue>>, VakintError> {
    reduce_requests_with(vakint, settings, requests, progress, false)
}

/// [`reduce_requests`] with every value factored: `reduced` is the rendered
/// reduction of a channel and `spectator` its coefficient, so that a caller
/// can expand both apart ([`ReducedPropagatorValue`]). The product of the two
/// is the value of [`reduce_requests`].
pub(crate) fn reduce_requests_factored(
    vakint: &Vakint,
    settings: &VakintSettings,
    requests: &[PropagatorRequest],
    progress: Option<&super::monitor::CallerProgress>,
) -> Result<BTreeMap<(usize, Atom), Vec<ReducedPropagatorValue>>, VakintError> {
    reduce_requests_with(vakint, settings, requests, progress, true)
}

fn reduce_requests_with(
    vakint: &Vakint,
    settings: &VakintSettings,
    requests: &[PropagatorRequest],
    progress: Option<&super::monitor::CallerProgress>,
    factored: bool,
) -> Result<BTreeMap<(usize, Atom), Vec<ReducedPropagatorValue>>, VakintError> {
    match fold_mode() {
        FoldMode::Off => reduce_requests_unfolded(vakint, settings, requests, progress, factored),
        mode => reduce_requests_folded(vakint, settings, requests, progress, mode, factored),
    }
}

/// How [`reduce_requests`] forms its reductions (`VAKINT_P_FOLD`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum FoldMode {
    /// One reduction per (loops, head, scale power, non-rational coefficient).
    Off,
    /// One reduction per (loops, head, scale power, spectator monomial,
    /// epsilon basis function), with rational inputs; the basis functions are
    /// multiplied in after the reduction.
    Basis,
    /// One reduction per (loops, head, scale power, spectator monomial), with
    /// the epsilon-rational parts in the inputs as rational functions of `d`.
    Coefficient,
}

pub(crate) fn fold_mode() -> FoldMode {
    static MODE: std::sync::LazyLock<FoldMode> =
        std::sync::LazyLock::new(|| match std::env::var("VAKINT_P_FOLD").as_deref() {
            Ok("off" | "0" | "false") => FoldMode::Off,
            Ok("basis") => FoldMode::Basis,
            Ok("coefficient") => FoldMode::Coefficient,
            _ => FoldMode::DEFAULT,
        });
    *MODE
}

impl FoldMode {
    /// The mode when `VAKINT_P_FOLD` is unset (or `on`).
    pub(crate) const DEFAULT: FoldMode = FoldMode::Coefficient;
}

/// What the last [`reduce_requests`] call on this thread did (diagnostics).
#[derive(Clone, Debug, Default)]
pub(crate) struct PReduceStats {
    pub channels: usize,
    pub reductions: usize,
    pub pairs: usize,
    pub split_seconds: f64,
    pub lower_seconds: f64,
    pub prepare_seconds: f64,
    pub library_seconds: f64,
    pub render_seconds: f64,
}

thread_local! {
    static STATS: std::cell::RefCell<PReduceStats> = std::cell::RefCell::new(PReduceStats::default());
}

/// The statistics of the last [`reduce_requests`] call on this thread.
pub(crate) fn last_stats() -> PReduceStats {
    STATS.with(|stats| stats.borrow().clone())
}

fn with_stats(update: impl FnOnce(&mut PReduceStats)) {
    STATS.with(|stats| update(&mut stats.borrow_mut()));
}

/// [`reduce_requests`] without folding.
pub(crate) fn reduce_requests_unfolded(
    vakint: &Vakint,
    settings: &VakintSettings,
    requests: &[PropagatorRequest],
    progress: Option<&super::monitor::CallerProgress>,
    factored: bool,
) -> Result<BTreeMap<(usize, Atom), Vec<ReducedPropagatorValue>>, VakintError> {
    with_stats(|stats| *stats = PReduceStats::default());
    let started = std::time::Instant::now();
    let handle = vakint.p_library()?;
    let mut forms = BTreeMap::new();
    let mut channels = PendingChannels::new();
    for request in requests {
        split_request(&handle, request, &mut forms, &mut channels)?;
    }
    with_stats(|stats| {
        stats.split_seconds = started.elapsed().as_secs_f64();
        stats.channels = channels.len();
    });
    let reducer = ChannelReducer {
        handle: &handle,
        options: request_options(settings),
        step_cache: if channels.is_empty() {
            None
        } else {
            vakint.ordered_step_cache(settings, &handle.cache)?
        },
        epsilon: vk_symbol!(&settings.epsilon_symbol),
    };
    let total = channels.len();
    let mut values = BTreeMap::<(usize, Atom), Vec<ReducedPropagatorValue>>::new();
    for (index, (channel, ChannelParts { parts, integrals })) in channels.into_iter().enumerate() {
        if let Some(progress) = progress {
            progress.cohorts(index, total);
        }
        let lowering = std::time::Instant::now();
        let inputs = lower_channel(&handle, parts, &mut forms)?;
        with_stats(|stats| {
            stats.lower_seconds += lowering.elapsed().as_secs_f64();
            stats.pairs += inputs.len();
        });
        let Some(reduced) = reduce_channel(&reducer, &channel, inputs, &integrals)? else {
            continue;
        };
        let Channel {
            loops,
            head,
            scale_power,
            coefficient,
        } = channel;
        let value = if factored {
            ReducedPropagatorValue {
                scale_power,
                reduced,
                spectator: Some(coefficient),
            }
        } else {
            ReducedPropagatorValue {
                scale_power,
                reduced: reduced * &coefficient,
                spectator: None,
            }
        };
        values.entry((loops, head)).or_default().push(value);
    }
    if let Some(progress) = progress {
        progress.cohorts(total, total);
    }
    Ok(values)
}

/// A term `value · basis(ε) · monomial` of a channel coefficient: `monomial`
/// holds the factors free of `ε` (and, as a fallback, any factor that is not
/// a rational function of `ε` alone), `basis` the product of the factors that
/// are rational functions of `ε` alone, `value` the rational number.
fn spectator_terms(coefficient: &Atom, epsilon: Symbol) -> Vec<(Atom, Atom, Rational)> {
    use symbolica::atom::AtomView;
    use symbolica::coefficient::CoefficientView;
    fn epsilon_rational(view: AtomView, epsilon: Symbol) -> bool {
        match view {
            AtomView::Num(_) => parts::rational(view).is_some(),
            AtomView::Var(var) => var.get_symbol() == epsilon,
            AtomView::Add(sum) => sum.iter().all(|term| epsilon_rational(term, epsilon)),
            AtomView::Mul(product) => product.iter().all(|f| epsilon_rational(f, epsilon)),
            AtomView::Pow(power) => {
                let (base, exponent) = power.get_base_exp();
                crate::get_integer_from_atom(exponent).is_some() && epsilon_rational(base, epsilon)
            }
            _ => false,
        }
    }
    let expanded = coefficient.expand();
    let mut terms = Vec::new();
    let mut add = |term: AtomView| {
        let mut value = Rational::one();
        let mut basis = Atom::num(1);
        let mut monomial = Atom::num(1);
        let factors = match term {
            AtomView::Mul(product) => product.iter().collect::<Vec<_>>(),
            other => vec![other],
        };
        for factor in factors {
            if let Some(number) = parts::rational(factor) {
                value *= &fraction(number);
            } else if let AtomView::Num(number) = factor
                && let CoefficientView::Natural(0, _, imaginary, denominator) =
                    number.get_coeff_view()
                && imaginary != 0
            {
                // A purely imaginary number: its rational magnitude goes to
                // the value and the unit to the monomial.
                let magnitude = Rational::from((imaginary, denominator));
                monomial *= factor.to_owned() / Atom::num(magnitude.clone());
                value *= &magnitude;
            } else if !factor.contains_symbol(epsilon) {
                monomial *= factor.to_owned();
            } else if epsilon_rational(factor, epsilon) {
                basis *= factor.to_owned();
            } else {
                monomial *= factor.to_owned();
            }
        }
        if !value.is_zero() {
            terms.push((monomial, basis, value));
        }
    };
    match expanded.as_view() {
        AtomView::Add(sum) => sum.iter().for_each(&mut add),
        other => add(other),
    }
    terms
}

/// A rational function of `ε` alone ([`spectator_terms`]) as a coefficient of
/// the library's context at `ε = (4 - d)/2`.
fn epsilon_coefficient(
    context: &CoefficientContext,
    view: symbolica::atom::AtomView,
    epsilon: Symbol,
    epsilon_at_d: &Coefficient,
    limits: ExactAlgebraLimits,
) -> Result<Coefficient, RustRedEvaluationError> {
    use symbolica::atom::AtomView;
    let native = RustRedEvaluationError::from_native;
    match view {
        AtomView::Num(_) => {
            let value = parts::rational(view).ok_or_else(|| invalid(format!("number {view}")))?;
            coefficient(context, &fraction(value), limits)
        }
        AtomView::Var(var) if var.get_symbol() == epsilon => Ok(epsilon_at_d.clone()),
        AtomView::Add(sum) => {
            let mut total = context.zero();
            for term in sum.iter() {
                let term = epsilon_coefficient(context, term, epsilon, epsilon_at_d, limits)?;
                total = context.try_add(&total, &term, limits).map_err(native)?;
            }
            Ok(total)
        }
        AtomView::Mul(product) => {
            let mut total = context.one();
            for factor in product.iter() {
                let factor = epsilon_coefficient(context, factor, epsilon, epsilon_at_d, limits)?;
                total = context.try_mul(&total, &factor, limits).map_err(native)?;
            }
            Ok(total)
        }
        AtomView::Pow(power) => {
            let (base, exponent) = power.get_base_exp();
            let exponent = crate::get_integer_from_atom(exponent)
                .ok_or_else(|| invalid(format!("power {view}")))?;
            let base = epsilon_coefficient(context, base, epsilon, epsilon_at_d, limits)?;
            let mut total = context.one();
            for _ in 0..exponent.unsigned_abs() {
                total = context.try_mul(&total, &base, limits).map_err(native)?;
            }
            if exponent < 0 {
                total = context
                    .try_div(&context.one(), &total, limits)
                    .map_err(native)?;
            }
            Ok(total)
        }
        _ => Err(invalid(format!(
            "not a rational function of epsilon: {view}"
        ))),
    }
}

/// The lowered inputs of one folded channel, by epsilon basis function.
#[derive(Default)]
struct FoldedInputs {
    by_basis: BTreeMap<Atom, BTreeMap<(PresentationId, Vec<i64>), Rational>>,
    integrals: BTreeSet<String>,
}

/// [`reduce_requests`] with the channels folded: every channel coefficient is
/// written as `Σ value · basis(ε) · monomial` ([`spectator_terms`]), and the
/// channel key keeps only the monomial. The value of a folded channel is
/// `Σ_basis basis(ε) · R(Σ value · inputs) · monomial`, by linearity the sum
/// of the values of the channels it folds. In [`FoldMode::Basis`] every basis
/// function is reduced apart with rational inputs and multiplied in at the
/// render level; in [`FoldMode::Coefficient`] the basis functions enter the
/// inputs as rational functions of `d` and each folded channel is reduced once.
pub(crate) fn reduce_requests_folded(
    vakint: &Vakint,
    settings: &VakintSettings,
    requests: &[PropagatorRequest],
    progress: Option<&super::monitor::CallerProgress>,
    mode: FoldMode,
    factored: bool,
) -> Result<BTreeMap<(usize, Atom), Vec<ReducedPropagatorValue>>, VakintError> {
    with_stats(|stats| *stats = PReduceStats::default());
    let started = std::time::Instant::now();
    let epsilon = vk_symbol!(&settings.epsilon_symbol);
    let handle = vakint.p_library()?;
    let mut forms = BTreeMap::new();
    let mut pending = PendingChannels::new();
    for request in requests {
        split_request(&handle, request, &mut forms, &mut pending)?;
    }
    with_stats(|stats| {
        stats.split_seconds = started.elapsed().as_secs_f64();
        stats.channels = pending.len();
    });
    let lowering = std::time::Instant::now();
    // Lower every channel once and distribute its lowered inputs, scaled by
    // the value of each term of its coefficient, over the folded channels and
    // basis functions. The integral keys are numbered once, and every folded
    // basis accumulates into a dense vector indexed by key number.
    let mut key_ids = HashMap::<(PresentationId, Vec<i64>), usize>::new();
    let mut keys = Vec::<(PresentationId, Vec<i64>)>::new();
    let mut buckets = Vec::<Vec<Rational>>::new();
    let mut targets = BTreeMap::<Channel, (BTreeMap<Atom, usize>, BTreeSet<String>)>::new();
    for (channel, ChannelParts { parts, integrals }) in pending {
        let lowered = lower_channel(&handle, parts, &mut forms)?;
        if lowered.is_empty() {
            continue;
        }
        let lowered = lowered
            .into_iter()
            .filter(|(_, value)| !value.is_zero())
            .map(|(key, value)| {
                let next = keys.len();
                let id = *key_ids.entry(key.clone()).or_insert_with(|| {
                    keys.push(key);
                    next
                });
                (id, value)
            })
            .collect::<Vec<_>>();
        for (monomial, basis, value) in spectator_terms(&channel.coefficient, epsilon) {
            let target = targets
                .entry(Channel {
                    loops: channel.loops,
                    head: channel.head.clone(),
                    scale_power: channel.scale_power,
                    coefficient: monomial,
                })
                .or_default();
            target.1.extend(integrals.iter().cloned());
            let next = buckets.len();
            let bucket = *target.0.entry(basis).or_insert_with(|| {
                buckets.push(Vec::new());
                next
            });
            let bucket = &mut buckets[bucket];
            if bucket.len() < keys.len() {
                bucket.resize(keys.len(), Rational::zero());
            }
            let unit = value == Rational::one();
            for (id, lowered_value) in &lowered {
                if unit {
                    bucket[*id] += lowered_value;
                } else {
                    bucket[*id] += &(lowered_value * &value);
                }
            }
        }
    }
    let mut folded = BTreeMap::<Channel, FoldedInputs>::new();
    for (channel, (by_basis, integrals)) in targets {
        let target = folded.entry(channel).or_default();
        target.integrals = integrals;
        for (basis, bucket) in by_basis {
            let inputs = std::mem::take(&mut buckets[bucket])
                .into_iter()
                .enumerate()
                .filter(|(_, value)| !value.is_zero())
                .map(|(id, value)| (keys[id].clone(), value))
                .collect::<BTreeMap<_, _>>();
            target.by_basis.insert(basis, inputs);
        }
    }
    with_stats(|stats| stats.lower_seconds = lowering.elapsed().as_secs_f64());
    let reducer = ChannelReducer {
        handle: &handle,
        options: request_options(settings),
        step_cache: if folded.is_empty() {
            None
        } else {
            vakint.ordered_step_cache(settings, &handle.cache)?
        },
        epsilon,
    };
    let library = &handle.library;
    let context = library.coefficient_context();
    let exact = reducer.options.limits.exact_algebra;
    let native = RustRedEvaluationError::from_native;
    let epsilon_at_d = {
        let d = context
            .parameter(&handle.dimension.to_canonical_string())
            .or_else(|| context.parameter("d"))
            .ok_or_else(|| invalid("the p library has no dimension parameter"))?;
        let four_minus_d = context
            .try_sub(&context.integer(4), &d, exact)
            .map_err(native)?;
        context
            .try_div(&four_minus_d, &context.integer(2), exact)
            .map_err(native)?
    };
    let total = folded.len();
    let mut values = BTreeMap::<(usize, Atom), Vec<ReducedPropagatorValue>>::new();
    for (
        index,
        (
            channel,
            FoldedInputs {
                by_basis,
                integrals,
            },
        ),
    ) in folded.into_iter().enumerate()
    {
        if let Some(progress) = progress {
            progress.cohorts(index, total);
        }
        let reduced = match mode {
            FoldMode::Coefficient => {
                let preparing = std::time::Instant::now();
                let mut inputs = BTreeMap::<(PresentationId, Vec<i64>), Coefficient>::new();
                for (basis, lowered) in by_basis {
                    let weight = epsilon_coefficient(
                        context,
                        basis.as_view(),
                        epsilon,
                        &epsilon_at_d,
                        exact,
                    )?;
                    for (key, value) in lowered {
                        if value.is_zero() {
                            continue;
                        }
                        let term = context
                            .try_mul(&weight, &coefficient(context, &value, exact)?, exact)
                            .map_err(native)?;
                        match inputs.get_mut(&key) {
                            Some(sum) => {
                                *sum = context.try_add(sum, &term, exact).map_err(native)?
                            }
                            None => {
                                inputs.insert(key, term);
                            }
                        }
                    }
                }
                with_stats(|stats| {
                    stats.prepare_seconds += preparing.elapsed().as_secs_f64();
                    stats.pairs += inputs.len();
                });
                let inputs = inputs
                    .into_iter()
                    .filter(|(_, value)| !value.is_zero())
                    .map(|((presentation, powers), value)| {
                        let key = IntegralKey::try_new(powers)
                            .map_err(|error| invalid(format!("integral key: {error}")))?;
                        Ok((presentation, key, value))
                    })
                    .collect::<Result<Vec<_>, RustRedEvaluationError>>()?;
                if inputs.is_empty() {
                    continue;
                }
                let reduction = reduce_inputs(&reducer, &channel, inputs, &integrals)?;
                let rendering = std::time::Instant::now();
                let terms = render_terms(&handle, &reduction, epsilon)?;
                let reduced = collect_terms(terms);
                with_stats(|stats| stats.render_seconds += rendering.elapsed().as_secs_f64());
                reduced
            }
            _ => {
                let mut combined = BTreeMap::<Atom, Atom>::new();
                for (basis, lowered) in by_basis {
                    with_stats(|stats| stats.pairs += lowered.len());
                    let inputs = lowered
                        .into_iter()
                        .filter(|(_, value)| !value.is_zero())
                        .map(|((presentation, powers), value)| {
                            let key = IntegralKey::try_new(powers)
                                .map_err(|error| invalid(format!("integral key: {error}")))?;
                            Ok((presentation, key, coefficient(context, &value, exact)?))
                        })
                        .collect::<Result<Vec<_>, RustRedEvaluationError>>()?;
                    if inputs.is_empty() {
                        continue;
                    }
                    let reduction = reduce_inputs(&reducer, &channel, inputs, &integrals)?;
                    let rendering = std::time::Instant::now();
                    for (monomial, value) in render_terms(&handle, &reduction, epsilon)? {
                        *combined.entry(monomial).or_default() += value * &basis;
                    }
                    with_stats(|stats| stats.render_seconds += rendering.elapsed().as_secs_f64());
                }
                let rendering = std::time::Instant::now();
                let reduced = collect_terms(combined);
                with_stats(|stats| stats.render_seconds += rendering.elapsed().as_secs_f64());
                reduced
            }
        };
        if reduced.is_zero() {
            continue;
        }
        let Channel {
            loops,
            head,
            scale_power,
            coefficient,
        } = channel;
        let value = if factored {
            ReducedPropagatorValue {
                scale_power,
                reduced,
                spectator: Some(coefficient),
            }
        } else {
            ReducedPropagatorValue {
                scale_power,
                reduced: reduced * &coefficient,
                spectator: None,
            }
        };
        values.entry((loops, head)).or_default().push(value);
    }
    if let Some(progress) = progress {
        progress.cohorts(total, total);
    }
    Ok(values)
}

/// One library reduction of `inputs` for `channel` (the checks of
/// [`reduce_channel`]), timed into the statistics.
fn reduce_inputs(
    reducer: &ChannelReducer,
    channel: &Channel,
    inputs: Vec<(PresentationId, IntegralKey, Coefficient)>,
    integrals: &BTreeSet<String>,
) -> Result<PReduction, VakintError> {
    let reducing = std::time::Instant::now();
    let reduction = {
        let _math = MathSpan::new();
        reducer
            .handle
            .library
            .reduce(p_reduction_request(
                &reducer.options,
                inputs,
                reducer.step_cache.clone(),
            ))
            .map_err(|error| reduction_error(error, integrals))?
    };
    with_stats(|stats| {
        stats.library_seconds += reducing.elapsed().as_secs_f64();
        stats.reductions += 1;
    });
    if reduction.dimension != channel.loops as i64 {
        return Err(invalid(format!(
            "a {}-loop request returned a reduction of dimension {}",
            channel.loops, reduction.dimension
        ))
        .into());
    }
    Ok(reduction)
}

/// `Σ coefficient.together() · monomial` of [`render_terms`]: [`render`].
fn collect_terms(terms: BTreeMap<Atom, Atom>) -> Atom {
    terms
        .into_iter()
        .fold(Atom::Zero, |reduced, (monomial, coefficient)| {
            reduced + coefficient.together() * monomial
        })
}

/// [`reduce_requests`] as it was before the lazy lowering: every channel is
/// lowered into one table before any is reduced. The reference of the tests.
#[cfg(test)]
pub(crate) fn reduce_requests_eager(
    vakint: &Vakint,
    settings: &VakintSettings,
    requests: &[PropagatorRequest],
    progress: Option<&super::monitor::CallerProgress>,
) -> Result<BTreeMap<(usize, Atom), Vec<ReducedPropagatorValue>>, VakintError> {
    let handle = vakint.p_library()?;
    let mut forms = BTreeMap::new();
    let mut channels = Channels::new();
    for request in requests {
        add_request(&handle, request, &mut forms, &mut channels)?;
    }
    let reducer = ChannelReducer {
        handle: &handle,
        options: request_options(settings),
        step_cache: if channels.is_empty() {
            None
        } else {
            vakint.ordered_step_cache(settings, &handle.cache)?
        },
        epsilon: vk_symbol!(&settings.epsilon_symbol),
    };
    let total = channels.len();
    let mut values = BTreeMap::<(usize, Atom), Vec<ReducedPropagatorValue>>::new();
    for (index, (channel, ChannelInputs { inputs, integrals })) in channels.into_iter().enumerate()
    {
        if let Some(progress) = progress {
            progress.cohorts(index, total);
        }
        let Some(reduced) = reduce_channel(&reducer, &channel, inputs, &integrals)? else {
            continue;
        };
        values
            .entry((channel.loops, channel.head))
            .or_default()
            .push(ReducedPropagatorValue {
                scale_power: channel.scale_power,
                reduced: reduced * &channel.coefficient,
                spectator: None,
            });
    }
    if let Some(progress) = progress {
        progress.cohorts(total, total);
    }
    Ok(values)
}

/// Reduce and finalize a batch of matched p contributions in `mode` (module
/// documentation). The result is final: the caller adds it to the expression
/// without RustRed's vacuum master substitution.
pub(crate) fn evaluate_ordered_batch(
    vakint: &Vakint,
    settings: &VakintSettings,
    contributions: &[MatchedContribution],
    mode: RustRedMasterMode,
    progress: Option<&super::monitor::CallerProgress>,
) -> Result<Atom, VakintError> {
    let requests = contributions
        .iter()
        .map(|contribution| {
            contribution
                .matched
                .propagator_request(contribution.numerator.as_view())
                .ok_or_else(|| {
                    VakintError::from(invalid("a contribution without a p-integral witness"))
                })
        })
        .collect::<Result<Vec<_>, VakintError>>()?;
    let mut result = Atom::Zero;
    for ((loops, head), values) in reduce_requests(vakint, settings, &requests, progress)? {
        let finalized = finalize::finalize(&values, loops, &head, settings, mode, "RustRed")?;
        result += simplify_real(finalized.as_view());
    }
    Ok(result)
}

#[cfg(test)]
#[path = "p_integral_tests.rs"]
mod tests;
