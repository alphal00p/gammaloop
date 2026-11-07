//! The RustRed p backend: massless propagator-type integrals reduced by the
//! engine's RustRed p library program (M2 plan, sections 3.3-3.5, M2-V3;
//! `docs/p_integrals.md`, sections 5-6).
//!
//! A batch of matched p contributions is split into the homogeneous parts of
//! the shared numerator split (`p_integral::parts`, the one Forcer renders),
//! grouped into channels keyed by (loop count, external head `p(i)`, scale
//! power `w`, non-rational coefficient), lowered onto the slots of the family's
//! top presentation (lines, then ISPs; a numerator scalar product is the
//! presentation's own affine form in its slots, Minkowski with `q·q = -1`) and
//! reduced with ONE [`PReductionRequest`] per channel onto the masters:
//! `PReduction` is `(εG(1,1))^-L Σ c·I_M(x)` at `Q² = 1` (convention C3), which
//! is exactly vakint's `reduced = (-1)^(Σa+r)·C_G`. The Wick sign of a part
//! (`PPart::wick_sign`) belongs to the Forcer path only and is never applied
//! here (critique-phys SF3). The terms are rendered as
//! `Σ coefficient(d = 4-2ε) · Π PGamma(m)^e · (1 | PMaster(name))` (convention
//! C4) and finalized by the shared p-finalizer with the backend name `RustRed`
//! (normalization, `(Q²)^w`, Laurent expansion and the `Oep` guard).

use std::collections::BTreeMap;

use rustred::algebra::{Coefficient, CoefficientContext, ExactAlgebraLimits};
use rustred::family::IntegralKey;
use rustred::propagator::PresentationId;
use rustred::propagator::library::{PBasis, PReduction, PReductionRequest, PUncoveredPolicy};
use rustred::propagator::presentation::{AffineSlotForm, Presentation};
use rustred::reduction::{OrderedSchedulerOptions, OrderedWeightArithmetic, ReductionLimits};
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
pub(super) fn supports(topology: &Topology) -> bool {
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
pub(super) fn has_p_data(vakint: &Vakint, loops: usize) -> bool {
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

/// One reduction request: the inputs of one loop count, external head, scale
/// power and non-rational coefficient.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
struct Channel {
    loops: usize,
    head: Atom,
    scale_power: i64,
    coefficient: Atom,
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
/// `PBasis::Master(name)` as `PMaster(name)` for a tabulated Forcer master. A
/// residual terminal or an unknown master fails closed, and so does a
/// nonzero condition that vanishes identically at `d = 4 - 2ε`.
fn render(
    handle: &PLibraryHandle,
    reduction: &PReduction,
    epsilon: Symbol,
) -> Result<Atom, RustRedEvaluationError> {
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
    let mut reduced = Atom::Zero;
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
        reduced += at_dimension(&term.coefficient) * tag * basis;
    }
    Ok(reduced)
}

/// The lowered inputs of a batch, by channel.
type Channels = BTreeMap<Channel, BTreeMap<(PresentationId, Vec<i64>), Rational>>;

/// Lower one backend request (`docs/p_integrals.md`, section 5) into
/// `channels`; a scaleless request adds nothing.
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
    for part in split {
        let lowered = lower_part(presentation, top, request.powers, &part.polynomial, forms)?;
        let inputs = channels
            .entry(Channel {
                loops: request.family.loops,
                head: head.clone(),
                scale_power: part.scale_power,
                coefficient: part.coefficient,
            })
            .or_default();
        for (key, value) in lowered {
            let entry = inputs
                .entry((top.presentation, key))
                .or_insert_with(Rational::zero);
            *entry += &value;
        }
    }
    Ok(())
}

/// The reduced values of a batch of requests, by loop count and external
/// head: one [`PReductionRequest`] per channel (module documentation).
fn reduce_requests(
    vakint: &Vakint,
    settings: &VakintSettings,
    requests: &[PropagatorRequest],
    progress: Option<&super::monitor::CallerProgress>,
) -> Result<BTreeMap<(usize, Atom), Vec<ReducedPropagatorValue>>, VakintError> {
    let handle = vakint.p_library()?;
    let library = &handle.library;
    let options = request_options(settings);
    let exact = options.limits.exact_algebra;
    let mut forms = BTreeMap::new();
    let mut channels = Channels::new();
    for request in requests {
        add_request(&handle, request, &mut forms, &mut channels)?;
    }
    let epsilon = vk_symbol!(&settings.epsilon_symbol);
    let step_cache = if channels.is_empty() {
        None
    } else {
        vakint.ordered_step_cache(settings, &handle.cache)?
    };
    let total = channels.len();
    let mut values = BTreeMap::<(usize, Atom), Vec<ReducedPropagatorValue>>::new();
    for (index, (channel, inputs)) in channels.into_iter().enumerate() {
        if let Some(progress) = progress {
            progress.cohorts(index, total);
        }
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
            continue;
        }
        let reduction = {
            let _math = MathSpan::new();
            library
                .reduce(PReductionRequest {
                    inputs,
                    arithmetic: options.arithmetic,
                    scheduler: options.scheduler,
                    step_cache: step_cache.clone(),
                    uncovered: PUncoveredPolicy::Reject,
                    onto_masters: true,
                    limits: options.limits,
                    inspector: None,
                })
                .map_err(RustRedEvaluationError::from_native)?
        };
        if reduction.dimension != channel.loops as i64 {
            return Err(invalid(format!(
                "a {}-loop request returned a reduction of dimension {}",
                channel.loops, reduction.dimension
            ))
            .into());
        }
        let reduced = render(&handle, &reduction, epsilon)? * &channel.coefficient;
        values
            .entry((channel.loops, channel.head))
            .or_default()
            .push(ReducedPropagatorValue {
                scale_power: channel.scale_power,
                reduced,
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
pub(super) fn evaluate_ordered_batch(
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
