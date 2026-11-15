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
//!
//! Entry points (crate-visible): [`supports`] and [`has_p_data`] for the
//! dispatch, [`evaluate_ordered_batch`] for a batch of matched contributions,
//! and [`reduce_requests`], the request level below the finalizer (backend
//! requests in, reduced values per scale power out), for crate tests that
//! compare the RustRed p path with Forcer's frozen values per scale power
//! (M2-V5 Tier A).

use std::collections::{BTreeMap, BTreeSet};
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
use crate::p_integral::master_catalogue::{PMasterCatalogue, PMasterSource};
use crate::p_integral::parts::{self, Fraction, PDot, PMonomial};
use crate::p_integral::{
    P_FAMILIES, P_INTEGRAL_LOOP_COUNTS, PropagatorRequest, ReducedPropagatorValue, external_head,
    finalize, masters, p_gamma_symbol, p_master_symbol, phases,
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
/// the settings (default `modular_polynomial`, the `Z_p(d)` replay, the
/// fastest up to four loops), the arithmetic of five-loop requests
/// ([`super::five_loop_ordered_weight_arithmetic`]: `modular_sampled` under
/// the default selection), the Ordered scheduler on `n_cores` threads and the
/// native reduction limits.
pub(super) struct PRequestOptions {
    pub arithmetic: OrderedWeightArithmetic,
    pub five_loop_arithmetic: OrderedWeightArithmetic,
    pub scheduler: OrderedSchedulerOptions,
    pub limits: ReductionLimits,
    /// vk5l WP4: the coverage policy of every request ([`uncovered_policy`]).
    pub uncovered: PUncoveredPolicy,
    /// The master mode of the settings' RustRed backend (request dumps).
    pub master_mode: Option<RustRedMasterMode>,
}

impl PRequestOptions {
    /// The arithmetic of a request of `loops` loops: five-loop p requests take
    /// `five_loop_arithmetic` (measured, notes arith-5l: `modular_sampled`
    /// 1.06-1.65x faster than `modular_polynomial` on 5L cells and the only
    /// one whose replay of a heavy 5L top finishes), one to four loops
    /// `arithmetic` (`modular_polynomial` 1.08-1.27x faster at 4L). Every
    /// arithmetic returns the same exact result.
    pub fn arithmetic_for(&self, loops: usize) -> OrderedWeightArithmetic {
        if loops == 5 {
            self.five_loop_arithmetic
        } else {
            self.arithmetic
        }
    }
}

pub(super) fn request_options(settings: &VakintSettings) -> PRequestOptions {
    let master_mode = super::request_dump::master_mode(settings);
    PRequestOptions {
        arithmetic: settings.rustred_ordered_weight_arithmetic(),
        five_loop_arithmetic: super::five_loop_ordered_weight_arithmetic(settings),
        scheduler: super::ordered_scheduler(settings),
        limits: settings.rustred_application_options.reduction_limits,
        uncovered: uncovered_policy(settings.rustred_p_uncovered, master_mode),
        master_mode,
    }
}

/// vk5l WP4 (plan VK5L §3.2, decision Q9): the coverage policy of the p
/// requests. [`crate::RustRedPUncoveredPolicy::Reject`] (the default): a key
/// the library cannot reduce is a typed error. `Retain` mirrors the vacuum
/// path's retained unresolved masters: such a key stays in the result as a
/// `PResidual(node, k..)` symbol with a warning naming the refusal, but only
/// under `master_mode=retain` (the vacuum path, too, refuses an unresolved
/// key once masters must be substituted); `Symbolic` and `Numerical` keep
/// `Reject`, so the caller sees the typed refusal rather than an unmapped
/// terminal. A Retain request runs without the step cache (RustRed: a retained
/// terminal's id is local to the process).
fn uncovered_policy(
    setting: crate::RustRedPUncoveredPolicy,
    master_mode: Option<RustRedMasterMode>,
) -> PUncoveredPolicy {
    match (setting, master_mode) {
        (crate::RustRedPUncoveredPolicy::Retain, Some(RustRedMasterMode::Retain)) => {
            PUncoveredPolicy::Retain
        }
        _ => PUncoveredPolicy::Reject,
    }
}

/// The request [`reduce_requests`] sends for one channel: `inputs` with the
/// `options` of the settings, the library's step cache, the coverage policy
/// of the settings ([`uncovered_policy`], `Reject` by default), reduced onto
/// the masters (module documentation).
fn p_reduction_request<'a>(
    options: &PRequestOptions,
    loops: usize,
    inputs: Vec<(PresentationId, IntegralKey, Coefficient)>,
    step_cache: Option<Arc<OrderedStepCache>>,
) -> PReductionRequest<'a> {
    PReductionRequest {
        inputs,
        arithmetic: options.arithmetic_for(loops),
        scheduler: options.scheduler,
        step_cache,
        uncovered: options.uncovered,
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
#[derive(Default)]
struct ChannelInputs {
    inputs: BTreeMap<(PresentationId, Vec<i64>), Rational>,
    integrals: BTreeSet<String>,
}

/// At most this many vakint integrals are named in a reduction error.
const NAMED_INTEGRALS: usize = 8;

/// The vakint integrals of a request, at most [`NAMED_INTEGRALS`], then a count.
fn named_integrals(integrals: &BTreeSet<String>) -> String {
    let mut named = integrals
        .iter()
        .take(NAMED_INTEGRALS)
        .cloned()
        .collect::<Vec<_>>();
    if integrals.len() > NAMED_INTEGRALS {
        named.push(format!("{} more", integrals.len() - NAMED_INTEGRALS));
    }
    named.join(", ")
}

/// Highest arity RustRed compiles by default (`RUSTRED_RUNTIME_ARITIES` unset:
/// `1..=16`); five-loop p families have 20 slots.
const DEFAULT_RUNTIME_ARITY: usize = 16;

/// The RustRed build a reproduction command of a key of `arity` slots needs:
/// empty up to [`DEFAULT_RUNTIME_ARITY`], else the arity-20 build (the
/// default build refuses a 20-slot chart).
fn reproduction_build(arity: usize) -> String {
    if arity <= DEFAULT_RUNTIME_ARITY {
        String::new()
    } else {
        format!(
            " (the reproduction command needs a RustRed build that compiles arity {arity}: \
             RUSTRED_RUNTIME_ARITIES=1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,{arity} with \
             --features rustred-app/runtime-arity-selection)"
        )
    }
}

/// A failed p reduction. A key the library does not cover (`Uncovered`,
/// `MissingOwner`) names the vakint integrals of the channel next to RustRed's
/// owner class and key. vk5l WP4: the typed refusals of boxed and open nodes
/// (`OutsideBox`, `OutsideCertifiedSet`, `NodeNotClosed`: node, key, violated
/// bound, original request and reproduction command) pass through unchanged
/// as [`RustRedEvaluationError::PRefusal`] with the vakint integrals and the
/// build the reproduction needs. Any other error maps as for the vacuum
/// programs.
fn reduction_error(error: PError, integrals: &BTreeSet<String>) -> RustRedEvaluationError {
    match &error {
        PError::OutsideBox { key, .. }
        | PError::OutsideCertifiedSet { key, .. }
        | PError::NodeNotClosed { key, .. } => {
            let build = reproduction_build(key.powers().len());
            RustRedEvaluationError::PRefusal {
                error: Box::new(error),
                integrals: named_integrals(integrals),
                build,
            }
        }
        PError::Uncovered { .. } | PError::MissingOwner { .. } => {
            RustRedEvaluationError::Reduction {
                detail: format!(
                    "{error}, reducing the p-integrals {} with the RustRed p library",
                    named_integrals(integrals)
                ),
            }
        }
        _ => RustRedEvaluationError::from_native(error),
    }
}

/// vk5l WP4: one p request (a channel or a joint group of `loops` loops):
/// reduced with the settings' `options`, written first when
/// `VAKINT_DUMP_ORDERED_REQUEST` is set (the request, then its result:
/// [`super::request_dump::PRequestDump`]; with
/// `VAKINT_DUMP_ORDERED_REQUEST_ONLY=1` the reduction is skipped and the
/// evaluation fails after the first written request), its failure mapped by
/// [`reduction_error`].
pub(super) fn reduce_one(
    library: &rustred::propagator::library::PLibrary,
    options: &PRequestOptions,
    loops: usize,
    inputs: Vec<(PresentationId, IntegralKey, Coefficient)>,
    step_cache: Option<Arc<OrderedStepCache>>,
    integrals: &BTreeSet<String>,
    label: &str,
) -> Result<PReduction, RustRedEvaluationError> {
    let dump = super::request_dump::PRequestDump::begin(super::request_dump::PDumpRequest {
        library,
        loops,
        label,
        arithmetic: options.arithmetic_for(loops),
        uncovered: options.uncovered,
        master_mode: options.master_mode,
        inputs: &inputs,
        integrals,
    })?;
    if let Some(dump) = &dump
        && dump.only()
    {
        return Err(super::request_dump::skipped(1));
    }
    let result = library.reduce(p_reduction_request(options, loops, inputs, step_cache));
    if let Ok(reduction) = &result
        && reduction.relation_cache != Default::default()
    {
        log::info!(
            "RustRed p request {label}: {}{}",
            super::p_relations::describe(&reduction.relation_cache),
            if reduction.outside_certified.is_empty() {
                String::new()
            } else {
                format!(" outside_certified_keys={}", reduction.outside_certified.len())
            }
        );
    }
    if let Some(dump) = dump {
        dump.finish(library, &result)?;
    }
    result.map_err(|error| reduction_error(error, integrals))
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
/// `PBasis::Master(name)` as `PMaster(name)` for a tabulated Forcer master and
/// `PBasis::Residual` as `PResidual(node, k..)` (t162-export);
/// the coefficient of each tag and master is one simplified rational function
/// of `ε` (`together`). An unknown master fails closed,
/// and so does a nonzero condition that vanishes identically at
/// `d = 4 - 2ε`.
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
            // t162-export: a residual terminal no basis map evaluates (a boxed
            // node's unknown integral, a five-loop LP master) stays a symbol;
            // the finalizer refuses it unless master_mode=retain.
            PBasis::Residual(terminal) => handle.residual_atom(*terminal)?,
        };
        let mut tag = Atom::num(1);
        for &(m, e) in term.tag.exponents() {
            tag *= p_gamma_symbol()
                .call((i64::from(m),))
                .pow(Atom::num(i64::from(e)));
        }
        *coefficients.entry(tag * basis).or_default() += at_dimension(&term.coefficient);
    }
    Ok(coefficients
        .into_iter()
        .fold(Atom::Zero, |reduced, (monomial, coefficient)| {
            reduced + coefficient.together() * monomial
        }))
}

/// The engine inputs of one channel: nonzero coefficients in the library's
/// context.
fn channel_inputs(
    library: &rustred::propagator::library::PLibrary,
    inputs: BTreeMap<(PresentationId, Vec<i64>), Rational>,
    exact: ExactAlgebraLimits,
) -> Result<Vec<(PresentationId, IntegralKey, Coefficient)>, RustRedEvaluationError> {
    inputs
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
        .collect()
}

/// Whether a batch of `channels` requests on `threads` workers runs its
/// requests side by side, each on one thread, instead of one after the
/// other on all threads: when there are at least two requests per worker
/// (many small requests pay a thread pool and a window dispatch each).
/// `VAKINT_P_ACROSS_CHANNELS=on|off` overrides the choice (benchmarks).
fn across_channels(threads: usize, channels: usize) -> bool {
    static MODE: std::sync::OnceLock<Option<bool>> = std::sync::OnceLock::new();
    let forced = *MODE.get_or_init(
        || match std::env::var("VAKINT_P_ACROSS_CHANNELS").as_deref() {
            Ok("on" | "1") => Some(true),
            Ok("off" | "0") => Some(false),
            _ => None,
        },
    );
    threads > 1 && channels > 1 && forced.unwrap_or(channels >= 2 * threads)
}

/// The lowered inputs of a batch, by channel.
type Channels = BTreeMap<Channel, ChannelInputs>;

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
    let integral = format!(
        "{}({})",
        request.family.name,
        request
            .powers
            .iter()
            .map(i64::to_string)
            .collect::<Vec<_>>()
            .join(",")
    );
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

/// The reduced values of a batch of backend requests (`docs/p_integrals.md`,
/// section 5), by loop count and external head: one [`PReductionRequest`] per
/// channel (module documentation), the values of one head NOT yet finalized
/// (`reduced` per scale power, comparable with Forcer's `combine_parts` in
/// the tag basis). The request-level entry of the RustRed p path below
/// [`evaluate_ordered_batch`], for crate tests (M2-V5 Tier A).
pub(crate) fn reduce_requests(
    vakint: &Vakint,
    settings: &VakintSettings,
    requests: &[PropagatorRequest],
    progress: Option<&super::monitor::CallerProgress>,
) -> Result<BTreeMap<(usize, Atom), Vec<ReducedPropagatorValue>>, VakintError> {
    // Before the library is loaded or used: a boxed node binds the relation
    // cache roots on its first use (p_relations.rs).
    super::p_relations::install(settings);
    let handle = vakint.p_library()?;
    let library = &handle.library;
    let options = request_options(settings);
    let exact = options.limits.exact_algebra;
    let mut forms = BTreeMap::new();
    let mut channels = Channels::new();
    let started = phases::start();
    for request in requests {
        add_request(&handle, request, &mut forms, &mut channels)?;
    }
    phases::end("p.split_and_lower", started);
    let epsilon = vk_symbol!(&settings.epsilon_symbol);
    let step_cache = phases::time("p.attach_step_cache", || {
        if channels.is_empty() {
            Ok(None)
        } else {
            vakint.ordered_step_cache(settings, &handle.cache)
        }
    })?;
    let (mut convert, mut reduce, mut rendering) = (
        phases::Accumulator::default(),
        phases::Accumulator::default(),
        phases::Accumulator::default(),
    );
    if phases::enabled() {
        let union = channels
            .values()
            .flat_map(|inputs| inputs.inputs.keys())
            .collect::<BTreeSet<_>>();
        let presentations = union.iter().map(|(id, _)| *id).collect::<BTreeSet<_>>();
        let integrals = channels
            .values()
            .flat_map(|inputs| inputs.integrals.iter())
            .collect::<BTreeSet<_>>();
        eprintln!(
            "vakint-phase-info p_batch requests={} channels={} inputs={} distinct_keys={} presentations={} vakint_integrals={}",
            requests.len(),
            channels.len(),
            channels
                .values()
                .map(|inputs| inputs.inputs.len())
                .sum::<usize>(),
            union.len(),
            presentations.len(),
            integrals.len()
        );
    }
    let mut values = BTreeMap::<(usize, Atom), Vec<ReducedPropagatorValue>>::new();
    let joint_mode = joint::mode();
    if joint_mode != joint::JointMode::Off && !channels.is_empty() {
        let groups = joint::groups(library, &channels, epsilon, exact)?;
        let check = joint_mode == joint::JointMode::Check;
        let reductions = reduce.time(|| {
            joint::reduce_groups(
                library,
                &options,
                &step_cache,
                &groups,
                &channels,
                check,
                exact,
            )
        })?;
        if !check {
            for (group, reduction) in groups.iter().zip(&reductions) {
                for (channel, _) in &group.channels {
                    channels.remove(channel);
                }
                let Some(reduction) = reduction else {
                    continue; // the summed inputs cancel: an exact zero
                };
                if reduction.dimension != group.loops as i64 {
                    return Err(invalid(format!(
                        "a {}-loop request returned a reduction of dimension {}",
                        group.loops, reduction.dimension
                    ))
                    .into());
                }
                let reduced =
                    rendering.time(|| render(&handle, reduction, epsilon))? * &group.monomial;
                values
                    .entry((group.loops, group.head.clone()))
                    .or_default()
                    .push(ReducedPropagatorValue {
                        scale_power: group.scale_power,
                        reduced,
                    });
            }
        }
    }
    let total = channels.len();
    let threads = settings.n_cores.max(1);
    if across_channels(threads, total) {
        // Many small requests (physics numerators: one per distinct
        // coefficient): each on one thread, the requests spread over the
        // workers, results consumed in channel order. A request returns the
        // same terms on any thread count and whatever its step cache holds,
        // so the values do not depend on the scheduling.
        let jobs = channels
            .into_iter()
            .map(|(channel, ChannelInputs { inputs, integrals })| {
                let inputs = convert.time(|| channel_inputs(library, inputs, exact))?;
                Ok((channel, inputs, integrals))
            })
            .collect::<Result<Vec<_>, RustRedEvaluationError>>()?
            .into_iter()
            .filter(|(_, inputs, _)| !inputs.is_empty())
            .collect::<Vec<_>>();
        let single = PRequestOptions {
            scheduler: OrderedSchedulerOptions {
                threads: 1,
                ..options.scheduler
            },
            ..options
        };
        let reductions = reduce.time(|| {
            let _math = MathSpan::new();
            super::parallel::map_ordered_by_cost(
                threads,
                &jobs,
                |(_, inputs, _)| inputs.len(),
                |index, (channel, inputs, integrals)| {
                    let started = phases::per_request().then(std::time::Instant::now);
                    let reduction = reduce_one(
                        library,
                        &single,
                        channel.loops,
                        inputs.clone(),
                        step_cache.clone(),
                        integrals,
                        &format!("channel {index} scale_power={}", channel.scale_power),
                    );
                    if let (Some(started), Ok(reduction)) = (started, &reduction) {
                        eprintln!(
                            "vakint-p-request index={index} inputs={} terms={} seconds={:.6}",
                            inputs.len(),
                            reduction.terms.len(),
                            started.elapsed().as_secs_f64()
                        );
                    }
                    reduction
                },
            )
        });
        let reductions = reductions.into_iter().collect::<Result<Vec<_>, _>>()?;
        for ((channel, _, _), reduction) in jobs.iter().zip(&reductions) {
            if reduction.dimension != channel.loops as i64 {
                return Err(invalid(format!(
                    "a {}-loop request returned a reduction of dimension {}",
                    channel.loops, reduction.dimension
                ))
                .into());
            }
        }
        let rendered = rendering.time(|| {
            super::parallel::map_ordered(threads, &reductions, |_, reduction| {
                render(&handle, reduction, epsilon)
            })
        });
        for ((channel, _, _), reduced) in jobs.into_iter().zip(rendered) {
            values
                .entry((channel.loops, channel.head))
                .or_default()
                .push(ReducedPropagatorValue {
                    scale_power: channel.scale_power,
                    reduced: reduced? * &channel.coefficient,
                });
        }
        if let Some(progress) = progress {
            progress.cohorts(total, total);
        }
        convert.report("p.input_coefficients");
        reduce.report("p.rustred_reduce");
        rendering.report("p.render");
        return Ok(values);
    }
    for (index, (channel, ChannelInputs { inputs, integrals })) in channels.into_iter().enumerate()
    {
        if let Some(progress) = progress {
            progress.cohorts(index, total);
        }
        let inputs = convert.time(|| channel_inputs(library, inputs, exact))?;
        if inputs.is_empty() {
            continue;
        }
        let input_count = inputs.len();
        let before = reduce.total;
        let reduction = reduce.time(|| {
            let _math = MathSpan::new();
            reduce_one(
                library,
                &options,
                channel.loops,
                inputs,
                step_cache.clone(),
                &integrals,
                &format!("channel {index} scale_power={}", channel.scale_power),
            )
        })?;
        if phases::per_request() {
            eprintln!(
                "vakint-p-request index={index} inputs={input_count} terms={} seconds={:.6}",
                reduction.terms.len(),
                (reduce.total - before).as_secs_f64()
            );
        }
        if reduction.dimension != channel.loops as i64 {
            return Err(invalid(format!(
                "a {}-loop request returned a reduction of dimension {}",
                channel.loops, reduction.dimension
            ))
            .into());
        }
        let reduced =
            rendering.time(|| render(&handle, &reduction, epsilon))? * &channel.coefficient;
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
    convert.report("p.input_coefficients");
    reduce.report("p.rustred_reduce");
    rendering.report("p.render");
    Ok(values)
}

/// The p master catalogue that finalizes the `loops`-loop channels of
/// `vakint`. The five-loop dispatch renders residual terminals as
/// `PResidual(node, k..)`; with a catalogue the finalizer maps them onto its
/// masters and applies the catalogue guards (`p_integral::master_catalogue`);
/// without one `Symbolic` and `Numerical` refuse them with
/// `PMasterCatalogueError::NoCatalogue` and `Retain` returns them. vk5l WP4:
/// the engine's catalogue (`Vakint::set_p_master_catalogue`) for five loops;
/// one to four loops keep the Forcer table. No catalogue ships yet (decision
/// Q3: Retain-only first); WP6 sets the catalogue of the loaded p area.
fn p_master_catalogue(vakint: &Vakint, loops: usize) -> Option<&PMasterCatalogue> {
    if loops >= 5 {
        vakint.p_master_catalogue()
    } else {
        None
    }
}

/// The vakint integrals of `requests` by loop count and external head
/// (`FAMILY(a_1,..)`), named by the finalizer's guards.
fn request_names(
    requests: &[PropagatorRequest],
) -> Result<BTreeMap<(usize, Atom), BTreeSet<String>>, VakintError> {
    let mut names = BTreeMap::<(usize, Atom), BTreeSet<String>>::new();
    for request in requests {
        let (head, _) = external_head(&request.external)?;
        let powers = request
            .powers
            .iter()
            .map(i64::to_string)
            .collect::<Vec<_>>()
            .join(",");
        names
            .entry((request.family.loops, head))
            .or_default()
            .insert(format!("{}({powers})", request.family.name));
    }
    Ok(names)
}

/// Finalize the reduced values of a batch (`reduce_requests`) in `mode`: per
/// loop count and head, the shared p finalizer with the master source of
/// that loop count ([`p_master_catalogue`]) and the request's vakint
/// integrals in its messages.
pub(crate) fn finalize_reduced(
    vakint: &Vakint,
    settings: &VakintSettings,
    reduced: BTreeMap<(usize, Atom), Vec<ReducedPropagatorValue>>,
    names: &BTreeMap<(usize, Atom), BTreeSet<String>>,
    mode: RustRedMasterMode,
) -> Result<Atom, VakintError> {
    let mut result = Atom::Zero;
    for ((loops, head), values) in reduced {
        let catalogue = p_master_catalogue(vakint, loops);
        let request = match names.get(&(loops, head.clone())) {
            Some(integrals) if !integrals.is_empty() => format!(
                "the {loops}-loop p-integrals {} of external {head}",
                named_integrals(integrals)
            ),
            _ => format!("the {loops}-loop p-integrals of external {head}"),
        };
        let finalized = phases::time("p.finalize", || {
            finalize::finalize_with(
                &values,
                loops,
                &head,
                settings,
                mode,
                "RustRed",
                PMasterSource::new(catalogue, &request),
            )
        })?;
        result += phases::time("p.final_simplify_real", || {
            simplify_real(finalized.as_view())
        });
    }
    Ok(result)
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
    let started = phases::start();
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
    phases::end("p.propagator_requests", started);
    let names = request_names(&requests)?;
    let reduced = reduce_requests(vakint, settings, &requests, progress)?;
    finalize_reduced(vakint, settings, reduced, &names, mode)
}

#[path = "p_joint.rs"]
mod joint;

#[cfg(test)]
#[path = "p_integral_tests.rs"]
mod tests;
