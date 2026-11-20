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
pub(crate) struct Channel {
    pub(crate) loops: usize,
    pub(crate) head: Atom,
    pub(crate) scale_power: i64,
    pub(crate) coefficient: Atom,
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
            super::p_relations::count_refusal();
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
    if let Ok(reduction) = &result {
        super::p_relations::count_table_steps(reduction.statistics.work.rule_applications() as u64);
    }
    if let Ok(reduction) = &result
        && reduction.relation_cache != Default::default()
    {
        log::info!(
            "RustRed p request {label}: {}{}",
            super::p_relations::describe(&reduction.relation_cache),
            if reduction.outside_certified.is_empty() {
                String::new()
            } else {
                format!(
                    " outside_certified_keys={}",
                    reduction.outside_certified.len()
                )
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

/// vk5l WP4: the lowering of one unit loop-dot monomial `Π dot^e` of a
/// sector, `Σ c · Π D_s^m`: packed slot exponents ([`PACKED_BITS`] per
/// slot) when they fit, else exponent vectors. Coefficients accumulate in a
/// hash map of packed exponents (profile of a six-loop census request: 52 %
/// of the time went into a `BTreeMap` of key vectors).
#[derive(Debug)]
enum UnitLowering {
    Packed(Vec<(u128, Rational)>),
    Unpacked(Vec<(Vec<u32>, Rational)>),
}

/// The lowering of the unit monomial `dots` (sorted) of a sector.
fn lower_unit(
    top: PTop,
    dots: &[(PDot, i64)],
    forms: &BTreeMap<PDot, AffineSlotForm>,
) -> UnitLowering {
    let degree = dots.iter().map(|(_, exponent)| *exponent).sum();
    if packable(top, degree) {
        UnitLowering::Packed(lower_unit_packed(dots, forms))
    } else {
        let mut product = BTreeMap::from([(vec![0u32; top.arity], Rational::one())]);
        for &(dot, exponent) in dots {
            for _ in 0..exponent {
                product = multiply(&product, &forms[&dot]);
            }
        }
        UnitLowering::Unpacked(product.into_iter().collect())
    }
}

/// Bits per slot exponent in a packed slot monomial.
const PACKED_BITS: u32 = 6;

/// Whether the slot monomials of `top` up to `degree` fit [`PACKED_BITS`].
fn packable(top: PTop, degree: i64) -> bool {
    top.arity as u32 * PACKED_BITS <= 128 && degree < (1 << PACKED_BITS)
}

/// `Π form(dot)^e` over the slots of `top`, packed (slot `s` exponent in
/// bits `PACKED_BITS·s ..`); the caller checks [`packable`].
fn lower_unit_packed(
    dots: &[(PDot, i64)],
    forms: &BTreeMap<PDot, AffineSlotForm>,
) -> Vec<(u128, Rational)> {
    let mut product = HashMap::<u128, Rational>::from([(0u128, Rational::one())]);
    for (dot, exponent) in dots {
        let form = &forms[dot];
        for _ in 0..*exponent {
            let mut next =
                HashMap::<u128, Rational>::with_capacity(product.len() * (form.terms.len() + 1));
            for (monomial, coefficient) in &product {
                if !form.constant.is_zero() {
                    *next.entry(*monomial).or_insert_with(Rational::zero) +=
                        &(coefficient * &form.constant);
                }
                for (slot, value) in &form.terms {
                    let shifted = monomial + (1u128 << (PACKED_BITS * u32::from(*slot)));
                    *next.entry(shifted).or_insert_with(Rational::zero) += &(coefficient * value);
                }
            }
            next.retain(|_, value| !value.is_zero());
            product = next;
        }
    }
    product.into_iter().collect()
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
    let mut packed = HashMap::<u128, Rational>::new();
    let mut unpacked = HashMap::<Vec<u32>, Rational>::new();
    for (value, monomial) in polynomial {
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
        let mut dots = monomial.0.clone();
        dots.sort();
        let unit = lower_unit(top, &dots, forms);
        let value = fraction(*value);
        match &unit {
            UnitLowering::Packed(terms) => {
                for (monomial, coefficient) in terms {
                    *packed.entry(*monomial).or_insert_with(Rational::zero) +=
                        &(coefficient * &value);
                }
            }
            UnitLowering::Unpacked(terms) => {
                for (monomial, coefficient) in terms {
                    *unpacked
                        .entry(monomial.clone())
                        .or_insert_with(Rational::zero) += &(coefficient * &value);
                }
            }
        }
    }
    let key_of = |exponent: &dyn Fn(usize) -> u32| {
        (0..top.arity)
            .map(|slot| {
                let base = if slot < top.line_count {
                    powers[slot].max(0)
                } else {
                    0
                };
                base - i64::from(exponent(slot))
            })
            .collect::<Vec<_>>()
    };
    let mask = (1u128 << PACKED_BITS) - 1;
    let mut lowered = BTreeMap::<Vec<i64>, Rational>::new();
    for (monomial, coefficient) in packed {
        if coefficient.is_zero() {
            continue;
        }
        let key = key_of(&|slot| ((monomial >> (PACKED_BITS * slot as u32)) & mask) as u32);
        *lowered.entry(key).or_insert_with(Rational::zero) += &coefficient;
    }
    for (monomial, coefficient) in unpacked {
        if coefficient.is_zero() {
            continue;
        }
        let key = key_of(&|slot| monomial[slot]);
        *lowered.entry(key).or_insert_with(Rational::zero) += &coefficient;
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
    Ok(collect_terms(render_terms(handle, reduction, epsilon)?))
}

/// vk5l-physrun: the Retain value of a raw `p5l_library` result text (lines
/// `Unit GammaTag([(m, e), ..]) : <coefficient in d>` and
/// `PResidual(NODE, [k_1, ..]) GammaTag([..]) : <coefficient>`; other lines
/// are ignored), built exactly as [`render_terms`] renders a reduction:
/// `Π PGamma(m)^e · basis · coefficient(d = 4 - 2ε)` with `basis = 1` for
/// `Unit` and `PResidual(NODE, k_1, ..)` for a residual terminal. The master
/// substitution is NOT done here: pass the value to
/// [`Vakint::refinalise_p_retained`] (the finaliser of [`finalize_reduced`]).
pub fn retained_from_p5l_output(text: &str, epsilon: Symbol) -> Result<Atom, VakintError> {
    let bad = |line: &str, what: &str| {
        VakintError::InvalidIntegralFormat(format!("p5l output line {:.80}: {what}", line))
    };
    let dimension = Atom::num(4) - Atom::num(2) * Atom::var(epsilon);
    let d = vk_symbol!("d");
    let mut value = Atom::Zero;
    let mut terms = 0usize;
    for line in text.lines() {
        let line = line.trim();
        if !(line.starts_with("Unit GammaTag(") || line.starts_with("PResidual(")) {
            continue;
        }
        let (head, coefficient) = line
            .split_once(" : ")
            .ok_or_else(|| bad(line, "no ' : '"))?;
        let (basis_text, tag_text) = head
            .rsplit_once(" GammaTag([")
            .ok_or_else(|| bad(line, "no GammaTag"))?;
        let tag_text = tag_text
            .strip_suffix("])")
            .ok_or_else(|| bad(line, "unterminated GammaTag"))?;
        let basis = if basis_text == "Unit" {
            Atom::num(1)
        } else {
            let inner = basis_text
                .strip_prefix("PResidual(")
                .and_then(|rest| rest.strip_suffix("])"))
                .ok_or_else(|| bad(line, "malformed PResidual"))?;
            let (node, key) = inner
                .split_once(", [")
                .ok_or_else(|| bad(line, "malformed PResidual key"))?;
            let mut builder =
                symbolica::atom::FunctionBuilder::new(crate::p_integral::p_residual_symbol())
                    .add_arg(Atom::var(vk_symbol!(node.trim())));
            for power in key.split(',').filter(|p| !p.trim().is_empty()) {
                let power: i64 = power.trim().parse().map_err(|_| bad(line, "a key power"))?;
                builder = builder.add_arg(Atom::num(power));
            }
            builder.finish()
        };
        let mut tag = Atom::num(1);
        for pair in tag_text.split("), (").filter(|p| !p.trim().is_empty()) {
            let pair = pair.trim_matches(|c| c == '(' || c == ')');
            let (m, e) = pair
                .split_once(',')
                .ok_or_else(|| bad(line, "a Gamma tag pair"))?;
            let m: i64 = m.trim().parse().map_err(|_| bad(line, "a Gamma tag m"))?;
            let e: i64 = e
                .trim()
                .parse()
                .map_err(|_| bad(line, "a Gamma tag exponent"))?;
            tag *= crate::p_integral::p_gamma_symbol()
                .call((m,))
                .pow(Atom::num(e));
        }
        let coefficient = crate::utils::vakint_macros::vk_parse!(coefficient)
            .map_err(|error| bad(line, &format!("coefficient: {error}")))?
            .replace(Atom::var(d).to_pattern())
            .with(dimension.to_pattern());
        value += tag * basis * coefficient;
        terms += 1;
    }
    if terms == 0 {
        return Err(VakintError::InvalidIntegralFormat(
            "p5l output: no result term".into(),
        ));
    }
    Ok(value)
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
pub(crate) fn integral_name(request: &PropagatorRequest) -> String {
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

/// Lowered inputs by integral key.
type LoweredInputs = BTreeMap<(PresentationId, Vec<i64>), Rational>;

/// One request split, lowered and folded (the value of [`RequestMemo`]): its
/// inputs by folded channel (spectator monomial) and epsilon basis function,
/// the number of channels its split gave, and its vakint integral. The
/// folded inputs of a batch are the sums of those of its requests.
pub(crate) struct FoldedRequest {
    targets: BTreeMap<Channel, BTreeMap<Atom, LoweredInputs>>,
    split_channels: usize,
    integral: String,
    pairs: usize,
}

/// The identity of a request: what its split and lowering depend on (the
/// library is the memo's owner).
#[derive(Clone, PartialEq, Eq, Hash)]
struct RequestKey {
    /// The fold separates the rational functions of this symbol.
    epsilon: Symbol,
    family: &'static str,
    powers: Vec<i64>,
    external: Atom,
    numerator: Atom,
    /// lane prune: the direct lowering drops zero-sector keys in this mode.
    zero_prune: super::p_zero::ZeroPrune,
}

/// Estimated resident bytes of one lowered input pair (key, rational,
/// map node) for [`RequestMemo`]'s bound.
const MEMO_BYTES_PER_PAIR: usize = 256;

/// vk5l WP4: the request memo of a p library. A request met again (same
/// family, powers, external, canonical numerator and epsilon symbol) reuses
/// its split, lowered and folded inputs instead of redoing them (on the beta
/// team's heaviest four-loop request the split and the lowering are 91 % of
/// the call); the sum over the batch, the reduction (with the step cache)
/// and the rendering follow as for a new request. Least recently used entries go first beyond the bound
/// of [`crate::OrderedCacheSettings::max_resident_bytes`] (`None`: no bound);
/// `VAKINT_P_REQUEST_MEMO=off` disables it.
#[derive(Default)]
pub(crate) struct RequestMemo {
    entries: HashMap<RequestKey, (Arc<FoldedRequest>, u64)>,
    clock: u64,
    pairs: usize,
    pub hits: u64,
    pub misses: u64,
}

impl RequestMemo {
    fn get(&mut self, key: &RequestKey) -> Option<Arc<FoldedRequest>> {
        self.clock += 1;
        let clock = self.clock;
        let found = self.entries.get_mut(key).map(|(value, used)| {
            *used = clock;
            Arc::clone(value)
        });
        if found.is_some() {
            self.hits += 1;
        } else {
            self.misses += 1;
        }
        found
    }

    fn insert(&mut self, key: RequestKey, value: Arc<FoldedRequest>, max_pairs: Option<usize>) {
        if max_pairs.is_some_and(|max| value.pairs > max) {
            return;
        }
        self.clock += 1;
        self.pairs += value.pairs;
        if let Some((old, _)) = self.entries.insert(key, (value, self.clock)) {
            self.pairs -= old.pairs;
        }
        while let Some(max) = max_pairs
            && self.pairs > max
        {
            let Some(oldest) = self
                .entries
                .iter()
                .min_by_key(|(_, (_, used))| *used)
                .map(|(key, _)| key.clone())
            else {
                break;
            };
            if let Some((old, _)) = self.entries.remove(&oldest) {
                self.pairs -= old.pairs;
            }
        }
    }
}

/// vk5l-physrun: `VAKINT_P_DIRECT_REQUESTS=off` keeps plain requests on the
/// ordinary split and lowering.
fn direct_requests_enabled() -> bool {
    static ENABLED: std::sync::LazyLock<bool> = std::sync::LazyLock::new(|| {
        !matches!(
            std::env::var("VAKINT_P_DIRECT_REQUESTS").as_deref(),
            Ok("off" | "0" | "false")
        )
    });
    *ENABLED
}

fn request_memo_enabled() -> bool {
    static ENABLED: std::sync::LazyLock<bool> = std::sync::LazyLock::new(|| {
        !matches!(
            std::env::var("VAKINT_P_REQUEST_MEMO").as_deref(),
            Ok("off" | "0" | "false")
        )
    });
    *ENABLED
}

/// Split, lower and fold one request through the library's [`RequestMemo`];
/// split and lowering (with the fold) seconds go to the statistics.
fn fold_request(
    handle: &PLibraryHandle,
    settings: &VakintSettings,
    request: &PropagatorRequest,
    forms: &mut BTreeMap<PresentationId, BTreeMap<PDot, AffineSlotForm>>,
    epsilon: Symbol,
) -> Result<Arc<FoldedRequest>, VakintError> {
    let _timing = crate::timing::lower();
    let key = request_memo_enabled().then(|| RequestKey {
        epsilon,
        family: request.family.name,
        powers: request.powers.to_vec(),
        external: request.external.clone(),
        numerator: request.canonical_numerator.to_owned(),
        zero_prune: super::p_zero::zero_prune(),
    });
    if let Some(key) = &key
        && let Some(found) = handle.request_memo().get(key)
    {
        with_stats(|stats| stats.memo_hits += 1);
        return Ok(found);
    }
    // vk5l-physrun (rule 17): every request first through the direct reader
    // (Horner over its loop-factor trie) and the Horner slot substitution; the
    // folded inputs are those of the split and lowering below, which run when
    // the direct lowering declines. VAKINT_P_DIRECT_REQUESTS=off skips it.
    if direct_requests_enabled() {
        let started = std::time::Instant::now();
        let mut direct_stats = super::p_lowering::DirectStats::default();
        let lowered =
            super::p_lowering::lower_request(settings, handle, request, &mut direct_stats)?;
        // Helper VH lane refpath: which requests the direct lowering takes
        // (`None`) or declines (the reason), for an active gate observer.
        super::p_gate::observe_direct(lowered.as_ref().err().map(|declined| declined.0));
        if let Ok(folded) = lowered {
            let split_channels = folded.len();
            let mut pairs = 0;
            let targets = folded
                .into_iter()
                .map(|(channel, inputs)| {
                    pairs += inputs
                        .by_basis
                        .values()
                        .map(|keys| keys.len())
                        .sum::<usize>();
                    (channel, inputs.by_basis)
                })
                .collect();
            with_stats(|stats| {
                stats.split_seconds += direct_stats.read_seconds + direct_stats.classify_seconds;
                stats.lower_seconds += started.elapsed().as_secs_f64()
                    - direct_stats.read_seconds
                    - direct_stats.classify_seconds;
            });
            let value = Arc::new(FoldedRequest {
                targets,
                split_channels,
                integral: integral_name(request),
                pairs,
            });
            if let Some(key) = key {
                let max_pairs = settings
                    .ordered_cache
                    .max_resident_bytes
                    .map(|bytes| bytes / MEMO_BYTES_PER_PAIR);
                handle
                    .request_memo()
                    .insert(key, Arc::clone(&value), max_pairs);
            }
            return Ok(value);
        }
    }
    let splitting = std::time::Instant::now();
    let mut pending = PendingChannels::new();
    split_request(handle, request, forms, &mut pending)?;
    with_stats(|stats| stats.split_seconds += splitting.elapsed().as_secs_f64());
    let lowering = std::time::Instant::now();
    let split_channels = pending.len();
    // Fold before lowering (vk5l WP4): every channel's coefficient is
    // `Σ value · basis(ε) · monomial` ([`spectator_terms`]); each loop-dot
    // monomial of its parts gets the weight `value · fraction` in its folded
    // target and basis. Each distinct monomial of the request is then lowered
    // once ([`lower_unit`]) and added, times its weights, into packed
    // accumulators per target and basis; nothing is kept per monomial. The
    // channels of a request share their monomials (six-loop census request:
    // 7 channels, one input monomial whose routing gives thousands of
    // loop-dot monomials of up to 98,234 keys), and rational sums are exact:
    // the folded inputs are those of lowering every channel.
    type MonomialKey = (usize, Vec<(PDot, i64)>);
    let mut sectors = Vec::<(PTop, &[i64])>::new();
    let mut sector_ids = HashMap::<(PresentationId, Vec<i64>), usize>::new();
    let mut weights = HashMap::<MonomialKey, Vec<(usize, Rational)>>::new();
    let mut slots = 0usize;
    let mut targets = BTreeMap::<Channel, BTreeMap<Atom, usize>>::new();
    let mut spectators = HashMap::<Atom, Vec<(Atom, Atom, Rational)>>::new();
    for (channel, ChannelParts { parts, .. }) in pending {
        let terms = spectators
            .entry(channel.coefficient.clone())
            .or_insert_with(|| spectator_terms(&channel.coefficient, epsilon))
            .clone();
        let mut monomials = Vec::<(MonomialKey, Rational)>::new();
        for PendingPart {
            top,
            powers,
            polynomial,
        } in parts
        {
            let next = sectors.len();
            let sector = *sector_ids
                .entry((top.presentation, powers.to_vec()))
                .or_insert_with(|| {
                    sectors.push((top, powers));
                    next
                });
            for (value, monomial) in polynomial {
                let mut dots = monomial.0;
                dots.sort();
                monomials.push(((sector, dots), fraction(value)));
            }
        }
        if monomials.is_empty() {
            continue;
        }
        for (monomial, basis, value) in terms {
            let target = targets
                .entry(Channel {
                    loops: channel.loops,
                    head: channel.head.clone(),
                    scale_power: channel.scale_power,
                    coefficient: monomial,
                })
                .or_default();
            let slot = *target.entry(basis).or_insert_with(|| {
                slots += 1;
                slots - 1
            });
            let unit = value == Rational::one();
            for (key, fraction) in &monomials {
                let weight = if unit {
                    fraction.clone()
                } else {
                    fraction * &value
                };
                weights.entry(key.clone()).or_default().push((slot, weight));
            }
        }
    }
    // Lower each distinct monomial once into the accumulators of its slots.
    let mut packed = vec![HashMap::<(usize, u128), Rational>::new(); slots];
    let mut unpacked = vec![HashMap::<(usize, Vec<u32>), Rational>::new(); slots];
    for ((sector, dots), mut by_slot) in weights {
        // Weights of one slot summed first.
        by_slot.sort_by_key(|(slot, _)| *slot);
        let mut merged = Vec::<(usize, Rational)>::new();
        for (slot, weight) in by_slot {
            match merged.last_mut() {
                Some((last, sum)) if *last == slot => *sum += &weight,
                _ => merged.push((slot, weight)),
            }
        }
        merged.retain(|(_, weight)| !weight.is_zero());
        if merged.is_empty() {
            continue;
        }
        let (top, _) = sectors[sector];
        let forms = forms.entry(top.presentation).or_default();
        let presentation = handle.library.presentation(top.presentation);
        for &(dot, _) in &dots {
            if let std::collections::btree_map::Entry::Vacant(entry) = forms.entry(dot) {
                let (a, b) = dot_rows(dot, usize::from(presentation.loops));
                entry.insert(
                    presentation
                        .express_scalar_product(&a, &b)
                        .map_err(|error| invalid(format!("{dot:?}: {error}")))?,
                );
            }
        }
        match lower_unit(top, &dots, forms) {
            UnitLowering::Packed(terms) => {
                for (slot, weight) in &merged {
                    let accumulator = &mut packed[*slot];
                    for (monomial, coefficient) in &terms {
                        *accumulator
                            .entry((sector, *monomial))
                            .or_insert_with(Rational::zero) += &(coefficient * weight);
                    }
                }
            }
            UnitLowering::Unpacked(terms) => {
                for (slot, weight) in &merged {
                    let accumulator = &mut unpacked[*slot];
                    for (monomial, coefficient) in &terms {
                        *accumulator
                            .entry((sector, monomial.clone()))
                            .or_insert_with(Rational::zero) += &(coefficient * weight);
                    }
                }
            }
        }
    }
    let mask = (1u128 << PACKED_BITS) - 1;
    let input_key = |sector: usize, exponent: &dyn Fn(usize) -> u32| {
        let (top, powers) = sectors[sector];
        let key = (0..top.arity)
            .map(|slot| {
                let base = if slot < top.line_count {
                    powers[slot].max(0)
                } else {
                    0
                };
                base - i64::from(exponent(slot))
            })
            .collect::<Vec<_>>();
        (top.presentation, key)
    };
    let mut pairs = 0;
    let targets = targets
        .into_iter()
        .map(|(channel, by_basis)| {
            let by_basis = by_basis
                .into_iter()
                .map(|(basis, slot)| {
                    let mut inputs = LoweredInputs::new();
                    for ((sector, monomial), value) in std::mem::take(&mut packed[slot]) {
                        if value.is_zero() {
                            continue;
                        }
                        let entry = input_key(sector, &|s| {
                            ((monomial >> (PACKED_BITS * s as u32)) & mask) as u32
                        });
                        *inputs.entry(entry).or_insert_with(Rational::zero) += &value;
                    }
                    for ((sector, monomial), value) in std::mem::take(&mut unpacked[slot]) {
                        if value.is_zero() {
                            continue;
                        }
                        let entry = input_key(sector, &|s| monomial[s]);
                        *inputs.entry(entry).or_insert_with(Rational::zero) += &value;
                    }
                    inputs.retain(|_, value| !value.is_zero());
                    pairs += inputs.len();
                    (basis, inputs)
                })
                .collect();
            (channel, by_basis)
        })
        .collect();
    with_stats(|stats| stats.lower_seconds += lowering.elapsed().as_secs_f64());
    let value = Arc::new(FoldedRequest {
        targets,
        split_channels,
        integral: integral_name(request),
        pairs,
    });
    if let Some(key) = key {
        let max_pairs = settings
            .ordered_cache
            .max_resident_bytes
            .map(|bytes| bytes / MEMO_BYTES_PER_PAIR);
        handle
            .request_memo()
            .insert(key, Arc::clone(&value), max_pairs);
    }
    Ok(value)
}

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
    let _timing = crate::timing::lower();
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
    let _timing = crate::timing::lower();
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
    let _timing = crate::timing::reduce();
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
    // Before the library is loaded or used: a boxed node binds the relation
    // cache roots on its first use (p_relations.rs).
    let package_relations = vakint.p_library()?.relations_root.clone();
    super::p_relations::install(settings, package_relations.as_deref());
    let dispatch_before = super::p_relations::dispatch_counters();
    let values = match fold_mode() {
        FoldMode::Off => reduce_requests_unfolded(vakint, settings, requests, progress, factored),
        mode => reduce_requests_folded(vakint, settings, requests, progress, mode, factored),
    }?;
    let dispatch = super::p_relations::dispatch_counters().since(&dispatch_before);
    with_stats(|stats| stats.dispatch = dispatch);
    if super::p_gate::active() {
        super::p_gate::observe_stats(&last_stats());
    }
    if phases::enabled() {
        eprintln!("vakint-p-dispatch {dispatch:?}");
        let stats = last_stats();
        eprintln!(
            "vakint-p-reduce requests={} channels={} reductions={} pairs={} split={:.6} lower={:.6} \
             prepare={:.6} library={:.6} render={:.6} memo_hits={} fold={:?}",
            requests.len(),
            stats.channels,
            stats.reductions,
            stats.pairs,
            stats.split_seconds,
            stats.lower_seconds,
            stats.prepare_seconds,
            stats.library_seconds,
            stats.render_seconds,
            stats.memo_hits,
            fold_mode()
        );
    }
    Ok(values)
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
    /// vk5l-physrun: `VAKINT_P_CHANNEL_REPLAY=1` reduces the coefficient-folded
    /// channels of one loop count together (RustRed channel replay).
    pub(crate) fn channel_replay() -> bool {
        static ON: std::sync::LazyLock<bool> = std::sync::LazyLock::new(|| {
            std::env::var("VAKINT_P_CHANNEL_REPLAY").as_deref() == Ok("1")
        });
        *ON
    }

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
    /// Requests served by the request memo ([`RequestMemo`]).
    pub memo_hits: usize,
    /// vk5l-physrun level-1 checkpoints reused / written.
    pub checkpoints_resumed: usize,
    pub checkpoints_written: usize,
    /// vk5l-physrun channel replay decisions (one line per loop count).
    pub channel_replay: Vec<String>,
    /// Census: split parts, their loop-dot monomials, lowered (channel, key)
    /// entries, folded channels and distinct integral keys.
    pub parts: usize,
    pub monomials: usize,
    pub lowered: usize,
    pub folded: usize,
    pub keys: usize,
    /// How the keys were served, by kind (vk5l WP4, three kinds).
    pub dispatch: super::p_relations::PDispatchCounters,
}

impl PReduceStats {
    /// The census record of one reduction call.
    pub(crate) fn census_json(&self, requests: usize, seconds: f64) -> serde_json::Value {
        serde_json::json!({
            "requests": requests,
            "channels": self.channels,
            "parts": self.parts,
            "monomials": self.monomials,
            "lowered": self.lowered,
            "folded": self.folded,
            "keys": self.keys,
            "reductions": self.reductions,
            "pairs": self.pairs,
            "split_seconds": self.split_seconds,
            "lower_seconds": self.lower_seconds,
            "prepare_seconds": self.prepare_seconds,
            "library_seconds": self.library_seconds,
            "render_seconds": self.render_seconds,
            "seconds": seconds,
        })
    }
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
pub(crate) fn spectator_terms(coefficient: &Atom, epsilon: Symbol) -> Vec<(Atom, Atom, Rational)> {
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
#[derive(Clone, Debug, Default, PartialEq)]
pub(crate) struct FoldedInputs {
    pub(crate) by_basis: BTreeMap<Atom, BTreeMap<(PresentationId, Vec<i64>), Rational>>,
    pub(crate) integrals: BTreeSet<String>,
}

/// [`reduce_requests`] with the channels folded: every channel coefficient is
/// written as `Σ value · basis(ε) · monomial` ([`spectator_terms`]), and the
/// channel key keeps only the monomial. The value of a folded channel is
/// `Σ_basis basis(ε) · R(Σ value · inputs) · monomial`, by linearity the sum
/// of the values of the channels it folds. In [`FoldMode::Basis`] every basis
/// function is reduced apart with rational inputs and multiplied in at the
/// render level; in [`FoldMode::Coefficient`] the basis functions enter the
/// inputs as rational functions of `d` and each folded channel is reduced once.
/// vk5l-physrun: `VAKINT_P_KEY_DUMP=<file>` appends every reduced key of a
/// folded channel as `loops<TAB>presentation<TAB>powers` (whole-input key
/// statistics: distinct keys per closure and across closures).
/// The inputs of a coefficient-folded channel: per key, `Σ_basis value ·
/// basis(ε)` at `d` ([`epsilon_coefficient`]), the keys in ascending order.
///
/// vk5l-physrun: the rational values of a key are gathered first (one entry
/// per basis, in ascending basis order), and the exact rational function is
/// built once per distinct list of (basis, value) pairs, with the same
/// products and sums in the same order as one key at a time. Physical
/// closures repeat the same small rational values over millions of keys.
/// A lowered key (presentation, powers) and its coefficient.
type KeyCoefficients = Vec<((PresentationId, Vec<i64>), Coefficient)>;

fn coefficient_inputs(
    context: &CoefficientContext,
    by_basis: BTreeMap<Atom, BTreeMap<(PresentationId, Vec<i64>), Rational>>,
    epsilon: Symbol,
    epsilon_at_d: &Coefficient,
    exact: ExactAlgebraLimits,
) -> Result<KeyCoefficients, VakintError> {
    let native = RustRedEvaluationError::from_native;
    let mut weights = Vec::with_capacity(by_basis.len());
    let mut parts = ahash::AHashMap::<(PresentationId, Vec<i64>), Vec<(u32, Rational)>>::new();
    for (index, (basis, lowered)) in by_basis.into_iter().enumerate() {
        weights.push(epsilon_coefficient(
            context,
            basis.as_view(),
            epsilon,
            epsilon_at_d,
            exact,
        )?);
        let index =
            u32::try_from(index).map_err(|_| invalid("too many epsilon basis functions"))?;
        for (key, value) in lowered {
            if !value.is_zero() {
                parts.entry(key).or_default().push((index, value));
            }
        }
    }
    let mut built = ahash::AHashMap::<Vec<(u32, Rational)>, Coefficient>::new();
    let mut inputs = Vec::with_capacity(parts.len());
    for (key, list) in parts {
        let value = match built.get(&list) {
            Some(value) => value.clone(),
            None => {
                let mut sum: Option<Coefficient> = None;
                for (index, value) in &list {
                    let term = context
                        .try_mul(
                            &weights[*index as usize],
                            &coefficient(context, value, exact)?,
                            exact,
                        )
                        .map_err(native)?;
                    sum = Some(match sum {
                        Some(sum) => context.try_add(&sum, &term, exact).map_err(native)?,
                        None => term,
                    });
                }
                let value = sum.ok_or_else(|| invalid("a folded key without a nonzero value"))?;
                built.insert(list, value.clone());
                value
            }
        };
        inputs.push((key, value));
    }
    inputs.sort_unstable_by(|left, right| left.0.cmp(&right.0));
    Ok(inputs)
}

fn key_dump<'a>(loops: usize, keys: impl Iterator<Item = &'a (PresentationId, Vec<i64>)>) {
    static FILE: std::sync::LazyLock<Option<std::sync::Mutex<std::fs::File>>> =
        std::sync::LazyLock::new(|| {
            let path = std::env::var_os("VAKINT_P_KEY_DUMP")?;
            std::fs::OpenOptions::new()
                .create(true)
                .append(true)
                .open(path)
                .ok()
                .map(std::sync::Mutex::new)
        });
    let Some(file) = FILE.as_ref() else {
        return;
    };
    let mut text = String::new();
    for (presentation, powers) in keys {
        use std::fmt::Write as _;
        let _ = writeln!(text, "{loops}\t{}\t{powers:?}", presentation.0);
    }
    use std::io::Write as _;
    let _ = file.lock().unwrap().write_all(text.as_bytes());
}

pub(crate) fn reduce_requests_folded(
    vakint: &Vakint,
    settings: &VakintSettings,
    requests: &[PropagatorRequest],
    progress: Option<&super::monitor::CallerProgress>,
    mode: FoldMode,
    factored: bool,
) -> Result<BTreeMap<(usize, Atom), Vec<ReducedPropagatorValue>>, VakintError> {
    with_stats(|stats| *stats = PReduceStats::default());
    let handle = vakint.p_library()?;
    let folded = fold_requests(&handle, settings, requests)?;
    // Helper VH lane refpath: the gate observer (a no-op unless a caller on
    // this thread observes the folded inputs, `super::p_gate::observe`).
    super::p_gate::observe_folded(&handle, vk_symbol!(&settings.epsilon_symbol), &folded);
    reduce_folded_channels(vakint, settings, &handle, folded, progress, mode, factored)
}

/// The folded channel inputs of a batch of plain requests: what
/// [`reduce_folded_channels`] hands to RustRed (one library reduction per
/// folded channel). Every request is split, lowered and folded
/// ([`fold_request`], or taken from the request memo) and the requests'
/// folded inputs are summed. Split from [`reduce_requests_folded`] as a pure
/// refactor (helper VH lane refpath, reference path R2 of the numerator gate):
/// the statistics it records and the map it returns are those of the inline
/// fold it replaces.
pub(crate) fn fold_requests(
    handle: &PLibraryHandle,
    settings: &VakintSettings,
    requests: &[PropagatorRequest],
) -> Result<BTreeMap<Channel, FoldedInputs>, VakintError> {
    let _timing = crate::timing::lower();
    let epsilon = vk_symbol!(&settings.epsilon_symbol);
    let mut forms = BTreeMap::new();
    // Every request split, lowered and folded (or taken from the request
    // memo); the batch's folded inputs are the sums of the requests' (rational
    // sums are exact: the order of the requests does not change them).
    let folded_requests = requests
        .iter()
        .map(|request| fold_request(handle, settings, request, &mut forms, epsilon))
        .collect::<Result<Vec<_>, VakintError>>()?;
    with_stats(|stats| {
        stats.channels = folded_requests
            .iter()
            .map(|request| request.split_channels)
            .sum()
    });
    let summing = std::time::Instant::now();
    let mut folded = BTreeMap::<Channel, FoldedInputs>::new();
    for request in &folded_requests {
        for (channel, by_basis) in &request.targets {
            let target = folded.entry(channel.clone()).or_default();
            target.integrals.insert(request.integral.clone());
            for (basis, inputs) in by_basis {
                match target.by_basis.get_mut(basis) {
                    None => {
                        target.by_basis.insert(basis.clone(), inputs.clone());
                    }
                    Some(sum) => {
                        for (key, value) in inputs {
                            *sum.entry(key.clone()).or_insert_with(Rational::zero) += value;
                        }
                    }
                }
            }
        }
    }
    drop(folded_requests);
    with_stats(|stats| {
        stats.lower_seconds += summing.elapsed().as_secs_f64();
        stats.folded = folded.len();
    });
    Ok(folded)
}

/// The reductions of folded channels ([`reduce_requests_folded`] after the
/// fold): one library reduction per folded channel, rendered and returned by
/// loop count and external head. Shared with the direct lowering
/// ([`super::p_lowering`]), which builds the same folded channels without
/// requests.
pub(crate) fn reduce_folded_channels(
    vakint: &Vakint,
    settings: &VakintSettings,
    handle: &PLibraryHandle,
    folded: BTreeMap<Channel, FoldedInputs>,
    progress: Option<&super::monitor::CallerProgress>,
    mode: FoldMode,
    factored: bool,
) -> Result<BTreeMap<(usize, Atom), Vec<ReducedPropagatorValue>>, VakintError> {
    let _timing = crate::timing::reduce();
    let epsilon = vk_symbol!(&settings.epsilon_symbol);
    let reducer = ChannelReducer {
        handle,
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
    // vk5l-physrun channel replay (VAKINT_P_CHANNEL_REPLAY=1): the channels of
    // one loop count are reduced together by PLibrary::reduce_channels (one
    // forward pass over a random mix, one replay per channel; RustRed falls
    // back to one reduction per channel on any failed zero check).
    let checkpoints = Checkpoints::new(
        settings,
        format!(
            "{:?}\n{:?}\n{:?}\n{:?}\n{:?}\n{}",
            handle.cache,
            handle.relations_root,
            checkpoint_arithmetic(&reducer.options),
            reducer.options.uncovered,
            reducer.options.master_mode,
            settings.epsilon_symbol
        ),
    );
    if mode == FoldMode::Coefficient && (FoldMode::channel_replay() || checkpoints.is_some()) {
        let mut prepared = Vec::new();
        for (
            channel,
            FoldedInputs {
                by_basis,
                integrals,
            },
        ) in folded
        {
            let preparing = std::time::Instant::now();
            let inputs = coefficient_inputs(context, by_basis, epsilon, &epsilon_at_d, exact)?;
            with_stats(|stats| {
                stats.prepare_seconds += preparing.elapsed().as_secs_f64();
                stats.pairs += inputs.len();
            });
            key_dump(channel.loops, inputs.iter().map(|(key, _)| key));
            let inputs = inputs
                .into_iter()
                .filter(|(_, value)| !value.is_zero())
                .map(|((presentation, powers), value)| {
                    let key = IntegralKey::try_new(powers)
                        .map_err(|error| invalid(format!("integral key: {error}")))?;
                    Ok((presentation, key, value))
                })
                .collect::<Result<Vec<_>, RustRedEvaluationError>>()?;
            if !inputs.is_empty() {
                prepared.push((channel, inputs, integrals));
            }
        }
        // Level-1 checkpoints: a channel whose reduced value is stored is not reduced again.
        let keys: Vec<Option<String>> = prepared
            .iter()
            .map(|(channel, inputs, _)| {
                checkpoints.as_ref().map(|c| {
                    c.key(
                        "channel",
                        &mut std::iter::once(format!("{channel:?}")).chain(
                            inputs
                                .iter()
                                .map(|(p, k, v)| format!("{} {:?} {}", p.0, k.powers(), v)),
                        ),
                    )
                })
            })
            .collect();
        let mut resumed: Vec<Option<Atom>> = keys
            .iter()
            .map(|key| key.as_ref().and_then(|key| checkpoints.as_ref()?.load(key)))
            .collect();
        // Channel replay only pays when channels share their trace: record
        // together only the channels whose most complex input (dots + numerator
        // rank) reaches the class maximum; the others are reduced on their own
        // (their trace is a small part of the shared one). Measured predictor,
        // written to the log and the run record (stats.channel_replay).
        let complexity = |inputs: &[(PresentationId, IntegralKey, Coefficient)]| -> i64 {
            inputs
                .iter()
                .map(|(_, key, _)| {
                    key.powers()
                        .iter()
                        .map(|&p| {
                            if p > 1 {
                                p - 1
                            } else if p < 0 {
                                -p
                            } else {
                                0
                            }
                        })
                        .sum::<i64>()
                })
                .max()
                .unwrap_or(0)
        };
        let mut groups = BTreeMap::<(usize, bool), Vec<usize>>::new();
        let mut by_loops = BTreeMap::<usize, Vec<(usize, i64)>>::new();
        for (at, (channel, inputs, _)) in prepared.iter().enumerate() {
            if resumed[at].is_none() {
                by_loops
                    .entry(channel.loops)
                    .or_default()
                    .push((at, complexity(inputs)));
            }
        }
        for (loops, members) in by_loops {
            let top = members.iter().map(|(_, w)| *w).max().unwrap_or(0);
            let shared: Vec<usize> = members
                .iter()
                .filter(|(_, w)| *w == top)
                .map(|(at, _)| *at)
                .collect();
            let replay = FoldMode::channel_replay() && shared.len() >= 2;
            let note = format!(
                "loops {loops}: channel complexities {:?}, top {top}, recorded together {:?}",
                members.iter().map(|(_, w)| *w).collect::<Vec<_>>(),
                if replay { shared.clone() } else { Vec::new() }
            );
            eprintln!("vakint-channel-replay {note}");
            with_stats(|stats| stats.channel_replay.push(note));
            for (at, _) in members {
                let together = replay && shared.contains(&at);
                groups.entry((loops, together)).or_default().push(at);
            }
        }
        let mut reductions: Vec<Option<PReduction>> = (0..prepared.len()).map(|_| None).collect();
        let input_sizes: Vec<usize> = prepared.iter().map(|(_, inputs, _)| inputs.len()).collect();
        // VAKINT_P_CHANNEL_LIMIT=<k>: reduce at most k channels in this run
        // (smallest first), then stop with an error; checkpoints are kept.
        let channel_limit: Option<usize> = std::env::var("VAKINT_P_CHANNEL_LIMIT")
            .ok()
            .and_then(|v| v.parse().ok());
        let mut reduced_here = 0usize;
        let mut inputs_of: Vec<Option<Vec<(PresentationId, IntegralKey, Coefficient)>>> = prepared
            .iter_mut()
            .map(|(_, inputs, _)| Some(std::mem::take(inputs)))
            .collect();
        for ((loops, together), members) in groups {
            if members.len() == 1 || !together {
                let mut members = members;
                members.sort_by_key(|&at| input_sizes[at]);
                for &at in &members {
                    if channel_limit.is_some_and(|limit| reduced_here >= limit) {
                        return Err(invalid(format!(
                            "channel limit reached after {reduced_here} channels (checkpoints kept)"
                        ))
                        .into());
                    }
                    let inputs = inputs_of[at].take().unwrap_or_default();
                    let clock = std::time::Instant::now();
                    let reduction =
                        reduce_inputs(&reducer, &prepared[at].0, inputs, &prepared[at].2)?;
                    reduced_here += 1;
                    eprintln!(
                        "vakint-channel loops {loops} channel {at}: {} input keys reduced in {:.1} s (peak RSS {:.1} GB)",
                        input_sizes[at],
                        clock.elapsed().as_secs_f64(),
                        std::fs::read_to_string("/proc/self/status")
                            .ok()
                            .and_then(|t| t
                                .lines()
                                .find(|l| l.starts_with("VmHWM:"))
                                .and_then(|l| l.split_whitespace().nth(1)?.parse::<f64>().ok()))
                            .unwrap_or(0.0)
                            / 1048576.0
                    );
                    if let (Some(c), Some(key)) = (&checkpoints, &keys[at]) {
                        let terms = render_terms(handle, &reduction, epsilon)?;
                        let reduced = collect_terms(terms);
                        c.store(key, &reduced);
                        resumed[at] = Some(reduced);
                    } else {
                        reductions[at] = Some(reduction);
                    }
                }
                continue;
            }
            if members.len() == 1 {
                let at = members[0];
                let inputs = inputs_of[at].take().unwrap_or_default();
                reductions[at] = Some(reduce_inputs(
                    &reducer,
                    &prepared[at].0,
                    inputs,
                    &prepared[at].2,
                )?);
                continue;
            }
            let integrals: BTreeSet<String> = members
                .iter()
                .flat_map(|&at| prepared[at].2.iter().cloned())
                .collect();
            let channels: Vec<_> = members
                .iter()
                .map(|&at| inputs_of[at].take().unwrap_or_default())
                .collect();
            let reducing = std::time::Instant::now();
            let reduced = {
                let _math = MathSpan::new();
                handle
                    .library
                    .reduce_channels(
                        channels,
                        p_reduction_request(
                            &reducer.options,
                            loops,
                            Vec::new(),
                            reducer.step_cache.clone(),
                        ),
                    )
                    .map_err(|error| reduction_error(error, &integrals))?
            };
            with_stats(|stats| {
                stats.library_seconds += reducing.elapsed().as_secs_f64();
                stats.reductions += members.len();
            });
            for (&at, reduction) in members.iter().zip(reduced) {
                super::p_relations::count_table_steps(
                    reduction.statistics.work.rule_applications() as u64,
                );
                if reduction.dimension != loops as i64 {
                    return Err(invalid(format!(
                        "a {loops}-loop request returned a reduction of dimension {}",
                        reduction.dimension
                    ))
                    .into());
                }
                reductions[at] = Some(reduction);
            }
        }
        for (at, ((channel, _, _), reduction)) in prepared.into_iter().zip(reductions).enumerate() {
            let reduced = if let Some(value) = resumed[at].take() {
                value
            } else {
                let Some(reduction) = reduction else { continue };
                let rendering = std::time::Instant::now();
                let terms = render_terms(handle, &reduction, epsilon)?;
                let reduced = collect_terms(terms);
                with_stats(|stats| stats.render_seconds += rendering.elapsed().as_secs_f64());
                if let (Some(c), Some(key)) = (&checkpoints, &keys[at]) {
                    c.store(key, &reduced);
                }
                reduced
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
        return Ok(values);
    }
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
                let inputs = coefficient_inputs(context, by_basis, epsilon, &epsilon_at_d, exact)?;
                with_stats(|stats| {
                    stats.prepare_seconds += preparing.elapsed().as_secs_f64();
                    stats.pairs += inputs.len();
                });
                key_dump(channel.loops, inputs.iter().map(|(key, _)| key));
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
                let terms = render_terms(handle, &reduction, epsilon)?;
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
                    for (monomial, value) in render_terms(handle, &reduction, epsilon)? {
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
        reduce_one(
            &reducer.handle.library,
            &reducer.options,
            channel.loops,
            inputs,
            reducer.step_cache.clone(),
            integrals,
            &format!(
                "{}L scale_power={} coefficient={}",
                channel.loops, channel.scale_power, channel.coefficient
            ),
        )?
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

/// The p master catalogue that finalizes `loops`-loop values of `vakint`;
/// `residual`: whether the values carry a residual terminal
/// `PResidual(node, k..)` ([`carries_residual`]). With a catalogue the
/// finalizer maps residual terminals onto its masters and applies the
/// catalogue guards (`p_integral::master_catalogue`); without one `Symbolic`
/// and `Numerical` refuse them with `PMasterCatalogueError::NoCatalogue` and
/// `Retain` returns them. The engine's catalogue
/// ([`Vakint::p_master_catalogue`]: the one set, `VAKINT_P_MASTER_CATALOGUE`,
/// the loaded package's, the built-in one) serves every five-loop value, and
/// every value of fewer loops that carries a residual terminal: a library
/// whose masters are residual terminals at one to four loops too (a 1-5
/// package, loaded, listed or injected) needs the catalogue's terminal rows,
/// which map them onto the Forcer masters (values from the Forcer table,
/// `PMasterSource::expansion`). A one- to four-loop value without one (every
/// value of a 1-4 library, whose basis layer maps onto `PMaster`) keeps the
/// Forcer table alone: there is nothing a catalogue could resolve.
fn p_master_catalogue(
    vakint: &Vakint,
    loops: usize,
    residual: bool,
) -> Option<Arc<PMasterCatalogue>> {
    if loops >= 5 || residual {
        vakint.p_master_catalogue()
    } else {
        None
    }
}

/// Whether `expression` contains a residual terminal `PResidual(node, k..)`.
fn carries_residual(expression: symbolica::atom::AtomView<'_>) -> bool {
    let mut found = false;
    expression.visitor(&mut |part| {
        if let symbolica::atom::AtomView::Fun(function) = part
            && function.get_symbol() == crate::p_integral::p_residual_symbol()
        {
            found = true;
        }
        !found
    });
    found
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
/// that loop count and values ([`p_master_catalogue`]) and the request's
/// vakint integrals in its messages.
pub(crate) fn finalize_reduced(
    vakint: &Vakint,
    settings: &VakintSettings,
    reduced: BTreeMap<(usize, Atom), Vec<ReducedPropagatorValue>>,
    names: &BTreeMap<(usize, Atom), BTreeSet<String>>,
    mode: RustRedMasterMode,
) -> Result<Atom, VakintError> {
    let _timing = crate::timing::reduce();
    let _step = crate::timing::finalize();
    let mut result = Atom::Zero;
    for ((loops, head), values) in reduced {
        let residual = values.iter().any(|value| {
            carries_residual(value.reduced.as_view())
                || value
                    .spectator
                    .as_ref()
                    .is_some_and(|spectator| carries_residual(spectator.as_view()))
        });
        let catalogue = p_master_catalogue(vakint, loops, residual);
        // lane vk5l-mastermap: the pole-bounded master basis reads its data
        // from the five-loop package (only when selected).
        let masters_directory = (loops >= 5
            && settings.p_master_basis == crate::p_integral::PMasterBasis::PoleBounded)
            .then(|| vakint.p_masters_directory())
            .flatten();
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
                PMasterSource::new(catalogue.as_deref(), &request)
                    .with_masters_directory(masters_directory.as_deref()),
            )
        })?;
        result += phases::time("p.final_simplify_real", || {
            simplify_real(finalized.as_view())
        });
    }
    Ok(result)
}

impl Vakint {
    /// Helper VH, lane l5-stage1 (2026-11-18): the value in the master mode of
    /// `settings` (its RustRed option; `Symbolic` without one) of a p-integral
    /// sum from its finalized `Retain` value `retained` (`loops` loops), with
    /// this engine's p master catalogue and WITHOUT a reduction: the
    /// finalisation of [`finalize_reduced`] applied to the retained master-basis
    /// expression ([`finalize::refinalise_with`]). Use: bring the reductions of
    /// an earlier binary with the same reduction code under this binary's
    /// finaliser; the caller proves the result equal to a full evaluation.
    pub fn refinalise_p_retained(
        &self,
        settings: &VakintSettings,
        retained: symbolica::atom::AtomView,
        loops: usize,
    ) -> Result<Atom, VakintError> {
        let mode = settings
            .evaluation_order
            .0
            .iter()
            .find_map(|method| match method {
                EvaluationMethod::RustRed(options) => Some(options.master_mode),
                _ => None,
            })
            .unwrap_or(RustRedMasterMode::Symbolic);
        let catalogue = p_master_catalogue(self, loops, carries_residual(retained));
        let masters_directory = (loops >= 5
            && settings.p_master_basis == crate::p_integral::PMasterBasis::PoleBounded)
            .then(|| self.p_masters_directory())
            .flatten();
        let request = format!("the re-finalisation of a Retain value of {loops}-loop p-integrals");
        let retained = retained.to_owned();
        let finalized = phases::time("p.refinalise", || {
            finalize::refinalise_with(
                &retained,
                loops,
                settings,
                mode,
                "RustRed",
                PMasterSource::new(catalogue.as_deref(), &request)
                    .with_masters_directory(masters_directory.as_deref()),
            )
        })?;
        Ok(simplify_real(finalized.as_view()))
    }

    /// The process-wide number of Laurent-expanded terms so far in which a
    /// tabulated master's unknown tail lies inside the needed orders
    /// (laurent-depth diagnostic of `p_integral::laurent`; such a term needs
    /// the exact products of 7c75215ea and is refused unless its unknown
    /// orders cancel).
    pub fn p_laurent_tail_in_range_terms() -> u64 {
        crate::p_integral::laurent::tail_in_range_terms()
    }
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
    let finalizing = std::time::Instant::now();
    let result = finalize_reduced(vakint, settings, reduced, &names, mode);
    super::p_gate::observe_finalize(finalizing.elapsed());
    result
}

#[cfg(test)]
#[path = "p_integral_tests.rs"]
mod tests;

/// vk5l-physrun level-1 checkpoints ([`crate::VakintSettings::p_checkpoint_dir`]).
/// Bump when the stored form or the meaning of a unit changes.
pub(crate) const P_CHECKPOINT_VERSION: &str = "vk5l-physrun.p-checkpoint.v1";

/// The arithmetic part of the checkpoint salt: the five-loop arithmetic with
/// its worker count fixed to one. The worker count changes how the replay is
/// scheduled, never a stored value, so a unit reduced on one thread is reused
/// on sixteen and vice versa. At one thread this is the arithmetic itself, so
/// the keys of checkpoints written before this change (all on one thread) are
/// unchanged.
fn checkpoint_arithmetic(options: &PRequestOptions) -> OrderedWeightArithmetic {
    options.arithmetic_for(5).with_threads(1)
}

/// One checkpoint directory with the salt every unit key includes (version,
/// library identity, relation root, arithmetic, uncovered policy, master
/// mode, epsilon symbol).
pub(crate) struct Checkpoints {
    dir: std::path::PathBuf,
    salt: String,
}

impl Checkpoints {
    pub(crate) fn new(settings: &VakintSettings, salt: String) -> Option<Self> {
        let dir = settings.p_checkpoint_dir.clone()?;
        std::fs::create_dir_all(&dir).ok()?;
        Some(Self {
            dir,
            salt: format!("{P_CHECKPOINT_VERSION}\n{salt}"),
        })
    }

    /// The key of a unit: blake3 of the salt, the unit kind and its parts.
    pub(crate) fn key(&self, kind: &str, parts: &mut dyn Iterator<Item = String>) -> String {
        let mut hasher = blake3::Hasher::new();
        hasher.update(self.salt.as_bytes());
        hasher.update(kind.as_bytes());
        for part in parts {
            hasher.update(b"\x1f");
            hasher.update(part.as_bytes());
        }
        format!("{kind}-{}", hasher.finalize().to_hex())
    }

    pub(crate) fn load(&self, key: &str) -> Option<Atom> {
        let text = std::fs::read_to_string(self.dir.join(format!("{key}.txt"))).ok()?;
        let value = crate::utils::vakint_macros::vk_parse!(&text).ok()?;
        with_stats(|stats| stats.checkpoints_resumed += 1);
        Some(value)
    }

    pub(crate) fn store(&self, key: &str, value: &Atom) {
        let path = self.dir.join(format!("{key}.txt"));
        let tmp = self.dir.join(format!(".{key}.tmp"));
        if std::fs::write(&tmp, value.to_canonical_string()).is_ok() {
            let _ = std::fs::rename(&tmp, &path);
            with_stats(|stats| stats.checkpoints_written += 1);
        }
    }
}
