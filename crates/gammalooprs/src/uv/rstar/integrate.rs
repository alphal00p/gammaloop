//! Integration of the remaining, Taylor-expanded single-mass vacuum graphs.
//!
//! The existing vacuum tensor projector acts before the marked line is cut.
//! Its external tensor coefficients accompany the scalar numerator through the
//! exact p-integral reduction and the final Laurent expansion.

use std::{
    collections::BTreeSet,
    time::{Duration, Instant},
};

use eyre::{Result, WrapErr, ensure, eyre};
use symbolica::{
    atom::{Atom, AtomCore, AtomView},
    function,
    id::Replacement,
    poly::series::SeriesDepth,
};
use vakint::{
    EvaluationMethod, RustRedMasterMode, TensorReductionMethod, Vakint, VakintExpression,
    VakintSettings,
    p_integral::{DirectClosureInput, MassiveClosure},
    symbols::S,
    vakint_parse,
};

use crate::uv::approx::integrated::Integrated;

/// A closure of today's path captured in remainder-boundary check mode
/// (DESIGN 7): its `VakintTerm.integral`, closing line, direct input (the
/// projected numerator, unrouted) and v1 record sequence.
#[derive(Clone, Debug)]
pub(super) struct CapturedClosure {
    pub(super) integral: Atom,
    pub(super) closure: MassiveClosure,
    pub(super) input: DirectClosureInput,
    /// The v1 record of this closure (the boundary record's `ref_closure`,
    /// DESIGN 6b).
    #[allow(dead_code)]
    pub(super) sequence: Option<usize>,
    /// False when the census dump skipped its evaluation.
    pub(super) evaluated: bool,
}

/// Everything captured from one [`integrate_single_mass`] call.
#[derive(Default)]
pub(super) struct Capture {
    pub(super) closures: Vec<CapturedClosure>,
    /// The value before `finish_tensor` (`result.expand()`).
    pub(super) value: Option<Atom>,
}

thread_local! {
    static CAPTURE: std::cell::RefCell<Option<Capture>> = const { std::cell::RefCell::new(None) };
}

/// Run `f` (today's outer integration of one remainder) capturing its
/// closures and value. Inert outside: one thread-local borrow per closure.
pub(super) fn capturing_closures<T>(f: impl FnOnce() -> T) -> (T, Capture) {
    // Restored on unwinding too: a capture left set would keep a clone of
    // every later closure input of the thread.
    let scope = super::ScopedSlot::enter(&CAPTURE, Capture::default());
    let result = f();
    (result, scope.exit().unwrap_or_default())
}

/// Whether a capture is active on this thread.
#[cfg(test)]
pub(super) fn capture_active() -> bool {
    CAPTURE.with(|slot| slot.borrow().is_some())
}

fn capture(make: impl FnOnce() -> CapturedClosure) {
    CAPTURE.with(|slot| {
        if let Some(capture) = slot.borrow_mut().as_mut() {
            capture.closures.push(make());
        }
    });
}

fn capture_value(value: &Atom) {
    CAPTURE.with(|slot| {
        if let Some(capture) = slot.borrow_mut().as_mut() {
            capture.value = Some(value.clone());
        }
    });
}

thread_local! {
    /// Closure calls and integrated terms of this thread's
    /// [`integrate_single_mass`] calls (counters for `IR_HARD_CHILD`).
    static SINGLE_MASS_COUNTS: std::cell::Cell<(usize, usize)> = const { std::cell::Cell::new((0, 0)) };
}

/// This thread's running totals `(closure calls, integrated terms)` of
/// [`integrate_single_mass`]; differences bracket one caller's integrations.
pub(crate) fn single_mass_counts() -> (usize, usize) {
    SINGLE_MASS_COUNTS.with(|cell| cell.get())
}

/// Integrate explicit vacuum topologies with exactly one massive denominator.
///
/// The mass and power must already be those of the selected surviving line.
/// Massless vacuum terms vanish; local terms (`integral = 1`) are retained.
/// Proper UV counterterms are supplied by the ordinary massive R operation.
/// The marked momentum needs a unit coefficient in the supplied loop basis;
/// the remaining massless graph must be covered by Vakint's p-integral backend.
/// The requested Laurent order applies after multiplying insertion coefficients:
/// their poles increase the depth needed from the remaining integral.
pub fn integrate_single_mass(
    vakint: &Vakint,
    settings: &VakintSettings,
    input: &VakintExpression,
) -> Result<Atom> {
    let started = Instant::now();
    let mut timing = SingleMassTiming::default();
    let mut tensor_settings = settings.clone();
    tensor_settings.tensor_reduction_method = TensorReductionMethod::FeynKit;
    tensor_settings.use_dot_product_notation = true;
    let mut result = Atom::Zero;
    let mut closures = Vec::new();
    let mut closure_bytes = 0_usize;
    // `GAMMALOOP_LOWERING=legacy` builds every closure input here, as before.
    let direct_lowering = !crate::uv::approx::streamed::legacy("direct_lowering");
    let mut direct = Vec::new();
    let mut direct_bytes = 0_usize;
    // The projector image of every distinct product of tensor factors.
    let mut tensor_images = std::collections::HashMap::<(Atom, Atom), Atom>::new();
    // Census: (summands, bytes) before and after the tensor projector and of
    // the closure inputs.
    let census = vakint::census::enabled();
    let mut census_source = (0_usize, 0_usize);
    let mut census_tensor = (0_usize, 0_usize);
    let mut census_closure = (0_usize, 0_usize, 0_usize);
    let mut census_topologies = BTreeSet::new();
    for source in &input.0 {
        if source.numerator.is_zero() {
            timing.zero_terms += 1;
            continue;
        }
        if source.integral == Atom::one() {
            let local_started = Instant::now();
            ensure!(
                !source.numerator.contains_symbol(S.k),
                "a local R* insertion cannot contain an unintegrated loop momentum"
            );
            result += &source.numerator;
            timing.local += local_started.elapsed();
            timing.local_terms += 1;
            continue;
        }
        let topology_started = Instant::now();
        let Some(cut) = cut_frame(&source.integral, &source.numerator)? else {
            timing.topology_cut += topology_started.elapsed();
            timing.scaleless_terms += 1;
            continue;
        };
        timing.topology_cut += topology_started.elapsed();
        timing.remaining_loops_max = timing.remaining_loops_max.max(cut.loops);
        // A massive self-loop times a massless vacuum component is scaleless.
        if cut.scaleless {
            timing.scaleless_terms += 1;
            continue;
        }
        let tensor_started = Instant::now();
        let mut term = source.clone();
        let factored = if direct_lowering && !crate::uv::approx::streamed::legacy("tensor_factors")
        {
            tensor_reduce_by_factors(vakint, &tensor_settings, source, &mut tensor_images)?
        } else {
            None
        };
        if let Some(factored) = factored {
            if crate::uv::approx::streamed::check("tensor_factors") {
                let mut whole = source.clone();
                whole
                    .tensor_reduce(vakint, &tensor_settings)
                    .wrap_err("while applying the existing vacuum tensor projector")?;
                let whole = Integrated::restore_lorentz_dimensions(&whole.numerator);
                let equal = (&factored - &whole).expand().is_zero();
                eprintln!(
                    "STREAM_CHECK stage=tensor_factors bytes={} equal={equal}",
                    source.numerator.as_view().get_byte_size()
                );
            }
            term.numerator = factored;
        } else {
            term.tensor_reduce(vakint, &tensor_settings).wrap_err(
                "while applying the existing vacuum tensor projector to an R* remainder",
            )?;
            // Lorentz representation metadata is constant during the following
            // Laurent expansion; only scalar evanescent coefficients may expand.
            term.numerator = Integrated::restore_lorentz_dimensions(&term.numerator);
        }
        crate::uv::approx::streamed::dump("closure_numerator", term.numerator.as_view());
        crate::uv::approx::streamed::dump("closure_topology", cut.integrand.as_view());
        crate::uv::approx::streamed::dump(
            "closure_routing",
            Atom::mul_many(
                cut.routing_pairs
                    .iter()
                    .map(|(head, image)| function!(symbolica::symbol!("rr"), head, image))
                    .collect::<Vec<_>>()
                    .iter()
                    .map(|a| a.as_view()),
            )
            .as_view(),
        );
        if census {
            census_source.0 += vakint::census::summands(source.numerator.as_view());
            census_source.1 += source.numerator.as_view().get_byte_size();
            census_tensor.0 += vakint::census::summands(term.numerator.as_view());
            census_tensor.1 += term.numerator.as_view().get_byte_size();
            census_topologies.insert(cut.integrand.clone());
        }
        timing.tensor += tensor_started.elapsed();
        timing.integrated_terms += 1;
        if settings.verify_numerator_identification {
            ensure!(
                super::pint6_dump::settings().is_none(),
                "PINT6_CLOSURE_DUMP does not cover the strict split branch"
            );
            // Strict matching admits only authenticated scalar symbols. Keep
            // arbitrary tensor spectators outside that check, as before.
            let strict_started = Instant::now();
            let integrated = cut.evaluate_split(vakint, settings, &term.numerator, &mut timing)?;
            // Inclusive branch time, overlapping the detailed phase totals.
            timing.strict_branch += strict_started.elapsed();
            let combine_started = Instant::now();
            result += integrated;
            timing.combine += combine_started.elapsed();
            continue;
        }
        // The p-integral owner separates homogeneous loop monomials from their
        // tensor spectators. Keeping their coefficients here lets its Laurent
        // expansion include insertion poles before truncating the full product.
        let preparation_started = Instant::now();
        let mut dumped = None;
        if let Some(dump) = super::pint6_dump::settings() {
            // Local export patch: before the direct-lowering branch, so both
            // closure paths are recorded (and skipped at >= the skip loops).
            let numerator = term.numerator.replace_multiple(&cut.routing).expand();
            let (p_loops, sequence) = super::pint6_dump::write(
                dump,
                settings,
                &cut.closure,
                &numerator,
                &cut.integrand,
                cut.loops,
            )?;
            dumped = Some((p_loops >= dump.skip_min_p_loops, sequence));
        }
        // Remainder-boundary check mode: after the record, before the skip.
        capture(|| CapturedClosure {
            integral: source.integral.clone(),
            closure: cut.closure.clone(),
            input: DirectClosureInput {
                topology: cut.integrand.clone(),
                routing: cut.routing_pairs.clone(),
                numerator: term.numerator.clone(),
            },
            sequence: dumped.map(|(_, sequence)| sequence),
            evaluated: !dumped.is_some_and(|(skipped, _)| skipped),
        });
        if dumped.is_some_and(|(skipped, _)| skipped) {
            // Export only: the request is recorded and not evaluated.
            timing.closure_preparation += preparation_started.elapsed();
            continue;
        }
        if direct_lowering && matches!(cut.integrand.as_view(), AtomView::Fun(_)) {
            // The routing and the expansion are applied by the p-integral
            // owner's direct lowering, onto the canonical family at once.
            let input_bytes = term.numerator.as_view().get_byte_size();
            timing.closure_preparation += preparation_started.elapsed();
            if census {
                census_closure.1 += input_bytes;
                census_closure.2 += 1;
            }
            if direct_bytes.saturating_add(input_bytes) > 16 << 20 {
                flush_direct(vakint, settings, &mut direct, &mut timing, &mut result)?;
                direct_bytes = 0;
            }
            direct_bytes = direct_bytes.saturating_add(input_bytes);
            direct.push((
                cut.closure,
                DirectClosureInput {
                    topology: cut.integrand,
                    routing: cut.routing_pairs,
                    numerator: term.numerator,
                },
            ));
            if direct.len() == 32 || direct_bytes >= 16 << 20 {
                flush_direct(vakint, settings, &mut direct, &mut timing, &mut result)?;
                direct_bytes = 0;
            }
            continue;
        }
        let preparation_started = Instant::now();
        let numerator = term.numerator.replace_multiple(&cut.routing).expand();
        let input = numerator * &cut.integrand;
        timing.closure_preparation += preparation_started.elapsed();
        let input_bytes = input.as_view().get_byte_size();
        if census {
            census_closure.0 += vakint::census::summands(numerator_census(&input));
            census_closure.1 += input_bytes;
            census_closure.2 += 1;
        }
        if closure_bytes.saturating_add(input_bytes) > 16 << 20 {
            flush_closures(vakint, settings, &mut closures, &mut timing, &mut result)?;
            closure_bytes = 0;
        }
        closure_bytes = closure_bytes.saturating_add(input_bytes);
        closures.push((cut.closure, input));
        if closures.len() == 32 || closure_bytes >= 16 << 20 {
            flush_closures(vakint, settings, &mut closures, &mut timing, &mut result)?;
            closure_bytes = 0;
        }
    }
    flush_closures(vakint, settings, &mut closures, &mut timing, &mut result)?;
    flush_direct(vakint, settings, &mut direct, &mut timing, &mut result)?;
    let combine_started = Instant::now();
    let result = result.expand();
    timing.combine += combine_started.elapsed();
    SINGLE_MASS_COUNTS.with(|cell| {
        let (closures, terms) = cell.get();
        cell.set((
            closures + timing.closure_calls,
            terms + timing.integrated_terms,
        ));
    });
    capture_value(&result);
    if census {
        vakint::census::record(
            "single_mass",
            serde_json::json!({
                "source_terms": input.0.len(),
                "integrated_terms": timing.integrated_terms,
                "source_summands": census_source.0,
                "source_bytes": census_source.1,
                "tensor_summands": census_tensor.0,
                "tensor_bytes": census_tensor.1,
                "closure_summands": census_closure.0,
                "closure_bytes": census_closure.1,
                "closures": census_closure.2,
                "distinct_cut_topologies": census_topologies.len(),
                "result_summands": vakint::census::summands(result.as_view()),
                "result_bytes": result.as_view().get_byte_size(),
                "tensor_seconds": timing.tensor.as_secs_f64(),
                "closure_preparation_seconds": timing.closure_preparation.as_secs_f64(),
                "closure_backend_seconds": timing.closure_backend.as_secs_f64(),
                "seconds": started.elapsed().as_secs_f64(),
            }),
        );
    }
    tracing::info!(
        uv = true,
        rstar = true,
        timing = true,
        source_terms = input.0.len(),
        zero_terms = timing.zero_terms,
        local_terms = timing.local_terms,
        scaleless_terms = timing.scaleless_terms,
        integrated_terms = timing.integrated_terms,
        remaining_loops_max = timing.remaining_loops_max,
        p_integral_loops_max = timing.remaining_loops_max.saturating_sub(1),
        closure_calls = timing.closure_calls,
        strict_split_terms = timing.strict_split_terms,
        strict_numerator_identification = settings.verify_numerator_identification,
        topology_cut_seconds = timing.topology_cut.as_secs_f64(),
        native_tensor_seconds = timing.tensor.as_secs_f64(),
        closure_preparation_seconds = timing.closure_preparation.as_secs_f64(),
        closure_backend_seconds = timing.closure_backend.as_secs_f64(),
        result_combination_seconds = timing.combine.as_secs_f64(),
        local_seconds = timing.local.as_secs_f64(),
        strict_branch_seconds = timing.strict_branch.as_secs_f64(),
        seconds = started.elapsed().as_secs_f64(),
        "Integrated single-mass R* remainder"
    );
    Ok(result)
}

/// Whether a numerator factor is a closed scalar of the vacuum tensor
/// projector: numbers, symbols and products, sums and powers of
/// `dot(k(i),k(j))`, `dot(p(i),p(j))` and functions without vectors.
fn closed_scalar(factor: AtomView<'_>) -> bool {
    match factor {
        AtomView::Num(_) | AtomView::Var(_) => true,
        AtomView::Add(sum) => sum.iter().all(closed_scalar),
        AtomView::Mul(product) => product.iter().all(closed_scalar),
        AtomView::Pow(power) => {
            let (base, exponent) = power.get_base_exp();
            closed_scalar(base) && closed_scalar(exponent)
        }
        AtomView::Fun(function) => {
            let symbol = function.get_symbol();
            if symbol == S.dot && function.get_nargs() == 2 {
                fn bare_vector(argument: AtomView<'_>) -> bool {
                    matches!(argument, AtomView::Fun(vector)
                        if (vector.get_symbol() == S.k || vector.get_symbol() == S.p)
                            && vector.get_nargs() == 1
                            && matches!(vector.get(0), AtomView::Num(_)))
                }
                // A vector or a rational linear combination of one kind of
                // vectors (loop or external): a sum of closed scalar products.
                let bare = |argument: AtomView<'_>| match argument {
                    AtomView::Add(sum) => {
                        let terms = sum.iter().collect::<Vec<_>>();
                        let single = |term: AtomView<'_>| match term {
                            AtomView::Mul(product) => {
                                let factors = product.iter().collect::<Vec<_>>();
                                matches!(factors.as_slice(),
                                    [AtomView::Num(_), v] | [v, AtomView::Num(_)] if bare_vector(*v))
                            }
                            other => bare_vector(other),
                        };
                        terms.iter().all(|t| single(*t))
                            && (terms.iter().all(|t| t.contains_symbol(S.k))
                                || terms.iter().all(|t| !t.contains_symbol(S.k)))
                    }
                    other => bare_vector(other),
                };
                let (a, b) = (function.get(0), function.get(1));
                // A loop momentum contracted with an external vector carries
                // integrated rank.
                return bare(a) && bare(b) && (a.contains_symbol(S.k) == b.contains_symbol(S.k));
            }
            !factor.contains_symbol(S.k)
                && !factor.contains_symbol(S.p)
                && !factor.contains_symbol(S.g)
                && !factor.contains_symbol(S.metric)
                && !factor.contains_symbol(S.dot)
                && function.iter().all(closed_scalar)
        }
    }
}

/// [`VakintTerm::tensor_reduce`] (with the dimension restoration) of a vacuum
/// term, applied once per distinct product of the factors of a summand that
/// are not closed scalars (open indices, loop momenta contracted with
/// external vectors, metrics). The projector is linear in the closed scalars
/// (`∫ T(k) S(k_i·k_j) = P[T] S`), so every summand `S·T` becomes `S·P[T]`;
/// the result is left unexpanded for the direct lowering.
pub(super) fn tensor_reduce_by_factors(
    vakint: &Vakint,
    settings: &VakintSettings,
    source: &vakint::VakintTerm,
    images: &mut std::collections::HashMap<(Atom, Atom), Atom>,
) -> Result<Option<Atom>> {
    let summands = match source.numerator.as_view() {
        AtomView::Add(sum) => sum.iter().collect::<Vec<_>>(),
        single => vec![single],
    };
    // Summands grouped by their tensor product: Σ_T (Σ S) · P[T], so that
    // every projector image appears once in the result.
    let mut order = Vec::<Atom>::new();
    let mut groups = std::collections::HashMap::<Atom, Vec<Atom>>::new();
    let mut scalar_terms = Vec::new();
    let mut any_tensor = false;
    for summand in summands {
        let factors = match summand {
            AtomView::Mul(product) => product.iter().collect::<Vec<_>>(),
            single => vec![single],
        };
        let (scalars, tensors): (Vec<_>, Vec<_>) = factors
            .into_iter()
            .partition(|factor| closed_scalar(*factor));
        if tensors.is_empty() {
            scalar_terms.push(summand.to_owned());
            continue;
        }
        any_tensor = true;
        let tensor = Atom::mul_many(tensors.iter().copied());
        let scalar = Atom::mul_many(scalars.iter().copied());
        match groups.get_mut(&tensor) {
            Some(list) => list.push(scalar),
            None => {
                order.push(tensor.clone());
                groups.insert(tensor, vec![scalar]);
            }
        }
    }
    if !any_tensor {
        // A numerator of closed scalars is its own projection.
        return Ok(Some(source.numerator.clone()));
    }
    let mut terms = scalar_terms;
    for tensor in order {
        let scalars = groups.remove(&tensor).unwrap_or_default();
        let key = (source.integral.clone(), tensor);
        let image = match images.get(&key) {
            Some(image) => image.clone(),
            None => {
                let mut term = vakint::VakintTerm {
                    integral: source.integral.clone(),
                    numerator: key.1.clone(),
                    vectors: source.vectors.clone(),
                };
                term.tensor_reduce(vakint, settings).wrap_err(
                    "while applying the existing vacuum tensor projector to an R* remainder",
                )?;
                let image = Integrated::restore_lorentz_dimensions(&term.numerator);
                images.insert(key, image.clone());
                image
            }
        };
        if image.is_zero() {
            continue;
        }
        let scalar = Atom::add_many(scalars.iter());
        if !scalar.is_zero() {
            terms.push(scalar * image);
        }
    }
    Ok(Some(Atom::add_many(terms.iter())))
}

/// The numerator sum of a closure input `numerator * topo(...)`.
fn numerator_census(input: &Atom) -> AtomView<'_> {
    if let AtomView::Mul(product) = input.as_view() {
        for factor in product.iter() {
            if matches!(factor, AtomView::Add(_)) {
                return factor;
            }
        }
    }
    input.as_view()
}

fn flush_direct(
    vakint: &Vakint,
    settings: &VakintSettings,
    closures: &mut Vec<(MassiveClosure, DirectClosureInput)>,
    timing: &mut SingleMassTiming,
    result: &mut Atom,
) -> Result<()> {
    if !closures.is_empty() {
        let backend_started = Instant::now();
        let integrated = MassiveClosure::evaluate_direct_batch(vakint, settings, closures)
            .wrap_err("while integrating an R* remainder through a massless p-integral")?;
        timing.closure_backend += backend_started.elapsed();
        timing.closure_calls += closures.len();
        let combine_started = Instant::now();
        for value in integrated {
            *result += value;
        }
        timing.combine += combine_started.elapsed();
        closures.clear();
    }
    Ok(())
}

fn flush_closures(
    vakint: &Vakint,
    settings: &VakintSettings,
    closures: &mut Vec<(MassiveClosure, Atom)>,
    timing: &mut SingleMassTiming,
    result: &mut Atom,
) -> Result<()> {
    if !closures.is_empty() {
        let backend_started = Instant::now();
        let integrated = MassiveClosure::evaluate_batch(vakint, settings, closures)
            .wrap_err("while integrating an R* remainder through a massless p-integral")?;
        timing.closure_backend += backend_started.elapsed();
        timing.closure_calls += closures.len();
        let combine_started = Instant::now();
        for value in integrated {
            *result += value;
        }
        timing.combine += combine_started.elapsed();
        closures.clear();
    }
    Ok(())
}

#[derive(Default)]
struct SingleMassTiming {
    topology_cut: Duration,
    tensor: Duration,
    closure_preparation: Duration,
    closure_backend: Duration,
    combine: Duration,
    local: Duration,
    strict_branch: Duration,
    zero_terms: usize,
    local_terms: usize,
    scaleless_terms: usize,
    integrated_terms: usize,
    remaining_loops_max: usize,
    closure_calls: usize,
    strict_split_terms: usize,
}

struct VacuumLine {
    id: Atom,
    ends: Atom,
    momentum: Atom,
    mass_squared: Atom,
    power: i64,
}

impl VacuumLine {
    fn parse(integral: &Atom) -> Result<Vec<Self>> {
        let AtomView::Fun(topology) = integral.as_view() else {
            return Err(eyre!("R* integration requires an explicit topo(prop(...))"));
        };
        ensure!(
            topology.get_symbol() == S.topo && topology.get_nargs() == 1,
            "R* integration requires an explicit topo(prop(...))"
        );
        let argument = topology.iter().next().unwrap();
        let factors = match argument {
            AtomView::Mul(product) => product.iter().collect::<Vec<_>>(),
            _ => vec![argument],
        };
        factors
            .into_iter()
            .map(|factor| {
                let AtomView::Fun(prop) = factor else {
                    return Err(eyre!("invalid R* propagator {factor}"));
                };
                ensure!(
                    prop.get_symbol() == S.prop && prop.get_nargs() == 5,
                    "invalid R* propagator {factor}"
                );
                let arguments = prop.iter().collect::<Vec<_>>();
                ensure!(
                    !arguments[2].contains_symbol(S.p),
                    "the R* remainder must be Taylor-expanded to a vacuum before integration"
                );
                Ok(Self {
                    id: arguments[0].to_owned(),
                    ends: arguments[1].to_owned(),
                    momentum: arguments[2].to_owned(),
                    mass_squared: arguments[3].to_owned(),
                    power: arguments[4]
                        .try_into()
                        .map_err(|_| eyre!("an R* propagator must have integer power"))?,
                })
            })
            .collect()
    }
}

/// The cut of a vacuum term at its unique massive line ([`cut_frame`]).
pub(super) struct VacuumCut {
    /// The closing line: the cut head `p(hard)`, the mass squared and power
    /// of the massive line.
    pub(super) closure: MassiveClosure,
    pub(super) routing: Vec<Replacement>,
    /// `routing` as (vacuum loop momentum `k(i)`, its image in the cut basis).
    pub(super) routing_pairs: Vec<(Atom, Atom)>,
    /// `topo(prop(..))` of the other lines in the cut momenta (or the tree
    /// factor of a one-loop vacuum term).
    pub(super) integrand: Atom,
    pub(super) scaleless: bool,
    /// Vacuum loops (the p-integral has one less).
    pub(super) loops: usize,
}

/// The cut of the vacuum term `integral` (`topo(prop(..))`) at its unique
/// massive line with positive power, with the closing head `p(n)` for the
/// smallest `n >= 1` absent from `probe` (the numerator); `None` when no line
/// is massive (the term is scaleless). Extracted from
/// [`integrate_single_mass`].
pub(super) fn cut_frame(integral: &Atom, probe: &Atom) -> Result<Option<VacuumCut>> {
    let lines = VacuumLine::parse(integral)?;
    let massive = lines
        .iter()
        .filter(|line| !line.mass_squared.is_zero() && line.power > 0)
        .collect::<Vec<_>>();
    if massive.is_empty() {
        return Ok(None);
    }
    ensure!(
        massive.len() == 1,
        "an R* remainder needs exactly one massive propagator, found {}",
        massive.len()
    );
    ensure!(
        lines
            .iter()
            .all(|line| { line.mass_squared.is_zero() || std::ptr::eq(line, massive[0]) }),
        "an R* remainder contains another massive propagator with nonpositive power"
    );
    Ok(Some(VacuumCut::new(&lines, massive[0], probe)?))
}

impl VacuumCut {
    fn evaluate_split(
        &self,
        vakint: &Vakint,
        settings: &VakintSettings,
        numerator: &Atom,
        timing: &mut SingleMassTiming,
    ) -> Result<Atom> {
        let preparation_started = Instant::now();
        let mut invariants = BTreeSet::new();
        numerator.visitor(&mut |part| {
            if let AtomView::Fun(function) = part
                && function.get_symbol() == S.dot
                && part.contains_symbol(S.k)
            {
                invariants.insert(part.to_owned());
                return false;
            }
            true
        });
        let coefficients =
            numerator.coefficient_list::<i16>(&invariants.into_iter().collect::<Vec<_>>());
        timing.closure_preparation += preparation_started.elapsed();
        timing.strict_split_terms += coefficients.len();
        let mut result = Atom::Zero;
        for (scalar, coefficient) in coefficients {
            let preparation_started = Instant::now();
            ensure!(
                !coefficient.contains_symbol(S.k),
                "the vacuum tensor projector left an uncontracted loop momentum: {coefficient}"
            );
            let scalar = scalar.replace_multiple(&self.routing).expand() * &self.integrand;
            timing.closure_preparation += preparation_started.elapsed();
            let integrated = self
                .evaluate(vakint, settings, &scalar, &coefficient, timing)
                .wrap_err("while integrating an R* remainder through a massless p-integral")?;
            let combine_started = Instant::now();
            result += integrated;
            timing.combine += combine_started.elapsed();
        }
        Ok(result)
    }

    fn evaluate(
        &self,
        vakint: &Vakint,
        settings: &VakintSettings,
        scalar: &Atom,
        coefficient: &Atom,
        timing: &mut SingleMassTiming,
    ) -> Result<Atom> {
        let preparation_started = Instant::now();
        let retained = settings
            .evaluation_order
            .0
            .iter()
            .find_map(|method| match method {
                EvaluationMethod::Forcer(options) => Some(options.master_mode),
                _ => None,
            })
            == Some(RustRedMasterMode::Retain);
        if retained {
            timing.closure_preparation += preparation_started.elapsed();
            let backend_started = Instant::now();
            let integrated = self.closure.evaluate(vakint, settings, scalar.as_view())?;
            timing.closure_backend += backend_started.elapsed();
            timing.closure_calls += 1;
            let combine_started = Instant::now();
            let result = coefficient * integrated;
            timing.combine += combine_started.elapsed();
            return Ok(result);
        }
        let epsilon = vakint_parse!(&settings.epsilon_symbol)
            .map_err(|error| eyre!(error))?
            .get_symbol()
            .ok_or_else(|| eyre!("Vakint epsilon must be a symbol"))?;
        let Some(valuation) = insertion_coefficient_valuation(coefficient, epsilon)? else {
            timing.closure_preparation += preparation_started.elapsed();
            return Ok(Atom::Zero);
        };
        // A pole in the already integrated insertion promotes higher orders of
        // the remaining loop integral. Keep tensor spectators outside Vakint,
        // but supply precisely the extra depth their Laurent coefficient needs.
        let mut closure_settings = settings.clone();
        closure_settings.number_of_terms_in_epsilon_expansion = settings
            .number_of_terms_in_epsilon_expansion
            .checked_sub(valuation.min(0))
            .ok_or_else(|| eyre!("the R* insertion requires an overflowing Laurent depth"))?;
        timing.closure_preparation += preparation_started.elapsed();
        let backend_started = Instant::now();
        let integrated = self
            .closure
            .evaluate(vakint, &closure_settings, scalar.as_view())?;
        timing.closure_backend += backend_started.elapsed();
        timing.closure_calls += 1;
        let combine_started = Instant::now();
        let target = settings.number_of_terms_in_epsilon_expansion - self.loops as i64 - 1;
        let result = (coefficient * integrated)
            .series(epsilon, 0, target)
            .wrap_err("could not truncate the complete R* insertion and integral product")?
            .to_atom();
        timing.combine += combine_started.elapsed();
        Ok(result)
    }

    fn new(lines: &[VacuumLine], massive: &VacuumLine, numerator: &Atom) -> Result<Self> {
        let mut heads = BTreeSet::new();
        for line in lines {
            line.momentum.visitor(&mut |part| {
                if let AtomView::Fun(function) = part
                    && function.get_symbol() == S.k
                    && function.get_nargs() == 1
                {
                    heads.insert(part.to_owned());
                    return false;
                }
                true
            });
        }
        let heads = heads.into_iter().collect::<Vec<_>>();
        for line in lines {
            for (head, coefficient) in line.momentum.coefficient_list::<i16>(&heads) {
                ensure!(
                    heads.contains(&head) && i64::try_from(coefficient.as_view()).is_ok(),
                    "R* vacuum routing requires integer linear loop momenta: {}",
                    line.momentum
                );
            }
        }
        let (pivot, sign) = heads
            .iter()
            .find_map(|head| {
                let coefficient = massive.momentum.coefficient(head.as_view());
                let sign = i64::try_from(coefficient.as_view()).ok()?;
                (sign == 1 || sign == -1).then_some((head, sign))
            })
            .ok_or_else(|| {
                eyre!(
                    "the marked R* momentum needs a unit loop coefficient for a measure-preserving cut: {}",
                    massive.momentum
                )
            })?;
        // The closing momentum is integrated, whereas every pre-existing p(i)
        // remains an external spectator, including indexed and compact forms.
        let mut external_ids = BTreeSet::new();
        numerator.visitor(&mut |part| {
            if let AtomView::Fun(function) = part
                && function.get_symbol() == S.p
                && let Some(index) = function.iter().next()
            {
                external_ids.insert(index.to_owned());
                return false;
            }
            true
        });
        let hard_index = (1_i64..)
            .find(|index| !external_ids.contains(&Atom::num(*index)))
            .expect("a finite numerator leaves unused momentum identifiers");
        let hard = function!(S.p, hard_index);
        let mut routing_pairs = heads
            .iter()
            .filter(|head| *head != pivot)
            .enumerate()
            .map(|(index, head)| (head.clone(), function!(S.k, index + 1)))
            .collect::<Vec<_>>();
        let mut routing = routing_pairs
            .iter()
            .map(|(head, image)| Replacement::new(head.to_pattern(), image.clone()))
            .collect::<Vec<_>>();
        let rest = (&massive.momentum - Atom::num(sign) * pivot)
            .expand()
            .replace_multiple(&routing);
        let pivot_image = Atom::num(sign) * (&hard - rest);
        routing.push(Replacement::new(pivot.to_pattern(), pivot_image.clone()));
        routing_pairs.push((pivot.clone(), pivot_image));
        let mut propagators = Atom::one();
        let mut radial = Atom::one();
        let mut has_loop = false;
        let mut has_hard = false;
        for line in lines.iter().filter(|line| !std::ptr::eq(*line, massive)) {
            let momentum = line.momentum.replace_multiple(&routing).expand();
            ensure!(
                !momentum.is_zero(),
                "R* vacuum routing produced a zero-momentum propagator"
            );
            has_loop |= momentum.contains_symbol(S.k);
            has_hard |= momentum.contains_symbol(S.p);
            propagators *= function!(S.prop, &line.id, &line.ends, &momentum, 0, line.power);
            // Only used if the cut is a tree, where every momentum is ±p.
            // Multi-loop cuts leave tree contraction to Vakint's frontend.
            if !momentum.contains_symbol(S.k) {
                let coefficient = momentum.coefficient(hard.as_view());
                let square = coefficient.pow(2) * function!(S.dot, &hard, &hard);
                radial *= square.pow(-line.power);
            }
        }
        Ok(Self {
            closure: MassiveClosure {
                momentum: hard,
                mass_squared: massive.mass_squared.clone(),
                power: massive.power,
            },
            routing,
            routing_pairs,
            integrand: if has_loop {
                function!(S.topo, propagators)
            } else {
                radial
            },
            scaleless: has_loop && !has_hard,
            loops: heads.len(),
        })
    }
}

/// The Laurent order of an R* insertion coefficient (a vacuum tensor
/// reduction coefficient times the integrated insertion), `None` when it
/// vanishes. IRH-DEPTH (regularity review of MATRIX-R45, finding 2): the
/// coefficient is a rational function of a symbolic `d`, and an
/// intermediate reduction coefficient can have a `1/(d-4)` (every rank-ten
/// orthogonal Weingarten coefficient has one), which a series in `ε` with
/// `d` symbolic does not see; the closure integral would then be requested
/// one order short. The order is therefore the lower of the order as
/// written and the order after `d = 4-2ε` (never above today's).
pub(super) fn insertion_coefficient_valuation(
    coefficient: &Atom,
    epsilon: symbolica::atom::Symbol,
) -> Result<Option<i64>> {
    let order = |expression: &Atom| -> Result<Option<i64>> {
        let series = expression
            .series(epsilon, 0, SeriesDepth::relative(1))
            .wrap_err("could not determine the Laurent order of an R* insertion coefficient")?;
        if series.is_zero() {
            return Ok(None);
        }
        let exponent = series.get_trailing_exponent();
        ensure!(
            exponent.is_integer(),
            "an R* insertion coefficient needs integer epsilon powers"
        );
        Ok(Some(exponent.numerator().to_i64().ok_or_else(|| {
            eyre!("the R* insertion's epsilon power is too large")
        })?))
    };
    let Some(written) = order(coefficient)? else {
        return Ok(None);
    };
    if epsilon != crate::utils::GS.dim_epsilon
        || !coefficient.contains_symbol(crate::utils::GS.dim)
    {
        return Ok(Some(written));
    }
    let regulated = Integrated::dimensionally_regularized(coefficient);
    // A series that cannot be taken after the substitution (an ε inside a
    // function argument) keeps today's order.
    Ok(Some(match order(&regulated) {
        Ok(Some(regulated)) => written.min(regulated),
        Ok(None) | Err(_) => written,
    }))
}

#[cfg(test)]
mod tests {
    use std::sync::LazyLock;

    use super::*;
    use symbolica::{coefficient::CoefficientView, prelude::RealLike};
    use vakint::{EvaluationOrder, ForcerOptions, LoopNormalizationFactor, VakintTerm};

    static VAKINT: LazyLock<Vakint> = LazyLock::new(|| Vakint::new().unwrap());

    fn parse(input: &str) -> Atom {
        vakint_parse!(input).unwrap()
    }

    fn settings() -> VakintSettings {
        VakintSettings {
            integral_normalization_factor: LoopNormalizationFactor::FMFTandMATAD,
            number_of_terms_in_epsilon_expansion: 4,
            n_cores: 1,
            evaluation_order: EvaluationOrder(vec![EvaluationMethod::Forcer(ForcerOptions {
                master_mode: RustRedMasterMode::Symbolic,
            })]),
            ..Default::default()
        }
    }

    fn integrate(input: &str) -> Atom {
        integrate_single_mass(
            &VAKINT,
            &settings(),
            &VakintExpression::try_from(parse(input)).unwrap(),
        )
        .unwrap()
    }

    fn assert_equal(actual: Atom, expected: Atom) {
        assert!(
            (&actual - &expected).together().is_zero(),
            "actual: {actual}\nexpected: {expected}"
        );
    }

    fn assert_numerically_equal(actual: Atom, expected: Atom) {
        let difference = (actual - expected).expand();
        let mut spectators = BTreeSet::new();
        difference.visitor(&mut |part| {
            if matches!(part, AtomView::Var(_) | AtomView::Fun(_)) {
                spectators.insert(part.to_owned());
                return false;
            }
            true
        });
        for (_, coefficient) in
            difference.coefficient_list::<i16>(&spectators.into_iter().collect::<Vec<_>>())
        {
            if coefficient.is_zero() {
                continue;
            }
            let AtomView::Num(number) = coefficient.as_view() else {
                panic!("non-numeric residual coefficient: {coefficient}")
            };
            let CoefficientView::Float(real, imaginary) = number.get_coeff_view() else {
                panic!("nonzero exact residual coefficient: {coefficient}")
            };
            // Default precision is thirty decimal digits. Summing before or
            // after numerical materialization can differ in the last digits.
            assert!(real.to_float().to_f64().abs() < 1e-24, "{difference}");
            assert!(imaginary.to_float().to_f64().abs() < 1e-24, "{difference}");
        }
    }

    // The former integration path, retained only as an independent ordering
    // oracle: integrate each loop monomial before multiplying its spectator,
    // explicitly promoting the integral depth for spectator poles.
    fn split_reference(input: &VakintExpression, settings: &VakintSettings) -> Result<Atom> {
        let mut tensor_settings = settings.clone();
        tensor_settings.tensor_reduction_method = TensorReductionMethod::FeynKit;
        tensor_settings.use_dot_product_notation = true;
        let mut result = Atom::Zero;
        for source in &input.0 {
            if source.numerator.is_zero() {
                continue;
            }
            if source.integral == Atom::one() {
                result += &source.numerator;
                continue;
            }
            let lines = VacuumLine::parse(&source.integral)?;
            let Some(massive) = lines
                .iter()
                .find(|line| !line.mass_squared.is_zero() && line.power > 0)
            else {
                continue;
            };
            let cut = VacuumCut::new(&lines, massive, &source.numerator)?;
            if cut.scaleless {
                continue;
            }
            let mut term = source.clone();
            term.tensor_reduce(&VAKINT, &tensor_settings)?;
            term.numerator = Integrated::restore_lorentz_dimensions(&term.numerator);
            result += cut.evaluate_split(
                &VAKINT,
                settings,
                &term.numerator,
                &mut SingleMassTiming::default(),
            )?;
        }
        Ok(result.expand())
    }

    #[test]
    fn complete_numerator_matches_split_integration_with_tensor_spectators_and_poles() {
        let topology = parse(
            "topo(prop(1,edge(1,2),k(1),0,1)*prop(2,edge(2,1),k(1)-k(2),0,1)*prop(3,edge(2,1),k(2),1,2))",
        );
        let polynomial = parse(
            "1/ε+dot(k(1),k(1))/(3-2*ε)+dot(k(1),k(2))/ε^2+dot(p(1),p(1))*dot(k(2),k(2))/(ε*(1-ε))",
        );
        let spectators = [
            Atom::one(),
            parse("g(71,72)"),
            parse("spenso::gamma(spenso::bis(4,a),spenso::bis(4,b),spenso::mink(4,mu))"),
            function!(idenso::color::CS.cas, 2, parse("spenso::coad(8)")).pow(2),
        ];
        for mode in [
            RustRedMasterMode::Symbolic,
            RustRedMasterMode::Retain,
            RustRedMasterMode::Numerical,
        ] {
            let settings = VakintSettings {
                number_of_terms_in_epsilon_expansion: 3,
                evaluation_order: EvaluationOrder(vec![EvaluationMethod::Forcer(ForcerOptions {
                    master_mode: mode,
                })]),
                ..settings()
            };
            for spectator in &spectators {
                let input =
                    VakintExpression::try_from(spectator * &polynomial * &topology).unwrap();
                let actual = integrate_single_mass(&VAKINT, &settings, &input).unwrap();
                assert!(!actual.contains_symbol(symbolica::atom::Symbol::DERIVATIVE));
                let expected = split_reference(&input, &settings).unwrap();
                if mode == RustRedMasterMode::Numerical {
                    assert_numerically_equal(actual, expected);
                } else {
                    assert_equal(actual, expected);
                }
            }
        }
        let mut cancelled = VakintExpression::try_from(topology).unwrap();
        cancelled.0[0].numerator = &polynomial - &polynomial;
        assert!(
            integrate_single_mass(&VAKINT, &settings(), &cancelled)
                .unwrap()
                .is_zero()
        );
    }

    #[test]
    fn cut_momentum_avoids_indexed_and_compact_external_identifiers() {
        let topology = parse("topo(prop(1,edge(1,1),k(1),1,2))");
        let lines = VacuumLine::parse(&topology).unwrap();
        let numerator = parse(
            "p(1,71)*p(2,spenso::mink(4,72))*spenso::dot(p(3,spenso::mink(4)),p(3,spenso::mink(4)))",
        );
        let cut = VacuumCut::new(&lines, &lines[0], &numerator).unwrap();
        assert_eq!(cut.closure.momentum, parse("p(4)"));
    }

    #[test]
    fn strict_numerator_verification_keeps_scalar_spectators_outside_matching() {
        let input = VakintExpression::try_from(parse(
            "c*(1/ε+dot(k(1),k(2))/(3-2*ε))*topo(prop(1,edge(1,2),k(1),0,1)*prop(2,edge(2,1),k(1)-k(2),0,1)*prop(3,edge(2,1),k(2),1,2))",
        ))
        .unwrap();
        let expected = integrate_single_mass(&VAKINT, &settings(), &input).unwrap();
        let strict = VakintSettings {
            verify_numerator_identification: true,
            ..settings()
        };
        assert_equal(
            integrate_single_mass(&VAKINT, &strict, &input).unwrap(),
            expected,
        );
    }

    #[test]
    fn one_loop_cut_keeps_the_tree_denominator_and_marked_power() {
        let actual = integrate("topo(prop(1,edge(1,2),k(7),mass2,2)*prop(2,edge(2,1),k(7),0,1))");
        let expected = MassiveClosure {
            momentum: parse("p(1)"),
            mass_squared: parse("mass2"),
            power: 2,
        }
        .evaluate(&VAKINT, &settings(), parse("1/dot(p(1),p(1))").as_view())
        .unwrap();
        assert_equal(actual, expected);
    }

    #[test]
    fn tensor_prepass_retains_external_coefficients_even_with_the_cut_momentum_id() {
        let actual = integrate("spectator*dot(k(3),p(1))^2*topo(prop(1,edge(1,1),k(3),1,3))");
        let expected = MassiveClosure {
            momentum: parse("p(1)"),
            mass_squared: Atom::one(),
            power: 3,
        }
        .evaluate(&VAKINT, &settings(), parse("dot(p(1),p(1))").as_view())
        .unwrap()
            * parse("spectator*dot(p(1),p(1))/(4-2*ε)");
        assert_equal(
            actual,
            expected
                .series(parse("ε").get_symbol().unwrap(), 0, 2)
                .unwrap()
                .to_atom(),
        );
    }

    #[test]
    fn insertion_poles_promote_the_required_higher_orders_of_the_closure() {
        for (terms, numerator, expected) in [
            (2, "1/ε", "1/ε^2+𝜋^2/12"),
            (1, "1/ε^2", "1/ε^3+𝜋^2/(12*ε)"),
            (1, "1/(ε^2*(4-2*ε))", "1/(4*ε^3)+1/(8*ε^2)+(1/16+𝜋^2/48)/ε"),
            (
                2,
                "dot(p(1),p(1))/(ε*(1-ε))",
                "dot(p(1),p(1))*(1/ε^2+1/ε+1+𝜋^2/12)",
            ),
        ] {
            let settings = VakintSettings {
                number_of_terms_in_epsilon_expansion: terms,
                ..settings()
            };
            let input = VakintExpression::try_from(
                parse(numerator) * parse("topo(prop(1,edge(1,1),k(1),1,2))"),
            )
            .unwrap();
            let actual = integrate_single_mass(&VAKINT, &settings, &input)
                .unwrap()
                .replace(parse("mursq").to_pattern())
                .with(Atom::one());
            assert_equal(actual, parse(expected));
        }
    }

    #[test]
    fn massless_remainders_vanish_but_local_insertions_survive() {
        let mut input =
            VakintExpression::try_from(parse("local_mass*topo(prop(1,edge(1,1),k(1),0,2))"))
                .unwrap();
        input.0.push(VakintTerm {
            integral: Atom::one(),
            numerator: parse("local_mass/ε"),
            vectors: Vec::new(),
        });
        assert_equal(
            integrate_single_mass(&VAKINT, &settings(), &input).unwrap(),
            parse("local_mass/ε"),
        );
    }

    #[test]
    fn a_marked_linear_combination_routes_with_unit_jacobian() {
        let canonical = parse(
            "topo(prop(1,edge(1,2),k(1),0,1)*prop(2,edge(2,1),k(1)-k(2),0,1)*prop(3,edge(2,1),k(2),1,2))",
        );
        let rerouted = canonical.replace_multiple(&[
            Replacement::new(parse("k(1)").to_pattern(), parse("k(5)")),
            Replacement::new(parse("k(2)").to_pattern(), parse("k(5)-k(8)")),
        ]);
        let evaluate = |input| {
            integrate_single_mass(
                &VAKINT,
                &settings(),
                &VakintExpression::try_from(input).unwrap(),
            )
            .unwrap()
        };
        let actual = evaluate(canonical);
        assert_equal(evaluate(rerouted), actual.clone());
        let epsilon = parse("ε").get_symbol().unwrap();
        let series = actual.series(epsilon, 0, 0).unwrap();
        let (_, leading_pole) = series.terms().find(|(power, _)| *power == -2).unwrap();
        assert_equal(leading_pole.clone(), parse("1/2"));
    }

    #[test]
    fn a_pure_hard_bridge_is_contracted_by_the_p_integral_frontend() {
        let actual = integrate(
            "topo(prop(1,edge(1,2),k(1),0,1)*prop(2,edge(2,1),k(1)-k(2),0,1)*prop(3,edge(2,3),k(2),1,2)*prop(4,edge(3,1),k(2),0,1))",
        );
        let expected = MassiveClosure {
            momentum: parse("p(1)"),
            mass_squared: Atom::one(),
            power: 2,
        }
        .evaluate(
            &VAKINT,
            &settings(),
            parse(
                "topo(prop(1,edge(1,2),k(1),0,1)*prop(2,edge(2,1),k(1)-p(1),0,1))/dot(p(1),p(1))",
            )
            .as_view(),
        )
        .unwrap();
        assert_equal(actual, expected);
    }

    #[test]
    fn two_massive_propagators_are_rejected() {
        let input = VakintExpression::try_from(parse(
            "topo(prop(1,edge(1,2),k(1),1,1)*prop(2,edge(2,1),k(1),1,1))",
        ))
        .unwrap();
        assert!(integrate_single_mass(&VAKINT, &settings(), &input).is_err());
    }
}
