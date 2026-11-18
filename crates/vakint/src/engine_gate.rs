//! The correctness gate and the paired timing harness of the numerator engine
//! ([`crate::rustred_evaluation::p_contract`]) on FeynKit DOT graphs of the
//! background-field `B -> B` self-energy (helper VH, numerator assignment,
//! lane engine).
//!
//! * [`run_graph`]: one graph through the engine: `from_dot`, the colour split
//!   of [`FactorizedNumerator::from_feynkit`] (one value per colour branch),
//!   per branch the prefactor in the reference path's conventions
//!   ([`refpath_prefactor`]: colour invariants named `CA, CF, TF, NA`,
//!   couplings expanded with `gs = 1`, the beta normalisation `1/(i q²)` and
//!   `1/(D - 1)`, `n_f` per closed quark loop), `contract_and_lower`, the
//!   branches summed, the gate view (zero-sector keys dropped), and optionally
//!   the reduction and finalisation through the reference's own entry
//!   (`reduce_folded_channels`).
//! * [`reference_branch_view`]: the reference path R2 (lane refpath's
//!   [`crate::refpath::contract`] → `from_diagram` → `Vakint::evaluate`,
//!   observed) on ONE branch of the engine's input, colour coefficient and
//!   graph weight left out on both sides: the per-branch gate.
//!
//! The graph-level gate compares [`run_graph`]'s view with lane refpath's
//! `gate_view` of the same graph and projection (live, [`crate::refpath::run_graph`],
//! or its stored dump), key by key and exactly ([`p_gate::compare_gate_views`]).

use std::sync::Arc;
use std::time::Instant;

use eyre::{Result, eyre};
use feynkit_graph::FeynmanDiagram;
use feynkit_kinematics::Kinematics;
use feynkit_model::Model;
use symbolica::atom::{Atom, AtomCore};
use symbolica::{function, symbol};

use crate::diagram_integral::DiagramIntegralOptions;
use crate::p_integral::factorized::{
    FactorizedNumerator, FromFeynkitOptions, Projector, ProjectorChoice, names,
};
use crate::refpath::{self, ColourBranch, Projection, SYMBOLS};
use crate::rustred_evaluation::p_contract::{self, ContractOptions, ContractedInputs};
use crate::rustred_evaluation::p_gate::{self, GateView, ZeroSectorCount};
use crate::symbols::S;
use crate::utils::vakint_macros::vk_symbol;
use crate::{Vakint, VakintExpression, VakintSettings};

/// The projections the engine runs one by one (`Both` = transverse and
/// longitudinal).
pub fn projections(projection: Projection) -> Vec<Projection> {
    match projection {
        Projection::Both => vec![Projection::Transverse, Projection::Longitudinal],
        other => vec![other],
    }
}

/// The engine's input values (one per non-zero colour branch) of a diagram
/// with the projector of `projection` (longitudinal: `Q(x,μ₀) Q(x,μ₁)` with
/// `p2_power = -1`).
pub fn branch_values(
    diagram: &FeynmanDiagram,
    model: &Model,
    projection: Projection,
) -> Result<Vec<FactorizedNumerator>> {
    let choice = match projection {
        Projection::Metric => ProjectorChoice::Metric,
        Projection::Transverse => ProjectorChoice::Transverse,
        Projection::Longitudinal => ProjectorChoice::Metric,
        Projection::Both => return Err(eyre!("one projection at a time")),
    };
    let options = FromFeynkitOptions {
        projector: choice,
        ..Default::default()
    };
    let mut values = FactorizedNumerator::from_feynkit(diagram, model, &options)?;
    if projection == Projection::Longitudinal {
        for value in &mut values {
            let [a, b] = value.ports.as_slice() else {
                return Err(eyre!("a two-point graph needs two ports"));
            };
            let Some(Some(x)) = value.external_edges.first() else {
                return Err(eyre!("no independent external edge"));
            };
            let q = symbol!(names::MOMENTUM);
            value.projector = Some(Projector {
                tensor: function!(q, *x as i64, a) * function!(q, *x as i64, b),
                p2_power: -1,
            });
        }
    }
    Ok(values)
}

/// `q² = dot(p(1), p(1))` in vakint's p notation (the head of R2).
fn q_squared() -> Atom {
    let p = function!(S.p, 1);
    function!(S.dot, &p, &p)
}

/// The beta normalisation of a projection in the value's dimension `D`:
/// metric and transverse `1/((D - 1) i q²)`, longitudinal `1/(i q²)`.
pub fn normalisation(value: &FactorizedNumerator, projection: Projection) -> Atom {
    let d = Atom::var(value.dimension);
    let i_q2 = Atom::i() * q_squared();
    match projection {
        Projection::Longitudinal => Atom::num(1) / i_q2,
        _ => Atom::num(1) / ((d - Atom::num(1)) * i_q2),
    }
}

/// Couplings expanded with the model's values and `gs = 1` (the reference's
/// convention).
pub fn expand_couplings(model: &Model, value: &Atom) -> Atom {
    model
        .expand_couplings(value)
        .replace(Atom::var(symbol!("UFO::gs")))
        .with(Atom::num(1))
}

/// The input type's colour coefficient (idenso invariants over `coad(dA)`,
/// `cof(Nc)`) in the reference's names (`CA`, `CF`, `TF`, `NA`, higher
/// invariants as plain symbols, `N_R = T_F N_A / C_F`).
pub fn refpath_colour(colour: &Atom) -> Atom {
    let s = &*SYMBOLS;
    let renamed = colour
        .replace(Atom::var(symbol!("spenso::dA")))
        .with(Atom::var(s.na))
        .replace(Atom::var(idenso::color::CS.nc))
        .with(Atom::var(s.nr));
    refpath::named_invariants(&renamed)
}

/// The whole prefactor of one branch in the reference's conventions: colour,
/// couplings (`gs = 1`), graph weight (symmetry factor and closed-loop sign),
/// `n_f` per closed quark loop, and the projection's normalisation.
pub fn refpath_prefactor(
    model: &Model,
    value: &FactorizedNumerator,
    projection: Projection,
) -> Atom {
    let s = &*SYMBOLS;
    let prefactor = &value.prefactor;
    let scalars = expand_couplings(model, &prefactor.scalars);
    &prefactor.overall
        * Atom::num(prefactor.loop_sign)
        * Atom::var(s.nf).pow(Atom::num(prefactor.nf_power as i64))
        * scalars
        * refpath_colour(&prefactor.colour)
        * normalisation(value, projection)
}

/// The prefactor of the per-branch gate: couplings (`gs = 1`) and the
/// projection's normalisation only (no colour, no graph weight).
pub fn branch_prefactor(
    model: &Model,
    value: &FactorizedNumerator,
    projection: Projection,
) -> Atom {
    expand_couplings(model, &value.prefactor.scalars) * normalisation(value, projection)
}

/// Seconds of one engine run of a graph (all projections).
#[derive(Clone, Debug, Default)]
pub struct EngineSeconds {
    /// `FeynmanDiagram::from_dot`.
    pub parse: f64,
    /// `from_feynkit` (colour split, value construction and validation).
    pub colour: f64,
    /// The engine: match on denominators; contraction + lowering
    /// (preparation, traces, contraction, fold); of which traces.
    pub matching: f64,
    pub lower: f64,
    pub trace: f64,
    /// The reduction (`reduce_folded_channels`: vakint prepare + library +
    /// render) and the finalisation.
    pub reduce: f64,
    pub prepare: f64,
    pub library: f64,
    pub finalise: f64,
}

/// One projection of an engine run.
pub struct EngineProjected {
    pub projection: &'static str,
    pub inputs: ContractedInputs,
    pub view: GateView,
    pub zero: ZeroSectorCount,
    /// The finalised value (when reduced).
    pub value: Option<Atom>,
    pub reductions: usize,
}

/// One graph through the engine.
pub struct EngineRun {
    pub name: String,
    pub loops: usize,
    pub branches: usize,
    pub zero_reason: Option<String>,
    pub seconds: EngineSeconds,
    pub results: Vec<EngineProjected>,
}

/// Options of [`run_graph`].
#[derive(Clone, Debug)]
pub struct EngineRunOptions {
    pub projection: Projection,
    pub contract: ContractOptions,
    /// Reduce and finalise (end-to-end value).
    pub reduce: bool,
    /// Laurent terms = loops + `extra_terms`.
    pub extra_terms: i64,
}

impl Default for EngineRunOptions {
    fn default() -> Self {
        Self {
            projection: Projection::Metric,
            contract: ContractOptions::default(),
            reduce: true,
            extra_terms: 1,
        }
    }
}

/// Run the engine on one DOT graph (module documentation).
pub fn run_graph(
    vakint: &Vakint,
    model: &Arc<Model>,
    dot: &str,
    options: &EngineRunOptions,
) -> Result<EngineRun> {
    let started = Instant::now();
    let diagram = FeynmanDiagram::from_dot(Arc::clone(model), dot)?;
    let mut seconds = EngineSeconds {
        parse: started.elapsed().as_secs_f64(),
        ..Default::default()
    };
    let loops = diagram.loop_count();
    let mut run = EngineRun {
        name: diagram.name().to_string(),
        loops,
        branches: 0,
        zero_reason: None,
        seconds: EngineSeconds::default(),
        results: Vec::new(),
    };
    if refpath::scaleless_contact(&diagram) {
        run.zero_reason = Some("no internal edge carries the external momentum: scaleless".into());
        run.seconds = seconds;
        return Ok(run);
    }
    let settings = refpath::vakint_settings(loops as i64 + options.extra_terms);
    for projection in projections(options.projection) {
        let started = Instant::now();
        let values = branch_values(&diagram, model, projection)?;
        seconds.colour += started.elapsed().as_secs_f64();
        run.branches = values.len();
        if values.is_empty() {
            run.zero_reason = Some("every colour branch vanishes".into());
            break;
        }
        let mut total = ContractedInputs::default();
        for value in &values {
            let prefactor = refpath_prefactor(model, value, projection);
            let inputs = p_contract::contract_and_lower(
                vakint,
                &settings,
                value,
                &prefactor,
                &options.contract,
            )?;
            total.absorb(inputs);
        }
        seconds.matching += total.stats.match_seconds;
        seconds.lower += total.stats.lower_seconds();
        seconds.trace += total.stats.trace_seconds;
        let (view, zero) = total.gate_view(vakint, &settings)?;
        let (value, reductions) = if options.reduce {
            let (value, reduce) = total.reduce(vakint, &settings)?;
            seconds.reduce += reduce.reduce;
            seconds.prepare += reduce.prepare;
            seconds.library += reduce.library;
            seconds.finalise += reduce.finalise;
            (Some(value), reduce.reductions)
        } else {
            (None, 0)
        };
        run.results.push(EngineProjected {
            projection: projection.name(),
            inputs: total,
            view,
            zero,
            value,
            reductions,
        });
    }
    run.seconds = seconds;
    Ok(run)
}

/// The reference path R2 on ONE value (one colour branch): lane refpath's
/// contraction ([`refpath::contract`]) of the branch's factors with the
/// projection, then `from_diagram` and `Vakint::evaluate` with the gate
/// observer; colour coefficient and graph weight are left out (the
/// counterpart of [`branch_prefactor`]). Returns the view without zero-sector
/// keys and the seconds of the contraction and of the evaluation.
pub fn reference_branch_view(
    vakint: &Vakint,
    diagram: &FeynmanDiagram,
    value: &FactorizedNumerator,
    projection: Projection,
    extra_terms: i64,
) -> Result<(GateView, ZeroSectorCount, f64, f64)> {
    let s = &*SYMBOLS;
    let mut residuals = vec![value.prefactor.scalars.clone()];
    residuals.extend(value.vertices.iter().map(|v| v.factor.clone()));
    residuals.extend(value.edges.iter().map(|e| e.factor.clone()));
    let branch = ColourBranch {
        coefficient: Atom::num(1),
        residuals,
    };
    let dimension = Atom::var(s.d);
    let kinematics = Kinematics::in_dimension(&dimension)?;
    let family = diagram.propagator_family(&kinematics)?;
    let [momentum] = family.external_momenta() else {
        return Err(eyre!("a two-point graph needs one external momentum"));
    };
    let momentum = momentum.clone();
    let q2 = Kinematics::in_dimension(&dimension)?
        .with_momenta([momentum.clone()])?
        .scalar_product(&momentum, &momentum)?;
    let started = Instant::now();
    let numerators = refpath::contract(diagram, &[branch], &momentum, &q2, projection)?;
    let contract_seconds = started.elapsed().as_secs_f64();
    let settings = refpath::vakint_settings(diagram.loop_count() as i64 + extra_terms);
    let epsilon = Atom::var(vk_symbol!(&settings.epsilon_symbol));
    let in_epsilon = Atom::num(4) - Atom::num(2) * &epsilon;
    let started = Instant::now();
    let mut view = GateView::new();
    let mut zero = ZeroSectorCount::default();
    for (_, numerator) in numerators {
        if numerator.is_zero() {
            continue;
        }
        let input = Atom::from(VakintExpression::from_diagram(
            diagram,
            &family,
            &numerator,
            &DiagramIntegralOptions::default(),
        )?);
        let input = input.replace(Atom::var(s.d)).with(in_epsilon.clone());
        let (result, observation) = p_gate::observe(|| vakint.evaluate(&settings, input.as_view()));
        let _value = result?;
        let (part, count) = observation.view_without_zero_sectors();
        view.extend(part);
        zero.total += count.total;
        zero.zero += count.zero;
        zero.unknown += count.unknown;
    }
    Ok((
        view,
        zero,
        contract_seconds,
        started.elapsed().as_secs_f64(),
    ))
}

/// The per-branch gate of a graph: for every branch of every projection, the
/// engine's view with [`branch_prefactor`] against [`reference_branch_view`].
/// Returns, per (projection, branch), the comparison.
pub fn branch_gate(
    vakint: &Vakint,
    model: &Arc<Model>,
    dot: &str,
    projection: Projection,
    contract: &ContractOptions,
    extra_terms: i64,
) -> Result<Vec<(&'static str, usize, p_gate::GateDifference, usize)>> {
    let diagram = FeynmanDiagram::from_dot(Arc::clone(model), dot)?;
    if refpath::scaleless_contact(&diagram) {
        return Ok(Vec::new());
    }
    let settings: VakintSettings =
        refpath::vakint_settings(diagram.loop_count() as i64 + extra_terms);
    let mut out = Vec::new();
    for projection in projections(projection) {
        for (index, value) in branch_values(&diagram, model, projection)?
            .iter()
            .enumerate()
        {
            let prefactor = branch_prefactor(model, value, projection);
            let inputs =
                p_contract::contract_and_lower(vakint, &settings, value, &prefactor, contract)?;
            let (engine, _) = inputs.gate_view(vakint, &settings)?;
            let (reference, _, _, _) =
                reference_branch_view(vakint, &diagram, value, projection, extra_terms)?;
            let keys = reference.len();
            out.push((
                projection.name(),
                index,
                p_gate::compare_gate_views(&reference, &engine),
                keys,
            ));
        }
    }
    Ok(out)
}
