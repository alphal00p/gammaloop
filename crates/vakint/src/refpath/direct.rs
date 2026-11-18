//! R3: the reference chain with a DIRECT HAND-OVER from idenso to the lowering
//! (helper VH, numerator assignment, lane idenso-audit). New entry point; R2
//! ([`super::run_graph`]) is unchanged.
//!
//! The tensor algebra stays idenso's (the same `contract_branch` chain as R2),
//! but
//! 1. the product of a branch is contracted in the EDGE basis: the edge
//!    momenta `Q(e)` are not routed before the contraction, so every momentum
//!    is one vector and idenso's output is a polynomial in `dot(Q(e), Q(f))`;
//! 2. the denominators are matched FIRST, on the graph's topology alone, with
//!    the two steps of R2's front end (`canonicalize` of the user's topology,
//!    then the canonical form's match), so the keys are R2's;
//! 3. the scalar is read by vakint's direct reader without expanding it
//!    (`p_lowering::lower_closure`: every `dot` of edge momenta becomes a
//!    bilinear form in the canonical scalar products through the edge
//!    images; zero-sector pruning as configured), and the folded channels go
//!    to R2's own reduction and finalisation.
//!
//! Anything this path does not reproduce exactly declines (`Ok(None)` from
//! [`lower_graph`]) and the caller keeps R2.
use std::collections::BTreeMap;
use std::time::Instant;

use eyre::Result;
use feynkit_graph::FeynmanDiagram;
use idenso::W_;
use spenso::network::tags::SPENSO_TAG;
use spenso::structure::representation::{Minkowski, RepName};
use symbolica::atom::{Atom, AtomCore, AtomView};

use super::{ColourBranch, Projection, SYMBOLS};
use crate::rustred_evaluation::p_gate::{GateView, ZeroSectorCount};
use crate::rustred_evaluation::p_handover::{
    self, Matched, PEdge, PLine, VectorRef, match_lines,
};
use crate::rustred_evaluation::p_integral::{Channel, FoldedInputs};
use crate::rustred_evaluation::p_lowering::DirectStats;
use crate::utils::vakint_macros::vk_symbol;
use crate::{Vakint, VakintSettings};

/// Seconds of R3's own stages (one graph, all projections).
#[derive(Clone, Debug, Default)]
pub struct DirectSeconds {
    /// idenso, edge basis (whole contraction incl. couplings and weight).
    pub contract: f64,
    /// Topology, two-step denominator match, edge images.
    pub matching: f64,
    /// Renaming the scalar into the reader's momenta, `D -> 4 - 2ε`.
    pub handover: f64,
    /// `lower_closure` (read, classify, slot substitution, fold).
    pub lower: f64,
    pub reduce: f64,
    pub finalise: f64,
}

/// One projection of R3.
pub struct DirectProjection {
    pub projection: &'static str,
    pub(crate) folded: BTreeMap<Channel, FoldedInputs>,
    pub value: Option<Atom>,
    pub(crate) stats: DirectStats,
}

/// What R3 produced for one graph, or why it declined.
pub enum DirectOutcome {
    /// Exact zero without reduction (scaleless after the match, zero numerator).
    Zero(&'static str),
    Lowered(Vec<DirectProjection>),
    Declined(String),
}

/// The colour terms of a graph for R3: R2's branches ([`super::colour_branches`])
/// or, with `merge`, the MERGED colour terms (lane share's A1, the beta
/// thread's merged records): the coefficient tensor of the multi-row locals is
/// rank-factorised exactly ([`merge_colour_channels`]) and every multi-row
/// local of a merged term carries `Σ_k w_k residual_k`. The sum over the terms
/// of `coefficient × Π residuals` is the sum over R2's branches (an identity of
/// the factorisation). Falls back to R2's branches when the merge refuses or
/// does not reduce the count. Returns the terms and R2's branch total.
pub fn colour_terms(diagram: &FeynmanDiagram, merge: bool) -> Result<(Vec<ColourBranch>, usize)> {
    use crate::p_integral::factorized::{
        ColourMode, colour_coefficient, merge_colour_channels, split_local_colour,
    };
    use feynkit_graph::NumeratorInterface;
    let (branches, total) = super::colour_branches(diagram)?;
    if !merge || branches.len() < 2 {
        return Ok((branches, total));
    }
    let mut rows = Vec::new();
    for (_, vertex) in diagram.vertices() {
        rows.push(split_local_colour(&NumeratorInterface::export(&vertex.numerator).0)?);
    }
    for (_, _, edge) in diagram.edges() {
        rows.push(split_local_colour(&NumeratorInterface::export(&edge.numerator).0)?);
    }
    let multi = rows
        .iter()
        .enumerate()
        .filter(|(_, local)| local.len() > 1)
        .map(|(position, _)| position)
        .collect::<Vec<_>>();
    if multi.is_empty() {
        return Ok((branches, total));
    }
    let dims = multi.iter().map(|p| rows[*p].len()).collect::<Vec<_>>();
    // Every channel of the multi-row locals (the single-row locals are common).
    let mut channels = Vec::new();
    let mut choice = vec![0usize; multi.len()];
    loop {
        let mut colour = Atom::num(1);
        for (position, local) in rows.iter().enumerate() {
            let k = multi.iter().position(|p| *p == position).map_or(0, |axis| choice[axis]);
            colour = colour * &local[k].colour;
        }
        let coefficient = colour_coefficient(&colour, ColourMode::Invariants)?;
        if !coefficient.is_zero() {
            channels.push((choice.clone(), coefficient));
        }
        let mut axis = 0;
        while axis < choice.len() {
            choice[axis] += 1;
            if choice[axis] < dims[axis] {
                break;
            }
            choice[axis] = 0;
            axis += 1;
        }
        if axis == choice.len() {
            break;
        }
    }
    let Some(merged) = merge_colour_channels(&dims, &channels) else {
        return Ok((branches, total));
    };
    if merged.len() >= branches.len() {
        return Ok((branches, total));
    }
    let adjoint = Atom::var(symbolica::symbol!("spenso::dA"));
    let fundamental = Atom::var(idenso::color::CS.nc);
    let mut out = Vec::with_capacity(merged.len());
    for term in merged {
        let coefficient = super::map_invariants(&term.coefficient, &adjoint, &fundamental)
            .expand()
            .together();
        if coefficient.is_zero() {
            continue;
        }
        let mut residuals = Vec::with_capacity(rows.len());
        for (position, local) in rows.iter().enumerate() {
            let Some(axis) = multi.iter().position(|p| *p == position) else {
                residuals.push(local[0].residual.clone());
                continue;
            };
            let mut sum = Atom::Zero;
            for (k, weight) in term.weights[axis].iter().enumerate() {
                if !weight.is_zero() {
                    sum += Atom::num(weight.clone()) * &local[k].residual;
                }
            }
            residuals.push(sum);
        }
        out.push(ColourBranch {
            coefficient,
            residuals,
        });
    }
    Ok((out, total))
}

/// R2's `contract_branch` (refpath.rs), with an option: `dots` asks idenso's
/// frontier for dot notation directly (`AlgebraContraction::Dots`) instead of
/// `Fully` followed by a separate `to_dots` rewrite. The scalar is the same
/// polynomial; its atom form may differ (the reader does not care).
fn contract_branch_edges(
    product: &Atom,
    momentum: &Atom,
    metric_only: bool,
    dots: bool,
) -> Result<(Atom, Atom)> {
    if !dots {
        return super::contract_branch(product, momentum, metric_only);
    }
    use idenso::dirac::GammaSimplifySettings;
    use idenso::tensor::{AlgebraContraction, AlgebraSettings, SymbolicTensor};
    use idenso::{CookMode, CookSettings, IndexTooling};
    use spenso::network::library::symbolic::ETS;
    use spenso::structure::abstract_index::AbstractIndex;
    use spenso::structure::representation::LibrarySlot;
    use spenso::structure::slot::IsAbstractSlot;
    let s = &*SYMBOLS;
    let cooking = CookSettings::indices()
        .with_mode(CookMode::ReversibleEncoding)
        .with_representation_payloads(true, true);
    let cooked = cooking.try_cook(product.as_view())?;
    let ports = cooked
        .list_dangling::<AbstractIndex>()?
        .iter()
        .map(|port| LibrarySlot::<AbstractIndex>::try_from(port.as_view()))
        .collect::<std::result::Result<Vec<_>, _>>()?;
    eyre::ensure!(
        ports.len() == 2
            && ports
                .iter()
                .all(|slot| slot.rep_name() == (Minkowski {}).into()),
        "a branch must have two dangling Lorentz ports, found {ports:?}"
    );
    let (left, right) = (ports[0].to_atom(), ports[1].to_atom());
    let metric = ETS.metric(&left, &right);
    let projector = if metric_only {
        metric
    } else {
        Atom::var(s.metric_tag) * metric
            + Atom::var(s.qq_tag) * super::indexed(momentum, &left)? * super::indexed(momentum, &right)?
    };
    let simplified = SymbolicTensor::infer_indexed(&cooked * &projector)?.simplify_algebra(
        &AlgebraSettings {
            gamma: Some(GammaSimplifySettings::canonical()),
            epsilon: true,
            contract: AlgebraContraction::Dots,
            ..Default::default()
        },
    )?;
    let scalar = cooking.uncook(simplified.expression().as_view());
    if metric_only {
        return Ok((scalar, Atom::Zero));
    }
    let g = scalar
        .replace(Atom::var(s.qq_tag))
        .with(Atom::Zero)
        .replace(Atom::var(s.metric_tag))
        .with(Atom::num(1));
    let qq = scalar
        .replace(Atom::var(s.metric_tag))
        .with(Atom::Zero)
        .replace(Atom::var(s.qq_tag))
        .with(Atom::num(1));
    Ok((g, qq))
}

/// R2's contraction of every branch, in the edge basis (no routing).
pub fn contract_edges(
    diagram: &FeynmanDiagram,
    branches: &[ColourBranch],
    momentum: &Atom,
    q_squared: &Atom,
    projection: Projection,
    dots: bool,
) -> Result<Vec<(&'static str, Atom)>> {
    let s = &*SYMBOLS;
    let d = Atom::var(s.d);
    let metric_only = projection == Projection::Metric;
    let mut g_sum = Atom::Zero;
    let mut qq_sum = Atom::Zero;
    for branch in branches {
        let product = branch
            .residuals
            .iter()
            .fold(Atom::num(1), |product, residual| product * residual);
        let continued = product
            .replace(Minkowski {}.to_symbolic([Atom::var(W_.d_), Atom::var(W_.a___)]))
            .with(Minkowski {}.to_symbolic([d.clone(), Atom::var(W_.a___)]));
        let (g, qq) = contract_branch_edges(&continued, momentum, metric_only, dots)?;
        g_sum += &branch.coefficient * g;
        qq_sum += &branch.coefficient * qq;
    }
    let couplings = |value: Atom| {
        diagram
            .model()
            .expand_couplings(&value)
            .replace(Atom::var(symbolica::symbol!("UFO::gs")))
            .with(Atom::num(1))
    };
    let i_q2 = Atom::i() * q_squared;
    let mut out = Vec::new();
    let transverse = || (&g_sum - &qq_sum / q_squared) / (&d - Atom::num(1)) / &i_q2;
    let longitudinal = || &qq_sum / q_squared / &i_q2;
    match projection {
        Projection::Metric => out.push(("metric", couplings(&g_sum / (&d - Atom::num(1)) / &i_q2))),
        Projection::Transverse => out.push(("transverse", couplings(transverse()))),
        Projection::Longitudinal => out.push(("longitudinal", couplings(longitudinal()))),
        Projection::Both => {
            out.push(("transverse", couplings(transverse())));
            out.push(("longitudinal", couplings(longitudinal())));
        }
    }
    Ok(out)
}

/// The lines and edges of a FeynKit two-point graph in the hand-over's
/// coordinates (loop momenta `k(i+1)` = FeynKit's `K(i)`, `p(1)` flowing in
/// through the first leg, as vakint's FeynKit ingress), and how the scalar's
/// vectors are identified: `Q(e, ..)` -> edge `e`, the independent external
/// momentum -> `orientation · p(1)`.
pub struct Lines {
    pub loops: usize,
    pub lines: Vec<PLine>,
    pub edges: Vec<PEdge>,
    orientation: i64,
    independent: (symbolica::atom::Symbol, Atom),
}

pub fn lines_of(
    diagram: &FeynmanDiagram,
    independent: &Atom,
) -> std::result::Result<Lines, String> {
    let basis = diagram.loop_momentum_basis();
    let loops = basis.loop_edges.len();
    let mut edges = diagram.edges().collect::<Vec<_>>();
    edges.sort_by_key(|(id, _, _)| *id);
    let (legs, lines): (Vec<_>, Vec<_>) = edges.into_iter().partition(|(_, ends, edge)| {
        edge.external.is_some() || ends.source.is_none() || ends.target.is_none()
    });
    if legs.len() != 2 || lines.is_empty() {
        return Err("not a two-point graph".into());
    }
    let independent_mask = basis
        .external_edges
        .iter()
        .map(|edge| i64::from(!basis.dependent_externals.contains(edge)))
        .collect::<Vec<_>>();
    let coefficients = |id: &feynkit_graph::EdgeId| -> Option<(Vec<i64>, i64)> {
        let (l, x) = basis.edge_signatures.get(id)?.integer_coefficients();
        let along = x
            .iter()
            .zip(&independent_mask)
            .map(|(c, m)| *c as i64 * m)
            .sum::<i64>();
        Some((l.iter().map(|c| *c as i64).collect(), along))
    };
    let (first, first_ends, _) = &legs[0];
    let (_, along) = coefficients(first).ok_or("a leg without routing")?;
    let into = if first_ends.source.is_none() { along } else { -along };
    let orientation = match into {
        1 => 1,
        -1 => -1,
        _ => return Err("the first leg does not carry ±P".into()),
    };
    let mut out_lines = Vec::new();
    for (id, ends, edge) in &lines {
        let (l, c) = coefficients(id).ok_or("a line without routing")?;
        let (Some(source), Some(target)) = (ends.source, ends.target) else {
            return Err("an internal line without two endpoints".into());
        };
        out_lines.push(PLine {
            id: id.0 as i64,
            source: source.0 as i64,
            target: target.0 as i64,
            loops: l,
            head: c * orientation,
            power: edge.propagator_power as i64,
        });
    }
    let mut out_edges = Vec::new();
    for (id, _, _) in legs.iter().chain(&lines) {
        let (l, c) = coefficients(id).ok_or("an edge without routing")?;
        out_edges.push(PEdge {
            id: id.0 as i64,
            loops: l,
            head: c * orientation,
        });
    }
    let independent = match independent.as_view() {
        AtomView::Fun(f) if f.get_nargs() >= 1 => {
            (f.get_symbol(), f.iter().next().unwrap().to_owned())
        }
        _ => return Err(format!("an unexpected external momentum {independent}")),
    };
    Ok(Lines {
        loops,
        lines: out_lines,
        edges: out_edges,
        orientation,
        independent,
    })
}

impl Lines {
    /// Identify a compact vector of idenso's output.
    fn vector(&self, view: AtomView<'_>) -> Option<VectorRef> {
        let AtomView::Fun(f) = view else {
            return None;
        };
        let first = f.iter().next()?;
        if f.get_symbol() == feynkit_graph::symbols::momentum() {
            return i64::try_from(first).ok().map(VectorRef::Edge);
        }
        (f.get_symbol() == self.independent.0 && first == self.independent.1.as_view())
            .then_some(VectorRef::Head(self.orientation))
    }
}

/// R3 on one graph: contraction (edge basis), match, hand-over, lowering and,
/// unless `reduce` is false, R2's reduction and finalisation per projection.
/// The numerators are those of R2 times the same weight.
#[allow(clippy::too_many_arguments)]
pub fn lower_graph(
    vakint: &Vakint,
    settings: &VakintSettings,
    diagram: &FeynmanDiagram,
    branches: &[ColourBranch],
    momentum: &Atom,
    q_squared: &Atom,
    weight: &Atom,
    projection: Projection,
    reduce: bool,
    dots: bool,
    seconds: &mut DirectSeconds,
) -> Result<DirectOutcome> {
    let started = Instant::now();
    let numerators = contract_edges(diagram, branches, momentum, q_squared, projection, dots)?;
    seconds.contract += started.elapsed().as_secs_f64();
    let started = Instant::now();
    let lines = match lines_of(diagram, momentum) {
        Ok(lines) => lines,
        Err(reason) => return Ok(DirectOutcome::Declined(reason)),
    };
    let matched = match match_lines(vakint, settings, lines.loops, &lines.lines, &lines.edges)? {
        Matched::Ready(matched) => matched,
        Matched::Scaleless => return Ok(DirectOutcome::Zero("scaleless after the match")),
        Matched::Declined(reason) => return Ok(DirectOutcome::Declined(reason)),
    };
    seconds.matching += started.elapsed().as_secs_f64();
    let epsilon = Atom::var(vk_symbol!(&settings.epsilon_symbol));
    let in_epsilon = Atom::num(4) - Atom::num(2) * &epsilon;
    let handle = vakint.p_library().map_err(crate::VakintError::from)?;
    let forbidden = [
        feynkit_graph::symbols::momentum(),
        feynkit_graph::symbols::loop_momentum(),
        feynkit_graph::symbols::external_momentum(),
    ];
    let mut out = Vec::new();
    for (name, numerator) in numerators {
        let numerator = numerator * weight;
        if numerator.is_zero() {
            continue;
        }
        let started = Instant::now();
        let scalar = match matched.rename(&numerator, SPENSO_TAG.dot, |v| lines.vector(v), &forbidden) {
            Ok(scalar) => scalar.replace(Atom::var(SYMBOLS.d)).with(in_epsilon.clone()),
            Err(reason) => return Ok(DirectOutcome::Declined(reason)),
        };
        seconds.handover += started.elapsed().as_secs_f64();
        let started = Instant::now();
        let mut stats = DirectStats::default();
        let folded = match matched.lower(settings, &handle, scalar.as_view(), &mut stats)? {
            Ok(folded) => folded,
            Err(declined) => return Ok(DirectOutcome::Declined(declined.0.into())),
        };
        seconds.lower += started.elapsed().as_secs_f64();
        let mut value = None;
        if reduce && !folded.is_empty() {
            let started = Instant::now();
            value = Some(matched.reduce(vakint, settings, &handle, folded.clone())?);
            seconds.reduce += started.elapsed().as_secs_f64();
        }
        out.push(DirectProjection {
            projection: name,
            folded,
            value,
            stats,
        });
    }
    if out.is_empty() {
        return Ok(DirectOutcome::Zero("zero numerator"));
    }
    Ok(DirectOutcome::Lowered(out))
}

/// The gate view of one R3 projection without its zero-sector keys (the form
/// of [`super::gate_view`]).
pub fn gate_view(
    vakint: &Vakint,
    settings: &VakintSettings,
    projection: &DirectProjection,
) -> Result<(GateView, ZeroSectorCount)> {
    let handle = vakint.p_library().map_err(crate::VakintError::from)?;
    Ok(p_handover::gate_view(settings, &handle, &projection.folded))
}

/// Exact equality of two Laurent expansions (`(a - b).expand() == 0`).
pub fn values_equal(a: &Atom, b: &Atom) -> bool {
    (a - b).expand().is_zero()
}
