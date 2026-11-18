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

/// Choices of R3's idenso call (defaults: edge basis, `Fully` + `to_dots`).
#[derive(Clone, Copy, Debug, Default)]
pub struct DirectOptions {
    /// `AlgebraContraction::Dots` instead of `Fully` followed by `to_dots`.
    pub dots: bool,
    /// Route the edge momenta into FeynKit's loop momenta BEFORE the
    /// contraction, exactly as R2 (idenso's input and output are then R2's);
    /// the hand-over then reads `dot(K(i), ..)` through the loop momenta's
    /// images instead of the edges'.
    pub routed: bool,
    /// Route with the spanning tree of the fewest loop + root momentum terms
    /// over the lines (exact enumeration; the beta thread's `compact-free`
    /// key) instead of FeynKit's stored basis. Implies routing before idenso.
    pub compact: bool,
}

/// Pseudo-edge id of FeynKit's loop momentum `K(i)` in the hand-over tables.
const LOOP_EDGE_BASE: i64 = 900_000;

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
    /// Choosing the compact basis (option `compact`), and its routing size
    /// against FeynKit's stored basis, with the number of spanning trees.
    pub routing: f64,
    pub routing_size: usize,
    pub stored_routing_size: usize,
    pub trees: usize,
}

/// One projection of R3.
pub struct DirectProjection {
    pub projection: &'static str,
    pub(crate) folded: BTreeMap<Channel, FoldedInputs>,
    pub value: Option<Atom>,
    pub(crate) stats: DirectStats,
    distinct_dots: usize,
}

/// Counters of one projection's lowering (the direct reader's statistics) and
/// of the handed-over scalar.
#[derive(Clone, Copy, Debug, Default)]
pub struct LoweringCounts {
    /// Terms of the scalar the reader visited.
    pub numerator_terms: usize,
    /// Distinct `vakint::dot` atoms in the handed-over scalar.
    pub distinct_dots: usize,
    /// Canonical terms after reading, after merging, distinct canonical loop
    /// monomials, lowered entries, folded entries.
    pub canonical_terms: usize,
    pub merged_terms: usize,
    pub distinct_monomials: usize,
    pub lowered: usize,
    pub folded: usize,
    /// Products formed by the slot substitution (a term times a slot form's
    /// term: the substitution's work).
    pub slot_products: usize,
    pub pruned_early: usize,
}

impl DirectProjection {
    pub fn counts(&self) -> LoweringCounts {
        LoweringCounts {
            numerator_terms: self.stats.numerator_terms,
            distinct_dots: self.distinct_dots,
            canonical_terms: self.stats.canonical_terms,
            merged_terms: self.stats.merged_terms,
            distinct_monomials: self.stats.distinct_monomials,
            lowered: self.stats.lowered,
            folded: self.stats.folded,
            slot_products: self.stats.slot_products,
            pruned_early: self.stats.pruned_early,
        }
    }
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

/// R2's `contract_branch` (refpath.rs): the same idenso calls and the same
/// scalar, without R2's final `list_dangling` check (the hand-over refuses
/// leftover notation itself). Option `dots` asks idenso's frontier for dot
/// notation directly (`AlgebraContraction::Dots`) instead of `Fully` followed
/// by a separate `to_dots` rewrite: the same polynomial, possibly another atom
/// form (the reader does not care).
fn contract_branch_edges(
    product: &Atom,
    momentum: &Atom,
    metric_only: bool,
    dots: bool,
) -> Result<(Atom, Atom)> {
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
            contract: if dots {
                AlgebraContraction::Dots
            } else {
                AlgebraContraction::Fully
            },
            ..Default::default()
        },
    )?;
    let reduced = if dots {
        simplified
    } else {
        simplified.to_dots()?
    };
    // No `list_dangling` check of the scalar here (R2 builds a network for
    // it): the hand-over's rename refuses any leftover tensor notation.
    let scalar = cooking.uncook(reduced.expression().as_view());
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
    options: DirectOptions,
    compact: Option<(&Lines, &CompactBasis)>,
) -> Result<Vec<(&'static str, Atom)>> {
    let s = &*SYMBOLS;
    let basis = diagram.loop_momentum_basis();
    let d = Atom::var(s.d);
    let metric_only = projection == Projection::Metric;
    let mut g_sum = Atom::Zero;
    let mut qq_sum = Atom::Zero;
    // The couplings are expanded (and gs = 1) in every residual BEFORE the
    // contraction, not on the contracted scalar as R2 does: a scalar
    // substitution commutes with the contraction, and the residuals are small.
    let couplings = |value: &Atom| {
        diagram
            .model()
            .expand_couplings(value)
            .replace(Atom::var(symbolica::symbol!("UFO::gs")))
            .with(Atom::num(1))
    };
    for branch in branches {
        let product = branch
            .residuals
            .iter()
            .fold(Atom::num(1), |product, residual| product * couplings(residual));
        let product = match compact {
            Some((lines, chosen)) => lines.route_compact(&product, chosen).map_err(|e| eyre::eyre!(e))?,
            None if options.routed => basis.route_expression(&product),
            None => product,
        };
        let continued = product
            .replace(Minkowski {}.to_symbolic([Atom::var(W_.d_), Atom::var(W_.a___)]))
            .with(Minkowski {}.to_symbolic([d.clone(), Atom::var(W_.a___)]));
        let (g, qq) = contract_branch_edges(&continued, momentum, metric_only, options.dots)?;
        g_sum += &branch.coefficient * g;
        qq_sum += &branch.coefficient * qq;
    }
    let couplings = |value: Atom| value;
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
    /// With a compact basis: the edge whose momentum is `K(i)`.
    loop_edges: Option<Vec<i64>>,
}

/// A loop-momentum basis chosen by [`Lines::compact_basis`]: the cotree edges
/// (their momenta are the new `K(0..L)`) and every edge's momentum in it.
pub struct CompactBasis {
    pub cotree: Vec<i64>,
    /// edge id -> (coefficients over the new `K(i)`, coefficient of FeynKit's
    /// independent external momentum).
    pub routes: std::collections::BTreeMap<i64, (Vec<i64>, i64)>,
    /// Sum over the lines of the loop + root momentum terms (the key), and the
    /// same for FeynKit's stored basis.
    pub size: usize,
    pub stored_size: usize,
    pub trees: usize,
}

/// `A x = b` over the rationals for a small integer `A` (columns given),
/// `None` if singular or the solution is not integral.
fn solve_integral(columns: &[Vec<i64>], b: &[i64]) -> Option<Vec<i64>> {
    let n = b.len();
    if columns.len() != n {
        return None;
    }
    // augmented matrix of rationals (num, den) as i128
    let mut m: Vec<Vec<(i128, i128)>> = (0..n)
        .map(|r| {
            let mut row: Vec<(i128, i128)> =
                columns.iter().map(|c| (c[r] as i128, 1i128)).collect();
            row.push((b[r] as i128, 1));
            row
        })
        .collect();
    fn norm((a, b): (i128, i128)) -> (i128, i128) {
        fn gcd(mut a: i128, mut b: i128) -> i128 {
            while b != 0 {
                (a, b) = (b, a % b);
            }
            a.abs()
        }
        let g = gcd(a, b).max(1);
        let (a, b) = (a / g, b / g);
        if b < 0 { (-a, -b) } else { (a, b) }
    }
    for col in 0..n {
        let pivot = (col..n).find(|&r| m[r][col].0 != 0)?;
        m.swap(col, pivot);
        let (pn, pd) = m[col][col];
        for k in col..=n {
            let (a, b) = m[col][k];
            m[col][k] = norm((a * pd, b * pn));
        }
        for r in 0..n {
            if r == col || m[r][col].0 == 0 {
                continue;
            }
            let (fn_, fd) = m[r][col];
            for k in col..=n {
                let (a, b) = m[r][k];
                let (c, d) = m[col][k];
                // a/b - (fn/fd)(c/d)
                m[r][k] = norm((a * fd * d - fn_ * c * b, b * fd * d));
            }
        }
    }
    (0..n)
        .map(|r| {
            let (a, b) = m[r][n];
            (b == 1 || a % b == 0).then_some((a / b) as i64)
        })
        .collect()
}

impl Lines {
    /// The spanning tree of the lines with the fewest loop + root momentum
    /// terms (by enumeration; ties: the smallest sorted cotree), and every
    /// edge's momentum in its basis. Exact: momenta are re-expressed, the
    /// images through the match are unchanged.
    pub fn compact_basis(&self) -> std::result::Result<CompactBasis, String> {
        let vertices = self
            .lines
            .iter()
            .flat_map(|l| [l.source, l.target])
            .collect::<std::collections::BTreeSet<_>>()
            .into_iter()
            .collect::<Vec<_>>();
        let index = |v: i64| vertices.binary_search(&v).unwrap();
        let n = vertices.len();
        if self.lines.len() + 1 < n || self.lines.len() + 1 - n != self.loops {
            return Err("the lines do not form a connected graph with the loop count".into());
        }
        let vector = |loops: &[i64], head: i64| -> (Vec<i64>, i64) {
            (loops.to_vec(), head * self.orientation)
        };
        let all_edges = self
            .edges
            .iter()
            .filter(|e| e.id < LOOP_EDGE_BASE)
            .map(|e| (e.id, vector(&e.loops, e.head)))
            .collect::<Vec<_>>();
        let lines = self.lines.iter().map(|l| (l.id, vector(&l.loops, l.head))).collect::<Vec<_>>();
        let size_of = |routes: &[(Vec<i64>, i64)]| {
            routes
                .iter()
                .map(|(s, x)| s.iter().filter(|c| **c != 0).count() + usize::from(*x != 0))
                .sum::<usize>()
        };
        let stored_size = size_of(&lines.iter().map(|(_, v)| v.clone()).collect::<Vec<_>>());
        let mut best: Option<(usize, Vec<i64>, Vec<usize>)> = None;
        let mut trees = 0;
        let e = self.lines.len();
        let mut choose = (0..n - 1).collect::<Vec<usize>>();
        loop {
            // union-find on the chosen tree lines
            let mut parent = (0..n).collect::<Vec<usize>>();
            fn find(p: &mut [usize], mut x: usize) -> usize {
                while p[x] != x {
                    p[x] = p[p[x]];
                    x = p[x];
                }
                x
            }
            let mut ok = true;
            for &k in &choose {
                let (a, b) = (find(&mut parent, index(self.lines[k].source)), find(&mut parent, index(self.lines[k].target)));
                if a == b {
                    ok = false;
                    break;
                }
                parent[a] = b;
            }
            if ok {
                trees += 1;
                let cotree = (0..e).filter(|k| !choose.contains(k)).collect::<Vec<_>>();
                let columns = cotree.iter().map(|&k| lines[k].1 .0.clone()).collect::<Vec<_>>();
                let mut routes = Vec::with_capacity(e);
                let mut fine = true;
                for (_, (l, x)) in &lines {
                    match solve_integral(&columns, l) {
                        Some(sigma) => {
                            let xs = x - sigma.iter().zip(&cotree).map(|(s, &k)| s * lines[k].1 .1).sum::<i64>();
                            routes.push((sigma, xs));
                        }
                        None => {
                            fine = false;
                            break;
                        }
                    }
                }
                if fine {
                    let size = size_of(&routes);
                    let ids = cotree.iter().map(|&k| lines[k].0).collect::<Vec<_>>();
                    let better = match &best {
                        None => true,
                        Some((bs, bids, _)) => size < *bs || (size == *bs && ids < *bids),
                    };
                    if better {
                        best = Some((size, ids, cotree));
                    }
                }
            }
            // next combination of n-1 out of e
            let k = n - 1;
            let mut i = k;
            while i > 0 && choose[i - 1] == e - k + i - 1 {
                i -= 1;
            }
            if i == 0 {
                break;
            }
            choose[i - 1] += 1;
            for j in i..k {
                choose[j] = choose[j - 1] + 1;
            }
        }
        let (size, cotree_ids, cotree) = best.ok_or("no spanning tree")?;
        let columns = cotree.iter().map(|&k| lines[k].1 .0.clone()).collect::<Vec<_>>();
        let mut routes = std::collections::BTreeMap::new();
        for (id, (l, x)) in &all_edges {
            let sigma = solve_integral(&columns, l).ok_or("an edge outside the cotree span")?;
            let xs = x - sigma.iter().zip(&cotree).map(|(s, &k)| s * lines[k].1 .1).sum::<i64>();
            routes.insert(*id, (sigma, xs));
        }
        Ok(CompactBasis {
            cotree: cotree_ids,
            routes,
            size,
            stored_size,
            trees,
        })
    }

    /// Route `product` in a compact basis: `Q(e, args..)` -> `Σ σ_i K(i, args..)
    /// + x P(args..)` (FeynKit's symbols; numerical multiples collected as
    /// `route_expression` does).
    pub fn route_compact(&self, product: &Atom, basis: &CompactBasis) -> std::result::Result<Atom, String> {
        let q = feynkit_graph::symbols::momentum();
        let k = feynkit_graph::symbols::loop_momentum();
        let mut failure = None;
        let routed = product.replace_map(|view, _, out| {
            let AtomView::Fun(f) = view else {
                return;
            };
            if f.get_symbol() != q || f.get_nargs() < 1 {
                return;
            }
            let mut arguments = f.iter();
            let Some(edge) = arguments.next().and_then(|a| i64::try_from(a).ok()) else {
                failure = Some("a symbolic edge id".to_string());
                return;
            };
            let rest = arguments.collect::<Vec<_>>();
            let Some((sigma, x)) = basis.routes.get(&edge) else {
                failure = Some(format!("an edge {edge} without a route"));
                return;
            };
            let mut sum = Atom::Zero;
            for (i, c) in sigma.iter().enumerate() {
                if *c != 0 {
                    sum += Atom::num(*c)
                        * symbolica::atom::FunctionBuilder::new(k)
                            .add_arg(Atom::num(i as i64))
                            .add_args(&rest)
                            .finish();
                }
            }
            if *x != 0 {
                sum += Atom::num(*x)
                    * symbolica::atom::FunctionBuilder::new(self.independent.0)
                        .add_arg(&self.independent.1)
                        .add_args(&rest)
                        .finish();
            }
            **out = sum;
        });
        match failure {
            Some(error) => Err(error),
            None => Ok(routed.expand_num()),
        }
    }
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
    let mut out_edges = (0..loops)
        .map(|i| PEdge {
            id: LOOP_EDGE_BASE + i as i64,
            loops: (0..loops).map(|j| i64::from(i == j)).collect(),
            head: 0,
        })
        .collect::<Vec<_>>();
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
        loop_edges: None,
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
        if f.get_symbol() == feynkit_graph::symbols::loop_momentum() {
            let i = i64::try_from(first).ok()?;
            return match &self.loop_edges {
                Some(cotree) => usize::try_from(i).ok().and_then(|i| cotree.get(i)).map(|&e| VectorRef::Edge(e)),
                None => (0..LOOP_EDGE_BASE)
                    .contains(&i)
                    .then_some(VectorRef::Edge(LOOP_EDGE_BASE + i)),
            };
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
    options: DirectOptions,
    seconds: &mut DirectSeconds,
) -> Result<DirectOutcome> {
    let mut lines = match lines_of(diagram, momentum) {
        Ok(lines) => lines,
        Err(reason) => return Ok(DirectOutcome::Declined(reason)),
    };
    let chosen = if options.compact {
        let started = Instant::now();
        let chosen = match lines.compact_basis() {
            Ok(chosen) => chosen,
            Err(reason) => return Ok(DirectOutcome::Declined(reason)),
        };
        lines.loop_edges = Some(chosen.cotree.clone());
        seconds.routing += started.elapsed().as_secs_f64();
        seconds.routing_size += chosen.size;
        seconds.stored_routing_size += chosen.stored_size;
        seconds.trees += chosen.trees;
        Some(chosen)
    } else {
        None
    };
    let started = Instant::now();
    let numerators = contract_edges(
        diagram,
        branches,
        momentum,
        q_squared,
        projection,
        options,
        chosen.as_ref().map(|chosen| (&lines, chosen)),
    )?;
    seconds.contract += started.elapsed().as_secs_f64();
    let started = Instant::now();
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
        let mut dots = std::collections::BTreeSet::new();
        scalar.as_view().visitor(&mut |node| {
            if let AtomView::Fun(f) = node
                && f.get_symbol() == crate::symbols::S.dot
            {
                dots.insert(node.to_owned());
                return false;
            }
            true
        });
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
            distinct_dots: dots.len(),
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
