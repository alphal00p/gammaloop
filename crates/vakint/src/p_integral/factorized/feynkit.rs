//! [`FactorizedNumerator::from_feynkit`]: one value per non-zero colour branch
//! of a FeynGen [`FeynmanDiagram`] (e.g. loaded from the beta thread's DOT
//! files with `FeynmanDiagram::from_dot`).
//!
//! Steps (each documented in N/notes/D1-INPUT-TYPE.md):
//! 1. every edge and vertex numerator is exported to source notation
//!    (`NumeratorInterface::export`: `hedge(h,1)` / `edge(e,k)` /
//!    `vertex(v,k)` indices, not the cooked in-memory ones);
//! 2. `mink(4, ..)` → `mink(D, ..)`, once, here (beta pitfall (v));
//! 3. the colour split ([`split_local_colour`]) of every local, the Cartesian
//!    product of rows, one colour coefficient per branch
//!    ([`colour_coefficient`], closing the two external adjoint legs with
//!    `δ_ab / N_A`); zero branches dropped;
//! 4. per branch: index-free, momentum-free factors of each residual (couplings,
//!    `𝑖`, numbers) move to [`Prefactor::scalars`];
//! 5. `overall_factor` → [`Prefactor::overall`] (evaluated) and
//!    [`Prefactor::loop_sign`] (its `InternalFermionLoopSign`); n_f power =
//!    closed fermion loops; the value is validated.

use std::collections::BTreeMap;

use feynkit_graph::{DiagramEndpoint, DiagramHalfEdge, EdgeId, FeynmanDiagram, NumeratorInterface};
use feynkit_model::Model;
use spenso::network::library::symbolic::ETS;
use symbolica::atom::{Atom, AtomCore, AtomView, Symbol};

use super::{
    ColourBranchInfo, ColourMode, ColourRow, Endpoint, FactorizedEdge, FactorizedError,
    FactorizedKind, FactorizedNumerator, FactorizedVertex, LineKind, Prefactor, Projector, Result,
    colour_coefficient, contains_tensor_content, merge_colour_channels, names, split_local_colour,
};

/// Which projector [`FactorizedNumerator::from_feynkit`] attaches to the two
/// external Lorentz ports.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum ProjectorChoice {
    /// No projector: the value keeps its open Lorentz ports.
    None,
    /// `g(port₀, port₁)`: the beta's root metric closure.
    #[default]
    Metric,
    /// `p² g(port₀, port₁) − Q(x, port₀) Q(x, port₁)` with `p2_power = −1`
    /// (the transverse projector; `x` = the first independent external edge).
    Transverse,
}

/// Options of [`FactorizedNumerator::from_feynkit`].
#[derive(Clone, Debug)]
pub struct FromFeynkitOptions {
    /// The `D` of `mink(D, ..)`; default `gammalooprs::dim` (gammaloop's
    /// `GS.dim`).
    pub dimension: Symbol,
    pub colour: ColourMode,
    pub projector: ProjectorChoice,
    /// Move index-free, momentum-free factors of each local into
    /// [`Prefactor::scalars`] (default true).
    pub pull_scalars: bool,
    /// MERGED COLOUR TERMS (helper VH lane share, experiment A1; default
    /// false = one value per non-zero colour channel): the channels are
    /// rank-factorised exactly ([`merge_colour_channels`]) and one value is
    /// returned per merged term, its multi-row locals carrying the rational
    /// combination `Σ_k w_k residual_k` of their rows; used only when it gives
    /// fewer values than channels (as the beta crate). The sum over the values
    /// of `prefactor.total() × contraction` is the same graph.
    pub merge_colour: bool,
}

impl Default for FromFeynkitOptions {
    fn default() -> Self {
        Self {
            dimension: FactorizedNumerator::default_dimension(),
            colour: ColourMode::default(),
            projector: ProjectorChoice::default(),
            pull_scalars: true,
            merge_colour: false,
        }
    }
}

fn ctor(detail: impl Into<String>) -> FactorizedError {
    FactorizedError::Constructor {
        detail: detail.into(),
    }
}

/// `mink(4, i)` → `mink(D, i)`.
pub(crate) fn to_d_dim(atom: &Atom, dimension: Symbol) -> Atom {
    let d = Atom::var(dimension);
    atom.replace_map(|part, _ctx, out| {
        let AtomView::Fun(f) = part else { return };
        if f.get_symbol().get_name() != names::MINK || f.get_nargs() != 2 {
            return;
        }
        let args = f.iter().collect::<Vec<_>>();
        if i64::try_from(args[0]) != Ok(4) {
            return;
        }
        **out = symbolica::atom::FunctionBuilder::new(f.get_symbol())
            .add_arg(&d)
            .add_arg(args[1])
            .finish();
    })
}

/// Split `residual` into (scalar part, tensor part).
fn pull_scalars(residual: &Atom) -> (Atom, Atom) {
    match residual.as_view() {
        AtomView::Mul(m) => {
            let (mut scalar, mut tensor) = (Atom::num(1), Atom::num(1));
            for part in m.iter() {
                if contains_tensor_content(part) {
                    tensor *= part.to_owned();
                } else {
                    scalar *= part.to_owned();
                }
            }
            (scalar, tensor)
        }
        view if !contains_tensor_content(view) => (residual.clone(), Atom::num(1)),
        _ => (Atom::num(1), residual.clone()),
    }
}

/// `overall_factor` → (evaluated weight without the loop sign, loop sign).
fn split_overall(overall: &Atom) -> Result<(Atom, i64)> {
    let factors = match overall.as_view() {
        AtomView::Mul(m) => m.iter().collect::<Vec<_>>(),
        other => vec![other],
    };
    let mut rest = Atom::num(1);
    let mut sign = 1i64;
    for factor in factors {
        if let AtomView::Fun(f) = factor
            && f.get_symbol().get_stripped_name() == "InternalFermionLoopSign"
        {
            let value = f
                .iter()
                .next()
                .and_then(|a| i64::try_from(a).ok())
                .filter(|s| *s == 1 || *s == -1)
                .ok_or_else(|| ctor(format!("bad loop sign {}", factor.to_owned())))?;
            sign *= value;
        } else {
            rest *= factor.to_owned();
        }
    }
    let rest = feynkit_graph::expressions::evaluate_overall_factor(rest.as_view());
    if contains_tensor_content(rest.as_view()) {
        return Err(ctor(format!("overall factor {rest} is not a scalar")));
    }
    Ok((rest, sign))
}

/// Per value of a graph: the row choice (informational), the colour
/// coefficient and the `(scalar, tensor)` of every local.
type Term = (Vec<usize>, Atom, Vec<(Atom, Atom)>);

/// The merged colour terms of a graph (option `merge_colour`): `rows` per
/// local, `branches` = the non-zero channels (row choice per local,
/// coefficient). The multi-row locals of a term carry `Σ_k w_k residual_k`
/// over the rows with non-zero weight; their scalars (pulled when `pull`)
/// must be rational multiples of each other (the ratio joins the weight).
/// `None`: nothing to merge, no fewer terms than channels, or scalars that
/// are not rational multiples (the caller keeps the channels).
fn merged_terms(
    rows: &[Vec<ColourRow>],
    branches: &[(Vec<usize>, Atom)],
    pull: bool,
) -> Option<Vec<Term>> {
    let multi = rows
        .iter()
        .enumerate()
        .filter(|(_, local_rows)| local_rows.len() > 1)
        .map(|(position, _)| position)
        .collect::<Vec<_>>();
    if multi.is_empty() || branches.len() < 2 {
        return None;
    }
    let dims = multi.iter().map(|p| rows[*p].len()).collect::<Vec<_>>();
    let channels = branches
        .iter()
        .map(|(choice, coefficient)| {
            (
                multi.iter().map(|p| choice[*p]).collect::<Vec<_>>(),
                coefficient.clone(),
            )
        })
        .collect::<Vec<_>>();
    let debug = std::env::var_os("VAKINT_SHARE_MERGE_DEBUG").is_some();
    let Some(merged) = merge_colour_channels(&dims, &channels) else {
        if debug {
            eprintln!(
                "[merge] refused: a coefficient is not an exact number combination: {channels:?}"
            );
        }
        return None;
    };
    if merged.len() >= branches.len() {
        if debug {
            eprintln!(
                "[merge] {} terms for {} channels (dims {dims:?}): kept",
                merged.len(),
                branches.len()
            );
        }
        return None;
    }
    let split = |residual: &Atom| {
        if pull {
            pull_scalars(residual)
        } else {
            (Atom::num(1), residual.clone())
        }
    };
    let mut out = Vec::with_capacity(merged.len());
    for term in merged {
        let mut choice = Vec::with_capacity(rows.len());
        let mut parts = Vec::with_capacity(rows.len());
        for (position, local_rows) in rows.iter().enumerate() {
            let Some(axis) = multi.iter().position(|p| *p == position) else {
                choice.push(0);
                parts.push(split(&local_rows[0].residual));
                continue;
            };
            let mut scalar: Option<Atom> = None;
            let mut first = None;
            let mut sum = Vec::new();
            for (k, weight) in term.weights[axis].iter().enumerate() {
                if weight.is_zero() {
                    continue;
                }
                first.get_or_insert(k);
                let (s, tensor) = split(&local_rows[k].residual);
                let ratio = match &scalar {
                    None => {
                        scalar = Some(s);
                        symbolica::domains::rational::Rational::one()
                    }
                    Some(reference) => {
                        let ratio = (&s / reference).expand();
                        let Some(ratio) = super::colour::rational_number(ratio.as_view()) else {
                            if debug {
                                eprintln!(
                                    "[merge] refused: row scalars {s} / {reference} = {ratio}"
                                );
                            }
                            return None;
                        };
                        ratio
                    }
                };
                sum.push(Atom::num(weight * &ratio) * tensor);
            }
            choice.push(first?);
            parts.push((scalar?, Atom::add_many(sum.iter())));
        }
        out.push((choice, term.coefficient, parts));
    }
    Some(out)
}

/// The local factors of the diagram, vertices (by id) then edges (by id).
enum Local {
    Vertex(usize),
    Edge(usize),
}

impl FactorizedNumerator {
    /// One value per non-zero colour branch of `diagram` (see the module
    /// documentation of `factorized::feynkit` for the steps). `model` must be
    /// the diagram's model (same fingerprint).
    pub fn from_feynkit(
        diagram: &FeynmanDiagram,
        model: &Model,
        options: &FromFeynkitOptions,
    ) -> Result<Vec<Self>> {
        crate::Vakint::initialize_vakint_symbols();
        if diagram.model().fingerprint() != model.fingerprint() {
            return Err(ctor("the model is not the diagram's model"));
        }
        if !diagram.projector().expression().is_one() {
            return Err(ctor(
                "the diagram stores a projector; v1 attaches its own (FromFeynkitOptions::projector)",
            ));
        }
        let export = |tensor: &idenso::tensor::SymbolicTensor<
            spenso::structure::partial::PartialStructure,
        >| NumeratorInterface::export(tensor).0;
        let numerator_prefactor = export(diagram.numerator_prefactor());
        if contains_tensor_content(numerator_prefactor.as_view()) {
            return Err(ctor("numerator_prefactor is not a scalar"));
        }
        let lmb = diagram.loop_momentum_basis();
        let loops = lmb.loop_edges.len();
        let independent = lmb
            .external_edges
            .iter()
            .filter(|e| !lmb.dependent_externals.contains(e))
            .copied()
            .collect::<Vec<_>>();
        // Half-edge ids of every edge end.
        let hedge = |edge: EdgeId, endpoint| {
            diagram
                .half_edge_id(DiagramHalfEdge { edge, endpoint })
                .map(|h| h.0)
        };
        let mut edges = Vec::new();
        let mut locals = Vec::<(Local, Atom)>::new();
        for (id, vertex) in diagram.vertices() {
            let numerator = to_d_dim(&export(&vertex.numerator), options.dimension);
            locals.push((Local::Vertex(id.0), numerator));
        }
        let mut edge_locals = Vec::new();
        for (id, ends, edge) in diagram.edges() {
            if edge.is_dummy {
                return Err(ctor(format!("edge {}: dummy edges are not in v1", id.0)));
            }
            let particle = model
                .particle_by_id(edge.particle)
                .map_err(|e| ctor(e.to_string()))?;
            let line = if particle.ghost_number != 0 {
                LineKind::Ghost
            } else {
                match particle.spin {
                    1 => LineKind::Scalar,
                    2 => LineKind::Fermion,
                    3 => LineKind::Vector,
                    spin => return Err(ctor(format!("edge {}: spin {spin}", id.0))),
                }
            };
            if line.is_directed_loop() && !edge.directed {
                return Err(ctor(format!("edge {}: undirected {}", id.0, line.as_str())));
            }
            let mass = model
                .parameter_by_id(particle.mass)
                .map_err(|e| ctor(e.to_string()))?;
            let mass_squared = if mass.name == "ZERO" {
                Atom::num(0)
            } else {
                Atom::var(symbolica::symbol!(&format!("UFO::{}", mass.name))).pow(Atom::num(2))
            };
            let signature = lmb
                .edge_signatures
                .get(&id)
                .ok_or_else(|| ctor(format!("edge {}: no momentum signature", id.0)))?;
            let loop_part = signature.loops.integer_coefficients();
            let external_part = signature.external.integer_coefficients();
            if loop_part.len() != loops || external_part.len() != lmb.external_edges.len() {
                return Err(ctor(format!("edge {}: signature lengths", id.0)));
            }
            let mut momentum = loop_part.iter().map(|c| *c as i64).collect::<Vec<_>>();
            for (position, external) in lmb.external_edges.iter().enumerate() {
                if lmb.dependent_externals.contains(external) {
                    if external_part[position] != 0 {
                        return Err(ctor(format!(
                            "edge {}: a dependent external momentum in its signature",
                            id.0
                        )));
                    }
                } else {
                    momentum.push(external_part[position] as i64);
                }
            }
            let source_hedge = hedge(id, DiagramEndpoint::Source);
            let target_hedge = hedge(id, DiagramEndpoint::Target);
            edges.push(FactorizedEdge {
                id: id.0,
                source: ends
                    .source
                    .map_or(Endpoint::External, |v| Endpoint::Vertex(v.0)),
                target: ends
                    .target
                    .map_or(Endpoint::External, |v| Endpoint::Vertex(v.0)),
                source_hedge: ends.source.and(source_hedge),
                target_hedge: ends.target.and(target_hedge),
                particle: particle.name.clone(),
                line,
                power: i64::from(edge.propagator_power),
                mass_squared,
                momentum,
                factor: Atom::num(1),
            });
            let numerator = to_d_dim(&export(&edge.numerator), options.dimension);
            edge_locals.push((Local::Edge(id.0), numerator));
        }
        locals.sort_by_key(|(l, _)| match l {
            Local::Vertex(v) => (0, *v),
            Local::Edge(e) => (1, *e),
        });
        edge_locals.sort_by_key(|(l, _)| match l {
            Local::Edge(e) => *e,
            Local::Vertex(v) => *v,
        });
        locals.extend(edge_locals);
        edges.sort_by_key(|e| e.id);
        let rows = locals
            .iter()
            .map(|(_, local)| split_local_colour(local))
            .collect::<Result<Vec<_>>>()?;
        // Branches: Cartesian product of the rows (row choices per local).
        let mut choices = vec![Vec::<usize>::new()];
        for local_rows in &rows {
            choices = choices
                .into_iter()
                .flat_map(|c| {
                    (0..local_rows.len()).map(move |r| {
                        let mut next = c.clone();
                        next.push(r);
                        next
                    })
                })
                .collect();
        }
        let mut branches = Vec::new();
        for choice in choices {
            let colour = choice
                .iter()
                .zip(&rows)
                .fold(Atom::num(1), |acc, (r, local_rows)| {
                    acc * &local_rows[*r].colour
                });
            let coefficient = colour_coefficient(&colour, options.colour)?;
            if !coefficient.is_zero() {
                branches.push((choice, coefficient));
            }
        }
        // The graph data shared by every branch.
        let (overall, loop_sign) = split_overall(diagram.overall_factor())?;
        let vertex_rules = diagram
            .vertices()
            .map(|(id, vertex)| {
                let rule = match vertex.interaction {
                    Some(rule) => model
                        .vertex_rule_by_id(rule)
                        .map_err(|e| ctor(e.to_string()))?
                        .name
                        .clone(),
                    None => "none".into(),
                };
                Ok((id.0, rule))
            })
            .collect::<Result<BTreeMap<_, _>>>()?;
        let loop_edges = lmb.loop_edges.iter().map(|e| Some(e.0)).collect::<Vec<_>>();
        let external_edges = independent.iter().map(|e| Some(e.0)).collect::<Vec<_>>();
        // Per value: the row choice (informational), the colour coefficient
        // and the (scalar, tensor) of every local.
        let unmerged = |branches: Vec<(Vec<usize>, Atom)>| {
            branches
                .into_iter()
                .map(|(choice, coefficient)| {
                    let parts = choice
                        .iter()
                        .zip(&rows)
                        .map(|(r, local_rows)| {
                            let residual = &local_rows[*r].residual;
                            if options.pull_scalars {
                                pull_scalars(residual)
                            } else {
                                (Atom::num(1), residual.clone())
                            }
                        })
                        .collect::<Vec<_>>();
                    (choice, coefficient, parts)
                })
                .collect::<Vec<_>>()
        };
        let terms = if options.merge_colour {
            match merged_terms(&rows, &branches, options.pull_scalars) {
                Some(merged) => merged,
                None => unmerged(branches),
            }
        } else {
            unmerged(branches)
        };
        let count = terms.len();
        let mut values = Vec::with_capacity(count);
        for (index, (choice, coefficient, parts)) in terms.into_iter().enumerate() {
            let mut scalars = numerator_prefactor.clone();
            let mut edges = edges.clone();
            let mut vertices = Vec::new();
            for ((local, _), (scalar, tensor)) in locals.iter().zip(parts) {
                scalars *= scalar;
                match local {
                    Local::Vertex(v) => vertices.push(FactorizedVertex {
                        id: *v,
                        rule: vertex_rules[v].clone(),
                        factor: tensor,
                    }),
                    Local::Edge(e) => {
                        edges.iter_mut().find(|edge| edge.id == *e).unwrap().factor = tensor
                    }
                }
            }
            vertices.sort_by_key(|v| v.id);
            let mut value = Self {
                label: format!("graph={}", diagram.name()),
                kind: FactorizedKind::Direct,
                dimension: options.dimension,
                loop_edges: loop_edges.clone(),
                external_edges: external_edges.clone(),
                edges,
                vertices,
                ports: Vec::new(),
                projector: None,
                prefactor: Prefactor {
                    overall: overall.clone(),
                    loop_sign,
                    nf_power: 0,
                    scalars,
                    colour: coefficient,
                },
                colour_branch: Some(ColourBranchInfo {
                    index,
                    count,
                    rows: choice,
                }),
                insertions: Vec::new(),
            };
            value.ports = value.open_lorentz_ports();
            value.projector = value.projector_for(options.projector)?;
            value.prefactor.nf_power = value
                .closed_loops()?
                .iter()
                .filter(|l| l.line == LineKind::Fermion)
                .count() as u32;
            value.validate()?;
            values.push(value);
        }
        Ok(values)
    }

    /// The Lorentz slots `mink(D, hedge(h, k))` of external half-edges that
    /// occur in the factors, sorted by canonical string.
    fn open_lorentz_ports(&self) -> Vec<Atom> {
        let external_hedges = self
            .edges
            .iter()
            .filter(|e| e.is_external())
            .flat_map(|e| e.source_hedge.into_iter().chain(e.target_hedge))
            .collect::<std::collections::BTreeSet<_>>();
        let mut ports = BTreeMap::new();
        let factors = self
            .edges
            .iter()
            .map(|e| &e.factor)
            .chain(self.vertices.iter().map(|v| &v.factor));
        for factor in factors {
            factor.visitor(&mut |part| {
                if let Some(slot) = super::as_slot(part)
                    && slot.rep == names::MINK
                    && slot.hedge().is_some_and(|h| external_hedges.contains(&h))
                {
                    ports.insert(slot.key(), slot.atom);
                }
                true
            });
        }
        ports.into_values().collect()
    }

    fn projector_for(&self, choice: ProjectorChoice) -> Result<Option<Projector>> {
        if choice == ProjectorChoice::None {
            return Ok(None);
        }
        let [a, b] = self.ports.as_slice() else {
            return Err(ctor(format!(
                "a two-port projector needs two Lorentz ports, found {}",
                self.ports.len()
            )));
        };
        let metric = ETS.metric(a.clone(), b.clone());
        match choice {
            ProjectorChoice::None => unreachable!(),
            ProjectorChoice::Metric => Ok(Some(Projector {
                tensor: metric,
                p2_power: 0,
            })),
            ProjectorChoice::Transverse => {
                let Some(Some(x)) = self.external_edges.first() else {
                    return Err(ctor("the transverse projector needs an external edge"));
                };
                let q = feynkit_graph::momentum_symbol();
                let AtomView::Fun(port) = a.as_view() else {
                    unreachable!("ports are slots")
                };
                let dummy = symbolica::atom::FunctionBuilder::new(port.get_symbol())
                    .add_arg(Atom::var(self.dimension))
                    .add_arg(symbolica::function!(
                        feynkit_graph::symbols::edge_index(),
                        *x as i64,
                        1
                    ))
                    .finish();
                let qq = symbolica::function!(q, *x as i64, &dummy)
                    * symbolica::function!(q, *x as i64, &dummy);
                let pp =
                    symbolica::function!(q, *x as i64, a) * symbolica::function!(q, *x as i64, b);
                Ok(Some(Projector {
                    tensor: qq * metric - pp,
                    p2_power: -1,
                }))
            }
        }
    }
}
