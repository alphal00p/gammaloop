//! Bare poles of a fully massive, single-scale vacuum graph from the ordinary
//! R operation:
//!
//! `K G = -Delta(G) - K sum_{S proper, non-empty UV spinney} Delta(S) * G/S`.
//!
//! The proper children `Delta(S)` are the MS counterterms of the existing
//! massive R operation (Wood forest, `PolePart`, all-mass Taylor kernel around
//! `m_uv`), evaluated with the PHYSICAL masses of the root, exactly as
//! [`RStarPlan::evaluate_uv_children`] does for the massless R* roots. The
//! quotients `G/S` keep their physical masses and are integrated as ordinary
//! single-scale massive vacuum integrals by Vakint; no Taylor expansion, no
//! auxiliary mass and no p-integral enters them.
//!
//! The overall counterterm `Delta(G)` of a root with `omega >= 0` is a pure
//! number times `m^omega`. Two independent routes are provided:
//! - [`RStarPlan::massive_vacuum_counterterm`]: the ordinary massive R
//!   operation on the root itself (needs `L`-loop `m_uv` tadpole poles);
//! - [`VacuumTopology::overall_counterterm_rstar`]: the mass expansion
//!   `Delta(G) = (m^2)^(omega/2) sum_{|a| = omega/2} prod_e binom(p_e+a_e-1, a_e)
//!   Delta(G_0 with powers p+a)` over LOG-divergent massless vacuum graphs,
//!   each by [`RStarPlan::counterterm`] (p-integrals with at most `L-1` loops).
//!   A line of power `q` is a chain of `q` lines through two-point vertices
//!   with unit numerators.

use std::time::Instant;

use eyre::{Result, WrapErr, ensure, eyre};
use idenso::{
    CookMode, CookSettings,
    dirac::GammaSimplifySettings,
    tensor::{AlgebraSettings, SymbolicTensor},
};
use linnet::half_edge::subgraph::{InternalSubGraph, SuBitGraph, SubSetLike, SubSetOps};
use spenso::{
    network::tags::SPENSO_TAG,
    structure::representation::{BaseRepName, Minkowski, RepName},
};
use symbolica::{
    atom::{Atom, AtomCore},
    function,
};
use vakint::{Vakint, VakintSettings};

use crate::{
    graph::{FeynmanGraph, Graph, LMBext, parse::IntoFinalizedRuntimeGraph},
    model::Model,
    utils::{GS, W_},
    uv::{
        ApproximationType, RenormalizationPrescriptionSettings, UVgenerationSettings,
        UltravioletGraph, approx::integrated::to_vakint_integrand, uv_graph::UVE,
    },
};

use super::{GraphRegion, RStarPlan, expression::RStarIntegrand};

/// The ordinary-R pole decomposition of a massive vacuum root.
#[derive(Clone, Debug)]
pub struct MassiveVacuumSubdivergences {
    /// `-K sum_S Delta(S) G/S`, summed over all proper non-empty spinneys.
    pub poles: Atom,
    /// One entry per spinney: its label, child loop counts, quotient loops,
    /// the un-truncated contribution `-Delta(S) G/S` (through the requested
    /// depth) and the wall time of its quotient integration.
    pub terms: Vec<MassiveVacuumTerm>,
    /// Wall time of the proper UV children (massive R operation).
    pub children_ms: u128,
    /// The proper children `Delta(S)` (support label, loops, polynomial).
    pub children: Vec<(String, usize, Vec<usize>, Atom)>,
}

#[derive(Clone, Debug)]
pub struct MassiveVacuumTerm {
    pub label: String,
    pub child_loops: Vec<usize>,
    pub quotient_loops: usize,
    pub contribution: Atom,
    pub elapsed_ms: u128,
}

/// UV settings with the MS (pole-part) prescription everywhere, as R* uses.
pub fn massive_vacuum_settings(settings: &UVgenerationSettings) -> UVgenerationSettings {
    UVgenerationSettings {
        subtract_uv: true,
        renormalization_prescription: RenormalizationPrescriptionSettings {
            log_divergent: ApproximationType::PolePart,
            massive_power_divergent: ApproximationType::PolePart,
            massless_power_divergent: ApproximationType::PolePart,
            overrides: Vec::new(),
        },
        ..settings.clone()
    }
}

impl RStarPlan {
    /// `-K sum_{S proper, non-empty} Delta(S) G/S` for a massive vacuum root.
    ///
    /// `mass` is the common physical mass of the root's lines (the mass
    /// itself, not its square). The proper children `Delta(S)` are those of
    /// the ordinary massive R operation (`evaluate_uv_children`, physical
    /// masses on every line). Each quotient is requested through `epsilon^-1`
    /// of the complete product with its inserted poles, plus `extra_terms`
    /// further orders (for depth checks). The overall counterterm of the root
    /// is NOT included (see [`Self::massive_vacuum_counterterm`]).
    pub fn massive_vacuum_subdivergence_poles(
        &self,
        graph: &Graph,
        model: &Model,
        mass: &Atom,
        vakint: &Vakint,
        uv_settings: &UVgenerationSettings,
        extra_terms: i64,
    ) -> Result<MassiveVacuumSubdivergences> {
        ensure!(
            graph.external_filter::<SuBitGraph>().is_empty(),
            "a massive vacuum root has no external legs"
        );
        ensure!(
            graph
                .iter_edges_of(&self.support)
                .all(|(_, _, data)| data.data.mass_atom(model) == *mass),
            "every line of a massive vacuum root must carry the given mass"
        );
        let uv_settings = massive_vacuum_settings(uv_settings);
        let started = Instant::now();
        let children = if self.counterterms.is_empty() {
            Default::default()
        } else {
            self.evaluate_uv_children(graph, model, mass, vakint, &uv_settings)
                .wrap_err("while computing the proper UV counterterms with massive R")?
        };
        let children_ms = started.elapsed().as_millis();
        let base = uv_settings.vakint.true_settings();
        let epsilon = GS.dim_epsilon;
        let mut total = Atom::Zero;
        let mut terms = Vec::new();
        for (index, term) in self
            .ordinary_uv
            .iter()
            .filter(|term| !term.counterterms.is_empty())
            .enumerate()
        {
            let started = Instant::now();
            let insertions = term
                .counterterms
                .iter()
                .map(|call| {
                    children
                        .get(&call.argument.support)
                        .cloned()
                        .ok_or_else(|| eyre!("missing a completed proper UV counterterm"))
                })
                .collect::<Result<Vec<_>>>()?;
            let mut settings = base.clone();
            settings.number_of_terms_in_epsilon_expansion =
                i64::try_from(term.remainder.loops)? + extra_terms;
            let contribution = -integrate_massive_vacuum(
                graph,
                model,
                &term.remainder,
                &insertions,
                vakint,
                &uv_settings,
                &settings,
                10_000 + index,
            )
            .wrap_err_with(|| {
                format!(
                    "while integrating the massive quotient {}",
                    term.remainder.remaining.string_label()
                )
            })?;
            total += &contribution;
            terms.push(MassiveVacuumTerm {
                label: term
                    .counterterms
                    .iter()
                    .map(|call| call.argument.support.string_label())
                    .collect::<Vec<_>>()
                    .join("+"),
                child_loops: term
                    .counterterms
                    .iter()
                    .map(|call| call.argument.loops)
                    .collect(),
                quotient_loops: term.remainder.loops,
                contribution,
                elapsed_ms: started.elapsed().as_millis(),
            });
        }
        let children = self
            .counterterms
            .iter()
            .filter_map(|call| {
                children.get(&call.argument.support).map(|value| {
                    (
                        call.argument.support.string_label(),
                        call.argument.loops,
                        graph
                            .iter_edges_of(&call.argument.support)
                            .map(|(_, edge, _)| usize::from(edge))
                            .collect(),
                        value.clone(),
                    )
                })
            })
            .collect();
        Ok(MassiveVacuumSubdivergences {
            poles: total.series(epsilon, 0, -1)?.to_atom().expand(),
            terms,
            children_ms,
            children,
        })
    }

    /// `Delta(G)` of the massive root by the ordinary massive R operation:
    /// the forest of the complete catalogue, root included, with `L`-loop
    /// `m_uv` tadpoles (so only up to the loop count of Vakint's vacuum data).
    /// Returns the MS counterterm `-K R' G`.
    pub fn massive_vacuum_counterterm(
        &self,
        graph: &Graph,
        model: &Model,
        vakint: &Vakint,
        uv_settings: &UVgenerationSettings,
    ) -> Result<Atom> {
        if graph.compute_dod(&self.support) < 0 {
            return Ok(Atom::Zero);
        }
        let settings = UVgenerationSettings {
            generate_integrated: true,
            add_marker: false,
            ..massive_vacuum_settings(uv_settings)
        };
        let mut forest = self.uv_forest(graph, &settings);
        forest.integrate(graph, model, vakint, &settings)?;
        Ok(forest
            .counterterm_for_support(&self.support)?
            .replace(GS.m_uv_expansion)
            .with(GS.m_uv_vacuum)
            .expand())
    }
}

/// Poles (through `epsilon^-1`, plus `extra_terms`) of the bare massive vacuum
/// graph, evaluated directly by Vakint with the same conventions as the
/// quotients. This is the reference of the identity.
pub fn massive_vacuum_direct(
    graph: &Graph,
    model: &Model,
    vakint: &Vakint,
    uv_settings: &UVgenerationSettings,
    extra_terms: i64,
) -> Result<Atom> {
    let support = graph
        .no_dummy()
        .subtract(&graph.external_filter::<SuBitGraph>());
    let loops = graph.n_loops(&support);
    let region = GraphRegion {
        boundary: graph.dummy_less_full_crown(&support),
        remaining: support.clone(),
        support,
        insertions: Vec::new(),
        loops,
    };
    let uv_settings = massive_vacuum_settings(uv_settings);
    let mut settings = uv_settings.vakint.true_settings();
    settings.number_of_terms_in_epsilon_expansion = i64::try_from(loops)? + extra_terms;
    integrate_massive_vacuum(
        graph,
        model,
        &region,
        &[],
        vakint,
        &uv_settings,
        &settings,
        9_999,
    )
}

/// Integrate a (possibly decorated) massive vacuum quotient with its physical
/// masses. Insertions follow `region.insertions`; their polynomials use the
/// original boundary momentum labels, which are routed here.
#[allow(clippy::too_many_arguments)]
pub(super) fn integrate_massive_vacuum(
    graph: &Graph,
    model: &Model,
    region: &GraphRegion,
    insertions: &[Atom],
    vakint: &Vakint,
    uv_settings: &UVgenerationSettings,
    settings: &VakintSettings,
    unique_id: usize,
) -> Result<Atom> {
    ensure!(
        insertions.len() == region.insertions.len(),
        "each contraction requires an insertion polynomial"
    );
    let mut contracted: SuBitGraph = graph.empty_subgraph();
    for slot in &region.insertions {
        contracted.union_with(&slot.support);
    }
    ensure!(
        region.support.subtract(&contracted) == region.remaining,
        "remainder does not match its contracted supports"
    );
    let mut expression = graph
        .numerator(&region.remaining, &contracted)
        .to_d_dim(GS.dim)
        .get_single_atom()?;
    for insertion in insertions {
        expression *= insertion;
    }
    if region.loops == 0 {
        return RStarIntegrand::finish_tensor(expression, unique_id);
    }
    let contracted_internal = InternalSubGraph::try_new(contracted, &graph.underlying)
        .ok_or_else(|| eyre!("an insertion contains an unpaired half-edge"))?;
    let basis = graph
        .shrunken_sub_lmb(
            &region.support,
            &contracted_internal,
            region.boundary.clone(),
        )
        .wrap_err("could not route a massive vacuum quotient")?;
    expression = expression.replace_multiple(graph.normal_emr_replacement(
        &region.remaining,
        &basis,
        &[W_.x___],
        |_| true,
    ));
    let minkowski = Minkowski {}.new_rep(GS.dim).to_symbolic([]);
    for (pair, edge, data) in graph.iter_edges_of(&region.remaining) {
        if !pair.is_paired() {
            continue;
        }
        let mass = data.data.mass_atom(model);
        let hard = basis.loop_atom::<Atom>(edge, GS.emr_mom, &[], true);
        ensure!(
            !hard.is_zero(),
            "a remaining vacuum propagator has no loop momentum"
        );
        let momentum = basis.loop_atom(edge, GS.emr_mom, std::slice::from_ref(&minkowski), true)
            + basis.ext_atom(edge, GS.emr_mom, std::slice::from_ref(&minkowski), true);
        // The third argument is the mass itself; the Vakint bridge squares it.
        expression /= function!(
            GS.den,
            usize::from(edge),
            hard,
            &mass,
            function!(SPENSO_TAG.dot, &momentum, &momentum) - mass.pow(2)
        );
    }
    let mut converter = graph.clone();
    for slot in &region.insertions {
        let nodes = converter
            .iter_nodes_of(&slot.support)
            .map(|(node, _, _)| node)
            .collect::<Vec<_>>();
        let first = nodes
            .first()
            .ok_or_else(|| eyre!("a contracted support contains no vertices"))?;
        let vertex = converter.underlying[*first].clone();
        converter.underlying.identify_nodes(&nodes, vertex);
    }
    converter.underlying.forget_identification_history();
    let expression = simplify(&expression)?;
    let mut input = to_vakint_integrand(
        &expression,
        &converter.underlying,
        &region.remaining,
        &converter.empty_subgraph::<SuBitGraph>(),
        &uv_settings.vakint,
        false,
    )?;
    input.canonicalize(settings, &vakint.topologies, false)?;
    let minkowski = Minkowski::selfless_symbol();
    let dimension = Atom::num(4) - Atom::num(2) * Atom::var(GS.dim_epsilon);
    for term in &mut input.0 {
        term.numerator = term
            .numerator
            .replace(function!(minkowski, GS.dim, W_.x___))
            .with(function!(minkowski, &dimension, W_.x___));
    }
    input.tensor_reduce(vakint, settings)?;
    input.evaluate_integral(vakint, settings)?;
    let result: Atom = input.into();
    RStarIntegrand::finish_tensor(result, unique_id)
}

fn simplify(expression: &Atom) -> Result<Atom> {
    let cooking = CookSettings::indices()
        .with_mode(CookMode::ReversibleEncoding)
        .with_representation_payloads(true, true);
    let simplified = SymbolicTensor::infer_indexed(cooking.try_cook(expression.as_view())?)?
        .undo_trace()?
        .undo_chain()?
        .simplify_algebra(&AlgebraSettings {
            gamma: Some(GammaSimplifySettings::canonical()),
            epsilon: true,
            ..Default::default()
        })?
        .to_dots()?;
    Ok(cooking.uncook(simplified.expression().as_view()))
}

/// One term of the mass expansion of `Delta(G)`: extra powers per line, its
/// binomial weight, the massless R* counterterm and its wall time (ms).
pub type DistributionTerm = (Vec<u32>, u64, Atom, u128);

/// A scalar vacuum topology: lines `(a, b, power)` between vertices.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct VacuumTopology {
    pub name: String,
    pub edges: Vec<(usize, usize, u32)>,
    /// Irreducible numerator factors `(q^2 - m^2)`, each `q` a signed sum of
    /// line momenta: `[(line index, sign)]` (line orientation as in `edges`).
    pub numerators: Vec<Vec<(usize, i64)>>,
}

/// One monomial of the expanded numerator: a scalar coefficient (with its
/// mass powers) and, per line, the numerator factors attached to it (each a
/// rank-1 `Q(line)^mu_f` contracted with the partner line of factor `f`).
#[derive(Clone, Debug)]
pub struct NumeratorMonomial {
    pub coefficient: Atom,
    pub attachments: Vec<Vec<usize>>,
}

impl VacuumTopology {
    pub fn loops(&self) -> usize {
        let vertices = self
            .edges
            .iter()
            .flat_map(|&(a, b, _)| [a, b])
            .collect::<std::collections::BTreeSet<_>>()
            .len();
        self.edges.len() + 1 - vertices
    }

    /// `4 L - 2 sum p` (scalar numerator).
    pub fn superficial_degree(&self) -> i64 {
        4 * self.loops() as i64 - 2 * self.edges.iter().map(|e| e.2 as i64).sum::<i64>()
            + self.numerator_degree()
    }

    /// Runtime DOT: a line of power `q` is a chain of `q` lines through
    /// two-point vertices; every vertex and line has a unit numerator.
    pub fn dot(&self) -> String {
        let mut lines: Vec<(String, String)> = Vec::new();
        let mut chain = 0;
        for &(a, b, power) in &self.edges {
            let mut from = format!("v{a}");
            for step in 0..power {
                let to = if step + 1 == power {
                    format!("v{b}")
                } else {
                    chain += 1;
                    format!("w{chain}")
                };
                lines.push((from.clone(), to.clone()));
                from = to;
            }
        }
        let mut names = lines
            .iter()
            .flat_map(|(a, b)| [a.clone(), b.clone()])
            .collect::<Vec<_>>();
        names.sort();
        names.dedup();
        let index = |name: &String| names.iter().position(|n| n == name).unwrap();
        let mut components = (0..names.len()).collect::<Vec<_>>();
        let mut slots = vec![0usize; names.len()];
        let mut loops = 0;
        let mut dot = format!(
            "digraph {} {{\ngraph [projector=1]; node [num=1];\nedge [particle=\"scalar_1\" num=1];\n",
            self.name.replace(['-', '.'], "_")
        );
        for (id, (a, b)) in lines.iter().enumerate() {
            let (ia, ib) = (index(a), index(b));
            let source = slots[ia];
            slots[ia] += 1;
            let sink = slots[ib];
            slots[ib] += 1;
            let basis = if components[ia] == components[ib] {
                loops += 1;
                format!("lmb_id={}", loops - 1)
            } else {
                let from = components[ib];
                let into = components[ia];
                for component in &mut components {
                    if *component == from {
                        *component = into;
                    }
                }
                String::new()
            };
            dot.push_str(&format!(
                "{a} -> {b} [id={id} {basis} source=\"{{ufo_order:{source}}}\" sink=\"{{ufo_order:{sink}}}\"];\n"
            ));
        }
        dot.push('}');
        dot
    }

    /// The runtime graph with every line of the given mass (`Atom::Zero`
    /// for the massless graph).
    pub fn graph(&self, model: &Model, mass: &Atom) -> Result<Graph> {
        let mut graph: Graph = self
            .dot()
            .into_finalized_runtime_graph(model)
            .wrap_err_with(|| format!("while building the runtime graph {}", self.name))?;
        let edges = graph
            .iter_edges()
            .map(|(_, edge, _)| edge)
            .collect::<Vec<_>>();
        for edge in edges {
            let data = &mut graph.underlying[edge];
            data.particle = data.particle.clone().override_mass(Some(mass.clone()));
        }
        Ok(graph)
    }

    /// Expand `prod_f (q_f^2 - m^2)` with `q_f = sum s_i Q_i` into monomials
    /// `prod_f [s_i s_j Q_i.Q_j or -m^2]`. Exact by linearity of R' and K;
    /// each monomial is an ordinary graph with edge-attached rank-1 numerators,
    /// so gammaloop's power counting and Taylor operation see every factor on
    /// exactly the subgraphs containing its line.
    pub fn numerator_monomials(&self, mass: &Atom) -> Vec<NumeratorMonomial> {
        let mut out = vec![NumeratorMonomial {
            coefficient: Atom::one(),
            attachments: vec![Vec::new(); self.edges.len()],
        }];
        for (f, factor) in self.numerators.iter().enumerate() {
            let mut next = Vec::new();
            for monomial in &out {
                for (a, &(i, si)) in factor.iter().enumerate() {
                    for (b, &(j, sj)) in factor.iter().enumerate() {
                        if b < a {
                            continue;
                        }
                        let weight = if a == b { si * sj } else { 2 * si * sj };
                        let mut attachments = monomial.attachments.clone();
                        attachments[i].push(f);
                        attachments[j].push(f);
                        next.push(NumeratorMonomial {
                            coefficient: &monomial.coefficient * Atom::num(weight),
                            attachments,
                        });
                    }
                }
                next.push(NumeratorMonomial {
                    coefficient: -(&monomial.coefficient * mass.pow(Atom::num(2))),
                    attachments: monomial.attachments.clone(),
                });
            }
            out = next;
        }
        out
    }

    /// The runtime graph of one numerator monomial: rank-1 momenta of the
    /// first chain segment of each line, `Q(e)^{mu_f}` with the factor's index.
    pub fn monomial_graph(
        &self,
        model: &Model,
        mass: &Atom,
        attachments: &[Vec<usize>],
    ) -> Result<Graph> {
        let mut graph = self.graph(model, mass)?;
        let mut first = 0usize;
        for (line, &(_, _, power)) in self.edges.iter().enumerate() {
            if let Some(factors) = attachments.get(line).filter(|f| !f.is_empty()) {
                let mut numerator = Atom::one();
                for f in factors {
                    let index: Atom =
                        crate::numerator::aind::Aind::Edge(4000 + *f as u16, 1).into();
                    numerator *= function!(
                        GS.emr_mom,
                        first,
                        Minkowski {}.new_rep(GS.dim).to_symbolic([index])
                    );
                }
                let data =
                    &mut graph.underlying[linnet::half_edge::involution::EdgeIndex::from(first)];
                data.num = crate::graph::Autogen::explicit(numerator);
                data.dod = crate::graph::Autogen::explicit(-2 + factors.len() as i32);
            }
            first += power as usize;
        }
        Ok(graph)
    }

    /// Total numerator degree (each factor counts 2).
    pub fn numerator_degree(&self) -> i64 {
        2 * self.numerators.len() as i64
    }

    /// All distributions `a` with `|a| = n` and their Minkowski expansion
    /// weight `prod_e binom(p_e + a_e - 1, a_e)` of
    /// `1/(k^2-m^2)^p = sum_a binom(p+a-1, a) (m^2)^a / (k^2)^(p+a)`.
    pub fn power_distributions(&self, n: u32) -> Vec<(Vec<u32>, u64)> {
        fn binomial(n: u64, k: u64) -> u64 {
            (0..k).fold(1, |acc, i| acc * (n - i) / (i + 1))
        }
        fn recurse(
            edges: &[(usize, usize, u32)],
            left: u32,
            prefix: &mut Vec<u32>,
            out: &mut Vec<(Vec<u32>, u64)>,
        ) {
            if prefix.len() == edges.len() {
                if left == 0 {
                    let weight = edges
                        .iter()
                        .zip(prefix.iter())
                        .map(|(e, &a)| binomial((e.2 + a - 1) as u64, a as u64))
                        .product();
                    out.push((prefix.clone(), weight));
                }
                return;
            }
            for a in 0..=left {
                prefix.push(a);
                recurse(edges, left - a, prefix, out);
                prefix.pop();
            }
        }
        let mut out = Vec::new();
        recurse(&self.edges, n, &mut Vec::new(), &mut out);
        out
    }

    /// `Delta(G)` of the massive root by route (ii): for every numerator
    /// monomial (coefficient `c`, rank `r`, see [`Self::numerator_monomials`])
    /// the mass expansion `c (m^2)^(w/2) sum_{|a| = w/2} prod_e binom(p_e+a_e-1, a_e)
    /// Delta(G_0 with powers p+a and the monomial's edge numerators)` with
    /// `w = 4L - 2 sum p + r`, each massless graph LOG-divergent and evaluated
    /// by [`RStarPlan::counterterm`]. Returns `(Delta, per-graph log)`.
    pub fn overall_counterterm_rstar(
        &self,
        model: &Model,
        mass: &Atom,
        vakint: &Vakint,
        uv_settings: &UVgenerationSettings,
    ) -> Result<(Atom, Vec<DistributionTerm>)> {
        let mut total = Atom::Zero;
        let mut log = Vec::new();
        for monomial in self.numerator_monomials(mass) {
            let rank: usize = monomial.attachments.iter().map(Vec::len).sum();
            let omega = self.superficial_degree() - self.numerator_degree() + rank as i64;
            if omega < 0 {
                continue;
            }
            ensure!(
                omega % 2 == 0,
                "odd superficial degree of a vacuum monomial"
            );
            let n = u32::try_from(omega / 2)?;
            let mut monomial_total = Atom::Zero;
            for (a, weight) in self.power_distributions(n) {
                let started = Instant::now();
                let dotted = VacuumTopology {
                    numerators: Vec::new(),
                    name: format!(
                        "{}_a{}",
                        self.name,
                        a.iter().map(u32::to_string).collect::<String>()
                    ),
                    edges: self
                        .edges
                        .iter()
                        .zip(&a)
                        .map(|(&(x, y, p), &extra)| (x, y, p + extra))
                        .collect(),
                };
                let massless = dotted.monomial_graph(model, &Atom::Zero, &monomial.attachments)?;
                let plan = RStarPlan::with_automatic_edge(&massless, model, uv_settings)?;
                let delta = plan
                    .counterterm(&massless, model, vakint, uv_settings)
                    .wrap_err_with(|| {
                        format!("while computing the R* counterterm of {}", dotted.name)
                    })?;
                monomial_total += Atom::num(weight as i64) * &delta;
                log.push((a, weight, delta, started.elapsed().as_millis()));
            }
            total += &monomial.coefficient * monomial_total * mass.pow(Atom::num(omega));
        }
        Ok((total.expand(), log))
    }
}
