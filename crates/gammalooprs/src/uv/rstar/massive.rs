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

use std::{collections::BTreeMap, time::Instant};

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
        UltravioletGraph,
        approx::{CutStructure, integrated::to_vakint_integrand},
        hedge_poset::Wood,
        uv_graph::UVE,
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
    pub children: Vec<(String, usize, Atom)>,
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
    /// itself, not its square). It is used for the proper children exactly as
    /// `evaluate_uv_children` uses its reference mass. Each quotient is
    /// requested through `epsilon^-1` of the complete product with its
    /// inserted poles, plus `extra_terms` further orders (for depth checks).
    /// The overall counterterm of the root is NOT included.
    pub fn massive_vacuum_subdivergence_poles(
        &self,
        graph: &Graph,
        model: &Model,
        mass: &Atom,
        vakint: &Vakint,
        uv_settings: &UVgenerationSettings,
        extra_terms: i64,
        forest_extra_terms: i64,
        standalone_children: bool,
    ) -> Result<MassiveVacuumSubdivergences> {
        ensure!(
            graph.external_filter::<SuBitGraph>().is_empty(),
            "a massive vacuum root has no external legs"
        );
        let uv_settings = massive_vacuum_settings(uv_settings);
        let started = Instant::now();
        let _ = mass;
        let children = if self.counterterms.is_empty() {
            Default::default()
        } else {
            let supports = self
                .counterterms
                .iter()
                .map(|call| call.argument.support.clone())
                .collect::<Vec<_>>();
            if standalone_children {
                let mut values = BTreeMap::new();
                for support in &supports {
                    values.insert(
                        support.clone(),
                        standalone_child_counterterm(
                            graph,
                            model,
                            support,
                            vakint,
                            &uv_settings,
                            forest_extra_terms,
                        )
                        .wrap_err("while computing a standalone child counterterm")?,
                    );
                }
                values
            } else {
                self.massive_forest_counterterms(
                    graph,
                    model,
                    vakint,
                    &uv_settings,
                    false,
                    forest_extra_terms,
                    &supports,
                )
                .wrap_err("while computing the proper UV counterterms with massive R")?
            }
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

    /// `Delta(G)` of the massive root by the ordinary massive R operation
    /// (route (i)): the Wood forest of the complete catalogue, root included,
    /// with `L`-loop `m_uv` tadpoles. Returns the MS counterterm `-K R' G`.
    pub fn massive_vacuum_counterterm(
        &self,
        graph: &Graph,
        model: &Model,
        vakint: &Vakint,
        uv_settings: &UVgenerationSettings,
        forest_extra_terms: i64,
    ) -> Result<Atom> {
        if graph.compute_dod(&self.support) < 0 {
            return Ok(Atom::Zero);
        }
        let support = self.support.clone();
        let mut values = self.massive_forest_counterterms(
            graph,
            model,
            vakint,
            &massive_vacuum_settings(uv_settings),
            true,
            forest_extra_terms,
            std::slice::from_ref(&support),
        )?;
        Ok(values.remove(&support).unwrap_or_default().expand())
    }

    /// The ordinary massive R counterterms (Wood forest, MS, all-mass Taylor
    /// kernel around `m_uv`) of the requested supports, with the lines' own
    /// (physical) masses. Unlike `evaluate_uv_children`, the forest's Vakint
    /// depth is `max_loops + forest_extra_terms` (the existing forest uses
    /// `max_loops + 1`, needed for MUV finite parts only): every node then
    /// asks its `l`-loop tadpoles through `eps^(max_loops - l - 1 + extra)`,
    /// which is exactly what the outermost pole projection needs.
    #[allow(clippy::too_many_arguments)]
    pub fn massive_forest_counterterms(
        &self,
        graph: &Graph,
        model: &Model,
        vakint: &Vakint,
        uv_settings: &UVgenerationSettings,
        include_root: bool,
        forest_extra_terms: i64,
        supports: &[SuBitGraph],
    ) -> Result<BTreeMap<SuBitGraph, Atom>> {
        let settings = UVgenerationSettings {
            generate_integrated: true,
            add_marker: false,
            ..massive_vacuum_settings(uv_settings)
        };
        let spinneys = self
            .uv_catalogue()
            .iter()
            .filter(|spinney| include_root || spinney.filter() != &self.support)
            .cloned()
            .collect::<Vec<_>>();
        let max_loops = spinneys
            .iter()
            .map(|spinney| graph.n_loops(spinney.filter()))
            .max()
            .unwrap_or(0);
        let mut wood = Wood::from_spinneys(
            spinneys,
            graph,
            CutStructure::empty(graph),
            &settings.vakint,
        );
        wood.vakint_settings.number_of_terms_in_epsilon_expansion =
            i64::try_from(max_loops)? + forest_extra_terms;
        let mut forest = wood.unfold();
        forest.integrate(graph, model, vakint, &settings)?;
        supports
            .iter()
            .map(|support| {
                Ok((
                    support.clone(),
                    forest
                        .counterterm_for_support(support)?
                        .replace(GS.m_uv_expansion)
                        .with(GS.m_uv_vacuum),
                ))
            })
            .collect()
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
        .shrunken_sub_lmb(&region.support, &contracted_internal, region.boundary.clone())
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

/// A scalar vacuum topology: lines `(a, b, power)` between vertices.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct VacuumTopology {
    pub name: String,
    pub edges: Vec<(usize, usize, u32)>,
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

    /// `Delta(G)` of the massive root by route (ii): massless log-divergent
    /// dotted graphs through `RStarPlan::counterterm`, times
    /// `(mass^2)^(omega/2)`. Returns `(Delta, per-distribution log)`.
    pub fn overall_counterterm_rstar(
        &self,
        model: &Model,
        mass: &Atom,
        vakint: &Vakint,
        uv_settings: &UVgenerationSettings,
    ) -> Result<(Atom, Vec<(Vec<u32>, u64, Atom, u128)>)> {
        let omega = self.superficial_degree();
        if omega < 0 {
            return Ok((Atom::Zero, Vec::new()));
        }
        ensure!(omega % 2 == 0, "odd superficial degree of a scalar vacuum graph");
        let n = u32::try_from(omega / 2)?;
        let mut total = Atom::Zero;
        let mut log = Vec::new();
        for (a, weight) in self.power_distributions(n) {
            let started = Instant::now();
            let dotted = VacuumTopology {
                name: format!("{}_a{}", self.name, a.iter().map(u32::to_string).collect::<String>()),
                edges: self
                    .edges
                    .iter()
                    .zip(&a)
                    .map(|(&(x, y, p), &extra)| (x, y, p + extra))
                    .collect(),
            };
            debug_assert_eq!(dotted.superficial_degree(), 0);
            let massless = dotted.graph(model, &Atom::Zero)?;
            let plan = RStarPlan::with_automatic_edge(&massless, model, uv_settings)?;
            let delta = plan
                .counterterm(&massless, model, vakint, uv_settings)
                .wrap_err_with(|| format!("while computing the R* counterterm of {}", dotted.name))?;
            total += Atom::num(weight as i64) * &delta;
            log.push((a, weight, delta, started.elapsed().as_millis()));
        }
        Ok(((total * mass.pow(Atom::num(omega))).expand(), log))
    }
}

/// `Delta(S)` of a CONNECTED proper child `S` of `graph`, computed on a
/// standalone runtime graph: the lines of `S` (own masses) plus one external
/// leg per boundary half-edge of `S`, all incoming. The leg momenta are then
/// renamed to the momenta entering `S` through those half-edges: `+Q(e)` at
/// the sink end of line `e`, `-Q(e)` at its source end (a self-loop at a
/// vertex of `S` contributes `+Q(e) - Q(e) = 0`). The counterterm of `S`
/// depends on `S` and its incoming momenta only, so this is the same object
/// the in-graph forest computes, without its boundary bookkeeping.
#[allow(clippy::too_many_arguments)]
pub fn standalone_child_counterterm(
    graph: &Graph,
    model: &Model,
    support: &SuBitGraph,
    vakint: &Vakint,
    uv_settings: &UVgenerationSettings,
    forest_extra_terms: i64,
) -> Result<Atom> {
    use linnet::half_edge::{involution::{Flow, HedgePair}, subgraph::SubSetLike as _};
    let mut nodes = Vec::new();
    let mut internal = Vec::new();
    for (pair, edge, data) in graph.iter_edges_of(support) {
        let HedgePair::Paired { source, sink } = pair else {
            return Err(eyre!("a child contains an unpaired half-edge"));
        };
        let (a, b) = (graph.node_id(source), graph.node_id(sink));
        for node in [a, b] {
            if !nodes.contains(&node) {
                nodes.push(node);
            }
        }
        internal.push((edge, a, b, data.data.mass_atom(model)));
    }
    let crown = graph.dummy_less_full_crown(support);
    let mut boundary = Vec::new();
    for hedge in crown.included_iter() {
        let edge = graph[&hedge];
        let node = graph.node_id(hedge);
        ensure!(nodes.contains(&node), "a boundary half-edge outside the child");
        let sign = match graph.flow(hedge) {
            Flow::Sink => 1,
            Flow::Source => -1,
        };
        boundary.push((edge, node, sign));
    }
    ensure!(!boundary.is_empty(), "a proper child of a vacuum graph has a boundary");
    let index = |node| nodes.iter().position(|n| *n == node).unwrap();
    let mut components = (0..nodes.len()).collect::<Vec<_>>();
    let mut slots = vec![0usize; nodes.len()];
    let mut loops = 0;
    let mut dot = String::from(
        "digraph child {\ngraph [projector=1]; node [num=1];\nedge [particle=\"scalar_1\" num=1];\next [style=invis];\n",
    );
    for (id, (_, a, b, _)) in internal.iter().enumerate() {
        let (ia, ib) = (index(*a), index(*b));
        let source = slots[ia];
        slots[ia] += 1;
        let sink = slots[ib];
        slots[ib] += 1;
        let basis = if components[ia] == components[ib] {
            loops += 1;
            format!("lmb_id={}", loops - 1)
        } else {
            let (from, into) = (components[ib], components[ia]);
            for component in &mut components {
                if *component == from {
                    *component = into;
                }
            }
            String::new()
        };
        dot.push_str(&format!(
            "v{ia} -> v{ib} [id={id} {basis} source=\"{{ufo_order:{source}}}\" sink=\"{{ufo_order:{sink}}}\"];\n"
        ));
    }
    let first_external = internal.len();
    for (j, (_, node, _)) in boundary.iter().enumerate() {
        let i = index(*node);
        let sink = slots[i];
        slots[i] += 1;
        dot.push_str(&format!(
            "ext -> v{i} [id={} sink=\"{{ufo_order:{sink}}}\"];\n",
            first_external + j
        ));
    }
    dot.push('}');
    let mut child: Graph = dot
        .into_finalized_runtime_graph(model)
        .wrap_err("while building a standalone child graph")?;
    for (id, (_, _, _, mass)) in internal.iter().enumerate() {
        let data = &mut child.underlying[linnet::half_edge::involution::EdgeIndex::from(id)];
        data.particle = data.particle.clone().override_mass(Some(mass.clone()));
    }
    let plan = RStarPlan::new(
        &child,
        model,
        linnet::half_edge::involution::EdgeIndex::from(0),
        uv_settings,
    )?;
    let support = plan.support.clone();
    let mut values = plan.massive_forest_counterterms(
        &child,
        model,
        vakint,
        uv_settings,
        true,
        forest_extra_terms,
        std::slice::from_ref(&support),
    )?;
    let value = values.remove(&support).unwrap_or_default();
    let renames = boundary
        .iter()
        .enumerate()
        .map(|(j, (edge, _, sign))| {
            symbolica::id::Replacement::new(
                function!(GS.emr_mom, first_external + j, W_.x___).to_pattern(),
                Atom::num(*sign) * function!(GS.emr_mom, usize::from(*edge), W_.x___),
            )
        })
        .collect::<Vec<_>>();
    Ok(value.replace_multiple(&renames).expand())
}
