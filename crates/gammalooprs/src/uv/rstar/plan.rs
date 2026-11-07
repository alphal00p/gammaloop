//! Scoped graph preparation for the UV-counterterm form of the R* operation.
//!
//! All supports and boundary half-edges refer to the original typed graph. In
//! particular, contracting a region never renumbers an insertion's tensor slots.

use std::collections::{BTreeMap, BTreeSet};

use color_eyre::Result;
use eyre::{ensure, eyre};
use linnet::half_edge::{
    HedgeGraph, PowersetIterator,
    involution::{EdgeIndex, HedgePair},
    subgraph::{Inclusion, ModifySubSet, SuBitGraph, SubSetLike, SubSetOps},
};
use symbolica::atom::{Atom, AtomCore};
use vakint::Vakint;

use crate::{
    debug_tags,
    graph::{FeynmanGraph, Graph},
    model::Model,
    utils::GS,
    uv::{
        ApproximationType, RenormalizationPrescriptionSettings, Spinney, UVgenerationSettings,
        UltravioletGraph,
        approx::CutStructure,
        hedge_poset::{Forests, Wood},
    },
};

/// A local polynomial replaces this support, with its original open indices.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub struct InsertionSlot {
    pub support: SuBitGraph,
    pub boundary: SuBitGraph,
}

/// An explicit, possibly decorated, counterterm argument or remaining graph.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub struct GraphRegion {
    /// Original support before the listed contractions.
    pub support: SuBitGraph,
    /// Propagators whose loop integrations have not been consumed.
    pub remaining: SuBitGraph,
    pub insertions: Vec<InsertionSlot>,
    pub boundary: SuBitGraph,
    pub loops: usize,
}

impl GraphRegion {
    fn new(graph: &Graph, support: &SuBitGraph, contractions: &[SuBitGraph]) -> Self {
        let mut remaining = support.clone();
        let insertions = contractions
            .iter()
            .map(|contracted| {
                remaining.subtract_with(contracted);
                InsertionSlot {
                    support: contracted.clone(),
                    boundary: graph.dummy_less_full_crown(contracted),
                }
            })
            .collect();
        let quotient = Self::quotient_topology(graph, contractions);
        Self {
            support: support.clone(),
            loops: quotient.cyclotomatic_number(&remaining),
            remaining,
            insertions,
            boundary: graph.dummy_less_full_crown(support),
        }
    }

    fn quotient_topology(graph: &Graph, contractions: &[SuBitGraph]) -> HedgeGraph<(), (), ()> {
        let mut quotient = graph.underlying.just_structure();
        for contracted in contractions {
            let nodes = quotient
                .iter_nodes_of(contracted)
                .map(|(node, _, _)| node)
                .collect::<Vec<_>>();
            if !nodes.is_empty() {
                // Identify, rather than delete, so original half-edge labels survive.
                quotient.identify_nodes(&nodes, ());
            }
        }
        quotient
    }

    fn is_one_particle_irreducible(&self, graph: &Graph) -> bool {
        // A disconnected ancestor must not appear connected merely because its
        // contracted hard component has become an isolated, invisible vertex.
        if self.loops == 0 || !graph.is_connected(&self.support) {
            return false;
        }
        let quotient = self.contracted_topology(graph);
        quotient.is_connected(&self.remaining) && quotient.bridges_of(&self.remaining).is_empty()
    }

    pub(crate) fn contracted_topology(&self, graph: &Graph) -> HedgeGraph<(), (), ()> {
        let contractions = self
            .insertions
            .iter()
            .map(|slot| slot.support.clone())
            .collect::<Vec<_>>();
        Self::quotient_topology(graph, &contractions)
    }

    pub fn contains_edge(&self, graph: &Graph, edge: EdgeIndex) -> bool {
        self.remaining.includes(&graph[&edge].1)
    }

    /// A swallowed mass kills an unintegrated vacuum loop, not a local term.
    pub(crate) fn is_scaleless(&self, graph: &Graph, marked_edge: EdgeIndex) -> bool {
        self.loops > 0 && !self.contains_edge(graph, marked_edge)
    }

    fn largest_cut_component(&self, graph: &Graph, marked_edge: EdgeIndex) -> usize {
        if self.loops == 0 {
            return 0;
        }
        let mut cut = self.remaining.clone();
        cut.sub(graph[&marked_edge].1);
        let quotient = self.contracted_topology(graph);
        let cyclic = quotient.non_bridges_of(&cut);
        quotient
            .connected_components(&cyclic)
            .iter()
            .map(|component| quotient.cyclotomatic_number(component))
            .max()
            .unwrap_or(0)
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub enum CountertermScope {
    /// Child coefficients belong inside this prepared amplitude's own -K.
    PreparedAmplitude,
    /// Evaluate each inserted monomial; its coefficient stays outside this -K.
    InsertionMonomial,
}

/// One UV counterterm call, including its own integration and pole projection.
///
/// The subtraction sign -K belongs to the call, independently of an IR forest's
/// combinatorial sign. No projector from an enclosing graph is included here.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub struct CountertermCall {
    pub argument: GraphRegion,
    /// Upper bound from the original typed graph. Insertions can lower it.
    pub taylor_order: i32,
    pub scope: CountertermScope,
}

impl CountertermCall {
    fn original(graph: &Graph, support: &SuBitGraph, taylor_order: i32) -> Self {
        Self {
            argument: GraphRegion::new(graph, support, &[]),
            taylor_order,
            scope: CountertermScope::PreparedAmplitude,
        }
    }

    fn quotient(graph: &Graph, support: &SuBitGraph, contracted: &SuBitGraph) -> Self {
        Self {
            argument: GraphRegion::new(graph, support, std::slice::from_ref(contracted)),
            taylor_order: graph.local_dod(support),
            scope: CountertermScope::InsertionMonomial,
        }
    }
}

/// One spinney contribution: independent child calls inserted into a remainder.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct UvPreparationTerm {
    pub remainder: GraphRegion,
    pub counterterms: Vec<CountertermCall>,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct IrForestTerm {
    /// -(-1)^|F|; each step also carries its own UV counterterm sign.
    pub sign: i8,
    /// Smaller ancestors act first. The last call contracts the whole quotient.
    /// Coefficients of the current insertion remain outside each next call.
    pub steps: Vec<CountertermCall>,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct IrCorrection {
    pub hard: SuBitGraph,
    /// Total momentum degree, including forced external derivative factors.
    pub taylor_order: i32,
    /// Full R_U preparation, including the hard component's overall counterterm.
    /// Remaining hard loops stay unintegrated until the root Taylor operation.
    pub preparation: Vec<UvPreparationTerm>,
    pub forest: Vec<IrForestTerm>,
}

/// Structural work before algebraic zeros, tensor moments and backend reduction.
/// These counts predict preparation size, not elapsed integration time.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct RStarCost {
    /// Largest edge-dependent reduction: a cut loop component or an IR UV map.
    /// Proper UV child loop counts do not distinguish the marked lines.
    pub max_reduction_loops: usize,
    pub ordinary_terms: usize,
    pub ir_terms: usize,
    /// UV operations across surviving IR summands, before cache reuse.
    pub quotient_calls: usize,
    /// Degree-weighted hard preparations; tensor ranks can increase actual work.
    pub hard_taylor_work: usize,
    /// Cutting outside triangles preserves potential propagator reductions.
    /// This does not certify that a backend triangle rule can be applied.
    pub cuts_triangle: bool,
}

impl RStarCost {
    pub fn terms(&self) -> usize {
        self.ordinary_terms + self.ir_terms
    }

    fn selection_key(&self) -> (usize, usize, usize, usize, bool) {
        (
            self.max_reduction_loops,
            self.terms(),
            self.quotient_calls,
            self.hard_taylor_work,
            self.cuts_triangle,
        )
    }

    fn estimate(
        graph: &Graph,
        support: &SuBitGraph,
        marked_edge: EdgeIndex,
        ordinary: &[UvPreparationTerm],
        corrections: &[IrCorrection],
    ) -> Self {
        let mut cost = Self {
            max_reduction_loops: ordinary
                .iter()
                .filter(|term| !term.remainder.is_scaleless(graph, marked_edge))
                .map(|term| term.remainder.largest_cut_component(graph, marked_edge))
                .max()
                .unwrap_or(0),
            ordinary_terms: ordinary
                .iter()
                .filter(|term| !term.remainder.is_scaleless(graph, marked_edge))
                .count(),
            cuts_triangle: Self::cuts_triangle(graph, support, marked_edge),
            ..Self::default()
        };
        for correction in corrections
            .iter()
            .filter(|correction| correction.hard.includes(&graph[&marked_edge].1))
        {
            let preparations = correction
                .preparation
                .iter()
                .filter(|term| !term.remainder.is_scaleless(graph, marked_edge))
                .collect::<Vec<_>>();
            let surviving = preparations.len();
            if surviving > 0 {
                for term in preparations {
                    cost.max_reduction_loops = cost
                        .max_reduction_loops
                        .max(term.remainder.largest_cut_component(graph, marked_edge));
                }
                for call in correction.forest.iter().flat_map(|term| &term.steps) {
                    cost.max_reduction_loops = cost.max_reduction_loops.max(call.argument.loops);
                }
            }
            cost.ir_terms += surviving * correction.forest.len();
            cost.quotient_calls += surviving
                * correction
                    .forest
                    .iter()
                    .map(|term| term.steps.len())
                    .sum::<usize>();
            cost.hard_taylor_work += surviving * (correction.taylor_order.max(0) as usize + 1);
        }
        cost
    }

    fn cuts_triangle(graph: &Graph, support: &SuBitGraph, edge: EdgeIndex) -> bool {
        let HedgePair::Paired { source, sink } = graph[&edge].1 else {
            return false;
        };
        let source = graph.node_id(source);
        let sink = graph.node_id(sink);
        if source == sink {
            return false;
        }
        let neighbors = |node| {
            graph
                .iter_edges_of(support)
                .filter_map(|(pair, _, _)| {
                    let HedgePair::Paired { source: a, sink: b } = pair else {
                        return None;
                    };
                    let a = graph.node_id(a);
                    let b = graph.node_id(b);
                    let other = if a == node {
                        b
                    } else if b == node {
                        a
                    } else {
                        return None;
                    };
                    (other != source && other != sink).then_some(other)
                })
                .collect::<BTreeSet<_>>()
        };
        !neighbors(source).is_disjoint(&neighbors(sink))
    }
}

/// Topology of -K I_rem T_compensated [R'_U + sum_gamma C^U B_gamma] G.
///
/// This catalogue is independent of Feynman-rule evaluation. Hard coefficients
/// remain outside every quotient call, while the root -K encloses the complete
/// assembled expression. The ordinary and hard preparations deliberately retain
/// contractions that swallow the marked line: only an actual remaining scaleless
/// loop integral may be discarded, never an already local mixed contribution.
pub struct RStarPlan {
    pub support: SuBitGraph,
    pub marked_edge: EdgeIndex,
    /// Total degree of the local polynomial, including external vertex momenta.
    /// This excludes external projectors, as does the graph numerator.
    pub taylor_order: i32,
    pub ordinary_uv: Vec<UvPreparationTerm>,
    pub ir: Vec<IrCorrection>,
    /// Undecorated UV children, ordered from fewer loops to more loops.
    /// Decorated quotient dependencies have their explicit order in `forest`.
    pub counterterms: Vec<CountertermCall>,
    catalogue: Vec<Spinney>,
}

impl RStarPlan {
    pub fn new(
        graph: &Graph,
        model: &Model,
        marked_edge: EdgeIndex,
        settings: &UVgenerationSettings,
    ) -> Result<Self> {
        Self::build(graph, model, Some(marked_edge), settings)
    }

    /// Prefer smaller loop reductions, then fewer surviving graph terms,
    /// quotient calls and hard Taylor work, then a cut outside every triangle. Remaining ties
    /// use the original edge index. Classification is shared
    /// across candidates, and each hard component's forest is generated once.
    /// This graph-only heuristic does not predict numerator zeros or Forcer's
    /// reduction difficulty; an explicit mark remains available through `new`.
    pub fn with_automatic_edge(
        graph: &Graph,
        model: &Model,
        settings: &UVgenerationSettings,
    ) -> Result<Self> {
        Self::build(graph, model, None, settings)
    }

    pub fn estimated_cost(&self, graph: &Graph) -> RStarCost {
        if graph.compute_dod(&self.support) < 0 {
            return RStarCost::default();
        }
        RStarCost::estimate(
            graph,
            &self.support,
            self.marked_edge,
            &self.ordinary_uv,
            &self.ir,
        )
    }

    fn build(
        graph: &Graph,
        model: &Model,
        marked_edge: Option<EdgeIndex>,
        settings: &UVgenerationSettings,
    ) -> Result<Self> {
        let support = graph
            .no_dummy()
            .subtract(&graph.external_filter::<SuBitGraph>());
        ensure!(
            marked_edge.is_none_or(|marked| graph
                .iter_edges_of(&support)
                .any(|(pair, edge, _)| { pair.is_paired() && edge == marked })),
            "the R* marked edge must be an internal, non-dummy propagator"
        );
        ensure!(
            graph.is_connected(&support)
                && graph.n_loops(&support) > 0
                && graph.bridges_of(&support).is_empty(),
            "an R* counterterm task requires a connected 1PI graph with loops"
        );
        let overall_divergent = graph.compute_dod(&support) >= 0;
        // R* scales the complete boundary numerator. Forced external momenta
        // count toward this Taylor degree even though they do not contribute
        // to loop UV power counting (for example the antighost derivative).
        let taylor_order = graph.local_dod(&support);

        // R* computes the MS UV map. Reuse the established spinney classification
        // and routing, without importing MUV finite terms or user exclusions.
        let ms_settings = UVgenerationSettings {
            subtract_uv: true,
            renormalization_prescription: RenormalizationPrescriptionSettings {
                log_divergent: ApproximationType::PolePart,
                massive_power_divergent: ApproximationType::PolePart,
                massless_power_divergent: ApproximationType::PolePart,
                overrides: Vec::new(),
            },
            ..settings.clone()
        };
        let mut catalogue =
            graph.classified_spinneys(&support, model, &ms_settings, &graph.loop_momentum_basis);
        catalogue.sort_by(|left, right| {
            left.filter()
                .n_included()
                .cmp(&right.filter().n_included())
                .then_with(|| left.filter().cmp(right.filter()))
        });
        ensure!(
            !overall_divergent || catalogue.iter().any(|spinney| spinney.filter() == &support),
            "the R* graph has no UV spinney with a compatible loop momentum basis"
        );

        let ordinary_uv = Self::preparation(graph, &support, &catalogue, false);
        // An overall convergent graph can still have subdivergence poles.
        // Its local UV map vanishes, so only the physical UV preparation is
        // needed for reconstructing its bare pole part.
        let mut ir = if overall_divergent {
            Self::ir_corrections(graph, &support, marked_edge, &catalogue)?
        } else {
            Vec::new()
        };
        let marked_edge = if let Some(marked_edge) = marked_edge {
            marked_edge
        } else {
            graph
                .iter_edges_of(&support)
                .filter(|(pair, _, _)| pair.is_paired())
                .map(|(_, edge, _)| {
                    let cost = if overall_divergent {
                        RStarCost::estimate(graph, &support, edge, &ordinary_uv, &ir)
                    } else {
                        RStarCost::default()
                    };
                    debug_tags!(#uv, #integrated, #rstar;
                        stage = "massive_edge_candidate",
                        marked_edge = usize::from(edge),
                        max_reduction_loops = cost.max_reduction_loops,
                        ordinary_terms = cost.ordinary_terms,
                        ir_terms = cost.ir_terms,
                        quotient_calls = cost.quotient_calls,
                        hard_taylor_work = cost.hard_taylor_work,
                        cuts_triangle = cost.cuts_triangle,
                        "Estimated R* work for an internal massive line"
                    );
                    (cost.selection_key(), edge)
                })
                .min()
                .ok_or_else(|| eyre!("the R* graph has no internal propagator"))?
                .1
        };
        ir.retain(|correction| correction.hard.includes(&graph[&marked_edge].1));
        let mut counterterms = ordinary_uv
            .iter()
            .chain(ir.iter().flat_map(|correction| &correction.preparation))
            .flat_map(|term| term.counterterms.iter().cloned())
            .collect::<BTreeSet<_>>()
            .into_iter()
            .collect::<Vec<_>>();
        counterterms.sort_by(|left, right| {
            left.argument
                .loops
                .cmp(&right.argument.loops)
                .then_with(|| left.cmp(right))
        });
        Ok(Self {
            support,
            marked_edge,
            taylor_order,
            ordinary_uv,
            ir,
            counterterms,
            catalogue,
        })
    }

    /// Reuse the existing UV forest for evaluations of the original UV map.
    /// Its all-mass Taylor kernel is not the R* root's single-line kernel.
    pub fn uv_forest(&self, graph: &Graph, settings: &UVgenerationSettings) -> Forests {
        Wood::from_spinneys(
            self.catalogue.clone(),
            graph,
            CutStructure::empty(graph),
            &settings.vakint,
        )
        .unfold()
    }

    /// Evaluate only proper UV children with the existing massive R operation.
    /// `reference_mass` is a mass, rather than a mass squared. The returned
    /// polynomials retain its powers and their original open tensor slots.
    /// The N-loop root is excluded before forest construction, so this never
    /// computes the massive vacuum integral that the R* rearrangement replaces.
    pub fn evaluate_uv_children(
        &self,
        graph: &Graph,
        model: &Model,
        reference_mass: &Atom,
        vakint: &Vakint,
        settings: &UVgenerationSettings,
    ) -> Result<BTreeMap<SuBitGraph, Atom>> {
        let mut family = graph.clone();
        family.underlying[self.marked_edge].particle = family[self.marked_edge]
            .particle
            .clone()
            .override_mass(Some(reference_mass.clone()));
        let settings = UVgenerationSettings {
            generate_integrated: true,
            // Diagnostic wrappers are not part of an inserted polynomial.
            add_marker: false,
            ..settings.clone()
        };
        let mut forest = Wood::from_spinneys(
            self.catalogue
                .iter()
                .filter(|spinney| spinney.filter() != &self.support)
                .cloned(),
            &family,
            CutStructure::empty(&family),
            &settings.vakint,
        )
        .unfold();
        forest.integrate(&family, model, vakint, &settings)?;
        self.counterterms
            .iter()
            .map(|call| {
                let support = &call.argument.support;
                let polynomial = forest
                    .counterterm_for_support(support)?
                    .replace(GS.m_uv_expansion)
                    .with(GS.m_uv_vacuum);
                Ok((support.clone(), polynomial))
            })
            .collect()
    }

    fn preparation(
        graph: &Graph,
        support: &SuBitGraph,
        catalogue: &[Spinney],
        include_overall: bool,
    ) -> Vec<UvPreparationTerm> {
        catalogue
            .iter()
            .filter(|spinney| {
                support.includes(spinney.filter())
                    && (include_overall || spinney.filter() != support)
            })
            .map(|spinney| {
                let mut components = graph.connected_components(spinney.filter());
                components.sort();
                UvPreparationTerm {
                    remainder: GraphRegion::new(graph, support, &components),
                    counterterms: components
                        .iter()
                        .map(|component| {
                            CountertermCall::original(
                                graph,
                                component,
                                graph.compute_dod(component),
                            )
                        })
                        .collect(),
                }
            })
            .collect()
    }

    /// The IR kernel is a sum over strict connected ancestor chains, separate
    /// from ordinary UV spinneys. Candidates are tested on the quotient, not
    /// restricted to the UV catalogue of the original graph.
    fn ir_corrections(
        graph: &Graph,
        support: &SuBitGraph,
        marked_edge: Option<EdgeIndex>,
        catalogue: &[Spinney],
    ) -> Result<Vec<IrCorrection>> {
        let mut corrections = Vec::new();
        for hard in catalogue.iter().filter(|spinney| {
            spinney.n_components() == 1
                && spinney.filter() != support
                && marked_edge.is_none_or(|edge| spinney.filter().includes(&graph[&edge].1))
                && graph.bridges_of(spinney.filter()).is_empty()
        }) {
            let mut forest = Vec::new();
            Self::ancestor_chains(graph, support, hard.filter(), &mut Vec::new(), &mut forest)?;
            corrections.push(IrCorrection {
                hard: hard.filter().clone(),
                taylor_order: graph.local_dod(hard.filter()),
                preparation: Self::preparation(graph, hard.filter(), catalogue, true),
                forest,
            });
        }
        Ok(corrections)
    }

    fn ancestor_chains(
        graph: &Graph,
        full: &SuBitGraph,
        contracted: &SuBitGraph,
        prefix: &mut Vec<CountertermCall>,
        terms: &mut Vec<IrForestTerm>,
    ) -> Result<()> {
        let overall = CountertermCall::quotient(graph, full, contracted);
        if overall.argument.is_one_particle_irreducible(graph) {
            let mut steps = prefix.clone();
            steps.push(overall);
            terms.push(IrForestTerm {
                sign: if prefix.len().is_multiple_of(2) {
                    -1
                } else {
                    1
                },
                steps,
            });
        }

        let available = graph
            .iter_edges_of(&full.subtract(contracted))
            .filter_map(|(pair, _, _)| pair.is_paired().then_some(pair))
            .collect::<Vec<_>>();
        ensure!(
            available.len() < usize::BITS as usize,
            "the R* ancestor powerset has too many edges to enumerate"
        );
        let edge_count = u8::try_from(available.len())
            .map_err(|_| eyre!("the R* ancestor edge count exceeds the powerset limit"))?;
        for subset in PowersetIterator::<usize>::new(edge_count).skip(1) {
            let mut ancestor = contracted.clone();
            for edge in subset.included_iter() {
                ancestor.add(available[edge]);
            }
            // Contracting G here would leave a repeated overall c-operation.
            // An insertion has total degree at most local_dod(contracted),
            // including inherited masses and forced external momenta. Thus
            // the decorated quotient's degree is at most local_dod(ancestor).
            // A negative bound makes this UV map, and its entire branch, zero.
            // Use the full local degree, not original loop power counting:
            // forced external momenta also belong to the inserted polynomial.
            if &ancestor == full || graph.local_dod(&ancestor) < 0 {
                continue;
            }
            let call = CountertermCall::quotient(graph, &ancestor, contracted);
            if !call.argument.is_one_particle_irreducible(graph) {
                continue;
            }
            prefix.push(call);
            Self::ancestor_chains(graph, full, &ancestor, prefix, terms)?;
            prefix.pop();
        }
        Ok(())
    }
}

#[cfg(test)]
pub(super) mod tests {
    use super::*;
    use crate::{
        finalized_runtime_dot, graph::parse::IntoFinalizedRuntimeGraph,
        initialisation::test_initialise, utils::load_generic_model,
    };

    fn edges(graph: &Graph, selected: &[usize]) -> SuBitGraph {
        let mut support = graph.empty_subgraph::<SuBitGraph>();
        for edge in selected {
            support.add(graph[&EdgeIndex::from(*edge)].1);
        }
        support
    }

    pub(in crate::uv::rstar) fn logarithmic_graph() -> Result<Graph> {
        let mut graph: Graph = finalized_runtime_dot!(
            digraph G {
                graph [projector=1];
                node [num=1];
                edge [particle="scalar_1" num=1];
                A -> B [id=0 source="{ufo_order:0}" sink="{ufo_order:0}"];
                A -> B [id=1 lmb_id=0 source="{ufo_order:1}" sink="{ufo_order:1}"];
                A -> C [id=2 source="{ufo_order:2}" sink="{ufo_order:0}"];
                C -> B [id=3 lmb_id=1 source="{ufo_order:1}" sink="{ufo_order:2}"];
            }, "scalars"
        )?;
        let edges = graph
            .iter_edges()
            .map(|(_, edge, _)| edge)
            .collect::<Vec<_>>();
        for edge in edges {
            graph.underlying[edge].particle =
                graph[edge].particle.clone().override_mass(Some(Atom::Zero));
        }
        Ok(graph)
    }

    #[test]
    fn hard_preparation_keeps_local_mixed_term_and_original_boundary() -> Result<()> {
        test_initialise()?;
        let graph = logarithmic_graph()?;
        let model = load_generic_model("scalars");
        let plan = RStarPlan::new(&graph, &model, EdgeIndex::from(0), &Default::default())?;
        assert_eq!(plan.taylor_order, 0);
        assert_eq!(plan.ordinary_uv.len(), 2);
        assert_eq!(plan.ir.len(), 1);
        let correction = &plan.ir[0];
        assert_eq!(correction.hard, edges(&graph, &[0, 1]));
        assert_eq!(correction.taylor_order, 0);
        assert_eq!(correction.preparation.len(), 2);
        let local = correction
            .preparation
            .iter()
            .find(|term| term.remainder.loops == 0)
            .expect("full R_U includes the local hard counterterm");
        assert!(local.remainder.remaining.is_empty());
        assert_eq!(local.counterterms.len(), 1);
        assert_eq!(local.remainder.insertions[0].support, correction.hard);
        assert_eq!(
            local.remainder.insertions[0].boundary,
            graph.dummy_less_full_crown(&correction.hard)
        );
        assert_eq!(correction.forest.len(), 1);
        assert_eq!(correction.forest[0].sign, -1);
        let steps = &correction.forest[0].steps;
        assert_eq!(steps.len(), 1);
        assert_eq!(steps[0].scope, CountertermScope::InsertionMonomial);
        assert_eq!(steps[0].argument.remaining, edges(&graph, &[2, 3]));
        assert_eq!(steps[0].argument.loops, 1);
        // Before contraction these propagators form a path, not a 1PI graph.
        assert!(!graph.bridges_of(&steps[0].argument.remaining).is_empty());
        assert!(steps[0].argument.is_one_particle_irreducible(&graph));
        assert!(
            !plan.ordinary_uv[1]
                .remainder
                .contains_edge(&graph, plan.marked_edge)
        );

        let alternate = RStarPlan::new(&graph, &model, EdgeIndex::from(2), &Default::default())?;
        assert!(alternate.ir.is_empty());
        assert!(
            alternate.ordinary_uv[1]
                .remainder
                .contains_edge(&graph, alternate.marked_edge)
        );
        assert_eq!(plan.counterterms.len(), 1);
        assert_eq!(plan.counterterms[0].argument.loops, 1);
        assert!(plan.uv_forest(&graph, &Default::default()).graph.n_nodes() > 1);
        Ok(())
    }

    #[test]
    fn automatic_edge_avoids_ir_forests_and_counts_local_mixed_terms() -> Result<()> {
        test_initialise()?;
        let graph = logarithmic_graph()?;
        let model = load_generic_model("scalars");
        let settings = UVgenerationSettings::default();
        let automatic = RStarPlan::with_automatic_edge(&graph, &model, &settings)?;
        assert_eq!(automatic.marked_edge, EdgeIndex::from(2));
        assert!(automatic.ir.is_empty());
        assert_eq!(automatic.estimated_cost(&graph).terms(), 2);
        for (_, edge, _) in graph.iter_edges_of(&automatic.support) {
            let explicit = RStarPlan::new(&graph, &model, edge, &settings)?;
            let cost = explicit.estimated_cost(&graph);
            assert!(automatic.estimated_cost(&graph).selection_key() <= cost.selection_key());
            if edge == EdgeIndex::from(0) {
                // One ordinary term is scaleless, but BOTH hard-preparation
                // terms survive: the bubble and its local counterterm.
                assert_eq!(cost.ordinary_terms, 1);
                assert_eq!(cost.ir_terms, 2);
                assert_eq!(cost.quotient_calls, 2);
            }
            if edge == automatic.marked_edge {
                assert_eq!(explicit.ordinary_uv, automatic.ordinary_uv);
                assert_eq!(explicit.ir, automatic.ir);
                assert_eq!(explicit.counterterms, automatic.counterterms);
            }
        }
        Ok(())
    }

    #[test]
    fn automatic_edge_ranks_unavoidable_ir_corrections_deterministically() -> Result<()> {
        test_initialise()?;
        let graph = finalized_runtime_dot!(
            digraph G {
                graph [projector=1];
                node [num=1];
                edge [particle="scalar_1" num=1];
                A -> B [id=0 source="{ufo_order:0}" sink="{ufo_order:0}"];
                A -> B [id=1 lmb_id=0 source="{ufo_order:1}" sink="{ufo_order:1}"];
                A -> B [id=2 lmb_id=1 source="{ufo_order:2}" sink="{ufo_order:2}"];
            }, "scalars"
        )?;
        let model = load_generic_model("scalars");
        let settings = UVgenerationSettings::default();
        let automatic = RStarPlan::with_automatic_edge(&graph, &model, &settings)?;
        assert_eq!(automatic.marked_edge, EdgeIndex::from(0));
        assert!(!automatic.ir.is_empty());
        for (_, edge, _) in graph.iter_edges_of(&automatic.support) {
            let explicit = RStarPlan::new(&graph, &model, edge, &settings)?;
            assert_eq!(
                explicit.estimated_cost(&graph),
                automatic.estimated_cost(&graph)
            );
        }
        Ok(())
    }

    #[test]
    fn triangle_tie_break_ignores_orientation_parallel_edges_and_self_loops() -> Result<()> {
        test_initialise()?;
        let graph = finalized_runtime_dot!(
            digraph G {
                graph [projector=1]; node [num=1];
                edge [particle="scalar_1" num=1];
                A -> B [id=0 source="{ufo_order:0}" sink="{ufo_order:0}"];
                B -> A [id=1 lmb_id=0 source="{ufo_order:1}" sink="{ufo_order:1}"];
                A -> C [id=2 source="{ufo_order:2}" sink="{ufo_order:0}"];
                C -> B [id=3 lmb_id=1 source="{ufo_order:1}" sink="{ufo_order:2}"];
                A -> A [id=4 lmb_id=2 source="{ufo_order:3}" sink="{ufo_order:4}"];
            }, "scalars"
        )?;
        let full = edges(&graph, &[0, 1, 2, 3, 4]);
        for edge in [0, 1, 2, 3].map(EdgeIndex::from) {
            assert!(RStarCost::cuts_triangle(&graph, &full, edge));
        }
        assert!(!RStarCost::cuts_triangle(&graph, &full, EdgeIndex::from(4)));
        let bubble = edges(&graph, &[0, 1, 4]);
        assert!(!RStarCost::cuts_triangle(
            &graph,
            &bubble,
            EdgeIndex::from(0)
        ));
        Ok(())
    }

    #[test]
    fn restricted_ancestors_are_coefficientwise_and_smallest_first() -> Result<()> {
        test_initialise()?;
        // Four parallel edges produce two strict, incomparable one-edge
        // ancestors of the hard bubble. Each remainder becomes a self-loop.
        let graph = finalized_runtime_dot!(
            digraph G {
                graph [projector=1];
                node [num=1];
                edge [particle="scalar_1" num=1];
                A -> B [id=0 source="{ufo_order:0}" sink="{ufo_order:0}"];
                A -> B [id=1 lmb_id=0 source="{ufo_order:1}" sink="{ufo_order:1}"];
                A -> B [id=2 lmb_id=1 source="{ufo_order:2}" sink="{ufo_order:2}"];
                A -> B [id=3 lmb_id=2 source="{ufo_order:3}" sink="{ufo_order:3}"];
            }, "scalars"
        )?;
        let full = edges(&graph, &[0, 1, 2, 3]);
        let hard = edges(&graph, &[0, 1]);
        let mut terms = Vec::new();
        RStarPlan::ancestor_chains(&graph, &full, &hard, &mut Vec::new(), &mut terms)?;
        assert_eq!(terms.len(), 3);
        assert_eq!(terms[0].sign, -1);
        assert_eq!(terms[0].steps.len(), 1);
        for term in &terms[1..] {
            assert_eq!(term.sign, 1);
            assert_eq!(term.steps.len(), 2);
            let first = &term.steps[0];
            let second = &term.steps[1];
            assert_eq!(first.argument.insertions[0].support, hard);
            assert_eq!(first.argument.remaining.n_included(), 2);
            assert_eq!(
                second.argument.insertions[0].support,
                first.argument.support
            );
            assert_eq!(second.argument.support, full);
            assert_eq!(second.argument.remaining.n_included(), 2);
            assert!(term.steps.iter().all(|call| {
                call.scope == CountertermScope::InsertionMonomial && call.argument.loops > 0
            }));
        }
        Ok(())
    }

    #[test]
    fn ir_forests_skip_negative_degree_ancestors_before_expansion() -> Result<()> {
        test_initialise()?;
        let graph = finalized_runtime_dot!(
            digraph G {
                graph [projector=1]; node [num=1];
                edge [particle="scalar_1" num=1];
                A -> B [id=0 source="{ufo_order:0}" sink="{ufo_order:0}"];
                A -> B [id=1 lmb_id=0 source="{ufo_order:1}" sink="{ufo_order:1}"];
                A -> C [id=2 source="{ufo_order:2}" sink="{ufo_order:0}"];
                C -> D [id=3 source="{ufo_order:1}" sink="{ufo_order:0}"];
                D -> B [id=4 lmb_id=1 source="{ufo_order:1}" sink="{ufo_order:2}"];
                C -> B [id=5 lmb_id=2 source="{ufo_order:2}" sink="{ufo_order:3}"];
            }, "scalars"
        )?;
        let hard = edges(&graph, &[0, 1]);
        let convergent = edges(&graph, &[0, 1, 2, 3, 4]);
        let divergent = edges(&graph, &[0, 1, 2, 5]);
        let full = edges(&graph, &[0, 1, 2, 3, 4, 5]);
        assert_eq!(graph.local_dod(&convergent), -2);
        assert_eq!(graph.local_dod(&divergent), 0);
        assert!(
            GraphRegion::new(&graph, &convergent, std::slice::from_ref(&hard))
                .is_one_particle_irreducible(&graph)
        );
        let mut terms = Vec::new();
        RStarPlan::ancestor_chains(&graph, &full, &hard, &mut Vec::new(), &mut terms)?;
        let calls = terms
            .iter()
            .flat_map(|term| &term.steps)
            .collect::<Vec<_>>();
        assert!(calls.iter().all(|call| call.argument.support != convergent));
        assert!(calls.iter().any(|call| call.argument.support == divergent));
        assert!(calls.iter().all(|call| call.taylor_order >= 0));
        Ok(())
    }

    #[test]
    fn uv_spinneys_reject_a_convergent_component_despite_nonnegative_total_degree() -> Result<()> {
        test_initialise()?;
        // A quark bubble with nested gluon chords. The vertex-disjoint
        // one-loop qqBB box and quark self-energy have degrees -1 and +1.
        // Linear quark numerators retain precisely this power counting.
        let graph = finalized_runtime_dot!(
            digraph G {
                graph [projector=1];
                node [num=1];
                edge [particle="g" num=1 dod=-2];
                ext [style=invis];
                ext -> v1 [id=0 sink="{ufo_order:2}"];
                v0 -> ext [id=1 source="{ufo_order:2}"];
                v1 -> v0 [id=2 particle="d" num="Q(2,spenso::cind(1))" dod=-1 lmb_id=0 source="{ufo_order:0}" sink="{ufo_order:1}"];
                v0 -> v5 [id=3 particle="d" num="Q(3,spenso::cind(1))" dod=-1 source="{ufo_order:0}" sink="{ufo_order:1}"];
                v4 -> v1 [id=4 particle="d" num="Q(4,spenso::cind(1))" dod=-1 source="{ufo_order:0}" sink="{ufo_order:1}"];
                v2 -> v3 [id=5 particle="d" num="Q(5,spenso::cind(1))" dod=-1 lmb_id=1 source="{ufo_order:0}" sink="{ufo_order:1}"];
                v2 -> v3 [id=6 source="{ufo_order:2}" sink="{ufo_order:2}"];
                v5 -> v2 [id=7 particle="d" num="Q(7,spenso::cind(1))" dod=-1 source="{ufo_order:0}" sink="{ufo_order:1}"];
                v3 -> v4 [id=8 particle="d" num="Q(8,spenso::cind(1))" dod=-1 lmb_id=2 source="{ufo_order:0}" sink="{ufo_order:1}"];
                v4 -> v5 [id=9 source="{ufo_order:2}" sink="{ufo_order:2}"];
            }, "sm"
        )?;
        let convergent = edges(&graph, &[2, 3, 4, 9]);
        let divergent = edges(&graph, &[5, 6]);
        let mixed = convergent.union(&divergent);
        let support = edges(&graph, &[2, 3, 4, 5, 6, 7, 8, 9]);
        assert_eq!(graph.compute_dod(&convergent), -1);
        assert_eq!(graph.compute_dod(&divergent), 1);
        assert_eq!(graph.compute_dod(&mixed), 0);
        assert_eq!(graph.connected_components(&mixed).len(), 2);
        let candidate = graph
            .spinneys(&support)
            .into_iter()
            .find(|spinney| spinney.filter == mixed)
            .expect("the summed degree admits the mixed cycle union before classification");
        assert!(Spinney::new(candidate, &graph, &graph.loop_momentum_basis).is_none());

        let model = load_generic_model("sm");
        let plan = RStarPlan::new(&graph, &model, EdgeIndex::from(7), &Default::default())?;
        assert!(
            plan.catalogue
                .iter()
                .all(|spinney| spinney.filter() != &mixed)
        );
        assert!(
            plan.counterterms
                .iter()
                .all(|call| { call.taylor_order >= 0 && call.argument.support != convergent })
        );
        assert!(
            plan.counterterms
                .iter()
                .any(|call| call.argument.support == divergent)
        );
        assert!(
            plan.counterterms
                .iter()
                .any(|call| { call.argument.support == edges(&graph, &[5, 6, 7, 8, 9]) })
        );
        Ok(())
    }

    #[test]
    fn disjoint_uv_children_keep_independent_projections() -> Result<()> {
        test_initialise()?;
        let graph = finalized_runtime_dot!(
            digraph G {
                graph [projector=1];
                node [num=1];
                edge [particle="scalar_1" num=1];
                A -> B [id=0 source="{ufo_order:0}" sink="{ufo_order:0}"];
                A -> B [id=1 lmb_id=0 source="{ufo_order:1}" sink="{ufo_order:1}"];
                C -> D [id=2 source="{ufo_order:0}" sink="{ufo_order:0}"];
                C -> D [id=3 lmb_id=1 source="{ufo_order:1}" sink="{ufo_order:1}"];
                B -> C [id=4 source="{ufo_order:2}" sink="{ufo_order:2}"];
                D -> A [id=5 lmb_id=2 source="{ufo_order:2}" sink="{ufo_order:2}"];
            }, "scalars"
        )?;
        let model = load_generic_model("scalars");
        let plan = RStarPlan::new(&graph, &model, EdgeIndex::from(0), &Default::default())?;
        let disjoint = plan
            .ordinary_uv
            .iter()
            .find(|term| term.counterterms.len() == 2)
            .expect("the two vertex-disjoint bubbles form a spinney");
        assert_eq!(disjoint.remainder.insertions.len(), 2);
        assert_eq!(disjoint.remainder.loops, 1);
        assert!(disjoint.counterterms.iter().all(|call| {
            call.argument.insertions.is_empty()
                && call.scope == CountertermScope::PreparedAmplitude
                && call.argument.loops == 1
        }));
        let disconnected = edges(&graph, &[0, 1, 2, 3]);
        assert!(plan.ir.iter().flat_map(|term| &term.forest).all(|term| {
            term.steps
                .iter()
                .all(|call| call.argument.support != disconnected)
        }));
        assert!(
            plan.counterterms
                .windows(2)
                .all(|calls| { calls[0].argument.loops <= calls[1].argument.loops })
        );
        let automatic = RStarPlan::with_automatic_edge(&graph, &model, &Default::default())?;
        assert_eq!(automatic.marked_edge, EdgeIndex::from(4));
        // Cutting a joining line separates two one-loop p-integrals. Cutting
        // inside a bubble leaves a two-loop block and a two-loop IR UV map.
        assert_eq!(automatic.estimated_cost(&graph).max_reduction_loops, 1);
        assert_eq!(plan.estimated_cost(&graph).max_reduction_loops, 2);
        Ok(())
    }

    #[test]
    fn unfolded_example_keeps_original_uv_catalogue_for_both_marked_lines() -> Result<()> {
        test_initialise()?;
        // The five-edge topology in the unfolded formula. Linear vertex
        // numerators give its original cubic-theory bounds without imposing a
        // projector or performing any tensor contraction in the graph plan.
        let graph = finalized_runtime_dot!(
            digraph G {
                graph [projector=1];
                edge [particle="scalar_1" num=1];
                e [style=invis];
                A [num="Q(0,spenso::cind(1))" dod=1];
                B [num="Q(0,spenso::cind(2))" dod=1];
                C [num="Q(2,spenso::cind(3))+2*Q(3,spenso::cind(3))" dod=1];
                D [num="Q(2,spenso::cind(4))+3*Q(3,spenso::cind(4))" dod=1];
                A -> B [id=0 source="{ufo_order:0}" sink="{ufo_order:0}"];
                A -> C [id=1 source="{ufo_order:1}" sink="{ufo_order:0}"];
                C -> D [id=2 source="{ufo_order:1}" sink="{ufo_order:0}"];
                C -> D [id=3 lmb_id=0 source="{ufo_order:2}" sink="{ufo_order:1}"];
                D -> B [id=4 lmb_id=1 source="{ufo_order:2}" sink="{ufo_order:1}"];
                e -> A [id=5 sink="{ufo_order:2}"];
                B -> e [id=6 source="{ufo_order:2}"];
            }, "scalars"
        )?;
        let model = load_generic_model("scalars");
        let outside = RStarPlan::new(&graph, &model, EdgeIndex::from(1), &Default::default())?;
        let inside = RStarPlan::new(&graph, &model, EdgeIndex::from(2), &Default::default())?;
        assert_eq!(outside.taylor_order, 2);
        assert_eq!(outside.ordinary_uv, inside.ordinary_uv);
        assert_eq!(outside.ordinary_uv.len(), 4);
        let hard_supports = |plan: &RStarPlan| {
            plan.ir
                .iter()
                .map(|term| term.hard.clone())
                .collect::<BTreeSet<_>>()
        };
        assert_eq!(
            hard_supports(&outside),
            BTreeSet::from([edges(&graph, &[0, 1, 3, 4]), edges(&graph, &[0, 1, 2, 4])])
        );
        assert_eq!(
            hard_supports(&inside),
            BTreeSet::from([edges(&graph, &[2, 3]), edges(&graph, &[0, 1, 2, 4])])
        );
        let bubble = inside
            .ir
            .iter()
            .find(|term| term.hard == edges(&graph, &[2, 3]))
            .unwrap();
        assert_eq!(bubble.taylor_order, 2);
        assert_eq!(bubble.forest.len(), 1);
        assert_eq!(bubble.forest[0].steps.len(), 1);
        assert_eq!(bubble.preparation.len(), 2);
        Ok(())
    }

    #[test]
    fn uv_children_preserve_reference_mass_before_parent_taylor() -> Result<()> {
        test_initialise()?;
        let mut graph: Graph = finalized_runtime_dot!(
            digraph G {
                graph [projector=0];
                node [num=1];
                edge [particle="scalar_1" num=1];
                A -> A [id=0 lmb_id=0 source="{ufo_order:0}" sink="{ufo_order:1}"];
                A -> B [id=1 source="{ufo_order:2}" sink="{ufo_order:0}"];
                A -> B [id=2 lmb_id=1 source="{ufo_order:3}" sink="{ufo_order:1}"];
            }, "scalars"
        )?;
        for edge in [0, 1, 2].map(EdgeIndex::from) {
            graph.underlying[edge].particle =
                graph[edge].particle.clone().override_mass(Some(Atom::Zero));
        }
        let model = load_generic_model("scalars");
        let settings = UVgenerationSettings::default();
        let plan = RStarPlan::new(&graph, &model, EdgeIndex::from(0), &settings)?;
        let mass = symbolica::symbol!("rstar_plan_reference_mass"; Scalar);
        let children = plan.evaluate_uv_children(
            &graph,
            &model,
            &Atom::var(mass),
            crate::utils::vakint()?,
            &settings,
        )?;
        assert!(!children.contains_key(&plan.support));
        let tadpole = &children[&edges(&graph, &[0])];
        assert!(
            !tadpole.is_zero(),
            "the enclosing zero projector must not act"
        );
        assert!(tadpole.replace(mass).with(0).expand().is_zero());
        assert!(
            (tadpole - tadpole.replace(mass).with(1) * Atom::var(mass).pow(2))
                .expand()
                .is_zero(),
            "the complete UV polynomial retains its inherited mass squared"
        );
        Ok(())
    }
}
