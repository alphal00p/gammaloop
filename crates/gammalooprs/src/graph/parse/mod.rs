use feynkit_graph::DOD;
use std::{
    collections::{BTreeMap, BTreeSet},
    ops::Deref,
    path::{Path, PathBuf},
};

use crate::{
    cff::surface::SurfaceCache,
    graph::{
        FinalizedCut, FinalizedTopologyThresholdCandidate, GraphGroup, GroupId, LoopMomentumBasis,
        attribute_warnings::warn_about_unknown_attributes, edge::EdgeExtraData,
    },
    integrands::process::ParamBuilder,
    model::{Model, ParticleId, ParticleIdGammaLoopExt},
    momentum::sample::LoopIndex,
    numerator::{GlobalPrefactor, aind::Aind},
    processes::DotExportSettings,
    utils::GS,
    uv::UltravioletGraph,
};
use feynkit_graph::{
    DiagramEndpoint as FeynkitDiagramEndpoint, DiagramHalfEdge as FeynkitDiagramHalfEdge,
    FeynmanDiagram,
};
use idenso::{color::ColorSimplifySettings, tensor::SymbolicTensor};
use spenso::{
    shadowing::{IntoAtom, TensorCollectFilter, symbolica_utils::LogPrint},
    structure::{
        IndexedStructure,
        representation::{Minkowski, Representation},
        slot::{IsAbstractSlot, Slot},
    },
};

use color_eyre::{Result, Section};

use eyre::{Context, Ok, eyre};
use itertools::Itertools;
use linnet::{
    half_edge::{
        HedgeGraph, NodeIndex,
        involution::{EdgeData, EdgeIndex, Flow, Hedge, HedgePair},
        nodestore::NodeStorageVec,
        subgraph::{ModifySubSet, OrientedCut, SuBitGraph, SubSetOps},
        swap::Swap,
    },
    parser::{DotEdgeData, DotGraph, DotHedgeData, DotVertexData, GraphSet, HedgeParseError},
    permutation::Permutation,
};
use symbolica::atom::{Atom, AtomCore, AtomView};
use symbolica::id::ConditionResult;
use tracing::instrument;
use tracing::{debug, warn};
use typed_index_collections::TiVec;

use super::{
    Autogen, Edge, Graph, HedgeData, LMBext, Vertex,
    edge::{EdgeMass, ParseEdge},
    global::ParseData,
    hedge_data::{NumIndices, NumeratorIndexSource, ParseHedgeData},
    vertex::ParseVertex,
};

/// Extract oriented particles from hedges, filtering out dummy edges
pub fn extract_oriented_particles_from_vertex_hedges<I, V>(
    graph: &HedgeGraph<ParseEdge, V, ParseHedgeData>,
    hedges: I,
    model: &Model,
) -> Vec<ParticleId>
where
    I: Iterator<Item = Hedge>,
{
    hedges
        .filter_map(|h| {
            let eid = graph[&h];
            if graph[eid].is_dummy {
                return None;
            }
            let particle = graph[eid].particle.particle()?;
            Some(if graph.flow(h) != Flow::Sink {
                particle.antiparticle(model)
            } else {
                particle
            })
        })
        .collect()
}

// Type aliases for cleaner code
type NumGraph = HedgeGraph<ParseEdge, ParseVertex, HedgeData>;
type UnderlyingGraph = HedgeGraph<Edge, Vertex, HedgeData>;

pub(crate) use feynkit_graph::NumeratorInterface;

pub mod string_utils;
pub use string_utils::{FromStripedStr, StripParse, ToQuoted};

#[derive(Clone, Debug)]
pub struct ParseGraph {
    pub global_data: ParseData,
    pub graph: HedgeGraph<ParseEdge, ParseVertex, ParseHedgeData>,
}

impl Deref for ParseGraph {
    type Target = HedgeGraph<ParseEdge, ParseVertex, ParseHedgeData>;

    fn deref(&self) -> &Self::Target {
        &self.graph
    }
}

impl ParseGraph {
    pub fn n_anticommutating_loops(&self, model: &Model) -> usize {
        let anticommutating: SuBitGraph = self
            .graph
            .from_filter(|a| a.particle.is_anticommutating(model));

        self.graph.cyclotomatic_number(&anticommutating)
    }
    pub fn n_external_anticommutating_loops(&mut self, model: &Model) -> Result<usize> {
        let internal = self.n_anticommutating_loops(model);
        self.graph
            .sew(
                |_, ae, _, be| {
                    if let (Some(a), Some(b)) = (ae.data.is_cut, be.data.is_cut) {
                        a == b
                    } else {
                        false
                    }
                },
                |af, ae, bf, be| match (af, bf) {
                    (Flow::Sink, Flow::Source) => (Flow::Sink, ae),
                    (Flow::Source, Flow::Sink) => (Flow::Source, be),
                    _ => panic!("Cannot sew hedges with flow {:?} and {:?}", af, bf),
                },
            )
            .map_err(|e| eyre::eyre!("Graph sewing failed: {:?}", e))?;

        Ok(self.n_anticommutating_loops(model) - internal)
    }

    pub fn debug_dot(&self) -> String {
        DotGraph::from(self).debug_dot()
    }
    /// Return the explicit UFO slot recorded for every half-edge.
    ///
    /// A finalized runtime artifact must retain the rule-leg assignment made
    /// by FeynKit. Inferring it from particles here would reintroduce a second
    /// physics-generation path and is therefore rejected.
    pub(crate) fn hedge_order(&self) -> Result<Vec<u8>> {
        (0..self.n_hedges())
            .map(|index| {
                self.graph[Hedge(index)].ufo_order.ok_or_else(|| {
                    eyre!(
                        "finalized runtime DOT graph '{}' is missing ufo_order for half-edge {index}",
                        self.global_data.name
                    )
                })
            })
            .collect()
    }

    /// Attach runtime metadata to the finalized half-edge storage. The retained
    /// tensor fragments are transferred directly when runtime values are built.
    fn from_feynkit_diagram(diagram: &FeynmanDiagram) -> Result<Self> {
        let underlying = diagram.underlying();
        let graph = underlying.map_data_ref_result(
            |_, _, vertex| {
                Ok(ParseVertex {
                    name: Some(vertex.name.clone()),

                    vertex_rule: vertex.interaction,
                    num: None,
                    num_interface: None,
                    dod: None,
                })
            },
            |_, edge_id, pair, data| {
                let mut edge =
                    ParseEdge::new(data.data.particle).with_label(format!("e{}", edge_id.0));
                edge.is_dummy = data.data.is_dummy;
                edge.lmb_id = diagram
                    .loop_momentum_basis()
                    .loop_edges
                    .iter()
                    .position(|candidate| candidate.0 == edge_id.0)
                    .map(LoopIndex);
                if let (HedgePair::Paired { sink, .. }, Some(_)) = (pair, &data.data.external) {
                    edge.is_cut = Some(sink);
                    edge.initial_state_connection = true;
                }
                Ok(EdgeData::new(edge, data.orientation))
            },
            |(hedge, _)| {
                let edge = &underlying[underlying[&hedge]];
                let slot = match underlying.flow(hedge) {
                    Flow::Source => edge.source_slot(),
                    Flow::Sink => edge.target_slot(),
                };
                Ok(ParseHedgeData {
                    ufo_order: Some(
                        u8::try_from(slot.0)
                            .map_err(|_| eyre!("finalized UFO slot {} exceeds u8", slot.0))?,
                    ),
                })
            },
        )?;
        Ok(Self {
            global_data: ParseData {
                name: diagram.name().to_owned(),
                overall_factor: diagram.overall_factor().clone(),
                ..Default::default()
            },
            graph,
        })
    }

    pub(crate) fn from_parsed(graph: DotGraph, model: &Model) -> Result<Self> {
        warn_about_unknown_attributes(&graph);
        if graph
            .global_data
            .statements
            .contains_key("canonical_cuts_required")
        {
            return Err(eyre!(
                "cross-section runtime DOT does not carry canonical physical cuts or topology-threshold candidates; import the canonical FeynmanDiagram DOT artifact through FeynKit instead"
            ));
        }
        let global_data = graph.global_data.try_into()?;
        let graph = graph
            .graph
            .map_data_ref_result(
                |_, _, v| Ok(v),
                ParseEdge::parse(model),
                ParseHedgeData::parse(),
            )?
            .map_data_ref_result(
                ParseVertex::parse(model, &global_data),
                |_, _, _, e| Ok(e.map(Clone::clone)),
                |(_, h)| Ok(h.clone()),
            )?;

        Ok(Self { graph, global_data })
    }
}

/// Helper struct to hold initial data extracted from ParseGraph
struct InitialGraphData {
    overall_factor: Atom,
    global_prefactor: GlobalPrefactor,
    additional_params: Vec<Atom>,
    add_polarizations: bool,
    group_id: Option<GroupId>,
    is_group_master: bool,
    name: String,
    threshold_counterterms: Autogen<super::threshold_counterterms::ThresholdCountertermSpec>,
}

/// Result of processing cut edges
struct CutProcessingResult {
    lmb_ids: BTreeMap<LoopIndex, EdgeIndex>,
    xs_ext_id: BTreeMap<Hedge, (EdgeIndex, Hedge)>,
    initial_hedges: SuBitGraph,
    full_cut: SuBitGraph,
}

impl CutProcessingResult {
    fn permute(&mut self, graph: &mut NumGraph) -> Result<Permutation> {
        let (h_perm, edge_perm): (Vec<_>, Vec<_>) = self
            .xs_ext_id
            .iter()
            .enumerate()
            .map(|(target_pos, (_, (edge_idx, h_id)))| {
                ((h_id.0, target_pos), (edge_idx.0, target_pos))
            })
            .unzip();

        let per = Permutation::from_mappings(edge_perm, graph.n_edges()).unwrap();
        let perh = Permutation::from_mappings(h_perm, graph.n_hedges()).unwrap();

        debug!("Before: {}", graph.dot(&self.initial_hedges));
        <HedgeGraph<_, _, _> as Swap<Hedge>>::permute(graph, &perh);
        let trans = perh.transpositions();

        for (i, j) in trans.into_iter().rev() {
            self.full_cut.swap(Hedge(i), Hedge(j));
            // self.initial_hedges.swap(i, j);// initial hedges is already assuming permuted hedges
        }

        debug!("Before after: {}", graph.dot(&self.initial_hedges));
        <HedgeGraph<_, _, _> as Swap<EdgeIndex>>::permute(graph, &per);
        for edge in self.lmb_ids.values_mut() {
            *edge = EdgeIndex(per[edge.0]);
        }

        debug!(" after: {}", graph.dot(&self.initial_hedges));
        Ok(per)
    }
}

fn display_graph_source_path(path: &Path) -> PathBuf {
    if path.is_absolute() {
        return path.canonicalize().unwrap_or_else(|_| path.to_path_buf());
    }

    let joined = std::env::current_dir()
        .map(|cwd| cwd.join(path))
        .unwrap_or_else(|_| path.to_path_buf());
    joined.canonicalize().unwrap_or(joined)
}

impl Graph {
    pub fn dot_serialize(&self, settings: &DotExportSettings) -> String {
        let mut out = String::new();
        self.dot_serialize_fmt(&mut out, settings).unwrap();
        out
    }

    pub(crate) fn dot_serialize_io(
        &self,
        writer: &mut impl std::io::Write,
        settings: &DotExportSettings,
    ) -> Result<(), std::io::Error> {
        let g = self
            .to_dot_graph_with_settings(settings)
            .map_err(std::io::Error::other)?;
        g.write_io(writer)
    }

    #[allow(dead_code)]
    pub(crate) fn dot_split_serialize_io(
        &self,
        writer: &mut impl std::io::Write,
    ) -> Result<(), std::io::Error> {
        let g = self.to_split_dotgraph();
        g.write_io(writer)
    }

    pub fn dot_serialize_fmt(
        &self,
        writer: &mut impl std::fmt::Write,
        settings: &DotExportSettings,
    ) -> Result<(), std::fmt::Error> {
        let g = self
            .to_dot_graph_with_settings(settings)
            .map_err(|_| std::fmt::Error)?;
        g.write_fmt(writer)
    }

    pub(crate) fn from_parsed_with_validation(graph: ParseGraph, model: &Model) -> Result<Self> {
        let res = Self::from_parsed(graph, model)?;
        res.validate_full_numerator_tensor_network()
            .with_context(|| {
                format!(
                    "Failed to validate full numerator tensor network for graph {}",
                    res.name
                )
            })?;
        Ok(res)
    }

    /// Enrich a finalized FeynKit diagram with GammaLoop runtime caches.
    ///
    /// This is deliberately a mechanical downward conversion: it copies the
    /// canonical topology, symbolic fragments, factors, projector, and routing.
    /// It never canonicalizes again, resolves a vertex rule, generates a
    /// numerator, or chooses a second loop-momentum basis.
    /// Transfer finalized tensor fragments and global factors directly. The
    /// metadata map preserves the native frame without exporting or re-admitting
    /// its already checked numerator values.
    pub(crate) fn from_feynkit(
        diagram: &FeynmanDiagram,
        group_id: Option<GroupId>,
        is_group_master: bool,
    ) -> Result<Self> {
        diagram.validate().with_context(|| {
            format!(
                "FeynKit diagram '{}' is not finalized consistently",
                diagram.name()
            )
        })?;

        let model = diagram.model();
        let mut parsed = ParseGraph::from_feynkit_diagram(diagram)?;
        parsed.global_data.group_id = group_id;
        parsed.global_data.is_group_master = is_group_master;
        let (mut initial_data, graph) =
            Self::extract_initial_data(&parsed, model, NumeratorIndexSource::Feynkit)?;
        initial_data.global_prefactor = GlobalPrefactor {
            num: diagram.numerator_prefactor().clone().try_into_indexed()?,
            projector: diagram.projector().clone().try_into_indexed()?,
        };
        initial_data.add_polarizations = false;

        let filter = |half_edges: &[FeynkitDiagramHalfEdge]| -> Result<SuBitGraph> {
            let mut filter: SuBitGraph = graph.empty_subgraph();
            for half_edge in half_edges {
                let pair = diagram.underlying()[&EdgeIndex(half_edge.edge.0)].1;
                let flow = match half_edge.endpoint {
                    FeynkitDiagramEndpoint::Source => Flow::Source,
                    FeynkitDiagramEndpoint::Target => Flow::Sink,
                };
                let hedge = match pair {
                    HedgePair::Paired { source, sink } | HedgePair::Split { source, sink, .. } => {
                        if flow == Flow::Source {
                            source
                        } else {
                            sink
                        }
                    }
                    HedgePair::Unpaired {
                        hedge,
                        flow: attached,
                    } if attached == flow => hedge,
                    _ => {
                        return Err(eyre!(
                            "finalized selection refers to absent endpoint {half_edge:?}"
                        ));
                    }
                };
                filter.add(hedge);
            }
            Ok(filter)
        };
        let mut initial_hedges: SuBitGraph = graph.empty_subgraph();
        for (pair, _, edge) in diagram.underlying().iter_edges() {
            if edge.data.external.is_some()
                && let HedgePair::Paired { sink, .. } = pair
            {
                initial_hedges.add(sink);
            }
        }
        let initial_state_cut = OrientedCut::from_underlying_strict(initial_hedges, &graph)?;
        let finalized_cuts = diagram
            .cuts()
            .iter()
            .map(|cut| {
                Ok(FinalizedCut {
                    cut: OrientedCut::from_underlying_strict(filter(&cut.cut)?, &graph)?,
                    left: filter(&cut.left.half_edges)?,
                    right: filter(&cut.right.half_edges)?,
                })
            })
            .collect::<Result<Vec<_>>>()?;
        let finalized_topology_threshold_candidates = diagram
            .topology_threshold_candidates()
            .iter()
            .map(|candidate| {
                Ok(FinalizedTopologyThresholdCandidate {
                    cut: OrientedCut::from_underlying_strict(filter(&candidate.cut)?, &graph)?,
                    left: filter(&candidate.left)?,
                    right: filter(&candidate.right)?,
                })
            })
            .collect::<Result<Vec<_>>>()?;
        let loop_momentum_basis: LoopMomentumBasis = diagram
            .loop_momentum_basis()
            .to_routing(diagram.underlying())
            .into();
        let global_prefactor = initial_data.global_prefactor;
        let polarizations = global_prefactor.polarizations();
        let param_builder = ParamBuilder::new(
            &(&polarizations, &graph),
            model,
            &loop_momentum_basis,
            initial_data.additional_params.clone(),
        )?;

        let underlying = Self::build_underlying_graph(
            graph,
            model,
            &param_builder,
            &initial_data.name,
            |vertex, _| {
                Ok(diagram.underlying()[vertex]
                    .numerator
                    .clone()
                    .try_into_indexed()?)
            },
            |edge, _| {
                Ok(diagram.underlying()[edge]
                    .numerator
                    .clone()
                    .try_into_indexed()?)
            },
        )
        .with_context(|| {
            format!(
                "failed to build GammaLoop runtime storage for finalized FeynKit diagram {}",
                initial_data.name
            )
        })?;

        let mut full_without_initials = underlying.full_filter();
        full_without_initials.subtract_with(&initial_state_cut.left);
        let mut tree_edges = underlying.bridges_of(&full_without_initials);
        tree_edges.union_with(&initial_state_cut.left);
        let mut result = Graph {
            model: diagram.model_arc(),
            overall_factor: initial_data.overall_factor,
            polarizations,
            global_prefactor,
            tree_edges,
            name: initial_data.name,
            loop_momentum_basis,
            initial_state_cut,
            underlying,
            surface_cache: SurfaceCache::default(),
            group_id: initial_data.group_id,
            is_group_master: initial_data.is_group_master,
            param_builder,
            finalized_cuts,
            finalized_topology_threshold_candidates,
            threshold_counterterms: initial_data.threshold_counterterms,
        };
        result.param_builder = ParamBuilder::new(
            &result,
            model,
            &result.loop_momentum_basis,
            initial_data.additional_params,
        )?;
        let runtime = result.numerator(&result.full_filter(), &result.empty_subgraph())?;
        let runtime = runtime.state.expr;
        let source = diagram.numerator().clone().try_into_indexed()?;
        if !runtime.same_value(&source) {
            // Keep graph numerators factorized, including at this validation
            // boundary. Only an exact coefficient proof certifies a rewrite.
            if runtime.coefficients_equal(&source, TensorCollectFilter::<0>::Tensors)?
                != ConditionResult::True
            {
                return Err(eyre!(
                    "GammaLoop runtime conversion of FeynKit diagram '{}' could not certify its finalized numerator",
                    diagram.name()
                ));
            }
        }
        result.validate_spatial_momentum_convention()?;
        Ok(result)
    }

    #[instrument(skip_all, fields(graph= %graph.debug_dot(),name = %graph.global_data.name.as_str()))]
    pub(crate) fn from_parsed(graph: ParseGraph, model: &Model) -> Result<Self> {
        if graph.global_data.projectors.is_none() {
            return Err(eyre!(
                "finalized DOT graph '{}' must provide an explicit projector (use `1` when no projector is required)",
                graph.global_data.name
            ));
        }
        for (vertex, _, data) in graph.graph.iter_nodes() {
            if data.num.is_none() {
                return Err(eyre!(
                    "finalized DOT graph '{}' is missing the numerator for vertex {vertex}",
                    graph.global_data.name
                ));
            }
        }
        for (_, edge, data) in graph.graph.iter_edges() {
            if data.data.num.is_none() {
                return Err(eyre!(
                    "finalized DOT graph '{}' is missing the numerator for edge {edge}",
                    graph.global_data.name
                ));
            }
            if data.data.is_cut.is_some() && !data.data.initial_state_connection {
                return Err(eyre!(
                    "finalized runtime DOT graph '{}' contains cross-section sewing metadata but no canonical physical cuts; import the canonical FeynmanDiagram DOT artifact through FeynKit instead",
                    graph.global_data.name
                ));
            }
        }

        let (mut initial_data, mut graph) =
            Self::extract_initial_data(&graph, model, NumeratorIndexSource::Literal)?;

        // Sew the graph based on cut edges
        graph
            .sew(
                |_, ae, _, be| {
                    if let (Some(a), Some(b)) = (ae.data.is_cut, be.data.is_cut) {
                        a == b
                    } else {
                        false
                    }
                },
                |af, ae, bf, be| match (af, bf) {
                    (Flow::Sink, Flow::Source) => (Flow::Sink, ae),
                    (Flow::Source, Flow::Sink) => (Flow::Source, be),
                    _ => panic!("Cannot sew hedges with flow {:?} and {:?}", af, bf),
                },
            )
            .map_err(|e| eyre::eyre!("Graph sewing failed: {:?}", e))?;

        let mut cut_result = Self::process_cut_edges(&graph)?;

        let edge_permutation = cut_result.permute(&mut graph)?;
        initial_data
            .threshold_counterterms
            .remap_edges(edge_permutation.map())?;

        let initial_state_cut =
            OrientedCut::from_underlying_strict(cut_result.initial_hedges, &graph)?;

        debug!("Initial state cut: {}", graph.dot(&initial_state_cut.left));
        debug_assert!(!initial_data.add_polarizations);
        let global_prefactor = initial_data.global_prefactor;
        let polarizations = global_prefactor.polarizations();
        let loop_momentum_basis = Self::materialize_explicit_loop_momentum_basis(
            &graph,
            &cut_result.full_cut,
            &cut_result.lmb_ids,
            &cut_result.xs_ext_id,
        )
        .with_context(|| format!("Failed to build lmb for graph {}", initial_data.name))?;
        let param_builder = ParamBuilder::new(
            &(&polarizations, &graph),
            model,
            &loop_momentum_basis,
            initial_data.additional_params.clone(),
        )?;

        let underlying = Self::build_underlying_graph(
            graph,
            model,
            &param_builder,
            &initial_data.name,
            |vertex, data| {
                let expression = data.num.clone().ok_or_else(|| {
                    eyre!("finalized graph is missing the numerator for vertex {vertex}")
                })?;
                Ok(
                    NumeratorInterface::admit(expression, data.num_interface.clone())?
                        .try_into_indexed()?,
                )
            },
            |edge, data| {
                let expression = data.num.clone().ok_or_else(|| {
                    eyre!("finalized graph is missing the numerator for edge {edge}")
                })?;
                Ok(
                    NumeratorInterface::admit(expression, data.num_interface.clone())?
                        .try_into_indexed()?,
                )
            },
        )
        .with_context(|| format!("Failed to build underlying graph {}", initial_data.name))?;

        let mut full_without_initials = underlying.full_filter();
        full_without_initials.subtract_with(&initial_state_cut.left);
        let mut tree_edges = underlying.bridges_of(&full_without_initials);
        tree_edges.union_with(&initial_state_cut.left);

        let mut g = Graph {
            model: std::sync::Arc::new(model.clone()),
            overall_factor: initial_data.overall_factor,
            polarizations: global_prefactor.polarizations(),
            global_prefactor,
            tree_edges,
            name: initial_data.name,
            loop_momentum_basis,
            initial_state_cut,
            underlying,
            surface_cache: SurfaceCache::default(),
            group_id: initial_data.group_id,
            is_group_master: initial_data.is_group_master,
            param_builder,
            finalized_cuts: Vec::new(),
            finalized_topology_threshold_candidates: Vec::new(),
            threshold_counterterms: initial_data.threshold_counterterms,
        };

        g.threshold_counterterms.validate_for_graph(&g)?;

        let external_momentum_edge_order = g.external_momentum_edge_order();
        g.loop_momentum_basis
            .canonicalize_external_order(&external_momentum_edge_order);

        let updated_param_builder_with_lmb = ParamBuilder::new(
            &g,
            model,
            &g.loop_momentum_basis,
            initial_data.additional_params,
        )?;

        debug!(
            "Updated param builder with LMB: {}\n{}",
            g.loop_momentum_basis,
            updated_param_builder_with_lmb.table(),
        );

        g.param_builder = updated_param_builder_with_lmb;

        debug!("{}", g.debug_dot());

        g.validate_spatial_momentum_convention().with_context(|| {
            format!(
                "Failed to validate the spatial-momentum convention for graph {}",
                g.name
            )
        })?;

        Ok(g)
    }

    fn validate_spatial_momentum_convention(&self) -> Result<()> {
        let n_edges = self.underlying.n_edges();
        let full_num = self
            .numerator(&self.full_filter(), &self.empty_subgraph())
            .get_single_atom()
            .unwrap()
            * &self.global_prefactor.num
            * &self.global_prefactor.projector
            * &self.overall_factor;
        let mut invalid_q3 = None;
        full_num.visitor(&mut |view| {
            if invalid_q3.is_none()
                && let AtomView::Fun(q3) = view
                && q3.get_symbol() == GS.emr_vec
            {
                let arguments = q3.iter().collect::<Vec<_>>();
                let invalid_reason = match arguments.as_slice() {
                    [edge, structure] => {
                        let edge_id = i64::try_from(*edge)
                            .ok()
                            .and_then(|edge| usize::try_from(edge).ok());
                        let is_represented = matches!(structure, AtomView::Fun(rep)
                            if rep.get_nargs() == 1)
                            && Representation::<Minkowski>::try_from(*structure)
                                .ok()
                                .and_then(|rep| usize::try_from(rep).ok())
                                == Some(4);
                        let is_indexed = matches!(structure, AtomView::Fun(rep)
                            if rep.get_nargs() == 2)
                            && Slot::<Minkowski, Aind>::try_from(*structure)
                                .ok()
                                .and_then(|slot| usize::try_from(slot.rep()).ok())
                                == Some(4);

                        if edge_id.is_none() {
                            Some(format!(
                                "found non-canonical edge argument {}",
                                edge.log_print(None)
                            ))
                        } else if edge_id.is_some_and(|edge| edge >= n_edges) {
                            Some(format!(
                                "found edge ID {} outside the valid range 0..{}",
                                edge.log_print(None),
                                n_edges,
                            ))
                        } else {
                            (!is_represented && !is_indexed).then(|| {
                                format!(
                                    "found non-canonical second argument {}",
                                    structure.log_print(None)
                                )
                            })
                        }
                    }
                    _ => Some(format!("found arity {}", arguments.len())),
                };
                if let Some(reason) = invalid_reason {
                    invalid_q3 = Some((reason, view.into_atom()));
                }
            }
            true
        });
        if let Some((reason, q3)) = invalid_q3 {
            return Err(eyre!(
                "Spatial momentum Q3 in graph '{}' must have exactly two arguments, beginning with a non-negative edge ID: Q3(edge,mink(4)) for a represented vector or Q3(edge,mink(4,index)) for an indexed tensor-network image; explicit components are stored as Q(edge,cind(i)); {reason} in {}",
                self.name,
                q3.log_print(None),
            ));
        }

        let mut invalid_q = None;
        full_num.visitor(&mut |view| {
            if invalid_q.is_none()
                && let AtomView::Fun(q) = view
                && q.get_symbol() == GS.emr_mom
            {
                let Some(edge) = q.iter().next() else {
                    invalid_q = Some(("found no edge argument".to_owned(), view.into_atom()));
                    return true;
                };
                let edge_id = i64::try_from(edge)
                    .ok()
                    .and_then(|edge| usize::try_from(edge).ok());
                let invalid_reason = if edge_id.is_none() {
                    Some(format!(
                        "found non-canonical edge argument {}",
                        edge.log_print(None)
                    ))
                } else {
                    edge_id.filter(|edge| *edge >= n_edges).map(|_| {
                        format!(
                            "found edge ID {} outside the valid range 0..{}",
                            edge.log_print(None),
                            n_edges,
                        )
                    })
                };
                if let Some(reason) = invalid_reason {
                    invalid_q = Some((reason, view.into_atom()));
                }
            }
            true
        });
        if let Some((reason, q)) = invalid_q {
            return Err(eyre!(
                "Momentum Q in graph '{}' must begin with a valid non-negative graph-edge ID; {reason} in {}",
                self.name,
                q.log_print(None),
            ));
        }

        Ok(())
    }

    pub(crate) fn validate_full_numerator_tensor_network(&self) -> Result<()> {
        self.validate_spatial_momentum_convention()?;
        let full_num = self
            .numerator(&self.full_filter(), &self.empty_subgraph())
            .get_single_atom()
            .unwrap()
            * &self.global_prefactor.num
            * &self.global_prefactor.projector
            * &self.overall_factor;
        let color_simplified = full_num
            .as_view()
            .simplify_color_with(ColorSimplifySettings::default().with_cof_dimension_invariants());
        if !full_num.is_zero() && color_simplified.is_zero() {
            warn!(
                "Full numerator for graph '{}' becomes zero after color algebra. The graph/projector color structure likely annihilates the amplitude.",
                self.name
            );
        }
        let net = full_num
            .parse_to_symbolic_net::<Aind>(&ParseSettings::default())
            .map_err(Report::from)?;
        let dangling = net.graph.dangling_indices();
        if !dangling.is_empty() {
            return Err(eyre!(
                "Full numerator still has dangling tensor indices: \n{}",
                dangling
                    .iter()
                    .map(|slot| format!(
                        "{}:{}",
                        slot.to_atom().log_print(None),
                        slot.to_atom().to_plain_string()
                    ))
                    .join(",\n")
            ));
        }

        Ok(())
    }

    fn extract_initial_data(
        parse_graph: &ParseGraph,
        model: &Model,
    ) -> Result<(InitialGraphData, NumGraph)> {
        let hedge_order = parse_graph.hedge_order(model)?;
        let global_data = &parse_graph.global_data;

        let initial_data = InitialGraphData {
            additional_params: global_data.parameters.clone(),
            overall_factor: global_data.overall_factor.clone(),
            global_prefactor: GlobalPrefactor {
                num: global_data.num.clone(),
                projector: global_data.projectors.clone().unwrap_or(Atom::one()),
            },
            add_polarizations: global_data.projectors.is_none(),
            group_id: global_data.group_id,
            is_group_master: global_data.is_group_master,
            name: global_data.name.clone(),
            threshold_counterterms: Autogen::from_option_or_generate(
                global_data.threshold_counterterms.clone(),
                Default::default,
            ),
        };

        let num_graph = parse_graph.graph.map_data_ref(
            |_, _, v| v.clone(),
            |_, _, _, e| e.map(|e| e.clone()),
            |h, hd| HedgeData {
                num_indices: NumIndices::parse(parse_graph)(h, hd),
                ufo_order: Autogen::from_option_or_generate(hd.ufo_order, || hedge_order[h.0]),
            },
        );

        Ok((initial_data, num_graph))
    }

    fn process_cut_edges(graph: &NumGraph) -> Result<CutProcessingResult> {
        let mut lmb_ids: BTreeMap<LoopIndex, EdgeIndex> = BTreeMap::new();
        let mut xs_ext_id: BTreeMap<Hedge, (EdgeIndex, Hedge)> = BTreeMap::new();
        let mut full_cut: SuBitGraph = graph.full_filter();

        for (p, eid, e) in graph.iter_edges() {
            let HedgePair::Paired { sink, .. } = p else {
                if e.data.is_cut.is_some() {
                    //As we have already sewn the graph, all cut edges must be paired, failure to do so would indicate a bug
                    return Err(eyre!("Cut edge must be paired"));
                } else {
                    continue;
                }
            };

            if let Some(lmb_id) = e.data.lmb_id {
                if let Some(old_value) = lmb_ids.insert(lmb_id, eid) {
                    return Err(eyre!(
                        "lmb_id {lmb_id:?} already exists with value {old_value:?}",
                    ));
                }
                debug!("Cutting {eid} for lmb_id{lmb_id}");
                full_cut.sub(p);
            } else if let Some(h) = e.data.is_cut {
                if let Some(old_value) = xs_ext_id.insert(h, (eid, sink)) {
                    return Err(eyre!("h {h:?} already exists with value {old_value:?}",));
                }
                full_cut.sub(p);
            }
        }

        // debug!("Graph now:{}", graph.dot(full_cut));

        let mut initial_hedges: SuBitGraph = graph.empty_subgraph();
        for (target_pos, _) in xs_ext_id.iter().enumerate() {
            initial_hedges.add(Hedge(target_pos));
        }

        Ok(CutProcessingResult {
            full_cut,
            lmb_ids,
            xs_ext_id,
            initial_hedges,
        })
    }

    fn generate_numerators(graph: &NumGraph, model: &Model) -> Result<NumeratorData> {
        let mut color_edge: EdgeVec<_> = vec![Atom::num(1); graph.n_edges()].into();
        let mut spin_edge: EdgeVec<_> = vec![Atom::i(); graph.n_edges()].into();

        for (p, eid, e) in graph.iter_edges() {
            if let HedgePair::Paired { source, sink } = p {
                let prop = e
                    .data
                    .particle
                    .particle()
                    .map(|p| model.get_propagator_for_particle(&p.name).numerator.clone())
                    .unwrap_or(Atom::num(1));

                color_edge[eid] = graph[source].color_kronekers(&graph[sink]);

                let spin_slots = [
                    &graph[source].num_indices.spin_indices.edge_indices,
                    &graph[sink].num_indices.spin_indices.edge_indices,
                ];

                let momenta = [(Flow::Source, graph[&source]), (Flow::Sink, graph[&sink])];
                let spin_nume = UFO.reindex_spin(&spin_slots, &momenta, prop, |i| {
                    Aind::Edge(usize::from(eid) as u16, i as u16)
                })?;

                spin_edge[eid] = spin_nume;
            }
        }

        let (color_vertex, spin_vertex) = Self::generate_vertex_numerators(graph)?;

        Ok(NumeratorData {
            color_edge,
            spin_edge,
            color_vertex,
            spin_vertex,
        })
    }

    fn setup_global_prefactor_and_params<'a, A: Into<AtomOrView<'a>>, P: IntoIterator<Item = A>>(
        mut global_prefactor: GlobalPrefactor,
        add_polarizations: bool,
        params: P,
        initial_state_cut: &OrientedCut,
        graph: &NumGraph,
        model: &Model,
    ) -> Result<(GlobalPrefactor, ParamBuilder)> {
        debug!("Initial state cut: {}", graph.dot(&initial_state_cut.left));

        if add_polarizations {
            let external_edges = initial_state_cut
                .left
                .union(initial_state_cut)
                .union(&graph.external_filter::<SuBitGraph>());
            let polarizations = graph.generate_polarizations_of(&external_edges);
            global_prefactor.projector *= polarizations;
        }

        let polarizations = global_prefactor.polarizations();
        let param_builder =
            ParamBuilder::new(&(&polarizations, graph), model, &graph.lmb(), params)?;

        Ok((global_prefactor, param_builder))
    }

    #[allow(clippy::type_complexity)]
    fn generate_vertex_numerators(
        graph: &NumGraph,
    ) -> Result<(
        Vec<Option<ParamTensor<OrderedStructure<Euclidean, Aind>>>>,
        Vec<Option<ParamTensor<OrderedStructure<Euclidean, Aind>>>>,
    )> {
        let mut color_vertex: Vec<Option<ParamTensor<OrderedStructure<Euclidean, Aind>>>> =
            vec![None; graph.n_nodes()];
        let mut spin_vertex = color_vertex.clone();

        for (ni, c, v) in graph.iter_nodes() {
            let mut color_slots = vec![];
            let mut spin_slots = vec![];
            let mut order = vec![];
            let mut momenta = vec![];
            for h in c {
                color_slots.push(&graph[h].num_indices.color_indices.vertex_indices);
                spin_slots.push(&graph[h].num_indices.spin_indices.vertex_indices);
                order.push(graph[h].ufo_order.value);
                momenta.push((graph.flow(h), graph[&h]));
            }

            let perm = Permutation::sort(&order);
            perm.apply_slice_in_place(&mut color_slots);
            perm.apply_slice_in_place(&mut spin_slots);
            perm.apply_slice_in_place(&mut momenta);

            let Some(vertex_rule) = &v.vertex_rule else {
                continue;
            };

            let [mut color_structure, couplings, mut spin_structure] =
                vertex_rule.tensors(ni.aind(1), ni.aind(0));

            spin_structure.map_data_ref_mut_result(|a| {
                *a =
                    UFO.reindex_spin(&spin_slots, &momenta, (*a).clone(), |u| ni.aind(u as u16))?;

                Ok(())
            })?;

            // couplings.map_data_mut(|a| *a = UFO.normalize_complex((*a).clone()));

            color_structure.map_data_ref_mut_result(|a| {
                *a = UFO.reindex_color(&color_slots, (*a).clone(), |u| ni.aind(u as u16))?;
                Ok(())
            })?;

            spin_vertex[ni.0] = Some(spin_structure.contract(&couplings).unwrap());
            color_vertex[ni.0] = Some(color_structure);
        }

        Ok((color_vertex, spin_vertex))
    }

    fn build_underlying_graph(
        graph: NumGraph,
        initial_state_cut: &OrientedCut,
        numerators: &NumeratorData,
        model: &Model,
        param_builder: &ParamBuilder,
    ) -> Result<UnderlyingGraph> {
        let mut loop_edge_filter = graph.full_filter();
        loop_edge_filter.subtract_with(&graph.bridges_of(&loop_edge_filter));
        let mut vertex_color_nums = numerators.color_vertex.clone();
        let mut vertex_spin_nums = numerators.spin_vertex.clone();
        let intermediate: UnderlyingGraph = graph.map_result(
            |_, i, v| {
                let num = match v.num {
                    Some(num) => Autogen::explicit(num),
                    None => Autogen::generated(
                        vertex_spin_nums[i.0]
                            .take()
                            .unwrap()
                            .contract(&vertex_color_nums[i.0].take().unwrap())
                            .unwrap()
                            .scalar()
                            .unwrap(),
                    ),
                };

                let dod = match v.dod {
                    Some(dod) => Autogen::explicit(dod),
                    None => Autogen::generated(if num.autogenerated {
                        v.vertex_rule.as_ref().map(|vr| vr.dod).unwrap_or(0)
                    } else {
                        num.all_dod()
                    }),
                };

                Ok(Vertex {
                    name: Autogen::from_option_or_generate(v.name, || i.to_string()),
                    num,
                    dod,
                    vertex_rule: v.vertex_rule,
                })
            },
            |_, _, p, eid, ed| {
                let e = ed.data;
                if e.particle.is_fermion() && !e.particle.is_self_antiparticle()&&  e.particle.orientation() != ed.orientation {
                    return Err(eyre!(
                        "Edge orientation {:?} does not match particle orientation {:?} for edge {},{}",
                        ed.orientation,
                        e.particle.orientation(),
                        eid,
                        e
                    ));
                    }

                let mass = EdgeMass::from_atom(e.mass_atom(), model, param_builder)?;

                let num = match e.num {
                    Some(num) => Autogen::explicit(num),
                    None => Autogen::generated(if initial_state_cut.left.intersects(&p) {
                        numerators.color_edge[eid].clone()
                    } else {
                        &numerators.color_edge[eid] * &numerators.spin_edge[eid]
                    }),
                };

                let dod = match e.dod {
                    Some(dod) => Autogen::explicit(dod),
                    None => Autogen::generated(if num.autogenerated {
                        if let Some(particle) = e.particle.particle() {
                            model.get_propagator_for_particle(&particle.name).dod
                        } else {
                            -2
                        }
                    } else {
                        num.edge_dod(eid) -2
                    }),
                };

                Ok(EdgeData::new(
                    Edge {
                        mass,
                        is_dummy: e.is_dummy,
                        name: Autogen::from_option_or_generate(e.name, || eid.to_string()),
                        particle: e.particle,
                        num,
                        dod,
                        extra_data: EdgeExtraData {
                            momtrop_edge_power: e.momtrop_edge_power,
                            vakint_edge_power: e.vakint_edge_power,
                        }
                    },
                    ed.orientation,
                ))
            },
            |_, h| Ok(h),
        )?;

        Ok(intermediate)
    }

    /// Materialize signatures from the exact loop edges selected by `lmb_id`.
    ///
    /// The spanning-forest routine only propagates momenta through that fixed
    /// complement. Missing, extra, or substituted loop edges are rejected, so
    /// the DOT runtime import cannot silently choose a second basis.
    fn materialize_explicit_loop_momentum_basis(
        graph: &NumGraph,
        full_cut: &SuBitGraph,
        lmb_ids: &BTreeMap<LoopIndex, EdgeIndex>,
        xs_ext_id: &BTreeMap<Hedge, (EdgeIndex, Hedge)>,
    ) -> Result<LoopMomentumBasis> {
        debug!("{}", graph.dot(full_cut));

        let mut full = graph.full_filter();
        for (pair, _, edge) in graph.iter_edges() {
            if edge.data.is_dummy {
                full.sub(pair);
            }
        }
        let total_loops = graph.cyclotomatic_number(&full);
        let explicit_loops = total_loops.checked_sub(xs_ext_id.len()).ok_or_else(|| {
            eyre!(
                "graph has {total_loops} loops but {} cut edges were marked as external",
                xs_ext_id.len()
            )
        })?;
        let expected_ids = (0..explicit_loops).map(LoopIndex).collect_vec();
        let actual_ids = lmb_ids.keys().copied().collect_vec();
        if actual_ids != expected_ids {
            return Err(eyre!(
                "finalized DOT graph must label exactly {explicit_loops} loop edges with contiguous lmb_id values 0..{explicit_loops}; found {actual_ids:?}"
            ));
        }

        let external = graph.internal_crown(&full);
        // Dummy edges cannot carry physical momentum or enter the spanning
        // forest used to materialize the explicitly selected loop coordinates.
        let forest_guide = full_cut.intersection(&full);
        let mut loop_momentum_basis = graph.lmb_impl(&full, &forest_guide, external)?;

        for e in 0..xs_ext_id.len() {
            let (l, _) = loop_momentum_basis
                .loop_edges
                .iter()
                .find_position(|a| *a == &EdgeIndex(e))
                .ok_or_else(|| {
                    eyre!("cut edge {e} is not a loop edge in the explicit momentum basis")
                })?;

            loop_momentum_basis.put_loop_to_ext(LoopIndex(l));
        }

        let materialized_edges = loop_momentum_basis
            .loop_edges
            .iter()
            .copied()
            .collect::<BTreeSet<_>>();
        let selected_edges = lmb_ids.values().copied().collect::<BTreeSet<_>>();
        if materialized_edges != selected_edges {
            return Err(eyre!(
                "lmb_id edges {selected_edges:?} do not form a loop-momentum basis; materialized edges were {materialized_edges:?}"
            ));
        }

        // Put the explicitly labelled edges in their requested order.
        for (target, edge) in lmb_ids {
            let current = loop_momentum_basis
                .loop_edges
                .iter()
                .position(|candidate| candidate == edge)
                .ok_or_else(|| {
                    eyre!(
                        "explicit loop edge {edge} disappeared while ordering the finalized basis"
                    )
                })?;
            if current != target.0 {
                loop_momentum_basis.swap_loops(LoopIndex(current), *target);
            }
        }

        Ok(loop_momentum_basis)
    }

    /// Import a fully finalized amplitude runtime artifact.
    ///
    /// The artifact must contain explicit numerator fragments, projector, UFO
    /// half-edge slots, and loop-momentum-basis IDs. Cross-section DOT must be
    /// imported as a canonical [`FeynmanDiagram`] so its typed cuts survive.
    pub fn from_finalized_runtime_dot(graph: DotGraph, model: &Model) -> Result<Self> {
        Self::from_parsed(ParseGraph::from_parsed(graph, model)?, model)
    }

    /// Import finalized amplitude runtime artifacts from one DOT file.
    pub fn from_finalized_runtime_file<P>(p: P, model: &Model) -> Result<Vec<Self>>
    where
        P: AsRef<Path>,
    {
        Self::from_finalized_runtime_path(p, model)
    }

    /// Import finalized amplitude runtime artifacts from a file or directory.
    pub fn from_finalized_runtime_path<P>(p: P, model: &Model) -> Result<Vec<Self>>
    where
        P: AsRef<Path>,
    {
        let path = p.as_ref();

        if path.is_dir() {
            // Load all .dot files from directory
            let mut all_graphs = Vec::new();
            let entries = std::fs::read_dir(path)
                .with_context(|| format!("Failed to read directory: {}", path.display()))?;

            let mut dot_files = Vec::new();
            for entry in entries {
                let entry = entry?;
                let file_path = entry.path();
                if file_path.is_file() && file_path.extension().is_some_and(|ext| ext == "dot") {
                    dot_files.push(file_path);
                }
            }

            // Sort files for consistent ordering
            dot_files.sort();

            for dot_file in dot_files {
                let graphs = Self::from_single_finalized_runtime_file(&dot_file, model)?;
                all_graphs.extend(graphs);
            }

            if all_graphs.is_empty() {
                return Err(eyre!(
                    "No .dot files found in directory: {}",
                    path.display()
                ));
            }

            Ok(all_graphs)
        } else {
            // Load single file
            Self::from_single_finalized_runtime_file(path, model)
        }
    }

    fn from_single_finalized_runtime_file<P>(p: P, model: &Model) -> Result<Vec<Self>>
    where
        P: AsRef<Path>,
    {
        let hedge_graph_set: GraphSet<
            DotEdgeData,
            DotVertexData,
            DotHedgeData,
            linnet::parser::GlobalData,
            NodeStorageVec<DotVertexData>,
        > = GraphSet::from_file(p.as_ref()).map_err(|a| match a {
            HedgeParseError::GraphFromFile(e) => match e.as_ref() {
                dot_parser::ast::GraphFromFileError::FileError(e) => eyre!(e.to_string())
                    .with_note(|| {
                        format!(
                            "Tried to access the file at: {}",
                            display_graph_source_path(p.as_ref()).display()
                        )
                    }),
                dot_parser::ast::GraphFromFileError::ParseError(e) => {
                    eyre!("Dot parsing error: {}", e)
                }
                dot_parser::ast::GraphFromFileError::PestParseError(e) => {
                    eyre!(e.to_string())
                }
            },
            HedgeParseError::ParseError(i) => color_eyre::Report::from(i),
            _ => {
                eyre!("Hedge parse error")
            }
        })?;
        Self::from_finalized_runtime_graph_set(hedge_graph_set, model)
    }

    /// Import finalized amplitude runtime artifacts from a DOT string.
    pub fn from_finalized_runtime_string<Str: AsRef<str>>(
        s: Str,
        model: &Model,
    ) -> Result<Vec<Self>> {
        let hedge_graph_set: GraphSet<
            DotEdgeData,
            DotVertexData,
            DotHedgeData,
            linnet::parser::GlobalData,
            NodeStorageVec<DotVertexData>,
        > = GraphSet::from_string(s).map_err(|a| match a {
            HedgeParseError::GraphFromFile(e) => match e.as_ref() {
                dot_parser::ast::GraphFromFileError::FileError(e) => {
                    eyre!(e.to_string())
                }
                dot_parser::ast::GraphFromFileError::ParseError(e) => {
                    eyre!("Dot parsing error: {}", e)
                }
                dot_parser::ast::GraphFromFileError::PestParseError(e) => {
                    eyre!(e.to_string())
                }
            },
            HedgeParseError::ParseError(i) => color_eyre::Report::from(i),
            _ => {
                eyre!("Hedge parse error")
            }
        })?;

        Self::from_finalized_runtime_graph_set(hedge_graph_set, model)
    }

    fn from_finalized_runtime_graph_set(
        set: GraphSet<
            DotEdgeData,
            DotVertexData,
            DotHedgeData,
            linnet::parser::GlobalData,
            NodeStorageVec<DotVertexData>,
        >,
        model: &Model,
    ) -> Result<Vec<Self>> {
        let mut graphs = Vec::new();

        for (graph, global_data) in set.set.into_iter().zip(set.global_data) {
            let graph = DotGraph { global_data, graph };
            debug!("Parsing: \n{}", graph.debug_dot());
            graphs.push(Graph::from_parsed(
                ParseGraph::from_parsed(graph, model)?,
                model,
            )?);
        }
        Ok(graphs)
    }
}

pub mod serialization;

/// completes and extract the user defined group structure on a lis of graphs
pub(crate) fn complete_group_parsing(graphs: &mut [Graph]) -> Result<TiVec<GroupId, GraphGroup>> {
    // validate the input
    let defined_group_ids = graphs
        .iter()
        .filter_map(|graph| graph.group_id)
        .sorted()
        .dedup()
        .collect_vec();

    let expected_group_ids = (0..defined_group_ids.len()).map(GroupId).collect_vec();

    if defined_group_ids != expected_group_ids {
        return Err(eyre!(
            "invalid group ids, group ids must start at 0 and contain no gaps"
        ));
    }
    // now set the remaining group ids
    let mut current_group_id = defined_group_ids.len();
    for graph in graphs.iter_mut() {
        if graph.group_id.is_none() {
            graph.group_id = Some(GroupId(current_group_id));
            graph.is_group_master = true;
            current_group_id += 1;
        }
    }

    let num_groups = current_group_id;

    // build the groups
    (0..num_groups)
        .map(|group_id| {
            let group_id = GroupId(group_id);
            let graphs_in_group = graphs
                .iter()
                .enumerate()
                .filter_map(|(i, g)| {
                    if g.group_id == Some(group_id) {
                        Some(i)
                    } else {
                        None
                    }
                })
                .collect_vec();

            // the special case of a single graph in the group is easy
            if graphs_in_group.len() == 1 {
                graphs[graphs_in_group[0]].is_group_master = true;
                Ok(GraphGroup {
                    master: graphs_in_group[0],
                    remaining: vec![],
                })
            } else {
                // see if a master is defined
                let master = graphs_in_group
                    .iter()
                    .find(|&&i| graphs[i].is_group_master)
                    .copied();

                if let Some(master) = master {
                    // find the remaining graphs and make sure no other master is defined
                    let remaining = graphs_in_group
                        .into_iter()
                        .filter(|&i| i != master)
                        .collect_vec();

                    let duplicate_master = remaining.iter().any(|&i| graphs[i].is_group_master);

                    if duplicate_master {
                        return Err(eyre!(
                            "Multiple group masters defined for group {group_id:?}"
                        ));
                    }
                    Ok(GraphGroup { master, remaining })
                } else {
                    // no master defined, take the first graph as master
                    let master = graphs_in_group[0];
                    graphs[master].is_group_master = true;
                    Ok(GraphGroup {
                        master,
                        remaining: graphs_in_group[1..].to_vec(),
                    })
                }
            }
        })
        .collect::<Result<TiVec<GroupId, GraphGroup>>>()
}

pub mod from_dot;
pub use from_dot::*;
#[cfg(test)]
pub mod tests;
