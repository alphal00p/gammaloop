use feynkit_graph::DOD;
use std::collections::HashSet;
use std::ops::Deref;

use crate::cff::surface::SurfaceCache;
use crate::graph::edge::EdgeExtraData;
use crate::graph::{
    FinalizedCut, FinalizedTopologyThresholdCandidate, GraphGroup, GroupId, LoopMomentumBasis,
};
use crate::integrands::process::ParamBuilder;
use crate::model::{Model, ParticleIdGammaLoopExt};
use crate::numerator::GlobalPrefactor;
use crate::numerator::aind::Aind;
use crate::processes::DotExportSettings;
use crate::utils::GS;
use crate::uv::UltravioletGraph;
use feynkit_graph::{
    DiagramEndpoint as FeynkitDiagramEndpoint, DiagramHalfEdge as FeynkitDiagramHalfEdge,
    FeynmanDiagram,
};
use idenso::color::ColorSimplifySettings;
use idenso::tensor::SymbolicTensor;
use idenso::{CookMode, CookSettings};
use spenso::shadowing::symbolica_utils::LogPrint;
use spenso::shadowing::{IntoAtom, TensorCollectFilter};
use spenso::structure::IndexedStructure;
use spenso::structure::abstract_index::{AIND_SYMBOLS, AbstractIndex};
use spenso::structure::partial::{PartialIndex, PartialStructure, PartialStructureExt};
use spenso::structure::representation::{Minkowski, Representation};
use spenso::structure::slot::{IsAbstractSlot, Slot};

use color_eyre::Result;

use eyre::{Context, Ok, eyre};
use itertools::Itertools;
use linnet::half_edge::involution::{EdgeData, EdgeIndex, Flow, Hedge, HedgePair};
use linnet::half_edge::subgraph::{ModifySubSet, OrientedCut, SuBitGraph, SubSetOps};
use linnet::half_edge::{EdgeAccessors, HedgeGraph, NodeIndex};
use symbolica::atom::{Atom, AtomCore, AtomView};
use symbolica::id::ConditionResult;
use tracing::warn;
use typed_index_collections::TiVec;

use super::edge::{EdgeMass, ParseEdge};
use super::global::ParseData;
use super::hedge_data::{NumIndices, ParseHedgeData};
use super::vertex::ParseVertex;
use super::{Autogen, Edge, Graph, HedgeData, Vertex};

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
    /// Runtime payload overrides must preserve the finalized frame before its
    /// half-edge IDs are used to construct cuts and momentum carriers.
    fn validate_feynkit_frame(&self, diagram: &FeynmanDiagram) -> Result<()> {
        let source = diagram.underlying();
        let graph = &self.graph;
        let same = source.n_live() == graph.n_live()
            && source.node_bound() == graph.node_bound()
            && source.n_edges() == graph.n_edges()
            && source.n_hedges() == graph.n_hedges()
            && source.iter_node_ids().eq(graph.iter_node_ids())
            && source.iter_edges().all(|(pair, edge, data)| {
                graph[&edge].1 == pair && graph.orientation(edge) == data.orientation
            })
            && (0..source.n_hedges()).all(|index| {
                let hedge = Hedge(index);
                graph.node_id(hedge) == source.node_id(hedge)
                    && graph.flow(hedge) == source.flow(hedge)
                    && graph[&hedge] == source[&hedge]
            });
        if !same {
            return Err(eyre!(
                "GammaLoop runtime conversion of FeynKit diagram '{}' changed its finalized topology",
                diagram.name().unwrap_or("<unnamed>")
            ));
        }
        Ok(())
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
                    dod: None,
                })
            },
            |_, edge_id, _, data| {
                let mut edge =
                    ParseEdge::new(data.data.particle).with_label(format!("e{}", edge_id.0));
                edge.is_dummy = data.data.is_dummy;
                Ok(EdgeData::new(edge, data.orientation))
            },
            |(hedge, _)| {
                let edge = &underlying[underlying[&hedge]];
                let slot = match underlying.flow(hedge) {
                    Flow::Source => edge.source_slot(),
                    Flow::Sink => edge.target_slot(),
                };
                Ok(ParseHedgeData {
                    ufo_order: u8::try_from(slot.0)
                        .map_err(|_| eyre!("finalized UFO slot {} exceeds u8", slot.0))?,
                })
            },
        )?;
        Ok(Self {
            global_data: ParseData {
                name: diagram.name().unwrap_or_default().to_owned(),
                overall_factor: diagram.overall_factor().clone(),
                ..Default::default()
            },
            graph,
        })
    }
}

/// Helper struct to hold initial data extracted from ParseGraph
struct InitialGraphData {
    overall_factor: Atom,
    additional_params: Vec<Atom>,
    group_id: Option<GroupId>,
    is_group_master: bool,
    name: String,
    threshold_counterterms: Autogen<super::threshold_counterterms::ThresholdCountertermSpec>,
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

    /// Enrich a finalized FeynKit diagram with GammaLoop runtime caches.
    ///
    /// This is deliberately a mechanical downward conversion: it copies the
    /// canonical topology, symbolic fragments, factors, projector, and routing.
    /// It never canonicalizes again, resolves a vertex rule, generates a
    /// numerator, or chooses a second loop-momentum basis.
    pub(crate) fn from_feynkit(
        diagram: &FeynmanDiagram,
        group_id: Option<GroupId>,
        is_group_master: bool,
    ) -> Result<Self> {
        let mut parsed = ParseGraph::from_feynkit_diagram(diagram)?;
        parsed.global_data.group_id = group_id;
        parsed.global_data.is_group_master = is_group_master;
        Self::from_feynkit_parse_graph(diagram, parsed)
    }

    /// Build the runtime graph of `diagram` from its converted payloads `parsed`,
    /// which may carry GammaLoop's own overrides of names, masses and degrees.
    /// Finalized tensor fragments and global factors transfer directly: the
    /// metadata map preserves the native frame without exporting or re-admitting
    /// their already checked values.
    fn from_feynkit_parse_graph(diagram: &FeynmanDiagram, parsed: ParseGraph) -> Result<Self> {
        diagram.validate().with_context(|| {
            format!(
                "FeynKit diagram '{}' is not finalized consistently",
                diagram.name().unwrap_or("<unnamed>")
            )
        })?;

        parsed.validate_feynkit_frame(diagram)?;
        let model = diagram.model();
        let (initial_data, graph) = Self::extract_initial_data(&parsed, model)?;

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
        let mut loop_momentum_basis: LoopMomentumBasis = diagram
            .loop_momentum_basis()
            .to_routing(diagram.underlying())
            .into();
        // FeynKit orders external momenta by edge. GammaLoop orders them by half-edge,
        // like the bases it enumerates, and keeps dummy legs among them.
        loop_momentum_basis.canonicalize_external_order(&Self::external_momentum_edge_order_of(
            &graph,
            &initial_state_cut,
        ));
        let global_prefactor = GlobalPrefactor {
            num: diagram.numerator_prefactor().clone().try_into_indexed()?,
            projector: diagram.projector().clone().try_into_indexed()?,
        };
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
            |vertex| {
                Ok(diagram.underlying()[vertex]
                    .numerator
                    .clone()
                    .try_into_indexed()?)
            },
            |edge| {
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
        result.threshold_counterterms.validate_for_graph(&result)?;
        result.certify_feynkit_numerators(diagram)?;
        result.validate_spatial_momentum_convention()?;
        Ok(result)
    }

    /// Every fragment transfers unchanged, but runtime subgraph products omit
    /// dangling edges. Native validation forces physical external fragments to
    /// one and permits non-unit dummy fragments, whose product can cancel or
    /// meet a zero retained numerator: only then compare the full numerators.
    fn certify_feynkit_numerators(&self, diagram: &FeynmanDiagram) -> Result<()> {
        let omits_a_fragment = diagram.underlying().iter_edges().any(|(pair, _, data)| {
            matches!(pair, HedgePair::Unpaired { .. })
                && !(data.data.numerator.is_scalar() && data.data.numerator.expression().is_one())
        });
        if !omits_a_fragment {
            return Ok(());
        }
        let runtime = self
            .numerator(&self.full_filter(), &self.empty_subgraph())?
            .state
            .expr;
        let source = diagram.numerator()?.clone().try_into_indexed()?;
        // Keep graph numerators factorized, including at this validation
        // boundary. Only an exact coefficient proof certifies a rewrite.
        if !runtime.same_value(&source)
            && runtime.coefficients_equal(&source, TensorCollectFilter::<0>::Tensors)?
                != ConditionResult::True
        {
            return Err(eyre!(
                "GammaLoop runtime conversion of FeynKit diagram '{}' omits a non-unit dangling numerator",
                diagram.name().unwrap_or("<unnamed>")
            ));
        }
        Ok(())
    }

    fn validate_spatial_momentum_convention(&self) -> Result<()> {
        let factor =
            feynkit_graph::expressions::evaluate_overall_factor(self.overall_factor.as_view());
        if factor.is_zero() {
            return Ok(());
        }
        let full_num = self
            .numerator(&self.full_filter(), &self.empty_subgraph())?
            .state
            .expr
            .multiply(self.global_prefactor.num.clone())?
            .multiply(self.global_prefactor.projector.clone())?;
        if full_num.is_zero() {
            return Ok(());
        }
        // Keep the established tensor carrier while validating the remaining
        // scalar source factor before its admission. Both passes borrow payloads;
        // Q3 diagnostics retain priority across all factors without a raw mirror.
        Self::validate_spatial_momentum_expressions(
            &[full_num.expression(), &factor],
            self.underlying.n_edges(),
            &self.name,
        )
    }

    fn validate_spatial_momentum_expressions(
        expressions: &[&Atom],
        n_edges: usize,
        name: &str,
    ) -> Result<()> {
        let mut invalid_q3 = None;
        for expression in expressions {
            expression.visitor(&mut |view| {
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
        }
        if let Some((reason, q3)) = invalid_q3 {
            return Err(eyre!(
                "Spatial momentum Q3 in graph '{}' must have exactly two arguments, beginning with a non-negative edge ID: Q3(edge,mink(4)) for a represented vector or Q3(edge,mink(4,index)) for an indexed tensor-network image; explicit components are stored as Q(edge,cind(i)); {reason} in {}",
                name,
                q3.log_print(None),
            ));
        }

        let mut invalid_q = None;
        for expression in expressions {
            expression.visitor(&mut |view| {
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
        }
        if let Some((reason, q)) = invalid_q {
            return Err(eyre!(
                "Momentum Q in graph '{}' must begin with a valid non-negative graph-edge ID; {reason} in {}",
                name,
                q.log_print(None),
            ));
        }

        Ok(())
    }

    pub(crate) fn validate_full_numerator_tensor_network(&self) -> Result<()> {
        self.validate_spatial_momentum_convention()?;
        let numerator = self
            .numerator(&self.full_filter(), &self.empty_subgraph())?
            .state
            .expr;
        let prefactor = self.global_numerator()?;
        let full_num = numerator.multiply(prefactor)?;
        let color_simplified = full_num
            .simplify_algebra(&idenso::tensor::AlgebraSettings {
                color: Some(ColorSimplifySettings::default().with_cof_dimension_invariants()),
                ..Default::default()
            })?
            .contract(idenso::tensor::ContractSettings {
                collect_chains: false,
                collect_traces: false,
                ..Default::default()
            })?;
        if !full_num.is_zero() && color_simplified.expression().is_zero() {
            warn!(
                "Full numerator for graph '{}' becomes zero after color algebra. The graph/projector color structure likely annihilates the amplitude.",
                self.name
            );
        }
        Self::require_contracted_indices(&full_num.into(), &HashSet::new())
    }

    /// Fail when a full numerator leaves a tensor port uncontracted, other than the
    /// explicitly indexed slots in `open`.
    pub(crate) fn require_contracted_indices(
        full_numerator: &SymbolicTensor<PartialStructure>,
        open: &HashSet<Atom>,
    ) -> Result<()> {
        // Tensor structures keep their indices encoded; compare their source spelling.
        let cooking = CookSettings::indices()
            .with_mode(CookMode::ReversibleEncoding)
            .with_representation_payloads(true, true);
        let dangling = full_numerator
            .structure()
            .logical_slots()
            .into_iter()
            .filter_map(|slot| {
                let explicit = match slot.aind {
                    PartialIndex::Explicit(index) => Some(
                        cooking.uncook(
                            slot.rep()
                                .slot::<AbstractIndex, _>(index)
                                .to_atom()
                                .as_view(),
                        ),
                    ),
                    PartialIndex::Open(_) => None,
                }
                // Anonymous ports stay unresolved; products only rename them apart.
                .filter(|slot| !slot.contains_symbol(AIND_SYMBOLS.openind));
                match explicit {
                    Some(slot) => (!open.contains(&slot))
                        .then(|| format!("{}:{}", slot.log_print(None), slot.to_plain_string())),
                    None => Some(format!(
                        "unresolved port of {}: GammaLoop requires an explicit index",
                        slot.rep()
                    )),
                }
            })
            .collect_vec();
        if dangling.is_empty() {
            return Ok(());
        }
        Err(eyre!(
            "Full numerator still has dangling tensor indices: \n{}",
            dangling.join(",\n")
        ))
    }

    fn extract_initial_data(
        parse_graph: &ParseGraph,
        model: &Model,
    ) -> Result<(InitialGraphData, NumGraph)> {
        let global_data = &parse_graph.global_data;

        let initial_data = InitialGraphData {
            additional_params: global_data.parameters.clone(),
            overall_factor: global_data.overall_factor.clone(),
            group_id: global_data.group_id,
            is_group_master: global_data.is_group_master,
            name: global_data.name.clone(),
            threshold_counterterms: Autogen::from_option_or_generate(
                global_data.threshold_counterterms.clone(),
                Default::default,
            ),
        };

        let mut indices = NumIndices::parse(parse_graph, model);
        let num_graph = parse_graph.graph.map_data_ref_result(
            |_, _, v| Ok(v.clone()),
            |_, _, _, e| Ok(e.map(|e| e.clone())),
            |(h, hd)| {
                Ok(HedgeData {
                    num_indices: indices(h, hd)?,
                    ufo_order: Autogen::explicit(hd.ufo_order),
                })
            },
        )?;

        Ok((initial_data, num_graph))
    }

    fn build_underlying_graph(
        graph: NumGraph,
        model: &Model,
        param_builder: &ParamBuilder,
        mut vertex_numerator: impl FnMut(NodeIndex) -> Result<SymbolicTensor<IndexedStructure>>,
        mut edge_numerator: impl FnMut(EdgeIndex) -> Result<SymbolicTensor<IndexedStructure>>,
    ) -> Result<UnderlyingGraph> {
        let intermediate: UnderlyingGraph = graph.map_result(
            |_, i, v| {
                let numerator = vertex_numerator(i)?;
                let dod = match v.dod {
                    Some(dod) => Autogen::explicit(dod),
                    None => Autogen::generated(numerator.expression().all_dod(GS.emr_mom)?),
                };
                let num = Autogen::explicit(numerator);

                Ok(Vertex {
                    name: Autogen::from_option_or_generate(v.name, || i.to_string()),
                    num,
                    dod,
                    vertex_rule: v.vertex_rule,
                })
            },
            |_, _, _, eid, ed| {
                let e = ed.data;
                if e.particle.is_fermion(model)
                    && !e.particle.is_self_antiparticle(model)
                    && e.particle.orientation(model) != ed.orientation
                {
                    return Err(eyre!(
                        "Edge orientation {:?} does not match particle orientation {:?} for edge {eid} ({})",
                        ed.orientation,
                        e.particle.orientation(model),
                        e.particle.particle().resolve(model).name
                    ));
                }

                let mass = EdgeMass::from_atom(e.particle.mass_atom(model), model, param_builder)?;

                let numerator = edge_numerator(eid)?;
                let dod = match e.dod {
                    Some(dod) => Autogen::explicit(dod),
                    None => Autogen::generated(
                        numerator.expression().edge_dod(GS.emr_mom, usize::from(eid))? - 2,
                    ),
                };
                let num = Autogen::explicit(numerator);

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
