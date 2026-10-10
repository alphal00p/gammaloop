use feynkit_graph::DOD;
use std::{collections::HashSet, ops::Deref};

use crate::{
    cff::surface::SurfaceCache,
    graph::{
        FinalizedCut, FinalizedTopologyThresholdCandidate, GraphGroup, GroupId, LoopMomentumBasis,
        edge::EdgeExtraData,
    },
    integrands::process::ParamBuilder,
    model::{Model, ParticleIdGammaLoopExt},
    numerator::{GlobalPrefactor, aind::Aind},
    processes::DotExportSettings,
    utils::GS,
    uv::UltravioletGraph,
};
use feynkit_graph::{
    DiagramEndpoint as FeynkitDiagramEndpoint, DiagramHalfEdge as FeynkitDiagramHalfEdge,
    FeynmanDiagram,
};
use idenso::{
    CookMode, CookSettings, IndexTooling, color::ColorSimplifySettings, tensor::SymbolicTensor,
};
use spenso::shadowing::{TensorCollectFilter, symbolica_utils::LogPrint};

use color_eyre::Result;

use eyre::{Context, Ok, eyre};
use itertools::Itertools;
use linnet::half_edge::{
    HedgeGraph,
    involution::{EdgeData, EdgeIndex, Flow, HedgePair},
    subgraph::{ModifySubSet, OrientedCut, SuBitGraph, SubSetOps},
};
use symbolica::{atom::Atom, id::ConditionResult};
use tracing::warn;
use typed_index_collections::TiVec;

use super::{
    Autogen, Edge, Graph, HedgeData, Vertex,
    edge::{EdgeMass, ParseEdge},
    global::ParseData,
    hedge_data::{NumIndices, ParseHedgeData},
    vertex::ParseVertex,
};

// Type aliases for cleaner code
type NumGraph = HedgeGraph<ParseEdge, ParseVertex, HedgeData>;
type UnderlyingGraph = HedgeGraph<Edge, Vertex, HedgeData>;

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
    /// Attach runtime parse payloads to the finalized half-edge storage.
    /// The map retains every vertex, edge, half-edge, orientation and expression.
    fn from_feynkit_diagram(diagram: &FeynmanDiagram) -> Result<Self> {
        let underlying = diagram.underlying();
        let graph = underlying.map_data_ref_result(
            |_, _, vertex| {
                Ok(ParseVertex {
                    name: Some(vertex.name.clone()),

                    vertex_rule: vertex.interaction,
                    num: Some(vertex.numerator.clone()),
                    dod: None,
                })
            },
            |_, edge_id, _, data| {
                let mut edge = ParseEdge::new(data.data.particle)
                    .with_label(format!("feynkit_edge_{}", edge_id.0))
                    .with_num(data.data.numerator.clone());
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
                name: diagram.name().to_owned(),
                overall_factor: diagram.overall_factor().clone(),
                projector: diagram.projector().clone(),
                num: diagram.numerator_prefactor().clone(),
                ..Default::default()
            },
            graph,
        })
    }
}

/// Helper struct to hold initial data extracted from ParseGraph
struct InitialGraphData {
    overall_factor: Atom,
    global_prefactor: GlobalPrefactor,
    additional_params: Vec<Atom>,
    group_id: Option<GroupId>,
    is_group_master: bool,
    name: String,
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
    fn from_feynkit_parse_graph(diagram: &FeynmanDiagram, parsed: ParseGraph) -> Result<Self> {
        diagram.validate().with_context(|| {
            format!(
                "FeynKit diagram '{}' is not finalized consistently",
                diagram.name()
            )
        })?;

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
        let global_prefactor = initial_data.global_prefactor;
        let polarizations = global_prefactor.polarizations();
        let param_builder = ParamBuilder::new(
            &(&polarizations, &graph),
            model,
            &loop_momentum_basis,
            initial_data.additional_params.clone(),
        );

        let underlying =
            Self::build_underlying_graph(graph, model, &param_builder).with_context(|| {
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
        };
        result.param_builder = ParamBuilder::new(
            &result,
            model,
            &result.loop_momentum_basis,
            initial_data.additional_params,
        );
        let runtime_numerator = result
            .numerator(&result.full_filter(), &result.empty_subgraph())
            .get_single_atom()?;
        if runtime_numerator != diagram.numerator() {
            let cooking = CookSettings::indices()
                .with_mode(CookMode::ReversibleEncoding)
                .with_representation_payloads(true, true);
            let runtime = SymbolicTensor::infer(cooking.try_cook(runtime_numerator.as_view())?)?;
            let finalized =
                SymbolicTensor::infer(cooking.try_cook(diagram.numerator().as_view())?)?;
            // Keep graph numerators factorized, including at this validation
            // boundary. Only an exact coefficient proof certifies a rewrite.
            if runtime.coefficients_equal(&finalized, TensorCollectFilter::<0>::Tensors)?
                != ConditionResult::True
            {
                return Err(eyre!(
                    "GammaLoop runtime conversion of FeynKit diagram '{}' could not certify its finalized numerator",
                    diagram.name()
                ));
            }
        }
        Ok(result)
    }

    pub(crate) fn validate_full_numerator_tensor_network(&self) -> Result<()> {
        let full_num = self
            .numerator(&self.full_filter(), &self.empty_subgraph())
            .get_single_atom()
            .unwrap()
            * &self.global_prefactor.num
            * &self.global_prefactor.projector
            * &self.overall_factor;
        let cooking =
            idenso::CookSettings::indices().with_mode(idenso::CookMode::ReversibleEncoding);
        let color_simplified = idenso::tensor::SymbolicTensor::infer(
            cooking
                .try_cook(full_num.as_view())
                .map_err(|error| eyre!("{error}"))?,
        )?
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
        Self::require_contracted_indices(&full_num, &HashSet::new())
    }

    /// Fail when a full numerator leaves a tensor index other than those in `open`
    /// uncontracted.
    pub(crate) fn require_contracted_indices(
        full_numerator: &Atom,
        open: &HashSet<Atom>,
    ) -> Result<()> {
        let mut dangling = full_numerator.list_dangling::<Aind>()?;
        dangling.retain(|index| !open.contains(index));
        if dangling.is_empty() {
            return Ok(());
        }
        Err(eyre!(
            "Full numerator still has dangling tensor indices: \n{}",
            dangling
                .iter()
                .map(|index| format!("{}:{}", index.log_print(None), index.to_plain_string()))
                .join(",\n")
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
            global_prefactor: GlobalPrefactor {
                num: global_data.num.clone(),
                projector: global_data.projector.clone(),
            },
            group_id: global_data.group_id,
            is_group_master: global_data.is_group_master,
            name: global_data.name.clone(),
        };

        let num_graph = parse_graph.graph.map_data_ref(
            |_, _, v| v.clone(),
            |_, _, _, e| e.map(|e| e.clone()),
            |h, hd| HedgeData {
                num_indices: NumIndices::parse(parse_graph, model)(h, hd),
                ufo_order: Autogen::explicit(hd.ufo_order),
            },
        );

        Ok((initial_data, num_graph))
    }

    fn build_underlying_graph(
        graph: NumGraph,
        model: &Model,
        param_builder: &ParamBuilder,
    ) -> Result<UnderlyingGraph> {
        let intermediate: UnderlyingGraph = graph.map_result(
            |_, i, v| {
                let num = Autogen::explicit(v.num.ok_or_else(|| {
                    eyre!("finalized graph is missing the numerator for vertex {i}")
                })?);

                let dod = match v.dod {
                    Some(dod) => Autogen::explicit(dod),
                    None => Autogen::generated(num.all_dod(GS.emr_mom)?),
                };

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

                let num = Autogen::explicit(e.num.ok_or_else(|| {
                    eyre!("finalized graph is missing the numerator for edge {eid}")
                })?);

                let dod = match e.dod {
                    Some(dod) => Autogen::explicit(dod),
                    None => Autogen::generated(num.edge_dod(GS.emr_mom, usize::from(eid))? - 2),
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
