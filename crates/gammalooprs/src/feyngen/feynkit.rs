//! GammaLoop runtime enrichment for canonical FeynKit generation.
//!
//! FeynKit owns the complete deterministic physics pipeline. GammaLoop only
//! attaches its downstream runtime caches to finalized FeynKit diagrams.

use std::sync::Arc;

use feynkit_generator::{GenerationControl, GenerationReport, GenerationResult, Generator};
use feynkit_graph::{EdgeId, FeynmanDiagram};
use thiserror::Error;

use super::FeynGenError;
use crate::{
    graph::{Graph, GroupId},
    model::Model,
    processes::ProcessDefinition,
    settings::GlobalSettings,
};

/// Errors raised while enriching canonical FeynKit diagrams for GammaLoop.
#[derive(Debug, Error)]
pub enum FeynkitRuntimeError {
    #[error("generation interrupted by user")]
    Interrupted,
    #[error(transparent)]
    Generation(#[from] feynkit_generator::GenerationError),
    #[error("failed to convert FeynKit diagram '{diagram}' into a GammaLoop graph: {source}")]
    GraphConversion {
        diagram: String,
        #[source]
        source: color_eyre::Report,
    },
    #[error("forced cut edge '{edge}' must use the canonical eN edge-ID syntax")]
    InvalidForcedCutEdge { edge: String },
}

/// Downstream GammaLoop behavior for the canonical FeynKit diagram type.
///
/// This is intentionally an extension trait: GammaLoop does not wrap or
/// duplicate FeynKit's graph representation.
pub trait FeynmanDiagramGammaLoopExt {
    fn to_gamma_loop_graph(
        &self,
        group_id: Option<GroupId>,
        is_group_master: bool,
    ) -> Result<Graph, FeynkitRuntimeError>;
}

impl FeynmanDiagramGammaLoopExt for FeynmanDiagram {
    fn to_gamma_loop_graph(
        &self,
        group_id: Option<GroupId>,
        is_group_master: bool,
    ) -> Result<Graph, FeynkitRuntimeError> {
        Graph::from_feynkit(self, group_id, is_group_master).map_err(|source| {
            FeynkitRuntimeError::GraphConversion {
                diagram: self.name().to_owned(),
                source,
            }
        })
    }
}

/// Atomic GammaLoop enrichment for a canonical collapsed generation result.
pub trait GenerationResultGammaLoopExt {
    fn to_gamma_loop_graphs(&self) -> Result<Vec<Graph>, FeynkitRuntimeError>;
}

impl GenerationResultGammaLoopExt for GenerationResult {
    fn to_gamma_loop_graphs(&self) -> Result<Vec<Graph>, FeynkitRuntimeError> {
        self.validate_groups()?;
        let mut group_ids = vec![None; self.diagrams.len()];
        for (group_index, group) in self.groups.iter().enumerate() {
            group_ids[group.master] = Some(GroupId(group_index));
        }
        self.diagrams
            .iter()
            .enumerate()
            .map(|(index, diagram)| diagram.to_gamma_loop_graph(group_ids[index], true))
            .collect()
    }
}

impl ProcessDefinition {
    /// Run canonical FeynKit generation and attach GammaLoop runtime caches.
    #[tracing::instrument(skip_all)]
    pub fn generate_with_report(
        &self,
        model: &Model,
        settings: &GlobalSettings,
    ) -> Result<(Vec<Graph>, GenerationReport), FeynkitRuntimeError> {
        let result = (|| {
            let forced_cuts = settings
                .generation
                .force_cuts
                .iter()
                .map(|cut| {
                    cut.iter()
                        .map(|edge| {
                            edge.strip_prefix('e')
                                .filter(|index| {
                                    !index.is_empty()
                                        && index.bytes().all(|byte| byte.is_ascii_digit())
                                })
                                .and_then(|index| index.parse::<usize>().ok())
                                .map(EdgeId)
                                .ok_or_else(|| FeynkitRuntimeError::InvalidForcedCutEdge {
                                    edge: edge.clone(),
                                })
                        })
                        .collect::<Result<Vec<_>, _>>()
                })
                .collect::<Result<Vec<_>, _>>()?;
            let options = self
                .generation_options
                .clone()
                .threads(settings.n_cores.feyngen)
                .forced_cuts(forced_cuts)
                .cancellation_check(crate::is_interrupt_requested)
                .progress(|_| {
                    if crate::is_interrupt_requested() {
                        GenerationControl::Cancel
                    } else {
                        GenerationControl::Continue
                    }
                });
            let generated =
                Generator::new(Arc::new(model.clone())).generate(&self.process, &options)?;
            if !generated.report.completed || crate::is_interrupt_requested() {
                return Err(FeynkitRuntimeError::Interrupted);
            }
            let graphs = generated.to_gamma_loop_graphs()?;
            Ok((graphs, generated.report))
        })();
        crate::clear_interrupt_request();
        result
    }

    /// Generate GammaLoop runtime graphs through the canonical FeynKit engine.
    pub fn generate(
        &self,
        model: &Model,
        settings: &GlobalSettings,
    ) -> Result<Vec<Graph>, FeynGenError> {
        self.generate_with_report(model, settings)
            .map(|(graphs, _)| graphs)
            .map_err(Into::into)
    }
}

#[cfg(test)]
mod tests {
    use feynkit_generator::{GenerationOptions, Process};
    use feynkit_graph::{
        DiagramCut, DiagramCutSide, DiagramEdge, DiagramEndpoint, DiagramHalfEdge, DiagramVertex,
        EdgeId, ExternalState, FeynmanDiagram, LoopMomentumBasis,
    };
    use linnet::half_edge::{
        involution::{EdgeIndex, Hedge, HedgePair, Orientation},
        subgraph::{Inclusion, SuBitGraph, SubGraphLike, SubSetLike, SubSetOps},
    };
    use symbolica::{
        atom::{Atom, AtomCore},
        parser::ParseSettings,
    };

    use super::*;
    use crate::{
        graph::LMBext,
        momentum::{SignOrZero, sample::ExternalIndex},
        numerator::aind::NewAind,
        processes::{CrossSectionGraph, ProcessDefinition},
        settings::global::GenerationSettings,
        uv::UltravioletGraph,
    };

    fn atom(expression: &str) -> Atom {
        Atom::parse(expression, "feynkit_runtime_test", ParseSettings::default()).unwrap()
    }

    fn model() -> Arc<Model> {
        Arc::new(
            Model::from_json(
                r#"{
                    "name":"direct_runtime",
                    "restriction":null,
                    "orders":[],
                    "parameters":[
                        {"name":"ZERO","lhablock":null,"lhacode":null,"nature":"internal","parameter_type":"real","value":[0.0,0.0],"expression":null}
                    ],
                    "particles":[
                        {"pdg_code":25,"name":"phi","antiname":"phi","spin":1,"color":1,"mass":"ZERO","width":"ZERO","texname":"phi","antitexname":"phi","charge":0.0,"ghost_number":0,"lepton_number":0,"y_charge":0}
                    ],
                    "propagators":[
                        {"name":"phi_prop","particle":"phi","numerator":"model_edge_rule","denominator":"P^2"}
                    ],
                    "lorentz_structures":[
                        {"name":"L","spins":[1,1,1],"structure":"model_vertex_rule"}
                    ],
                    "couplings":[
                        {"name":"GC","expression":"model_coupling","orders":[],"value":null}
                    ],
                    "vertex_rules":[
                        {"name":"V","particles":["phi","phi","phi"],"color_structures":["1"],"lorentz_structures":["L"],"couplings":[["GC"]]}
                    ],
                    "functions":[],
                    "form_factors":[]
                }"#,
            )
            .unwrap(),
        )
    }

    fn diagram(model: &Arc<Model>, basis: Option<LoopMomentumBasis>) -> FeynmanDiagram {
        diagram_with_vertex_dummy(model, basis, false)
    }

    fn diagram_with_vertex_dummy(
        model: &Arc<Model>,
        basis: Option<LoopMomentumBasis>,
        vertex_dummy: bool,
    ) -> FeynmanDiagram {
        let particle = model.particle_id("phi").unwrap();
        let interaction = model.vertex_rule_id("V").unwrap();
        let left_numerator = if vertex_dummy {
            atom("FeynKit::VertexDummy(1,0)")
        } else {
            atom("kept_vertex_fragment")
        };
        let numerator = if vertex_dummy {
            atom("FeynKit::VertexDummy(1,0)*kept_vertex_fragment*kept_edge_fragment^2")
        } else {
            atom("kept_vertex_fragment^2*kept_edge_fragment^2")
        };
        let mut builder = FeynmanDiagram::builder(Arc::clone(model), "direct")
            .symmetry_factor(7)
            .overall_factor(atom("kept_overall"))
            .numerator(numerator)
            .numerator_prefactor(atom("kept_prefactor"))
            .projector(atom("kept_projector"));
        if let Some(basis) = basis {
            builder = builder.loop_momentum_basis(basis);
        }
        let incoming = builder.add_vertex(DiagramVertex::external(
            "incoming",
            0,
            ExternalState::Incoming,
        ));
        let left = builder.add_vertex(DiagramVertex {
            numerator: left_numerator,
            ..DiagramVertex::interaction("left", interaction)
        });
        let right = builder.add_vertex(DiagramVertex {
            numerator: atom("kept_vertex_fragment"),
            ..DiagramVertex::interaction("right", interaction)
        });
        let outgoing = builder.add_vertex(DiagramVertex::external(
            "outgoing",
            1,
            ExternalState::Outgoing,
        ));
        builder
            .add_edge(incoming, left, DiagramEdge::new(particle, false))
            .unwrap();
        for _ in 0..2 {
            let mut edge = DiagramEdge::new(particle, false);
            edge.numerator = atom("kept_edge_fragment");
            builder.add_edge(left, right, edge).unwrap();
        }
        builder
            .add_edge(right, outgoing, DiagramEdge::new(particle, false))
            .unwrap();
        builder.build().unwrap()
    }

    fn legacy_order_bubble(model: &Arc<Model>, basis: Option<LoopMomentumBasis>) -> FeynmanDiagram {
        let particle = model.particle_id("phi").unwrap();
        let interaction = model.vertex_rule_id("V").unwrap();
        let mut builder = FeynmanDiagram::builder(Arc::clone(model), "legacy_order_bubble")
            .numerator(atom("kept_vertex_fragment^4"));
        if let Some(basis) = basis {
            builder = builder.loop_momentum_basis(basis);
        }
        let incoming = builder.add_vertex(DiagramVertex::external_in_connection(
            "incoming",
            0,
            ExternalState::Incoming,
            0,
        ));
        let interaction_vertex = |name| DiagramVertex {
            numerator: atom("kept_vertex_fragment"),
            ..DiagramVertex::interaction(name, interaction)
        };
        let first = builder.add_vertex(interaction_vertex("first"));
        let second = builder.add_vertex(interaction_vertex("second"));
        let third = builder.add_vertex(interaction_vertex("third"));
        let fourth = builder.add_vertex(interaction_vertex("fourth"));
        let outgoing = builder.add_vertex(DiagramVertex::external_in_connection(
            "outgoing",
            1,
            ExternalState::Outgoing,
            0,
        ));

        for (expected, source, target) in [
            (0, incoming, first),
            (1, first, second),
            (2, second, third),
            (3, second, third),
            (4, third, fourth),
            (5, first, fourth),
            (6, fourth, outgoing),
        ] {
            assert_eq!(
                builder
                    .add_edge(source, target, DiagramEdge::new(particle, false))
                    .unwrap(),
                EdgeId(expected)
            );
        }
        let diagram = builder.build().unwrap();
        let half_edge = |edge, endpoint| DiagramHalfEdge {
            edge: EdgeId(edge),
            endpoint,
        };
        let left = vec![
            half_edge(0, DiagramEndpoint::Source),
            half_edge(0, DiagramEndpoint::Target),
            half_edge(1, DiagramEndpoint::Source),
            half_edge(5, DiagramEndpoint::Source),
        ];
        let right = (0..=6)
            .flat_map(|edge| {
                [
                    half_edge(edge, DiagramEndpoint::Source),
                    half_edge(edge, DiagramEndpoint::Target),
                ]
            })
            .filter(|half_edge| !left.contains(half_edge))
            .collect::<Vec<_>>();
        diagram
            .with_cut_partitions(vec![(left.clone(), right.clone())])
            .unwrap()
            .with_topology_threshold_partitions(vec![(left, right)])
            .unwrap()
    }

    fn two_connection_cross_section(model: &Arc<Model>) -> FeynmanDiagram {
        let particle = model.particle_id("phi").unwrap();
        let interaction = model.vertex_rule_id("V").unwrap();
        let mut builder = FeynmanDiagram::builder(Arc::clone(model), "two_connections");
        let incoming_0 = builder.add_vertex(DiagramVertex::external_in_connection(
            "incoming_0",
            0,
            ExternalState::Incoming,
            0,
        ));
        let incoming_1 = builder.add_vertex(DiagramVertex::external_in_connection(
            "incoming_1",
            1,
            ExternalState::Incoming,
            1,
        ));
        let outgoing_0 = builder.add_vertex(DiagramVertex::external_in_connection(
            "outgoing_0",
            2,
            ExternalState::Outgoing,
            0,
        ));
        let outgoing_1 = builder.add_vertex(DiagramVertex::external_in_connection(
            "outgoing_1",
            3,
            ExternalState::Outgoing,
            1,
        ));
        let internal = (0..4)
            .map(|index| {
                builder.add_vertex(DiagramVertex::interaction(format!("v{index}"), interaction))
            })
            .collect::<Vec<_>>();
        let scalar = || DiagramEdge::new(particle, false);

        // The second incoming edge points from its interaction vertex toward
        // the incoming external vertex. The sewing bridge must nevertheless
        // expose that connection as +P1 in the runtime coordinate frame.
        builder.add_edge(incoming_0, internal[0], scalar()).unwrap();
        builder.add_edge(internal[1], incoming_1, scalar()).unwrap();
        builder.add_edge(internal[2], outgoing_0, scalar()).unwrap();
        builder.add_edge(internal[3], outgoing_1, scalar()).unwrap();
        builder
            .add_edge(internal[0], internal[1], scalar())
            .unwrap();
        builder
            .add_edge(internal[1], internal[2], scalar())
            .unwrap();
        builder
            .add_edge(internal[2], internal[3], scalar())
            .unwrap();
        builder
            .add_edge(internal[3], internal[0], scalar())
            .unwrap();

        let half_edge = |edge, endpoint| DiagramHalfEdge {
            edge: EdgeId(edge),
            endpoint,
        };
        let left = vec![
            half_edge(0, DiagramEndpoint::Source),
            half_edge(0, DiagramEndpoint::Target),
            half_edge(1, DiagramEndpoint::Source),
            half_edge(1, DiagramEndpoint::Target),
            half_edge(4, DiagramEndpoint::Source),
            half_edge(4, DiagramEndpoint::Target),
            half_edge(5, DiagramEndpoint::Source),
            half_edge(7, DiagramEndpoint::Target),
        ];
        let right = (0..8)
            .flat_map(|edge| {
                [
                    half_edge(edge, DiagramEndpoint::Source),
                    half_edge(edge, DiagramEndpoint::Target),
                ]
            })
            .filter(|half_edge| !left.contains(half_edge))
            .collect::<Vec<_>>();
        let side = |half_edges| DiagramCutSide {
            half_edges,
            coupling_orders: Default::default(),
            loop_count: 0,
        };
        builder
            .cuts(vec![DiagramCut {
                cut: vec![
                    half_edge(5, DiagramEndpoint::Source),
                    half_edge(7, DiagramEndpoint::Target),
                ],
                left: side(left),
                right: side(right),
            }])
            .build()
            .unwrap()
    }

    #[test]
    fn direct_conversion_preserves_finalized_physics_and_routing() {
        let model = model();
        let draft = diagram(&model, None);
        let basis = draft
            .loop_momentum_bases()
            .unwrap()
            .into_iter()
            .find(|basis| basis.loop_edges == vec![EdgeId(2)])
            .expect("the bubble topology exposes edge 2 as a loop basis");
        let selected_loop_edge = basis.loop_edges[0];
        let diagram = diagram(&model, Some(basis));

        let runtime = diagram.to_gamma_loop_graph(None, true).unwrap();

        assert_eq!(diagram.symmetry_factor(), 7);
        assert_eq!(runtime.overall_factor, atom("kept_overall"));
        assert_eq!(runtime.global_prefactor.num, atom("kept_prefactor"));
        assert_eq!(runtime.global_prefactor.projector, atom("kept_projector"));
        assert!(
            runtime
                .underlying
                .iter_nodes()
                .any(|(_, _, vertex)| vertex.num.value == atom("kept_vertex_fragment"))
        );
        assert!(
            runtime
                .underlying
                .iter_edges()
                .any(|(_, _, edge)| { edge.data.num.value == atom("kept_edge_fragment") })
        );
        assert_eq!(
            runtime
                .numerator(&runtime.full_filter(), &runtime.empty_subgraph())
                .get_single_atom()
                .unwrap()
                .expand(),
            diagram.numerator().expand()
        );

        let routed_edge = runtime
            .edge_name_to_index(&format!("feynkit_edge_{}", selected_loop_edge.0))
            .unwrap();
        assert_eq!(
            runtime
                .loop_momentum_basis
                .loop_edges
                .iter()
                .copied()
                .collect::<Vec<_>>(),
            vec![routed_edge]
        );
    }

    #[test]
    fn direct_conversion_preserves_legacy_hedges_and_stable_edge_ids_separately() {
        let model = model();
        let draft = legacy_order_bubble(&model, None);
        let basis = draft
            .loop_momentum_bases()
            .unwrap()
            .into_iter()
            .find(|basis| basis.loop_edges == vec![EdgeId(1), EdgeId(2)])
            .expect("the two-loop bubble exposes the requested legacy basis");
        let runtime = legacy_order_bubble(&model, Some(basis))
            .to_gamma_loop_graph(None, true)
            .unwrap();

        for edge in 0..6 {
            assert_eq!(
                runtime
                    .edge_name_to_index(&format!("feynkit_edge_{edge}"))
                    .unwrap(),
                EdgeIndex(edge),
                "logical EdgeIndex must remain the stable FeynKit edge ID"
            );
        }
        let expected_pairs = [
            HedgePair::Paired {
                source: Hedge(11),
                sink: Hedge(0),
            },
            HedgePair::Paired {
                source: Hedge(1),
                sink: Hedge(2),
            },
            HedgePair::Paired {
                source: Hedge(5),
                sink: Hedge(6),
            },
            HedgePair::Paired {
                source: Hedge(7),
                sink: Hedge(8),
            },
            HedgePair::Paired {
                source: Hedge(9),
                sink: Hedge(10),
            },
            HedgePair::Paired {
                source: Hedge(3),
                sink: Hedge(4),
            },
        ];
        for (edge, expected) in expected_pairs.into_iter().enumerate() {
            assert_eq!(runtime.underlying[&EdgeIndex(edge)].1, expected);
        }
        assert_eq!(
            runtime
                .loop_momentum_basis
                .tree
                .included_iter()
                .collect::<Vec<_>>(),
            [3, 4, 7, 8, 9, 10].map(Hedge).to_vec()
        );
        assert_eq!(
            runtime
                .loop_momentum_basis
                .loop_edges
                .iter()
                .copied()
                .collect::<Vec<_>>(),
            vec![EdgeIndex(1), EdgeIndex(2)]
        );

        // Cross-section channel lookup keys each generated LMB by its sorted
        // loop-edge set but preserves the first basis's internal loop order.
        // The legacy topology order must therefore keep channel {e3, e5} in
        // coordinates (e5, e3), not (e3, e5).
        let full_filter = runtime.full_filter();
        let cut_graph = full_filter.subtract(&runtime.initial_state_cut.right);
        let externals: SuBitGraph = runtime.empty_subgraph();
        let channel_lmb = runtime
            .all_spanning_forests_of(&cut_graph)
            .into_iter()
            .map(|forest| {
                let mut lmb = runtime
                    .lmb_impl(&full_filter, &forest, externals.clone())
                    .unwrap();
                let mut external_loop_edges = lmb
                    .loop_edges
                    .iter()
                    .copied()
                    .filter(|edge| runtime.initial_state_cut.intersects(&runtime[edge].1))
                    .collect::<Vec<_>>();
                external_loop_edges.sort();
                for edge in external_loop_edges {
                    let loop_index = lmb.edge_signatures[edge]
                        .internal
                        .iter_enumerated()
                        .find_map(|(loop_index, sign)| {
                            (*sign != crate::momentum::SignOrZero::Zero).then_some(loop_index)
                        })
                        .expect("an initial-state basis edge carries one loop coordinate");
                    lmb.put_loop_to_ext(loop_index);
                }
                runtime.canonicalize_lmb_external_order(&mut lmb);
                lmb
            })
            .find(|lmb| {
                let mut loop_edges = lmb.loop_edges.iter().copied().collect::<Vec<_>>();
                loop_edges.sort();
                loop_edges == vec![EdgeIndex(3), EdgeIndex(5)]
            })
            .expect("the bubble exposes the legacy {e3, e5} channel");
        assert_eq!(channel_lmb.loop_edges.raw, vec![EdgeIndex(5), EdgeIndex(3)]);
    }

    #[test]
    fn sewn_initial_state_routing_uses_positive_external_momenta() {
        let model = model();
        let runtime = two_connection_cross_section(&model)
            .to_gamma_loop_graph(None, true)
            .unwrap();

        assert_eq!(
            runtime.loop_momentum_basis.ext_edges.raw,
            vec![EdgeIndex(0), EdgeIndex(1)]
        );
        assert_eq!(
            runtime.loop_momentum_basis.edge_signatures[EdgeIndex(0)]
                .external
                .iter()
                .copied()
                .collect::<Vec<_>>(),
            vec![SignOrZero::Plus, SignOrZero::Zero]
        );
        assert_eq!(
            runtime.loop_momentum_basis.edge_signatures[EdgeIndex(1)]
                .external
                .iter()
                .copied()
                .collect::<Vec<_>>(),
            vec![SignOrZero::Zero, SignOrZero::Plus]
        );

        let node = |name: &str| {
            runtime
                .underlying
                .iter_nodes()
                .find(|(_, _, vertex)| vertex.name.value == name)
                .map(|(node, _, _)| node)
                .unwrap_or_else(|| panic!("missing runtime vertex {name}"))
        };
        for (edge, incoming, outgoing) in [
            (EdgeIndex(0), node("v0"), node("v2")),
            (EdgeIndex(1), node("v1"), node("v3")),
        ] {
            let HedgePair::Paired { source, sink } = runtime.underlying[&edge].1 else {
                panic!("sewn connection {edge} is not paired");
            };
            assert_eq!(runtime.underlying.node_id(source), incoming);
            assert_eq!(runtime.underlying.node_id(sink), outgoing);
        }
    }

    #[test]
    fn direct_conversion_translates_vertex_dummy_by_typed_vertex_identity() {
        let model = model();
        let diagram = diagram_with_vertex_dummy(&model, None, true);
        let runtime = diagram.to_gamma_loop_graph(None, true).unwrap();

        let (left_node, _, left) = runtime
            .underlying
            .iter_nodes()
            .find(|(_, _, vertex)| vertex.name.value == "left")
            .expect("the left FeynKit vertex must survive runtime enrichment");
        assert_eq!(left.num.value, Atom::from(left_node.aind(0)));
    }

    #[test]
    fn direct_conversion_rejects_inconsistent_diagram_numerator() {
        let model = model();
        let diagram = diagram(&model, None);
        let mut builder = FeynmanDiagram::builder(Arc::clone(&model), diagram.name())
            .symmetry_factor(diagram.symmetry_factor())
            .overall_factor(diagram.overall_factor().clone())
            .numerator(atom("different_numerator"))
            .numerator_prefactor(diagram.numerator_prefactor().clone())
            .projector(diagram.projector().clone())
            .loop_momentum_basis(diagram.loop_momentum_basis().clone())
            .cuts(diagram.cuts().to_vec());
        for (_, vertex) in diagram.vertices() {
            builder.add_vertex(vertex.clone());
        }
        for (_, endpoints, edge) in diagram.edges() {
            builder
                .add_edge_with_slots(
                    endpoints.source,
                    endpoints.target,
                    edge.clone(),
                    edge.source_slot(),
                    edge.target_slot(),
                )
                .unwrap();
        }
        let inconsistent = builder.build().unwrap();
        let error = match inconsistent.to_gamma_loop_graph(None, true) {
            Ok(_) => panic!("an inconsistent numerator must be rejected"),
            Err(error) => error,
        };
        let FeynkitRuntimeError::GraphConversion { source, .. } = error else {
            panic!("unexpected conversion error variant: {error}");
        };
        let message = format!("{source:?}");
        assert!(
            message.contains("diagram-wide numerator differs"),
            "unexpected conversion error: {message}"
        );
    }

    #[test]
    fn canonical_cross_section_generation_survives_runtime_enrichment_and_consumption() {
        let model = model();
        let process = Process::cross_section(["phi"], ["phi", "phi"])
            .with_loop_count(1, 1)
            .unwrap();
        let options = GenerationOptions::default().threads(1).max_vertices(2);
        let generated = Generator::new(Arc::clone(&model))
            .generate(&process, &options)
            .unwrap();

        assert!(generated.report.completed);
        assert!(!generated.diagrams.is_empty());
        assert!(
            generated
                .diagrams
                .iter()
                .all(|diagram| !diagram.cuts().is_empty())
        );
        assert!(
            generated
                .diagrams
                .iter()
                .all(|diagram| !diagram.topology_threshold_candidates().is_empty())
        );

        let runtime_graphs = generated.to_gamma_loop_graphs().unwrap();
        assert_eq!(runtime_graphs.len(), generated.diagrams.len());
        let process_definition = ProcessDefinition {
            process,
            generation_options: options,
            folder_name: "canonical_cross_section_bridge".to_owned(),
            process_id: 0,
        };
        let settings = GenerationSettings::default();

        for (diagram, runtime) in generated.diagrams.iter().zip(runtime_graphs) {
            assert_eq!(runtime.finalized_cuts.len(), diagram.cuts().len());
            assert_eq!(
                runtime.finalized_topology_threshold_candidates.len(),
                diagram.topology_threshold_candidates().len()
            );

            let connections = diagram
                .vertices()
                .filter_map(|(_, vertex)| vertex.external.as_ref())
                .map(|external| external.connection)
                .collect::<std::collections::BTreeSet<_>>();
            assert_eq!(
                runtime.loop_momentum_basis.ext_edges.len(),
                connections.len()
            );
            assert_eq!(
                runtime.initial_state_cut.left.n_included(),
                connections.len()
            );
            for connection in connections {
                let runtime_external =
                    runtime.loop_momentum_basis.ext_edges[ExternalIndex(connection)];
                assert!(
                    runtime
                        .initial_state_cut
                        .as_subgraph()
                        .includes(&runtime.underlying[&runtime_external].1)
                );
            }

            assert!(
                runtime
                    .finalized_cuts
                    .windows(2)
                    .all(|cuts| cuts[0].cut <= cuts[1].cut),
                "runtime cuts must use the legacy canonical oriented-cut order"
            );
            for canonical in diagram.cuts() {
                let canonical_edges = canonical
                    .cut
                    .iter()
                    .map(|half_edge| {
                        runtime
                            .edge_name_to_index(&format!("feynkit_edge_{}", half_edge.edge.0))
                            .unwrap()
                    })
                    .collect::<std::collections::BTreeSet<_>>();
                let finalized = runtime
                    .finalized_cuts
                    .iter()
                    .find(|finalized| {
                        finalized
                            .cut
                            .left
                            .included_iter()
                            .map(|hedge| runtime.underlying[&hedge])
                            .collect::<std::collections::BTreeSet<_>>()
                            == canonical_edges
                    })
                    .expect("each canonical crossing edge set survives runtime enrichment");
                assert!(!finalized.left.intersects(&runtime.initial_state_cut.right));
                assert!(!finalized.right.intersects(&runtime.initial_state_cut.left));
                assert!(finalized.left.intersects(&runtime.initial_state_cut.left));
                assert!(finalized.right.intersects(&runtime.initial_state_cut.right));
                assert_eq!(
                    finalized.cut.left.n_included(),
                    canonical.cut.len(),
                    "each canonical crossing endpoint must remain one oriented runtime cut edge"
                );
                let mut globally_reversed = None;
                for half_edge in &canonical.cut {
                    let runtime_edge = runtime
                        .edge_name_to_index(&format!("feynkit_edge_{}", half_edge.edge.0))
                        .unwrap();
                    let expected = match half_edge.endpoint {
                        DiagramEndpoint::Source => Orientation::Default,
                        DiagramEndpoint::Target => Orientation::Reversed,
                    };
                    let actual = finalized
                        .cut
                        .get_from_pair(runtime.underlying[&runtime_edge].1);
                    let reversed = match (expected, actual) {
                        (Orientation::Default, Orientation::Default)
                        | (Orientation::Reversed, Orientation::Reversed) => false,
                        (Orientation::Default, Orientation::Reversed)
                        | (Orientation::Reversed, Orientation::Default) => true,
                        _ => panic!("a physical cut edge became unoriented"),
                    };
                    assert_eq!(
                        *globally_reversed.get_or_insert(reversed),
                        reversed,
                        "initial-state alignment must reverse every crossing edge together"
                    );
                }
            }

            for canonical_loop_edge in &diagram.loop_momentum_basis().loop_edges {
                let runtime_loop_edge = runtime
                    .edge_name_to_index(&format!("feynkit_edge_{}", canonical_loop_edge.0))
                    .unwrap();
                assert!(
                    runtime
                        .loop_momentum_basis
                        .loop_edges
                        .iter()
                        .any(|edge| *edge == runtime_loop_edge)
                );
            }

            for candidate in &runtime.finalized_topology_threshold_candidates {
                assert!(
                    !candidate
                        .cut
                        .left
                        .intersects(&runtime.initial_state_cut.left)
                );
                assert!(
                    !candidate
                        .cut
                        .left
                        .intersects(&runtime.initial_state_cut.right)
                );
                assert!(
                    !candidate
                        .cut
                        .right
                        .intersects(&runtime.initial_state_cut.left)
                );
                assert!(
                    !candidate
                        .cut
                        .right
                        .intersects(&runtime.initial_state_cut.right)
                );
                assert!(!candidate.left.intersects(&runtime.initial_state_cut.right));
                assert!(!candidate.right.intersects(&runtime.initial_state_cut.left));
                assert!(candidate.cut.nedges(&runtime.underlying) > 1);
            }

            let cross_section = CrossSectionGraph::new(runtime);
            assert!(cross_section.graph.finalized_cuts.is_empty());
            assert!(
                cross_section
                    .graph
                    .finalized_topology_threshold_candidates
                    .is_empty()
            );
            assert_eq!(cross_section.cuts.len(), diagram.cuts().len());
            assert_eq!(
                cross_section.topology_threshold_candidates.len(),
                diagram.topology_threshold_candidates().len()
            );
            assert_eq!(
                cross_section
                    .process_valid_cuts(&model, &process_definition, &settings)
                    .unwrap()
                    .len(),
                diagram.cuts().len()
            );
        }
    }
}
