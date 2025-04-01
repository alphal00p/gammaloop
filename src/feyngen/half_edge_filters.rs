use std::{fmt::Display, sync::Arc};

use bitvec::vec::BitVec;
use linnet::half_edge::{
    involution::{EdgeData, EdgeIndex, Flow, HedgePair, Orientation},
    nodestorage::NodeStorageOps,
    subgraph::{cut::PossiblyCutEdge, SubGraph, SubGraphOps},
    EdgeAccessors, HedgeGraph,
};
use symbolica::graph::Graph as SymbolicaGraph;

use crate::{
    graph::HedgeGraphExt,
    model::{Model, Particle},
};

use super::diagram_generator::{EdgeColor, NodeColorFunctions};

#[derive(Debug, Clone, Eq, PartialEq, Hash, PartialOrd, Ord, Copy)]
pub enum VertexType<V> {
    Internal(V),
    Left { id: EdgeIndex, vertex: V },
    Right { id: EdgeIndex, vertex: V },
}

#[derive(Debug, Clone, Eq, PartialEq, Hash, PartialOrd, Ord, Copy)]
pub enum CutState {
    Saturated,
    Cut,
    Vacuum,
}

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct FeynGenHedgeGraph<E, V> {
    pub graph: HedgeGraph<PossiblyCutEdge<E>, VertexType<V>>,
    pub state: CutState,
    // signature: BTreeSet<(isize, Orientation, EdgeIndex)>,
}

impl<E, V> FeynGenHedgeGraph<E, V> {
    pub fn cut(&mut self)
    where
        E: Clone,
    {
        // -> HedgeGraph<isize, DisCompVertex> {
        //

        self.graph.split_clone();
        self.state = CutState::Cut;
    }

    pub fn glue_external_hedges(&mut self) {
        if self.state == CutState::Cut {
            self.graph.glue_back_strict();
            self.state = CutState::Vacuum;
        }
    }

    /// Only on cut graphs
    pub fn only_default_orientations(&mut self) {
        let cut = self.graph.cut();

        for l in cut.iter_left_hedges() {
            let orientation = self.graph.orientation(l);
            let flow = self.graph.flow(l);

            match orientation {
                Orientation::Reversed => match flow {
                    Flow::Source => self.graph.set_flow(l, Flow::Sink),
                    Flow::Sink => self.graph.set_flow(l, Flow::Source),
                },
                Orientation::Undirected => {
                    if let Flow::Sink = flow {
                        self.graph.set_flow(l, Flow::Source);
                        self.graph[[&l]].cut(Flow::Source);
                    }
                    self.graph.set_orientation(l, Orientation::Default);
                }
                _ => {}
            }
        }

        for l in cut.iter_right_hedges() {
            let orientation = self.graph.orientation(l);
            let flow = self.graph.flow(l);

            match orientation {
                Orientation::Reversed => match flow {
                    Flow::Source => self.graph.set_flow(l, Flow::Sink),
                    Flow::Sink => self.graph.set_flow(l, Flow::Source),
                },
                Orientation::Undirected => {
                    if let Flow::Sink = flow {
                        self.graph.set_flow(l, Flow::Sink);
                        self.graph[[&l]].cut(Flow::Sink);
                    }
                    self.graph.set_orientation(l, Orientation::Default);
                }
                _ => {}
            }
        }
    }

    pub fn remove_external_nodes(&mut self)
    where
        E: Clone,
        V: Clone,
    {
        if self.state != CutState::Saturated {
            return;
        }
        let mut excised: BitVec = self.graph.empty_subgraph();

        for (n, d) in self.graph.iter_nodes() {
            if !matches!(d, VertexType::Internal(_)) {
                excised.union_with(&n.hairs)
            }
        }

        excised = excised.complement(&self.graph);

        let mut excised = self
            .graph
            .concretize(&excised)
            .map(|_, _, _, d| d.clone(), |_, _, _, e| e.map(|d| d.clone()));

        excised.align_underlying_to_superficial();
        self.graph = excised;
        self.only_default_orientations();
        self.state = CutState::Cut;
    }
}
impl<V> FeynGenHedgeGraph<Arc<Particle>, V> {
    pub fn number_of_fermion_loops(&self) -> usize {
        let mut fermions: BitVec = self.graph.empty_subgraph();

        for (p, _, d) in self.graph.iter_all_edges() {
            if d.data.edge_data().is_fermion() {
                #[allow(clippy::single_match)]
                match p {
                    HedgePair::Paired { source, sink } => {
                        fermions.set(source.0, true);
                        fermions.set(sink.0, true);
                    }
                    _ => {}
                }
            }
        }

        self.loop_count(&fermions)
    }

    pub fn loop_count<S: SubGraph>(&self, subgraph: &S) -> usize {
        let n_hedges = self.graph.count_internal_edges(subgraph);

        let n_nodes = self.graph.number_of_nodes_in_subgraph(subgraph);

        let n_components = self.graph.count_connected_components(subgraph);

        n_hedges + n_components - n_nodes
    }

    ///Mutates because it glues external hedges back together
    pub fn number_of_external_fermion_loops(&mut self) -> usize
    where
        V: Clone,
    {
        self.remove_external_nodes();
        let internal = self.number_of_fermion_loops();
        // println!("{}", self);

        self.glue_external_hedges();

        let all_fermion_loops = self.number_of_fermion_loops();

        all_fermion_loops - internal
    }

    pub fn from_feyn_gen_symbolica(
        graph: SymbolicaGraph<V, EdgeColor>,
        model: &Model,
        n_initials: usize,
    ) -> Self
    where
        V: NodeColorFunctions + Clone,
    {
        let mut graph = HedgeGraph::<EdgeColor, V>::from_sym(graph).map(
            |_, _, _, vertex| {
                let id = EdgeIndex::from(vertex.pairing_tag(n_initials) as usize);
                if vertex.is_incoming(n_initials) {
                    VertexType::Right { id, vertex }
                } else if vertex.is_outgoing(n_initials) {
                    VertexType::Left { id, vertex }
                } else {
                    VertexType::Internal(vertex)
                }
            },
            |inv, _, pair, d| {
                let mut particle = model.get_particle_from_pdg(d.data.pdg);
                let mut o = d.orientation;
                if !particle.is_fermion() {
                    o = Orientation::Undirected;
                } else if particle.is_antiparticle() {
                    o = Orientation::Reversed;
                    particle = particle.get_anti_particle(model);
                }
                let id = inv[pair.any_hedge()];
                EdgeData::new(PossiblyCutEdge::uncut(particle, id), o)
            },
        );

        graph = graph.map(
            |_, _, _, v| v,
            |_, node_data, pair, mut e| {
                if let HedgePair::Paired { source, sink } = pair {
                    let src = node_data.node_id_ref(source);
                    let sk = node_data.node_id_ref(sink);

                    let source_node = node_data.get_node_data(src);
                    let sink_node = node_data.get_node_data(sk);

                    match (source_node, sink_node) {
                        (VertexType::Left { id, .. }, _) => {
                            e.data.cut(Flow::Sink);
                            e.data.index = *id;
                        }
                        (VertexType::Right { id, .. }, _) => {
                            e.data.cut(Flow::Source);
                            e.data.index = *id;
                        }
                        (_, VertexType::Left { id, .. }) => {
                            e.data.cut(Flow::Sink);
                            e.data.index = *id;
                        }
                        (_, VertexType::Right { id, .. }) => {
                            e.data.cut(Flow::Source);
                            e.data.index = *id;
                        }
                        _ => {}
                    }
                };
                e
            },
        );
        Self {
            graph,
            state: CutState::Saturated,
        }
    }
}

impl<V> Display for FeynGenHedgeGraph<Arc<Particle>, V> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}",
            self.graph.dot_impl(
                &self.graph.full_filter(),
                "",
                &|a| {
                    match a.flow() {
                        Some(Flow::Sink) => {
                            Some(format!("label=\"sink{}{}\"", a.index, a.edge_data().name))
                        }
                        Some(Flow::Source) => {
                            Some(format!("label=\"source{}{}\"", a.index, a.edge_data().name))
                        }
                        None => Some(format!("label=\"{}\"", a.edge_data().name)),
                    }
                },
                &|a| Some(format!(
                    "label=\"{}\"",
                    match a {
                        VertexType::Internal(_) => "i".to_owned(),
                        VertexType::Right { id, .. } => format!("r{}", id),
                        VertexType::Left { id, .. } => format!("l{}", id),
                    }
                ))
            )
        )?;
        Ok(())
    }
}

#[cfg(test)]
pub mod test {
    use symbolica::graph::Graph as SymbolicaGraph;

    use crate::{
        feyngen::diagram_generator::{EdgeColor, NodeColorWithoutVertexRule},
        tests_from_pytest::load_generic_model,
    };

    use super::FeynGenHedgeGraph;

    #[test]
    fn fermion_loop_count() {
        let mut graph = SymbolicaGraph::new();

        let ext1 = graph.add_node(NodeColorWithoutVertexRule { external_tag: 1 });
        let ext2 = graph.add_node(NodeColorWithoutVertexRule { external_tag: 4 });

        let ext3 = graph.add_node(NodeColorWithoutVertexRule { external_tag: 3 });
        let ext4 = graph.add_node(NodeColorWithoutVertexRule { external_tag: 2 });

        let internal = NodeColorWithoutVertexRule { external_tag: 0 };
        let n1 = graph.add_node(internal);
        let n2 = graph.add_node(internal);
        let n3 = graph.add_node(internal);
        let n4 = graph.add_node(internal);
        let n5 = graph.add_node(internal);
        let n6 = graph.add_node(internal);

        graph
            .add_edge(ext1, n1, true, EdgeColor { pdg: 11 })
            .unwrap();
        graph.add_edge(n1, n2, true, EdgeColor { pdg: 11 }).unwrap();
        graph
            .add_edge(n2, n3, false, EdgeColor { pdg: 21 })
            .unwrap();
        graph.add_edge(n4, n3, true, EdgeColor { pdg: 11 }).unwrap();
        graph
            .add_edge(n4, n5, false, EdgeColor { pdg: 21 })
            .unwrap();
        graph.add_edge(n5, n6, true, EdgeColor { pdg: 11 }).unwrap();
        graph.add_edge(n6, n5, true, EdgeColor { pdg: 11 }).unwrap();
        graph
            .add_edge(n6, n1, false, EdgeColor { pdg: 21 })
            .unwrap();
        graph
            .add_edge(ext4, n4, true, EdgeColor { pdg: 11 })
            .unwrap();
        graph
            .add_edge(n2, ext2, true, EdgeColor { pdg: 11 })
            .unwrap();
        graph
            .add_edge(n3, ext3, true, EdgeColor { pdg: 11 })
            .unwrap();

        let model = load_generic_model("sm");
        let mut he_graph = FeynGenHedgeGraph::from_feyn_gen_symbolica(graph, &model, 2);

        let n_external_fermion_loops = he_graph.number_of_external_fermion_loops();

        assert_eq!(n_external_fermion_loops, 1);

        let mut graph = SymbolicaGraph::new();

        let ext1 = graph.add_node(NodeColorWithoutVertexRule { external_tag: 1 });
        let ext4 = graph.add_node(NodeColorWithoutVertexRule { external_tag: 2 });

        let ext2 = graph.add_node(NodeColorWithoutVertexRule { external_tag: 3 });
        let ext3 = graph.add_node(NodeColorWithoutVertexRule { external_tag: 4 });

        let internal = NodeColorWithoutVertexRule { external_tag: 0 };
        let n1 = graph.add_node(internal);
        let n2 = graph.add_node(internal);
        let n3 = graph.add_node(internal);
        let n4 = graph.add_node(internal);
        let n5 = graph.add_node(internal);
        let n6 = graph.add_node(internal);

        graph
            .add_edge(ext1, n1, true, EdgeColor { pdg: 11 })
            .unwrap();
        graph.add_edge(n1, n2, true, EdgeColor { pdg: 11 }).unwrap();
        graph
            .add_edge(n2, n3, false, EdgeColor { pdg: 21 })
            .unwrap();
        graph.add_edge(n4, n3, true, EdgeColor { pdg: 11 }).unwrap();
        graph
            .add_edge(n4, n5, false, EdgeColor { pdg: 21 })
            .unwrap();
        graph.add_edge(n5, n6, true, EdgeColor { pdg: 11 }).unwrap();
        graph.add_edge(n6, n5, true, EdgeColor { pdg: 11 }).unwrap();
        graph
            .add_edge(n6, n1, false, EdgeColor { pdg: 21 })
            .unwrap();
        graph
            .add_edge(ext4, n4, true, EdgeColor { pdg: 11 })
            .unwrap();
        graph
            .add_edge(n2, ext2, true, EdgeColor { pdg: 11 })
            .unwrap();
        graph
            .add_edge(n3, ext3, true, EdgeColor { pdg: 11 })
            .unwrap();

        let mut he_graph = FeynGenHedgeGraph::from_feyn_gen_symbolica(graph, &model, 2);

        let n_external_fermion_loops = he_graph.number_of_external_fermion_loops();

        assert_eq!(n_external_fermion_loops, 2)
        // SymbolicaGraph<NodeColorWithVertexRule, EdgeColor>,
    }

    #[test]
    fn external_gluon() {
        let mut graph = SymbolicaGraph::new();

        let ext1 = graph.add_node(NodeColorWithoutVertexRule { external_tag: 1 });
        let ext2 = graph.add_node(NodeColorWithoutVertexRule { external_tag: 3 });

        let ext3 = graph.add_node(NodeColorWithoutVertexRule { external_tag: 4 });
        let ext4 = graph.add_node(NodeColorWithoutVertexRule { external_tag: 2 });

        let internal = NodeColorWithoutVertexRule { external_tag: 0 };
        let n1 = graph.add_node(internal);
        let n2 = graph.add_node(internal);
        let n3 = graph.add_node(internal);
        let n4 = graph.add_node(internal);
        let n5 = graph.add_node(internal);
        let n6 = graph.add_node(internal);

        graph
            .add_edge(ext1, n1, true, EdgeColor { pdg: 21 })
            .unwrap();
        graph.add_edge(n1, n2, true, EdgeColor { pdg: 21 }).unwrap();
        graph
            .add_edge(n2, n3, false, EdgeColor { pdg: 21 })
            .unwrap();
        graph.add_edge(n4, n3, true, EdgeColor { pdg: 11 }).unwrap();
        graph
            .add_edge(n4, n5, false, EdgeColor { pdg: 21 })
            .unwrap();
        graph.add_edge(n5, n6, true, EdgeColor { pdg: 11 }).unwrap();
        graph.add_edge(n6, n5, true, EdgeColor { pdg: 11 }).unwrap();
        graph
            .add_edge(n6, n1, false, EdgeColor { pdg: 21 })
            .unwrap();
        graph
            .add_edge(ext4, n4, true, EdgeColor { pdg: 11 })
            .unwrap();
        graph
            .add_edge(ext2, n2, true, EdgeColor { pdg: 21 })
            .unwrap();
        graph
            .add_edge(n3, ext3, true, EdgeColor { pdg: 11 })
            .unwrap();

        let model = load_generic_model("sm");
        let mut he_graph = FeynGenHedgeGraph::from_feyn_gen_symbolica(graph, &model, 2);
        he_graph.number_of_external_fermion_loops();
        // SymbolicaGraph<NodeColorWithVertexRule, EdgeColor>,
    }
}
