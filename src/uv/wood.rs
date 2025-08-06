use std::{cmp::Ordering, collections::VecDeque, hash::Hash, ops::Deref};

use crate::{
    cff::{
        expression::{GraphOrientation, OrientationData, OrientationID},
        generation::{generate_uv_cff, ShiftRewrite},
    },
    model::ArcParticle,
    momentum::Sign,
    new_graph::{Edge, LMBext, LoopMomentumBasis, Vertex},
    numerator::aind::Aind,
    utils::{sign_atom, GS, W_},
};
use ahash::AHashSet;
use bitvec::vec::BitVec;
use eyre::eyre;
use idenso::metric::MS;
use log::debug;
use pathfinding::prelude::BfsReachable;
use serde::{Deserialize, Serialize};
use slotmap::SecondaryMap;
use spenso::{
    algebra::ScalarMul,
    network::parsing::ShadowedStructure,
    shadowing::symbolica_utils::SerializableAtom,
    structure::{
        dimension::Dimension,
        representation::{Minkowski, RepName},
        HasStructure, NamedStructure, OrderedStructure, ToSymbolic,
    },
    tensors::parametric::{
        atomcore::{PatternReplacement, TensorAtomMaps},
        ParamTensor,
    },
};
use symbolica::{
    atom::{Atom, AtomCore, FunctionBuilder, Symbol},
    function,
    id::Replacement,
    parse,
    printer::PrintOptions,
    symbol,
};

use linnet::half_edge::{
    involution::{EdgeIndex, HedgePair, SignOrZero},
    subgraph::{Inclusion, InternalSubGraph, SubGraph, SubGraphOps},
    HedgeGraph,
};

use typed_index_collections::TiVec;
use vakint::{
    vakint_symbol, EvaluationOrder, LoopNormalizationFactor, Vakint, VakintExpression,
    VakintSettings,
};
// use vakint::{EvaluationOrder, LoopNormalizationFactor, Vakint, VakintSettings};

use crate::{
    graph::{BareEdge, BareGraph, BareVertex},
    model::normalise_complex,
};

use super::{
    approx::Approximation,
    poset::{DagNode, PosetNode, DAG},
    Forest, Poset, UltravioletGraph,
};

pub struct Wood {
    poset: Poset<InternalSubGraph, ()>,
    additional_unions: SecondaryMap<PosetNode, Vec<PosetNode>>,
}

impl Wood {
    pub fn n_spinneys(&self) -> usize {
        self.poset.n_nodes()
    }

    pub fn from_spinneys<E, V, H, I: IntoIterator<Item = InternalSubGraph>>(
        s: I,
        graph: impl AsRef<HedgeGraph<E, V, H>>,
    ) -> Self {
        let mut poset = Poset::from_iter(s.into_iter().map(|s| (s, ())));
        let ref_graph = graph.as_ref();

        poset.invert();
        poset.compute_topological_order();

        let mut unions = SecondaryMap::new();

        for (i, sg) in poset.nodes.iter() {
            let cs = ref_graph.connected_components(&sg.data);

            if cs.len() > 1 {
                // sg is a disjoint union of spinneys (at the level of half-edges) (strongly disjoint)
                let mut union = vec![];

                for &c in sg.parents.iter() {
                    let mut is_in = 0;
                    for comp in &cs {
                        let comp =
                            InternalSubGraph::cleaned_filter_optimist(comp.clone(), ref_graph);
                        if comp == poset.nodes[c].data {
                            // find the components in the wood that this union is made of
                            union.push(c);
                            is_in += 1;
                        }
                    }

                    if is_in > 1 {
                        panic!("is in too many components")
                    }
                }

                unions.insert(i, union.clone());
                // sg.children = union;
            }
        }

        // let coverset = poset.to_cover_set();
        Wood {
            poset,
            additional_unions: unions,
        }
    }

    fn unfold_bfs<E, V, H, G>(
        &self,
        graph: &G,
        lmb: &LoopMomentumBasis,
        dag: &mut DAG<Approximation, DagNode, ()>,
        unions: &mut SecondaryMap<PosetNode, Option<Vec<(PosetNode, Option<DagNode>)>>>,
        root: PosetNode,
    ) -> DagNode
    where
        G: UltravioletGraph + AsRef<HedgeGraph<E, V, H>>,
    {
        // let graph = graph.as_ref();
        let mut search_front = VecDeque::new();

        let tree_root = dag.add_node(Approximation::new(
            self.poset.data(root).clone(),
            graph,
            lmb,
        ));
        search_front.push_front((root, tree_root));

        while let Some((node, parent)) = search_front.pop_front() {
            for c in &self.poset.nodes[node].children {
                if let Some(tagged_union) = unions.get_mut(*c) {
                    // Is this node a disjoint union of spinneys
                    if let Some(mut union) = tagged_union.take() {
                        let mut all_supplied = true;
                        for (p, d) in &mut union {
                            if self.poset.nodes[*p].data == self.poset.nodes[node].data {
                                *d = Some(parent);
                            }
                            if d.is_none() {
                                all_supplied = false;
                            }
                        }
                        if all_supplied {
                            let child = dag.add_node(Approximation::new(
                                self.poset.data(*c).clone(),
                                graph,
                                lmb,
                            ));
                            for (_, d) in union {
                                dag.add_edge(d.unwrap(), child);
                            }
                            search_front.push_front((*c, child));
                        } else {
                            *tagged_union = Some(union);
                        }
                    }
                } else {
                    let child =
                        dag.add_node(Approximation::new(self.poset.data(*c).clone(), graph, lmb));
                    dag.add_edge(parent, child);
                    search_front.push_front((*c, child));
                }
            }
        }
        tree_root
    }

    pub fn unfold<E, V, H, G>(&self, graph: &G, lmb: &LoopMomentumBasis) -> Forest
    where
        G: UltravioletGraph + AsRef<HedgeGraph<E, V, H>>,
    {
        let mut dag: DAG<Approximation, DagNode, ()> = DAG::new();

        let root = self.poset.minimum().unwrap();

        let mut unions = SecondaryMap::new();

        for (p, u) in self.additional_unions.iter() {
            let union: Vec<(PosetNode, Option<DagNode>)> = u.iter().map(|i| (*i, None)).collect();
            unions.insert(p, Some(union));
        }

        let _ = self.unfold_bfs(graph, lmb, &mut dag, &mut unions, root);

        Forest { dag }
    }

    pub fn dot(&self, graph: &impl UltravioletGraph) -> String {
        let shift = self.poset.shift();
        self.poset.to_dot_impl(&|n| {
            format!(
                "label={}, dod={},topo_order = {}",
                n.dot_id(shift),
                graph.dod(&n.data),
                // graph.as_ref().count_internal_edges(&n.data),
                n.order.unwrap()
            )
        })
    }

    pub fn dot_spinneys<E, V, H, G>(&self, graph: &G)
    where
        G: UltravioletGraph + AsRef<HedgeGraph<E, V, H>>,
    {
        for s in self.poset.node_values() {
            println!(
                "found {} loop spinney with dod {}:{} ",
                graph.n_loops(s),
                graph.dod(s),
                graph.as_ref().dot(s)
            );
        }
    }
}

impl Wood {
    pub fn show_graphs<H, G>(&self, graph: &G) -> String
    where
        G: UltravioletGraph + AsRef<HedgeGraph<Edge, Vertex, H>>,
    {
        let mut out = String::new();
        out.push_str("Poset structure:\n");
        out.push_str(&self.poset.dot_structure());

        out.push_str("Graphs:\n");
        for (k, n) in self.poset.nodes.iter() {
            out.push_str(&graph.as_ref().dot_impl(
                &n.data,
                format!(
                    "dod={};nodeid ={};\n",
                    graph.dod(&n.data),
                    self.poset.dot_id(k)
                ),
                &|_h| None,
                &|e| Some(format!("dod={}", e.dod)),
                &|n| Some(format!("dod={}", n.dod)),
            ));
            out.push('\n');
        }
        out
    }
}
