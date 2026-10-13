use std::fmt;

use crate::utils::vakint_macros::vk_symbol;
use ahash::{HashMap, HashSet};
use symbolica::{
    atom::{Atom, AtomView},
    graph,
};

use crate::{VakintError, get_individual_momenta_from_atom, get_node_ids, get_prop_with_id};

#[derive(Debug, Clone)]
pub struct Edge {
    pub id: usize,
    pub left_node_id: usize,
    pub right_node_id: usize,
    pub momentum: Atom,
    pub mass: Atom,
}

impl fmt::Display for Edge {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "(#{}|{}->{}|{}|{})",
            self.id, self.left_node_id, self.right_node_id, self.momentum, self.mass
        )
    }
}

#[derive(Debug, Clone)]
pub struct Node {
    pub id: usize,
    pub edges: Vec<(usize, EdgeDirection)>,
}

impl fmt::Display for Node {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "(#{}|[{}])",
            self.id,
            self.edges
                .iter()
                .map(|(e_id, dir)| format!("{}@{}", dir, e_id))
                .collect::<Vec<_>>()
                .join(",")
        )
    }
}

#[derive(Debug, Clone)]
pub enum EdgeDirection {
    Incoming,
    Outgoing,
}

impl EdgeDirection {
    pub fn is_incoming(&self) -> bool {
        match self {
            EdgeDirection::Incoming => true,
            EdgeDirection::Outgoing => false,
        }
    }
    #[allow(unused)]
    pub fn is_outgoing(&self) -> bool {
        match self {
            EdgeDirection::Incoming => false,
            EdgeDirection::Outgoing => true,
        }
    }
}

#[derive(Debug, Clone, Default)]
pub struct Graph {
    pub edges: HashMap<usize, Edge>,
    pub nodes: HashMap<usize, Node>,
}

impl fmt::Display for Graph {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let mut sorted_edges = self.edges.iter().collect::<Vec<_>>();
        sorted_edges.sort_by(|(e1_id, _), (e2_id, _)| e1_id.partial_cmp(e2_id).unwrap());
        let mut sorted_nodes = self.nodes.iter().collect::<Vec<_>>();
        sorted_nodes.sort_by(|(n1_id, _), (n2_id, _)| n1_id.partial_cmp(n2_id).unwrap());
        write!(
            f,
            "Edges: {}\nNodes: {}",
            sorted_edges
                .iter()
                .map(|(_e_id, e)| format!("{}", e))
                .collect::<Vec<_>>()
                .join(" "),
            sorted_nodes
                .iter()
                .map(|(_n_id, n)| format!("{}", n))
                .collect::<Vec<_>>()
                .join(" ")
        )
    }
}

impl fmt::Display for EdgeDirection {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            EdgeDirection::Incoming => write!(f, "IN"),
            EdgeDirection::Outgoing => write!(f, "OUT"),
        }
    }
}

impl Graph {
    pub fn new_from_atom(
        integral_expression: AtomView,
        tot_n_props: usize,
    ) -> Result<Graph, VakintError> {
        let mut graph = Graph::default();
        for i_prop in 1..=tot_n_props {
            if let Some(m) = get_prop_with_id(integral_expression, i_prop) {
                // Format check for the momenta
                let momentum = m.get(&vk_symbol!("q_")).unwrap().to_owned();

                let (left_node_id, right_node_id) = get_node_ids(&m)?;

                graph.edges.insert(
                    i_prop,
                    Edge {
                        id: i_prop,
                        momentum,
                        left_node_id,
                        right_node_id,
                        mass: m.get(&vk_symbol!("mUVsq_")).unwrap().to_owned(),
                    },
                );
            }
        }

        for (&e_id, edge) in graph.edges.iter() {
            if let Some(node) = graph.nodes.get_mut(&edge.left_node_id) {
                node.edges.push((e_id, EdgeDirection::Outgoing));
            } else {
                graph.nodes.insert(
                    edge.left_node_id,
                    Node {
                        id: edge.left_node_id,
                        edges: vec![(e_id, EdgeDirection::Outgoing)],
                    },
                );
            }
            if let Some(node) = graph.nodes.get_mut(&edge.right_node_id) {
                node.edges.push((e_id, EdgeDirection::Incoming));
            } else {
                graph.nodes.insert(
                    edge.right_node_id,
                    Node {
                        id: edge.right_node_id,
                        edges: vec![(e_id, EdgeDirection::Incoming)],
                    },
                );
            }
        }

        for (n_id, nodes) in graph.nodes.iter() {
            if nodes.edges.len() <= 1 {
                return Err(VakintError::MalformedGraph(format!(
                    "Node {} is connected to only {} edges, this cannot be for a vaccuum graph. Graph:\n{}",
                    n_id,
                    nodes.edges.len(),
                    integral_expression
                )));
            }
        }

        Ok(graph)
    }

    pub fn to_dot(&self) -> String {
        format!(
            "digraph G {{\n{}\n}}",
            self.edges
                .values()
                .map(|e| format!(
                    "  {} -> {} [label=\"{}|{}\"]",
                    e.left_node_id, e.right_node_id, e.id, e.momentum
                ))
                .collect::<Vec<_>>()
                .join("\n")
        )
    }

    pub fn to_symbolica_graph(&self, directed: bool) -> graph::Graph<graph::Empty, Atom> {
        let mut g = graph::Graph::new();
        let mut nodes = self.nodes.values().collect::<Vec<_>>();
        nodes.sort_by_key(|node| node.id);
        let symb_nodes: HashMap<usize, usize> =
            HashMap::from_iter(nodes.into_iter().map(|n| (n.id, g.add_node(graph::Empty))));
        let mut edges = self.edges.values().collect::<Vec<_>>();
        edges.sort_by_key(|edge| edge.id);
        for e in edges {
            g.add_edge(
                *symb_nodes.get(&e.right_node_id).unwrap(),
                *symb_nodes.get(&e.left_node_id).unwrap(),
                directed,
                e.mass.to_owned(),
            )
            .unwrap();
        }
        g
    }

    pub fn contract_edges(&mut self, edges_to_contract: &HashSet<usize>) {
        // Merge whole equivalence classes before touching incidence. Independent
        // renames lose unions such as (3,1) and (3,2) after the first rename.
        let mut representatives: HashMap<usize, usize> =
            self.nodes.keys().map(|id| (*id, *id)).collect();
        let mut contractions = edges_to_contract.iter().copied().collect::<Vec<_>>();
        contractions.sort_unstable();
        for id in contractions {
            let edge = &self.edges[&id];
            let left = representatives[&edge.left_node_id];
            let right = representatives[&edge.right_node_id];
            let (keep, remove) = (left.min(right), left.max(right));
            for representative in representatives.values_mut() {
                if *representative == remove {
                    *representative = keep;
                }
            }
        }
        self.edges.retain(|id, _| !edges_to_contract.contains(id));
        self.nodes.retain(|id, _| representatives[id] == *id);
        for node in self.nodes.values_mut() {
            node.edges.clear();
        }
        for edge in self.edges.values_mut() {
            edge.left_node_id = representatives[&edge.left_node_id];
            edge.right_node_id = representatives[&edge.right_node_id];
            self.nodes
                .get_mut(&edge.left_node_id)
                .unwrap()
                .edges
                .push((edge.id, EdgeDirection::Outgoing));
            self.nodes
                .get_mut(&edge.right_node_id)
                .unwrap()
                .edges
                .push((edge.id, EdgeDirection::Incoming));
        }
    }

    pub fn add_one_contraction(
        &self,
        current_contractions: &Vec<HashSet<usize>>,
        n_loops: usize,
        edge_order: &[usize],
    ) -> Vec<HashSet<usize>> {
        let mut res: Vec<HashSet<usize>> = vec![];
        let mut next_possible_contractions = vec![];

        //println!("current_contractions: {:?}", current_contractions);

        // Generate possible unique next contractions
        for c in current_contractions {
            for i_e in edge_order.iter() {
                if c.contains(i_e) {
                    continue;
                }
                let mut new_c = c.clone();
                new_c.insert(*i_e);
                if !next_possible_contractions.contains(&new_c) {
                    next_possible_contractions.push(new_c)
                }
            }
        }
        // println!(
        //     "next_possible_contractions: {:?}",
        //     next_possible_contractions
        // );

        // Filter for unique ones
        let mut canonized_graphs_identified: Vec<graph::Graph<graph::Empty, Atom>> = vec![];
        for c in next_possible_contractions {
            let mut contracted_g = self.clone();
            contracted_g.contract_edges(&c);
            let contracted_g_sym = contracted_g.to_symbolica_graph(false);
            if contracted_g_sym.num_loops() != n_loops {
                continue;
            }
            let canonized_graph = contracted_g_sym.canonize();
            if !canonized_graphs_identified.contains(&canonized_graph.graph) {
                res.push(c);
                canonized_graphs_identified.push(canonized_graph.graph);
            }
        }
        //println!("Resulting next single contractions: {:?}", res);
        res
    }

    pub fn find_unique_contractions(&self) -> Vec<Vec<usize>> {
        let mut unique_contractions: Vec<HashSet<usize>> = vec![];
        let mut next_contractions: Vec<HashSet<usize>> = vec![HashSet::default()];

        // Try and shrink edges with higher IDs and the more complex composite of loop momenta first
        let mut edge_exploration_order = self.edges.keys().cloned().collect::<Vec<_>>();
        edge_exploration_order.sort_by(|e_a, e_b| {
            let key_a = (
                get_individual_momenta_from_atom(self.edges.get(e_a).unwrap().momentum.as_view())
                    .unwrap()
                    .len(),
                *e_a,
            );
            let key_b = (
                get_individual_momenta_from_atom(self.edges.get(e_b).unwrap().momentum.as_view())
                    .unwrap()
                    .len(),
                *e_b,
            );
            key_b.partial_cmp(&key_a).unwrap()
        });

        let master_graph = self.to_symbolica_graph(false);
        loop {
            next_contractions = self.add_one_contraction(
                &next_contractions,
                master_graph.num_loops(),
                &edge_exploration_order,
            );

            if next_contractions.is_empty() {
                break;
            } else {
                unique_contractions.extend(next_contractions.clone());
            }
        }
        let mut res: Vec<Vec<usize>> = unique_contractions
            .iter()
            .map(|c| {
                let mut c_v = c.iter().cloned().collect::<Vec<_>>();
                c_v.sort();
                c_v
            })
            .collect();
        res.sort_by(|a, b| {
            let key_a = (a.len(), a);
            let key_b = (b.len(), b);
            key_a.partial_cmp(&key_b).unwrap()
        });
        res
    }

    /// Independent edges are the complement of a native spanning forest.
    /// Every connected component contributes; self loops are always independent.
    pub fn get_one_lmb(&self) -> Result<Vec<i64>, VakintError> {
        let native = self.to_symbolica_graph(false);
        let mut visited = vec![false; native.nodes().len()];
        let mut tree_edges = HashSet::default();
        for node in 0..visited.len() {
            if visited[node] {
                continue;
            }
            for (index, info) in native.get_spanning_tree(node).nodes.iter().enumerate() {
                visited[index] |= info.position.is_some();
                if let Some(edge) = info.edge_id {
                    tree_edges.insert(edge);
                }
            }
        }
        let mut ids = self.edges.keys().copied().collect::<Vec<_>>();
        ids.sort_unstable();
        ids.into_iter()
            .enumerate()
            .filter(|(index, _)| !tree_edges.contains(index))
            .map(|(_, id)| {
                i64::try_from(id).map_err(|_| {
                    VakintError::MalformedGraph(format!(
                        "edge ID {id} does not fit a loop-basis index"
                    ))
                })
            })
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::utils::vakint_macros::vk_parse;

    #[test]
    fn contractions_join_shared_endpoint_equivalence_classes() {
        crate::Vakint::initialize_vakint_symbols();
        let expression = vk_parse!("topo(prop(1,edge(3,1),k(1),m2,1)*prop(2,edge(3,2),-k(1),m2,1)*prop(3,edge(1,2),k(1),m2,1))").unwrap();
        let mut graph = Graph::new_from_atom(expression.as_view(), 3).unwrap();
        graph.contract_edges(&[1, 2].into_iter().collect());
        assert_eq!(graph.nodes.len(), 1);
        let edge = &graph.edges[&3];
        assert_eq!((edge.left_node_id, edge.right_node_id), (1, 1));
        assert_eq!(graph.nodes[&1].edges.len(), 2);
        assert_eq!(graph.get_one_lmb().unwrap(), vec![3]);
    }

    #[test]
    fn native_spanning_forest_retains_disconnected_self_loops() {
        crate::Vakint::initialize_vakint_symbols();
        let expression = vk_parse!("topo(prop(1,edge(1,1),k(1),m2,1)*prop(2,edge(2,2),k(2),m2,1)*prop(3,edge(3,3),k(3),m2,1)*prop(4,edge(4,4),k(4),m2,1)*prop(5,edge(5,5),k(5),m2,1))").unwrap();
        let graph = Graph::new_from_atom(expression.as_view(), 5).unwrap();
        assert_eq!(graph.get_one_lmb().unwrap(), vec![1, 2, 3, 4, 5]);
    }
}
