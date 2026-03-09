use std::fmt;

use crate::utils::vakint_macros::{vk_parse, vk_symbol};
use ahash::{HashMap, HashSet};
use symbolica::{
    atom::{Atom, AtomCore, AtomView, SliceType},
    graph,
    id::Match,
};

use crate::{
    get_individual_momenta_from_atom, get_node_ids, get_prop_with_id, utils::replace_until_stable,
    VakintError,
};

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
                return Err(VakintError::MalformedGraph(format!("Node {} is connected to only {} edges, this cannot be for a vaccuum graph. Graph:\n{}",
                    n_id, nodes.edges.len(), integral_expression
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
        let symb_nodes: HashMap<usize, usize> = HashMap::from_iter(
            self.nodes
                .values()
                .map(|n| (n.id, g.add_node(graph::Empty))),
        );
        for e in self.edges.values() {
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
        let mut nodes_to_merge = vec![];
        for contracted_edge in edges_to_contract {
            let edge = self.edges.get(contracted_edge).unwrap();
            if edge.left_node_id > edge.right_node_id {
                nodes_to_merge.push((edge.left_node_id, edge.right_node_id));
                self.nodes.remove(&edge.left_node_id);
            } else {
                nodes_to_merge.push((edge.right_node_id, edge.left_node_id));
                self.nodes.remove(&edge.right_node_id);
            }
            self.edges.remove(contracted_edge);
        }

        'relabel_nodes_in_edges: loop {
            let mut did_relabel = false;
            for (old_node_id, new_node_id) in nodes_to_merge.iter() {
                for edge in self.edges.values_mut() {
                    if edge.left_node_id == *old_node_id {
                        edge.left_node_id = *new_node_id;
                        did_relabel = true;
                    }
                    if edge.right_node_id == *old_node_id {
                        edge.right_node_id = *new_node_id;
                        did_relabel = true;
                    }
                }
            }
            if !did_relabel {
                break 'relabel_nodes_in_edges;
            }
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

    // See "counting cycles II" of Ben's blog post at: https://symbolica.io/posts/pattern_matching/
    pub fn get_one_lmb(&self) -> Result<Vec<i64>, VakintError> {
        let mut atom_graph = Atom::num(1);
        let mut node_values = self.nodes.values().collect::<Vec<_>>();
        node_values.sort_by_key(|n| n.id);
        for node in node_values {
            let mut sorted_edges = node.edges.clone();
            sorted_edges.sort_by_key(|(e, _dir)| *e);
            atom_graph *= vk_parse!(&format!(
                "v(n({}),{})",
                node.id,
                sorted_edges
                    .iter()
                    .map(|(e, _dir)| format!("{}", e))
                    .collect::<Vec<_>>()
                    .join(",")
            ))
            .unwrap();
        }
        atom_graph = replace_until_stable(
            atom_graph.as_view(),
            &vk_parse!("v(l1___,x_,r1___)*v(l2___,x_,r2___)")
                .unwrap()
                .to_pattern(),
            &vk_parse!("v(l1___,r1___,l2___,r2___)")
                .unwrap()
                .to_pattern(),
            None,
            None,
        );

        let mut remaining_edges: Vec<i64> = vec![];
        if let Some(m) = atom_graph
            .pattern_match(&vk_parse!("v(args__)").unwrap().to_pattern(), None, None)
            .next_detailed()
        {
            if let Some(Match::Multiple(SliceType::Arg, aviews)) =
                m.match_stack.get(vk_symbol!("args__"))
            {
                for a_id in aviews {
                    let remaining_edge: i64 = match a_id.to_owned().try_into() {
                        Ok(res) => res,
                        Err(_) => continue,
                    };
                    if !remaining_edges.contains(&remaining_edge) {
                        remaining_edges.push(remaining_edge);
                    }
                }
            } else {
                return Err(VakintError::MalformedGraph(format!(
                    "Could not construct an LMB for graph:\n{}",
                    self
                )));
            }
        } else {
            return Err(VakintError::MalformedGraph(format!(
                "Could not construct an LMB for graph:\n{}",
                self
            )));
        }
        Ok(remaining_edges)
    }
}
