use indicatif::ParallelProgressIterator;
use indicatif::ProgressBar;
use indicatif::{ParallelProgressIterator, ProgressStyle};
use log::warn;
use rayon::iter::IntoParallelIterator;
use rayon::iter::ParallelIterator;
use rayon::prelude::*;
use rayon::ThreadPoolBuilder;
use spenso::structure::representation::ExtendibleReps;
use std::collections::hash_map::Entry;
use std::collections::HashSet;
use std::sync::{Arc, Mutex};
use std::time::Instant;
use symbolica::atom::AtomView;
use symbolica::atom::Symbol;
use symbolica::domains::finite_field::PrimeIteratorU64;
use symbolica::fun;
use symbolica::symb;

use ahash::AHashMap;
use ahash::AHashSet;
use ahash::HashMap;
use colored::Colorize;
use log::debug;
use log::info;
use smartstring::{LazyCompact, SmartString};
use symbolica::atom::AtomCore;

use super::NumeratorAwareGraphGroupingOption;
use super::SelfEnergyFilterOptions;
use super::SnailFilterOptions;
use super::TadpolesFilterOptions;
use super::{FeynGenError, FeynGenOptions};

use crate::graph::EdgeType;
use crate::model::ColorStructure;
use crate::model::Particle;
use crate::model::VertexRule;
use crate::momentum::SignOrZero;
use crate::numerator::Color;
use crate::numerator::GlobalPrefactor;
use crate::numerator::Numerator;
use crate::numerator::SymbolicExpression;
use crate::utils;
use crate::{
    feyngen::{FeynGenFilter, GenerationType},
    graph::BareGraph,
    model::Model,
};
use itertools::Itertools;
use linnet::half_edge::subgraph::InternalSubGraph;
use linnet::half_edge::subgraph::OrientedCut;
use linnet::half_edge::HedgeGraph;
use linnet::half_edge::NodeIndex;
use linnet::half_edge::Orientation;
use symbolica::{atom::Atom, graph::Graph as SymbolicaGraph};

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct NodeColorWithVertexRule {
    pub external_tag: i32,
    pub vertex_rule: Arc<VertexRule>,
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Default)]
pub struct NodeColorWithoutVertexRule {
    pub external_tag: i32,
}

#[derive(Debug, Copy, Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Default)]
pub struct EdgeColor {
    pub pdg: isize,
}

impl std::fmt::Display for EdgeColor {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.pdg)
    }
}

impl EdgeColor {
    pub fn from_particle(particle: Arc<Particle>) -> Self {
        Self {
            pdg: particle.pdg_code,
        }
    }
}

pub trait NodeColorFunctions: Sized {
    fn get_external_tag(&self) -> i32;
    fn set_external_tag(&mut self, external_tag: i32);

    fn coupling_orders(&self) -> AHashMap<SmartString<LazyCompact>, usize> {
        AHashMap::default()
    }

    fn get_sign(&self, n_initial_states: usize) -> SignOrZero {
        if self.get_external_tag() > 0 {
            if self.get_external_tag() <= n_initial_states as i32 {
                SignOrZero::Plus
            } else {
                SignOrZero::Minus
            }
        } else {
            SignOrZero::Zero
        }
    }

    fn passes_amplitude_filter(
        _amplitude_subgraph: &InternalSubGraph,
        _graph: &HedgeGraph<Arc<Particle>, Self>,
        _amp_couplings: Option<&std::collections::HashMap<String, usize, ahash::RandomState>>,
    ) -> bool {
        true
    }
}

impl NodeColorFunctions for NodeColorWithVertexRule {
    fn get_external_tag(&self) -> i32 {
        self.external_tag
    }
    fn set_external_tag(&mut self, external_tag: i32) {
        self.external_tag = external_tag;
    }

    fn coupling_orders(&self) -> AHashMap<SmartString<LazyCompact>, usize> {
        let mut coupling_orders = AHashMap::default();
        let vr = self.vertex_rule.clone();
        if vr.name == "external" {
            return coupling_orders;
        }
        for (k, v) in vr.coupling_orders() {
            *coupling_orders.entry(k).or_insert(0) += v;
        }
        coupling_orders
    }

    fn passes_amplitude_filter(
        amplitude_subgraph: &InternalSubGraph,
        graph: &HedgeGraph<Arc<Particle>, Self>,
        amp_couplings: Option<&std::collections::HashMap<String, usize, ahash::RandomState>>,
    ) -> bool {
        if let Some(amp_couplings) = amp_couplings {
            let mut coupling_orders = AHashMap::default();
            for (_, s) in graph.iter_node_data(amplitude_subgraph) {
                if s.get_external_tag() != 0 {
                    for (k, v) in s.coupling_orders() {
                        *coupling_orders.entry(k).or_insert(0) += v;
                    }
                }
            }

            amp_couplings.iter().all(|(k, v)| {
                coupling_orders
                    .get(&SmartString::from(k))
                    .map_or(0 == *v, |o| *o == *v)
            })
        } else {
            true
        }
    }
}

impl NodeColorFunctions for NodeColorWithoutVertexRule {
    fn get_external_tag(&self) -> i32 {
        self.external_tag
    }
    fn set_external_tag(&mut self, external_tag: i32) {
        self.external_tag = external_tag;
    }
}

impl std::fmt::Display for NodeColorWithoutVertexRule {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.external_tag)
    }
}

impl std::fmt::Display for NodeColorWithVertexRule {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "({}|{})",
            self.external_tag,
            self.vertex_rule.name.clone()
        )
    }
}

#[derive(Debug, Clone)]
pub struct FeynGen {
    pub options: FeynGenOptions,
}

impl FeynGen {
    pub fn new(options: FeynGenOptions) -> Self {
        Self { options }
    }

    #[allow(clippy::type_complexity)]
    pub fn assign_node_colors(
        model: &Model,
        graph: &SymbolicaGraph<NodeColorWithoutVertexRule, EdgeColor>,
        node_colors: &HashMap<
            Vec<(Option<bool>, SmartString<LazyCompact>)>,
            Vec<SmartString<LazyCompact>>,
        >,
    ) -> Result<Vec<(SymbolicaGraph<NodeColorWithVertexRule, EdgeColor>, usize)>, FeynGenError>
    {
        // println!("graph = {}", graph.to_dot());
        let mut colored_nodes: Vec<Vec<NodeColorWithVertexRule>> = vec![vec![]];
        let edges = graph.edges();
        for (i_n, node) in graph.nodes().iter().enumerate() {
            let mut node_edges = vec![];
            for e in node.edges.iter() {
                let orientation: Option<bool> = if edges[*e].directed {
                    Some(edges[*e].vertices.0 == i_n)
                } else {
                    None
                };
                let particle = model.get_particle_from_pdg(edges[*e].data.pdg);
                node_edges.push((orientation, particle.name.clone()));
                // self-loop edge must be counted twice
                if edges[*e].vertices.0 == edges[*e].vertices.1 {
                    node_edges.push((orientation.map(|o| !o), particle.name.clone()));
                }
            }
            node_edges.sort();
            let dummy_external_vertex_rule = Arc::new(VertexRule {
                name: "external".into(),
                couplings: vec![],
                lorentz_structures: vec![],
                particles: vec![],
                color_structures: ColorStructure::new(vec![]),
            });
            let colors = if node_edges.len() == 1 {
                &vec![dummy_external_vertex_rule]
            } else if let Some(cs) = node_colors.get(&node_edges) {
                &cs.iter()
                    .map(|c| model.get_vertex_rule(c))
                    .collect::<Vec<_>>()
            } else {
                return Err(FeynGenError::GenericError(format!(
                    "Could not find node colors for node edges {:?}",
                    node_edges
                )));
            };
            let mut new_colored_nodes: Vec<Vec<NodeColorWithVertexRule>> = vec![];
            for current_colors in colored_nodes.iter_mut() {
                for color in colors {
                    current_colors.push(NodeColorWithVertexRule {
                        external_tag: node.data.external_tag,
                        vertex_rule: color.clone(),
                    });
                    new_colored_nodes.push(current_colors.clone());
                }
            }

            colored_nodes = new_colored_nodes;
        }

        let mut colored_graphs = vec![];
        for color in colored_nodes.clone() {
            let mut g: SymbolicaGraph<_, _> = SymbolicaGraph::new();
            for colored_node in color {
                g.add_node(colored_node);
            }
            for edge in graph.edges() {
                _ = g.add_edge(
                    edge.vertices.0,
                    edge.vertices.1,
                    edge.directed,
                    edge.data.clone(),
                );
            }
            colored_graphs.push(g);
        }
        Ok(FeynGen::group_isomorphic_graphs(&colored_graphs))
    }

    #[allow(clippy::type_complexity)]
    pub fn group_isomorphic_graphs<
        NodeColor: Clone + PartialEq + Eq + PartialOrd + Ord + std::hash::Hash,
    >(
        graphs: &[SymbolicaGraph<NodeColor, EdgeColor>],
    ) -> Vec<(SymbolicaGraph<NodeColor, EdgeColor>, usize)> {
        if graphs.len() == 1 {
            return vec![(graphs[0].clone(), 1)];
        }
        let mut iso_buckets: HashMap<SymbolicaGraph<NodeColor, EdgeColor>, usize> =
            HashMap::default();
        for g in graphs.iter() {
            let canonized_g = g.canonize();
            *iso_buckets.entry(canonized_g.graph).or_insert_with(|| 1) += 1;
        }
        let mapped_graph = iso_buckets
            .iter()
            .map(|(g, count)| (g.clone(), *count))
            .collect();

        mapped_graph
    }

    pub fn contains_particles(
        graph: &SymbolicaGraph<NodeColorWithoutVertexRule, EdgeColor>,
        particles: &[isize],
    ) -> bool {
        let mut particles_stack = particles.to_vec();
        for edge in graph.edges().iter() {
            if let Some(pos) = particles_stack.iter().position(|&x| x == edge.data.pdg) {
                particles_stack.swap_remove(pos);
            }
            if particles_stack.is_empty() {
                return true;
            }
        }
        false
    }

    pub fn find_edge_position(
        graph: &SymbolicaGraph<NodeColorWithoutVertexRule, EdgeColor>,
        edge_vertices: (usize, usize),
        oriented: bool,
    ) -> Option<usize> {
        for (i_e, edge) in graph.edges().iter().enumerate() {
            if edge.vertices == edge_vertices
                || (!oriented && edge.vertices == (edge_vertices.1, edge_vertices.0))
            {
                return Some(i_e);
            }
        }
        None
    }

    pub fn veto_special_topologies(
        model: &Model,
        graph: &SymbolicaGraph<NodeColorWithoutVertexRule, EdgeColor>,
        veto_self_energy: Option<&SelfEnergyFilterOptions>,
        veto_tadpole: Option<&TadpolesFilterOptions>,
        veto_snails: Option<&SnailFilterOptions>,
    ) -> bool {
        if graph.nodes().iter().any(|n| n.data.external_tag < 0) {
            panic!("External tag must be positive, but found negative as obtained when performing external state symmetrization");
        }
        let debug = false;
        if debug {
            debug!(
                "\n\n>> Vetoing special topologies for the following {}-loop graph:\n{}",
                graph.num_loops(),
                graph.to_dot()
            );
            debug!(
                "Graph edges raw representation:\n{}",
                graph
                    .edges()
                    .iter()
                    .enumerate()
                    .map(|(i_e, e)| format!(
                        "[ #{} vertices: {:?}, e.directed: {}, e.data: {:?} ]",
                        i_e, e.vertices, e.directed, e.data
                    ))
                    .join("\n")
            );
            debug!(
                "Graph nodes raw representation:\n{}",
                graph
                    .nodes()
                    .iter()
                    .enumerate()
                    .map(|(i_n, n)| format!(
                        "[ #{} edges: {:?}, n.data: {:?} ]",
                        i_n, n.edges, n.data
                    ))
                    .join("\n")
            );
        }
        if veto_self_energy.is_none() && veto_tadpole.is_none() && veto_snails.is_none() {
            return false;
        }
        let graph_nodes: &[symbolica::graph::Node<NodeColorWithoutVertexRule>] = graph.nodes();
        let graph_edges: &[symbolica::graph::Edge<EdgeColor>] = graph.edges();

        let max_external = graph
            .nodes()
            .iter()
            .filter(|n| n.data.external_tag > 0)
            .map(|n| n.data.external_tag)
            .max()
            .unwrap_or(0) as usize;
        if max_external == 0 {
            // Do not implement any veto for vacuum graphs
            return false;
        }
        let max_external_node_position = graph
            .nodes()
            .iter()
            .position(|n| n.data.external_tag == (max_external as i32))
            .unwrap();
        if debug {
            debug!(
                "Spanning tree root position: max_external={},max_external_node_position={}",
                max_external, max_external_node_position
            );
        }

        let mut external_partices: Vec<Arc<Particle>> =
            vec![model.particles[0].clone(); max_external];
        for e in graph_edges {
            if graph_nodes[e.vertices.0].data.external_tag != 0 {
                external_partices[(graph_nodes[e.vertices.0].data.external_tag - 1) as usize] =
                    model.get_particle_from_pdg(e.data.pdg);
            } else if graph_nodes[e.vertices.1].data.external_tag != 0 {
                external_partices[(graph_nodes[e.vertices.1].data.external_tag - 1) as usize] =
                    model.get_particle_from_pdg(e.data.pdg);
            }
        }

        let mut spanning_tree = graph.get_spanning_tree(max_external_node_position);
        spanning_tree.chain_decomposition();

        if debug {
            debug!(
                "Spanning tree: order={:?}\n{}",
                spanning_tree.order,
                spanning_tree
                    .nodes
                    .iter()
                    .enumerate()
                    .map(|(i_n, n)| format!(
                        "[ #{} position: {:?}, parent: {:?}, back_edges: {:?}, external: {}, chain_id: {:?} ]",
                        i_n, n.position, n.parent, n.back_edges, n.external, n.chain_id
                    ))
                    .join("\n")
            );
        }
        let mut node_children: Vec<Vec<usize>> = vec![vec![]; spanning_tree.nodes.len()];
        for (i_n, node) in spanning_tree.nodes.iter().enumerate() {
            node_children[node.parent].push(i_n);
        }
        if debug {
            debug!("node_children={:?}", node_children);
        }
        let mut external_momenta_routing: Vec<Vec<usize>> = vec![vec![]; spanning_tree.nodes.len()];
        for (i_n, node) in graph.nodes().iter().enumerate() {
            if (node.edges.len() != 1)
                || node.data.external_tag == (max_external_node_position as i32)
            {
                continue;
            }
            let external_index = node.data.external_tag as usize;
            external_momenta_routing[i_n].push(external_index);
            let mut next_node = spanning_tree.nodes[i_n].parent;
            // println!("max_external={},max_external_node_position={},sink_node_position_in_spanning_tree={}", max_external, max_external_node_position, sink_node_position_in_spanning_tree);
            // println!("external_momenta_routing={:?}", external_momenta_routing);
            // println!("spanning_tree={:?}", spanning_tree);
            // println!("Starting from node #{} = {:?}", i_n, node);
            while next_node != max_external_node_position {
                external_momenta_routing[next_node].push(external_index);
                next_node = spanning_tree.nodes[next_node].parent;
            }
        }
        for route in external_momenta_routing.iter_mut() {
            if route.len() == max_external - 1 {
                *route = vec![max_external];
            }
        }
        if debug {
            debug!("external_momenta_routing={:?}", external_momenta_routing);
        }

        // See https://arxiv.org/pdf/1209.0700 for information on the logic of this algorithm

        // Tuple format: (external_leg_id, back_edge_start_node_position, back_edge_position_in_list, chain_id)
        let mut self_energy_attachments: HashSet<(usize, usize, usize, usize)> = HashSet::default();
        // Tuple format: (back_edge_start_node_position, back_edge_position_in_list, chain_id)
        let mut vacuum_attachments: HashSet<(usize, usize, usize)> = HashSet::default();
        // Tuple format: (back_edge_start_node_position, back_edge_position_in_list, chain_id)
        let mut self_loops: HashSet<(usize, usize, usize)> = HashSet::default();
        for &i_n in &spanning_tree.order {
            let node = &spanning_tree.nodes[i_n];
            for (i_back_edge, &back_edge) in node.back_edges.iter().enumerate() {
                let i_chain = i_n;
                if back_edge == i_n {
                    self_loops.insert((i_n, i_back_edge, i_chain));
                    continue;
                }
                let mut self_energy_external_leg_id: Option<usize> = None;
                let mut curr_chain_node = back_edge;
                let mut is_valid_chain = true;
                'follow_chain: loop {
                    if curr_chain_node == i_n {
                        break 'follow_chain;
                    }
                    let moms = &external_momenta_routing[curr_chain_node];
                    if moms.len() == 1 {
                        if let Some(se_leg) = self_energy_external_leg_id {
                            if se_leg != moms[0] {
                                is_valid_chain = false;
                                break 'follow_chain;
                            }
                        }
                        self_energy_external_leg_id = Some(moms[0]);
                        for child in node_children[curr_chain_node].iter() {
                            if (!external_momenta_routing[*child].is_empty())
                                && external_momenta_routing[*child][0] != moms[0]
                            {
                                is_valid_chain = false;
                                break 'follow_chain;
                            }
                        }
                    } else if !moms.is_empty() {
                        is_valid_chain = false;
                        break 'follow_chain;
                    }
                    if let Some(chain_id) = spanning_tree.nodes[curr_chain_node].chain_id {
                        if chain_id != i_chain {
                            is_valid_chain = false;
                            break 'follow_chain;
                        }
                    } else {
                        is_valid_chain = false;
                        break 'follow_chain;
                    }
                    curr_chain_node = spanning_tree.nodes[curr_chain_node].parent;
                }
                if is_valid_chain {
                    if let Some(leg_id) = self_energy_external_leg_id {
                        // Make sure the attachment point of the self-energy does not receive any other external momenta
                        if external_momenta_routing[i_n].len() == 1
                            && external_momenta_routing[i_n][0] == leg_id
                        {
                            // Also For 1 -> 1 processes, we must also verify that it is not the whole graph
                            // Also For 1 -> 1 processes, we must also verify that it is not the whole graph
                            if max_external > 2
                                || spanning_tree
                                    .nodes
                                    .iter()
                                    .filter(|n| (!n.external) && n.chain_id.is_none())
                                    .count()
                                    > 1
                            {
                                self_energy_attachments.insert((leg_id, i_n, i_back_edge, i_chain));
                            }
                        }
                    } else {
                        vacuum_attachments.insert((i_n, i_back_edge, i_chain));
                    }
                }
            }
        }
        if debug {
            debug!("self_energy_attachments={:?}", self_energy_attachments);
            debug!("vacuum_attachments={:?}", vacuum_attachments);
            debug!("self_loops={:?}", self_loops);
        }
        let mut tree_bridge_node_indices: HashSet<usize> = HashSet::default();
        for (i_n, node) in spanning_tree.nodes.iter().enumerate() {
            if node.chain_id.is_none()
                && !node.external
                && !external_momenta_routing[i_n].is_empty()
                && !spanning_tree.nodes[node.parent]
                    .back_edges
                    .iter()
                    .any(|&end| i_n == end)
            {
                tree_bridge_node_indices.insert(i_n);
            }
        }
        if debug {
            debug!("bridge_node_positions={:?}", tree_bridge_node_indices);
        }

        // For self-energies we must confirm that they are self-energies by checking if the back edge start node is a bridge
        for (leg_id, back_edge_start_node_index, _back_edge_position_in_list, _chain_id) in
            self_energy_attachments.iter()
        {
            if tree_bridge_node_indices.contains(back_edge_start_node_index) {
                if let Some(veto_self_energy_options) = veto_self_energy {
                    if debug {
                        debug!(
                            "Vetoing self-energy for leg_id={}, back_edge_start_node_index={}, back_edge_position_in_list={}, chain_id={}, with options:\n{:?}",
                            leg_id, back_edge_start_node_index, _back_edge_position_in_list, _chain_id, veto_self_energy_options
                        );
                    }
                    if veto_self_energy_options.veto_only_scaleless_self_energy {
                        panic!(
                            "Option to only remove scaleless self-energies is not implemented yet"
                        );
                    } else {
                        #[allow(clippy::unnecessary_unwrap)]
                        if external_partices[leg_id - 1].is_massive() {
                            if veto_self_energy_options.veto_self_energy_of_massive_lines {
                                return true;
                            }
                        } else {
                            #[allow(clippy::unnecessary_unwrap)]
                            if veto_self_energy_options.veto_self_energy_of_massless_lines {
                                return true;
                            }
                        }
                    }
                }
            }
        }
        // For vaccuum attachments, we must differentiate if there are snails (start node is a tree bridge) or tadpoles (start node *is not* a tree bridge)
        for (back_edge_start_node_index, _back_edge_position_in_list, _chain_id) in
            vacuum_attachments.iter().chain(self_loops.iter())
        {
            let mut first_tree_attachment_node_index = *back_edge_start_node_index;
            while external_momenta_routing[first_tree_attachment_node_index].is_empty() {
                first_tree_attachment_node_index =
                    spanning_tree.nodes[first_tree_attachment_node_index].parent;
            }
            let attachment_edge = &graph_edges[FeynGen::find_edge_position(
                graph,
                (
                    first_tree_attachment_node_index,
                    spanning_tree.nodes[first_tree_attachment_node_index].parent,
                ),
                false,
            )
            .unwrap_or_else(|| {
                panic!("Could not find edge between bridge node parent and grandparent")
            })];
            let attachment_particle = model.get_particle_from_pdg(attachment_edge.data.pdg);

            if !tree_bridge_node_indices.contains(back_edge_start_node_index) {
                // Tadpole
                if let Some(veto_tadpole_options) = veto_tadpole {
                    if debug {
                        debug!(
                            "Vetoing tadpole for back_edge_start_node_index={}, back_edge_position_in_list={}, chain_id={}, with options:\n{:?}",
                            back_edge_start_node_index, _back_edge_position_in_list, _chain_id, veto_tadpole_options
                        );
                    }
                    if veto_tadpole_options.veto_only_scaleless_tadpoles {
                        panic!(
                            "Option to only remove scaleless self-energies is not implemented yet"
                        );
                    } else {
                        #[allow(clippy::unnecessary_unwrap)]
                        if attachment_particle.is_massive() {
                            if veto_tadpole_options.veto_tadpoles_attached_to_massive_lines {
                                return true;
                            }
                        } else {
                            #[allow(clippy::unnecessary_unwrap)]
                            if veto_tadpole_options.veto_tadpoles_attached_to_massless_lines {
                                return true;
                            }
                        }
                    }
                }
            } else {
                #[allow(clippy::unnecessary_unwrap)]
                // Snail
                if let Some(veto_snails_options) = veto_snails {
                    if debug {
                        debug!(
                            "Vetoing snail for back_edge_start_node_index={}, back_edge_position_in_list={}, chain_id={}, with options:\n{:?}",
                            back_edge_start_node_index, _back_edge_position_in_list, _chain_id, veto_snails_options
                        );
                    }
                    #[allow(clippy::unnecessary_unwrap)]
                    if veto_snails_options.veto_only_scaleless_snails {
                        panic!(
                            "Option to only remove scaleless self-energies is not implemented yet"
                        );
                    } else {
                        #[allow(clippy::unnecessary_unwrap)]
                        if attachment_particle.is_massive() {
                            if veto_snails_options.veto_snails_attached_to_massive_lines {
                                return true;
                            }
                        } else {
                            #[allow(clippy::unnecessary_unwrap)]
                            if veto_snails_options.veto_snails_attached_to_massless_lines {
                                return true;
                            }
                        }
                    }
                }
            }
        }

        if debug {
            debug!(">> No special topology veto applied to this graph");
        }
        false
    }

    pub(super) fn unresolved_cut_content(&self, model: &Model) -> (usize, AHashSet<Arc<Particle>>) {
        if let Some(p) = self.options.cross_section_filters.get_perturbative_orders() {
            let mut unresolved = AHashSet::new();
            for k in p.keys() {
                let k: SmartString<_> = k.clone().into();
                if let Some(p) = model.unresolved_particles.get(&k) {
                    unresolved = unresolved.union(p).cloned().collect();
                }
            }
            (p.values().sum(), unresolved)
        } else {
            (0, AHashSet::new())
        }
    }

    pub fn contains_cut<NodeColor>(
        &self,
        model: &Model,
        graph: &SymbolicaGraph<NodeColor, EdgeColor>,
        n_unresolved: usize,
        unresolved_type: &AHashSet<Arc<Particle>>,
    ) -> bool
    where
        NodeColor: NodeColorFunctions + Clone,
    {
        #[allow(clippy::too_many_arguments)]
        fn is_valid_cut<NodeColor: NodeColorFunctions>(
            cut: &(InternalSubGraph, OrientedCut, InternalSubGraph),
            s_set: &AHashSet<NodeIndex>,
            model: &Model,
            n_unresolved: usize,
            unresolved_type: &AHashSet<Arc<Particle>>,
            particle_content: &[Arc<Particle>],
            amp_couplings: Option<&std::collections::HashMap<String, usize, ahash::RandomState>>,
            amp_loop_count: Option<(usize, usize)>,
            graph: &HedgeGraph<Arc<Particle>, NodeColor>,
        ) -> bool {
            if is_s_channel(cut, s_set, graph) {
                let mut cut_content: Vec<_> = cut
                    .1
                    .iter_edges_relative(graph)
                    .map(|(o, d)| {
                        if matches!(o, Orientation::Reversed) {
                            d.data.as_ref().unwrap().get_anti_particle(model)
                        } else {
                            d.data.as_ref().unwrap().clone()
                        }
                    })
                    .collect();

                for p in particle_content.iter() {
                    if let Some(pos) = cut_content.iter().position(|c| c == p) {
                        cut_content.swap_remove(pos);
                    } else {
                        return false;
                    }
                }

                if cut_content.len() > n_unresolved {
                    return false;
                }

                for p in cut_content.iter() {
                    if !unresolved_type.contains(p) {
                        return false;
                    }
                }

                if let Some((min_loop, max_loop)) = amp_loop_count {
                    let left_loop = graph.cyclotomatic_number(&cut.0);
                    let right_loop = graph.cyclotomatic_number(&cut.2);

                    let sum = left_loop + right_loop;
                    if sum < min_loop || sum > max_loop {
                        return false;
                    }
                }
                NodeColor::passes_amplitude_filter(&cut.0, graph, amp_couplings)
                    && NodeColor::passes_amplitude_filter(&cut.2, graph, amp_couplings)
            } else {
                false
            }
        }

        fn is_s_channel<NodeColor>(
            cut: &(InternalSubGraph, OrientedCut, InternalSubGraph),
            s_set: &AHashSet<NodeIndex>,
            graph: &HedgeGraph<Arc<Particle>, NodeColor>,
        ) -> bool {
            let nodes: AHashSet<_> = graph
                .iter_node_data(&cut.0)
                .map(|(i, _)| graph.id_from_hairs(i).unwrap())
                .collect();

            s_set.is_subset(&nodes)
        }

        let particle_content = self
            .options
            .final_pdgs
            .iter()
            .map(|&pdg| model.get_particle_from_pdg(pdg as isize))
            .collect::<Vec<_>>();
        let amp_couplings = self.options.amplitude_filters.get_coupling_orders();
        let amp_loop_count = self.options.amplitude_filters.get_loop_count_range();
        let n_particles = self.options.initial_pdgs.len();
        let he_graph = HedgeGraph::from(graph.clone()).map(
            |node_color| node_color,
            |_, d| d.map(|d| model.get_particle_from_pdg(d.pdg)),
        );

        let mut s_set = AHashSet::new();
        let mut t_set = vec![];

        for (n, f) in he_graph.iter_nodes() {
            let id = he_graph.id_from_hairs(n).unwrap();
            match f.get_sign(n_particles) {
                SignOrZero::Plus => {
                    s_set.insert(id);
                }
                SignOrZero::Minus => {
                    t_set.push(id);
                }
                _ => {}
            }
        }
        if let (Some(&s), Some(&t)) = (s_set.iter().next(), t_set.first()) {
            let cuts = he_graph.all_cuts(s, t);

            let pass_cut_filter = cuts.iter().any(|c| {
                is_valid_cut(
                    c,
                    &s_set,
                    model,
                    n_unresolved,
                    unresolved_type,
                    &particle_content,
                    amp_couplings,
                    amp_loop_count,
                    &he_graph,
                )
            });

            pass_cut_filter
        } else {
            true //TODO still check the amplitude level filters in the case where there is no initial-state specified
        }
    }

    // This fast cut checker does not enumerate all cuts, but rather checks if the graph can contain a cut with the right particles
    // It also does not consider amplitude-level filters
    pub fn contains_cut_fast<NodeColor: NodeColorFunctions>(
        &self,
        model: &Model,
        graph: &SymbolicaGraph<NodeColor, EdgeColor>,
        n_unresolved: usize,
        unresolved_type: &AHashSet<Arc<Particle>>,
        particles: &[isize],
    ) -> bool {
        let n_initial_states = self.options.initial_pdgs.len();
        // Quick filter to check if the diagram can possibly contain a cut with the right particles and splitting the graph in two parts,
        // as expected for a final-state cut of a forward scattering graph.

        // TODO: replace this with s_and_t_cut from hedge, so as to have proper handling of particle vs anti-particles.
        //let he_graph = HedgeGraph::from(graph.clone());
        //he_graph.all_s_t_cuts(s, t, regions)

        // for now normalize it all to particles
        let mut cut_map: std::collections::HashMap<
            isize,
            std::vec::Vec<usize>,
            ahash::RandomState,
        > = HashMap::default();
        for &p in particles.iter().collect::<HashSet<_>>() {
            let particle = model.get_particle_from_pdg(p);
            if particle.is_antiparticle() {
                cut_map.insert(particle.get_anti_particle(model).pdg_code, vec![]);
            } else {
                cut_map.insert(p, vec![]);
            }
        }
        let mut unresolved_set: Vec<usize> = vec![];
        for (i_e, edge) in graph.edges().iter().enumerate() {
            // filter out external edges
            if graph.nodes()[edge.vertices.0].data.get_external_tag() != 0
                || graph.nodes()[edge.vertices.1].data.get_external_tag() != 0
            {
                continue;
            }
            let particle = model.get_particle_from_pdg(edge.data.pdg);
            let e = if particle.is_antiparticle() {
                cut_map.get_mut(&particle.get_anti_particle(model).pdg_code)
            } else {
                cut_map.get_mut(&particle.pdg_code)
            };
            if let Some(cut_entry) = e {
                cut_entry.push(i_e);
            }
            if unresolved_type.contains(&particle) {
                unresolved_set.push(i_e);
            }
        }
        // Generate unique, unordered combinations, for each possible unresolved multiplicity
        let mut unique_combinations: Vec<
            std::collections::HashSet<Vec<usize>, ahash::RandomState>,
        > = vec![];
        for unresolved_multiplicity in 0..=n_unresolved {
            let mut unique_combinations_for_this_multiplicity: std::collections::HashSet<
                Vec<usize>,
                ahash::RandomState,
            > = HashSet::default();
            let mut particles_for_this_multiplicity =
                particles.iter().map(Some).collect::<Vec<_>>();
            for _ in 0..unresolved_multiplicity {
                particles_for_this_multiplicity.push(None);
            }
            for combination in particles_for_this_multiplicity
                .iter()
                .map(|p| {
                    if let Some(particle) = p {
                        cut_map.get(particle).unwrap().iter()
                    } else {
                        unresolved_set.iter()
                    }
                })
                .multi_cartesian_product()
            {
                // Sort the combination and insert into HashSet
                let mut sorted_combination: Vec<_> = combination.into_iter().cloned().collect();
                if sorted_combination.iter().collect::<HashSet<_>>().len()
                    != sorted_combination.len()
                {
                    continue;
                }
                sorted_combination.sort_unstable();
                unique_combinations_for_this_multiplicity.insert(sorted_combination);
            }
            unique_combinations.push(unique_combinations_for_this_multiplicity);
        }

        let mut adj_list: HashMap<usize, Vec<(usize, usize)>> = HashMap::default();
        for (i_e, e) in graph.edges().iter().enumerate() {
            adj_list
                .entry(e.vertices.0)
                .or_default()
                .push((i_e, e.vertices.1));
            adj_list
                .entry(e.vertices.1)
                .or_default()
                .push((i_e, e.vertices.0));
        }

        let left_initial_state_positions = (1..=n_initial_states)
            .map(|leg_id| {
                graph
                    .nodes()
                    .iter()
                    .position(|n| n.data.get_external_tag() == (leg_id as i32))
                    .unwrap()
            })
            .collect::<Vec<_>>();
        let right_initial_state_positions = (n_initial_states + 1..=2 * n_initial_states)
            .map(|leg_id| {
                graph
                    .nodes()
                    .iter()
                    .position(|n| n.data.get_external_tag() == (leg_id as i32))
                    .unwrap()
            })
            .collect::<Vec<_>>();

        fn are_incoming_connected_to_outgoing<NodeColor>(
            cut: &[usize],
            graph: &SymbolicaGraph<NodeColor, EdgeColor>,
            adj_list: &HashMap<usize, Vec<(usize, usize)>>,
            left_nodes: &[usize],
            right_nodes: &[usize],
        ) -> bool {
            for left_node in left_nodes.iter() {
                for right_node in right_nodes.iter() {
                    let mut visited: Vec<bool> = vec![false; graph.nodes().len()];
                    let mut stack: Vec<usize> = vec![*left_node];
                    while let Some(node) = stack.pop() {
                        if node == *right_node {
                            return true;
                        }
                        visited[node] = true;
                        for (i_e, neighbor) in adj_list[&node].iter() {
                            if !cut.contains(i_e) && !visited[*neighbor] {
                                stack.push(*neighbor);
                            }
                        }
                    }
                }
            }
            false
        }
        for unique_combination_with_fixed_multiplicity in unique_combinations {
            'cutloop: for cut in unique_combination_with_fixed_multiplicity {
                if are_incoming_connected_to_outgoing(
                    &cut,
                    graph,
                    &adj_list,
                    &left_initial_state_positions,
                    &right_initial_state_positions,
                ) {
                    continue 'cutloop;
                }
                // Make sure the cut is minimal and that no subcut is a valid cut.
                for subcut in (1..=(cut.len() - 1))
                    .flat_map(|offset| cut.iter().combinations(cut.len() - offset))
                {
                    if !are_incoming_connected_to_outgoing(
                        subcut.iter().map(|&id| *id).collect::<Vec<_>>().as_slice(),
                        graph,
                        &adj_list,
                        &left_initial_state_positions,
                        &right_initial_state_positions,
                    ) {
                        continue 'cutloop;
                    }
                }
                return true;
            }
        }

        false
    }

    fn group_isomorphic_graphs_after_node_color_change<NC>(
        graphs: &HashMap<SymbolicaGraph<NC, EdgeColor>, Atom>,
        node_colors_to_change: &HashMap<i32, i32>,
        pool: &rayon::ThreadPool,
        progress_bar_style: &ProgressStyle,
    ) -> HashMap<SymbolicaGraph<NC, EdgeColor>, Atom>
    where
        NC: NodeColorFunctions + Clone + PartialOrd + Ord + Eq + std::hash::Hash + Send + Sync,
    {
        #[allow(clippy::type_complexity)]
        let iso_buckets: Arc<
            Mutex<
                HashMap<
                    SymbolicaGraph<NC, EdgeColor>,
                    (usize, (SymbolicaGraph<NC, EdgeColor>, Vec<usize>, Atom)),
                >,
            >,
        > = Arc::new(Mutex::new(HashMap::default()));

        let iso_buckets_clone = iso_buckets.clone();
        let bar = ProgressBar::new(graphs.len() as u64);
        bar.set_style(progress_bar_style.clone());
        bar.set_message("Grouping isomorphic topologies including external leg symmetrization...");
        pool.install(|| {
            graphs
                .par_iter()
                .progress_with(bar.clone())
                .for_each(|(g, symmetry_factor)| {
                    let mut g_node_color_modified = g.clone();
                    let mut modifications: Vec<(usize, i32)> = vec![];
                    for (i_n, node) in g.nodes().iter().enumerate() {
                        for (src_node_color, trgt_node_color) in node_colors_to_change {
                            if node.data.get_external_tag() == *src_node_color {
                                modifications.push((i_n, *trgt_node_color));
                            }
                        }
                    }
                    for (i_n, new_color) in modifications {
                        let mut node_data = g.nodes()[i_n].data.clone();
                        node_data.set_external_tag(new_color);
                        g_node_color_modified.set_node_data(i_n, node_data);
                    }
                    let canonized_g = g_node_color_modified.canonize();
                    let ext_ordering = g
                        .nodes()
                        .iter()
                        .filter_map(|n| {
                            if n.data.get_external_tag() > 0 {
                                Some(*n.edges.first().unwrap())
                            } else {
                                None
                            }
                        })
                        .collect::<Vec<_>>();

                    {
                        let mut iso_buckets_guard = iso_buckets_clone.lock().unwrap();

                        if let Some((count, (other_g, external_ordering, sym_fact))) =
                            iso_buckets_guard.get_mut(&canonized_g.graph)
                        {
                            // Make sure to pick a representative with a canonical external ordering so that additional flavour grouping afterwards is possible
                            if ext_ordering < *external_ordering {
                                *external_ordering = ext_ordering;
                                *other_g = g.clone();
                            }
                            *count += 1;

                            assert!(sym_fact == symmetry_factor);
                        } else {
                            iso_buckets_guard.insert(
                                canonized_g.graph,
                                (1, (g.clone(), ext_ordering, symmetry_factor.clone())),
                            );
                        }
                    }
                });
        });
        bar.finish_and_clear();

        let iso_buckets_guard = iso_buckets.lock().unwrap();
        iso_buckets_guard
            .iter()
            .map(
                |(_canonized_g, (count, (g, _external_ordering, symmetry_factor)))| {
                    (g.clone(), symmetry_factor * Atom::new_num(*count as i64))
                },
            )
            .collect()
    }

    pub fn canonize_external_momenta_assignment(
        &self,
        model: &Model,
        graph: &mut SymbolicaGraph<NodeColorWithVertexRule, EdgeColor>,
    ) {
        let mut initial_pdgs = self
            .options
            .initial_pdgs
            .iter()
            .enumerate()
            .map(|(i_n, pdg)| (pdg, i_n + 1))
            .collect::<Vec<_>>();
        let mut final_pdgs = if matches!(self.options.generation_type, GenerationType::CrossSection)
        {
            self.options
                .initial_pdgs
                .iter()
                .enumerate()
                .map(|(i_n, pdg)| (pdg, i_n + 1 + initial_pdgs.len()))
                .collect::<Vec<_>>()
        } else {
            self.options
                .final_pdgs
                .iter()
                .enumerate()
                .map(|(i_n, pdg)| (pdg, i_n + 1 + initial_pdgs.len()))
                .collect::<Vec<_>>()
        };

        let mut all_pdgs = initial_pdgs.clone();
        all_pdgs.extend(final_pdgs.clone());

        let mut new_node_data = vec![];
        let mut new_edge_data = vec![];
        for (i_e, e) in graph.edges().iter().enumerate() {
            // All edges supposed to be incoming at this stage
            assert!(graph.nodes()[e.vertices.1].data.external_tag == 0);
            let p = model.get_particle_from_pdg(e.data.pdg);
            let external_tag = graph.nodes()[e.vertices.0].data.external_tag;
            let is_initial_state = external_tag <= self.options.initial_pdgs.len() as i32;
            let container = if is_initial_state {
                &mut initial_pdgs
            } else {
                &mut final_pdgs
            };
            if external_tag != 0 {
                let pdg_code = if is_initial_state {
                    p.pdg_code
                } else {
                    p.get_anti_particle(model).pdg_code
                };
                if self.options.symmetrize_left_right_states {
                    let matched_external_pos = all_pdgs
                        .iter()
                        .position(|(&pdg, _i_ext)| pdg == pdg_code as i64)
                        .unwrap();
                    let mut new_data = graph.nodes()[e.vertices.0].data.clone();
                    let new_external_tag = all_pdgs[matched_external_pos].1 as i32;
                    // If we swapped initial and final state assignment, then we must also flip the pdg code of the corresponding half-edges
                    if (new_external_tag > initial_pdgs.len() as i32
                        && new_data.external_tag <= initial_pdgs.len() as i32)
                        || (new_external_tag <= initial_pdgs.len() as i32
                            && new_data.external_tag > initial_pdgs.len() as i32)
                    {
                        let mut e_data = e.data.clone();
                        e_data.pdg = model
                            .get_particle_from_pdg(e_data.pdg)
                            .get_anti_particle(model)
                            .pdg_code;
                        new_edge_data.push((i_e, e_data));
                    }
                    new_data.set_external_tag(new_external_tag);
                    new_node_data.push((e.vertices.0, new_data));
                    all_pdgs.remove(matched_external_pos);
                } else {
                    let matched_external_pos = container
                        .iter()
                        .position(|(&pdg, _i_ext)| pdg == pdg_code as i64)
                        .unwrap();
                    let mut new_data = graph.nodes()[e.vertices.0].data.clone();
                    new_data.set_external_tag(container[matched_external_pos].1 as i32);
                    new_node_data.push((e.vertices.0, new_data));
                    container.remove(matched_external_pos);
                }
            }
            /*
            let external_tag = graph.nodes()[e.vertices.1].data.external_tag;
            if external_tag != 0 {
                // let pdg_code = if is_initial_state {
                //     p.get_anti_particle(model).pdg_code
                // } else {
                //     p.pdg_code
                // };
                let pdg_code = p.get_anti_particle(model).pdg_code;
                if self.options.symmetrize_left_right_states {
                    let matched_external_pos = all_pdgs
                        .iter()
                        .position(|(&pdg, _i_ext)| pdg == pdg_code as i64)
                        .unwrap();
                    let mut new_data = graph.nodes()[e.vertices.0].data.clone();
                    new_data.set_external_tag(all_pdgs[matched_external_pos].1 as i32);
                    new_node_data.push((e.vertices.0, new_data));
                    all_pdgs.remove(matched_external_pos);
                } else {
                    let matched_external_pos = container
                        .iter()
                        .position(|(&pdg, _i_ext)| pdg == pdg_code as i64)
                        .unwrap();
                    let mut new_data = graph.nodes()[e.vertices.0].data.clone();
                    new_data.set_external_tag(container[matched_external_pos].1 as i32);
                    new_node_data.push((e.vertices.0, new_data));
                    container.remove(matched_external_pos);
                }
            }
             */
        }

        for (node_pos, new_data) in new_node_data {
            graph.set_node_data(node_pos, new_data);
        }
        for (edge_pos, new_data) in new_edge_data {
            graph.set_edge_data(edge_pos, new_data);
        }
    }

    /// This function canonizes the edge and vertex ordering of a graph based on the skeletton graph with only propagator mass as edge color.
    /// This is useful to then allow for further grouping of isomorphic graphs, incl numerator.
    #[allow(clippy::type_complexity)]
    pub fn canonicalize_edge_and_vertex_ordering(
        &self,
        model: &Model,
        input_graph: &SymbolicaGraph<NodeColorWithVertexRule, EdgeColor>,
        node_colors_for_external_symmetrization: &HashMap<i32, i32>,
        numerator_aware_isomorphism_grouping: &NumeratorAwareGraphGroupingOption,
        manually_canonalize_initial_final_swap: bool,
    ) -> (
        SymbolicaGraph<NodeColorWithoutVertexRule, std::string::String>,
        SymbolicaGraph<NodeColorWithVertexRule, EdgeColor>,
    ) {
        // println!("INPUT GRAPH:\n{}", input_graph.to_dot());
        let mut canonized_skelettons = vec![];
        let external_node_positions = input_graph
            .nodes()
            .iter()
            .enumerate()
            .filter_map(|(i_n, n)| {
                if n.data.external_tag > 0 {
                    Some((n.data.external_tag, i_n))
                } else {
                    None
                }
            })
            .collect::<HashMap<_, _>>();
        let mut external_remappings = vec![AHashMap::<usize, usize>::default()];
        if manually_canonalize_initial_final_swap {
            let mut remapping = AHashMap::default();
            for i_ext in 1..=self.options.initial_pdgs.len() {
                remapping.insert(
                    *external_node_positions.get(&(i_ext as i32)).unwrap(),
                    *external_node_positions
                        .get(&((i_ext + self.options.initial_pdgs.len()) as i32))
                        .unwrap(),
                );
                remapping.insert(
                    *external_node_positions
                        .get(&((i_ext + self.options.initial_pdgs.len()) as i32))
                        .unwrap(),
                    *external_node_positions.get(&(i_ext as i32)).unwrap(),
                );
            }
            external_remappings.push(remapping);
        }
        for remapping in external_remappings {
            // Make sure to canonize the edge ordering based on the skeletton graph with only propagator mass as edge color
            let mut skeletton_graph = SymbolicaGraph::new();
            for node in input_graph.nodes() {
                skeletton_graph.add_node(NodeColorWithoutVertexRule {
                    external_tag: *node_colors_for_external_symmetrization
                        .get(&node.data.external_tag)
                        .unwrap_or(&node.data.external_tag),
                });
            }
            let color_according_to_mass: bool = if let Some(grouping_options) =
                numerator_aware_isomorphism_grouping.get_options()
            {
                grouping_options.differentiate_particle_masses_only
            } else {
                true
            };
            for edge in input_graph.edges() {
                // We must maintain the directionality of the external edges to make sure that antifermions are not mapped to fermions
                // when grouping is done based on the skeletton graph and only the mass of the propagators
                // Also external edges must always retain their PDG since exchanging the momentum of e.g. a massless electron and muon
                // would not be a valid isomorphism w.r.t the process definition
                let is_edge_external = [edge.vertices.0, edge.vertices.1]
                    .iter()
                    .any(|e| input_graph.nodes()[*e].data.external_tag != 0);
                let mut remapped_edge_vertices: (usize, usize) = (
                    *remapping.get(&edge.vertices.0).unwrap_or(&edge.vertices.0),
                    *remapping.get(&edge.vertices.1).unwrap_or(&edge.vertices.1),
                );
                // Force all externals incoming in the canonicalization
                let particle = if input_graph.nodes()[edge.vertices.1].data.external_tag > 0 {
                    remapped_edge_vertices = (remapped_edge_vertices.1, remapped_edge_vertices.0);
                    model
                        .get_particle_from_pdg(edge.data.pdg)
                        .get_anti_particle(model)
                } else {
                    model.get_particle_from_pdg(edge.data.pdg)
                };
                skeletton_graph
                    .add_edge(
                        remapped_edge_vertices.0,
                        remapped_edge_vertices.1,
                        is_edge_external, //edge.directed && is_edge_external,
                        if color_according_to_mass && !is_edge_external {
                            particle.mass.to_string()
                        } else {
                            particle.name.to_string()
                        },
                    )
                    .unwrap();
            }
            canonized_skelettons.push((remapping, skeletton_graph.canonize()));
        }
        canonized_skelettons.sort_by(|a, b| a.1.graph.partial_cmp(&b.1.graph).unwrap());
        let (selected_external_remapping, canonized_skeletton) =
            canonized_skelettons.first().unwrap();

        let mut can_graph_node_pos_to_input_graph_node_pos: Vec<usize> =
            vec![0; input_graph.nodes().len()];
        for (input_graph_node_position, node_order) in
            canonized_skeletton.vertex_map.iter().enumerate()
        {
            can_graph_node_pos_to_input_graph_node_pos[*node_order] = *selected_external_remapping
                .get(&input_graph_node_position)
                .unwrap_or(&input_graph_node_position);
        }
        let mut input_graph_node_pos_to_can_graph_node_pos: Vec<usize> =
            vec![0; input_graph.nodes().len()];
        for (can_graph_node_pos, input_graph_node_pos) in can_graph_node_pos_to_input_graph_node_pos
            .iter()
            .enumerate()
        {
            input_graph_node_pos_to_can_graph_node_pos[*input_graph_node_pos] = can_graph_node_pos;
        }

        // Sort nodes according to the canonized skeletton graph
        // This will also ensure that EMR variables line up
        let mut sorted_g: SymbolicaGraph<NodeColorWithVertexRule, EdgeColor> =
            SymbolicaGraph::new();
        for &node_order in can_graph_node_pos_to_input_graph_node_pos.iter() {
            let node_data = input_graph.nodes()[*selected_external_remapping
                .get(&node_order)
                .unwrap_or(&node_order)]
            .data
            .clone();
            sorted_g.add_node(node_data);
        }

        let input_graph_nodes = input_graph.nodes();
        let mut reordered_edges = input_graph
            .edges()
            .iter()
            .map(|e| {
                // Set all external edges incoming and internal edges going from lower to higher node order
                let is_external = input_graph_nodes[e.vertices.0].data.external_tag > 0
                    || input_graph_nodes[e.vertices.1].data.external_tag > 0;
                let is_external_outgoing = input_graph_nodes[e.vertices.1].data.external_tag > 0;
                let is_flipped = is_external_outgoing
                    || (!is_external
                        && input_graph_node_pos_to_can_graph_node_pos[e.vertices.0]
                            > input_graph_node_pos_to_can_graph_node_pos[e.vertices.1]);

                // let is_flipped = if input_graph_nodes[e.vertices.0].data.external_tag > 0 {
                //     false
                // } else if input_graph_nodes[e.vertices.1].data.external_tag > 0 {
                //     // let (node_data_in, node_data_out) = (
                //     //     input_graph_nodes[e.vertices.0].data.clone(),
                //     //     input_graph_nodes[e.vertices.1].data.clone(),
                //     // );
                //     // sorted_g.set_node_data(e.vertices.0, node_data_out);
                //     // sorted_g.set_node_data(e.vertices.1, node_data_in);
                //     true
                // } else {
                //     canonized_skeletton.vertex_map[e.vertices.0]
                //         > canonized_skeletton.vertex_map[e.vertices.1]
                // };
                /*
                let is_flipped = if input_graph_nodes[e.vertices.0].data.external_tag > 0 {
                    let (node_data_in, node_data_out) = (
                        input_graph_nodes[e.vertices.0].data.clone(),
                        input_graph_nodes[e.vertices.1].data.clone(),
                    );
                    if input_graph_nodes[e.vertices.0].data.external_tag
                        > self.options.initial_pdgs.len() as i32
                    {
                        sorted_g.set_node_data(e.vertices.0, node_data_out);
                        sorted_g.set_node_data(e.vertices.1, node_data_in);
                        true
                    } else {
                        false
                    }
                } else if input_graph_nodes[e.vertices.1].data.external_tag > 0 {
                    let (node_data_in, node_data_out) = (
                        input_graph_nodes[e.vertices.0].data.clone(),
                        input_graph_nodes[e.vertices.1].data.clone(),
                    );
                    if input_graph_nodes[e.vertices.0].data.external_tag
                        <= self.options.initial_pdgs.len() as i32
                    {
                        sorted_g.set_node_data(e.vertices.0, node_data_out);
                        sorted_g.set_node_data(e.vertices.1, node_data_in);
                        true
                    } else {
                        false
                    }
                } else {
                    canonized_skeletton.vertex_map[e.vertices.0]
                        > canonized_skeletton.vertex_map[e.vertices.1]
                };
                */
                /*
                let p_is_antiparticle = !p.is_self_antiparticle() && p.is_antiparticle();
                #[allow(clippy::collapsible_else_if)]
                let is_flipped = if input_graph_nodes[e.vertices.0].data.external_tag > 0 {
                    if input_graph_nodes[e.vertices.0].data.external_tag
                        <= self.options.initial_pdgs.len() as i32
                    {
                        false // p_is_antiparticle
                    } else {
                        true // !p_is_antiparticle
                    }
                } else if input_graph_nodes[e.vertices.1].data.external_tag > 0 {
                    if input_graph_nodes[e.vertices.1].data.external_tag
                        > self.options.initial_pdgs.len() as i32
                    {
                        false // p_is_antiparticle
                    } else {
                        true // !p_is_antiparticle
                    }
                } else {
                    canonized_skeletton.vertex_map[e.vertices.0]
                        > canonized_skeletton.vertex_map[e.vertices.1]
                };
                */
                /*
                let is_flipped = if input_graph_nodes[e.vertices.0].data.external_tag > 0
                    || input_graph_nodes[e.vertices.1].data.external_tag > 0
                {
                    canonized_skeletton.vertex_map[e.vertices.0]
                        > canonized_skeletton.vertex_map[e.vertices.1]
                } else {
                    false
                };
                 */

                let mut particle = if is_flipped {
                    model
                        .get_particle_from_pdg(e.data.pdg)
                        .get_anti_particle(model)
                } else {
                    model.get_particle_from_pdg(e.data.pdg)
                };
                // Apply the switch from particle to anti-particle if the canonilazation swapped initial and final states.
                if !selected_external_remapping.is_empty() && is_external {
                    particle = particle.get_anti_particle(model);
                }
                (
                    (e, is_flipped),
                    // This key will serve to give a unique ordering of the edges
                    (
                        if !is_flipped {
                            input_graph_node_pos_to_can_graph_node_pos[e.vertices.0]
                        } else {
                            input_graph_node_pos_to_can_graph_node_pos[e.vertices.1]
                        },
                        if !is_flipped {
                            input_graph_node_pos_to_can_graph_node_pos[e.vertices.1]
                        } else {
                            input_graph_node_pos_to_can_graph_node_pos[e.vertices.0]
                        },
                        particle.pdg_code,
                    ),
                )
            })
            .collect::<Vec<_>>();
        reordered_edges.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap());

        for &((e, _is_flipped), (v_in, v_out, pdg)) in reordered_edges.iter() {
            // We must canonize the fermion flow as well, so we force the fermion flow to always go from the lower to the higher node order,
            // and adjust the edge colour (particle vs anti-particle) accordingly.
            // This is necessary to ensure that the meaning of the edge momentum representation (which always aligns the monentum with the fermion flow)
            // is preserved
            sorted_g
                .add_edge(
                    v_in,
                    v_out,
                    e.directed,
                    EdgeColor::from_particle(model.get_particle_from_pdg(pdg)),
                )
                .unwrap();
        }

        // In order for the external assignment momenta to line up between the canonized versions of these graphs
        self.canonize_external_momenta_assignment(model, &mut sorted_g);

        (canonized_skeletton.graph.to_owned(), sorted_g)
    }

    #[allow(clippy::too_many_arguments)]
    pub fn generate(
        &self,
        model: &Model,
        numerator_aware_isomorphism_grouping: &NumeratorAwareGraphGroupingOption,
        filter_self_loop: bool,
        graph_prefix: String,
        selected_graphs: Option<Vec<String>>,
        vetoed_graphs: Option<Vec<String>>,
        loop_momentum_bases: Option<HashMap<String, Vec<String>>>,
        numerator_global_prefactor: GlobalPrefactor,
        num_threads: Option<usize>,
    ) -> Result<Vec<BareGraph>, FeynGenError> {
        let progress_bar_style = ProgressStyle::with_template(
            "[{elapsed_precise} | ETA: {eta_precise}] {bar:40.cyan/blue} {pos:>7}/{len:7} ({percent}%) {msg}",
        )
        .unwrap();

        // Build a custom ThreadPool. If `Some(threads)` is given, use that;
        // otherwise, default to the number of logical CPUs on the system.
        let pool = match num_threads {
            Some(threads) => ThreadPoolBuilder::new()
                .num_threads(threads)
                .build()
                .expect("Failed to build custom Rayon ThreadPool"),
            None => ThreadPoolBuilder::new()
                .build()
                .expect("Failed to build default Rayon ThreadPool"),
        };

        debug!(
            "Generating Feynman diagrams over {} cores for model {} and process:\n{}",
            pool.current_num_threads(),
            model.name,
            self.options
        );

        let filters = match self.options.generation_type {
            GenerationType::Amplitude => &self.options.amplitude_filters,
            GenerationType::CrossSection => &self.options.cross_section_filters,
        };

        #[allow(clippy::type_complexity)]
        let mut vertex_signatures: HashMap<
            Vec<(Option<bool>, SmartString<LazyCompact>)>,
            Vec<SmartString<LazyCompact>>,
        > = HashMap::default();
        'add_vertex_rules: for vertex_rule in model.vertex_rules.iter() {
            let mut oriented_particles = vec![];
            for p in vertex_rule.particles.iter() {
                if p.is_self_antiparticle() {
                    oriented_particles.push((None, p.name.clone()));
                } else if p.is_antiparticle() {
                    oriented_particles.push((Some(true), p.get_anti_particle(model).name.clone()));
                } else {
                    oriented_particles.push((Some(false), p.name.clone()));
                }
                if let Some(vetoed_particles) = filters.get_particle_vetos() {
                    if vetoed_particles.contains(&(p.pdg_code as i64))
                        || vetoed_particles.contains(&(p.get_anti_particle(model).pdg_code as i64))
                    {
                        continue 'add_vertex_rules;
                    }
                }
            }
            vertex_signatures
                .entry(oriented_particles)
                .or_default()
                .push(vertex_rule.name.clone());
        }

        let mut external_edges = self
            .options
            .initial_pdgs
            .iter()
            .enumerate()
            .map(|(i_initial, pdg)| {
                let p = model.get_particle_from_pdg(*pdg as isize);
                (
                    NodeColorWithoutVertexRule {
                        external_tag: (i_initial + 1) as i32,
                    },
                    if p.is_self_antiparticle() {
                        (None, EdgeColor::from_particle(p))
                    } else if p.is_antiparticle() {
                        (
                            Some(false),
                            EdgeColor::from_particle(p.get_anti_particle(model)),
                        )
                    } else {
                        (Some(true), EdgeColor::from_particle(p))
                    },
                )
            })
            .collect::<Vec<_>>();

        let mut external_connections = vec![];
        match self.options.generation_type {
            GenerationType::Amplitude => {
                for i_initial in 1..=self.options.initial_pdgs.len() {
                    external_connections.push((Some(i_initial), None));
                }
                for i_final in (self.options.initial_pdgs.len() + 1)
                    ..=(self.options.initial_pdgs.len() + self.options.final_pdgs.len())
                {
                    external_connections.push((None, Some(i_final)));
                }
                external_edges.extend(
                    self.options
                        .final_pdgs
                        .iter()
                        .enumerate()
                        .map(|(i_final, pdg)| {
                            let p = model.get_particle_from_pdg(*pdg as isize);
                            (
                                NodeColorWithoutVertexRule {
                                    external_tag: (self.options.initial_pdgs.len() + i_final + 1)
                                        as i32,
                                },
                                if p.is_self_antiparticle() {
                                    (None, EdgeColor::from_particle(p))
                                } else if p.is_antiparticle() {
                                    (
                                        Some(true),
                                        EdgeColor::from_particle(p.get_anti_particle(model)),
                                    )
                                } else {
                                    (Some(false), EdgeColor::from_particle(p))
                                },
                            )
                        })
                        .collect::<Vec<_>>(),
                );
            }
            GenerationType::CrossSection => {
                let mut i_final = self.options.initial_pdgs.len();
                for (i_initial, pdg) in self.options.initial_pdgs.iter().enumerate() {
                    i_final += 1;
                    let p = model.get_particle_from_pdg(*pdg as isize);
                    external_connections.push((Some(i_initial + 1), Some(i_final)));
                    external_edges.push((
                        NodeColorWithoutVertexRule {
                            external_tag: i_final as i32,
                        },
                        if p.is_self_antiparticle() {
                            (None, EdgeColor::from_particle(p))
                        } else if p.is_antiparticle() {
                            (
                                Some(true),
                                EdgeColor::from_particle(p.get_anti_particle(model)),
                            )
                        } else {
                            (Some(false), EdgeColor::from_particle(p))
                        },
                    ));
                }
            }
        }

        // debug!("external_edges = {:?}", external_edges);
        // debug!("vertex_signatures = {:?}", vertex_signatures);

        // let external_edges_for_generation = external_edges
        //     .iter()
        //     .map(|(i, (orientation, name))| (i.clone(), (*orientation, name)))
        //     .collect::<Vec<_>>();
        let external_edges_for_generation = external_edges.clone();
        let vertex_signatures_for_generation = vertex_signatures
            .keys()
            .map(|v| {
                v.iter()
                    .map(|(orientation, p)| {
                        (
                            *orientation,
                            EdgeColor::from_particle(model.get_particle(p)),
                        )
                    })
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();
        debug!(
            "generation_external_edges = {:?}",
            external_edges_for_generation
        );
        debug!(
            "generation_vertex_signatures = {:?}",
            vertex_signatures_for_generation
        );

        debug!("feygen options: {:?}", self);

        // Record the start time
        let start = Instant::now();
        let mut last_step = start;
        info!(
            "{} | Δ={} | {:<95}",
            format!("{:<6}", utils::format_wdhms(0)).blue().bold(),
            format!("{:<6}", utils::format_wdhms(0)).blue(),
            "Starting Feynman graph generation with Symbolica..."
        );
        let mut graphs = SymbolicaGraph::generate(
            external_edges_for_generation.as_slice(),
            vertex_signatures_for_generation.as_slice(),
            None,
            Some(self.options.loop_count_range.1),
            filters.get_max_bridge(),
            !filter_self_loop,
        );

        // Immediately drop lower loop count contributions
        graphs.retain(|g, _| g.num_loops() >= self.options.loop_count_range.0);
        let mut step = Instant::now();
        info!(
            "{} | Δ={} | {:<95}{}",
            format!("{:<6}", utils::format_wdhms_from_duration(step - start))
                .blue()
                .bold(),
            format!("{:<6}", utils::format_wdhms_from_duration(step - last_step)).blue(),
            "Number of graphs resulting from Symbolica generation:",
            format!("{}", graphs.len()).green().bold()
        );
        last_step = step;

        let unoriented_final_state_particles = if self.options.generation_type
            == GenerationType::CrossSection
            && !self.options.final_pdgs.is_empty()
        {
            self.options
                .final_pdgs
                .iter()
                .map(|pdg| {
                    let p = model.get_particle_from_pdg(*pdg as isize);
                    if p.is_antiparticle() {
                        p.get_anti_particle(model).pdg_code
                    } else {
                        p.pdg_code
                    }
                })
                .collect::<Vec<_>>()
        } else {
            vec![]
        };

        if !unoriented_final_state_particles.is_empty() {
            let bar = ProgressBar::new(graphs.len() as u64);
            bar.set_style(progress_bar_style.clone());
            bar.set_message("Enforcing particle content...");
            pool.install(|| {
                graphs = graphs
                    .par_iter()
                    .progress_with(bar.clone())
                    .filter(|(g, _symmetry_factor)| {
                        FeynGen::contains_particles(g, unoriented_final_state_particles.as_slice())
                    })
                    .map(|(g, sf)| (g.clone(), sf.clone()))
                    .collect::<HashMap<_, _>>()
            });
            bar.finish_and_clear();

            step = Instant::now();
            info!(
                "{} | Δ={} | {:<95}{}",
                format!("{:<6}", utils::format_wdhms_from_duration(step - start))
                    .blue()
                    .bold(),
                format!("{:<6}", utils::format_wdhms_from_duration(step - last_step)).blue(),
                "Number of graphs retained after enforcing supergraph particle content:",
                format!("{}", graphs.len()).green()
            );
            last_step = step;
        }

        let tadpole_filter = filters
            .0
            .iter()
            .filter_map(|f| {
                if let FeynGenFilter::TadpolesFilter(opts) = f {
                    Some(opts)
                } else {
                    None
                }
            })
            .next();
        let external_self_energy_filter = filters
            .0
            .iter()
            .filter_map(|f| {
                if let FeynGenFilter::SelfEnergyFilter(opts) = f {
                    Some(opts)
                } else {
                    None
                }
            })
            .next();
        let zero_snails_filter = filters
            .0
            .iter()
            .filter_map(|f| {
                if let FeynGenFilter::ZeroSnailsFilter(opts) = f {
                    Some(opts)
                } else {
                    None
                }
            })
            .next();
        if tadpole_filter.is_some()
            || external_self_energy_filter.is_some()
            || zero_snails_filter.is_some()
        {
            let bar = ProgressBar::new(graphs.len() as u64);
            bar.set_style(progress_bar_style.clone());
            bar.set_message("Vetoing special topologies...");
            pool.install(|| {
                graphs = graphs
                    .par_iter()
                    .progress_with(bar.clone())
                    .filter(|(g, _symmetry_factor)| {
                        !FeynGen::veto_special_topologies(
                            model,
                            g,
                            external_self_energy_filter,
                            tadpole_filter,
                            zero_snails_filter,
                        )
                    })
                    .map(|(g, sf)| (g.clone(), sf.clone()))
                    .collect::<HashMap<_, _>>()
            });
            bar.finish_and_clear();
        }

        // Re-interpret the symmetry factor as a multiplier where the symbolica factor from symbolica appears in the denominator
        let graphs = graphs
            .iter()
            .map(|(g, symmetry_factor)| {
                (
                    g.clone(),
                    Atom::new_num(1) / Atom::new_num(symmetry_factor.clone()),
                )
            })
            .collect::<HashMap<_, _>>();

        step = Instant::now();
        info!(
            "{} | Δ={} | {:<95}{}",
            format!("{:<6}", utils::format_wdhms_from_duration(step - start))
                .blue()
                .bold(),
            format!("{:<6}", utils::format_wdhms_from_duration(step - last_step)).blue(),
            "Number of graphs retained after removal of vetoed topologies:",
            format!("{}", graphs.len()).green()
        );
        last_step = step;

        #[allow(clippy::type_complexity)]
        let mut node_colors: HashMap<
            Vec<(Option<bool>, SmartString<LazyCompact>)>,
            Vec<SmartString<LazyCompact>>,
        > = HashMap::default();
        for (v_legs, v_colors) in vertex_signatures.iter() {
            let mut sorted_ps = v_legs.clone();
            sorted_ps.sort();
            node_colors.insert(sorted_ps, v_colors.clone());
        }

        let mut processed_graphs = vec![];
        for (g, symmetry_factor) in graphs.iter() {
            for (colored_g, multiplicity) in FeynGen::assign_node_colors(model, g, &node_colors)? {
                processed_graphs.push((
                    colored_g.canonize().graph,
                    (Atom::new_num(multiplicity as i64) * symmetry_factor).to_canonical_string(),
                ));
            }
        }

        let bar = ProgressBar::new(graphs.len() as u64);
        bar.set_style(progress_bar_style.clone());
        bar.set_message("Assigning interactions to vertices...");
        let mut processed_graphs = vec![];
        let lists_of_entries = pool.install(|| {
            graphs
                .par_iter()
                .progress_with(bar.clone())
                .map(|(g, symmetry_factor)| {
                    FeynGen::assign_node_colors(model, g, &node_colors).map(|colored_graphs| {
                        colored_graphs
                            .iter()
                            .map(|(colored_g, multiplicity)| {
                                (
                                    colored_g.canonize().graph,
                                    (Atom::new_num(*multiplicity as i64) * symmetry_factor)
                                        .to_canonical_string(),
                                )
                            })
                            .collect::<Vec<_>>()
                    })
                })
                .collect::<Result<Vec<_>, FeynGenError>>()
        })?;
        bar.finish_and_clear();

        for entries in lists_of_entries {
            processed_graphs.extend(entries);
        }
        processed_graphs.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());

        step = Instant::now();
        info!(
            "{} | Δ={} | {:<95}{}",
            format!("{:<6}", utils::format_wdhms_from_duration(step - start))
                .blue()
                .bold(),
            format!("{:<6}", utils::format_wdhms_from_duration(step - last_step)).blue(),
            "Number of graphs after vertex info assignment:",
            format!("{}", processed_graphs.len()).green()
        );
        last_step = step;

        filters.apply_filters(&mut processed_graphs, &pool, &progress_bar_style)?;

        step = Instant::now();
        info!(
            "{} | Δ={} | {:<95}{}",
            format!("{:<6}", utils::format_wdhms_from_duration(step - start))
                .blue()
                .bold(),
            format!("{:<6}", utils::format_wdhms_from_duration(step - last_step)).blue(),
            "Number of graphs after all complete graphs filters are applied:",
            format!("{}", processed_graphs.len()).green()
        );
        last_step = step;

        let (n_unresolved, unresolved_type) = self.unresolved_cut_content(model);

        // The fast cutksoky filter is only fast for up to ~ 6 particles to check
        const MAX_FAST_CUTKOSKY_PARTICLES: usize = 6;
        let mut applied_fast_cutksosky_cut_filter = false;
        if self.options.generation_type == GenerationType::CrossSection
            && !self.options.final_pdgs.is_empty()
            && self.options.final_pdgs.len() + n_unresolved <= MAX_FAST_CUTKOSKY_PARTICLES
        {
            applied_fast_cutksosky_cut_filter = true;

            let bar = ProgressBar::new(processed_graphs.len() as u64);
            bar.set_style(progress_bar_style.clone());
            bar.set_message("Applying fast Cutkosky cut filter...");
            pool.install(|| {
                processed_graphs = processed_graphs
                    .par_iter()
                    .progress_with(bar.clone())
                    .filter(|(g, _symmetry_factor)| {
                        self.contains_cut_fast(
                            model,
                            g,
                            n_unresolved,
                            &unresolved_type,
                            unoriented_final_state_particles.as_slice(),
                        )
                    })
                    .map(|(g, sf)| (g.clone(), sf.clone()))
                    .collect::<Vec<_>>()
            });
            bar.finish_and_clear();

            step = Instant::now();
            info!(
                "{} | Δ={} | {:<95}{}",
                format!("{:<6}", utils::format_wdhms_from_duration(step - start))
                    .blue()
                    .bold(),
                format!("{:<6}", utils::format_wdhms_from_duration(step - last_step)).blue(),
                "Number of graphs after fast Cutkosky cut filter is applied:",
                format!("{}", processed_graphs.len()).green()
            );
            last_step = step;
        }

        // This secondary cutkosky cut filter is always necessary as the fast one can have false positives
        // even in the absence of constraints on amplitudes on either side of the cut
        if self.options.generation_type == GenerationType::CrossSection
            && !self.options.final_pdgs.is_empty()
            && (self
                .options
                .amplitude_filters
                .get_coupling_orders()
                .is_some()
                || self
                    .options
                    .amplitude_filters
                    .get_loop_count_range()
                    .is_some()
                || !applied_fast_cutksosky_cut_filter)
        {
            // if self.options.generation_type == GenerationType::CrossSection
            //     && !self.options.final_pdgs.is_empty()
            // {
            let bar = ProgressBar::new(processed_graphs.len() as u64);
            bar.set_style(progress_bar_style.clone());
            bar.set_message("Applying secondary exact Cutkosky cut filter...");

            pool.install(|| {
                processed_graphs = processed_graphs
                    .par_iter()
                    .progress_with(bar.clone())
                    .filter(|(g, _)| self.contains_cut(model, g, n_unresolved, &unresolved_type))
                    .map(|(g, sf)| (g.clone(), sf.clone()))
                    .collect::<Vec<_>>()
            });
            bar.finish_and_clear();

            step = Instant::now();
            info!(
                "{} | Δ={} | {:<95}{}",
                format!("{:<6}", utils::format_wdhms_from_duration(step - start))
                    .blue()
                    .bold(),
                format!("{:<6}", utils::format_wdhms_from_duration(step - last_step)).blue(),
                "Number of graphs after exact Cutkosky cut filter is applied:",
                format!("{}", processed_graphs.len()).green()
            );
            last_step = step;
        }

        // Now account for initial state symmetry by further grouping contributions
        let mut node_colors_for_external_symmetrization: HashMap<i32, i32> = HashMap::default();
        // Because of the interplay with the cutkosky cut filter and left-right canonization when using symmetrize_left_right_states
        // we must do to canonicalizations here and we will select the "smallest one"
        let mut node_colors_for_canonicalization_left_right: HashMap<i32, i32> = HashMap::default();
        if self.options.symmetrize_initial_states {
            for initial_color in 1..=self.options.initial_pdgs.len() {
                let trgt_color = if self.options.symmetrize_left_right_states {
                    -1
                } else {
                    -2
                };
                node_colors_for_external_symmetrization.insert(initial_color as i32, trgt_color);
                if self.options.generation_type == GenerationType::Amplitude {
                    node_colors_for_canonicalization_left_right
                        .insert(initial_color as i32, trgt_color);
                } else {
                    node_colors_for_canonicalization_left_right.insert(initial_color as i32, -2);
                }
            }
            if self.options.generation_type == GenerationType::CrossSection {
                for final_color in
                    self.options.initial_pdgs.len() + 1..=2 * self.options.initial_pdgs.len()
                {
                    let trgt_color = if self.options.symmetrize_left_right_states {
                        -1
                    } else {
                        -3
                    };
                    node_colors_for_external_symmetrization.insert(final_color as i32, trgt_color);
                    node_colors_for_canonicalization_left_right.insert(final_color as i32, -3);
                }
            }
        }
        if self.options.generation_type == GenerationType::Amplitude
            && self.options.symmetrize_final_states
        {
            for final_color in self.options.initial_pdgs.len() + 1
                ..=self.options.initial_pdgs.len() + self.options.final_pdgs.len()
            {
                let trgt_color = if self.options.symmetrize_left_right_states {
                    -1
                } else {
                    -3
                };
                node_colors_for_external_symmetrization.insert(final_color as i32, trgt_color);
                node_colors_for_canonicalization_left_right.insert(final_color as i32, trgt_color);
            }
        }
        // We must disable the naive grouping for cross-sections when symmetrizing the left-right states
        // because we would wrongly count permutations that could be disabled by the cutkosky cut filter
        if !node_colors_for_external_symmetrization.is_empty() {
            processed_graphs = FeynGen::group_isomorphic_graphs_after_node_color_change(
                &processed_graphs
                    .iter()
                    .map(|(g, m)| (g.clone(), Atom::parse(m).unwrap()))
                    .collect::<HashMap<_, _>>(),
                &node_colors_for_canonicalization_left_right,
                &pool,
                &progress_bar_style,
            )
            .iter()
            .map(|(g, m)| (g.clone(), m.to_canonical_string()))
            .collect::<Vec<_>>();
        }

        step = Instant::now();
        info!(
            "{} | Δ={} | {:<95}{}",
            format!("{:<6}", utils::format_wdhms_from_duration(step - start))
                .blue()
                .bold(),
            format!("{:<6}", utils::format_wdhms_from_duration(step - last_step)).blue(),
            "Number of graphs after symmetrization of external states:",
            format!("{}", processed_graphs.len()).green()
        );
        last_step = step;

        let bar = ProgressBar::new(graphs.len() as u64);
        bar.set_style(progress_bar_style.clone());
        bar.set_message("Canonalizing edge and ndoes ordering of selected graphs, including leg symmetrization...");
        let mut canonized_processed_graphs = pool.install(|| {
            processed_graphs
                .par_iter()
                .progress_with(bar.clone())
                .map(|(g, symmetry_factor)| {
                    let manually_canonalize_initial_final_swap = self.options.generation_type
                        == GenerationType::CrossSection
                        && self.options.symmetrize_left_right_states;
                    let (mut canonical_repr, mut sorted_g) = self
                        .canonicalize_edge_and_vertex_ordering(
                            model,
                            g,
                            &node_colors_for_canonicalization_left_right,
                            numerator_aware_isomorphism_grouping,
                            manually_canonalize_initial_final_swap,
                        );
                    // If we are symmetrizing the left-right states in the context of a cross-section, the canonaliztion above
                    // has canonized the choice of which nodes to assign to the initial and final state.
                    // We now do a second pass to canonalize the vertex ordering for that particular choice.
                    if manually_canonalize_initial_final_swap {
                        (canonical_repr, sorted_g) = self.canonicalize_edge_and_vertex_ordering(
                            model,
                            &sorted_g,
                            &node_colors_for_canonicalization_left_right,
                            numerator_aware_isomorphism_grouping,
                            false,
                        );
                    }
                    (canonical_repr, sorted_g, symmetry_factor.clone())
                })
                .collect::<Vec<_>>()
        });
        bar.finish_and_clear();

        canonized_processed_graphs.sort_by(|a, b| (a.1).partial_cmp(&b.1).unwrap());

        let n_zeroes = Arc::new(Mutex::new(0));
        let n_zeroes_clone = n_zeroes.clone();
        let n_groupings = Arc::new(Mutex::new(0));
        let n_groupings_clone = n_groupings.clone();
        #[allow(clippy::type_complexity)]
        let pooled_bare_graphs: Arc<
            Mutex<
                HashMap<
                    SymbolicaGraph<NodeColorWithoutVertexRule, std::string::String>,
                    Vec<(usize, Option<ProcessedNumeratorForComparison>, BareGraph)>,
                >,
            >,
        > = Arc::new(Mutex::new(HashMap::default()));
        let pooled_bare_graphs_clone = pooled_bare_graphs.clone();

        let bar = ProgressBar::new(canonized_processed_graphs.len() as u64);
        bar.set_style(progress_bar_style.clone());
        bar.set_message(format!(
            "Final numerator-aware processing of remaining graphs ({} found: {} | {} found: {})...",
            "#zeros".green(),
            format!("{}", 0).green().bold(),
            "#groupings".green(),
            format!("{}", 0).green().bold(),
        ));
        pool.install(|| {
            canonized_processed_graphs
                .par_iter()
                .progress_with(bar.clone())
                .enumerate()
                .map(|(i_g, (canonical_repr, g, symmetry_factor))| {
                    let graph_name = format!("{}{}", graph_prefix, i_g);
                    if let Some(selected_graphs) = &selected_graphs {
                        if !selected_graphs.contains(&graph_name) {
                            return Ok(())
                        }
                    }
                    if let Some(vetoed_graphs) = &vetoed_graphs {
                        if vetoed_graphs.contains(&graph_name) {
                            return Ok(())
                        }
                    }

                    let bare_graph = BareGraph::from_symbolica_graph(
                        model,
                        graph_name,
                        g,
                        symmetry_factor.clone(),
                        external_connections.clone(),
                        None,
                    )?;
                    // When disabling numerator-aware graph isomorphism, each graph is added separately
                    if matches!(
                        numerator_aware_isomorphism_grouping,
                        NumeratorAwareGraphGroupingOption::NoGrouping
                    ) {
                        {
                            let mut pooled_bare_graphs_lock = pooled_bare_graphs_clone.lock().unwrap();
                            match pooled_bare_graphs_lock.entry(canonical_repr.clone()) {
                                Entry::Vacant(entry) => {
                                    entry.insert(vec![(i_g, None, bare_graph)]);
                                }
                                Entry::Occupied(mut entry) => {
                                    entry.get_mut().push((i_g, None, bare_graph));
                                }
                            }
                        }
                    } else {
                        let numerator =
                            Numerator::default().from_graph(&bare_graph, &numerator_global_prefactor);
                        let numerator_color_simplified: Numerator<SymbolicExpression<Color>> =
                            numerator.color_simplify();
                        if numerator_color_simplified
                            .get_single_atom()
                            .unwrap()
                            .0
                            .is_zero()
                        {
                            {
                                let mut n_zeroes_value = n_zeroes_clone.lock().unwrap();
                                let n_groupings_value = n_groupings_clone.lock().unwrap();
                                *n_zeroes_value += 1;
                                bar.set_message(format!("Final numerator-aware processing of remaining graphs ({} found: {} | {} found: {})...",
                                    "#zeros".green(),
                                    format!("{}",n_zeroes_value).green().bold(),
                                    "#groupings".green(),
                                    format!("{}",n_groupings_value).green().bold(),
                                ));
                            }
                            return Ok(())
                        }
                        if matches!(
                            numerator_aware_isomorphism_grouping,
                            NumeratorAwareGraphGroupingOption::OnlyDetectZeroes
                        ) {
                            {
                                let mut pooled_bare_graphs_lock = pooled_bare_graphs_clone.lock().unwrap();

                                match pooled_bare_graphs_lock.entry(canonical_repr.clone()) {
                                    Entry::Vacant(entry) => {
                                        entry.insert(vec![(i_g, None, bare_graph)]);
                                    }
                                    Entry::Occupied(mut entry) => {
                                        entry.get_mut().push((i_g, None, bare_graph));
                                    }
                                }
                            }
                        } else {
                            // println!(
                            //     "I have:\n{}",
                            //     numerator_color_simplified
                            //         .clone()
                            //         .parse()
                            //         .contract(ContractionSettings::<Rational>::Normal)
                            //         .unwrap()
                            //         .state
                            //         .tensor
                            //         .iter_flat()
                            //         .map(|(id, d)| format!("{}: {}", id, d))
                            //         .collect::<Vec<_>>()
                            //         .join("\n")
                            // );
                            let numerator_data = Some(
                                ProcessedNumeratorForComparison::from_numerator_symbolic_expression(
                                    i_g,
                                    &bare_graph,
                                    numerator_color_simplified,
                                    numerator_aware_isomorphism_grouping,
                                )?,
                            );
                            // println!("Skeletton G#{}:\n{}", i_g, canonical_repr.to_dot());
                            {
                                let mut pooled_bare_graphs_lock = pooled_bare_graphs_clone.lock().unwrap();

                                match pooled_bare_graphs_lock.entry(canonical_repr.clone()) {
                                    Entry::Vacant(entry) => {
                                        entry.insert(vec![(i_g, numerator_data, bare_graph)]);
                                    }
                                    Entry::Occupied(mut entry) => {
                                        let mut found_match = false;
                                        for (_graph_id, other_numerator, other_graph) in entry.get_mut() {
                                            if let Some(ratio) = FeynGen::compare_numerator_tensors(
                                                numerator_aware_isomorphism_grouping,
                                                numerator_data.as_ref().unwrap(),
                                                other_numerator.as_ref().unwrap(),
                                            ) {
                                                {
                                                    let n_zeroes_value = n_zeroes_clone.lock().unwrap();
                                                    let mut n_groupings_value =
                                                        n_groupings_clone.lock().unwrap();
                                                    *n_groupings_value += 1;
                                                    bar.set_message(format!("Final numerator-aware processing of remaining graphs ({} found: {} | {} found: {})...",
                                                        "#zeros".green(),
                                                        format!("{}",n_zeroes_value).green().bold(),
                                                        "#groupings".green(),
                                                        format!("{}",n_groupings_value).green().bold(),
                                                    ));
                                                }
                                                found_match = true;
                                                other_graph.overall_factor =
                                                    (Atom::parse(other_graph.overall_factor.as_str())
                                                        .unwrap()
                                                        + ratio
                                                            * Atom::parse(
                                                                bare_graph.overall_factor.as_str(),
                                                            )
                                                            .unwrap())
                                                    .expand()
                                                    .to_canonical_string();
                                            }
                                        }
                                        if !found_match {
                                            entry.get_mut().push((i_g, numerator_data, bare_graph));
                                        }
                                    }
                                }
                            }
                        }
                    }
                    Ok(())

                }).collect::<Result<Vec<()>, FeynGenError>>()
        })?;
        bar.finish_and_clear();

        let mut n_cancellations: i32 = 0;
        let (mut bare_graphs, pooled_bare_graphs_len) = {
            let pooled_bare_graphs_lock = pooled_bare_graphs.lock().unwrap();

            (
                pooled_bare_graphs_lock
                    .values()
                    .flatten()
                    .filter_map(|(graph_id, _numerator, graph)| {
                        if Atom::parse(&graph.overall_factor)
                            .unwrap()
                            .expand()
                            .is_zero()
                        {
                            n_cancellations += 1;
                            None
                        } else {
                            Some((*graph_id, graph.clone()))
                        }
                    })
                    .collect::<Vec<_>>(),
                pooled_bare_graphs_lock.len(),
            )
        };
        bare_graphs.sort_by(|a: &(usize, BareGraph), b| (a.0).cmp(&b.0));

        let (n_zeroes_value, n_groupings_value) = {
            (
                *n_zeroes_clone.lock().unwrap(),
                *n_groupings_clone.lock().unwrap(),
            )
        };
        step = Instant::now();
        info!(
            "{} | Δ={} | {:<95}{} ({} isomorphically unique graph{}, {} zero{}, {} grouped and {} cancellation{})",
            format!("{:<6}", utils::format_wdhms_from_duration(step - start)).blue().bold(),
            format!(
                "{:<6}",
                utils::format_wdhms_from_duration(step - last_step)
            )
            .blue(),
            format!(
                "Number of graphs after numerator-aware grouping with strategy '{}':",
                numerator_aware_isomorphism_grouping
            ),
            format!("{}", bare_graphs.len()).green().bold(),
            format!("{}", pooled_bare_graphs_len).blue().bold(),
            if pooled_bare_graphs_len > 1 { "s" } else { "" },
            format!("{}", n_zeroes_value).blue().bold(),
            if n_zeroes_value > 1 { "s" } else { "" },
            format!("{}", n_groupings_value).blue().bold(),
            format!("{}", n_cancellations).blue().bold(),
            if n_cancellations > 1 { "s" } else { "" },
        );

        for (_graph_id, graph) in bare_graphs.iter_mut() {
            let forced_lmb = if let Some(lmbs) = loop_momentum_bases.as_ref() {
                let g_name = String::from(graph.name.clone());
                lmbs.get(&g_name).map(|lmb: &Vec<String>| {
                    lmb.iter().map(SmartString::<LazyCompact>::from).collect()
                })
            } else {
                None
            };
            if forced_lmb.is_some() {
                graph.set_loop_momentum_basis(&forced_lmb)?;
            }
        }
        // println!(
        //     "Graphs: [{}]",
        //     bare_graphs
        //         .iter()
        //         .map(|(_graph_id, graph)| format!("\"{}\"", graph.name.clone()))
        //         .collect::<Vec<_>>()
        //         .join(", "),
        // );

        Ok(bare_graphs
            .iter()
            .map(|(_graph_id, graph)| graph.clone())
            .collect::<Vec<_>>())
    }

    fn compare_numerator_tensors(
        numerator_aware_isomorphism_grouping: &NumeratorAwareGraphGroupingOption,
        numerator_a: &ProcessedNumeratorForComparison,
        numerator_b: &ProcessedNumeratorForComparison,
    ) -> Option<Atom> {
        fn analyze_ratios(ratios: &HashSet<Atom>) -> Option<Atom> {
            if ratios.len() > 1 {
                None
            } else {
                let ratio = ratios.iter().next().unwrap().to_owned();
                for head in ExtendibleReps::BUILTIN_SELFDUAL_NAMES
                    .iter()
                    .chain(ExtendibleReps::BUILTIN_DUALIZABLE_NAMES.iter())
                    .chain(["Q"].iter())
                    .map(Symbol::new)
                {
                    if ratio
                        .pattern_match(&fun!(head, symb!("args__")).to_pattern(), None, None)
                        .next()
                        .is_some()
                    {
                        return None;
                    }
                }
                Some(ratio)
            }
        }

        fn analyze_diff_and_sum(a: AtomView, b: AtomView) -> Option<Atom> {
            if (a - b).expand_num().collect_num().is_zero() {
                return Some(Atom::new_num(1));
            }
            if (a + b).expand_num().collect_num().is_zero() {
                return Some(Atom::new_num(-1));
            }
            None
        }

        fn compare_sample_points(
            sample_points_a: &[Vec<(Atom, Atom)>],
            sample_points_b: &[Vec<(Atom, Atom)>],
        ) -> bool {
            sample_points_a
                .iter()
                .zip(sample_points_b.iter())
                .all(|(sp_a, sp_b)| {
                    sp_a.len() == sp_b.len() && sp_a.iter().zip(sp_b.iter()).all(|(a, b)| a == b)
                })
        }
        // if numerator_a.diagram_id == 22 || numerator_b.diagram_id == 22 {
        //     println!(
        //         "STARTING COMPARISON of diag #{} and #{}",
        //         numerator_b.diagram_id, numerator_a.diagram_id,
        //     );
        // }
        match numerator_aware_isomorphism_grouping {
            NumeratorAwareGraphGroupingOption::NoGrouping
            | NumeratorAwareGraphGroupingOption::OnlyDetectZeroes => None,
            NumeratorAwareGraphGroupingOption::GroupIdenticalGraphUpToScalarRescaling(
                grouping_options,
            ) => {
                if grouping_options.test_canonized_numerator {
                    if let (Some(canonized_num_a), Some(canonized_num_b)) = (
                        numerator_a.canonized_numerator.as_ref(),
                        numerator_b.canonized_numerator.as_ref(),
                    ) {
                        let mut ratios = HashSet::<Atom>::default();
                        let r = if canonized_num_a == canonized_num_b {
                            Atom::new_num(1)
                        } else if *canonized_num_a == canonized_num_b * Atom::new_num(-1) {
                            Atom::new_num(-1)
                        } else {
                            (canonized_num_a.as_view() / canonized_num_b.as_view())
                                .expand_num()
                                .collect_num()
                        };
                        ratios.insert(r);
                        // println!(
                        //     "Canonalized numerator of diagiaram #{}: {}",
                        //     numerator_a.diagram_id, numerator_a.canonized_numerator
                        // );
                        // println!(
                        //     "Canonalized numerator of diagiaram #{}: {}",
                        //     numerator_b.diagram_id, numerator_b.canonized_numerator
                        // );
                        // println!(
                        //     "ratio from canonalized numerators: {:?}",
                        //     ratios
                        //         .iter()
                        //         .map(|av| av.to_canonical_string())
                        //         .collect::<Vec<_>>()
                        //         .join(",")
                        // );
                        if let Some(ratio) = analyze_ratios(&ratios) {
                            debug!(
                                "Combining graph #{} with #{} using canonized numerators, with ratio #{}/#{}={}",
                                numerator_b.diagram_id,
                                numerator_a.diagram_id,
                                numerator_a.diagram_id,
                                numerator_b.diagram_id,
                                ratio
                            );
                            return Some(ratio);
                        }
                    }
                }
                if grouping_options.number_of_numerical_samples > 0 {
                    if let (
                        (Some(sample_points_a), Some(evaluations_a)),
                        (Some(sample_points_b), Some(evaluations_b)),
                    ) = (
                        (
                            numerator_a.sample_points.as_ref(),
                            numerator_a.sample_evaluations.as_ref(),
                        ),
                        (
                            numerator_b.sample_points.as_ref(),
                            numerator_b.sample_evaluations.as_ref(),
                        ),
                    ) {
                        // Make sure the variables detected are the same to begin with
                        if compare_sample_points(sample_points_a, sample_points_b) {
                            // println!(
                            //     "Sample evaluations a:\n{}",
                            //     evaluations_a
                            //         .iter()
                            //         .map(|av| av.to_canonical_string())
                            //         .collect::<Vec<_>>()
                            //         .join("\n")
                            // );
                            // println!(
                            //     "Sample evaluations b:\n{}",
                            //     evaluations_b
                            //         .iter()
                            //         .map(|av| av.to_canonical_string())
                            //         .collect::<Vec<_>>()
                            //         .join("\n")
                            // );
                            let ratios = evaluations_a
                                .iter()
                                .zip(evaluations_b.iter())
                                .map(|(a, b)| (a / b).expand_num().collect_num())
                                .collect::<HashSet<_>>();
                            // if numerator_a.diagram_id == 22 || numerator_b.diagram_id == 22 {
                            //     println!(
                            //     "ratios from numerical evaluations when comparing diagram #{} and #{}:\n{}",
                            //     numerator_b.diagram_id,
                            //     numerator_a.diagram_id,
                            //     ratios
                            //         .iter()
                            //         .map(|av| av.to_canonical_string())
                            //         .collect::<Vec<_>>()
                            //         .join("\n")
                            // );
                            // }
                            if let Some(ratio) = analyze_ratios(&ratios) {
                                //     debug!(
                                //     "Combining graph #{} with #{} using numerical evaluation, with ratio #{}/#{}={}",
                                //     numerator_b.diagram_id,
                                //     numerator_a.diagram_id,
                                //     numerator_a.diagram_id,
                                //     numerator_b.diagram_id,
                                //     ratio
                                // );
                                return Some(ratio);
                            }
                        } else {
                            warn!("Skipping comparison of numerical samples between diagrams #{} and #{} because their variables differ", numerator_b.diagram_id, numerator_a.diagram_id);
                            warn!(
                                "Sample points A:\n{}",
                                sample_points_a
                                    .first()
                                    .unwrap()
                                    .iter()
                                    .map(|(a, b)| format!("{} -> {}", a, b))
                                    .join("\n")
                            );
                            warn!(
                                "Sample points B:\n{}",
                                sample_points_b
                                    .first()
                                    .unwrap()
                                    .iter()
                                    .map(|(a, b)| format!("{} -> {}", a, b))
                                    .join("\n")
                            );
                        }
                    }
                }
                None
            }
            NumeratorAwareGraphGroupingOption::GroupIdenticalGraphUpToSign(grouping_options) => {
                if grouping_options.test_canonized_numerator {
                    if let (Some(canonized_num_a), Some(canonized_num_b)) = (
                        numerator_a.canonized_numerator.as_ref(),
                        numerator_b.canonized_numerator.as_ref(),
                    ) {
                        if let Some(ratio) = analyze_diff_and_sum(
                            canonized_num_a.as_view(),
                            canonized_num_b.as_view(),
                        ) {
                            debug!(
                            "Combining graph #{} with #{} using canonized numerators, with ratio #{}/#{}={}",
                            numerator_b.diagram_id,
                            numerator_a.diagram_id,
                            numerator_a.diagram_id,
                            numerator_b.diagram_id,
                            ratio
                        );
                            return Some(ratio);
                        }
                    }
                }
                if grouping_options.number_of_numerical_samples > 0 {
                    if let (
                        (Some(sample_points_a), Some(evaluations_a)),
                        (Some(sample_points_b), Some(evaluations_b)),
                    ) = (
                        (
                            numerator_a.sample_points.as_ref(),
                            numerator_a.sample_evaluations.as_ref(),
                        ),
                        (
                            numerator_b.sample_points.as_ref(),
                            numerator_b.sample_evaluations.as_ref(),
                        ),
                    ) {
                        if compare_sample_points(sample_points_a, sample_points_b) {
                            let ratios = evaluations_a
                                .iter()
                                .zip(evaluations_b.iter())
                                .map(|(a, b)| analyze_diff_and_sum(a.as_view(), b.as_view()))
                                .collect::<HashSet<_>>();
                            if ratios.len() == 1 {
                                if let Some(ratio) = ratios.iter().next().unwrap().to_owned() {
                                    debug!(
                                    "Combining graph #{} with #{} using numerical evaluation, with ratio = {}",
                                    numerator_b.diagram_id,
                                    numerator_a.diagram_id,
                                    ratio
                                );
                                    return Some(ratio);
                                }
                            }
                        } else {
                            //println!("Skipping comparison of numerical samples between diagrams #{} and #{} because their variables differ", numerator_b.diagram_id, numerator_a.diagram_id);
                        }
                    }
                }
                None
            }
        }
    }

    /*
    pub fn remove_external_polarization_vectors(
        numerator: symbolica::atom::AtomView,
    ) -> (Atom, AHashMap<Atom, Vec<Atom>>) {
        let mut new_numerator = numerator.to_owned();
        let arc_mutex_external_indices = Arc::new(std::sync::Mutex::new(AHashMap::new()));
        for pol_head in crate::momentum::POLARIZATION_TYPES {
            let arc_mutex_external_indices_sent = arc_mutex_external_indices.clone();
            let pattern = fun!(symb!(pol_head.to_string()), symb!("args__")).to_pattern();
            new_numerator = numerator.replace_all(
                &pattern,
                symbolica::id::Pattern::Transformer(Box::new((
                    Some(Atom::new_num(1).to_pattern()),
                    vec![symbolica::transformer::Transformer::Map(Box::new(
                        move |a_in: symbolica::atom::AtomView, a_out: &mut Atom| {
                            if let symbolica::atom::AtomView::Fun(s) = a_in {
                                let mut args_iter = s.to_slice().iter();
                                // Skip the first argument
                                args_iter.next();
                                args_iter.for_each(|idx| {
                                    if let Some(m) = idx
                                        .pattern_match(
                                            &Atom::parse("idxType_(idx_)").unwrap().to_pattern(),
                                            None,
                                            None,
                                        )
                                        .next()
                                    {
                                        arc_mutex_external_indices_sent
                                            .lock()
                                            .unwrap()
                                            .entry(m[symb!("idxType_")])
                                            .or_insert_with(|| vec![])
                                            .push(m[symb!("idx_")]);
                                    }
                                });
                            }
                            *a_out = a_in.to_owned();
                            Ok(())
                        },
                    ))],
                ))),
                None,
                None,
            );
        }
        let external_indices = arc_mutex_external_indices.lock().unwrap().clone();
        for (idx_type, idxs) in external_indices.iter_mut() {
            idxs.sort();
        }
        (new_numerator, external_indices)
    }
     */
}

struct ProcessedNumeratorForComparison {
    diagram_id: usize,
    canonized_numerator: Option<Atom>,
    sample_points: Option<Vec<Vec<(Atom, Atom)>>>,
    sample_evaluations: Option<Vec<Atom>>,
}

impl ProcessedNumeratorForComparison {
    fn from_numerator_symbolic_expression(
        diagram_id: usize,
        bare_graph: &BareGraph,
        numerator: Numerator<SymbolicExpression<Color>>,
        numerator_aware_isomorphism_grouping: &NumeratorAwareGraphGroupingOption,
    ) -> Result<Self, FeynGenError> {
        // println!("----");
        // println!(
        //     "Numerator input for diagram        #{}: {}",
        //     diagram_id,
        //     numerator.get_single_atom().unwrap().0
        // );
        let default_processed_data = ProcessedNumeratorForComparison {
            diagram_id,
            canonized_numerator: None,
            sample_points: None,
            sample_evaluations: None,
        };
        let res = if let Some(group_options) = numerator_aware_isomorphism_grouping.get_options() {
            if group_options.test_canonized_numerator
                || group_options.number_of_numerical_samples > 0
            {
                let mut processed_numerator = numerator.clone();
                let mut lmb_replacements = bare_graph.generate_lmb_replacement_rules(
                    "Q(<i>,x<j>__)",
                    "K(<i>,x<j>__)",
                    "P(<i>,x<j>__)",
                );
                // // Flip the momentum direction of all antiparticle
                // for (i_edge, edge) in bare_graph.edges.iter().enumerate() {
                //     if edge.particle.is_antiparticle() {
                //         'find_rep: for rep in lmb_replacements.iter_mut() {
                //             if rep
                //                 .0
                //                 .pattern_match(
                //                     &Atom::parse(&format!("Q({},x__)", i_edge))
                //                         .unwrap()
                //                         .to_pattern(),
                //                     None,
                //                     None,
                //                 )
                //                 .next()
                //                 .is_some()
                //             {
                //                 rep.1 = (rep.1.clone() * 1).expand();
                //                 break 'find_rep;
                //             }
                //         }
                //     }
                // }
                // Force "final-state" momenta and pol vectors to be identical to external momenta
                for (i_ext, connection) in bare_graph.external_connections.iter().enumerate() {
                    if let (Some(left_edge_pos), Some(right_edge_pos)) = connection {
                        let right_edge = &bare_graph.edges[*right_edge_pos];
                        let left_edge = &bare_graph.edges[*left_edge_pos];
                        let connected_external_id = bare_graph.external_connections.len() + i_ext;
                        for rep in lmb_replacements.iter_mut() {
                            rep.1 = rep.1.replace_all(
                                &Atom::parse(&format!("P({},x__)", connected_external_id))
                                    .unwrap()
                                    .to_pattern(),
                                Atom::parse(&format!("P({},x__)", i_ext))
                                    .unwrap()
                                    .to_pattern(),
                                None,
                                None,
                            );
                        }
                        let left_edge_pol = match left_edge.edge_type {
                            EdgeType::Incoming => left_edge.particle.in_pol_symbol(),
                            EdgeType::Outgoing => left_edge.particle.out_pol_symbol(),
                            _ => unreachable!(),
                        };
                        let right_edge_pol = match right_edge.edge_type {
                            EdgeType::Incoming => right_edge.particle.in_pol_symbol(),
                            EdgeType::Outgoing => right_edge.particle.out_pol_symbol(),
                            _ => unreachable!(),
                        };
                        if let (Some(left_edge_pol), Some(right_edge_pol)) =
                            (left_edge_pol, right_edge_pol)
                        {
                            lmb_replacements.push((
                                Atom::parse(&format!(
                                    "{}({},x__)",
                                    right_edge_pol, connected_external_id
                                ))
                                .unwrap(),
                                Atom::parse(&format!("{}({},x__)", left_edge_pol, i_ext)).unwrap(),
                            ));
                            // lmb_replacements.push((
                            //     Atom::parse(&format!(
                            //         "{}({},xA__)*{}({},xB__)",
                            //         left_edge_pol, i_ext, right_edge_pol, connected_external_id,
                            //     ))
                            //     .unwrap(),
                            //     Atom::parse("Metric(xA__,xB__)").unwrap(),
                            // ));
                        }
                    }
                }
                //lmb_replacements.push((Atom::parse("MB").unwrap(), Atom::Zero));
                // let test: Atom = Atom::parse("-1/9*𝑖*ee^2*G^2*(MB*id(bis(4,2),bis(4,5))+γ(mink(4,25),bis(4,2),bis(4,5))*Q(2,mink(4,25)))*(MB*id(bis(4,3),bis(4,8))-γ(mink(4,27),bis(4,8),bis(4,3))*Q(4,mink(4,27)))*(MB*id(bis(4,4),bis(4,7))+γ(mink(4,28),bis(4,4),bis(4,7))*Q(5,mink(4,28)))*(MB*id(bis(4,6),bis(4,9))+γ(mink(4,29),bis(4,6),bis(4,9))*Q(6,mink(4,29)))*Metric(mink(4,0),mink(4,1))*Metric(mink(4,2),mink(4,3))*id(mink(4,0),mink(4,5))*id(mink(4,1),mink(4,4))*γ(mink(4,2),bis(4,3),bis(4,2))*γ(mink(4,3),bis(4,5),bis(4,4))*γ(mink(4,4),bis(4,7),bis(4,6))*γ(mink(4,5),bis(4,9),bis(4,8))").unwrap();
                // let reps = lmb_replacements
                //     .iter()
                //     .map(|(l, r)| Replacement::new(l.to_pattern(), r.to_pattern()))
                //     .collect::<Vec<_>>();
                // let test2 = test.replace_all_multiple(&reps);
                // println!("BEFORE: {}", test);
                // println!("AFTER: {}", test2);
                // println!(
                //     "BEFORE: {}",
                //     processed_numerator.get_single_atom().unwrap().0
                // );
                // println!(
                //     "REPLACEMENTS:\n{}",
                //     lmb_replacements
                //         .iter()
                //         .map(|(a, b)| format!("{} -> {}", a, b))
                //         .collect::<Vec<_>>()
                //         .join("\n")
                // );
                processed_numerator = processed_numerator.apply_reps(
                    lmb_replacements
                        .iter()
                        .map(|(a, b)| (a.as_view(), b.as_view()))
                        .collect::<Vec<_>>(),
                );
                // println!(
                //     "AFTER: {}",
                //     processed_numerator.get_single_atom().unwrap().0
                // );
                // processed_numerator = processed_numerator.apply_reps2(&lmb_replacements);
                // println!(
                //     "processed_numerator A:\n{}",
                //     processed_numerator.state.colorless
                // );
                // println!(
                //     "processed_numerator:\n{}",
                //     processed_numerator.get_single_atom().unwrap().0
                // );
                let canonized_numerator = if group_options.test_canonized_numerator {
                    Some(
                        processed_numerator
                            .canonize()
                            .unwrap()
                            .get_single_atom()
                            .unwrap()
                            .0,
                    )
                } else {
                    None
                };

                let (sample_points, sample_evaluations) =
                    if group_options.number_of_numerical_samples > 0 {
                        let parsed_numerator = processed_numerator.parse();

                        // println!(
                        //     "NUMERATOR DOT:\n{}",
                        //     parsed_numerator.state.net.rich_graph().dot()
                        // );
                        // panic!("stop");

                        // println!("----");

                        if let Some(grouping_options) =
                            numerator_aware_isomorphism_grouping.get_options()
                        {
                            let mut prime_iterator = PrimeIteratorU64::new(1);
                            prime_iterator.nth(grouping_options.numerical_sample_seed as usize);
                            let sample_points = (0..grouping_options.number_of_numerical_samples)
                                .map(|_| {
                                    parsed_numerator.random_concretize_reps(
                                        Some(&mut prime_iterator),
                                        grouping_options
                                            .fully_numerical_substitution_when_comparing_numerators,
                                    )
                                })
                                .collect::<Vec<_>>();
                            let sample_evaluations = sample_points
                                .iter()
                                .map(|sp| {
                                    let reps = sp
                                        .iter()
                                        .map(|(a, b)| (a.as_view(), b.as_view()))
                                        .collect::<Vec<_>>();
                                    parsed_numerator.evaluate_with_replacements(
                                        reps,
                                        grouping_options
                                            .fully_numerical_substitution_when_comparing_numerators,
                                    )
                                })
                                .collect::<Result<Vec<_>, _>>()?;
                            (Some(sample_points), Some(sample_evaluations))
                        } else {
                            (None, None)
                        }
                    } else {
                        (None, None)
                    };

                ProcessedNumeratorForComparison {
                    diagram_id,
                    canonized_numerator,
                    sample_points,
                    sample_evaluations,
                }
            } else {
                default_processed_data
            }
        } else {
            default_processed_data
        };
        Ok(res)
    }
}
