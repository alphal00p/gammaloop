use indicatif::ParallelProgressIterator;
use indicatif::ProgressBar;
use indicatif::ProgressStyle;
use rayon::iter::IntoParallelIterator;
use rayon::iter::ParallelIterator;
use std::collections::hash_map::Entry;
use std::collections::HashSet;
use std::str::FromStr;
use std::sync::Arc;
use std::time::Duration;

use ahash::AHashMap;
use ahash::AHashSet;
use ahash::HashMap;
use colored::Colorize;
use itertools::zip;
use log::debug;
use log::info;
use smartstring::{LazyCompact, SmartString};
use spenso::iterators::IteratableTensor;
use spenso::parametric::atomcore::TensorAtomMaps;
use spenso::parametric::ParamTensor;
use spenso::upgrading_arithmetic::FallibleAdd;
use spenso::upgrading_arithmetic::FallibleSub;
use symbolica::atom::AtomCore;
use symbolica::domains::rational::Rational;
use symbolica::id::Pattern;

use super::NumeratorAwareGraphGroupingOption;
use super::SelfEnergyFilterOptions;
use super::SnailFilterOptions;
use super::TadpolesFilterOptions;
use super::{FeynGenError, FeynGenOptions};
use crate::model::ColorStructure;
use crate::model::Particle;
use crate::model::VertexRule;
use crate::momentum::SignOrZero;
use crate::numerator::Color;
use crate::numerator::ContractionSettings;
use crate::numerator::GlobalPrefactor;
use crate::numerator::Numerator;
use crate::numerator::SymbolicExpression;
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
    pub fn assign_node_colors<'a>(
        model: &Model,
        graph: &SymbolicaGraph<NodeColorWithoutVertexRule, &'a str>,
        node_colors: &HashMap<
            Vec<(Option<bool>, SmartString<LazyCompact>)>,
            Vec<SmartString<LazyCompact>>,
        >,
    ) -> Result<Vec<(SymbolicaGraph<NodeColorWithVertexRule, &'a str>, usize)>, FeynGenError> {
        // println!("graph = {}", graph.to_dot());
        let mut colored_nodes: Vec<Vec<NodeColorWithVertexRule>> = vec![vec![]];
        let edges = graph.edges();
        for (i_n, node) in graph.nodes().iter().enumerate() {
            let mut node_edges = vec![];
            for e in node.edges.iter() {
                let orientation = if edges[*e].directed {
                    Some(edges[*e].vertices.0 == i_n)
                } else {
                    None
                };
                node_edges.push((orientation, edges[*e].data.into()));
                // self-loop edge must be counted twice
                if edges[*e].vertices.0 == edges[*e].vertices.1 {
                    node_edges.push((orientation.map(|o| !o), edges[*e].data.into()));
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
                _ = g.add_edge(edge.vertices.0, edge.vertices.1, edge.directed, edge.data);
            }
            colored_graphs.push(g);
        }
        Ok(FeynGen::group_isomorphic_graphs(&colored_graphs))
    }

    #[allow(clippy::type_complexity)]
    pub fn group_isomorphic_graphs<
        'a,
        NodeColor: Clone + PartialEq + Eq + PartialOrd + Ord + std::hash::Hash,
    >(
        graphs: &[SymbolicaGraph<NodeColor, &'a str>],
    ) -> Vec<(SymbolicaGraph<NodeColor, &'a str>, usize)> {
        if graphs.len() == 1 {
            return vec![(graphs[0].clone(), 1)];
        }
        let mut iso_buckets: HashMap<SymbolicaGraph<NodeColor, &'a str>, usize> =
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
        graph: &SymbolicaGraph<NodeColorWithoutVertexRule, &str>,
        particles: &[SmartString<LazyCompact>],
    ) -> bool {
        let mut particles_stack = particles
            .iter()
            .map(|p_name| p_name.as_str())
            .collect::<Vec<_>>();
        for edge in graph.edges().iter() {
            if let Some(pos) = particles_stack.iter().position(|&x| x == edge.data) {
                particles_stack.swap_remove(pos);
            }
            if particles_stack.is_empty() {
                return true;
            }
        }
        false
    }

    pub fn find_edge_position(
        graph: &SymbolicaGraph<NodeColorWithoutVertexRule, &str>,
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
        graph: &SymbolicaGraph<NodeColorWithoutVertexRule, &str>,
        veto_self_energy: Option<&SelfEnergyFilterOptions>,
        veto_tadpole: Option<&TadpolesFilterOptions>,
        veto_snails: Option<&SnailFilterOptions>,
    ) -> bool {
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
        let graph_edges: &[symbolica::graph::Edge<&str>] = graph.edges();

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
                    model.get_particle(&SmartString::<LazyCompact>::from(e.data));
            } else if graph_nodes[e.vertices.1].data.external_tag != 0 {
                external_partices[(graph_nodes[e.vertices.1].data.external_tag - 1) as usize] =
                    model.get_particle(&SmartString::<LazyCompact>::from(e.data));
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
            let attachment_particle =
                model.get_particle(&SmartString::<LazyCompact>::from(attachment_edge.data));

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
        graph: &SymbolicaGraph<NodeColor, &str>,
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
            |_, d| d.map(|d| model.get_particle(&SmartString::from_str(d).unwrap())),
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
        graph: &SymbolicaGraph<NodeColor, &str>,
        n_unresolved: usize,
        unresolved_type: &AHashSet<Arc<Particle>>,
        particles: &[SmartString<LazyCompact>],
    ) -> bool {
        let n_initial_states = self.options.initial_pdgs.len();
        // Quick filter to check if the diagram can possibly contain a cut with the right particles and splitting the graph in two parts,
        // as expected for a final-state cut of a forward scattering graph.

        // TODO: replace this with s_and_t_cut from hedge, so as to have proper handling of particle vs anti-particles.
        //let he_graph = HedgeGraph::from(graph.clone());
        //he_graph.all_s_t_cuts(s, t, regions)

        // for now normalize it all to particles
        let mut cut_map: std::collections::HashMap<
            SmartString<LazyCompact>,
            std::vec::Vec<usize>,
            ahash::RandomState,
        > = HashMap::default();
        for p in particles.iter().collect::<HashSet<_>>() {
            let particle = model.get_particle(p);
            if particle.is_antiparticle() {
                cut_map.insert(particle.get_anti_particle(model).name.clone(), vec![]);
            } else {
                cut_map.insert(p.clone(), vec![]);
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
            let particle = model.get_particle(&edge.data.into());
            let e = if particle.is_antiparticle() {
                cut_map.get_mut(&particle.get_anti_particle(model).name)
            } else {
                cut_map.get_mut(&particle.name)
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
                        cut_map.get(particle.as_str()).unwrap().iter()
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
            graph: &SymbolicaGraph<NodeColor, &str>,
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
                for subcut in (1..=(particles.len() - 1))
                    .flat_map(|offset| cut.iter().combinations(particles.len() - offset))
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

    fn group_isomorphic_graphs_after_node_color_change<'a, NC>(
        graphs: &HashMap<SymbolicaGraph<NC, &'a str>, Atom>,
        node_colors_to_change: &HashMap<i32, i32>,
    ) -> HashMap<SymbolicaGraph<NC, &'a str>, Atom>
    where
        NC: NodeColorFunctions + Clone + PartialOrd + Ord + Eq + std::hash::Hash,
    {
        #[allow(clippy::type_complexity)]
        let mut iso_buckets: HashMap<
            SymbolicaGraph<NC, &'a str>,
            (usize, (SymbolicaGraph<NC, &'a str>, Vec<usize>, Atom)),
        > = HashMap::default();

        for (g, symmetry_factor) in graphs.iter() {
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
            if let Some((count, (other_g, external_ordering, sym_fact))) =
                iso_buckets.get_mut(&canonized_g.graph)
            {
                // Make sure to pick a representative with a canonical external ordering so that additional flavour grouping afterwards is possible
                if ext_ordering < *external_ordering {
                    *external_ordering = ext_ordering;
                    *other_g = g.clone();
                }
                *count += 1;

                assert!(sym_fact == symmetry_factor);
            } else {
                iso_buckets.insert(
                    canonized_g.graph,
                    (1, (g.clone(), ext_ordering, symmetry_factor.clone())),
                );
            }
        }

        iso_buckets
            .iter()
            .map(
                |(_canonized_g, (count, (g, _external_ordering, symmetry_factor)))| {
                    (g.clone(), symmetry_factor * Atom::new_num(*count as i64))
                },
            )
            .collect()
    }

    /// This function canonizes the edge and vertex ordering of a graph based on the skeletton graph with only propagator mass as edge color.
    /// This is useful to then allow for further grouping of isomorphic graphs, incl numerator.
    #[allow(clippy::type_complexity)]
    pub fn canonized_edge_and_vertex_ordering<'a>(
        &self,
        model: &Model,
        input_graph: &SymbolicaGraph<NodeColorWithVertexRule, &'a str>,
    ) -> (
        SymbolicaGraph<NodeColorWithoutVertexRule, std::string::String>,
        SymbolicaGraph<NodeColorWithVertexRule, &'a str>,
    ) {
        // Make sure to canonize the edge ordering based on the skeletton graph with only propagator mass as edge color
        let mut skeletton_graph = SymbolicaGraph::new();
        for node in input_graph.nodes() {
            skeletton_graph.add_node(NodeColorWithoutVertexRule {
                external_tag: node.data.external_tag,
            });
        }
        for edge in input_graph.edges() {
            skeletton_graph
                .add_edge(
                    edge.vertices.0,
                    edge.vertices.1,
                    false,
                    model.get_particle(&edge.data.into()).mass.to_string(),
                )
                .unwrap();
        }
        let canonized_skeletton = skeletton_graph.canonize();

        let mut reordered_nodes: Vec<usize> = vec![0; input_graph.nodes().len()];
        for (input_graph_node_position, node_order) in
            canonized_skeletton.vertex_map.iter().enumerate()
        {
            reordered_nodes[*node_order] = input_graph_node_position;
        }
        let mut reordered_edges = input_graph
            .edges()
            .iter()
            .map(|e| {
                let unoriented_edge = if canonized_skeletton.vertex_map[e.vertices.0]
                    < canonized_skeletton.vertex_map[e.vertices.1]
                {
                    (
                        canonized_skeletton.vertex_map[e.vertices.0],
                        canonized_skeletton.vertex_map[e.vertices.1],
                    )
                } else {
                    (
                        canonized_skeletton.vertex_map[e.vertices.1],
                        canonized_skeletton.vertex_map[e.vertices.0],
                    )
                };
                (
                    e,
                    // This key will serve to give a unique ordering of the edges
                    (
                        unoriented_edge.0,
                        unoriented_edge.1,
                        model.get_particle(&e.data.into()).mass.to_string(),
                    ),
                )
            })
            .collect::<Vec<_>>();
        reordered_edges.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap());

        // Sort nodes according to the canonized skeletton graph
        let mut sorted_g = SymbolicaGraph::new();
        for node_order in reordered_nodes {
            sorted_g.add_node(input_graph.nodes()[node_order].data.clone());
        }
        for (e, _sorting_key) in reordered_edges.iter() {
            sorted_g
                .add_edge(
                    canonized_skeletton.vertex_map[e.vertices.0],
                    canonized_skeletton.vertex_map[e.vertices.1],
                    e.directed,
                    e.data,
                )
                .unwrap();
        }

        (canonized_skeletton.graph, sorted_g)
    }

    #[allow(clippy::too_many_arguments)]
    pub fn generate(
        &self,
        model: &Model,
        numerator_aware_isomorphism_grouping: NumeratorAwareGraphGroupingOption,
        filter_self_loop: bool,
        graph_prefix: String,
        selected_graphs: Option<Vec<String>>,
        vetoed_graphs: Option<Vec<String>>,
        loop_momentum_bases: Option<HashMap<String, Vec<String>>>,
    ) -> Result<Vec<BareGraph>, FeynGenError> {
        let progress_bar_style = ProgressStyle::with_template(
            "[{elapsed_precise} | ETA: {eta_precise}] {bar:40.cyan/blue} {pos:>7}/{len:7} ({percent}%) {msg}",
        )
        .unwrap();

        debug!(
            "Generating Feynman diagrams for model {} and process:\n{}",
            model.name, self.options
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
                        (None, p.name.clone())
                    } else if p.is_antiparticle() {
                        (Some(false), p.get_anti_particle(model).name.clone())
                    } else {
                        (Some(true), p.name.clone())
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
                                    (None, p.name.clone())
                                } else if p.is_antiparticle() {
                                    (Some(true), p.get_anti_particle(model).name.clone())
                                } else {
                                    (Some(false), p.name.clone())
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
                            (None, p.name.clone())
                        } else if p.is_antiparticle() {
                            (Some(true), p.get_anti_particle(model).name.clone())
                        } else {
                            (Some(false), p.name.clone())
                        },
                    ));
                }
            }
        }

        // debug!("external_edges = {:?}", external_edges);
        // debug!("vertex_signatures = {:?}", vertex_signatures);
        let external_edges_for_generation = external_edges
            .iter()
            .map(|(i, (orientation, name))| (i.clone(), (*orientation, name.as_str())))
            .collect::<Vec<_>>();
        let vertex_signatures_for_generation = vertex_signatures
            .keys()
            .map(|v| {
                v.iter()
                    .map(|(orientation, p)| (*orientation, p.as_str()))
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
        info!("Starting Feynman graph generation with Symbolica...");
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
        info!(
            "{:<95}{}",
            "Number of graphs resulting from Symbolica generation:",
            format!("{}", graphs.len()).green().bold()
        );

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
                        p.get_anti_particle(model).name.clone()
                    } else {
                        p.name.clone()
                    }
                })
                .collect::<Vec<_>>()
        } else {
            vec![]
        };

        if !unoriented_final_state_particles.is_empty() {
            graphs.retain(|g, _symmetry_factor| {
                FeynGen::contains_particles(g, unoriented_final_state_particles.as_slice())
            });
            info!(
                "{:<95}{}",
                "Number of graphs retained after enforcing supergraph particle content:",
                format!("{}", graphs.len()).green()
            );
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
            graphs.retain(|g, _symmetry_factor| {
                !FeynGen::veto_special_topologies(
                    model,
                    g,
                    external_self_energy_filter,
                    tadpole_filter,
                    zero_snails_filter,
                )
            });
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

        info!(
            "{:<95}{}",
            "Number of graphs retained after removal of vetoed topologies:",
            format!("{}", graphs.len()).green()
        );

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
        processed_graphs.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());

        info!(
            "{:<95}{}",
            "Number of graphs after vertex info assignment:",
            format!("{}", processed_graphs.len()).green()
        );

        filters.apply_filters(&mut processed_graphs)?;

        info!(
            "{:<95}{}",
            "Number of graphs after all complete graphs filters are applied:",
            format!("{}", processed_graphs.len()).green()
        );

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
            bar.enable_steady_tick(Duration::from_millis(100));
            processed_graphs.retain(|(g, _symmetry_factor)| {
                bar.inc(1);
                self.contains_cut_fast(
                    model,
                    g,
                    n_unresolved,
                    &unresolved_type,
                    unoriented_final_state_particles.as_slice(),
                )
            });
            bar.finish_and_clear();

            info!(
                "{:<95}{}",
                "Number of graphs after fast Cutkosky cut filter is applied:",
                format!("{}", processed_graphs.len()).green()
            );
        }

        // This secondary cutkosky cut filter is only necessary if some amplitude coupling order constraints where requested
        // because that one cannot be done until the nodes have been colored with the proper choice of vertex rule
        // Also the amplitude loop count restriction would not have been correctly handled when using the short-circuit.
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
            let bar = ProgressBar::new(processed_graphs.len() as u64);
            bar.set_style(progress_bar_style.clone());
            bar.set_message("Applying secondary exact Cutkosky cut filter...");
            bar.enable_steady_tick(Duration::from_millis(100));
            processed_graphs = processed_graphs
                .into_par_iter()
                .filter(|(g, _)| self.contains_cut(model, g, n_unresolved, &unresolved_type))
                .progress_with(bar)
                .collect();

            info!(
                "{:<95}{}",
                "Number of graphs after exact Cutkosky cut filter is applied:",
                format!("{}", processed_graphs.len()).green()
            );
        }

        // Now account for initial state symmetry by further grouping contributions
        let mut node_colors_to_change: HashMap<i32, i32> = HashMap::default();
        if self.options.symmetrize_initial_states {
            for initial_color in 1..=self.options.initial_pdgs.len() {
                node_colors_to_change.insert(
                    initial_color as i32,
                    if self.options.symmetrize_left_right_states {
                        0
                    } else {
                        -1
                    },
                );
            }
            if self.options.generation_type == GenerationType::CrossSection {
                for final_color in
                    self.options.initial_pdgs.len() + 1..=2 * self.options.initial_pdgs.len()
                {
                    node_colors_to_change.insert(
                        final_color as i32,
                        if self.options.symmetrize_left_right_states {
                            0
                        } else {
                            -2
                        },
                    );
                }
            }
        }
        if self.options.generation_type == GenerationType::Amplitude
            && self.options.symmetrize_final_states
        {
            for final_color in self.options.initial_pdgs.len() + 1
                ..=self.options.initial_pdgs.len() + self.options.final_pdgs.len()
            {
                node_colors_to_change.insert(
                    final_color as i32,
                    if self.options.symmetrize_left_right_states {
                        0
                    } else {
                        -2
                    },
                );
            }
        }
        if !node_colors_to_change.is_empty() {
            processed_graphs = FeynGen::group_isomorphic_graphs_after_node_color_change(
                &processed_graphs
                    .iter()
                    .map(|(g, m)| (g.clone(), Atom::parse(m).unwrap()))
                    .collect::<HashMap<_, _>>(),
                &node_colors_to_change,
            )
            .iter()
            .map(|(g, m)| (g.clone(), m.to_canonical_string()))
            .collect::<Vec<_>>();
        }

        info!(
            "{:<95}{}",
            "Number of graphs after symmetrization of external states:",
            format!("{}", processed_graphs.len()).green()
        );

        #[allow(clippy::type_complexity)]
        let mut pooled_bare_graphs: HashMap<
            SymbolicaGraph<NodeColorWithoutVertexRule, std::string::String>,
            Vec<(usize, Option<ProcessedNumeratorForComparison>, BareGraph)>,
        > = HashMap::default();

        let bar = ProgressBar::new(processed_graphs.len() as u64);
        bar.set_style(progress_bar_style.clone());
        bar.set_message("Final numerator-aware processing of remaining graphs...");
        bar.enable_steady_tick(Duration::from_millis(100));
        for (i_g, (g, symmetry_factor)) in processed_graphs.iter().enumerate() {
            bar.inc(1);
            // println!(
            //     "Computing canonical representation for graph #{}:\nedges: {:?}\n nodes: {:?}",
            //     i_g,
            //     g.edges()
            //         .iter()
            //         .enumerate()
            //         .map(|(i_e, e)| format!("{}|{}", i_e, e.data))
            //         .collect::<Vec<_>>(),
            //     g.nodes()
            //         .iter()
            //         .enumerate()
            //         .map(|(i_n, n)| format!(
            //             "{}|{}|{}|edges:{:?}",
            //             i_n, n.data.0, n.data.1, n.edges
            //         ))
            //         .collect::<Vec<_>>()
            // );
            let (canonical_repr, sorted_g) = self.canonized_edge_and_vertex_ordering(model, g);
            // println!(
            //     "Canonical representation for graph #{}:\nedges: {:?}\n nodes: {:?}",
            //     i_g,
            //     sorted_g
            //         .edges()
            //         .iter()
            //         .enumerate()
            //         .map(|(i_e, e)| format!("{}|{}", i_e, e.data))
            //         .collect::<Vec<_>>(),
            //     sorted_g
            //         .nodes()
            //         .iter()
            //         .enumerate()
            //         .map(|(i_n, n)| format!(
            //             "{}|{}|{}|edges:{:?}",
            //             i_n, n.data.0, n.data.1, n.edges
            //         ))
            //         .collect::<Vec<_>>()
            // );
            // println!("Canonical dot for graph #{}:", canonical_repr.to_dot());
            let graph_name = format!("{}{}", graph_prefix, i_g);
            if let Some(selected_graphs) = &selected_graphs {
                if !selected_graphs.contains(&graph_name) {
                    continue;
                }
            }
            if let Some(vetoed_graphs) = &vetoed_graphs {
                if vetoed_graphs.contains(&graph_name) {
                    continue;
                }
            }
            let forced_lmb = if let Some(lmbs) = loop_momentum_bases.as_ref() {
                lmbs.get(&graph_name)
                    .map(|lmb| lmb.iter().map(SmartString::<LazyCompact>::from).collect())
            } else {
                None
            };
            let bare_graph = BareGraph::from_symbolica_graph(
                model,
                graph_name,
                &sorted_g,
                symmetry_factor.clone(),
                external_connections.clone(),
                forced_lmb,
            )?;
            // When disabling numerator-aware graph isomorphism, each graph is added separately
            if numerator_aware_isomorphism_grouping == NumeratorAwareGraphGroupingOption::NoGrouping
            {
                match pooled_bare_graphs.entry(canonical_repr) {
                    Entry::Vacant(entry) => {
                        entry.insert(vec![(i_g, None, bare_graph)]);
                    }
                    Entry::Occupied(mut entry) => {
                        entry.get_mut().push((i_g, None, bare_graph));
                    }
                }
            } else {
                let numerator =
                    Numerator::default().from_graph(&bare_graph, &GlobalPrefactor::default());
                let numerator_color_simplified = numerator.color_simplify();
                if numerator_color_simplified
                    .get_single_atom()
                    .unwrap()
                    .0
                    .is_zero()
                {
                    continue;
                }
                if numerator_aware_isomorphism_grouping
                    == NumeratorAwareGraphGroupingOption::OnlyDetectZeroes
                {
                    match pooled_bare_graphs.entry(canonical_repr) {
                        Entry::Vacant(entry) => {
                            entry.insert(vec![(i_g, None, bare_graph)]);
                        }
                        Entry::Occupied(mut entry) => {
                            entry.get_mut().push((i_g, None, bare_graph));
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
                            &numerator_aware_isomorphism_grouping,
                            i_g,
                            numerator_color_simplified,
                        ),
                    );
                    // let numerator = Some(
                    //     numerator_color_simplified
                    //         .parse()
                    //         .contract(ContractionSettings::<Rational>::Normal)
                    //         .unwrap()
                    //         .state
                    //         .tensor
                    //         // TODO: A future implementation of "compare_tensors" should be able to avoid the expansion below.
                    //         .map_data(|d| d.expand()),
                    // );

                    match pooled_bare_graphs.entry(canonical_repr) {
                        Entry::Vacant(entry) => {
                            entry.insert(vec![(i_g, numerator_data, bare_graph)]);
                        }
                        Entry::Occupied(mut entry) => {
                            let mut found_match = false;
                            for (_graph_id, other_numerator, other_graph) in entry.get_mut() {
                                if let Some(ratio) = FeynGen::compare_numerator_tensors(
                                    &numerator_aware_isomorphism_grouping,
                                    numerator_data.as_ref().unwrap(),
                                    other_numerator.as_ref().unwrap(),
                                ) {
                                    found_match = true;
                                    other_graph.overall_factor =
                                        (Atom::parse(other_graph.overall_factor.as_str()).unwrap()
                                            + ratio)
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
        bar.finish_and_clear();

        let mut bare_graphs = pooled_bare_graphs
            .values()
            .flatten()
            .filter_map(|(graph_id, _numerator, graph)| {
                if Atom::parse(&graph.overall_factor)
                    .unwrap()
                    .expand()
                    .is_zero()
                {
                    None
                } else {
                    Some((*graph_id, graph))
                }
            })
            .collect::<Vec<_>>();
        bare_graphs.sort_by(|a, b| (a.0).cmp(&b.0));

        info!(
            "{:<95}{} ({} isomorphically unique graphs when ignoring numerators)",
            format!(
                "Number of graphs after numerator-aware grouping with strategy '{}':",
                numerator_aware_isomorphism_grouping
            ),
            format!("{}", bare_graphs.len()).green().bold(),
            pooled_bare_graphs.len()
        );

        Ok(bare_graphs
            .iter()
            .cloned()
            .map(|(_graph_id, graph)| graph.clone())
            .collect::<Vec<_>>())
    }

    fn compare_numerator_tensors(
        numerator_aware_isomorphism_grouping: &NumeratorAwareGraphGroupingOption,
        numerator_a: &ProcessedNumeratorForComparison,
        numerator_b: &ProcessedNumeratorForComparison,
    ) -> Option<Atom> {
        match numerator_aware_isomorphism_grouping {
            NumeratorAwareGraphGroupingOption::NoGrouping
            | NumeratorAwareGraphGroupingOption::OnlyDetectZeroes => None,
            NumeratorAwareGraphGroupingOption::GroupIdenticalGraphUpToScalarRescaling => {
                let ratios = zip(
                    numerator_a.tensor.iter_flat(),
                    numerator_b.tensor.iter_flat(),
                )
                .map(|((_idx_a, a), (_idx_b, b))| (a / b))
                .collect::<HashSet<_>>();
                // for (a_idx, a) in numerator_a.tensor.iter_flat() {
                //     println!("Numerator A #{} : {}", a_idx, a);
                // }
                // for (b_idx, b) in numerator_b.tensor.iter_flat() {
                //     println!("Numerator B #{} : {}", b_idx, b);
                // }
                // println!(
                //     "ratios: {:?}",
                //     ratios
                //         .iter()
                //         .map(|av| av.to_canonical_string())
                //         .collect::<Vec<_>>()
                //         .join(",")
                // );
                if ratios.len() > 1 {
                    None
                } else {
                    let ratio = ratios.iter().next().unwrap().to_owned();
                    // Make sure the ratio only consists of couplings
                    for pattern in ["cind(arg_)"] {
                        if ratio
                            .pattern_match(&Pattern::parse(pattern).unwrap(), None, None)
                            .next()
                            .is_some()
                        {
                            return None;
                        }
                    }
                    debug!(
                        "Combining graph #{} with #{}, with ratio #{}/#{}={}",
                        numerator_b.diagram_id,
                        numerator_a.diagram_id,
                        numerator_a.diagram_id,
                        numerator_b.diagram_id,
                        ratio
                    );
                    Some(ratio)
                }
            }
            NumeratorAwareGraphGroupingOption::GroupIdenticalGraphUpToSign => {
                // for (a_idx, a) in numerator_a.tensor.iter_flat() {
                //     println!("Numerator A #{} : {}", a_idx, a);
                // }
                // for (b_idx, b) in numerator_b.tensor.iter_flat() {
                //     println!("Numerator B #{} : {}", b_idx, b);
                // }
                let res = if let Some(diff) = numerator_a.tensor.sub_fallible(&numerator_b.tensor) {
                    if diff
                        .iter_flat()
                        .all(|(_idx, d)| d.expand_num().collect_num().is_zero())
                    {
                        Some(1)
                    } else {
                        let sum = numerator_a
                            .tensor
                            .add_fallible(&numerator_b.tensor)
                            .unwrap();
                        // for (sum_idx, sum) in sum.iter_flat() {
                        //     println!("Sum #{}: {}", sum_idx, sum);
                        // }
                        if sum
                            .iter_flat()
                            .all(|(_idx, d)| d.expand_num().collect_num().is_zero())
                        {
                            Some(-1)
                        } else {
                            None
                        }
                    }
                } else {
                    None
                };
                if let Some(r) = res {
                    debug!(
                        "Combining identical graph{} #{} with #{}",
                        if r == -1 {
                            " (up to a relative sign)"
                        } else {
                            ""
                        },
                        numerator_b.diagram_id,
                        numerator_a.diagram_id
                    );
                    Some(Atom::new_num(r))
                } else {
                    None
                }
            }
        }
    }
}

struct ProcessedNumeratorForComparison {
    diagram_id: usize,
    tensor: ParamTensor,
}

impl ProcessedNumeratorForComparison {
    fn from_numerator_symbolic_expression(
        numerator_aware_isomorphism_grouping: &NumeratorAwareGraphGroupingOption,
        diagram_id: usize,
        numerator: Numerator<SymbolicExpression<Color>>,
    ) -> Self {
        let mut tensor = numerator
            .parse()
            .contract::<Rational>(ContractionSettings::<Rational>::Normal)
            .unwrap()
            .state
            .tensor;
        match numerator_aware_isomorphism_grouping {
            NumeratorAwareGraphGroupingOption::NoGrouping
            | NumeratorAwareGraphGroupingOption::OnlyDetectZeroes => {}
            NumeratorAwareGraphGroupingOption::GroupIdenticalGraphUpToSign => {
                // d.expand() is much slower and not necessary as we just need to canonicalize the product, sadly it's necessary in complicated cases, i.e.
                // generate h > a a | h a b [{{3}}] --symmetrize_left_right_states -num_grouping group_identical_graphs_up_to_sign
                // should yield only one graph!
                // TODO: FIX so that d.expand_num().collect_num() works!
                tensor = tensor.expand_num().collect_num()
                // tensor = tensor.expand()

                // tensor = tensor.map_data(|d| d.expand_num().collect_num())
            }
            NumeratorAwareGraphGroupingOption::GroupIdenticalGraphUpToScalarRescaling => {
                // TOFIX: this does not work in more complicated cases. Need to use canonicalized tensors instead of spenso-expanded ones.
                tensor = tensor.expand_num().collect_num()
            }
        }
        ProcessedNumeratorForComparison { diagram_id, tensor }
    }
}
