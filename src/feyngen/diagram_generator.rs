use std::collections::HashSet;
use std::sync::Arc;

use ahash::HashMap;
use log::debug;
use smartstring::{LazyCompact, SmartString};

use super::SelfEnergyFilterOptions;
use super::SnailFilterOptions;
use super::TadpolesFilterOptions;
use super::{FeynGenError, FeynGenOptions};
use crate::model::Particle;
use crate::numerator::Numerator;
use crate::{
    feyngen::{FeynGenFilter, GenerationType},
    graph::BareGraph,
    model::Model,
};
use itertools::Itertools;
use symbolica::{atom::Atom, domains::integer::Integer, graph::Graph as SymbolicaGraph};

pub struct FeynGen {
    pub options: FeynGenOptions,
}

impl FeynGen {
    pub fn new(options: FeynGenOptions) -> Self {
        Self { options }
    }

    #[allow(clippy::type_complexity)]
    pub fn assign_node_colors<'a>(
        graph: &SymbolicaGraph<usize, &'a str>,
        node_colors: &HashMap<Vec<SmartString<LazyCompact>>, Vec<SmartString<LazyCompact>>>,
    ) -> Result<
        Vec<(
            SymbolicaGraph<(usize, SmartString<LazyCompact>), &'a str>,
            usize,
        )>,
        FeynGenError,
    > {
        // println!("graph = {}", graph.to_dot());
        let mut colored_nodes: Vec<Vec<(usize, SmartString<LazyCompact>)>> = vec![vec![]];
        let edges = graph.edges();
        for node in graph.nodes().iter() {
            let mut node_edges = vec![];
            for e in node.edges.iter() {
                node_edges.push(edges[*e].data.into());
                // self-loop edge must be counted twice
                if edges[*e].vertices.0 == edges[*e].vertices.1 {
                    node_edges.push(edges[*e].data.into());
                }
            }
            node_edges.sort();
            let colors = if node_edges.len() == 1 {
                &vec!["external".into()]
            } else if let Some(cs) = node_colors.get(&node_edges) {
                cs
            } else {
                return Err(FeynGenError::GenericError(format!(
                    "Could not find node colors for node edges {:?}",
                    node_edges
                )));
            };
            let mut new_colored_nodes: Vec<Vec<(usize, SmartString<LazyCompact>)>> = vec![];
            for current_colors in colored_nodes.iter_mut() {
                for color in colors {
                    current_colors.push((node.data, color.clone()));
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
    pub fn group_isomorphic_graphs<'a>(
        graphs: &[SymbolicaGraph<(usize, SmartString<LazyCompact>), &'a str>],
    ) -> Vec<(
        SymbolicaGraph<(usize, SmartString<LazyCompact>), &'a str>,
        usize,
    )> {
        if graphs.len() == 1 {
            return vec![(graphs[0].clone(), 1)];
        }
        let mut iso_buckets: HashMap<
            SymbolicaGraph<(usize, SmartString<LazyCompact>), &'a str>,
            usize,
        > = HashMap::default();
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
        graph: &SymbolicaGraph<usize, &str>,
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
        graph: &SymbolicaGraph<usize, &str>,
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
        graph: &SymbolicaGraph<usize, &str>,
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
        let graph_nodes: &[symbolica::graph::Node<usize>] = graph.nodes();
        let graph_edges: &[symbolica::graph::Edge<&str>] = graph.edges();

        let max_external = graph
            .nodes()
            .iter()
            .filter(|n| n.data > 0)
            .map(|n| n.data)
            .max()
            .unwrap_or(0);
        if max_external == 0 {
            // Do not implement any veto for vacuum graphs
            return false;
        }
        let max_external_node_position = graph
            .nodes()
            .iter()
            .position(|n| n.data == max_external)
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
            if graph_nodes[e.vertices.0].data != 0 {
                external_partices[graph_nodes[e.vertices.0].data - 1] =
                    model.get_particle(&SmartString::<LazyCompact>::from(e.data));
            } else if graph_nodes[e.vertices.1].data != 0 {
                external_partices[graph_nodes[e.vertices.1].data - 1] =
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
            if (node.edges.len() != 1) || node.data == max_external_node_position {
                continue;
            }
            let external_index = node.data;
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
                            if spanning_tree.nodes.iter().any(|n| {
                                !n.external && n.chain_id.is_none() && n.back_edges.is_empty()
                            }) {
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

    pub fn contains_cut(
        graph: &SymbolicaGraph<usize, &str>,
        n_initial_states: usize,
        particles: &[SmartString<LazyCompact>],
    ) -> bool {
        // Quick filter to check if the diagram can possibly contain a cut with the right particles and splitting the graph in two parts,
        // as expected for a final-state cut of a forward scattering graph.

        let mut cut_map: std::collections::HashMap<&str, std::vec::Vec<usize>, ahash::RandomState> =
            HashMap::default();
        for p in particles {
            cut_map.insert(p.as_str(), vec![]);
        }
        for (i_e, edge) in graph.edges().iter().enumerate() {
            if let Some(cut_entry) = cut_map.get_mut(&edge.data) {
                cut_entry.push(i_e);
            }
        }

        let cut_options_values: Vec<_> = cut_map.values().collect();

        // Generate unique, unordered combinations
        let mut unique_combinations: std::collections::HashSet<Vec<usize>, ahash::RandomState> =
            HashSet::default();
        for combination in cut_options_values
            .iter()
            .map(|v| v.iter())
            .multi_cartesian_product()
        {
            // Sort the combination and insert into HashSet
            let mut sorted_combination: Vec<_> = combination.into_iter().cloned().collect();
            sorted_combination.sort_unstable();
            unique_combinations.insert(sorted_combination);
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
            .map(|leg_id| graph.nodes().iter().position(|n| n.data == leg_id).unwrap())
            .collect::<Vec<_>>();
        let right_initial_state_positions = (n_initial_states + 1..=2 * n_initial_states)
            .map(|leg_id| graph.nodes().iter().position(|n| n.data == leg_id).unwrap())
            .collect::<Vec<_>>();

        'cutloop: for cut in cut_options_values {
            for left_node in left_initial_state_positions.iter() {
                for right_node in right_initial_state_positions.iter() {
                    let mut visited: Vec<bool> = vec![false; graph.nodes().len()];
                    let mut stack: Vec<usize> = vec![*left_node];
                    while let Some(node) = stack.pop() {
                        if node == *right_node {
                            continue 'cutloop;
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
            return true;
        }

        false
    }

    fn group_isomorphic_graphs_after_node_color_change<'a>(
        graphs: &HashMap<SymbolicaGraph<usize, &'a str>, Integer>,
        node_colors_to_change: &HashMap<usize, usize>,
    ) -> HashMap<SymbolicaGraph<usize, &'a str>, Integer> {
        #[allow(clippy::type_complexity)]
        let mut iso_buckets: HashMap<
            SymbolicaGraph<usize, &'a str>,
            (usize, (SymbolicaGraph<usize, &'a str>, Integer)),
        > = HashMap::default();

        for (g, symmetry_factor) in graphs.iter() {
            let mut g_node_color_modified = g.clone();
            let mut modifications: Vec<(usize, usize)> = vec![];
            for (i_n, node) in g.nodes().iter().enumerate() {
                for (src_node_color, trgt_node_color) in node_colors_to_change {
                    if node.data == *src_node_color {
                        modifications.push((i_n, *trgt_node_color));
                    }
                }
            }
            for (i_n, new_color) in modifications {
                g_node_color_modified.set_node_data(i_n, new_color);
            }
            let canonized_g = g_node_color_modified.canonize();
            (iso_buckets
                .entry(canonized_g.graph)
                .or_insert_with(|| (1, (g.clone(), symmetry_factor.clone()))))
            .0 += 1;
        }

        iso_buckets
            .iter()
            .map(|(_canonized_g, (count, (g, symmetry_factor)))| {
                (
                    g.clone(),
                    Integer::Natural(symmetry_factor.to_i64().unwrap() * (*count as i64)),
                )
            })
            .collect()
    }

    #[allow(clippy::too_many_arguments)]
    pub fn generate(
        &self,
        model: &Model,
        numerator_aware_isomorphism_grouping: bool,
        filter_self_loop: bool,
        graph_prefix: String,
        selected_graphs: Option<Vec<String>>,
        vetoed_graphs: Option<Vec<String>>,
        loop_momentum_bases: Option<HashMap<String, Vec<String>>>,
    ) -> Result<Vec<BareGraph>, FeynGenError> {
        debug!(
            "Generating Feynman diagrams for model {} and process:\n{}",
            model.name, self.options
        );

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
                if let Some(FeynGenFilter::ParticleVeto(vetoed_particles)) =
                    self.options.filters.get_particle_vetos()
                {
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
                    i_initial + 1,
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
                                self.options.initial_pdgs.len() + i_final + 1,
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
                        i_final,
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
            .map(|(i, (orientation, name))| (*i, (*orientation, name.as_str())))
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
        let mut graphs = SymbolicaGraph::generate(
            external_edges_for_generation.as_slice(),
            vertex_signatures_for_generation.as_slice(),
            None,
            Some(self.options.loop_count_range.1),
            self.options.filters.get_max_bridge(),
            !filter_self_loop,
        );
        // Immediately drop lower loop count contributions
        graphs.retain(|g, _| g.num_loops() >= self.options.loop_count_range.0);

        debug!("Symbolica generated {} graphs", graphs.len());

        if self.options.generation_type == GenerationType::CrossSection {
            let final_state_particles = self
                .options
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
                .collect::<Vec<_>>();
            graphs.retain(|g, _symmetry_factor| {
                FeynGen::contains_particles(g, final_state_particles.as_slice())
            });

            // Make sure that there can be a valid Cutkosky cut in each graph retained
            if !final_state_particles.is_empty() {
                graphs.retain(|g, _symmetry_factor| {
                    FeynGen::contains_cut(
                        g,
                        self.options.initial_pdgs.len(),
                        final_state_particles.as_slice(),
                    )
                });
            }
        }

        let tadpole_filter = self
            .options
            .filters
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
        let external_self_energy_filter = self
            .options
            .filters
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
        let zero_snails_filter = self
            .options
            .filters
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

        // Now account for initial state symmetry by further grouping contributions
        let mut node_colors_to_change: HashMap<usize, usize> = HashMap::default();
        if self.options.symmetrize_initial_states {
            for initial_color in 1..=self.options.initial_pdgs.len() {
                node_colors_to_change.insert(
                    initial_color,
                    if self.options.generation_type == GenerationType::CrossSection
                        && self.options.symmetrize_left_right_states
                    {
                        0
                    } else {
                        1
                    },
                );
            }
            if self.options.generation_type == GenerationType::CrossSection {
                for final_color in
                    self.options.initial_pdgs.len() + 1..=2 * self.options.initial_pdgs.len()
                {
                    node_colors_to_change.insert(
                        final_color,
                        if self.options.symmetrize_left_right_states {
                            0
                        } else {
                            2
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
                node_colors_to_change.insert(final_color, 0);
            }
        }
        if !node_colors_to_change.is_empty() {
            graphs = FeynGen::group_isomorphic_graphs_after_node_color_change(
                &graphs,
                &node_colors_to_change,
            );
        }

        debug!("Number of graphs retained: {}", graphs.len());

        let mut node_colors: HashMap<Vec<SmartString<LazyCompact>>, Vec<SmartString<LazyCompact>>> =
            HashMap::default();
        for (v_legs, v_colors) in vertex_signatures.iter() {
            let mut sorted_ps = v_legs.iter().map(|(_, p)| p.clone()).collect::<Vec<_>>();
            sorted_ps.sort();
            node_colors.insert(sorted_ps, v_colors.clone());
        }

        let mut processed_graphs = vec![];
        for (g, symmetry_factor) in graphs.iter() {
            for (colored_g, multiplicity) in FeynGen::assign_node_colors(g, &node_colors)? {
                processed_graphs.push((
                    colored_g.canonize().graph,
                    (Atom::new_num(multiplicity as i64)
                        / Atom::new_num(symmetry_factor.to_i64().unwrap()))
                    .to_canonical_string(),
                ));
            }
        }
        processed_graphs.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());

        debug!(
            "Number of graphs after node-color dressing: {}",
            processed_graphs.len()
        );

        self.options
            .filters
            .apply_filters(model, &mut processed_graphs)?;

        debug!(
            "Number of graphs after all filters are applied: {}",
            processed_graphs.len()
        );

        let mut bare_graphs: Vec<BareGraph> = vec![];

        for (i_g, (g, symmetry_factor)) in processed_graphs.iter().enumerate() {
            // Now generate multiple copies of the graph for each vertex and
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
                g,
                symmetry_factor.clone(),
                external_connections.clone(),
                forced_lmb,
            )?;
            bare_graphs.push(bare_graph);
        }

        debug!(
            "Total number of graphs generated by gammaLoop: {}",
            bare_graphs.len()
        );
        let previous_length = bare_graphs.len();
        if numerator_aware_isomorphism_grouping {
            bare_graphs.retain(|g| {
                let numerator = Numerator::default().from_graph(g, None);
                let numerator_color_simplified =
                    numerator.color_simplify().get_single_atom().unwrap().0;
                // println!(
                //     "numerator_color_simplified of graph {} = {}",
                //     g.name, numerator_color_simplified
                // );
                !numerator_color_simplified.is_zero()
            });
            debug!(
                "Total number of graphs remaining after graph isomorphism check including numerator: {} / {}",
                bare_graphs.len(), previous_length
            );
        }

        Ok(bare_graphs)
    }
}
