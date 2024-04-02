use ahash::HashSet;
use color_eyre::Report;
use log::debug;
use num::traits::Zero;
use serde::{Deserialize, Serialize};

use crate::{
    graph::{EdgeType, Graph},
    utils::FloatLike,
};

/// Dimensionality of space
pub const D: usize = 3; // we are always going to do 3-d representations, so I hardcode this.

/// Simplified version of the Graph struct, used to generate the tropical subgraph table
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TropicalGraph {
    dod: f64,
    is_edge_massive: Vec<bool>,
    pub topology: Vec<TropicalEdge>,
    num_massive_edges: usize,
    external_vertices: Vec<u8>,
    _edge_map: Vec<usize>, // maps the edges in this graph to the edges in the parent graph
    _inverse_edge_map: Vec<Option<usize>>, // maps the edges in the parent graph to the edges in this graph
    signature_matrix: Vec<Vec<isize>>,
}

impl TropicalGraph {
    /// Create the simplfied graph from the full Graph struct
    pub fn from_graph(graph: &Graph, tropical_edge_weights: &[f64]) -> Self {
        let dod = tropical_edge_weights.iter().sum::<f64>()
            - D as f64 / 2.0 * graph.loop_momentum_basis.basis.len() as f64;

        let mut topology = Vec::with_capacity(graph.edges.len());
        let mut is_edge_massive = Vec::with_capacity(graph.edges.len());
        let mut edge_map = Vec::with_capacity(graph.edges.len());
        let mut external_vertices = Vec::new();
        let mut num_massive_edges = 0;
        let mut _inverse_edge_map = vec![None; graph.edges.len()];

        for (tropical_edge_index, (index, edge)) in graph.get_loop_edges_iterator().enumerate() {
            let tropical_edge = TropicalEdge {
                edge_id: tropical_edge_index,
                left: edge.vertices[0] as u8,
                right: edge.vertices[1] as u8,
                weight: tropical_edge_weights[tropical_edge_index],
            };

            topology.push(tropical_edge);

            if let Some(mass) = edge.particle.mass.value {
                if mass.is_zero() {
                    is_edge_massive.push(false);
                } else {
                    is_edge_massive.push(true);
                    num_massive_edges += 1;
                }
            } else {
                is_edge_massive.push(false);
            }

            edge_map.push(index);
            _inverse_edge_map[index] = Some(tropical_edge_index);
        }

        // Candidates to be external vertices of the tree-stripped graph
        let mut external_vertices_pool = HashSet::default();

        for edge in graph
            .edges
            .iter()
            .filter(|e| e.edge_type == EdgeType::Incoming)
            .chain(
                graph
                    .edges
                    .iter()
                    .filter(|e| e.edge_type == EdgeType::Outgoing),
            )
        {
            external_vertices_pool.insert(edge.vertices[0]);
            external_vertices_pool.insert(edge.vertices[1]);
        }

        for (_, edge) in graph.get_tree_level_edges_iterator() {
            external_vertices_pool.insert(edge.vertices[0]);
            external_vertices_pool.insert(edge.vertices[1]);
        }

        // if an edge contains a vertex in the external vertex pool, this vertex must be labeled as external for the tropical sampling.
        for edge in &topology {
            if external_vertices_pool.contains(&(edge.left as usize)) {
                external_vertices.push(edge.left);
            }

            if external_vertices_pool.contains(&(edge.right as usize)) {
                external_vertices.push(edge.right);
            }
        }

        let edge_number = topology.len();
        let loop_number = graph.loop_momentum_basis.basis.len();

        let mut signature_matrix = vec![vec![0; loop_number]; edge_number];
        for e in 0..edge_number {
            for l in 0..loop_number {
                signature_matrix[e][l] =
                    graph.loop_momentum_basis.edge_signatures[edge_map[e]].0[l];
            }
        }

        Self {
            dod,
            is_edge_massive,
            topology,
            num_massive_edges,
            external_vertices,
            _edge_map: edge_map,
            _inverse_edge_map,
            signature_matrix,
        }
    }

    fn get_full_subgraph_id(&self) -> TropicalSubGraphId {
        TropicalSubGraphId::new(self.topology.len())
    }

    /// subgraph gamma of a parent graph G is mass-spanning if it contains all massive propagators of G.
    /// and momentum-spanning if it has a connected component that contains all external vertices of G.
    /// A mass-momentum-spanning subgraph is both mass-spanning and momentum-spanning.
    fn is_mass_momentum_spanning(&self, edges_in_subgraph: &[usize]) -> bool {
        let num_massive_edges = edges_in_subgraph
            .iter()
            .filter(|&&i| self.is_edge_massive[i])
            .count();

        let is_mass_spanning = num_massive_edges == self.num_massive_edges;

        // this does not check the connected component property, needs to be fixe

        let connected_compoenents = self.get_connected_components(edges_in_subgraph);

        let is_momentum_spanning = connected_compoenents.iter().any(|component| {
            let edges_in_connected_subgraph = component.contains_edges();
            self.external_vertices.iter().all(|&v| {
                edges_in_connected_subgraph
                    .iter()
                    .any(|&i| self.topology[i].contains_vertex(v))
            })
        });

        is_mass_spanning && is_momentum_spanning
    }

    /// Get all connected components of a graph, used to compute loop number of possible disconnected graph
    fn get_connected_components(&self, edges_in_subgraph: &[usize]) -> Vec<TropicalSubGraphId> {
        let num_edges_in_subgraph = edges_in_subgraph.len();

        if num_edges_in_subgraph == 0 {
            return vec![];
        }

        let mut visited_edges: HashSet<usize> = HashSet::default();
        let mut connected_components: Vec<HashSet<usize>> = vec![];
        let mut current_component: HashSet<usize> = HashSet::default();

        // start search in the first edge
        let mut current_edges = vec![edges_in_subgraph[0]];

        visited_edges.insert(current_edges[0]);
        current_component.insert(current_edges[0]);

        while num_edges_in_subgraph
            > connected_components
                .iter()
                .map(std::collections::HashSet::len)
                .sum::<usize>()
        {
            let neighbours = current_edges
                .iter()
                .flat_map(|&edge_id| {
                    self.get_neighbouring_edges_in_subgraph(edge_id, edges_in_subgraph)
                })
                .collect::<Vec<usize>>();

            let mut current_component_grown = false;
            for &neighbour in &neighbours {
                current_edges.push(neighbour);
                visited_edges.insert(neighbour);
                if current_component.insert(neighbour) {
                    current_component_grown = true;
                };
            }

            if !current_component_grown {
                connected_components.push(current_component);
                current_component = HashSet::default();
                current_edges.clear();
                for &edge_id in edges_in_subgraph {
                    if !visited_edges.contains(&edge_id) {
                        current_edges.push(edge_id);
                        visited_edges.insert(edge_id);
                        current_component.insert(edge_id);
                        break;
                    }
                }
            }
        }

        connected_components
            .into_iter()
            .map(|component| {
                TropicalSubGraphId::from_edge_list(
                    &component.into_iter().collect::<Vec<usize>>(),
                    self.topology.len(),
                )
            })
            .collect()
    }

    /// get the loop number of a potentially disconnected grraph
    fn get_loop_number(&self, edges_in_subgraph: &[usize]) -> u8 {
        if edges_in_subgraph.is_empty() {
            return 0;
        }

        let connected_components = self.get_connected_components(edges_in_subgraph);

        connected_components
            .iter()
            .map(|c| self.get_loop_number_of_connected_component(c))
            .sum()
    }

    /// Get the loop number of a connected graph, by Euler's formula
    fn get_loop_number_of_connected_component(&self, subgraph_id: &TropicalSubGraphId) -> u8 {
        let edges_in_connected_subgraph = subgraph_id.contains_edges();
        let mut vertices: HashSet<u8> = HashSet::default();

        for &edge in &edges_in_connected_subgraph {
            vertices.insert(self.topology[edge].left);
            vertices.insert(self.topology[edge].right);
        }

        let num_vertices = vertices.len();
        let num_edges = edges_in_connected_subgraph.len();
        1 + num_edges as u8 - num_vertices as u8
    }

    /// Sum all the weights of the edges in the subgraph
    fn compute_weight_sum(&self, edges_in_subgraph: &[usize]) -> f64 {
        edges_in_subgraph
            .iter()
            .map(|&i| self.topology[i].weight)
            .sum()
    }

    /// Get all the edges that share a vertex with a given edge
    fn get_neighbouring_edges_in_subgraph(
        &self,
        edge_id: usize,
        edges_in_subgraph: &[usize],
    ) -> Vec<usize> {
        edges_in_subgraph
            .iter()
            .filter(|&&i| self.are_neighbours(edge_id, i))
            .copied()
            .collect()
    }

    /// Check if two edges share a vertex.
    fn are_neighbours(&self, edge_id_1: usize, edge_id_2: usize) -> bool {
        self.topology[edge_id_1].contains_vertex(self.topology[edge_id_2].left)
            || self.topology[edge_id_1].contains_vertex(self.topology[edge_id_2].right)
    }

    /// Definition of the j-function for a subgraph, described in the tropical sampling papers.
    fn recursive_j_function_eval(
        subgraph_id: &TropicalSubGraphId,
        table: &mut Vec<OptionTropicalSubgraphTableEntry>,
    ) -> f64 {
        if subgraph_id.is_empty() {
            let j_function = 1.0;
            table[subgraph_id.id].j_function = Some(j_function);
            j_function
        } else {
            // this guards against infinite recursion, for some reason it occurs in the 4-loop ladder
            if let Some(j_function) = table[subgraph_id.id].j_function {
                return j_function;
            }

            let edges_in_subgraph = subgraph_id.contains_edges();
            let subgraphs = edges_in_subgraph.iter().map(|&e| subgraph_id.pop_edge(e));

            let j_function = subgraphs
                .map(|g| {
                    TropicalGraph::recursive_j_function_eval(&g, table)
                        / table[g.id].generalized_dod.unwrap()
                })
                .sum();
            table[subgraph_id.id].j_function = Some(j_function);
            j_function
        }
    }
}

/// The bits in the id represent whether the corresponding edge is in the subgraph
#[derive(Debug, Clone, Copy, PartialEq)]
struct TropicalSubGraphId {
    id: usize,
    num_edges: usize,
}

impl TropicalSubGraphId {
    fn new(num_edges: usize) -> Self {
        Self {
            id: (1 << num_edges) - 1,
            num_edges,
        }
    }

    fn from_id(id: usize, num_edges: usize) -> Self {
        Self { id, num_edges }
    }

    /// remove an edge from the subgraph
    fn pop_edge(&self, edge_id: usize) -> Self {
        Self {
            id: self.id ^ (1 << edge_id),
            num_edges: self.num_edges,
        }
    }

    /// The 0 id as all bits set to 0, so it contains no edges and thus represents the empty graph
    fn is_empty(&self) -> bool {
        self.id == 0
    }

    /// Check wheter the subgraph contains a specific edge
    fn has_edge(&self, edge_id: usize) -> bool {
        self.id & (1 << edge_id) != 0
    }

    /// Get the edges contained in the subgraph
    fn contains_edges(&self) -> Vec<usize> {
        (0..self.num_edges).filter(|&i| self.has_edge(i)).collect()
    }

    /// Check if the subgraph contains only one edge, this is a special case in the tropical sampling algorithm
    fn has_one_edge(&self) -> bool {
        self.id.count_ones() == 1
    }

    /// Create a subgraph id from a list of edges
    fn from_edge_list(edge_list: &[usize], num_edges: usize) -> Self {
        let mut id = 0;
        for &edge_id in edge_list {
            id |= 1 << edge_id;
        }
        Self { id, num_edges }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct TropicalEdge {
    edge_id: usize, // position in the list of edges
    left: u8,
    right: u8,
    pub weight: f64,
}

impl TropicalEdge {
    fn contains_vertex(&self, vertex: u8) -> bool {
        self.left == vertex || self.right == vertex
    }
}

/// Helper struct for generation, to store entries before computing the value of the j_function (Which requires the value of the other quantities to be present)
#[derive(Debug, Clone, Copy)]
struct OptionTropicalSubgraphTableEntry {
    loop_number: Option<u8>,
    mass_momentum_spanning: Option<bool>,
    j_function: Option<f64>,
    generalized_dod: Option<f64>,
}

impl OptionTropicalSubgraphTableEntry {
    fn all_none() -> Self {
        Self {
            loop_number: None,
            mass_momentum_spanning: None,
            j_function: None,
            generalized_dod: None,
        }
    }

    /// Check if all fields are set, and panic if it has failed
    fn to_entry(self) -> TropicalSubgraphTableEntry {
        TropicalSubgraphTableEntry {
            loop_number: self
                .loop_number
                .unwrap_or_else(|| panic!("loop number not set")),
            mass_momentum_spanning: self
                .mass_momentum_spanning
                .unwrap_or_else(|| panic!("mass-momentum spanning not set")),
            j_function: self
                .j_function
                .unwrap_or_else(|| panic!("j-function not set")),
            generalized_dod: self
                .generalized_dod
                .unwrap_or_else(|| panic!("generalized dod not set")),
        }
    }
}

/// Data that needs to be stored for each subgraph
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct TropicalSubgraphTableEntry {
    pub loop_number: u8,
    pub mass_momentum_spanning: bool,
    pub j_function: f64,
    pub generalized_dod: f64,
}

/// The list of data for all subgraphs, indexed using the TropicalSubGraphId
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TropicalSubgraphTable {
    pub table: Vec<TropicalSubgraphTableEntry>,
    pub tropical_graph: TropicalGraph,
}

impl TropicalSubgraphTable {
    /// Generate the tropical subgraph table from a graph
    pub fn generate_from_graph(
        graph: &Graph,
        tropical_edge_weights: &[f64],
    ) -> Result<Self, Report> {
        debug!("🌴🥥 Generating tropical subgraph table 🥥🌴");

        let tropical_graph = TropicalGraph::from_graph(graph, tropical_edge_weights);
        Self::generate_from_tropical(&tropical_graph)
    }

    fn generate_from_tropical(tropical_graph: &TropicalGraph) -> Result<Self, Report> {
        let num_edges = tropical_graph.topology.len();
        let powerset_size = 2usize.pow(num_edges as u32);

        // allocate the subgraph table
        let mut option_subgraph_table =
            vec![OptionTropicalSubgraphTableEntry::all_none(); powerset_size];

        // create iterator over all subgraphs
        let subgraph_iterator =
            (0..powerset_size).map(|i| TropicalSubGraphId::from_id(i, num_edges));

        let full_subgraph_id = tropical_graph.get_full_subgraph_id();

        for subgraph in subgraph_iterator {
            let edges_in_subgraph = subgraph.contains_edges();

            // check the mass-momentum spanning property
            let is_mass_momentum_spanning =
                tropical_graph.is_mass_momentum_spanning(&edges_in_subgraph);

            option_subgraph_table[subgraph.id].mass_momentum_spanning =
                Some(is_mass_momentum_spanning);

            // compute the sum of weights of all edges in the subgraph
            let weight_sum = tropical_graph.compute_weight_sum(&edges_in_subgraph);
            let loop_number = tropical_graph.get_loop_number(&edges_in_subgraph);

            let generalized_dod = if !subgraph.is_empty() {
                if is_mass_momentum_spanning {
                    weight_sum - loop_number as f64 * D as f64 / 2.0 - tropical_graph.dod
                } else {
                    weight_sum - loop_number as f64 * D as f64 / 2.0
                }
            } else {
                1.0
            };

            if generalized_dod <= 0.0 && !subgraph.is_empty() && subgraph != full_subgraph_id {
                return Err(Report::msg(format!(
                    "Generalized DoD: {generalized_dod} is negative for subgraph {subgraph:?}\n
                    loop number: {loop_number}, mass-momentum spanning: {is_mass_momentum_spanning}, weight sum: {weight_sum}"
                )));
            }

            option_subgraph_table[subgraph.id].loop_number = Some(loop_number);
            option_subgraph_table[subgraph.id].generalized_dod = Some(generalized_dod);
        }

        TropicalGraph::recursive_j_function_eval(&full_subgraph_id, &mut option_subgraph_table);

        Ok(TropicalSubgraphTable {
            table: option_subgraph_table
                .into_iter()
                .map(OptionTropicalSubgraphTableEntry::to_entry)
                .collect(),
            tropical_graph: tropical_graph.clone(),
        })
    }

    /// sample an edge from a subgraph, according to the relevant probability distribution, returns the edge and the subgraph without the edge for later use
    /// Panics if the uniform random number is greater than one or the probability distribution is incorrectly normalized.
    fn sample_edge<T: FloatLike>(
        &self,
        uniform: T,
        subgraph: &TropicalSubGraphId,
    ) -> (usize, TropicalSubGraphId) {
        let edges_in_subgraph = subgraph.contains_edges();
        let j = Into::<T>::into(self.table[subgraph.id].j_function);

        let mut cum_sum = T::zero();
        for &edge in &edges_in_subgraph {
            let graph_without_edge = subgraph.pop_edge(edge);
            let p_e = Into::<T>::into(self.table[graph_without_edge.id].j_function)
                / Into::<T>::into(j)
                / Into::<T>::into(self.table[graph_without_edge.id].generalized_dod);
            cum_sum += p_e;
            if cum_sum >= uniform {
                return (edge, graph_without_edge);
            }
        }

        panic!("Sampling could not sample edge, with uniform_random_number: {}, cumulative sum evaluated to: {}", uniform, cum_sum);
    }
}

pub mod tropical_parameterization {
    use color_eyre::Report;

    use itertools::{izip, Itertools};
    use lorentz_vector::LorentzVector;
    use statrs::function::gamma::gamma;

    use crate::{
        graph::Graph,
        linalg::SquareMatrix,
        utils::{box_muller, inverse_gamma_lr, FloatLike},
    };

    use super::{TropicalSubgraphTable, D};

    /// Fake random number generator, perhaps it is more readable to get rid of this.
    struct MimicRng<'a, T> {
        cache: &'a [T],
        counter: usize,
        tokens: Vec<&'a str>,
    }

    impl<'a, T: Copy> MimicRng<'a, T> {
        fn get_random_number(&mut self, token: Option<&'a str>) -> T {
            let random_number = self.cache[self.counter];
            self.counter += 1;
            if let Some(token) = token {
                self.tokens.push(token);
            }
            random_number
        }

        fn new(cache: &'a [T]) -> Self {
            Self {
                cache,
                counter: 0,
                tokens: Vec::new(),
            }
        }
    }

    struct PermatuhedralSamplingResult<T> {
        x: Vec<T>,
        u_trop: T,
        v_trop: T,
    }

    /// This function returns the feynman parameters for a given graph and sample point, it also computes u_trop and v_trop.
    /// A rescaling is performed for numerical stability, with this rescaling u_trop and v_trop always evaluate to 1.
    fn permatuhedral_sampling<T: FloatLike>(
        tropical_subgraph_table: &TropicalSubgraphTable,
        rng: &mut MimicRng<T>,
        debug: usize,
    ) -> PermatuhedralSamplingResult<T> {
        let mut kappa = T::one();
        let mut x_vec = vec![T::zero(); tropical_subgraph_table.tropical_graph.topology.len()];
        let mut u_trop = T::one();
        let mut v_trop = T::one();

        let mut graph = tropical_subgraph_table
            .tropical_graph
            .get_full_subgraph_id();

        while !graph.is_empty() {
            // this saves a random variable
            let (edge, graph_without_edge) = if graph.has_one_edge() {
                let edge = graph.contains_edges()[0];
                let graph_without_edge = graph.pop_edge(edge);
                (edge, graph_without_edge)
            } else {
                tropical_subgraph_table
                    .sample_edge(rng.get_random_number(Some("sample_edge")), &graph)
            };

            x_vec[edge] = kappa;

            if tropical_subgraph_table.table[graph.id].mass_momentum_spanning
                && !tropical_subgraph_table.table[graph_without_edge.id].mass_momentum_spanning
            {
                v_trop = x_vec[edge];
            }

            if tropical_subgraph_table.table[graph_without_edge.id].loop_number
                < tropical_subgraph_table.table[graph.id].loop_number
            {
                u_trop *= x_vec[edge];
            }

            // Terminate early, so we do not waste a random variable in the final step
            graph = graph_without_edge;
            if graph.is_empty() {
                break;
            }

            let xi = rng.get_random_number(Some("sample xi"));
            kappa *= xi.powf(
                Into::<T>::into(tropical_subgraph_table.table[graph.id].generalized_dod).inv(),
            );

            if debug > 4 {
                println!(
                    "generalized_dod: {}",
                    tropical_subgraph_table.table[graph.id].generalized_dod
                );
                println!("xi: {xi}");
                println!("kappa: {kappa}");
            }
        }

        let xi_trop = u_trop * v_trop;

        // perform rescaling for numerical stability
        let target = u_trop.powf(Into::<T>::into(-(D as f64 / 2.0)))
            * (u_trop / xi_trop).powf(Into::<T>::into(tropical_subgraph_table.tropical_graph.dod));

        let loop_number = tropical_subgraph_table.table.last().unwrap().loop_number;
        let scaling = target.powf(
            Into::<T>::into(
                D as f64 / 2.0 * loop_number as f64 + tropical_subgraph_table.tropical_graph.dod,
            )
            .inv(),
        );

        if debug > 2 {
            println!("scaling: {scaling}");
            println!("u_trop before scaling: {u_trop}");
            println!("xi_trop before scaling: {xi_trop}");
            println!("v_trop before scaling: {v_trop}");
        }

        x_vec.iter_mut().for_each(|x| *x *= scaling);
        u_trop = T::one();
        v_trop = T::one();

        PermatuhedralSamplingResult {
            x: x_vec,
            u_trop,
            v_trop,
        }
    }

    /// In this version the steps are reordered, perhaps useful to keep when we want to do learning on top of tropical sampling, but for now it is not used.
    /// It has been checked that it produces the same result as the normal version.
    fn _adapted_permatuhedral_sampling<T: FloatLike>(
        tropical_subgraph_table: &TropicalSubgraphTable,
        rng: &mut MimicRng<T>,
        edge_rng: &[T],
        _debug: usize,
    ) -> PermatuhedralSamplingResult<T> {
        let mut graph = tropical_subgraph_table
            .tropical_graph
            .get_full_subgraph_id();

        let mut graph_ordering = Vec::with_capacity(graph.num_edges + 1);
        graph_ordering.push((graph, None));

        while !graph.is_empty() {
            let (edge, graph_without_edge) = if graph.has_one_edge() {
                let edge = graph.contains_edges()[0];
                let graph_without_edge = graph.pop_edge(edge);
                (edge, graph_without_edge)
            } else {
                tropical_subgraph_table
                    .sample_edge(rng.get_random_number(Some("sample_edge")), &graph)
            };

            graph_ordering.push((graph_without_edge, Some(edge)));
            graph = graph_without_edge;
        }

        let mut kappa = T::one();
        let mut x_vec = vec![T::zero(); tropical_subgraph_table.tropical_graph.topology.len()];
        let mut u_trop = T::one();
        let mut v_trop = T::one();

        for i in 1..=graph.num_edges {
            let (graph, _) = graph_ordering[i - 1];
            if let (graph_without_edge, Some(edge)) = graph_ordering[i] {
                x_vec[edge] = kappa;

                if tropical_subgraph_table.table[graph.id].mass_momentum_spanning
                    && !tropical_subgraph_table.table[graph_without_edge.id].mass_momentum_spanning
                {
                    v_trop = x_vec[edge];
                }

                if tropical_subgraph_table.table[graph_without_edge.id].loop_number
                    < tropical_subgraph_table.table[graph.id].loop_number
                {
                    u_trop *= x_vec[edge];
                }

                if graph_without_edge.is_empty() {
                    break;
                }

                let graph_without_edge_dod =
                    tropical_subgraph_table.table[graph_without_edge.id].generalized_dod;

                let next_edge = graph_ordering[i + 1].1.unwrap();
                let xi = edge_rng[next_edge];
                kappa *= xi.powf(Into::<T>::into(graph_without_edge_dod).inv());
            } else {
                unreachable!()
            }
        }

        let xi_trop = u_trop * v_trop;

        // perform rescaling for numerical stability
        let target = u_trop.powf(Into::<T>::into(-(D as f64 / 2.0)))
            * (u_trop / xi_trop).powf(Into::<T>::into(tropical_subgraph_table.tropical_graph.dod));

        let loop_number = tropical_subgraph_table.table.last().unwrap().loop_number;
        let scaling = target.powf(
            Into::<T>::into(
                D as f64 / 2.0 * loop_number as f64 + tropical_subgraph_table.tropical_graph.dod,
            )
            .inv(),
        );

        x_vec.iter_mut().for_each(|x| *x *= scaling);
        u_trop = T::one();
        v_trop = T::one();

        PermatuhedralSamplingResult {
            x: x_vec,
            u_trop,
            v_trop,
        }
    }

    /// This function is called from the paramatrization stage of evaluate_sample
    pub fn generate_tropical_sample<T: FloatLike + Into<f64>>(
        x_space_point: &[T],
        external_momenta: &[LorentzVector<T>],
        graph: &Graph,
        debug: usize,
    ) -> Result<(Vec<LorentzVector<T>>, T), Report> {
        let tropical_subgraph_table = graph.derived_data.tropical_subgraph_table.as_ref().unwrap();
        let signature_matrix = &tropical_subgraph_table.tropical_graph.signature_matrix;
        let _num_edges = signature_matrix.len();
        let num_loops = signature_matrix[0].len();

        let mut rng = MimicRng::new(x_space_point);

        // perhaps this data should be stored in some cache in the integrand.
        let (edge_masses, edge_shifts): (Vec<T>, Vec<LorentzVector<T>>) = tropical_subgraph_table
            .tropical_graph
            .topology
            .iter()
            .enumerate()
            .map(|(index, _edge)| {
                let edge_mass = Into::<T>::into(
                    match graph.edges[tropical_subgraph_table.tropical_graph._edge_map[index]]
                        .particle
                        .mass
                        .value
                    {
                        Some(mass) => mass.re,
                        None => 0.0,
                    },
                );

                let external_signature = &graph.loop_momentum_basis.edge_signatures
                    [tropical_subgraph_table.tropical_graph._edge_map[index]]
                    .1;

                let edge_shift = external_momenta.iter().enumerate().fold(
                    LorentzVector::<T>::new(),
                    |acc, (i, external_momentum)| {
                        let external_spatial = LorentzVector::from_args(
                            T::zero(),
                            external_momentum.x,
                            external_momentum.y,
                            external_momentum.z,
                        );
                        acc + external_spatial * Into::<T>::into(external_signature[i] as f64)
                    },
                );

                (edge_mass, edge_shift)
            })
            .unzip();

        if debug > 1 {
            println!("edge shifts: {:?}", edge_shifts);
            println!("edge masses: {:?}", edge_masses);
        }

        let permatuhedral_result = permatuhedral_sampling(tropical_subgraph_table, &mut rng, debug);

        if debug > 1 {
            println!("feynman parameters: {:?}", permatuhedral_result.x);
            println!("u_trop: {}", permatuhedral_result.u_trop);
            println!("v_trop: {}", permatuhedral_result.v_trop);
            println!("signature matrix: {:?}", signature_matrix);
        }

        let l_matrix = compute_l_matrix(&permatuhedral_result.x, signature_matrix);
        let decomposed_l_matrix = l_matrix.decompose_for_tropical()?;

        if debug > 1 {
            println!("L matrix: {:?}", l_matrix);
        }

        let lambda = Into::<T>::into(inverse_gamma_lr(
            tropical_subgraph_table.tropical_graph.dod,
            Into::<f64>::into(rng.get_random_number(Some("sample lambda"))),
            50,
        ));

        if debug > 1 {
            println!("lambda: {}", lambda);
        }

        let q_vectors = sample_q_vectors(&mut rng, num_loops);
        let u_vectors = compute_u_vectors(&permatuhedral_result.x, signature_matrix, &edge_shifts);

        if debug > 1 {
            println!("u vectors: {:?}", u_vectors);
            println!("q vectors: {:?}", q_vectors);
        }

        if debug > 1 {
            println!("det L: {}", decomposed_l_matrix.determinant);
        }

        let v_polynomial = compute_v_polynomial(
            &permatuhedral_result.x,
            &u_vectors,
            &decomposed_l_matrix.inverse,
            &edge_shifts,
            &edge_masses,
            debug,
        );

        if debug > 1 {
            println!("v polynomial: {}", v_polynomial);
        }

        let loop_momenta = compute_loop_momenta(
            v_polynomial,
            lambda,
            &decomposed_l_matrix.q_transposed_inverse,
            &q_vectors,
            &decomposed_l_matrix.inverse,
            &u_vectors,
            debug,
        );

        if debug > 1 {
            println!("loop momenta: {:?}", loop_momenta);
        }

        let jacobian = compute_jacobian(
            tropical_subgraph_table,
            permatuhedral_result.u_trop,
            permatuhedral_result.v_trop,
            v_polynomial,
            decomposed_l_matrix.determinant,
            num_loops,
            debug,
        );

        if debug > 1 {
            println!("jacobian: {}", jacobian);
            println!("tokens: {:?}", rng.tokens);
        }

        Ok((loop_momenta, jacobian))
    }

    /// Compute the L x L matrix from the feynman parameters and the signature matrix
    #[inline]
    fn compute_l_matrix<T: FloatLike>(
        x_vec: &[T],
        signature_matrix: &[Vec<isize>],
    ) -> SquareMatrix<T> {
        let num_edges = signature_matrix.len();
        let num_loops = signature_matrix[0].len();

        let mut temp_l_matrix = SquareMatrix::new_zeros(num_loops);

        for i in 0..num_loops {
            for j in 0..num_loops {
                for e in 0..num_edges {
                    temp_l_matrix[(i, j)] += x_vec[e]
                        * Into::<T>::into((signature_matrix[e][i] * signature_matrix[e][j]) as f64);
                }
            }
        }

        temp_l_matrix
    }

    /// Sample Gaussian distributed vectors, using the Box-Muller transform
    #[inline]
    fn sample_q_vectors<T: FloatLike>(
        rng: &mut MimicRng<T>,
        num_loops: usize,
    ) -> Vec<LorentzVector<T>> {
        let token = Some("box muller");
        let num_variables = D * num_loops;

        let num_uniform_variables = num_variables + num_variables % 2;
        let gaussians = (0..num_uniform_variables / 2).flat_map(|_| {
            let (box_muller_1, box_muller_2) =
                box_muller(rng.get_random_number(token), rng.get_random_number(token));

            [box_muller_1, box_muller_2]
        });

        #[allow(clippy::useless_conversion)] // without the conversion I get an error
        (0..num_loops)
            .zip(gaussians.chunks(3).into_iter())
            .map(|(_, mut chunk)| {
                LorentzVector::from_args(
                    T::zero(),
                    chunk.next().unwrap(),
                    chunk.next().unwrap(),
                    chunk.next().unwrap(),
                )
            })
            .collect_vec()
    }

    /// Compute the vectors u, according to the formula in the notes
    #[inline]
    fn compute_u_vectors<T: FloatLike>(
        x_vec: &[T],
        signature_marix: &[Vec<isize>],
        edge_shifts: &[LorentzVector<T>],
    ) -> Vec<LorentzVector<T>> {
        let num_loops = signature_marix[0].len();
        let num_edges = signature_marix.len();

        (0..num_loops)
            .map(|l| {
                (0..num_edges).fold(LorentzVector::new(), |acc: LorentzVector<T>, e| {
                    acc + edge_shifts[e] * x_vec[e] * Into::<T>::into(signature_marix[e][l] as f64)
                })
            })
            .collect_vec()
    }

    /// Compute the polynomial v, according to the formula in the notes
    #[inline]
    fn compute_v_polynomial<T: FloatLike>(
        x_vec: &[T],
        u_vectors: &[LorentzVector<T>],
        inverse_l: &SquareMatrix<T>,
        edge_shifts: &[LorentzVector<T>],
        edge_masses: &[T],
        debug: usize,
    ) -> T {
        let num_loops = inverse_l.get_dim();

        if debug > 6 {
            println!("analyzing computation_of_polynomial");
            let mut term_1 = T::zero();
            for (&x_e, &mass, &shift) in izip!(x_vec, edge_masses, edge_shifts) {
                let quant = x_e * (mass * mass + shift.spatial_squared());
                println!("x_e: {x_e}, mass: {mass}, shift: {shift}");
                println!("quant: {quant}");
                term_1 += quant;
            }

            println!("term_1: {term_1}");

            let mut term_2 = T::zero();
            for (i, j) in (0..num_loops).cartesian_product(0..num_loops) {
                let quant = u_vectors[i].spatial_dot(&u_vectors[j]) * inverse_l[(i, j)];
                println!(
                    "u_i: {}, u_j: {}, inverse_l: {}",
                    u_vectors[i],
                    u_vectors[j],
                    inverse_l[(i, j)]
                );
                println!("quant: {quant}");
                term_2 += quant;
            }

            println!("term_2: {term_2}");
            let res = term_1 - term_2;
            println!("res: {res}");
            res
        } else {
            let term_1 = izip!(x_vec, edge_masses, edge_shifts)
                .map(|(&x_e, &mass, &shift)| x_e * (mass * mass + shift.spatial_squared()))
                .sum::<T>();

            let term_2 = (0..num_loops)
                .cartesian_product(0..num_loops)
                .map(|(i, j)| u_vectors[i].spatial_dot(&u_vectors[j]) * inverse_l[(i, j)])
                .sum::<T>();

            term_1 - term_2
        }
    }

    /// Compute the loop momenta, according to the formula in the notes
    #[inline]
    fn compute_loop_momenta<T: FloatLike>(
        v: T,
        lambda: T,
        q_t_inverse: &SquareMatrix<T>,
        q_vectors: &[LorentzVector<T>],
        l_inverse: &SquareMatrix<T>,
        u_vectors: &[LorentzVector<T>],
        debug: usize,
    ) -> Vec<LorentzVector<T>> {
        let num_loops = q_t_inverse.get_dim();
        let prefactor = (v / lambda * Into::<T>::into(0.5)).sqrt();

        if debug > 6 {
            for l in 0..num_loops {
                let mut weird_shift: LorentzVector<T> = LorentzVector::new();
                for j in 0..num_loops {
                    let term_in_weird_shift = u_vectors[j] * l_inverse[(l, j)];
                    println!("term in weird shift: {term_in_weird_shift}");
                    weird_shift += u_vectors[j] * l_inverse[(l, j)];
                }
                println!("weird shift: {weird_shift}");
            }
        }
        (0..num_loops)
            .map(|l| {
                q_vectors.iter().zip(u_vectors.iter()).enumerate().fold(
                    LorentzVector::new(),
                    |acc, (l_prime, (q, u))| {
                        acc + q * prefactor * q_t_inverse[(l, l_prime)]
                            - u * l_inverse[(l, l_prime)]
                    },
                )
            })
            .collect_vec()
    }

    /// Compute the Jacobian of the tropical sampling algorithm. This does not include the factors of positive powers of E that need to
    /// be multiplied with the three-dimensional representation. Perhaps it is more readable to include that here.
    #[inline]
    fn compute_jacobian<T: FloatLike>(
        tropical_subgraph_table: &TropicalSubgraphTable,
        u_trop: T,
        v_trop: T,
        v: T,
        u: T,
        num_loops: usize,
        debug: usize,
    ) -> T {
        let polynomial_ratio = (v_trop / v)
            .powf(Into::<T>::into(tropical_subgraph_table.tropical_graph.dod))
            * (u_trop / u).powf(Into::<T>::into(D as f64 / 2.0));

        let i_trop = Into::<T>::into(tropical_subgraph_table.table.last().unwrap().j_function);

        if debug > 1 {
            println!("i_trop: {i_trop}");
        }

        let gamma_omega = Into::<T>::into(gamma(tropical_subgraph_table.tropical_graph.dod));
        let denom = Into::<T>::into(
            tropical_subgraph_table
                .tropical_graph
                .topology
                .iter()
                .map(|e| gamma(e.weight))
                .product::<f64>(),
        );

        if debug > 1 {
            println!("gamma omega: {gamma_omega}");
            println!("gamma nu product: {denom}");
        }

        // we divide by (2pi)^L later, tropcial sampling already contains a part of this, so we undo that here.
        let pi_power = Into::<T>::into(std::f64::consts::PI.powf((D * num_loops) as f64 / 2.0));

        // we include the factor of 2^E here, so we don't need to include it in the evaluate function, we still need to add it for the tree-like edges
        let num_edges = tropical_subgraph_table.tropical_graph.topology.len();
        let two_to_the_e: usize = 1 << num_edges;

        i_trop * gamma_omega / denom * polynomial_ratio * pi_power
            / Into::<T>::into(two_to_the_e as f64)
    }
}

// some tests
#[cfg(test)]
mod tests {
    use crate::utils::assert_approx_eq;

    const TOLERANCE: f64 = 1e-14;

    use super::*;

    #[test]
    fn test_connected_components() {
        let topology = vec![
            TropicalEdge {
                edge_id: 0,
                left: 0,
                right: 1,
                weight: 1.0,
            },
            TropicalEdge {
                edge_id: 1,
                left: 1,
                right: 0,
                weight: 1.0,
            },
            TropicalEdge {
                edge_id: 2,
                left: 2,
                right: 3,
                weight: 1.0,
            },
            TropicalEdge {
                edge_id: 3,
                left: 3,
                right: 2,
                weight: 1.0,
            },
        ];

        let tropical_graph = TropicalGraph {
            dod: 0.0,
            is_edge_massive: vec![false; 4],
            topology,
            num_massive_edges: 0,
            external_vertices: vec![0, 1, 2, 3],
            _edge_map: vec![0, 1, 2, 3],
            _inverse_edge_map: vec![], // not needed for this test
            signature_matrix: vec![vec![0; 4]; 4], // not needed for this test
        };

        let components = tropical_graph.get_connected_components(&[0, 1, 2, 3]);
        assert_eq!(components.len(), 2);

        let components = tropical_graph.get_connected_components(&[0, 1]);
        assert_eq!(components.len(), 1);

        let components = tropical_graph.get_connected_components(&[0, 2]);
        assert_eq!(components.len(), 2);
    }

    #[test]
    fn test_loop_number() {
        let topology1 = vec![
            TropicalEdge {
                edge_id: 0,
                left: 0,
                right: 1,
                weight: 1.0,
            },
            TropicalEdge {
                edge_id: 1,
                left: 1,
                right: 0,
                weight: 1.0,
            },
            TropicalEdge {
                edge_id: 2,
                left: 2,
                right: 3,
                weight: 1.0,
            },
            TropicalEdge {
                edge_id: 3,
                left: 3,
                right: 2,
                weight: 1.0,
            },
        ];

        let tropical_graph1 = TropicalGraph {
            dod: 0.0,
            is_edge_massive: vec![false; 4],
            topology: topology1,
            num_massive_edges: 0,
            external_vertices: vec![0, 1, 2, 3],
            _edge_map: vec![0, 1, 2, 3],
            _inverse_edge_map: vec![], // not needed for this test
            signature_matrix: vec![vec![0; 4]; 4], // not needed for this test
        };

        let loop_number = tropical_graph1.get_loop_number(&[0, 1, 2, 3]);
        assert_eq!(loop_number, 2);

        let loop_number = tropical_graph1.get_loop_number(&[0, 1]);
        assert_eq!(loop_number, 1);

        let loop_number = tropical_graph1.get_loop_number(&[0, 2]);
        assert_eq!(loop_number, 0);
    }

    // tests compared against the output of feyntrop
    #[test]
    fn test_triangle() {
        let triangle_topology = vec![
            TropicalEdge {
                edge_id: 0,
                left: 0,
                right: 1,
                weight: 0.66,
            },
            TropicalEdge {
                edge_id: 1,
                left: 1,
                right: 2,
                weight: 0.66,
            },
            TropicalEdge {
                edge_id: 2,
                left: 2,
                right: 0,
                weight: 0.66,
            },
        ];

        let triangle_graph = TropicalGraph {
            dod: 3. * 0.66 - 3. / 2.,
            is_edge_massive: vec![false; 3],
            topology: triangle_topology,
            num_massive_edges: 0,
            external_vertices: vec![0, 1, 2],
            _edge_map: vec![0, 1, 2],
            _inverse_edge_map: vec![], // not needed for this test
            signature_matrix: vec![vec![0; 3]; 3], // not needed for this test
        };

        let subgraph_table = TropicalSubgraphTable::generate_from_tropical(&triangle_graph)
            .expect("Failed to generate subgraph table");

        assert_eq!(subgraph_table.table.len(), 8);

        assert_eq!(
            subgraph_table.table[0],
            TropicalSubgraphTableEntry {
                loop_number: 0,
                mass_momentum_spanning: false,
                j_function: 1.0,
                generalized_dod: 1.0,
            }
        );

        assert_eq!(
            subgraph_table.table[1],
            TropicalSubgraphTableEntry {
                loop_number: 0,
                mass_momentum_spanning: false,
                j_function: 1.0,
                generalized_dod: 0.66,
            }
        );

        assert_eq!(subgraph_table.table[2], subgraph_table.table[1]);
        assert_eq!(subgraph_table.table[4], subgraph_table.table[1]);

        for i in [3, 5, 6] {
            let table = subgraph_table.table[i];
            assert!(table.mass_momentum_spanning);
            assert_eq!(table.loop_number, 0);
            assert_approx_eq(table.generalized_dod, 0.84, TOLERANCE);
            assert_approx_eq(table.j_function, 3.030_303_030_303_03, TOLERANCE);
        }

        let final_table = subgraph_table.table[7];
        assert!(final_table.mass_momentum_spanning);
        assert_eq!(final_table.loop_number, 1);
        assert_approx_eq(final_table.generalized_dod, 0.0, TOLERANCE);
        assert_approx_eq(final_table.j_function, 10.822_510_822_510_82, TOLERANCE);
    }

    #[test]
    fn test_sunrise() {
        let sunrise_topology = vec![
            TropicalEdge {
                edge_id: 0,
                left: 0,
                right: 1,
                weight: 0.66,
            },
            TropicalEdge {
                edge_id: 1,
                left: 0,
                right: 1,
                weight: 0.66,
            },
            TropicalEdge {
                edge_id: 2,
                left: 0,
                right: 1,
                weight: 0.66,
            },
        ];

        let sunrise_graph = TropicalGraph {
            dod: 3. * 0.66 - 6. / 2.,
            is_edge_massive: vec![false; 3],
            topology: sunrise_topology,
            num_massive_edges: 0,
            external_vertices: vec![0, 1],
            _edge_map: vec![0, 1, 2],
            _inverse_edge_map: vec![], // not needed for this test
            signature_matrix: vec![vec![0; 3]; 3], // not needed for this test
        };

        let subgraph_table = TropicalSubgraphTable::generate_from_tropical(&sunrise_graph)
            .expect("Failed to generate subgraph table");

        assert_eq!(subgraph_table.table.len(), 8);

        assert_eq!(
            subgraph_table.table[0],
            TropicalSubgraphTableEntry {
                loop_number: 0,
                mass_momentum_spanning: false,
                j_function: 1.0,
                generalized_dod: 1.0,
            }
        );

        for i in [1, 2, 4] {
            let table = subgraph_table.table[i];
            assert!(table.mass_momentum_spanning);
            assert_eq!(table.loop_number, 0);
            assert_approx_eq(table.generalized_dod, 1.68, TOLERANCE);
            assert_approx_eq(table.j_function, 1.0, TOLERANCE);
        }

        for i in [3, 5, 6] {
            let table = subgraph_table.table[i];
            assert!(table.mass_momentum_spanning);
            assert_eq!(table.loop_number, 1);
            assert_approx_eq(table.generalized_dod, 0.84, TOLERANCE);
            assert_approx_eq(table.j_function, 1.190_476_190_476_19, TOLERANCE);
        }

        let final_table = subgraph_table.table[7];
        assert!(final_table.mass_momentum_spanning);
        assert_eq!(final_table.loop_number, 2);
        assert_approx_eq(final_table.generalized_dod, 0.0, TOLERANCE);
        assert_approx_eq(final_table.j_function, 4.251_700_680_272_108, TOLERANCE);
    }

    #[test]
    fn mercedes() {
        let weight = 11. / 14.;
        let externals = vec![0, 2];

        let mercedes_topology = vec![
            TropicalEdge {
                edge_id: 0,
                left: 0,
                right: 1,
                weight,
            },
            TropicalEdge {
                edge_id: 1,
                left: 1,
                right: 2,
                weight,
            },
            TropicalEdge {
                edge_id: 2,
                left: 2,
                right: 3,
                weight,
            },
            TropicalEdge {
                edge_id: 3,
                left: 3,
                right: 0,
                weight,
            },
            TropicalEdge {
                edge_id: 4,
                left: 1,
                right: 4,
                weight,
            },
            TropicalEdge {
                edge_id: 5,
                left: 2,
                right: 4,
                weight,
            },
            TropicalEdge {
                edge_id: 6,
                left: 3,
                right: 4,
                weight,
            },
        ];

        let gr = TropicalGraph {
            dod: 1.0,
            num_massive_edges: 0,
            is_edge_massive: vec![false; 7],
            topology: mercedes_topology,
            external_vertices: externals,
            _edge_map: vec![0, 1, 2, 3, 4, 5, 6, 7],
            _inverse_edge_map: vec![],
            signature_matrix: vec![vec![0; 3]; 3],
        };

        let subgraph_table = TropicalSubgraphTable::generate_from_tropical(&gr)
            .expect("Failed to generate subgraph table");

        let i_tr = subgraph_table.table.last().unwrap().j_function;

        assert_approx_eq(i_tr, 1_818.303_855_640_347_1, TOLERANCE);
    }
}
