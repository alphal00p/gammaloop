use ahash::HashSet;
use color_eyre::Report;

use crate::graph::Graph;

#[derive(Debug, Clone)]
struct TropicalGraph {
    dod: f64,
    is_edge_massive: Vec<bool>,
    topology: Vec<TropicalEdge>,
    num_massive_edges: usize,
    external_vertices: Vec<u8>,
}

impl TropicalGraph {
    pub fn from_graph(graph: &Graph, tropical_edge_weights: &[f64]) -> Self {
        todo!()
    }

    fn get_full_subgraph_id(&self) -> TropicalSubGraphId {
        TropicalSubGraphId::new(self.topology.len())
    }

    // subgraph gamma of a parent graph G is mass-spanning if it contains all massive propagators of G
    // and momentum-spanning if it contains all external vertices of G
    // a mass-momentum-spanning subgraph is both mass-spanning and momentum-spanning
    fn is_mass_momentum_spanning(&self, edges_in_subgraph: &[usize]) -> bool {
        let num_massive_edges = edges_in_subgraph
            .iter()
            .filter(|&&i| self.is_edge_massive[i])
            .count();

        let is_mass_spanning = num_massive_edges == self.num_massive_edges;

        let is_momentum_spanning = self.external_vertices.iter().all(|&v| {
            edges_in_subgraph
                .iter()
                .any(|&i| self.topology[i].contains_vertex(&v))
        });

        is_mass_spanning && is_momentum_spanning
    }

    fn get_connected_components(&self, edges_in_subgraph: &[usize]) -> Vec<TropicalSubGraphId> {
        let num_edges_in_subgraph = edges_in_subgraph.len();

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
                .map(|component| component.len())
                .sum::<usize>()
        {
            let neighbours = current_edges
                .iter()
                .flat_map(|edge_id| {
                    self.get_neighbouring_edges_in_subgraph(edge_id, edges_in_subgraph)
                })
                .collect::<Vec<usize>>();

            let mut current_component_grown = false;
            for neighbour in neighbours.iter() {
                current_edges.push(*neighbour);
                visited_edges.insert(*neighbour);
                if current_component.insert(*neighbour) {
                    current_component_grown = true;
                };
            }

            if !current_component_grown {
                connected_components.push(current_component);
                current_component = HashSet::default();
                current_edges.clear();
                for edge_id in edges_in_subgraph.iter() {
                    if !visited_edges.contains(edge_id) {
                        current_edges.push(*edge_id);
                        visited_edges.insert(*edge_id);
                        current_component.insert(*edge_id);
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

    fn get_loop_number_of_connected_component(&self, subgraph_id: &TropicalSubGraphId) -> u8 {
        let edges_in_connected_subgraph = subgraph_id.contains_edges();
        let mut vertices: HashSet<u8> = HashSet::default();

        for edge in edges_in_connected_subgraph.iter() {
            vertices.insert(self.topology[*edge].left);
            vertices.insert(self.topology[*edge].right);
        }

        let num_vertices = vertices.len();
        let num_edges = edges_in_connected_subgraph.len();
        1 + num_edges as u8 - num_vertices as u8
    }

    fn compute_weight_sum(&self, edges_in_subgraph: &[usize]) -> f64 {
        edges_in_subgraph
            .iter()
            .map(|&i| self.topology[i].weight)
            .sum()
    }

    fn get_neighbouring_edges_in_subgraph(
        &self,
        edge_id: &usize,
        edges_in_subgraph: &[usize],
    ) -> Vec<usize> {
        edges_in_subgraph
            .iter()
            .filter(|&&i| self.are_neighbours(edge_id, &i))
            .copied()
            .collect()
    }

    // check if two edges share a vertex
    fn are_neighbours(&self, edge_id_1: &usize, edge_id_2: &usize) -> bool {
        self.topology[*edge_id_1].contains_vertex(&self.topology[*edge_id_2].left)
            || self.topology[*edge_id_1].contains_vertex(&self.topology[*edge_id_2].right)
    }

    fn recursive_j_function_eval(
        subgraph_id: &TropicalSubGraphId,
        table: &mut Vec<OptionTropicalSubgraphTableEntry>,
    ) -> f64 {
        if subgraph_id.is_empty() {
            let j_function = 1.0;
            table[subgraph_id.id].j_function = Some(j_function);
            j_function
        } else {
            let edges_in_subgraph = subgraph_id.contains_edges();
            let subgraphs = edges_in_subgraph.iter().map(|e| subgraph_id.pop_edge(e));

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

    fn pop_edge(&self, edge_id: &usize) -> Self {
        Self {
            id: self.id ^ (1 << *edge_id),
            num_edges: self.num_edges,
        }
    }

    fn is_empty(&self) -> bool {
        self.id == 0
    }

    fn has_edge(&self, edge_id: &usize) -> bool {
        self.id & (1 << *edge_id) != 0
    }

    fn contains_edges(&self) -> Vec<usize> {
        (0..self.num_edges).filter(|i| self.has_edge(i)).collect()
    }

    fn from_edge_list(edge_list: &[usize], num_edges: usize) -> Self {
        let mut id = 0;
        for edge_id in edge_list.iter() {
            id |= 1 << edge_id;
        }
        Self { id, num_edges }
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
struct TropicalEdge {
    edge_id: usize, // position in the list of edges
    left: u8,
    right: u8,
    weight: f64,
}

impl TropicalEdge {
    fn contains_vertex(&self, vertex: &u8) -> bool {
        self.left == *vertex || self.right == *vertex
    }
}

// helper struct for generation
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

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TropicalSubgraphTableEntry {
    pub loop_number: u8,
    pub mass_momentum_spanning: bool,
    pub j_function: f64,
    pub generalized_dod: f64,
}

#[derive(Debug, Clone)]
pub struct TropicalSubgraphTable {
    pub table: Vec<TropicalSubgraphTableEntry>,
}

impl TropicalSubgraphTable {
    pub fn generate_from_graph(
        graph: &Graph,
        tropical_edge_weights: &[f64],
    ) -> Result<Self, Report> {
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
                    weight_sum - loop_number as f64 * 3.0 / 2.0 - tropical_graph.dod
                } else {
                    weight_sum - loop_number as f64 * 3.0 / 2.0
                }
            } else {
                1.0
            };

            if generalized_dod <= 0.0 && !subgraph.is_empty() && subgraph != full_subgraph_id {
                return Err(Report::msg(format!(
                    "Generalized DoD: {} is negative for subgraph {:?}\n
                    loop number: {}, mass-momentum spanning: {}, weight sum: {}",
                    generalized_dod, subgraph, loop_number, is_mass_momentum_spanning, weight_sum
                )));
            }

            option_subgraph_table[subgraph.id].loop_number = Some(loop_number);
            option_subgraph_table[subgraph.id].generalized_dod = Some(generalized_dod);
        }

        TropicalGraph::recursive_j_function_eval(&full_subgraph_id, &mut option_subgraph_table);

        Ok(TropicalSubgraphTable {
            table: option_subgraph_table
                .into_iter()
                .map(|e| e.to_entry())
                .collect(),
        })
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
            assert_approx_eq(table.j_function, 3.03030303030303, TOLERANCE);
        }

        let final_table = subgraph_table.table[7];
        assert!(final_table.mass_momentum_spanning);
        assert_eq!(final_table.loop_number, 1);
        assert_approx_eq(final_table.generalized_dod, 0.0, TOLERANCE);
        assert_approx_eq(final_table.j_function, 10.82251082251082, TOLERANCE);
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
            assert_approx_eq(table.j_function, 1.19047619047619, TOLERANCE);
        }

        let final_table = subgraph_table.table[7];
        assert!(final_table.mass_momentum_spanning);
        assert_eq!(final_table.loop_number, 2);
        assert_approx_eq(final_table.generalized_dod, 0.0, TOLERANCE);
        assert_approx_eq(final_table.j_function, 4.251700680272108, TOLERANCE);
    }
}
