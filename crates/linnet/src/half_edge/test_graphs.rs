use nodestore::NodeStorageVec;

use super::*;
use crate::half_edge::builder::HedgeGraphBuilder;

#[test]
fn test_petersen() -> TestResult {
    TestGraph::Petersen.test_all()
}

#[test]
fn test_mobius_ladder_6() -> TestResult {
    println!("{}", TestGraph::MoebiusLadder(6).build().0.base_dot());
    TestGraph::MoebiusLadder(6).test_all()
}

#[test]
fn test_complete_5() -> TestResult {
    println!("{}", TestGraph::Complete(5).build().0.base_dot());
    let complete = TestGraph::Complete(5).build().0;
    let trees = complete.all_spanning_forests_of(&complete.full_filter());

    // for t in &trees {
    //     println!("{}", complete.dot(t))
    // }

    assert_eq!(trees.len(), 5usize.pow(5 - 2));

    TestGraph::Complete(5).test_all()
}

#[test]
fn test_flower_4() -> TestResult {
    println!("{}", TestGraph::Flower(4).build().0.base_dot());
    TestGraph::Flower(4).test_all()
}

#[test]
fn test_cycle_5() -> TestResult {
    println!("{}", TestGraph::Cycle(5).build().0.base_dot());

    let cycle = TestGraph::Cycle(5).build();
    let trees = cycle.0.all_spanning_forests_of(&cycle.0.full_filter());

    // for t in &trees {
    //     println!("{}", cycle.0.dot(t))
    // }

    assert_eq!(trees.len(), 5);
    TestGraph::Cycle(5).test_all()
}

// #[test]
// fn test_complete_bipartite_3_4() -> TestResult {
//     TestGraph::CompleteBipartite(3, 4).test_all()
// }

#[test]
fn test_coxeter() -> TestResult {
    println!("{}", TestGraph::Coxeter.build().0.base_dot());
    TestGraph::Coxeter.test_all()
}

// #[test]
// fn test_prism_5() -> TestResult {
//     TestGraph::PrismGraph(5).test_all()
// }

// #[test]
// fn test_generalized_petersen_5_2() -> TestResult {
//     TestGraph::GeneralizedPetersen(5, 2).test_all()
// }

// Test Graph Structures
#[derive(Debug, Clone)]
#[allow(dead_code)]
enum TestGraph {
    Petersen,
    MoebiusLadder(usize),              // n rungs
    Complete(usize),                   // Kn
    CompleteBipartite(usize, usize),   // Km,n
    PrismGraph(usize),                 // n-prism
    GeneralizedPetersen(usize, usize), // GP(n,k)
    Flower(usize),                     // n petals
    Cycle(usize),                      // cycle of length n
    Heawood,
    MobiusKantor,
    Folkman,
    Coxeter,
    Leech,
    Tutte,
}

impl fmt::Display for TestGraph {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            TestGraph::Heawood => write!(f, "Heawood"),
            TestGraph::MobiusKantor => write!(f, "MobiusKantor"),
            TestGraph::Folkman => write!(f, "Folkman"),
            TestGraph::Coxeter => write!(f, "Coxeter"),
            TestGraph::Leech => write!(f, "Leech"),
            TestGraph::Tutte => write!(f, "Tutte"),
            TestGraph::Petersen => write!(f, "Petersen"),
            TestGraph::MoebiusLadder(n) => write!(f, "MoebiusLadder({n})"),
            TestGraph::Complete(n) => write!(f, "Complete({n})"),
            TestGraph::Flower(n) => write!(f, "Flower({n})"),
            TestGraph::Cycle(n) => write!(f, "Cycle({n})"),
            TestGraph::CompleteBipartite(n, m) => write!(f, "CompleteBipartite({n}, {m})"),
            TestGraph::PrismGraph(n) => write!(f, "PrismGraph({n})"),
            TestGraph::GeneralizedPetersen(n, k) => write!(f, "GeneralizedPetersen({n}, {k})"),
        }
    }
}

impl TestGraph {
    fn build(&self) -> (HedgeGraph<(), (), (), NodeStorageVec<()>>, GraphProperties) {
        let mut builder = HedgeGraphBuilder::new();

        match self {
            TestGraph::Heawood => {
                // We'll label 7 "points" as 0..6, 7 "lines" as 7..13
                let points: Vec<_> = (0..7).map(|_| builder.add_node(())).collect();
                let lines: Vec<_> = (0..7).map(|_| builder.add_node(())).collect();

                // The Fano-plane incidence: each "line" node connects to 3 distinct "point" nodes
                let line_incidence = [
                    (7, [0, 1, 3]),
                    (8, [0, 2, 6]),
                    (9, [0, 4, 5]),
                    (10, [1, 2, 4]),
                    (11, [1, 5, 6]),
                    (12, [2, 3, 5]),
                    (13, [3, 4, 6]),
                ];
                for (l, pts) in line_incidence {
                    let ln = lines[l - 7]; // lines go 7..13 → index 0..6
                    for &p in pts.iter() {
                        builder.add_edge(ln, points[p], (), false);
                    }
                }

                let graph = builder.build();
                let properties = GraphProperties {
                    n_nodes: 14,
                    n_edges: 21,
                    connected_components: 1,
                    min_degree: 3,
                    max_degree: 3,
                    girth: Some(6),
                    diameter: 3,
                    is_bipartite: true,
                    // Cyclomatic: E - V + 1 = 21 - 14 + 1 = 8
                    cyclotomatic_number: 8,
                    n_cycles: 8, // if you're using that to store E - V + 1
                };
                (graph, properties)
            }

            // ----------------------------------------------------------------
            // 2) Möbius–Kantor Graph
            //    - 16 vertices, 24 edges
            //    - 3-regular, bipartite, girth=6, diameter=4
            // ----------------------------------------------------------------
            TestGraph::MobiusKantor => {
                // We’ll use a standard labeling 0..15 with this adjacency list.
                // Each entry i → adjacency[i] is a triple (3-regular).
                let adjacency = vec![
                    vec![1, 3, 8],   // 0
                    vec![0, 2, 9],   // 1
                    vec![1, 3, 10],  // 2
                    vec![0, 2, 11],  // 3
                    vec![5, 7, 12],  // 4
                    vec![4, 6, 13],  // 5
                    vec![5, 7, 14],  // 6
                    vec![4, 6, 15],  // 7
                    vec![0, 10, 15], // 8
                    vec![1, 11, 12], // 9
                    vec![2, 8, 13],  // 10
                    vec![3, 9, 14],  // 11
                    vec![4, 9, 15],  // 12
                    vec![5, 10, 14], // 13
                    vec![6, 11, 13], // 14
                    vec![7, 8, 12],  // 15
                ];
                // Build edges
                let node_idx: Vec<_> = adjacency.iter().map(|_| builder.add_node(())).collect();

                // Now add edges
                for (i, nbrs) in adjacency.iter().enumerate() {
                    for &j in nbrs {
                        // add each edge only once, say if i<j
                        if i < j {
                            builder.add_edge(node_idx[i], node_idx[j], (), false);
                        }
                    }
                }

                let graph = builder.build();
                let properties = GraphProperties {
                    n_nodes: 16,
                    n_edges: 24,
                    connected_components: 1,
                    min_degree: 3,
                    max_degree: 3,
                    girth: Some(6),
                    diameter: 4,        // known property for Möbius–Kantor
                    is_bipartite: true, // it is bipartite
                    // E - V + 1 = 24 - 16 + 1 = 9
                    cyclotomatic_number: 9,
                    n_cycles: 9,
                };
                (graph, properties)
            }

            // ----------------------------------------------------------------
            // 3) Folkman Graph
            //    - 20 vertices, 40 edges
            //    - 4-regular, girth=4, diameter=2, not bipartite
            // ----------------------------------------------------------------
            TestGraph::Folkman => {
                // Adjacency from a standard labeling 0..19
                // Each line is the list of neighbors of vertex i
                let adjacency = vec![
                    vec![1, 2, 4, 5],    // 0
                    vec![0, 2, 6, 7],    // 1
                    vec![0, 1, 8, 9],    // 2
                    vec![4, 6, 8, 10],   // 3
                    vec![0, 3, 11, 12],  // 4
                    vec![0, 13, 14, 15], // 5
                    vec![1, 3, 16, 17],  // 6
                    vec![1, 8, 14, 16],  // 7
                    vec![2, 3, 7, 9],    // 8
                    vec![2, 8, 18, 19],  // 9
                    vec![3, 12, 13, 16], // 10
                    vec![4, 13, 14, 17], // 11
                    vec![4, 10, 15, 18], // 12
                    vec![5, 10, 11, 19], // 13
                    vec![5, 7, 11, 15],  // 14
                    vec![5, 12, 14, 19], // 15
                    vec![6, 7, 10, 17],  // 16
                    vec![6, 11, 16, 18], // 17
                    vec![9, 12, 17, 19], // 18
                    vec![9, 13, 15, 18], // 19
                ];
                let node_idx: Vec<_> = adjacency.iter().map(|_| builder.add_node(())).collect();

                // Add edges
                for (i, nbrs) in adjacency.iter().enumerate() {
                    for &j in nbrs {
                        if i < j {
                            builder.add_edge(node_idx[i], node_idx[j], (), false);
                        }
                    }
                }

                let graph = builder.build();
                let properties = GraphProperties {
                    n_nodes: 20,
                    n_edges: 40,
                    connected_components: 1,
                    min_degree: 4,
                    max_degree: 4,
                    girth: Some(4),
                    diameter: 2,
                    is_bipartite: false, // it has odd cycles (just no triangles)
                    // E - V + 1 = 40 - 20 + 1 = 21
                    cyclotomatic_number: 21,
                    n_cycles: 21,
                };
                (graph, properties)
            }

            // ----------------------------------------------------------------
            // 4) Coxeter Graph
            //    - 28 vertices, 42 edges
            //    - 3-regular, bipartite, girth=4, diameter=4
            // ----------------------------------------------------------------
            TestGraph::Coxeter => {
                // A standard 0..27 labeling with 3 neighbors each
                // (Distance-transitive, 3-regular, bipartite, girth=4, diameter=4).
                let adjacency = vec![
                    vec![1, 2, 3],    // 0
                    vec![0, 4, 5],    // 1
                    vec![0, 6, 7],    // 2
                    vec![0, 8, 9],    // 3
                    vec![1, 6, 10],   // 4
                    vec![1, 8, 11],   // 5
                    vec![2, 4, 12],   // 6
                    vec![2, 9, 13],   // 7
                    vec![3, 5, 14],   // 8
                    vec![3, 7, 15],   // 9
                    vec![4, 11, 16],  // 10
                    vec![5, 10, 17],  // 11
                    vec![6, 13, 18],  // 12
                    vec![7, 12, 19],  // 13
                    vec![8, 15, 20],  // 14
                    vec![9, 14, 21],  // 15
                    vec![10, 17, 22], // 16
                    vec![11, 16, 23], // 17
                    vec![12, 19, 24], // 18
                    vec![13, 18, 25], // 19
                    vec![14, 21, 26], // 20
                    vec![15, 20, 27], // 21
                    vec![16, 23, 26], // 22
                    vec![17, 22, 27], // 23
                    vec![18, 25, 26], // 24
                    vec![19, 24, 27], // 25
                    vec![20, 22, 24], // 26
                    vec![21, 23, 25], // 27
                ];
                let node_idx: Vec<_> = adjacency.iter().map(|_| builder.add_node(())).collect();
                for (i, nbrs) in adjacency.iter().enumerate() {
                    for &j in nbrs {
                        if i < j {
                            builder.add_edge(node_idx[i], node_idx[j], (), false);
                        }
                    }
                }

                let graph = builder.build();
                let properties = GraphProperties {
                    n_nodes: 28,
                    n_edges: 42,
                    connected_components: 1,
                    min_degree: 3,
                    max_degree: 3,
                    girth: Some(4),
                    diameter: 4,
                    is_bipartite: true,
                    // E - V + 1 = 42 - 28 + 1 = 15
                    cyclotomatic_number: 15,
                    n_cycles: 15,
                };
                (graph, properties)
            }

            TestGraph::Petersen => {
                // Build Petersen Graph
                let outer: Vec<_> = (0..5).map(|_| builder.add_node(())).collect();
                let inner: Vec<_> = (0..5).map(|_| builder.add_node(())).collect();

                // Outer pentagon
                for i in 0..5 {
                    builder.add_edge(outer[i], outer[(i + 1) % 5], (), false);
                }

                // Inner pentagram
                for i in 0..5 {
                    builder.add_edge(inner[i], inner[(i + 2) % 5], (), false);
                }

                // Spokes
                for i in 0..5 {
                    builder.add_edge(outer[i], inner[i], (), false);
                }

                let properties = GraphProperties {
                    n_nodes: 10,
                    n_edges: 15,
                    n_cycles: 6,
                    connected_components: 1,
                    min_degree: 3,
                    max_degree: 3,
                    girth: Some(5),
                    diameter: 2,
                    is_bipartite: false,
                    cyclotomatic_number: 6,
                };

                (builder.build(), properties)
            }
            TestGraph::MoebiusLadder(n) => {
                // Need at least 3 vertices per side for a valid Möbius ladder
                assert!(*n >= 3, "Möbius ladder needs at least 3 vertices per side");
                assert!(*n >= 3, "Möbius ladder needs at least 3 rungs");

                // Create 2n nodes
                let vs: Vec<_> = (0..2 * n).map(|_| builder.add_node(())).collect();

                // 1) Add the single cycle of length 2n
                for i in 0..(2 * n) {
                    builder.add_edge(vs[i], vs[(i + 1) % (2 * n)], (), false);
                }

                // 2) Add the n "rungs": i -- i+n
                for i in 0..*n {
                    builder.add_edge(vs[i], vs[i + n], (), false);
                }

                let properties = GraphProperties {
                    n_nodes: 2 * n,             // n vertices per side [1]
                    n_edges: 3 * n,             // n rungs + 2n side edges [1]
                    n_cycles: n + 1,            // n squares + 1 large cycle [2]
                    connected_components: 1,    // Graph is connected [1]
                    min_degree: 3,              // Each vertex has exactly 3 edges [1]
                    max_degree: 3,              // Cubic graph [1]
                    girth: Some(4),             // Smallest cycle is a square [3]
                    diameter: n / 2,            // Maximum distance between any two vertices [4]
                    is_bipartite: n % 2 == 0,   // Bipartite only when n is even [5]
                    cyclotomatic_number: n + 1, // |E| - |V| + 1 = 3n - 2n + 1 = n + 1 [6]
                };

                (builder.build(), properties)
            }
            TestGraph::Complete(n) => {
                let nodes: Vec<_> = (0..*n).map(|_| builder.add_node(())).collect();

                for i in 0..*n {
                    for j in (i + 1)..*n {
                        builder.add_edge(nodes[i], nodes[j], (), false);
                    }
                }

                let properties = GraphProperties {
                    n_nodes: *n,
                    n_edges: n * (n - 1) / 2,
                    n_cycles: if *n >= 3 {
                        n * (n - 1) * (n - 2) / 6
                    } else {
                        0
                    },
                    connected_components: 1,
                    min_degree: n - 1,
                    max_degree: n - 1,
                    girth: if *n >= 3 { Some(3) } else { None },
                    diameter: 1,
                    is_bipartite: *n <= 2,
                    cyclotomatic_number: (n - 1) * (n - 2) / 2,
                };

                (builder.build(), properties)
            }
            TestGraph::Flower(n) => {
                let center = builder.add_node(());
                let petals: Vec<_> = (0..*n).map(|_| builder.add_node(())).collect();

                for i in 0..*n {
                    builder.add_edge(center, petals[i], (), false);
                    builder.add_edge(petals[i], petals[(i + 1) % n], (), false);
                }

                let properties = GraphProperties {
                    n_nodes: n + 1,
                    n_edges: 2 * n,
                    n_cycles: *n,
                    connected_components: 1,
                    min_degree: 3,
                    max_degree: *n,
                    girth: Some(3),
                    diameter: 2,
                    is_bipartite: n % 2 == 0,
                    cyclotomatic_number: *n,
                };

                (builder.build(), properties)
            }

            TestGraph::Cycle(n) => {
                // Create n nodes in a cycle
                let nodes: Vec<_> = (0..*n).map(|_| builder.add_node(())).collect();

                // Connect nodes in a cycle: 0->1->2->...->n-1->0
                for i in 0..*n {
                    let j = (i + 1) % *n;
                    builder.add_edge(nodes[i], nodes[j], (), false);
                }

                let properties = GraphProperties {
                    n_nodes: *n,
                    n_edges: *n,
                    n_cycles: if *n >= 3 { 1 } else { 0 },
                    connected_components: 1,
                    min_degree: 2,
                    max_degree: 2,
                    girth: if *n >= 3 { Some(*n) } else { None },
                    diameter: *n / 2,
                    is_bipartite: *n % 2 == 0,
                    cyclotomatic_number: 1,
                };

                (builder.build(), properties)
            }
            _ => unimplemented!(),
        }
    }
}

use std::fmt;

type TestResult = Result<(), TestError>;

#[derive(Debug)]
#[allow(dead_code)]
enum TestError {
    WrongNodeCount {
        expected: usize,
        found: usize,
    },
    WrongEdgeCount {
        expected: usize,
        found: usize,
    },
    WrongComponentCount {
        expected: usize,
        found: usize,
    },
    WrongCyclomaticNumber {
        expected: usize,
        found: usize,
    },
    WrongDegree {
        kind: &'static str,
        expected: usize,
        found: usize,
    },
    InvalidCycle,
    Custom(String),
}

impl std::fmt::Display for TestError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            TestError::WrongNodeCount { expected, found } => write!(
                f,
                "Wrong number of nodes. Expected: {expected}, Found: {found}",
            ),
            TestError::WrongEdgeCount { expected, found } => write!(
                f,
                "Wrong number of edges. Expected: {expected}, Found: {found}",
            ),
            TestError::WrongComponentCount { expected, found } => write!(
                f,
                "Wrong number of components. Expected: {expected}, Found: {found}",
            ),
            TestError::WrongCyclomaticNumber { expected, found } => write!(
                f,
                "Wrong cyclotomatic number. Expected: {expected}, Found: {found}",
            ),
            TestError::WrongDegree {
                kind,
                expected,
                found,
            } => write!(
                f,
                "Wrong {kind} degree. Expected: {expected}, Found: {found}",
            ),
            TestError::InvalidCycle => write!(f, "Invalid cycle detected"),
            TestError::Custom(msg) => write!(f, "{msg}"),
        }
    }
}

// #[derive(Debug, Clone)]
// enum TestGraph {
//     Petersen,
//     MoebiusLadder(usize),
//     Complete(usize),
//     Flower(usize),
// }

#[derive(Debug)]
#[allow(dead_code)]
struct GraphProperties {
    n_nodes: usize,
    n_edges: usize,
    n_cycles: usize,
    connected_components: usize,
    min_degree: usize,
    max_degree: usize,
    girth: Option<usize>,
    diameter: usize,
    is_bipartite: bool,
    cyclotomatic_number: usize,
}

impl TestGraph {
    fn test_basic_properties(
        &self,
        graph: &HedgeGraph<(), (), (), NodeStorageVec<()>>,
        props: &GraphProperties,
    ) -> TestResult {
        if graph.n_nodes() != props.n_nodes {
            return Err(TestError::WrongNodeCount {
                expected: props.n_nodes,
                found: graph.n_nodes(),
            });
        }

        if graph.n_internals() != props.n_edges {
            return Err(TestError::WrongEdgeCount {
                expected: props.n_edges,
                found: graph.n_internals(),
            });
        }

        let components = graph.count_connected_components(&graph.full_filter());
        if components != props.connected_components {
            return Err(TestError::WrongComponentCount {
                expected: props.connected_components,
                found: components,
            });
        }

        Ok(())
    }

    fn test_cycle_properties(
        &self,
        graph: &HedgeGraph<(), (), (), NodeStorageVec<()>>,
        props: &GraphProperties,
    ) -> TestResult {
        let cyclomatic = graph.cyclotomatic_number(&graph.full_graph());
        if cyclomatic != props.cyclotomatic_number {
            return Err(TestError::WrongCyclomaticNumber {
                expected: props.cyclotomatic_number,
                found: cyclomatic,
            });
        }

        let cycles = graph.cycle_basis().0;
        for cycle in &cycles {
            if !Self::verify_cycle(graph, cycle) {
                return Err(TestError::InvalidCycle);
            }
        }

        Ok(())
    }

    fn test_degree_distribution(
        &self,
        graph: &HedgeGraph<(), (), (), NodeStorageVec<()>>,
        props: &GraphProperties,
    ) -> TestResult {
        let mut min_degree = usize::MAX;
        let mut max_degree = 0;

        for i in 0..graph.n_nodes() {
            let degree = graph.iter_crown(NodeIndex(i)).count();
            min_degree = min_degree.min(degree);
            max_degree = max_degree.max(degree);
        }

        if min_degree != props.min_degree {
            return Err(TestError::WrongDegree {
                kind: "minimum",
                expected: props.min_degree,
                found: min_degree,
            });
        }

        if max_degree != props.max_degree {
            return Err(TestError::WrongDegree {
                kind: "maximum",
                expected: props.max_degree,
                found: max_degree,
            });
        }

        Ok(())
    }

    fn verify_cycle(graph: &HedgeGraph<(), (), (), NodeStorageVec<()>>, cycle: &Cycle) -> bool {
        let mut degree_map = std::collections::HashMap::new();
        for hedge in cycle.filter.included_iter() {
            let node = graph.node_id(hedge);
            *degree_map.entry(node).or_insert(0) += 1;
        }
        degree_map.values().all(|&degree| degree % 2 == 0)
    }

    fn test_all(&self) -> TestResult {
        let (graph, properties) = self.build();

        self.test_basic_properties(&graph, &properties)?;
        self.test_cycle_properties(&graph, &properties)?;
        self.test_degree_distribution(&graph, &properties)?;

        Ok(())
    }
}
