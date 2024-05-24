use crate::{
    cff::{
        expression::OrientationExpression,
        hsurface::HsurfaceID,
        surface::{HybridSurface, HybridSurfaceID},
        tree::Tree,
    },
    graph::Graph,
};
use ahash::HashMap;
use color_eyre::Report;
use eyre::Result;
use itertools::Itertools;
use std::fmt::Debug;

use serde::{Deserialize, Serialize};

use log::debug;

use super::{
    cff_graph::{CFFGenerationGraph, VertexSet},
    esurface::{EsurfaceCollection, EsurfaceID},
    expression::{CFFExpression, CFFExpressionNode, CFFLimit, TermId},
    hsurface::HsurfaceCollection,
    tree::NodeId,
};

#[derive(Debug, Clone)]
enum GenerationData {
    Data {
        graph: CFFGenerationGraph,
        surface_id: Option<HybridSurfaceID>,
    },
    Pointer {
        term_id: usize,
        node_id: usize,
    },
}

fn forget_graphs(data: GenerationData) -> CFFExpressionNode {
    match data {
        GenerationData::Data {
            surface_id: esurface_id,
            ..
        } => CFFExpressionNode::Data(esurface_id.unwrap()),
        GenerationData::Pointer { term_id, node_id } => CFFExpressionNode::Pointer {
            term_id: Into::<TermId>::into(term_id),
            node_id: Into::<NodeId>::into(node_id),
        },
    }
}

impl GenerationData {
    fn insert_esurface(&mut self, surface_id: HybridSurfaceID) {
        match self {
            GenerationData::Data {
                surface_id: ref mut id,
                ..
            } => {
                *id = Some(surface_id);
            }
            GenerationData::Pointer { .. } => {}
        }
    }
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
struct CFFTreeNodePointer {
    term_id: usize,
    node_id: usize,
}

// Orientation is a bitset that represents the orientation of the edges in a graph.
// 0 means +, 1 means -, in this representation the original graph is represented by the number 0
#[derive(Debug, Clone, Copy)]
struct OrientationGenerator {
    identifier: usize,
    num_edges: usize,
}

impl OrientationGenerator {
    #[allow(unused)]
    fn default(num_edges: usize) -> Self {
        Self {
            identifier: 0,
            num_edges,
        }
    }
}

impl IntoIterator for OrientationGenerator {
    type Item = bool;
    type IntoIter = OrientationIterator;

    fn into_iter(self) -> Self::IntoIter {
        OrientationIterator {
            identifier: self.identifier,
            current_location: 0,
            num_edges: self.num_edges,
        }
    }
}

// OrientationIterator allows us to iterate over the edges in a graph, and
// view their orientation as a boolean
struct OrientationIterator {
    identifier: usize,
    current_location: usize,
    num_edges: usize,
}

impl Iterator for OrientationIterator {
    type Item = bool;
    fn next(&mut self) -> Option<Self::Item> {
        if self.current_location < self.num_edges {
            let result = self.identifier & (1 << self.current_location) == 0;
            self.current_location += 1;
            Some(result)
        } else {
            None
        }
    }
}

// This function returns an iterator over all possible orientations of a graph
fn iterate_possible_orientations(num_edges: usize) -> impl Iterator<Item = OrientationGenerator> {
    if num_edges > 64 {
        panic!("Maximum number of edges supported is currently 64")
    }

    let max_size = 2_usize.pow(num_edges as u32);
    (0..max_size).map(move |x| OrientationGenerator {
        identifier: x,
        num_edges,
    })
}

fn get_orientations(graph: &Graph) -> Vec<CFFGenerationGraph> {
    let num_virtual_edges = graph.get_virtual_edges_iterator().count();
    let possible_orientations = iterate_possible_orientations(num_virtual_edges);

    possible_orientations
        .map(|orientation_of_virtuals| {
            let orientation_of_virtuals = orientation_of_virtuals.into_iter().collect_vec();

            CFFGenerationGraph::new(graph, orientation_of_virtuals)
        })
        .collect()
}

pub fn generate_cff_expression(
    graph: &Graph,
    allow_hsurface: bool,
) -> Result<CFFExpression, Report> {
    // construct a hashmap that contains as keys all vertices that connect to external edges
    // and as values those external edges that it connects to

    let graphs = get_orientations(graph);
    debug!("generating cff for graph: {}", graph.name);
    debug!("number of orientations: {}", graphs.len());

    let graph_cff = generate_cff_from_orientations(graphs, None, None, allow_hsurface)?;

    Ok(graph_cff)
}

pub fn generate_cff_limit(
    left_dags: Vec<CFFGenerationGraph>,
    right_dags: Vec<CFFGenerationGraph>,
    esurfaces: &EsurfaceCollection,
) -> CFFLimit {
    assert_eq!(
        left_dags.len(),
        right_dags.len(),
        "number of left and right dags must match"
    );

    let left =
        generate_cff_from_orientations(left_dags, Some(esurfaces.clone()), None, true).unwrap();
    assert_eq!(
        left.esurfaces.len(),
        esurfaces.len(),
        "new esurfaces generated during factorisation"
    );

    let right =
        generate_cff_from_orientations(right_dags, Some(esurfaces.clone()), None, true).unwrap();
    assert_eq!(
        right.esurfaces.len(),
        esurfaces.len(),
        "new esurfaces generated during factorisation"
    );

    CFFLimit { left, right }
}

fn generate_cff_from_orientations(
    orientations_and_graphs: Vec<CFFGenerationGraph>,
    optional_esurface_cache: Option<EsurfaceCollection>,
    optional_hsurface_cache: Option<HsurfaceCollection>,
    allow_hsurface: bool,
) -> Result<CFFExpression, Report> {
    let esurface_cache = if let Some(cache) = optional_esurface_cache {
        cache
    } else {
        EsurfaceCollection::from_iter(std::iter::empty())
    };

    let hsurface_cache = if let Some(cache) = optional_hsurface_cache {
        cache
    } else {
        HsurfaceCollection::from_iter(std::iter::empty())
    };

    let graph_cache = HashMap::default();

    let mut generator_cache = GeneratorCache {
        graph_cache,
        esurface_cache,
        hsurface_cache,
        vertices_used: vec![],
        cache_hits: 0,
        non_cache_hits: 0,
    };

    // filter cyclic orientations beforehand
    let acyclic_orientations_and_graphs = orientations_and_graphs
        .into_iter()
        .filter(|graph| !graph.has_directed_cycle_initial())
        .collect_vec();

    debug!(
        "number of acyclic orientations: {}",
        acyclic_orientations_and_graphs.len()
    );

    let terms = acyclic_orientations_and_graphs
        .into_iter()
        .enumerate()
        .map(|(term_id, graph)| {
            let global_orientation = graph.global_orientation.clone();
            let tree = generate_tree_for_orientation(
                graph.clone(),
                term_id,
                &mut generator_cache,
                allow_hsurface,
            );
            let expression = tree.map(forget_graphs);

            OrientationExpression {
                expression,
                orientation: global_orientation,
                dag: graph,
            }
        })
        .collect_vec();

    debug!("number of cache hits: {}", generator_cache.cache_hits);
    debug!(
        "percentage of cache hits: {:.1}%",
        generator_cache.cache_hits as f64
            / (generator_cache.cache_hits + generator_cache.non_cache_hits) as f64
            * 100.0
    );

    #[cfg(test)]
    {
        println!("number of cache hits: {}", generator_cache.cache_hits);
        println!(
            "percentage of cache hits: {:.1}%",
            generator_cache.cache_hits as f64
                / (generator_cache.cache_hits + generator_cache.non_cache_hits) as f64
                * 100.0
        );
    }

    Ok(CFFExpression {
        orientations: terms.into(),
        esurfaces: generator_cache.esurface_cache,
        hsurfaces: generator_cache.hsurface_cache,
    })
}

struct GeneratorCache {
    graph_cache: HashMap<CFFGenerationGraph, (usize, usize)>,
    esurface_cache: EsurfaceCollection,
    hsurface_cache: HsurfaceCollection,
    vertices_used: Vec<VertexSet>,
    cache_hits: usize,
    non_cache_hits: usize,
}

fn generate_tree_for_orientation(
    graph: CFFGenerationGraph,
    term_id: usize,
    generator_cache: &mut GeneratorCache,
    allow_hsurface: bool,
) -> Tree<GenerationData> {
    let mut tree = Tree::from_root(GenerationData::Data {
        graph,
        surface_id: None,
    });

    while let Some(()) = advance_tree(&mut tree, term_id, generator_cache, allow_hsurface) {}

    tree
}

fn advance_tree(
    tree: &mut Tree<GenerationData>,
    term_id: usize,
    generator_cache: &mut GeneratorCache,
    allow_hsurface: bool,
) -> Option<()> {
    let bottom_layer = tree
        .get_bottom_layer()
        .into_iter()
        .filter(|&node_id| matches!(&tree.get_node(node_id).data, GenerationData::Data { .. }))
        .collect_vec(); // allocation needed because tree is mutable

    let (children_optional, new_surfaces_for_tree): (
        Vec<Option<Vec<CFFGenerationGraph>>>,
        Vec<HybridSurfaceID>,
    ) = bottom_layer
        .iter()
        .map(|&node_id| {
            let node = &tree.get_node(node_id);
            let graph = match &node.data {
                GenerationData::Data { graph, .. } => graph,
                GenerationData::Pointer { .. } => {
                    unreachable!("filtered")
                }
            };

            let (option_children, surface) =
                graph.generate_children(&mut generator_cache.vertices_used, allow_hsurface);

            let surface_id = match surface {
                HybridSurface::Esurface(esurface) => {
                    let option_esurface_id = generator_cache
                        .esurface_cache
                        .position(|val| val == &esurface);

                    let esurface_id = match option_esurface_id {
                        Some(esurface_id) => esurface_id,
                        None => {
                            generator_cache.esurface_cache.push(esurface);
                            Into::<EsurfaceID>::into(generator_cache.esurface_cache.len() - 1)
                        }
                    };

                    HybridSurfaceID::Esurface(esurface_id)
                }
                HybridSurface::Hsurface(hsurface) => {
                    let option_hsurface_id = generator_cache
                        .hsurface_cache
                        .position(|val| val == &hsurface);

                    let hsurface_id = match option_hsurface_id {
                        Some(hsurface_id) => hsurface_id,
                        None => {
                            generator_cache.hsurface_cache.push(hsurface);
                            Into::<HsurfaceID>::into(generator_cache.hsurface_cache.len() - 1)
                        }
                    };

                    HybridSurfaceID::Hsurface(hsurface_id)
                }
            };

            (option_children, surface_id)
        })
        .unzip();

    bottom_layer
        .iter()
        .zip(new_surfaces_for_tree)
        .for_each(|(&node_id, esurface_id)| {
            tree.apply_mut_closure(node_id, |data| data.insert_esurface(esurface_id))
        });

    let all_some = children_optional.iter().all(Option::is_some);
    let all_none = children_optional.iter().all(Option::is_none);

    assert!(
        all_some || all_none,
        "Some cff branches have finished earlier than others"
    );

    let children = if all_some && !all_none {
        children_optional
            .into_iter()
            .map(Option::unwrap)
            .collect_vec()
    } else {
        return None;
    };

    bottom_layer
        .iter()
        .zip(children)
        .for_each(|(&node_id, children)| {
            children.into_iter().for_each(|child| {
                let hashable_child = child.clone();
                if let Some((cff_expression_term_id, cff_expression_node_id)) =
                    generator_cache.graph_cache.get(&hashable_child)
                {
                    let new_pointer = GenerationData::Pointer {
                        term_id: *cff_expression_term_id,
                        node_id: *cff_expression_node_id,
                    };
                    tree.insert_node(node_id, new_pointer);
                    generator_cache.cache_hits += 1;
                } else {
                    let child_node_id = tree.get_num_nodes();
                    let child_node = GenerationData::Data {
                        graph: child,
                        surface_id: None,
                    };

                    generator_cache
                        .graph_cache
                        .insert(hashable_child, (term_id, child_node_id));
                    tree.insert_node(node_id, child_node);
                    generator_cache.non_cache_hits += 1;
                }
            });
        });

    Some(())
}

#[cfg(test)]
mod tests_cff {
    use lorentz_vector::LorentzVector;
    use num::traits::Inv;
    use utils::FloatLike;

    use crate::utils;

    use super::*;

    // helper function to do some quick tests
    #[allow(unused)]
    fn generate_orientations_for_testing(
        edges: Vec<(usize, usize)>,
        incoming_vertices: Vec<usize>,
    ) -> Vec<CFFGenerationGraph> {
        let num_edges = edges.len();

        iterate_possible_orientations(num_edges)
            .map(|or| {
                let orientation_vector = or.into_iter().collect_vec();
                let mut new_edges = edges.clone();
                for (edge_id, edge_orientation) in orientation_vector.iter().enumerate() {
                    if *edge_orientation {
                        new_edges[edge_id] = edges[edge_id];
                    } else {
                        let rotated_edge = (edges[edge_id].1, edges[edge_id].0);
                        new_edges[edge_id] = rotated_edge;
                    }
                }

                CFFGenerationGraph::from_vec(
                    new_edges,
                    incoming_vertices.clone(),
                    orientation_vector,
                )
            })
            .filter(|graph| !graph.has_directed_cycle_initial())
            .collect_vec()
    }

    #[allow(unused)]
    fn compute_one_loop_energy<T: FloatLike>(k: LorentzVector<T>, p: LorentzVector<T>, m: T) -> T {
        ((k + p).spatial_squared() + m * m).sqrt()
    }

    #[test]
    fn test_orientation_struct() {
        let orientations = iterate_possible_orientations(3).collect_vec();
        assert_eq!(orientations.len(), 8);

        let orientation1 = orientations[0].into_iter().collect_vec();
        assert_eq!(orientation1, vec![true, true, true]);

        let orientation2 = orientations[1].into_iter().collect_vec();
        assert_eq!(orientation2, vec![false, true, true]);

        let orientation3 = orientations[2].into_iter().collect_vec();
        assert_eq!(orientation3, vec![true, false, true]);

        let orientation4 = orientations[3].into_iter().collect_vec();
        assert_eq!(orientation4, vec![false, false, true]);

        let orientation5 = orientations[4].into_iter().collect_vec();
        assert_eq!(orientation5, vec![true, true, false]);

        let orientation6 = orientations[5].into_iter().collect_vec();
        assert_eq!(orientation6, vec![false, true, false]);

        let orientation7 = orientations[6].into_iter().collect_vec();
        assert_eq!(orientation7, vec![true, false, false]);

        let orientation8 = orientations[7].into_iter().collect_vec();
        assert_eq!(orientation8, vec![false, false, false]);
    }

    #[test]
    fn test_cff_generation_triangle() {
        let triangle = vec![(2, 0), (0, 1), (1, 2)];

        let incoming_vertices = vec![0, 1, 2];
        let orientations = generate_orientations_for_testing(triangle, incoming_vertices);
        assert_eq!(orientations.len(), 6);

        let cff = generate_cff_from_orientations(orientations, None, None, false).unwrap();
        assert_eq!(cff.esurfaces.len(), 6);

        let p1 = LorentzVector::from_args(1., 3., 4., 5.);
        let p2 = LorentzVector::from_args(1., 6., 7., 8.);
        let p3 = -p1 - p2;
        let zero = LorentzVector::from_args(0., 0., 0., 0.);
        let m = 0.;

        let k = LorentzVector::from_args(0., 1., 2., 3.);

        let virtual_energy_cache = [
            compute_one_loop_energy(k, zero, m),
            compute_one_loop_energy(k, p1, m),
            compute_one_loop_energy(k, p1 + p2, m),
        ];

        let external_energy_cache = [p1.t, p2.t, p3.t];

        // combine the virtual and external energies
        let mut energy_cache = external_energy_cache.to_vec();
        energy_cache.extend(virtual_energy_cache);

        let energy_prefactor = virtual_energy_cache
            .iter()
            .map(|e| (2. * e).inv())
            .product::<f64>();

        let cff_res: f64 =
            energy_prefactor * cff.evaluate(&energy_cache) * (2. * std::f64::consts::PI).powi(-3);

        let target_res = 6.333_549_225_536_17e-9_f64;
        let absolute_error: f64 = cff_res - target_res;
        let relative_error = absolute_error.abs() / cff_res.abs();

        assert!(
            relative_error.abs() < 1.0e-15,
            "relative error: {:+e} (ground truth: {:+e} vs reproduced: {:+e})",
            relative_error,
            target_res,
            cff_res
        );
    } //

    #[test]
    fn test_cff_test_double_triangle() {
        let double_triangle_edges = vec![(0, 1), (0, 2), (1, 2), (1, 3), (2, 3)];
        let incoming_vertices = vec![0, 3];

        let orientations =
            generate_orientations_for_testing(double_triangle_edges, incoming_vertices);

        let cff = generate_cff_from_orientations(orientations, None, None, false).unwrap();

        let q = LorentzVector::from_args(1., 2., 3., 4.);
        let zero = LorentzVector::from_args(0., 0., 0., 0.);

        let k = LorentzVector::from_args(3., 6., 23., 9.);
        let l = LorentzVector::from_args(0., 3., 12., 34.);

        let virtual_energy_cache = [
            compute_one_loop_energy(k, zero, 0.),
            compute_one_loop_energy(q - k, zero, 0.),
            compute_one_loop_energy(k - l, zero, 0.),
            compute_one_loop_energy(l, zero, 0.),
            compute_one_loop_energy(q - l, zero, 0.),
        ];

        let external_energy_cache = [q.t, -q.t];

        let mut energy_cache = external_energy_cache.to_vec();
        energy_cache.extend(virtual_energy_cache);

        let energy_prefactor = virtual_energy_cache
            .iter()
            .map(|e| (2. * e).inv())
            .product::<f64>();
        let cff_res = energy_prefactor * cff.evaluate(&energy_cache);

        let target = 1.0794792137096797e-13;
        let absolute_error = cff_res - target;
        let relative_error = absolute_error / cff_res;

        assert!(
            relative_error.abs() < 1.0e-15,
            "relative error: {:+e}, target: {:+e}, result: {:+e}",
            relative_error,
            target,
            cff_res
        );
    }

    #[test]
    fn test_cff_tbt() {
        let tbt_edges = vec![
            (0, 1),
            (2, 0),
            (1, 2),
            (1, 3),
            (2, 4),
            (3, 4),
            (3, 5),
            (5, 4),
        ];

        let incoming_vertices = vec![0, 5];

        let orientataions = generate_orientations_for_testing(tbt_edges, incoming_vertices);
        let cff = generate_cff_from_orientations(orientataions, None, None, false).unwrap();

        let q = LorentzVector::from_args(1.0, 2.0, 3.0, 4.0);
        let zero_vector = LorentzVector::from_args(0., 0., 0., 0.);

        let p0 = q;
        let p5 = -q;

        let k = LorentzVector::from_args(3., 6., 23., 9.);
        let l = LorentzVector::from_args(0., 3., 12., 34.);
        let m = LorentzVector::from_args(0., 7., 24., 1.);

        let mass = 0.;

        let energies_cache = [
            p0.t,
            p5.t,
            compute_one_loop_energy(k, zero_vector, mass),
            compute_one_loop_energy(k - q, zero_vector, mass),
            compute_one_loop_energy(k - l, zero_vector, mass),
            compute_one_loop_energy(l, zero_vector, mass),
            compute_one_loop_energy(q - l, zero_vector, mass),
            compute_one_loop_energy(l - m, zero_vector, mass),
            compute_one_loop_energy(m, zero_vector, mass),
            compute_one_loop_energy(m - q, zero_vector, mass),
        ];

        let virtual_energy_cache = energies_cache[2..].to_vec();

        let energy_prefactor = virtual_energy_cache
            .iter()
            .map(|e| (2. * e).inv())
            .product::<f64>();
        let res = cff.evaluate(&energies_cache) * energy_prefactor;

        let absolute_error = res - 1.2625322619777278e-21;
        let relative_error = absolute_error / res;
        assert!(
            relative_error.abs() < 1.0e-15,
            "relative error: {:+e}",
            relative_error
        );
    }

    #[test]
    #[ignore]
    fn fishnet2b2() {
        let edges = vec![
            (0, 1),
            (1, 2),
            (3, 4),
            (4, 5),
            (6, 7),
            (7, 8),
            (0, 3),
            (1, 4),
            (2, 5),
            (3, 6),
            (4, 7),
            (5, 8),
        ];

        let incoming_vertices = vec![0, 2, 6, 8];

        let orientations = generate_orientations_for_testing(edges, incoming_vertices);

        // get time before cff generation
        let start = std::time::Instant::now();

        let _cff = generate_cff_from_orientations(orientations, None, None, false).unwrap();

        let finish = std::time::Instant::now();
        println!("time to generate cff: {:?}", finish - start);
    }

    #[test]
    #[ignore]
    fn cube() {
        let edges = vec![
            (0, 1),
            (1, 3),
            (3, 2),
            (2, 0),
            (4, 5),
            (5, 7),
            (7, 6),
            (6, 4),
            (0, 4),
            (1, 5),
            (2, 6),
            (3, 7),
        ];

        let mut external_data = HashMap::default();
        for v in 0..8 {
            external_data.insert(v, vec![12 + v]);
        }

        let mut position_map = HashMap::default();
        for i in 0..edges.len() {
            position_map.insert(i, i);
        }

        let incoming_vertices = vec![0, 1, 2, 3, 4, 5, 6, 7];

        let orientations = generate_orientations_for_testing(edges, incoming_vertices);

        // get time before cff generation
        let _start = std::time::Instant::now();

        let _cff = generate_cff_from_orientations(orientations, None, None, false).unwrap();

        let _finish = std::time::Instant::now();
    }

    #[test]
    #[ignore]
    fn fishnet2b3() {
        let edges = vec![
            (0, 1),
            (1, 2),
            (3, 4),
            (4, 5),
            (6, 7),
            (7, 8),
            (9, 10),
            (10, 11),
            (0, 3),
            (3, 6),
            (6, 9),
            (1, 4),
            (4, 7),
            (7, 10),
            (2, 5),
            (5, 8),
            (8, 11),
        ];

        let incoming_vertices = vec![];

        let _energy_cache = [3.0; 17];

        let orientations = generate_orientations_for_testing(edges, incoming_vertices);
        let energy_cache = [3.0; 17];

        let cff = generate_cff_from_orientations(orientations, None, None, false).unwrap();

        let start = std::time::Instant::now();
        for _ in 0..100 {
            let _res = cff.evaluate(&energy_cache);
        }
        let finish = std::time::Instant::now();
        println!("time to evaluate cff: {:?}", (finish - start) / 100);
    }
}
