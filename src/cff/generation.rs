use crate::{
    cff::{
        esurface::add_external_shifts,
        expression::{CompiledCFFExpression, OrientationExpression},
        hsurface::HsurfaceID,
        surface::{HybridSurface, HybridSurfaceID},
        tree::Tree,
    },
    graph::BareGraph,
};
use ahash::HashMap;
use color_eyre::Report;
use color_eyre::Result;
use itertools::Itertools;
use std::fmt::Debug;

use serde::{Deserialize, Serialize};

use log::debug;

use super::{
    cff_graph::{CFFGenerationGraph, VertexSet},
    esurface::{Esurface, EsurfaceCollection, EsurfaceID, ExternalShift},
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

fn get_orientations(graph: &BareGraph) -> Vec<CFFGenerationGraph> {
    let num_virtual_edges = graph.get_virtual_edges_iterator().count();
    let possible_orientations = iterate_possible_orientations(num_virtual_edges);

    possible_orientations
        .map(|orientation_of_virtuals| {
            let orientation_of_virtuals = orientation_of_virtuals.into_iter().collect_vec();

            CFFGenerationGraph::new(graph, orientation_of_virtuals)
        })
        .collect()
}

pub fn generate_cff_expression(graph: &BareGraph) -> Result<CFFExpression, Report> {
    // construct a hashmap that contains as keys all vertices that connect to external edges
    // and as values those external edges that it connects to

    let graphs = get_orientations(graph);
    debug!("generating cff for graph: {}", graph.name);
    debug!("number of orientations: {}", graphs.len());

    let (dep_mom, dep_mom_expr) = graph.get_dep_mom_expr();

    let graph_cff =
        generate_cff_from_orientations(graphs, None, None, None, dep_mom, &dep_mom_expr)?;

    Ok(graph_cff)
}

pub fn generate_cff_limit(
    left_dags: Vec<CFFGenerationGraph>,
    right_dags: Vec<CFFGenerationGraph>,
    esurfaces: &EsurfaceCollection,
    limit_esurface: &Esurface,
    dep_mom: usize,
    dep_mom_expr: &ExternalShift,
    orientations_in_limit: (Vec<Vec<bool>>, Vec<TermId>),
) -> Result<CFFLimit, String> {
    assert_eq!(
        left_dags.len(),
        right_dags.len(),
        "number of left and right dags must match"
    );

    let left = generate_cff_from_orientations(
        left_dags,
        Some(esurfaces.clone()),
        None,
        Some(limit_esurface),
        dep_mom,
        dep_mom_expr,
    )
    .unwrap();
    let right = generate_cff_from_orientations(
        right_dags,
        Some(esurfaces.clone()),
        None,
        Some(limit_esurface),
        dep_mom,
        dep_mom_expr,
    )
    .unwrap();

    Ok(CFFLimit {
        left,
        right,
        orientations_in_limit,
    })
}

fn generate_cff_from_orientations(
    orientations_and_graphs: Vec<CFFGenerationGraph>,
    optional_esurface_cache: Option<EsurfaceCollection>,
    optional_hsurface_cache: Option<HsurfaceCollection>,
    rewrite_at_cache_growth: Option<&Esurface>,
    dep_mom: usize,
    dep_mom_expr: &ExternalShift,
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
                rewrite_at_cache_growth,
                dep_mom,
                dep_mom_expr,
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

    //#[cfg(test)]
    //{
    //    println!("number of cache hits: {}", generator_cache.cache_hits);
    //    println!(
    //        "percentage of cache hits: {:.1}%",
    //        generator_cache.cache_hits as f64
    //            / (generator_cache.cache_hits + generator_cache.non_cache_hits) as f64
    //            * 100.0
    //    );
    //}

    // let terms =vec![terms[0].clone(),terms[1].clone(),terms[2].clone(),terms[3].clone()];

    Ok(CFFExpression {
        orientations: terms.into(),
        esurfaces: generator_cache.esurface_cache,
        hsurfaces: generator_cache.hsurface_cache,
        compiled: CompiledCFFExpression::None,
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
    rewrite_at_cache_growth: Option<&Esurface>,
    dep_mom: usize,
    dep_mom_expr: &ExternalShift,
) -> Tree<GenerationData> {
    let mut tree = Tree::from_root(GenerationData::Data {
        graph,
        surface_id: None,
    });

    while let Some(()) = advance_tree(
        &mut tree,
        term_id,
        generator_cache,
        rewrite_at_cache_growth,
        dep_mom,
        dep_mom_expr,
    ) {}

    tree
}

fn advance_tree(
    tree: &mut Tree<GenerationData>,
    term_id: usize,
    generator_cache: &mut GeneratorCache,
    rewrite_at_cache_growth: Option<&Esurface>,
    dep_mom: usize,
    dep_mom_expr: &ExternalShift,
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
                graph.generate_children(&mut generator_cache.vertices_used);

            let surface_id = match surface {
                HybridSurface::Esurface(mut esurface) => {
                    esurface.canonicalize_shift(dep_mom, dep_mom_expr);
                    let option_esurface_id = generator_cache
                        .esurface_cache
                        .position(|val| *val == esurface);

                    let esurface_id = match option_esurface_id {
                        Some(esurface_id) => esurface_id,
                        None => {
                            if let Some(rewrite_esurface) = rewrite_at_cache_growth {
                                esurface
                                    .energies
                                    .retain(|e| !rewrite_esurface.energies.contains(e));

                                let negative_rewriter_esurface_shift = rewrite_esurface
                                    .external_shift
                                    .iter()
                                    .map(|(index, sign)| (*index, -sign))
                                    .collect_vec();

                                esurface.external_shift = add_external_shifts(
                                    &esurface.external_shift,
                                    &negative_rewriter_esurface_shift,
                                );

                                esurface.canonicalize_shift(dep_mom, dep_mom_expr);
                                let new_option_esurface_id = generator_cache
                                    .esurface_cache
                                    .position(|val| *val == esurface);

                                match new_option_esurface_id {
                                    Some(new_esurface_id) => new_esurface_id,
                                    None => panic!(
                                        "rewriting the esurface did not yield an existing esurface\n
                                        rewritten esurface: {:#?} \n
                                        using {:#?}\n ", esurface, rewrite_esurface,
                                    ),
                                }
                            } else {
                                generator_cache.esurface_cache.push(esurface);
                                Into::<EsurfaceID>::into(generator_cache.esurface_cache.len() - 1)
                            }
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
                HybridSurface::Unit(_) => HybridSurfaceID::Unit,
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
    use symbolica::{
        atom::{Atom, AtomCore},
        domains::float::{NumericalFloatLike, Real},
        id::Pattern,
    };
    use utils::FloatLike;

    use crate::{
        cff::cff_graph::CFFEdgeType,
        momentum::{FourMomentum, ThreeMomentum},
        tests::{self, load_default_settings},
        utils::{self, RefDefault, F},
    };

    use super::*;

    // helper function to do some quick tests
    #[allow(unused)]
    fn generate_orientations_for_testing(
        edges: Vec<(usize, usize)>,
        incoming_vertices: Vec<usize>,
    ) -> Vec<CFFGenerationGraph> {
        let num_edges = edges.len();
        let incoming_vertices = incoming_vertices
            .into_iter()
            .map(|v| (v, CFFEdgeType::External))
            .collect_vec();

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
    fn compute_one_loop_energy<T: FloatLike>(
        k: ThreeMomentum<F<T>>,
        p: ThreeMomentum<F<T>>,
        m: F<T>,
    ) -> F<T> {
        ((k + p).norm_squared() + &m * &m).sqrt()
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

        let dep_mom = 2;
        let dep_mom_expr = vec![(0, -1), (1, -1)];

        let cff =
            generate_cff_from_orientations(orientations, None, None, None, dep_mom, &dep_mom_expr)
                .unwrap();
        assert_eq!(
            cff.esurfaces.len(),
            6,
            "too many esurfaces: {:#?}",
            cff.esurfaces
        );

        let p1 = FourMomentum::from_args(F(1.), F(3.), F(4.), F(5.));
        let p2 = FourMomentum::from_args(F(1.), F(6.), F(7.), F(8.));
        let p3 = -p1 - p2;
        let zero = FourMomentum::from_args(F(0.), F(0.), F(0.), F(0.));
        let m = F(0.);

        let k = ThreeMomentum::new(F(1.), F(2.), F(3.));

        let virtual_energy_cache = [
            compute_one_loop_energy(k, zero.spatial, m),
            compute_one_loop_energy(k, p1.spatial, m),
            compute_one_loop_energy(k, p1.spatial + p2.spatial, m),
        ];

        let external_energy_cache = [p1.temporal.value, p2.temporal.value, p3.temporal.value];

        // combine the virtual and external energies
        let mut energy_cache = external_energy_cache.to_vec();
        energy_cache.extend(virtual_energy_cache);

        let energy_prefactor = virtual_energy_cache
            .iter()
            .map(|e| (F(2.) * e).inv())
            .reduce(|acc, x| acc * x)
            .unwrap();

        let settings = tests::load_default_settings();

        let cff_res: F<f64> = energy_prefactor
            * cff.eager_evaluate(&energy_cache, &settings)
            * F((2. * std::f64::consts::PI).powi(-3));

        let target_res = F(6.333_549_225_536_17e-9_f64);
        let absolute_error = cff_res - target_res;
        let relative_error = absolute_error.abs() / cff_res.abs();

        assert!(
            relative_error.abs() < F(1.0e-15),
            "relative error: {:+e} (ground truth: {:+e} vs reproduced: {:+e})",
            relative_error,
            target_res,
            cff_res
        );

        // test the generation for each possible limit
        for (esurface_id, _) in cff.esurfaces.iter_enumerated() {
            let expanded_limit = cff.expand_limit_to_atom(HybridSurfaceID::Esurface(esurface_id));

            let limit = cff
                .limit_for_esurface(esurface_id, dep_mom, &dep_mom_expr)
                .unwrap();
            let limit_atom = limit.limit_to_atom_with_rewrite(Some(&cff.esurfaces[esurface_id]));

            let p2_atom = Atom::parse("p2").unwrap();
            let rhs = Atom::parse("- p0 - p1").unwrap();

            let p2_pattern = Pattern::Literal(p2_atom);
            let rhs_pattern = Pattern::Literal(rhs);

            let conditions = None;
            let settings = None;

            let atom_limit_0 =
                expanded_limit.replace_all(&p2_pattern, &rhs_pattern, conditions, settings);

            let limit_atom =
                limit_atom.replace_all(&p2_pattern, &rhs_pattern, conditions, settings);

            let diff = (limit_atom - &atom_limit_0).expand();
            assert_eq!(diff, Atom::new());
        }
    } //

    #[test]
    fn test_cff_test_double_triangle() {
        let double_triangle_edges = vec![(0, 1), (0, 2), (1, 2), (1, 3), (2, 3)];
        let incoming_vertices = vec![0, 3];

        let orientations =
            generate_orientations_for_testing(double_triangle_edges, incoming_vertices);

        let dep_mom = 1;
        let dep_mom_expr = vec![(0, -1)];

        let cff =
            generate_cff_from_orientations(orientations, None, None, None, dep_mom, &dep_mom_expr)
                .unwrap();

        let q = FourMomentum::from_args(F(1.), F(2.), F(3.), F(4.));
        let zero = FourMomentum::from_args(F(0.), F(0.), F(0.), F(0.));

        let k = ThreeMomentum::new(F(6.), F(23.), F(9.));
        let l = ThreeMomentum::new(F(3.), F(12.), F(34.));

        let m = F::from_f64(0.);

        let virtual_energy_cache = [
            compute_one_loop_energy(k, zero.spatial, m),
            compute_one_loop_energy(q.spatial - k, zero.spatial, m),
            compute_one_loop_energy(k - l, zero.spatial, m),
            compute_one_loop_energy(l, zero.spatial, m),
            compute_one_loop_energy(q.spatial - l, zero.spatial, m),
        ];

        let external_energy_cache = [q.temporal.value, -q.temporal.value];
        let settings = tests::load_default_settings();

        let mut energy_cache = external_energy_cache.to_vec();
        energy_cache.extend(virtual_energy_cache);

        let energy_prefactor = virtual_energy_cache
            .iter()
            .map(|e| (F(2.) * e).inv())
            .reduce(|acc, x| acc * x)
            .unwrap();

        let cff_res = energy_prefactor * cff.eager_evaluate(&energy_cache, &settings);

        let target = F(1.0794792137096797e-13);
        let absolute_error = cff_res - target;
        let relative_error = absolute_error / cff_res;

        assert!(
            relative_error.abs() < F(1.0e-15),
            "relative error: {:+e}, target: {:+e}, result: {:+e}",
            relative_error,
            target,
            cff_res
        );

        // this does not work yet

        // for (esurface_id, esurface) in cff.esurfaces.iter_enumerated() {

        // let expanded_limit: RationalPolynomial<IntegerRing, u8> = cff.expand_limit_to_atom(HybridSurfaceID::Esurface(esurface_id)).to_rational_polynomial(&Z, &Z, None);

        // let factorised_limit = cff.limit_for_esurface(esurface_id, dep_mom, &dep_mom_expr).unwrap();
        // let factorised_limit_atom = factorised_limit.limit_to_atom_with_rewrite(Some(esurface)).to_rational_polynomial(&Z, &Z, None);

        // apply energy conservation
        // let diff = expanded_limit - factorised_limit_atom;
        // println!("diff: {}", diff);
        // can't test all, but probably works?
        // symbolica crash, probably works on newer version? can't change because everything is outdated, need to merge with main
        //}
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

        let dep_mom = 1;
        let dep_mom_expr = vec![(0, -1)];

        let orientataions = generate_orientations_for_testing(tbt_edges, incoming_vertices);
        let cff =
            generate_cff_from_orientations(orientataions, None, None, None, dep_mom, &dep_mom_expr)
                .unwrap();

        let q = FourMomentum::from_args(F(1.0), F(2.0), F(3.0), F(4.0));
        let zero_vector = q.default();

        let p0 = q;
        let p5 = -q;

        let k = ThreeMomentum::new(F(6.), F(23.), F(9.));
        let l = ThreeMomentum::new(F(3.), F(12.), F(34.));
        let m = ThreeMomentum::new(F(7.), F(24.), F(1.));

        let mass = F(0.);

        let energies_cache = [
            p0.temporal.value,
            p5.temporal.value,
            compute_one_loop_energy(k, zero_vector.spatial, mass),
            compute_one_loop_energy(k - q.spatial, zero_vector.spatial, mass),
            compute_one_loop_energy(k - l, zero_vector.spatial, mass),
            compute_one_loop_energy(l, zero_vector.spatial, mass),
            compute_one_loop_energy(q.spatial - l, zero_vector.spatial, mass),
            compute_one_loop_energy(l - m, zero_vector.spatial, mass),
            compute_one_loop_energy(m, zero_vector.spatial, mass),
            compute_one_loop_energy(m - q.spatial, zero_vector.spatial, mass),
        ];

        let virtual_energy_cache = energies_cache[2..].to_vec();

        let energy_prefactor = virtual_energy_cache
            .iter()
            .map(|e| (F(2.) * e).inv())
            .reduce(|acc, x| acc * x)
            .unwrap();

        let settings = tests::load_default_settings();

        let res = cff.eager_evaluate(&energies_cache, &settings) * energy_prefactor;

        let absolute_error = res - F(1.2625322619777278e-21);
        let relative_error = absolute_error / res;
        assert!(
            relative_error.abs() < F(1.0e-15),
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

        let dep_mom = 3;
        let dep_mom_expr = vec![(0, -1), (1, -1), (2, -1)];

        let orientations = generate_orientations_for_testing(edges, incoming_vertices);

        // get time before cff generation
        let start = std::time::Instant::now();

        let _cff =
            generate_cff_from_orientations(orientations, None, None, None, dep_mom, &dep_mom_expr)
                .unwrap();

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

        let dep_mom = 7;
        let dep_mom_expr = (0..7).map(|i| (i, -1)).collect();

        let incoming_vertices = vec![0, 1, 2, 3, 4, 5, 6, 7];

        let orientations = generate_orientations_for_testing(edges, incoming_vertices);

        // get time before cff generation
        let _start = std::time::Instant::now();

        let _cff =
            generate_cff_from_orientations(orientations, None, None, None, dep_mom, &dep_mom_expr)
                .unwrap();

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

        let dep_mom = 3;
        let dep_mom_expr = vec![(0, -1), (1, -1), (2, -1)];

        let start = std::time::Instant::now();
        let orientations = generate_orientations_for_testing(edges, incoming_vertices);
        let cff =
            generate_cff_from_orientations(orientations, None, None, None, dep_mom, &dep_mom_expr)
                .unwrap();
        let finish = std::time::Instant::now();
        println!("time to generate cff: {:?}", finish - start);

        let settings = load_default_settings();

        let energy_cache = [F(3.0); 17];
        let start = std::time::Instant::now();
        for _ in 0..100 {
            let _res = cff.evaluate(&energy_cache, &settings);
        }
        let finish = std::time::Instant::now();
        println!("time to evaluate cff: {:?}", (finish - start) / 100);
    }
}
