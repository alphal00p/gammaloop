use crate::{
    cff::{
        cut_expression::CFFCutExpression,
        esurface::add_external_shifts,
        expression::CompiledCFFExpression,
        hsurface::HsurfaceID,
        surface::{HybridSurface, HybridSurfaceID},
        tree::Tree,
    },
    disable,
    new_cs::{CrossSectionCut, CutId},
};
use ahash::HashMap;
use color_eyre::Report;
use color_eyre::Result;
use itertools::Itertools;
use linnet::half_edge::{
    hedgevec::HedgeVec, involution::HedgePair, subgraph::OrientedCut, HedgeGraph,
};
use linnet::half_edge::{
    involution::{EdgeIndex, Orientation},
    subgraph::InternalSubGraph,
};
use rayon::result;
use std::{
    fmt::Debug,
    iter::{Chain, Enumerate, Map},
    ops::Index,
    vec::IntoIter,
};
use symbolica::{
    atom::{Atom, AtomCore},
    id::{Pattern, Replacement},
};
use typed_index_collections::TiVec;

use serde::{Deserialize, Serialize};

use log::debug;

use super::{
    cff_graph::{CFFGenerationGraph, VertexSet},
    cut_expression::{
        CutOrientationExpression, OrientationData, OrientationExpression, OrientationID,
    },
    esurface::{Esurface, EsurfaceCollection, EsurfaceID, ExternalShift},
    expression::{CFFExpression, CFFLimit, TermId},
    hsurface::{self, Hsurface, HsurfaceCollection},
    surface::{HybridSurfaceRef, Surface, UnitSurface},
    tree::NodeId,
};

#[derive(Debug, Clone)]
struct GenerationData {
    graph: CFFGenerationGraph,
    surface_id: Option<HybridSurfaceID>,
}

#[derive(Debug, Clone)]
pub struct ShiftRewrite {
    pub dependent_momentum: EdgeIndex,
    pub dependent_momentum_expr: ExternalShift,
}

fn forget_graphs(data: GenerationData) -> HybridSurfaceID {
    data.surface_id.expect("corrupted expression tree")
}

impl GenerationData {
    fn insert_esurface(&mut self, surface_id: HybridSurfaceID) {
        self.surface_id = Some(surface_id);
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
    type Item = Orientation;
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
    type Item = Orientation;
    fn next(&mut self) -> Option<Self::Item> {
        if self.current_location < self.num_edges {
            let result_bool = self.identifier & (1 << self.current_location) == 0;
            let result = match result_bool {
                true => Orientation::Default,
                false => Orientation::Reversed,
            };

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

fn get_orientations<E, V>(graph: &HedgeGraph<E, V>) -> Vec<CFFGenerationGraph> {
    let internal_subgraph = InternalSubGraph::cleaned_filter_pessimist(graph.full_filter(), graph);
    let num_virtual_edges = graph.count_internal_edges(&internal_subgraph);
    let virtual_possible_orientations = iterate_possible_orientations(num_virtual_edges);

    virtual_possible_orientations
        .map(|orientation_of_virtuals| {
            let mut orientation_of_virtuals = orientation_of_virtuals.into_iter();

            let global_orientation = graph
                .new_hedgevec_from_iter(graph.iter_all_edges().map(|(hedge_pair, __, _)| {
                    match hedge_pair {
                        HedgePair::Unpaired { .. } => Orientation::Default,
                        HedgePair::Paired { .. } => orientation_of_virtuals
                            .next()
                            .expect(" unable to reconstruct orientation"),
                        HedgePair::Split { .. } => todo!(),
                    }
                }))
                .expect("unable to construct global orientation");

            assert!(
                orientation_of_virtuals.next().is_none(),
                "did not saturate virtual orientations when constructing global orientation"
            );

            CFFGenerationGraph::new_new(graph, global_orientation)
        })
        .collect_vec()
}

fn get_orientations_with_cut<E, V>(
    graph: &HedgeGraph<E, V>,
    oriented_cut: &OrientedCut,
) -> Vec<HedgeVec<Orientation>> {
    let internal_subgraph = InternalSubGraph::cleaned_filter_pessimist(graph.full_filter(), graph);
    let num_virtual_edges = graph.count_internal_edges(&internal_subgraph);

    let virtual_possible_orientations = iterate_possible_orientations(num_virtual_edges);

    let orientations_consistent_with_cut = virtual_possible_orientations
        .map(|orientation_of_virtuals| {
            // pad a virtual orientation with orientations of externals.
            let mut orientation_of_virtuals = orientation_of_virtuals.into_iter();

            let global_orientation = graph
                .new_hedgevec_from_iter(graph.iter_all_edges().map(|(hedge_pair, __, _)| {
                    match hedge_pair {
                        HedgePair::Unpaired { .. } => Orientation::Default,
                        HedgePair::Paired { .. } => orientation_of_virtuals
                            .next()
                            .expect(" unable to reconstruct orientation"),
                        HedgePair::Split { .. } => todo!(),
                    }
                }))
                .expect("unable to construct global orientation");

            assert!(
                orientation_of_virtuals.next().is_none(),
                "did not saturate virtual orientations when constructing global orientation"
            );

            global_orientation
        })
        .filter(|global_orientation| {
            // filter out orientations that are not consistent with the cut
            let edges_in_cut = graph.iter_edges(oriented_cut).map(|(_, id, _)| id);
            let orientation_of_edges_in_cut =
                oriented_cut.iter_edges_relative(graph).map(|(or, _)| or);

            edges_in_cut
                .zip(orientation_of_edges_in_cut)
                .all(|(edge_id, orientation)| global_orientation[edge_id] == orientation)
        })
        .filter(|global_orientation| {
            // filter out orientations that have a directed cycle
            let graph = CFFGenerationGraph::new_new(graph, global_orientation.clone());
            !graph.has_directed_cycle_initial()
        });

    orientations_consistent_with_cut.collect_vec()
}

fn get_possible_orientations_for_cut_list<E, V>(
    graph: &HedgeGraph<E, V>,
    cuts: &TiVec<CutId, CrossSectionCut>,
) -> TiVec<OrientationID, OrientationData> {
    let internal_subgraph = InternalSubGraph::cleaned_filter_pessimist(graph.full_filter(), graph);
    let num_virtual_edges = graph.count_internal_edges(&internal_subgraph);

    let virtual_possible_orientations = iterate_possible_orientations(num_virtual_edges);

    let global_orientations = virtual_possible_orientations.map(|orientation_of_virtuals| {
        // pad a virtual orientation with orientations of externals.
        let mut orientation_of_virtuals = orientation_of_virtuals.into_iter();

        let global_orientation = graph
            .new_hedgevec_from_iter(graph.iter_all_edges().map(|(hedge_pair, __, _)| {
                match hedge_pair {
                    HedgePair::Unpaired { .. } => Orientation::Default,
                    HedgePair::Paired { .. } => orientation_of_virtuals
                        .next()
                        .expect(" unable to reconstruct orientation"),
                    HedgePair::Split { .. } => todo!(),
                }
            }))
            .expect("unable to construct global orientation");

        assert!(
            orientation_of_virtuals.next().is_none(),
            "did not saturate virtual orientations when constructing global orientation"
        );

        global_orientation
    });

    // filter out orientations that are not dags
    let filter_non_dag = global_orientations.filter(|global_orientation| {
        let graph = CFFGenerationGraph::new_new(graph, global_orientation.clone());
        !graph.has_directed_cycle_initial()
    });

    // find the cuts that are consistent with the orientation
    // remove orientations with no consistent cuts
    let orientations = filter_non_dag
        .filter_map(|global_orientation| {
            let cuts_consistent_with_orientation = cuts
                .iter_enumerated()
                .filter(|(_cut_id, cut)| {
                    let edges_in_cut = graph.iter_edges(&cut.cut).map(|(_, id, _)| id);
                    let orientation_of_edges_in_cut =
                        cut.cut.iter_edges_relative(graph).map(|(or, _)| or);

                    edges_in_cut
                        .zip(orientation_of_edges_in_cut)
                        .all(|(edge_id, orientation)| global_orientation[edge_id] == orientation)
                })
                .map(|(cut_id, _)| cut_id)
                .collect_vec();

            if cuts_consistent_with_orientation.is_empty() {
                None
            } else {
                Some(OrientationData {
                    orientation: global_orientation,
                    cuts: cuts_consistent_with_orientation,
                })
            }
        })
        .collect();

    orientations
}

pub fn generate_cff_expression<E, V>(
    graph: &HedgeGraph<E, V>,
    canonize_esurface: &Option<ShiftRewrite>,
) -> Result<CFFExpression, Report> {
    let graphs = get_orientations(graph);
    debug!("number of orientations: {}", graphs.len());

    let graph_cff = generate_cff_from_orientations(graphs, None, None, None, canonize_esurface)?;

    Ok(graph_cff)
}

fn generate_cff_for_orientation<E, V>(
    graph: &HedgeGraph<E, V>,
    canonize_esurface: &Option<ShiftRewrite>,
    cache: &mut SurfaceCache,
    cuts: &TiVec<CutId, CrossSectionCut>,
    orientation_data: &OrientationData,
) -> Vec<CutOrientationExpression> {
    orientation_data
        .cuts
        .iter()
        .map(|cut_id| {
            let cut = &cuts[*cut_id];
            let left_diagram = CFFGenerationGraph::new_from_subgraph(
                graph,
                orientation_data.orientation.clone(),
                &cut.left,
            );
            let right_diagram = CFFGenerationGraph::new_from_subgraph(
                graph,
                orientation_data.orientation.clone(),
                &cut.right,
            );

            let left_tree =
                generate_tree_for_orientation(left_diagram, cache, None, canonize_esurface)
                    .map(forget_graphs);

            let right_tree =
                generate_tree_for_orientation(right_diagram, cache, None, canonize_esurface)
                    .map(forget_graphs);

            CutOrientationExpression {
                left: left_tree,
                right: right_tree,
            }
        })
        .collect()
}

pub fn generate_cff_with_cuts<E, V>(
    graph: &HedgeGraph<E, V>,
    canonize_esurface: &Option<ShiftRewrite>,
    cuts: &TiVec<CutId, CrossSectionCut>,
) -> CFFCutExpression {
    let orientations = get_possible_orientations_for_cut_list(graph, cuts);
    let mut surface_cache = SurfaceCache {
        esurface_cache: EsurfaceCollection::from_iter(std::iter::empty()),
        hsurface_cache: HsurfaceCollection::from_iter(std::iter::empty()),
    };

    let orientation_expressions = orientations
        .into_iter()
        .map(|orientation_data| {
            let cut_expressions = generate_cff_for_orientation(
                graph,
                canonize_esurface,
                &mut surface_cache,
                cuts,
                &orientation_data,
            );

            OrientationExpression {
                data: orientation_data,
                expressions: cut_expressions,
            }
        })
        .collect();

    CFFCutExpression {
        orientations: orientation_expressions,
        surfaces: surface_cache,
    }
}

pub fn generate_cff_limit(
    left_dags: Vec<CFFGenerationGraph>,
    right_dags: Vec<CFFGenerationGraph>,
    esurfaces: &EsurfaceCollection,
    limit_esurface: &Esurface,
    canonize_esurface: &Option<ShiftRewrite>,
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
        canonize_esurface,
    )
    .unwrap();
    let right = generate_cff_from_orientations(
        right_dags,
        Some(esurfaces.clone()),
        None,
        Some(limit_esurface),
        canonize_esurface,
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
    canonize_esurface: &Option<ShiftRewrite>,
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

    let mut generator_cache = SurfaceCache {
        esurface_cache,
        hsurface_cache,
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
        .map(|graph| {
            let global_orientation = graph.global_orientation.clone();
            let tree = generate_tree_for_orientation(
                graph.clone(),
                &mut generator_cache,
                rewrite_at_cache_growth,
                canonize_esurface,
            );
            let expression = tree.map(forget_graphs);

            crate::cff::expression::OrientationExpression {
                expression,
                orientation: global_orientation,
                dag: graph,
            }
        })
        .collect_vec();

    Ok(CFFExpression {
        orientations: terms.into(),
        esurfaces: generator_cache.esurface_cache,
        hsurfaces: generator_cache.hsurface_cache,
        compiled: CompiledCFFExpression::None,
    })
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct SurfaceCache {
    pub esurface_cache: EsurfaceCollection,
    pub hsurface_cache: HsurfaceCollection,
}

impl SurfaceCache {
    pub fn substitute_energies(&self, atom: &Atom) -> Atom {
        let replacement_rules = self.get_all_replacements();
        atom.replace_all_multiple(&replacement_rules)
    }

    pub fn iter_all_surfaces(
        &self,
    ) -> impl Iterator<Item = (HybridSurfaceID, HybridSurfaceRef)> + '_ {
        let esurface_id_iter = self.esurface_cache.iter_enumerated().map(|(id, esurface)| {
            (
                HybridSurfaceID::Esurface(id),
                HybridSurfaceRef::Esurface(esurface),
            )
        });

        let hsurface_id_iter = self.hsurface_cache.iter_enumerated().map(|(id, hsurface)| {
            (
                HybridSurfaceID::Hsurface(id),
                HybridSurfaceRef::Hsurface(hsurface),
            )
        });

        esurface_id_iter.chain(hsurface_id_iter)
    }

    pub fn get_all_replacements(&self) -> Vec<Replacement> {
        self.iter_all_surfaces()
            .map(|(id, surface)| {
                let id_atom = Pattern::from(Atom::from(id));
                let surface_atom = Pattern::from(Atom::from(surface));
                Replacement::new(id_atom, surface_atom)
            })
            .collect()
    }

    pub fn get_surface(&self, surface_id: HybridSurfaceID) -> HybridSurfaceRef {
        match surface_id {
            HybridSurfaceID::Esurface(id) => HybridSurfaceRef::Esurface(&self.esurface_cache[id]),
            HybridSurfaceID::Hsurface(id) => HybridSurfaceRef::Hsurface(&self.hsurface_cache[id]),
            HybridSurfaceID::Unit => HybridSurfaceRef::Unit(UnitSurface {}),
        }
    }
}

fn generate_tree_for_orientation(
    graph: CFFGenerationGraph,
    generator_cache: &mut SurfaceCache,
    rewrite_at_cache_growth: Option<&Esurface>,
    canonize_esurface: &Option<ShiftRewrite>,
) -> Tree<GenerationData> {
    let mut tree = Tree::from_root(GenerationData {
        graph,
        surface_id: None,
    });

    while let Some(()) = advance_tree(
        &mut tree,
        generator_cache,
        rewrite_at_cache_growth,
        canonize_esurface,
    ) {}

    tree
}

fn advance_tree(
    tree: &mut Tree<GenerationData>,
    generator_cache: &mut SurfaceCache,
    rewrite_at_cache_growth: Option<&Esurface>,
    canonize_esurface: &Option<ShiftRewrite>,
) -> Option<()> {
    let bottom_layer = tree.get_bottom_layer();

    let (children_optional, new_surfaces_for_tree): (
        Vec<Option<Vec<CFFGenerationGraph>>>,
        Vec<HybridSurfaceID>,
    ) = bottom_layer
        .iter()
        .map(|&node_id| {
            let node = &tree.get_node(node_id);
            let graph = &node.data.graph;

            let (option_children, surface) = graph.generate_children();

            let surface_id = match surface {
                HybridSurface::Esurface(mut esurface) => {
                    if let Some(shift_rewrite) = canonize_esurface {
                        esurface.canonicalize_shift(shift_rewrite);
                    }
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

                                if let Some(shift_rewrite) = canonize_esurface {
                                    esurface.canonicalize_shift(shift_rewrite);
                                }

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
                let child_node = GenerationData {
                    graph: child,
                    surface_id: None,
                };

                tree.insert_node(node_id, child_node);
            });
        });
    Some(())
}

#[cfg(test)]
mod tests_cff {
    use linnet::half_edge::{
        builder::HedgeGraphBuilder, involution::Flow, nodestorage::NodeStorageVec,
    };
    use symbolica::domains::float::{NumericalFloatLike, Real};
    use utils::FloatLike;

    use crate::{
        cff::cff_graph::CFFEdgeType,
        momentum::{FourMomentum, ThreeMomentum},
        tests::{self, load_default_settings},
        utils::{self, dummy_hedge_graph, RefDefault, F},
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
                    match edge_orientation {
                        Orientation::Default => {
                            new_edges[edge_id] = edges[edge_id];
                        }
                        Orientation::Reversed => {
                            let rotated_edge = (edges[edge_id].1, edges[edge_id].0);
                            new_edges[edge_id] = rotated_edge;
                        }
                        Orientation::Undirected => {
                            unreachable!("unexpected orientation")
                        }
                    }
                }

                CFFGenerationGraph::from_vec(new_edges, incoming_vertices.clone(), None)
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
        assert_eq!(
            orientation1,
            vec![
                Orientation::Default,
                Orientation::Default,
                Orientation::Default
            ]
        );

        let orientation2 = orientations[1].into_iter().collect_vec();
        assert_eq!(
            orientation2,
            vec![
                Orientation::Reversed,
                Orientation::Default,
                Orientation::Default
            ]
        );

        let orientation3 = orientations[2].into_iter().collect_vec();
        assert_eq!(
            orientation3,
            vec![
                Orientation::Default,
                Orientation::Reversed,
                Orientation::Default
            ]
        );

        let orientation4 = orientations[3].into_iter().collect_vec();
        assert_eq!(
            orientation4,
            vec![
                Orientation::Reversed,
                Orientation::Reversed,
                Orientation::Default
            ]
        );

        let orientation5 = orientations[4].into_iter().collect_vec();
        assert_eq!(
            orientation5,
            vec![
                Orientation::Default,
                Orientation::Default,
                Orientation::Reversed
            ]
        );

        let orientation6 = orientations[5].into_iter().collect_vec();
        assert_eq!(
            orientation6,
            vec![
                Orientation::Reversed,
                Orientation::Default,
                Orientation::Reversed
            ]
        );

        let orientation7 = orientations[6].into_iter().collect_vec();
        assert_eq!(
            orientation7,
            vec![
                Orientation::Default,
                Orientation::Reversed,
                Orientation::Reversed
            ]
        );

        let orientation8 = orientations[7].into_iter().collect_vec();
        assert_eq!(
            orientation8,
            vec![
                Orientation::Reversed,
                Orientation::Reversed,
                Orientation::Reversed
            ]
        );
    }

    #[test]
    fn test_cff_generation_triangle() {
        let triangle = vec![(2, 0), (0, 1), (1, 2)];

        let incoming_vertices = vec![0, 1, 2];
        let orientations = generate_orientations_for_testing(triangle, incoming_vertices);
        assert_eq!(orientations.len(), 6);

        let dep_mom = EdgeIndex::from(2);
        let dep_mom_expr = vec![(EdgeIndex::from(0), -1), (EdgeIndex::from(1), -1)];

        let shift_rewrite = Some(ShiftRewrite {
            dependent_momentum: dep_mom,
            dependent_momentum_expr: dep_mom_expr,
        });

        let cff =
            generate_cff_from_orientations(orientations, None, None, None, &shift_rewrite.clone())
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

        let energy_cache = dummy_hedge_graph(6)
            .new_hedgevec_from_iter(energy_cache)
            .unwrap();

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

        // test cff from hedge graph
        let mut triangle_hedge_graph_builder = HedgeGraphBuilder::new();

        let nodes = (0..3)
            .map(|_| triangle_hedge_graph_builder.add_node(()))
            .collect_vec();

        for node in nodes.clone() {
            triangle_hedge_graph_builder.add_external_edge(
                node,
                (),
                Orientation::Undirected,
                Flow::Sink,
            );
        }

        triangle_hedge_graph_builder.add_edge(nodes[2], nodes[0], (), Orientation::Undirected);
        triangle_hedge_graph_builder.add_edge(nodes[0], nodes[1], (), Orientation::Undirected);
        triangle_hedge_graph_builder.add_edge(nodes[1], nodes[2], (), Orientation::Undirected);

        let triangle_hedge_graph = triangle_hedge_graph_builder.build::<NodeStorageVec<()>>();

        let cff_hedge = generate_cff_expression(&triangle_hedge_graph, &shift_rewrite).unwrap();

        let cff_res: F<f64> = energy_prefactor
            * cff_hedge.eager_evaluate(&energy_cache, &settings)
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
    }

    #[test]
    fn test_cff_test_double_triangle() {
        let double_triangle_edges = vec![(0, 1), (0, 2), (1, 2), (1, 3), (2, 3)];
        let incoming_vertices = vec![0, 3];

        let orientations =
            generate_orientations_for_testing(double_triangle_edges, incoming_vertices);

        let dep_mom = EdgeIndex::from(1);
        let dep_mom_expr = vec![(EdgeIndex::from(0), -1)];

        let shift_rewrite = Some(ShiftRewrite {
            dependent_momentum: dep_mom,
            dependent_momentum_expr: dep_mom_expr,
        });

        let cff =
            generate_cff_from_orientations(orientations, None, None, None, &shift_rewrite).unwrap();

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

        let energy_cache = dummy_hedge_graph(energy_cache.len())
            .new_hedgevec_from_iter(energy_cache)
            .unwrap();

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

        let mut hedge_double_triangle_builder = HedgeGraphBuilder::new();
        let nodes = (0..4)
            .map(|_| hedge_double_triangle_builder.add_node(()))
            .collect_vec();

        hedge_double_triangle_builder.add_external_edge(
            nodes[0],
            (),
            Orientation::Undirected,
            Flow::Sink,
        );
        hedge_double_triangle_builder.add_external_edge(
            nodes[3],
            (),
            Orientation::Undirected,
            Flow::Sink,
        );

        hedge_double_triangle_builder.add_edge(nodes[0], nodes[1], (), Orientation::Undirected);
        hedge_double_triangle_builder.add_edge(nodes[0], nodes[2], (), Orientation::Undirected);
        hedge_double_triangle_builder.add_edge(nodes[1], nodes[2], (), Orientation::Undirected);
        hedge_double_triangle_builder.add_edge(nodes[1], nodes[3], (), Orientation::Undirected);
        hedge_double_triangle_builder.add_edge(nodes[2], nodes[3], (), Orientation::Undirected);

        let hedge_double_traingle = hedge_double_triangle_builder.build::<NodeStorageVec<()>>();
        let cff_hedge = generate_cff_expression(&hedge_double_traingle, &shift_rewrite).unwrap();

        let cff_res = energy_prefactor * cff_hedge.eager_evaluate(&energy_cache, &settings);

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

        let cuts = hedge_double_traingle.all_cuts(
            hedge_double_traingle[&nodes[3]].clone(),
            hedge_double_traingle[&nodes[0]].clone(),
        );
        let mut num_with_6_ors = 0;
        let mut num_with_4_ors = 0;
        assert_eq!(cuts.len(), 4);
        for (_, cut, _) in &cuts {
            let orientations = get_orientations_with_cut(&hedge_double_traingle, cut);
            if orientations.len() == 4 {
                num_with_4_ors += 1
            } else if orientations.len() == 6 {
                num_with_6_ors += 1
            }
        }

        assert_eq!(num_with_4_ors, 2);
        assert_eq!(num_with_6_ors, 2);
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

        let dep_mom = EdgeIndex::from(1);
        let dep_mom_expr = vec![(EdgeIndex::from(0), -1)];

        let shift_rewrite = Some(ShiftRewrite {
            dependent_momentum: dep_mom,
            dependent_momentum_expr: dep_mom_expr,
        });

        let orientataions = generate_orientations_for_testing(tbt_edges, incoming_vertices);
        let cff = generate_cff_from_orientations(orientataions, None, None, None, &shift_rewrite)
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

        let energies_cache = dummy_hedge_graph(energies_cache.len())
            .new_hedgevec_from_iter(energies_cache)
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

        let mut tbt_hedge_builder = HedgeGraphBuilder::new();
        let nodes = (0..6).map(|_| tbt_hedge_builder.add_node(())).collect_vec();
        tbt_hedge_builder.add_external_edge(nodes[0], (), Orientation::Undirected, Flow::Sink);
        tbt_hedge_builder.add_external_edge(nodes[5], (), Orientation::Undirected, Flow::Sink);

        tbt_hedge_builder.add_edge(nodes[0], nodes[1], (), Orientation::Undirected);
        tbt_hedge_builder.add_edge(nodes[2], nodes[0], (), Orientation::Undirected);
        tbt_hedge_builder.add_edge(nodes[1], nodes[2], (), Orientation::Undirected);
        tbt_hedge_builder.add_edge(nodes[1], nodes[3], (), Orientation::Undirected);
        tbt_hedge_builder.add_edge(nodes[2], nodes[4], (), Orientation::Undirected);
        tbt_hedge_builder.add_edge(nodes[3], nodes[4], (), Orientation::Undirected);
        tbt_hedge_builder.add_edge(nodes[3], nodes[5], (), Orientation::Undirected);
        tbt_hedge_builder.add_edge(nodes[5], nodes[4], (), Orientation::Undirected);

        let tbt_hedge = tbt_hedge_builder.build::<NodeStorageVec<()>>();
        let cff_hedge = generate_cff_expression(&tbt_hedge, &shift_rewrite).unwrap();
        let res = cff_hedge.eager_evaluate(&energies_cache, &settings) * energy_prefactor;

        let absolute_error = res - F(1.2625322619777278e-21);
        let relative_error = absolute_error / res;
        assert!(
            relative_error.abs() < F(1.0e-15),
            "relative error: {:+e}",
            relative_error
        );

        let cuts = tbt_hedge.all_cuts(tbt_hedge[&nodes[0]].clone(), tbt_hedge[&nodes[5]].clone());
        assert_eq!(cuts.len(), 9);
        let mut num_with_24 = 0;
        let mut num_with_16 = 0;
        let mut num_with_42 = 0;
        let mut num_with_36 = 0;
        for (_, cut, _) in cuts.iter() {
            let orientations = get_orientations_with_cut(&tbt_hedge, cut);
            if orientations.len() == 24 {
                num_with_24 += 1;
            }
            if orientations.len() == 16 {
                num_with_16 += 1;
            }
            if orientations.len() == 42 {
                num_with_42 += 1
            }
            if orientations.len() == 36 {
                num_with_36 += 1;
            }
        }

        assert_eq!(num_with_24, 4);
        assert_eq!(num_with_16, 2);
        assert_eq!(num_with_42, 2);
        assert_eq!(num_with_36, 1);
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

        let dep_mom = EdgeIndex::from(3);
        let dep_mom_expr = vec![
            (EdgeIndex::from(0), -1),
            (EdgeIndex::from(1), -1),
            (EdgeIndex::from(2), -1),
        ];

        let shift_rewrite = ShiftRewrite {
            dependent_momentum: dep_mom,
            dependent_momentum_expr: dep_mom_expr,
        };

        let orientations = generate_orientations_for_testing(edges, incoming_vertices);

        // get time before cff generation
        let start = std::time::Instant::now();

        let _cff =
            generate_cff_from_orientations(orientations, None, None, None, &Some(shift_rewrite))
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

        let dep_mom = EdgeIndex::from(7);
        let dep_mom_expr = (0..7).map(|i| (EdgeIndex::from(i), -1)).collect();

        let shift_rewrite = ShiftRewrite {
            dependent_momentum: dep_mom,
            dependent_momentum_expr: dep_mom_expr,
        };

        let incoming_vertices = vec![0, 1, 2, 3, 4, 5, 6, 7];

        let orientations = generate_orientations_for_testing(edges, incoming_vertices);

        // get time before cff generation
        let _start = std::time::Instant::now();

        let _cff =
            generate_cff_from_orientations(orientations, None, None, None, &Some(shift_rewrite))
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

        let dep_mom = EdgeIndex::from(3);
        let dep_mom_expr = vec![
            (EdgeIndex::from(0), -1),
            (EdgeIndex::from(1), -1),
            (EdgeIndex::from(2), -1),
        ];

        let shift_rewrite = ShiftRewrite {
            dependent_momentum: dep_mom,
            dependent_momentum_expr: dep_mom_expr,
        };

        let start = std::time::Instant::now();
        let orientations = generate_orientations_for_testing(edges, incoming_vertices);
        let cff =
            generate_cff_from_orientations(orientations, None, None, None, &Some(shift_rewrite))
                .unwrap();
        let finish = std::time::Instant::now();
        println!("time to generate cff: {:?}", finish - start);

        let settings = load_default_settings();

        let energy_cache = [F(3.0); 17];

        let energy_cache = dummy_hedge_graph(energy_cache.len())
            .new_hedgevec_from_iter(energy_cache)
            .unwrap();
        let start = std::time::Instant::now();
        for _ in 0..100 {
            let _res = cff.evaluate(&energy_cache, &settings);
        }
        let finish = std::time::Instant::now();
        println!("time to evaluate cff: {:?}", (finish - start) / 100);
    }
}
