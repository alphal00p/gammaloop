use std::collections::HashMap;

use crate::{
    cff::{
        hsurface::{Hsurface, HsurfaceID},
        surface::{HybridSurface, HybridSurfaceID, InfiniteSurface},
        tree::Tree,
    },
    graph::{Graph, LoopMomentumBasis, get_cff_inverse_energy_product_impl},
    processes::{CrossSectionCut, CutId},
    settings::global::OrientationPattern,
};
use ahash::HashSet;
use bincode::{Decode, Encode};
use color_eyre::Result;
use itertools::Itertools;
use linnet::half_edge::{
    HedgeGraph,
    involution::{EdgeIndex, EdgeVec, Orientation},
    subgraph::{SubGraphLike, SubSetOps},
};
use symbolica::{
    atom::{Atom, AtomCore},
    id::{Pattern, Replacement},
};
use typed_index_collections::TiVec;

use serde::{Deserialize, Serialize};

use super::{
    cff_graph::CFFGenerationGraph,
    esurface::{Esurface, EsurfaceCollection, EsurfaceID, ExternalShift},
    expression::{CFFExpression, OrientationID},
    hsurface::HsurfaceCollection,
    surface::{HybridSurfaceRef, UnitSurface},
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

#[cfg(test)]
fn generate_cff_expression<E, V, H>(
    graph: &HedgeGraph<E, V, H>,
    canonize_esurface: &Option<ShiftRewrite>,
    edges_in_initial_state_cut: &[EdgeIndex],
) -> Result<CFFExpression<OrientationID>> {
    let mut surface_cache = SurfaceCache::new();
    generate_cff_expression_from_subgraph(
        graph,
        &graph.full_graph(),
        canonize_esurface,
        &[],
        edges_in_initial_state_cut,
        &mut surface_cache,
    )
}

impl Graph {
    pub(crate) fn generate_cff(
        &mut self,
        contract_edges: &[EdgeIndex],
        canonize_esurface: &Option<ShiftRewrite>,
        orientation_pattern: &OrientationPattern,
    ) -> Result<CFFExpression<OrientationID>> {
        self.generate_cff_with_feynkit(contract_edges, canonize_esurface, orientation_pattern)
    }
}

pub fn generate_cff_expression_from_subgraph<E, V, H, S: SubGraphLike>(
    graph: &HedgeGraph<E, V, H>,
    subgraph: &S,
    canonize_esurface: &Option<ShiftRewrite>,
    reversed_dangling: &[EdgeIndex],
    edges_in_initial_state_cut: &[EdgeIndex],
    surface_cache: &mut SurfaceCache,
) -> Result<CFFExpression<OrientationID>> {
    super::feynkit::generate_cff_expression_from_subgraph_with_feynkit(
        graph,
        subgraph,
        canonize_esurface,
        reversed_dangling,
        edges_in_initial_state_cut,
        surface_cache,
    )
}

#[derive(Copy, Clone, Debug)]
pub struct ConstraintData<'a> {
    pub constraints: &'a [&'a Esurface],
    pub illegal_esurfaces: &'a [&'a Esurface],
}

#[derive(Copy, Clone, Debug)]
pub struct UvCffTopology<'a> {
    pub contract_edges: &'a [EdgeIndex],
    pub edges_in_initial_state_cut: &'a [EdgeIndex],
    pub orientation: &'a EdgeVec<Orientation>,
    pub cut_edges: &'a [EdgeIndex],
}

pub fn generate_uv_cff<E, V, H, S: SubGraphLike>(
    graph: &HedgeGraph<E, V, H>,
    subgraph: &S,
    canonize_esurface: &Option<ShiftRewrite>,
    topology: UvCffTopology<'_>,
    setup: PostProcessingSetup<'_>,
) -> Result<Atom> {
    let mut generation_graph =
        CFFGenerationGraph::new_from_subgraph(graph, topology.orientation.clone(), subgraph)?;

    for contracted_edge in topology.contract_edges {
        generation_graph = generation_graph.contract_edge(*contracted_edge);
    }

    generation_graph.remove_self_edges();

    if generation_graph.has_directed_cycle_initial() {
        return Ok(Atom::new());
    }

    let mut surface_cache = SurfaceCache {
        esurface_cache: EsurfaceCollection::from_iter(std::iter::empty()),
        hsurface_cache: HsurfaceCollection::from_iter(std::iter::empty()),
    };

    let generate_tree_for_orientation = generate_tree_for_orientation(
        generation_graph,
        &mut surface_cache,
        topology.edges_in_initial_state_cut,
        canonize_esurface,
    );

    let mut tree = generate_tree_for_orientation.map(forget_graphs);

    post_process(
        &mut tree,
        topology.orientation,
        subgraph,
        &surface_cache,
        setup,
    );

    let surface_cache_to_use = setup
        .rewrite_esurfaces
        .map_or(&surface_cache, |rewrite| rewrite.allowed_targets);

    let atom_tree = tree.to_atom_inv();
    let atom_tree_substituted =
        surface_cache_to_use.substitute_energies(&atom_tree, topology.cut_edges);
    let inverse_energies =
        get_cff_inverse_energy_product_impl(graph, subgraph, topology.contract_edges);

    Ok(atom_tree_substituted * &inverse_energies)
}

#[derive(Clone, Copy)]
pub struct PostProcessingSetup<'a> {
    pub constraint_data: Option<ConstraintData<'a>>,
    pub rewrite_esurfaces: Option<EsurfaceRewritingInstructions<'a>>,
}

#[derive(Clone, Copy)]
pub struct EsurfaceRewritingInstructions<'a> {
    pub allowed_targets: &'a SurfaceCache,
    pub graph: &'a Graph,
    pub cuts: &'a TiVec<CutId, CrossSectionCut>,
    pub subgraph_location: (Option<CutId>, Option<CutId>),
}

fn post_process<S: SubGraphLike>(
    tree: &mut Tree<HybridSurfaceID>,
    orientation: &EdgeVec<Orientation>,
    subgraph: &S,
    surface_cache: &SurfaceCache,
    setup: PostProcessingSetup<'_>,
) {
    if let Some(constraint_data) = setup.constraint_data {
        tree.map_mut(|surface_id| {
            let esurface_is_allowed = match surface_id {
                HybridSurfaceID::Esurface(esurface_id) => {
                    let esurface_to_compare = &surface_cache.esurface_cache[*esurface_id];
                    constraint_data
                        .illegal_esurfaces
                        .iter()
                        .all(|illegal_esurface| esurface_to_compare != *illegal_esurface)
                }
                HybridSurfaceID::Hsurface(hsurface_id) => {
                    let hsurface_to_compare = &surface_cache.hsurface_cache[*hsurface_id];
                    constraint_data
                        .illegal_esurfaces
                        .iter()
                        .all(|illegal_esurface| {
                            !hsurface_to_compare
                                .equality_under_energy_conservation(
                                    illegal_esurface,
                                    constraint_data.constraints,
                                )
                                .unwrap_or(
                                    hsurface_to_compare.equality_by_try_convert(illegal_esurface),
                                )
                        })
                }
                HybridSurfaceID::Unit => true,
                HybridSurfaceID::Infinite => true,
            };

            if !esurface_is_allowed {
                *surface_id = HybridSurfaceID::Infinite
            }
        });
    }

    if let Some(rewrite_esurfaces) = setup.rewrite_esurfaces {
        let hashset_of_appearing_ids = tree
            .iter_nodes()
            .map(|node| node.data)
            .collect::<HashSet<HybridSurfaceID>>();

        let mut id_map = HashMap::<HybridSurfaceID, HybridSurfaceID>::new();
        id_map.insert(HybridSurfaceID::Unit, HybridSurfaceID::Unit);
        id_map.insert(HybridSurfaceID::Infinite, HybridSurfaceID::Infinite);

        for appearing_id in hashset_of_appearing_ids.iter() {
            let surface_to_rewrite = surface_cache.get_surface(*appearing_id);

            match surface_to_rewrite {
                HybridSurfaceRef::Unit(_) => continue,
                HybridSurfaceRef::Infinite(_) => continue,
                HybridSurfaceRef::Esurface(esurface) => {
                    if let Some(esurface_id) = rewrite_esurfaces
                        .allowed_targets
                        .esurface_cache
                        .position(|allowed_esurface| allowed_esurface == esurface)
                    {
                        let new_id = HybridSurfaceID::Esurface(esurface_id);
                        id_map.insert(*appearing_id, new_id);
                    } else {
                        let complete_to_right =
                            if let Some(cut_id) = rewrite_esurfaces.subgraph_location.1 {
                                let edges_in_cut = rewrite_esurfaces
                                    .graph
                                    .iter_edges_of(&rewrite_esurfaces.cuts[cut_id].cut)
                                    .map(|(_, edge_id, _)| edge_id)
                                    .collect_vec();

                                edges_in_cut
                                    .iter()
                                    .all(|edge_id| esurface.energies.contains(edge_id))
                            } else {
                                false
                            };

                        let complete_to_left =
                            if let Some(cut_id) = rewrite_esurfaces.subgraph_location.0 {
                                let edges_in_cut = rewrite_esurfaces
                                    .graph
                                    .iter_edges_of(&rewrite_esurfaces.cuts[cut_id].cut)
                                    .map(|(_, edge_id, _)| edge_id)
                                    .collect_vec();

                                edges_in_cut
                                    .iter()
                                    .all(|edge_id| esurface.energies.contains(edge_id))
                            } else {
                                false
                            };

                        if complete_to_left && complete_to_right {
                            panic!("esurface has no connected component");
                        }

                        if !complete_to_left && !complete_to_right {
                            println!("esurface: {:#?}", esurface);
                            panic!("esurface cannot be rewritten to any allowed target");
                        }

                        let vertices_to_add = if complete_to_left {
                            let cut_id = rewrite_esurfaces.subgraph_location.0.unwrap();
                            &rewrite_esurfaces.cuts[cut_id].left
                        } else if complete_to_right {
                            let cut_id = rewrite_esurfaces.subgraph_location.1.unwrap();
                            &rewrite_esurfaces.cuts[cut_id].right
                        } else {
                            unreachable!()
                        };

                        let new_esurface_subgraph = esurface
                            .vertex_set
                            .subgraph(rewrite_esurfaces.graph)
                            .union(vertices_to_add);

                        let new_esurface = Esurface::new_from_subgraph(
                            &new_esurface_subgraph,
                            rewrite_esurfaces.graph,
                            orientation,
                        );

                        let new_esurface_id = rewrite_esurfaces
                            .allowed_targets
                            .esurface_cache
                            .position(|allowed_esurface| allowed_esurface == &new_esurface)
                            .expect("constructed esurface not in allowed targets");

                        let new_id = HybridSurfaceID::Esurface(new_esurface_id);
                        id_map.insert(*appearing_id, new_id);
                    }
                }
                HybridSurfaceRef::Hsurface(hsurface) => {
                    let complete_to_left =
                        if let Some(cut_id) = rewrite_esurfaces.subgraph_location.0 {
                            let edges_in_cut = rewrite_esurfaces
                                .graph
                                .iter_edges_of(&rewrite_esurfaces.cuts[cut_id].cut)
                                .map(|(_, edge_id, _)| edge_id)
                                .collect_vec();

                            hsurface
                                .negative_energies
                                .iter()
                                .all(|edge_id| edges_in_cut.contains(edge_id))
                        } else {
                            false
                        };

                    let complete_to_right =
                        if let Some(cut_id) = rewrite_esurfaces.subgraph_location.1 {
                            let edges_in_cut = rewrite_esurfaces
                                .graph
                                .iter_edges_of(&rewrite_esurfaces.cuts[cut_id].cut)
                                .map(|(_, edge_id, _)| edge_id)
                                .collect_vec();

                            hsurface
                                .negative_energies
                                .iter()
                                .all(|edge_id| edges_in_cut.contains(edge_id))
                        } else {
                            false
                        };

                    if complete_to_left && complete_to_right {
                        panic!(
                            "hsurface has no connected component supergraph, it cannot exist, but it does"
                        );
                    }

                    if !complete_to_left && !complete_to_right {
                        println!("hsurface: {:#?}", hsurface);
                        panic!("hsurface cannot be rewritten to any allowed target");
                    }

                    let vertices_to_add = if complete_to_left {
                        let cut_id = rewrite_esurfaces.subgraph_location.0.unwrap();
                        &rewrite_esurfaces.cuts[cut_id].left
                    } else if complete_to_right {
                        let cut_id = rewrite_esurfaces.subgraph_location.1.unwrap();
                        &rewrite_esurfaces.cuts[cut_id].right
                    } else {
                        unreachable!()
                    };

                    let new_esurface_subgraph = hsurface
                        .vertex_set
                        .subgraph(rewrite_esurfaces.graph)
                        .union(vertices_to_add);

                    let new_esurface = Esurface::new_from_subgraph(
                        &new_esurface_subgraph,
                        rewrite_esurfaces.graph,
                        orientation,
                    );

                    let new_esurface_id = rewrite_esurfaces
                        .allowed_targets
                        .esurface_cache
                        .position(|allowed_esurface| allowed_esurface == &new_esurface)
                        .unwrap_or_else(|| {
                            println!("for graph: {}", rewrite_esurfaces.graph.name.clone());
                            println!("dot: \n {}", rewrite_esurfaces.graph.debug_dot());
                            println!("subgraph: \n {}", rewrite_esurfaces.graph.dot(subgraph));

                            println!("from hsurface: {:?}", hsurface);
                            println!("constructed esurface: {:?}", new_esurface);
                            panic!("constructed esurface not in allowed targets");
                        });

                    let new_id = HybridSurfaceID::Esurface(new_esurface_id);
                    id_map.insert(*appearing_id, new_id);
                }
            }
        }

        tree.map_mut(|surface_id| *surface_id = id_map[surface_id]);
    }
}

#[derive(Clone, Debug, Serialize, Deserialize, Encode, Decode)]
pub struct SurfaceCache {
    #[bincode(with_serde)]
    pub esurface_cache: EsurfaceCollection, // Esurfaces of the supergraph
    #[bincode(with_serde)]
    pub hsurface_cache: HsurfaceCollection, // Anything else.
}

impl SurfaceCache {
    pub fn substitute_energies(&self, atom: &Atom, cut_edges: &[EdgeIndex]) -> Atom {
        let replacement_rules = self.get_all_replacements(cut_edges);
        atom.replace_multiple(&replacement_rules)
    }

    pub(crate) fn iter_all_surfaces(
        &'_ self,
    ) -> impl Iterator<Item = (HybridSurfaceID, HybridSurfaceRef<'_>)> + '_ {
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

    pub(crate) fn get_all_replacements(&self, cut_edges: &[EdgeIndex]) -> Vec<Replacement> {
        self.iter_all_surfaces()
            .map(|(id, surface)| {
                let id_atom = Pattern::from(Atom::from(id));
                let surface_atom = Pattern::from(surface.to_atom(cut_edges));
                Replacement::new(id_atom, surface_atom)
            })
            .collect()
    }

    pub(crate) fn get_all_replacements_in_lmb(
        &self,
        cut_edges: &[EdgeIndex],
        lmb: &LoopMomentumBasis,
    ) -> Vec<Replacement> {
        self.iter_all_surfaces()
            .map(|(id, surface)| {
                let id_atom = Pattern::from(Atom::from(id));
                let surface_atom = Pattern::from(surface.to_atom_in_lmb(cut_edges, lmb));
                Replacement::new(id_atom, surface_atom)
            })
            .collect()
    }

    #[allow(dead_code)]
    pub(crate) fn get_surface(&self, surface_id: HybridSurfaceID) -> HybridSurfaceRef<'_> {
        match surface_id {
            HybridSurfaceID::Esurface(id) => HybridSurfaceRef::Esurface(&self.esurface_cache[id]),
            HybridSurfaceID::Hsurface(id) => HybridSurfaceRef::Hsurface(&self.hsurface_cache[id]),
            HybridSurfaceID::Unit => HybridSurfaceRef::Unit(UnitSurface {}),
            HybridSurfaceID::Infinite => HybridSurfaceRef::Infinite(InfiniteSurface {}),
        }
    }

    #[allow(dead_code)]
    pub(crate) fn new() -> Self {
        Self {
            esurface_cache: EsurfaceCollection::from_iter(std::iter::empty()),
            hsurface_cache: HsurfaceCollection::from_iter(std::iter::empty()),
        }
    }
}

fn generate_tree_for_orientation(
    graph: CFFGenerationGraph,
    generator_cache: &mut SurfaceCache,
    edges_in_initial_state_cut: &[EdgeIndex],
    canonize_esurface: &Option<ShiftRewrite>,
) -> Tree<GenerationData> {
    let mut tree = Tree::from_root(GenerationData {
        graph,
        surface_id: None,
    });

    while let Some(()) = advance_tree(
        &mut tree,
        generator_cache,
        edges_in_initial_state_cut,
        canonize_esurface,
    ) {}

    tree
}

fn advance_tree(
    tree: &mut Tree<GenerationData>,
    generator_cache: &mut SurfaceCache,
    edges_in_initial_state_cut: &[EdgeIndex],
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

            // treat the edges in the initial state cut as true externals
            let surface = match surface {
                HybridSurface::Esurface(esurface) => {
                    let energies_to_be_moved = esurface
                        .energies
                        .iter()
                        .filter(|edge_id| edges_in_initial_state_cut.contains(edge_id))
                        .copied()
                        .collect_vec();

                    if energies_to_be_moved.is_empty() {
                        HybridSurface::Esurface(esurface)
                    } else {
                        let new_energies = esurface
                            .energies
                            .iter()
                            .filter(|edge_id| !energies_to_be_moved.contains(edge_id))
                            .copied()
                            .collect_vec();

                        let mut new_shift = esurface.external_shift.clone();
                        for energy_to_move in energies_to_be_moved.iter() {
                            new_shift.push((*energy_to_move, 1));
                        }

                        new_shift.sort_by_key(|(edge_id, _)| *edge_id);

                        HybridSurface::Esurface(Esurface {
                            energies: new_energies,
                            external_shift: new_shift,
                            vertex_set: esurface.vertex_set,
                        })
                    }
                }
                HybridSurface::Unit(unit) => HybridSurface::Unit(unit),
                HybridSurface::Infinite(infinite) => HybridSurface::Infinite(infinite),
                HybridSurface::Hsurface(hsurface) => {
                    let positive_energies_to_be_moved = hsurface
                        .positive_energies
                        .iter()
                        .filter(|edge_id| edges_in_initial_state_cut.contains(edge_id))
                        .copied()
                        .collect_vec();

                    let negative_energies_to_be_moved = hsurface
                        .negative_energies
                        .iter()
                        .filter(|edge_id| edges_in_initial_state_cut.contains(edge_id))
                        .copied()
                        .collect_vec();

                    if positive_energies_to_be_moved.is_empty()
                        && negative_energies_to_be_moved.is_empty()
                    {
                        HybridSurface::Hsurface(hsurface)
                    } else if !positive_energies_to_be_moved.is_empty()
                        && negative_energies_to_be_moved.is_empty()
                    {
                        let new_positive_energies = hsurface
                            .positive_energies
                            .iter()
                            .filter(|edge_id| !positive_energies_to_be_moved.contains(edge_id))
                            .copied()
                            .collect_vec();

                        let mut new_shift = hsurface.external_shift.clone();

                        for positive_energy_to_move in positive_energies_to_be_moved.iter() {
                            new_shift.push((*positive_energy_to_move, 1));
                        }

                        new_shift.sort_by_key(|(edge_id, _)| *edge_id);

                        HybridSurface::Hsurface(Hsurface {
                            positive_energies: new_positive_energies,
                            negative_energies: hsurface.negative_energies.clone(),
                            external_shift: new_shift,
                            vertex_set: hsurface.vertex_set,
                        })
                    } else if !negative_energies_to_be_moved.is_empty()
                        && positive_energies_to_be_moved.is_empty()
                    {
                        let new_negative_energies = hsurface
                            .negative_energies
                            .iter()
                            .filter(|edge_id| !negative_energies_to_be_moved.contains(edge_id))
                            .copied()
                            .collect_vec();

                        let mut new_shift = hsurface.external_shift.clone();

                        for negative_energy_to_move in negative_energies_to_be_moved.iter() {
                            new_shift.push((*negative_energy_to_move, -1));
                        }

                        new_shift.sort_by_key(|(edge_id, _)| *edge_id);

                        if new_negative_energies.is_empty() {
                            HybridSurface::Esurface(Esurface {
                                energies: hsurface.positive_energies.clone(),
                                external_shift: new_shift,
                                vertex_set: hsurface.vertex_set,
                            })
                        } else {
                            HybridSurface::Hsurface(Hsurface {
                                positive_energies: hsurface.positive_energies.clone(),
                                negative_energies: new_negative_energies,
                                external_shift: new_shift,
                                vertex_set: hsurface.vertex_set,
                            })
                        }
                    } else {
                        unreachable!()
                    }
                }
            };

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
                HybridSurface::Unit(_) => HybridSurfaceID::Unit,
                HybridSurface::Infinite(_) => HybridSurfaceID::Infinite,
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
    use std::ops::Range;

    use linnet::half_edge::{
        builder::HedgeGraphBuilder, involution::Flow, nodestore::NodeStorageVec,
    };
    use symbolica::{
        evaluate::{ExpressionEvaluator, FunctionMap, OptimizationSettings},
        parse, symbol,
    };
    use utils::FloatLike;

    use crate::{
        momentum::{FourMomentum, ThreeMomentum},
        settings::global::OrientationPattern,
        utils::{
            self, F, RefDefault, external_energy_atom_from_index, ose_atom_from_index,
            test_utils::dummy_hedge_graph,
        },
    };

    use super::*;

    // helper function to make a symbolica evaluator
    impl CFFExpression<OrientationID> {
        fn quick_symbolica_evaluator(
            &self,
            external_range: Range<usize>,
            virtual_range: Range<usize>,
        ) -> ExpressionEvaluator<F<f64>> {
            let expression_atom_no_energy_sub = self.to_atom(OrientationPattern::default());
            let num_energies = external_range.end.max(virtual_range.end);
            let mut expression_atom = self
                .surfaces
                .substitute_energies(&expression_atom_no_energy_sub, &[]);
            for edge_id in 0..num_energies {
                let edge_id = EdgeIndex::from(edge_id);
                expression_atom = expression_atom
                    .replace(ose_atom_from_index(edge_id))
                    .with(external_energy_atom_from_index(edge_id));
            }

            let params = (0..num_energies)
                .map(|i| external_energy_atom_from_index(EdgeIndex::from(i)))
                .collect_vec();

            let function_map = FunctionMap::new();

            let mut tree = expression_atom
                .as_view()
                .to_evaluation_tree(&function_map, &params)
                .unwrap();

            tree.horner_scheme();
            tree.common_subexpression_elimination();
            tree.linearize(&OptimizationSettings::default())
                .map_coeff(&|c| (&c.re).into())
        }
    }

    fn compute_one_loop_energy<T: FloatLike>(
        k: ThreeMomentum<F<T>>,
        p: ThreeMomentum<F<T>>,
        m: F<T>,
    ) -> F<T> {
        ((k + p).norm_squared() + &m * &m).sqrt()
    }

    fn proper_atom(graph: &HedgeGraph<(), ()>) -> Atom {
        let cff = generate_cff_expression(graph, &None, &[]).unwrap();

        let mut cff_atom = cff.to_atom(OrientationPattern::default());
        cff_atom = cff.surfaces.substitute_energies(&cff_atom, &[]);
        let inverse_energy_product =
            get_cff_inverse_energy_product_impl(graph, &graph.full_graph(), &[]);

        cff_atom *= inverse_energy_product;
        cff_atom
    }

    #[test]
    fn test_dot_trick_bubble() {
        let mut dotted_topology_builder = HedgeGraphBuilder::new();
        let dotted_nodes = (0..3)
            .map(|_| dotted_topology_builder.add_node(()))
            .collect_vec();

        dotted_topology_builder.add_edge(dotted_nodes[0], dotted_nodes[1], (), false);
        dotted_topology_builder.add_edge(dotted_nodes[1], dotted_nodes[2], (), false);
        dotted_topology_builder.add_edge(dotted_nodes[2], dotted_nodes[0], (), false);
        let dotted_topology = dotted_topology_builder.build();

        let mut dotted_cff_atom = proper_atom(&dotted_topology);
        dotted_cff_atom = dotted_cff_atom
            .replace(parse!("OSE(2)"))
            .with(parse!("OSE(1)"));

        let mut topology_builder = HedgeGraphBuilder::new();
        let nodes = (0..2).map(|_| topology_builder.add_node(())).collect_vec();

        topology_builder.add_edge(nodes[0], nodes[1], (), false);
        topology_builder.add_edge(nodes[0], nodes[1], (), false);
        let toplogy = topology_builder.build();
        let mut cff_atom = proper_atom(&toplogy);

        cff_atom = cff_atom.replace(parse!("OSE(1)")).with(parse!("OSE1"));
        cff_atom = cff_atom.derivative(symbol!("OSE1"));
        cff_atom = cff_atom.replace(parse!("OSE1")).with(parse!("OSE(1)")) / parse!("2*OSE(1)");

        let diff = (&cff_atom - &dotted_cff_atom).expand();

        println!("cff_atom: {}", cff_atom.expand());
        println!("dotted_cff_atom: {}", dotted_cff_atom.expand());

        println!("diff: {}", diff);
    }

    #[test]
    fn test_dot_trick_amg() {
        let mut dotted_topology_builder = HedgeGraphBuilder::new();
        let dotted_nodes = (0..5)
            .map(|_| dotted_topology_builder.add_node(()))
            .collect_vec();

        dotted_topology_builder.add_edge(dotted_nodes[0], dotted_nodes[3], (), false);
        dotted_topology_builder.add_edge(dotted_nodes[0], dotted_nodes[2], (), false);
        dotted_topology_builder.add_edge(dotted_nodes[0], dotted_nodes[1], (), false);
        dotted_topology_builder.add_edge(dotted_nodes[3], dotted_nodes[1], (), false);
        dotted_topology_builder.add_edge(dotted_nodes[1], dotted_nodes[2], (), false);
        dotted_topology_builder.add_edge(dotted_nodes[3], dotted_nodes[4], (), false);
        dotted_topology_builder.add_edge(dotted_nodes[4], dotted_nodes[2], (), false);
        let dotted_topology = dotted_topology_builder.build();

        let mut dotted_cff_atom = proper_atom(&dotted_topology);

        dotted_cff_atom = dotted_cff_atom
            .replace(parse!("OSE(6)"))
            .with(parse!("OSE(5)"));

        let mut topology_builder = HedgeGraphBuilder::new();
        let _nodes = (0..4).map(|_| topology_builder.add_node(())).collect_vec();

        topology_builder.add_edge(dotted_nodes[0], dotted_nodes[3], (), false);
        topology_builder.add_edge(dotted_nodes[0], dotted_nodes[2], (), false);
        topology_builder.add_edge(dotted_nodes[0], dotted_nodes[1], (), false);
        topology_builder.add_edge(dotted_nodes[3], dotted_nodes[1], (), false);
        topology_builder.add_edge(dotted_nodes[1], dotted_nodes[2], (), false);
        topology_builder.add_edge(dotted_nodes[3], dotted_nodes[2], (), false);
        let topology = topology_builder.build();

        let mut cff_atom = proper_atom(&topology);
        cff_atom = cff_atom.replace(parse!("OSE(5)")).with(parse!("OSE5"));
        cff_atom = cff_atom.derivative(symbol!("OSE5"));
        cff_atom = cff_atom.replace(parse!("OSE5")).with(parse!("OSE(5)")) / parse!("2*OSE(5)");

        let diff = (cff_atom + dotted_cff_atom).expand();
        //.replace(function!(GS.ose, W_.x_))
        //.with(parse!("E"))
        //.expand();

        println!("diff: {}", diff);
    }

    #[test]
    fn test_cff_generation_triangle() {
        let dep_mom = EdgeIndex::from(2);
        let dep_mom_expr = vec![(EdgeIndex::from(0), -1), (EdgeIndex::from(1), -1)];

        let shift_rewrite = Some(ShiftRewrite {
            dependent_momentum: dep_mom,
            dependent_momentum_expr: dep_mom_expr,
        });

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
            .new_edgevec_from_iter(energy_cache)
            .unwrap();

        let energy_prefactor = virtual_energy_cache
            .iter()
            .map(|e| (F(2.) * e).inv())
            .reduce(|acc, x| acc * x)
            .unwrap();

        // Test the production adapter from a hedge graph.
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

        let triangle_hedge_graph: HedgeGraph<(), (), ()> =
            triangle_hedge_graph_builder.build::<NodeStorageVec<()>>();

        let cff_hedge =
            generate_cff_expression(&triangle_hedge_graph, &shift_rewrite, &[]).unwrap();
        assert_eq!(
            cff_hedge.surfaces.esurface_cache.len(),
            6,
            "too many esurfaces: {:#?}",
            cff_hedge.surfaces.esurface_cache,
        );
        let mut cff_hedge_evaluator = cff_hedge.quick_symbolica_evaluator(0..3, 3..6);

        let cff_res: F<f64> = energy_prefactor
            * cff_hedge_evaluator.evaluate_single(energy_cache.as_ref())
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

    mod failing {
        use super::*;

        #[test]
        fn test_cff_test_double_triangle() {
            let dep_mom = EdgeIndex::from(1);
            let dep_mom_expr = vec![(EdgeIndex::from(0), -1)];

            let shift_rewrite = Some(ShiftRewrite {
                dependent_momentum: dep_mom,
                dependent_momentum_expr: dep_mom_expr,
            });

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

            let mut energy_cache = external_energy_cache.to_vec();
            energy_cache.extend(virtual_energy_cache);

            let energy_cache = dummy_hedge_graph(energy_cache.len())
                .new_edgevec_from_iter(energy_cache)
                .unwrap();

            let energy_prefactor = virtual_energy_cache
                .iter()
                .map(|e| (F(2.) * e).inv())
                .reduce(|acc, x| acc * x)
                .unwrap();

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

            let hedge_double_traingle: HedgeGraph<(), (), ()> =
                hedge_double_triangle_builder.build::<NodeStorageVec<()>>();
            let cff_hedge =
                generate_cff_expression(&hedge_double_traingle, &shift_rewrite, &[]).unwrap();
            let mut cff_hedge_evaluator = cff_hedge.quick_symbolica_evaluator(0..2, 2..7);
            let cff_res =
                energy_prefactor * cff_hedge_evaluator.evaluate_single(energy_cache.as_ref());

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
        }

        #[test]
        fn test_cff_tbt() {
            let dep_mom = EdgeIndex::from(1);
            let dep_mom_expr = vec![(EdgeIndex::from(0), -1)];

            let shift_rewrite = Some(ShiftRewrite {
                dependent_momentum: dep_mom,
                dependent_momentum_expr: dep_mom_expr,
            });

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
                .new_edgevec_from_iter(energies_cache)
                .unwrap();

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

            let tbt_hedge: HedgeGraph<(), (), ()> = tbt_hedge_builder.build::<NodeStorageVec<()>>();
            let cff_hedge = generate_cff_expression(&tbt_hedge, &shift_rewrite, &[]).unwrap();

            let mut cff_hedge_evaluator = cff_hedge.quick_symbolica_evaluator(0..2, 2..10);
            let res =
                cff_hedge_evaluator.evaluate_single(energies_cache.as_ref()) * energy_prefactor;

            let absolute_error = res - F(1.2625322619777278e-21);
            let relative_error = absolute_error / res;
            assert!(
                relative_error.abs() < F(1.0e-15),
                "relative error: {:+e}",
                relative_error
            );
        }
    }
}
