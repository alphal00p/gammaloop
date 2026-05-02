use std::collections::{BTreeMap, HashMap};

use crate::{
    cff::{
        expression::OrientationData,
        hsurface::{Hsurface, HsurfaceID},
        surface::{
            GammaLoopSurfaceCache, HybridSurface, HybridSurfaceID, LinearSurface, LinearSurfaceID,
            LinearSurfaceKind,
        },
        tree::Tree,
    },
    graph::{Graph, get_cff_inverse_energy_product_impl},
    processes::{CrossSectionCut, CutId},
};
use ahash::HashSet;
use color_eyre::Report;
use color_eyre::Result;
use itertools::Itertools;
use linnet::half_edge::{
    HedgeGraph,
    involution::{EdgeVec, Flow, HedgePair},
    subgraph::{OrientedCut, SubGraphLike, SubSetOps},
};
use linnet::half_edge::{
    involution::{EdgeIndex, Orientation},
    subgraph::InternalSubGraph,
};
use symbolica::atom::{Atom, AtomCore};
use three_dimensional_reps::{
    Generate3DExpressionOptions, NumeratorSamplingScaleMode, RepresentationMode, ThreeDGraphSource,
};
use typed_index_collections::TiVec;

use serde::{Deserialize, Serialize};

use tracing::debug;

use super::{
    cff_graph::CFFGenerationGraph,
    esurface::{Esurface, EsurfaceCollection, EsurfaceID, ExternalShift},
    expression::{CFFExpression, OrientationExpression, OrientationID},
    hsurface::HsurfaceCollection,
    surface::HybridSurfaceRef,
};

pub use super::surface::SurfaceCache;

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
#[allow(dead_code)]
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

#[cfg(test)]
fn get_orientations<E, V, H>(
    graph: &HedgeGraph<E, V, H>,
    dummy_edges: &[EdgeIndex],
) -> Vec<CFFGenerationGraph> {
    let internal_subgraph = InternalSubGraph::cleaned_filter_pessimist(graph.full_filter(), graph);
    let num_virtual_edges = graph.count_internal_edges(&internal_subgraph);
    let virtual_possible_orientations = iterate_possible_orientations(num_virtual_edges);

    virtual_possible_orientations
        .map(|orientation_of_virtuals| {
            let mut orientation_of_virtuals = orientation_of_virtuals.into_iter();

            let global_orientation = graph.new_edgevec(|_, __, hedge_pair| match hedge_pair {
                HedgePair::Unpaired { .. } => Orientation::Default,
                HedgePair::Paired { .. } => orientation_of_virtuals
                    .next()
                    .expect(" unable to reconstruct orientation"),
                HedgePair::Split { .. } => todo!(),
            });

            assert!(
                orientation_of_virtuals.next().is_none(),
                "did not saturate virtual orientations when constructing global orientation"
            );

            CFFGenerationGraph::new(graph, global_orientation, dummy_edges)
        })
        .collect_vec()
}

pub(crate) fn get_orientations_from_subgraph<E, V, H, S: SubGraphLike>(
    graph: &HedgeGraph<E, V, H>,
    subgraph: &S,
    reversed_dangling: &[EdgeIndex],
) -> Vec<CFFGenerationGraph> {
    let num_virtual_edges = graph.count_internal_edges(subgraph);
    let virtual_possible_orientations = iterate_possible_orientations(num_virtual_edges);

    virtual_possible_orientations
        .map(|orientation_of_virtuals| {
            let mut orientation_of_virtuals = orientation_of_virtuals.into_iter();

            let global_orientation = graph.new_edgevec(|_, edge_id, _| {
                if let Some((pair, _, _)) = graph
                    .iter_edges_of(subgraph)
                    .find(|(_pair, id, _)| *id == edge_id)
                {
                    match pair {
                        HedgePair::Paired { .. } => orientation_of_virtuals
                            .next()
                            .expect("orientation generation corrupted, not enough edges"),
                        HedgePair::Unpaired { .. } => Orientation::Default,
                        HedgePair::Split { .. } => {
                            if reversed_dangling.contains(&edge_id) {
                                Orientation::Reversed
                            } else {
                                Orientation::Default
                            }
                        }
                    }
                } else {
                    Orientation::Undirected
                }
            });

            CFFGenerationGraph::new_from_subgraph(graph, global_orientation, subgraph).unwrap()
        })
        .filter(|cff_graph| !cff_graph.has_directed_cycle_initial())
        .collect()
}

#[allow(unused)]
fn get_orientations_with_cut<E, V, H>(
    graph: &HedgeGraph<E, V, H>,
    oriented_cut: &OrientedCut,
) -> Vec<EdgeVec<Orientation>> {
    let internal_subgraph = InternalSubGraph::cleaned_filter_pessimist(graph.full_filter(), graph);
    let num_virtual_edges = graph.count_internal_edges(&internal_subgraph);

    let virtual_possible_orientations = iterate_possible_orientations(num_virtual_edges);

    let orientations_consistent_with_cut = virtual_possible_orientations
        .map(|orientation_of_virtuals| {
            // pad a virtual orientation with orientations of externals.
            let mut orientation_of_virtuals = orientation_of_virtuals.into_iter();

            let global_orientation = graph.new_edgevec(|_, __, hedge_pair| match hedge_pair {
                HedgePair::Unpaired { .. } => Orientation::Default,
                HedgePair::Paired { .. } => orientation_of_virtuals
                    .next()
                    .expect(" unable to reconstruct orientation"),
                HedgePair::Split { .. } => todo!(),
            });

            assert!(
                orientation_of_virtuals.next().is_none(),
                "did not saturate virtual orientations when constructing global orientation"
            );

            global_orientation
        })
        .filter(|global_orientation| {
            // filter out orientations that are not consistent with the cut
            let edges_in_cut = graph.iter_edges_of(oriented_cut).map(|(_, id, _)| id);
            let orientation_of_edges_in_cut = oriented_cut.iter_edges(graph).map(|(or, _)| or);

            edges_in_cut
                .zip(orientation_of_edges_in_cut)
                .all(|(edge_id, orientation)| global_orientation[edge_id] == orientation)
        })
        .filter(|global_orientation| {
            // filter out orientations that have a directed cycle
            let graph = CFFGenerationGraph::new(graph, global_orientation.clone(), &[]);
            !graph.has_directed_cycle_initial()
        });

    orientations_consistent_with_cut.collect_vec()
}

#[cfg(test)]
fn generate_cff_expression<E, V, H>(
    graph: &HedgeGraph<E, V, H>,
    canonize_esurface: &Option<ShiftRewrite>,
    edges_in_initial_state_cut: &[EdgeIndex],
    dummy_edges: &[EdgeIndex],
) -> Result<CFFExpression<OrientationID>> {
    let graphs = get_orientations(graph, dummy_edges);
    debug!("number of orientations: {}", graphs.len());
    let mut surface_cache = SurfaceCache {
        esurface_cache: EsurfaceCollection::from_iter(std::iter::empty()),
        hsurface_cache: HsurfaceCollection::from_iter(std::iter::empty()),
        linear_surface_cache: crate::cff::surface::LinearSurfaceCollection::from_iter(
            std::iter::empty(),
        ),
    };
    let graph_cff = generate_cff_from_orientations(
        graphs,
        &mut surface_cache,
        edges_in_initial_state_cut,
        canonize_esurface,
    )?;

    // patch the surface cache
    Ok(graph_cff)
}

impl Graph {
    pub(crate) fn generate_3d_expression_for_integrand(
        &mut self,
        contract_edges: &[EdgeIndex],
        canonize_esurface: &Option<ShiftRewrite>,
    ) -> Result<CFFExpression<OrientationID>> {
        let initial_state_cut_edges = self
            .iter_edges_of(&self.initial_state_cut)
            .map(|(_, edge_id, _)| edge_id)
            .collect_vec();
        if !initial_state_cut_edges.is_empty() {
            debug!(
                "using legacy CFF generation for graphs with initial-state cuts until generalized initial-cut parity is complete"
            );
            return self.generate_cff(contract_edges, canonize_esurface);
        }

        if !contract_edges.is_empty() {
            debug!(
                "using legacy CFF generation for contracted production subgraphs until generalized contracted-basis support is complete"
            );
            return self.generate_cff(contract_edges, canonize_esurface);
        }

        let options = Generate3DExpressionOptions {
            representation: RepresentationMode::Cff,
            energy_degree_bounds: Vec::new(),
            numerator_sampling_scale: NumeratorSamplingScaleMode::None,
        };
        let mut expression = match {
            let source = ContractedGraphSource {
                graph: self,
                contract_edges,
            };
            three_dimensional_reps::generate_3d_expression(&source, &options)
        } {
            Ok(expression) => expression,
            Err(error) => {
                debug!(
                    ?error,
                    "falling back to legacy CFF generation for a contracted graph not yet supported by generalized 3D expression generation"
                );
                return self.generate_cff(contract_edges, canonize_esurface);
            }
        };

        if !initial_state_cut_edges.is_empty() {
            expression.orientations = expression
                .orientations
                .into_iter()
                .filter(|orientation| {
                    initial_state_cut_edges.iter().all(|edge_id| {
                        orientation
                            .data
                            .orientation
                            .get(*edge_id)
                            .is_some_and(|orientation| *orientation == Orientation::Default)
                    })
                })
                .collect();
        }

        self.convert_generated_expression_surfaces(
            expression,
            canonize_esurface,
            &initial_state_cut_edges,
        )
    }

    #[allow(dead_code)]
    pub(crate) fn generate_cff(
        &mut self,
        contract_edges: &[EdgeIndex],
        canonize_esurface: &Option<ShiftRewrite>,
    ) -> Result<CFFExpression<OrientationID>> {
        let mut seed_graph = CFFGenerationGraph::new_from_graph(self);

        for edge in contract_edges {
            seed_graph = seed_graph.contract_edge(*edge);
        }

        let edges_in_initial_state_cut = self
            .iter_edges_of(&self.initial_state_cut)
            .map(|x| x.1)
            .collect_vec();

        let virtual_edges_of_contracted_graph = seed_graph.num_virtual_edges();

        let orientations = iterate_possible_orientations(virtual_edges_of_contracted_graph);

        let mut oriented_acyclic_graphs = vec![];

        for orientation in orientations {
            let mut orientation_iterator = orientation.into_iter();

            let global_orientation = self.new_edgevec(|_, edge_id, hedge_pair| {
                if hedge_pair.is_unpaired() || contract_edges.contains(&edge_id) {
                    Orientation::Undirected
                } else if edges_in_initial_state_cut.contains(&edge_id) {
                    Orientation::Default
                } else {
                    orientation_iterator
                        .next()
                        .expect("orientation generation corrupted, not enough edges")
                }
            });

            let mut cff_graph = seed_graph.clone();
            cff_graph.apply_orientation(global_orientation)?;

            if !cff_graph.has_directed_cycle_initial() {
                oriented_acyclic_graphs.push(cff_graph);
            }
        }

        generate_cff_from_orientations(
            oriented_acyclic_graphs,
            &mut self.surface_cache,
            &edges_in_initial_state_cut,
            canonize_esurface,
        )
    }

    fn convert_generated_expression_surfaces(
        &mut self,
        mut expression: three_dimensional_reps::ThreeDExpression<OrientationID>,
        canonize_esurface: &Option<ShiftRewrite>,
        initial_state_cut_edges: &[EdgeIndex],
    ) -> Result<CFFExpression<OrientationID>> {
        let mut linear_surface_map = BTreeMap::<LinearSurfaceID, HybridSurfaceID>::new();
        for (linear_surface_id, surface) in
            expression.surfaces.linear_surface_cache.iter_enumerated()
        {
            let converted = self.intern_generated_linear_surface(
                surface,
                canonize_esurface,
                initial_state_cut_edges,
            )?;
            linear_surface_map.insert(linear_surface_id, converted);
        }

        for orientation in expression.orientations.iter_mut() {
            orientation.for_each_denominator_tree_mut(|denominator| {
                denominator.map_mut(|surface_id| {
                    remap_generated_surface_id(surface_id, &linear_surface_map)
                });
            });
            for variant in &mut orientation.variants {
                for surface_id in &mut variant.numerator_surfaces {
                    remap_generated_surface_id(surface_id, &linear_surface_map);
                }
                // GammaLoop production CFF keeps inverse on-shell-energy factors
                // outside the CFF denominator tree. The generalized crate stores
                // them per variant, so strip them here to preserve the current
                // production convention.
                variant.half_edges.clear();
            }
        }

        Ok(CFFExpression {
            orientations: expression.orientations,
            surfaces: self.surface_cache.clone(),
        })
    }

    fn intern_generated_linear_surface(
        &mut self,
        surface: &LinearSurface,
        canonize_esurface: &Option<ShiftRewrite>,
        initial_state_cut_edges: &[EdgeIndex],
    ) -> Result<HybridSurfaceID> {
        if !surface.expression.uniform_scale_coeff.is_zero()
            || !surface.expression.constant.is_zero()
        {
            return Err(eyre::eyre!(
                "generalized CFF production cannot convert non-homogeneous linear surface {:?}",
                surface
            ));
        }

        let mut positive_energies = Vec::new();
        let mut negative_energies = Vec::new();
        let mut external_shift = Vec::new();
        collect_linear_surface_terms(
            &surface.expression.internal_terms,
            initial_state_cut_edges,
            &mut positive_energies,
            &mut negative_energies,
            &mut external_shift,
        )?;
        for (edge_id, coeff) in &surface.expression.external_terms {
            let coeff = atom_integer_coeff(coeff)?;
            external_shift.push((*edge_id, coeff));
        }
        positive_energies.sort();
        negative_energies.sort();
        external_shift.sort_by_key(|(edge_id, _)| *edge_id);

        match surface.kind {
            LinearSurfaceKind::Esurface => {
                if !negative_energies.is_empty() {
                    return Err(eyre::eyre!(
                        "generalized CFF production cannot convert E-surface with negative internal terms {:?}",
                        surface
                    ));
                }
                let mut esurface = Esurface {
                    energies: positive_energies,
                    external_shift,
                    vertex_set: crate::cff::cff_graph::VertexSet::dummy(),
                };
                if let Some(shift_rewrite) = canonize_esurface {
                    esurface.canonicalize_shift(shift_rewrite);
                }

                let esurface_id = self
                    .surface_cache
                    .esurface_cache
                    .position(|existing| *existing == esurface)
                    .unwrap_or_else(|| {
                        self.surface_cache.esurface_cache.push(esurface);
                        Into::<EsurfaceID>::into(self.surface_cache.esurface_cache.len() - 1)
                    });

                Ok(HybridSurfaceID::Esurface(esurface_id))
            }
            LinearSurfaceKind::Hsurface => {
                let hsurface = Hsurface {
                    positive_energies,
                    negative_energies,
                    external_shift,
                    vertex_set: crate::cff::cff_graph::VertexSet::dummy(),
                };
                let hsurface_id = self
                    .surface_cache
                    .hsurface_cache
                    .position(|existing| existing == &hsurface)
                    .unwrap_or_else(|| {
                        self.surface_cache.hsurface_cache.push(hsurface);
                        Into::<HsurfaceID>::into(self.surface_cache.hsurface_cache.len() - 1)
                    });

                Ok(HybridSurfaceID::Hsurface(hsurface_id))
            }
        }
    }
}

struct ContractedGraphSource<'a> {
    graph: &'a Graph,
    contract_edges: &'a [EdgeIndex],
}

impl ThreeDGraphSource for ContractedGraphSource<'_> {
    fn to_three_d_parsed_graph(
        &self,
    ) -> three_dimensional_reps::graph_io::Result<three_dimensional_reps::ParsedGraph> {
        use three_dimensional_reps::graph_io::{
            GraphIoError, ParsedGraph, ParsedGraphExternalEdge, ParsedGraphInternalEdge,
        };

        let mut parent = (0..self.graph.n_nodes()).collect_vec();
        for (pair, edge_id, edge_data) in self.graph.underlying.iter_edges() {
            if edge_data.data.is_dummy || !self.contract_edges.contains(&edge_id) {
                continue;
            }
            if let HedgePair::Paired { source, sink } = pair {
                union_parent(
                    &mut parent,
                    usize::from(self.graph.node_id(source)),
                    usize::from(self.graph.node_id(sink)),
                );
            }
        }

        let node_ids = self
            .graph
            .underlying
            .iter_nodes()
            .map(|(node_id, _, _)| node_id)
            .sorted()
            .collect::<Vec<_>>();
        let mut root_to_internal = BTreeMap::<usize, usize>::new();
        for node_id in &node_ids {
            let root = find_parent(&mut parent, usize::from(*node_id));
            if !root_to_internal.contains_key(&root) {
                root_to_internal.insert(root, root_to_internal.len());
            }
        }
        let mut node_to_internal = BTreeMap::new();
        for node_id in &node_ids {
            let root = find_parent(&mut parent, usize::from(*node_id));
            node_to_internal.insert(*node_id, root_to_internal[&root]);
        }

        let loop_names = self
            .graph
            .loop_momentum_basis
            .loop_edges
            .iter()
            .map(|edge_id| self.graph.underlying[*edge_id].name.value.clone())
            .collect::<Vec<_>>();
        let external_names = self
            .graph
            .loop_momentum_basis
            .ext_edges
            .iter()
            .map(|edge_id| self.graph.underlying[*edge_id].name.value.clone())
            .collect::<Vec<_>>();

        let mut internal_edges = Vec::new();
        let mut external_edges = Vec::new();
        let mut next_external_id = 10_000_000usize;

        for (pair, edge_index, edge_data) in self
            .graph
            .underlying
            .iter_edges()
            .sorted_by_key(|(_, edge_index, _)| *edge_index)
        {
            if edge_data.data.is_dummy {
                continue;
            }
            let signature = &self.graph.loop_momentum_basis.edge_signatures[edge_index];
            let momentum_signature = three_dimensional_reps::MomentumSignature {
                loop_signature: (&signature.internal).into_iter().map(sign_to_i32).collect(),
                external_signature: (&signature.external).into_iter().map(sign_to_i32).collect(),
            };
            let label = edge_data.data.name.value.clone();
            match pair {
                HedgePair::Paired { source, sink } => {
                    if self.contract_edges.contains(&edge_index) {
                        continue;
                    }
                    let tail = *node_to_internal
                        .get(&self.graph.node_id(source))
                        .ok_or_else(|| {
                            GraphIoError::Source(format!(
                                "missing contracted source node mapping for edge {}",
                                usize::from(edge_index)
                            ))
                        })?;
                    let head =
                        *node_to_internal
                            .get(&self.graph.node_id(sink))
                            .ok_or_else(|| {
                                GraphIoError::Source(format!(
                                    "missing contracted sink node mapping for edge {}",
                                    usize::from(edge_index)
                                ))
                            })?;
                    internal_edges.push(ParsedGraphInternalEdge {
                        edge_id: internal_edges.len(),
                        tail,
                        head,
                        label,
                        mass_key: Some(edge_data.data.particle.mass_atom().to_canonical_string()),
                        signature: momentum_signature,
                        had_pow: false,
                    });
                }
                HedgePair::Unpaired { hedge, flow } => {
                    let node = *node_to_internal
                        .get(&self.graph.node_id(hedge))
                        .ok_or_else(|| {
                            GraphIoError::Source(format!(
                                "missing contracted external node mapping for edge {}",
                                usize::from(edge_index)
                            ))
                        })?;
                    let (source, destination) = match flow {
                        Flow::Source => (Some(node), None),
                        Flow::Sink => (None, Some(node)),
                    };
                    external_edges.push(ParsedGraphExternalEdge {
                        edge_id: next_external_id,
                        source,
                        destination,
                        label,
                        external_coefficients: momentum_signature.external_signature,
                    });
                    next_external_id += 1;
                }
                HedgePair::Split { .. } => {
                    return Err(GraphIoError::Source(
                        "split edges are not supported when extracting contracted Graph input"
                            .to_string(),
                    ));
                }
            }
        }

        Ok(ParsedGraph {
            internal_edges,
            external_edges,
            loop_names,
            external_names,
            node_name_to_internal: root_to_internal
                .into_iter()
                .map(|(root, node)| (format!("n{root}"), node))
                .collect(),
        })
    }

    fn energy_edge_index_map(
        &self,
        _parsed: &three_dimensional_reps::ParsedGraph,
    ) -> Option<three_dimensional_reps::EnergyEdgeIndexMap> {
        let internal = self
            .graph
            .underlying
            .iter_edges()
            .sorted_by_key(|(_, edge_index, _)| *edge_index)
            .filter_map(|(pair, edge_index, edge_data)| {
                (pair.is_paired()
                    && !edge_data.data.is_dummy
                    && !self.contract_edges.contains(&edge_index))
                .then_some(usize::from(edge_index))
            })
            .enumerate()
            .collect::<BTreeMap<_, _>>();

        let external = self
            .graph
            .loop_momentum_basis
            .ext_edges
            .iter()
            .enumerate()
            .map(|(external_id, edge_id)| (external_id, usize::from(*edge_id)))
            .collect::<BTreeMap<_, _>>();

        Some(three_dimensional_reps::EnergyEdgeIndexMap {
            internal,
            external,
            orientation_edge_count: self.graph.underlying.n_edges(),
        })
    }
}

fn remap_generated_surface_id(
    surface_id: &mut HybridSurfaceID,
    linear_surface_map: &BTreeMap<LinearSurfaceID, HybridSurfaceID>,
) {
    if let HybridSurfaceID::Linear(linear_surface_id) = surface_id {
        *surface_id = *linear_surface_map
            .get(linear_surface_id)
            .expect("all generated linear surfaces should have been interned");
    }
}

fn collect_linear_surface_terms(
    terms: &[(EdgeIndex, Atom)],
    initial_state_cut_edges: &[EdgeIndex],
    positive_energies: &mut Vec<EdgeIndex>,
    negative_energies: &mut Vec<EdgeIndex>,
    external_shift: &mut Vec<(EdgeIndex, i64)>,
) -> Result<()> {
    for (edge_id, coeff) in terms {
        let coeff = atom_integer_coeff(coeff)?;
        if coeff == 0 {
            continue;
        }
        if initial_state_cut_edges.contains(edge_id) {
            // GammaLoop cut E-surfaces store initial-state energies on the
            // external-shift side with the opposite sign.
            external_shift.push((*edge_id, -coeff));
            continue;
        }

        let target = if coeff > 0 {
            &mut *positive_energies
        } else {
            &mut *negative_energies
        };
        for _ in 0..coeff.unsigned_abs() {
            target.push(*edge_id);
        }
    }
    Ok(())
}

fn atom_integer_coeff(coeff: &Atom) -> Result<i64> {
    let coeff_text = coeff.to_canonical_string();
    coeff_text
        .parse::<i64>()
        .map_err(|_| eyre::eyre!("expected integer linear-surface coefficient, found {coeff_text}"))
}

fn find_parent(parent: &mut [usize], node: usize) -> usize {
    let parent_node = parent[node];
    if parent_node == node {
        node
    } else {
        let root = find_parent(parent, parent_node);
        parent[node] = root;
        root
    }
}

fn union_parent(parent: &mut [usize], left: usize, right: usize) {
    let left_root = find_parent(parent, left);
    let right_root = find_parent(parent, right);
    if left_root != right_root {
        parent[right_root] = left_root;
    }
}

fn sign_to_i32(sign: crate::momentum::SignOrZero) -> i32 {
    match sign {
        crate::momentum::SignOrZero::Minus => -1,
        crate::momentum::SignOrZero::Zero => 0,
        crate::momentum::SignOrZero::Plus => 1,
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
    let graphs = get_orientations_from_subgraph(graph, subgraph, reversed_dangling);
    let cff = generate_cff_from_orientations(
        graphs,
        surface_cache,
        edges_in_initial_state_cut,
        canonize_esurface,
    )?;
    Ok(cff)
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
        linear_surface_cache: crate::cff::surface::LinearSurfaceCollection::from_iter(
            std::iter::empty(),
        ),
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
        surface_cache_to_use.substitute_energies_gs(&atom_tree, topology.cut_edges);
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
                HybridSurfaceID::Linear(_) => true,
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
                HybridSurfaceRef::Linear(_) => {
                    id_map.insert(*appearing_id, *appearing_id);
                }
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

fn generate_cff_from_orientations<O: From<usize> + Into<usize>>(
    orientations_and_graphs: Vec<CFFGenerationGraph>,
    generator_cache: &mut SurfaceCache,
    edges_in_initial_state_cut: &[EdgeIndex],
    canonize_esurface: &Option<ShiftRewrite>,
) -> Result<CFFExpression<O>, Report> {
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
                generator_cache,
                edges_in_initial_state_cut,
                canonize_esurface,
            );
            let expression = tree.map(forget_graphs);

            OrientationExpression::pure_cff(OrientationData::new(global_orientation), expression)
        })
        .collect_vec();

    Ok(CFFExpression {
        orientations: terms.into(),
        surfaces: generator_cache.clone(),
    })
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
                HybridSurface::Linear(surface) => HybridSurface::Linear(surface),
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
                HybridSurface::Linear(surface) => {
                    let option_surface_id = generator_cache
                        .linear_surface_cache
                        .position(|val| val == &surface);

                    let surface_id = match option_surface_id {
                        Some(surface_id) => surface_id,
                        None => {
                            generator_cache.linear_surface_cache.push(surface);
                            Into::<crate::cff::surface::LinearSurfaceID>::into(
                                generator_cache.linear_surface_cache.len() - 1,
                            )
                        }
                    };

                    HybridSurfaceID::Linear(surface_id)
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
    use std::{ops::Range, vec};

    use ahash::HashMap;

    use linnet::half_edge::{
        builder::HedgeGraphBuilder, involution::Flow, nodestore::NodeStorageVec,
    };
    use symbolica::{
        atom::AtomCore,
        evaluate::{ExpressionEvaluator, FunctionMap, OptimizationSettings},
        parse, symbol,
    };
    use utils::FloatLike;

    use crate::{
        cff::{cff_graph::CFFEdgeType, expression::GammaLoopThreeDExpression},
        dot,
        graph::{feynman_graph::FeynmanGraph, parse::IntoGraph},
        initialisation::test_initialise,
        momentum::{FourMomentum, ThreeMomentum},
        processes::AmplitudeGraph,
        settings::global::OrientationPattern,
        utils::{
            self, F, RefDefault, external_energy_atom_from_index, ose_atom_from_index,
            test_utils::dummy_hedge_graph,
        },
    };

    use super::*;

    trait QuickCffSymbolicaEvaluator {
        fn quick_symbolica_evaluator(
            &self,
            external_range: Range<usize>,
            virtual_range: Range<usize>,
        ) -> ExpressionEvaluator<F<f64>>;
    }

    // helper function to make a symbolica evaluator
    impl QuickCffSymbolicaEvaluator for CFFExpression<OrientationID> {
        fn quick_symbolica_evaluator(
            &self,
            external_range: Range<usize>,
            virtual_range: Range<usize>,
        ) -> ExpressionEvaluator<F<f64>> {
            let expression_atom_no_energy_sub = self.to_atom_gs(&OrientationPattern::default());
            let num_energies = external_range.end.max(virtual_range.end);
            let mut expression_atom = self
                .surfaces
                .substitute_energies_gs(&expression_atom_no_energy_sub, &[]);
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
                .to_evaluation_tree(&function_map, &params)
                .unwrap();

            tree.horner_scheme();
            tree.common_subexpression_elimination();
            tree.linearize(&OptimizationSettings::default())
                .map_coeff(&|c| (&c.re).into())
        }
    }

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

        let mut surface_cache = SurfaceCache::new();

        let _cff = generate_cff_from_orientations::<OrientationID>(
            orientations,
            &mut surface_cache,
            &[],
            &Some(shift_rewrite),
        )
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

        let mut surface_cache = SurfaceCache::new();

        let _cff = generate_cff_from_orientations::<OrientationID>(
            orientations,
            &mut surface_cache,
            &[],
            &Some(shift_rewrite),
        )
        .unwrap();

        let _finish = std::time::Instant::now();
    }

    fn proper_atom(graph: &HedgeGraph<(), ()>) -> Atom {
        let cff = generate_cff_expression(graph, &None, &[], &[]).unwrap();

        let mut cff_atom = cff.to_atom_gs(&OrientationPattern::default());
        cff_atom = cff.surfaces.substitute_energies_gs(&cff_atom, &[]);
        let inverse_energy_product =
            get_cff_inverse_energy_product_impl(graph, &graph.full_graph(), &[]);

        cff_atom *= inverse_energy_product;
        cff_atom
    }

    fn cff_atom_with_surface_substitutions(cff: &CFFExpression<OrientationID>) -> Atom {
        let atom = cff.to_atom_gs(&OrientationPattern::default());
        cff.surfaces.substitute_energies_gs(&atom, &[])
    }

    fn amplitude_production_contract_edges(graph: &Graph) -> Vec<EdgeIndex> {
        graph
            .iter_edges_of(
                &graph
                    .tree_edges
                    .subtract(&graph.initial_state_cut)
                    .subtract(&graph.external_filter::<linnet::half_edge::subgraph::SuBitGraph>()),
            )
            .map(|(_, edge_id, _)| edge_id)
            .collect_vec()
    }

    #[test]
    fn generated_cff_matches_legacy_cff_for_scalar_bubble() {
        test_initialise().unwrap();
        let graph: AmplitudeGraph = dot!(
            digraph bub {
                edge [particle=scalar_1]
                node [num=1]
                e [style=invis]
                e -> A:0 [id=3]
                B:1 -> e [id=2]
                A -> B [id=1]
                A -> B [id=0]
            },
            "scalars"
        )
        .unwrap();

        for (case, contract_edges) in [
            ("uncontracted", Vec::new()),
            (
                "amplitude production contractions",
                amplitude_production_contract_edges(&graph.graph),
            ),
        ] {
            let mut legacy_graph = graph.graph.clone();
            let legacy_canonize =
                legacy_graph.get_esurface_canonization(&legacy_graph.loop_momentum_basis);
            let legacy = legacy_graph
                .generate_cff(&contract_edges, &legacy_canonize)
                .unwrap();

            let mut generated_graph = graph.graph.clone();
            let generated_canonize =
                generated_graph.get_esurface_canonization(&generated_graph.loop_momentum_basis);
            let generated = generated_graph
                .generate_3d_expression_for_integrand(&contract_edges, &generated_canonize)
                .unwrap();

            let legacy_orientations = legacy
                .orientations
                .iter()
                .map(|orientation| orientation.data.to_string())
                .sorted()
                .collect_vec();
            let generated_orientations = generated
                .orientations
                .iter()
                .map(|orientation| orientation.data.to_string())
                .sorted()
                .collect_vec();
            assert_eq!(legacy_orientations, generated_orientations, "{case}");

            let diff = (cff_atom_with_surface_substitutions(&legacy)
                - cff_atom_with_surface_substitutions(&generated))
            .expand();
            assert!(
                diff.is_zero(),
                "generalized CFF path differs from legacy path for {case} by {diff}"
            );
        }
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

        let mut surface_cache = SurfaceCache::new();

        let cff = generate_cff_from_orientations(
            orientations,
            &mut surface_cache,
            &[],
            &shift_rewrite.clone(),
        )
        .unwrap();
        assert_eq!(
            cff.surfaces.esurface_cache.len(),
            6,
            "too many esurfaces: {:#?}",
            cff.surfaces.esurface_cache,
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
            .new_edgevec_from_iter(energy_cache)
            .unwrap();

        let energy_prefactor = virtual_energy_cache
            .iter()
            .map(|e| (F(2.) * e).inv())
            .reduce(|acc, x| acc * x)
            .unwrap();

        let mut evaluator = cff.quick_symbolica_evaluator(0..3, 3..6);

        let cff_res: F<f64> = energy_prefactor
            * evaluator.evaluate_single(energy_cache.clone().as_ref())
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

        let triangle_hedge_graph: HedgeGraph<(), (), ()> =
            triangle_hedge_graph_builder.build::<NodeStorageVec<()>>();

        let cff_hedge =
            generate_cff_expression(&triangle_hedge_graph, &shift_rewrite, &[], &[]).unwrap();
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

            let mut surface_cache = SurfaceCache::new();

            let cff = generate_cff_from_orientations(
                orientations,
                &mut surface_cache,
                &[],
                &shift_rewrite,
            )
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

            let mut evaluator = cff.quick_symbolica_evaluator(0..2, 2..7);

            let cff_res =
                energy_prefactor * evaluator.evaluate_single(energy_cache.clone().as_ref());

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

            let hedge_double_traingle: HedgeGraph<(), (), ()> =
                hedge_double_triangle_builder.build::<NodeStorageVec<()>>();
            let cff_hedge =
                generate_cff_expression(&hedge_double_traingle, &shift_rewrite, &[], &[]).unwrap();
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

            let node_3 = hedge_double_traingle.iter_crown(nodes[3]).into();
            let node_0 = hedge_double_traingle.iter_crown(nodes[0]).into();

            let cuts = hedge_double_traingle.all_cuts(node_3, node_0);
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
            let mut surface_cache = SurfaceCache::new();
            let orientataions = generate_orientations_for_testing(tbt_edges, incoming_vertices);
            let cff = generate_cff_from_orientations(
                orientataions,
                &mut surface_cache,
                &[],
                &shift_rewrite,
            )
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
                .new_edgevec_from_iter(energies_cache)
                .unwrap();

            let mut evaluator = cff.quick_symbolica_evaluator(0..2, 2..10);

            let res = evaluator.evaluate_single(energies_cache.clone().as_ref()) * energy_prefactor;

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

            let tbt_hedge: HedgeGraph<(), (), ()> = tbt_hedge_builder.build::<NodeStorageVec<()>>();
            let cff_hedge = generate_cff_expression(&tbt_hedge, &shift_rewrite, &[], &[]).unwrap();

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

            let node_0 = tbt_hedge.iter_crown(nodes[0]).into();
            let node_5 = tbt_hedge.iter_crown(nodes[5]).into();

            let cuts = tbt_hedge.all_cuts(node_0, node_5).clone();
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
    }

    mod slow {
        use super::*;

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
            let mut surface_cache = SurfaceCache::new();

            let cff = generate_cff_from_orientations(
                orientations,
                &mut surface_cache,
                &[],
                &Some(shift_rewrite),
            )
            .unwrap();
            let finish = std::time::Instant::now();
            println!("time to generate cff: {:?}", finish - start);

            let energy_cache = [F(3.0); 17];

            //let energy_cache = dummy_hedge_graph(energy_cache.len())
            //    .new_hedgevec_from_iter(energy_cache)
            //    .unwrap();

            let mut evaluator = cff.quick_symbolica_evaluator(0..4, 4..17);
            let start = std::time::Instant::now();
            for _ in 0..100 {
                let _res = evaluator.evaluate_single(&energy_cache);
            }
            let finish = std::time::Instant::now();
            println!("time to evaluate cff: {:?}", (finish - start) / 100);
        }
    }
}
