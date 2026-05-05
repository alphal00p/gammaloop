use std::collections::BTreeMap;

use crate::{
    cff::{
        VertexSet,
        hsurface::{Hsurface, HsurfaceID},
        surface::{HybridSurfaceID, LinearSurface, LinearSurfaceID, LinearSurfaceKind},
    },
    graph::{Graph, GraphThreeDSource},
    settings::global::{GenerationSettings, ThreeDRepresentation, UniformNumeratorSamplingScale},
};
use ahash::HashSet;
use color_eyre::Result;
use itertools::Itertools;
use linnet::half_edge::{
    involution::EdgeIndex,
    subgraph::{SuBitGraph, SubSetOps},
};
use symbolica::atom::{Atom, AtomCore};
use three_dimensional_reps::{
    Generate3DExpressionOptions, NumeratorSamplingScaleMode, RepresentationMode, ThreeDGraphSource,
    tree::{NodeId, Tree},
};

use tracing::debug;

use super::{
    esurface::{Esurface, EsurfaceID, ExternalShift},
    expression::{CFFExpression, OrientationID},
};

#[derive(Debug, Clone)]
pub struct ShiftRewrite {
    pub dependent_momentum: EdgeIndex,
    pub dependent_momentum_expr: ExternalShift,
}

#[derive(Debug, Clone, Copy)]
struct SurfaceMapEntry {
    surface_id: HybridSurfaceID,
    sign: i64,
}

impl Graph {
    pub(crate) fn generate_3d_expression_for_integrand(
        &mut self,
        contract_edges: &[EdgeIndex],
        canonize_esurface: &Option<ShiftRewrite>,
        options: &Generate3DExpressionOptions,
    ) -> Result<CFFExpression<OrientationID>> {
        let initial_state_cut_edges = self
            .iter_edges_of(&self.initial_state_cut)
            .map(|(_, edge_id, _)| edge_id)
            .collect_vec();
        let preserved_edges = self.normalized_preserved_4d_denominator_edges(options);
        let source_contract_edges = source_contract_edges_for_3d_expression(
            contract_edges,
            &initial_state_cut_edges,
            &preserved_edges,
        );
        if !source_contract_edges.is_empty() {
            debug!(
                "contracting {} internal edges before generalized 3D expression generation",
                source_contract_edges.len()
            );
        }
        if !preserved_edges.is_empty() {
            debug!(
                "preserving {} internal edges as four-dimensional residual denominators",
                preserved_edges.len()
            );
        }
        let mut local_options = options.clone();
        local_options.preserve_internal_edges_as_four_d_denominators = preserved_edges
            .iter()
            .map(|edge_id| usize::from(*edge_id))
            .collect();

        let source = GraphThreeDSource::new(self, &source_contract_edges);

        let expression = {
            three_dimensional_reps::generate_3d_expression(&source, &local_options).map_err(
                |error| {
                    let source_summary = source
                        .to_three_d_parsed_graph()
                        .map(|parsed| three_d_source_summary(&parsed))
                        .unwrap_or_else(|source_error| {
                            format!("failed to rebuild 3D source summary: {source_error}")
                        });
                    eyre::eyre!(
                        "generalized 3D expression generation failed for graph `{}` with representation {:?} and energy-degree bounds {:?}: {error}\n{source_summary}",
                        self.name,
                        local_options.representation,
                        local_options.energy_degree_bounds,
                    )
                },
            )
        }?;

        self.convert_generated_expression_surfaces(
            expression,
            local_options.representation,
            generated_cff_expression_uses_variant_half_edges(&local_options),
            canonize_esurface,
            &initial_state_cut_edges,
        )
    }

    pub(crate) fn production_3d_expression_options(
        &self,
        representation: ThreeDRepresentation,
        settings: &GenerationSettings,
    ) -> Result<Generate3DExpressionOptions> {
        let mut options = self.three_d_expression_options(
            representation_mode(representation),
            numerator_sampling_scale_mode(settings.uniform_numerator_sampling_scale),
        )?;
        if representation == ThreeDRepresentation::Cff && !settings.explicit_orientation_sum_only {
            // The orientation-localized CFF forest path substitutes the
            // numerator pointwise on each orientation and only needs the pure
            // CFF denominator expression. The bounded high-energy CFF root is
            // needed when production builds one explicitly summed 3D expression
            // up front.
            options.energy_degree_bounds.clear();
        }
        Ok(options)
    }

    pub fn three_d_expression_options(
        &self,
        representation: RepresentationMode,
        numerator_sampling_scale: NumeratorSamplingScaleMode,
    ) -> Result<Generate3DExpressionOptions> {
        Ok(Generate3DExpressionOptions {
            representation,
            energy_degree_bounds: self
                .automatic_numerator_energy_degree_bounds_for_3d_expression(representation)?,
            numerator_sampling_scale,
            include_cff_duplicate_signature_excess_sign: true,
            preserve_internal_edges_as_four_d_denominators: self
                .preserved_4d_denominator_edges_for_3d_expression(representation)
                .into_iter()
                .map(usize::from)
                .collect(),
        })
    }

    pub fn automatic_numerator_energy_degree_bounds_for_3d_expression(
        &self,
        _representation: RepresentationMode,
    ) -> Result<Vec<(usize, usize)>> {
        let excluded_edges = self
            .iter_edges_of(&self.initial_state_cut)
            .map(|(_, edge_id, _)| edge_id)
            .chain(self.external_tree_4d_denominator_edges())
            .collect_vec();
        Ok(self.automatic_numerator_energy_degree_bounds_excluding(excluded_edges)?)
    }

    pub fn preserved_4d_denominator_edges_for_3d_expression(
        &self,
        _representation: RepresentationMode,
    ) -> Vec<EdgeIndex> {
        self.external_tree_4d_denominator_edges()
    }

    fn normalized_preserved_4d_denominator_edges(
        &self,
        options: &Generate3DExpressionOptions,
    ) -> Vec<EdgeIndex> {
        let mut preserved_edges =
            self.preserved_4d_denominator_edges_for_3d_expression(options.representation);
        preserved_edges.extend(
            options
                .preserve_internal_edges_as_four_d_denominators
                .iter()
                .copied()
                .map(EdgeIndex),
        );
        preserved_edges.sort_unstable();
        preserved_edges.dedup();
        preserved_edges
    }

    pub(crate) fn external_tree_4d_denominator_edges(&self) -> Vec<EdgeIndex> {
        self.iter_edges_of(
            &self
                .tree_edges
                .subtract(&self.initial_state_cut)
                .subtract(&self.external_filter::<SuBitGraph>()),
        )
        .filter_map(|(_, edge_id, _)| {
            (!self.loop_momentum_basis.edge_signatures[edge_id].is_loop_dependent())
                .then_some(edge_id)
        })
        .collect_vec()
    }

    pub(crate) fn denominator_only_cff_3d_expression_options(&self) -> Generate3DExpressionOptions {
        Generate3DExpressionOptions {
            representation: RepresentationMode::Cff,
            energy_degree_bounds: Vec::new(),
            numerator_sampling_scale: NumeratorSamplingScaleMode::None,
            include_cff_duplicate_signature_excess_sign: false,
            preserve_internal_edges_as_four_d_denominators: Vec::new(),
        }
    }

    pub(crate) fn convert_generated_expression_surfaces(
        &mut self,
        mut expression: three_dimensional_reps::ThreeDExpression<OrientationID>,
        representation: RepresentationMode,
        use_generated_cff_half_edges: bool,
        canonize_esurface: &Option<ShiftRewrite>,
        initial_state_cut_edges: &[EdgeIndex],
    ) -> Result<CFFExpression<OrientationID>> {
        let mut linear_surface_map = BTreeMap::<LinearSurfaceID, SurfaceMapEntry>::new();
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
            let mut remapped_variants = Vec::with_capacity(orientation.variants.len());
            for mut variant in std::mem::take(&mut orientation.variants) {
                let signed_denominators =
                    remap_denominator_tree_surface_ids(&variant.denominator, &linear_surface_map);
                for surface_id in &mut variant.numerator_surfaces {
                    let numerator_sign =
                        remap_generated_surface_id(surface_id, &linear_surface_map);
                    if numerator_sign < 0 {
                        variant.prefactor *= Atom::num(-1);
                    }
                }
                if representation == RepresentationMode::Cff && !use_generated_cff_half_edges {
                    // GammaLoop's CFF evaluator convention keeps the
                    // on-shell-energy factors as one global product
                    // 1/prod(-2E_i) for the pure CFF denominator sector. The
                    // bounded higher-energy CFF algorithm has variant-local
                    // half-edge factors, so those must remain attached to the
                    // variants. LTD half-edge factors are residue Jacobians and
                    // always stay on their generated variants.
                    variant.half_edges.clear();
                }
                for (denominator_sign, denominator) in signed_denominators {
                    let mut signed_variant = variant.clone();
                    signed_variant.denominator = denominator;
                    if denominator_sign < 0 {
                        signed_variant.prefactor *= Atom::num(-1);
                    }
                    remapped_variants.push(signed_variant);
                }
            }
            orientation.variants = remapped_variants;
        }

        Ok(CFFExpression {
            orientations: expression.orientations,
            surfaces: self.surface_cache.clone(),
            residual_denominators: expression.residual_denominators,
        })
    }

    fn intern_generated_linear_surface(
        &mut self,
        surface: &LinearSurface,
        canonize_esurface: &Option<ShiftRewrite>,
        initial_state_cut_edges: &[EdgeIndex],
    ) -> Result<SurfaceMapEntry> {
        if surface.numerator_only {
            let surface_id =
                Into::<LinearSurfaceID>::into(self.surface_cache.linear_surface_cache.len());
            self.surface_cache
                .linear_surface_cache
                .push(surface.clone());
            return Ok(SurfaceMapEntry {
                surface_id: HybridSurfaceID::Linear(surface_id),
                sign: 1,
            });
        }

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

        let (surface_id, sign) = match surface.kind {
            LinearSurfaceKind::Esurface => {
                let mut sign = 1;
                if positive_energies.is_empty() && !negative_energies.is_empty() {
                    sign = -1;
                    positive_energies = negative_energies;
                    negative_energies = Vec::new();
                    for (_, coeff) in &mut external_shift {
                        *coeff = -*coeff;
                    }
                }
                if !negative_energies.is_empty() {
                    return Err(eyre::eyre!(
                        "generalized 3D production cannot convert E-surface with mixed-sign internal terms {:?}",
                        surface
                    ));
                }
                let mut esurface = Esurface {
                    energies: positive_energies,
                    external_shift,
                    vertex_set: VertexSet::dummy(),
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

                (HybridSurfaceID::Esurface(esurface_id), sign)
            }
            LinearSurfaceKind::Hsurface => {
                let hsurface = Hsurface {
                    positive_energies,
                    negative_energies,
                    external_shift,
                    vertex_set: VertexSet::dummy(),
                };
                let hsurface_id = self
                    .surface_cache
                    .hsurface_cache
                    .position(|existing| existing == &hsurface)
                    .unwrap_or_else(|| {
                        self.surface_cache.hsurface_cache.push(hsurface);
                        Into::<HsurfaceID>::into(self.surface_cache.hsurface_cache.len() - 1)
                    });

                (HybridSurfaceID::Hsurface(hsurface_id), 1)
            }
        };
        Ok(SurfaceMapEntry { surface_id, sign })
    }
}

pub(crate) fn generated_cff_expression_uses_variant_half_edges(
    options: &Generate3DExpressionOptions,
) -> bool {
    options.representation == RepresentationMode::Cff
        && options
            .energy_degree_bounds
            .iter()
            .any(|(_, degree)| *degree > 1)
}

fn three_d_source_summary(parsed: &three_dimensional_reps::ParsedGraph) -> String {
    let internal_edges = parsed
        .internal_edges
        .iter()
        .map(|edge| {
            format!(
                "edge {}: {} -> {}, loops={:?}, externals={:?}",
                edge.edge_id,
                edge.tail,
                edge.head,
                edge.signature.loop_signature,
                edge.signature.external_signature
            )
        })
        .join("; ");
    format!(
        "3D source has {} loop names {:?}, {} internal edges [{}]",
        parsed.loop_names.len(),
        parsed.loop_names,
        parsed.internal_edges.len(),
        internal_edges
    )
}

fn source_contract_edges_for_3d_expression(
    contract_edges: &[EdgeIndex],
    initial_state_cut_edges: &[EdgeIndex],
    preserved_edges: &[EdgeIndex],
) -> Vec<EdgeIndex> {
    let mut non_contractible_edges = initial_state_cut_edges
        .iter()
        .copied()
        .collect::<HashSet<_>>();
    non_contractible_edges.extend(preserved_edges.iter().copied());
    contract_edges
        .iter()
        .copied()
        .filter(|edge_id| !non_contractible_edges.contains(edge_id))
        .collect_vec()
}

fn remap_generated_surface_id(
    surface_id: &mut HybridSurfaceID,
    linear_surface_map: &BTreeMap<LinearSurfaceID, SurfaceMapEntry>,
) -> i64 {
    if let HybridSurfaceID::Linear(linear_surface_id) = surface_id {
        let entry = *linear_surface_map
            .get(linear_surface_id)
            .expect("all generated linear surfaces should have been interned");
        *surface_id = entry.surface_id;
        entry.sign
    } else {
        1
    }
}

fn remap_denominator_tree_surface_ids(
    denominator: &Tree<HybridSurfaceID>,
    linear_surface_map: &BTreeMap<LinearSurfaceID, SurfaceMapEntry>,
) -> Vec<(i64, Tree<HybridSurfaceID>)> {
    let node_signs = denominator
        .iter_nodes()
        .map(|node| {
            let sign = if let HybridSurfaceID::Linear(linear_surface_id) = node.data {
                linear_surface_map
                    .get(&linear_surface_id)
                    .map(|entry| entry.sign)
                    .unwrap_or(1)
            } else {
                1
            };
            (node.node_id, sign)
        })
        .collect::<BTreeMap<_, _>>();

    let mut chains_by_sign = BTreeMap::<i64, Vec<Vec<HybridSurfaceID>>>::new();
    for leaf in denominator.get_bottom_layer() {
        let mut sign = 1;
        let mut chain = Vec::new();
        let mut current = Some(leaf);
        while let Some(node_id) = current {
            sign *= node_signs.get(&node_id).copied().unwrap_or(1);
            let mut surface_id = denominator.get_node(node_id).data;
            remap_generated_surface_id(&mut surface_id, linear_surface_map);
            if surface_id != HybridSurfaceID::Unit {
                chain.push(surface_id);
            }
            current = denominator.get_node(node_id).parent;
        }
        chain.reverse();
        chains_by_sign.entry(sign).or_default().push(chain);
    }

    if chains_by_sign.is_empty() {
        return vec![(1, Tree::from_root(HybridSurfaceID::Unit))];
    }

    chains_by_sign
        .into_iter()
        .map(|(sign, chains)| (sign, denominator_tree_from_chains(&chains)))
        .collect()
}

fn denominator_tree_from_chains(chains: &[Vec<HybridSurfaceID>]) -> Tree<HybridSurfaceID> {
    if chains.is_empty() || chains.iter().all(Vec::is_empty) {
        return Tree::from_root(HybridSurfaceID::Unit);
    }

    let mut tree = Tree::from_root(HybridSurfaceID::Unit);
    for chain in chains {
        if chain.is_empty() {
            insert_terminal_unit_if_missing(&mut tree, NodeId::root());
            continue;
        }
        let mut parent = NodeId::root();
        for surface_id in chain {
            let existing_child = tree
                .get_node(parent)
                .children
                .iter()
                .copied()
                .find(|child| tree.get_node(*child).data == *surface_id);
            if let Some(child) = existing_child {
                parent = child;
                continue;
            }
            let child = NodeId(tree.get_num_nodes());
            tree.insert_node(parent, *surface_id);
            parent = child;
        }
        insert_terminal_unit_if_missing(&mut tree, parent);
    }
    tree
}

fn insert_terminal_unit_if_missing(tree: &mut Tree<HybridSurfaceID>, parent: NodeId) {
    let has_terminal_unit = tree
        .get_node(parent)
        .children
        .iter()
        .any(|child| tree.get_node(*child).data == HybridSurfaceID::Unit);
    if !has_terminal_unit {
        tree.insert_node(parent, HybridSurfaceID::Unit);
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

fn numerator_sampling_scale_mode(
    setting: UniformNumeratorSamplingScale,
) -> NumeratorSamplingScaleMode {
    match setting {
        UniformNumeratorSamplingScale::None => NumeratorSamplingScaleMode::None,
        UniformNumeratorSamplingScale::BeyondQuadratic => {
            NumeratorSamplingScaleMode::BeyondQuadratic
        }
        UniformNumeratorSamplingScale::All => NumeratorSamplingScaleMode::All,
    }
}

fn representation_mode(representation: ThreeDRepresentation) -> RepresentationMode {
    match representation {
        ThreeDRepresentation::Cff => RepresentationMode::Cff,
        ThreeDRepresentation::Ltd => RepresentationMode::Ltd,
    }
}
