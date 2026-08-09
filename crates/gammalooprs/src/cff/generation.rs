use std::collections::BTreeMap;

use crate::{
    cff::{
        VertexSet,
        hsurface::{Hsurface, HsurfaceID},
        surface::{
            HybridSurfaceID, LinearEnergyExpr, LinearSurface, LinearSurfaceID, LinearSurfaceKind,
        },
    },
    graph::{FeynmanGraph, Graph, GraphThreeDSource},
    settings::global::{GenerationSettings, UniformNumeratorSamplingScale},
};
use ahash::HashSet;
use color_eyre::Result;
use itertools::Itertools;
use linnet::half_edge::involution::EdgeIndex;
use linnet::num_traits::SignOrZero;
use symbolica::atom::{Atom, AtomCore};
use three_dimensional_reps::{
    CffEnergyFactorOwnership, Generate3DExpressionOptions, GeneratedThreeDExpression,
    NumeratorSamplingScaleMode, RepresentationMode, ThreeDGraphSource,
    tree::{NodeId, Tree},
};

use tracing::debug;

use super::{
    esurface::{Esurface, EsurfaceID, ExternalShift},
    expression::CFFExpression,
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

struct RemappedDenominatorTree {
    prefactor_sign: i64,
    surface_signs: BTreeMap<HybridSurfaceID, i64>,
    denominator: Tree<HybridSurfaceID>,
}

impl Graph {
    pub(crate) fn generate_3d_expression_for_4d_term(
        &self,
        source: &GraphThreeDSource<'_>,
        options: &Generate3DExpressionOptions,
        analysis_numerator: &Atom,
    ) -> Result<GeneratedThreeDExpression> {
        let mut source_options = options.clone();
        let initial_state_cut_edges = self
            .iter_edges_of(&self.initial_state_cut)
            .map(|(_, edge_id, _)| edge_id)
            .collect::<HashSet<_>>();
        let bridge_edges = self
            .iter_edges_of(&self.tree_edges)
            .map(|(_, edge_id, _)| edge_id)
            .filter(|edge_id| !initial_state_cut_edges.contains(edge_id))
            .collect::<HashSet<_>>();
        let contracted_source = source.contract_subgraph();
        for (_, edge_id, _) in self.iter_edges_of(&contracted_source) {
            if !bridge_edges.contains(&edge_id) {
                continue;
            }
            let Some(coordinates) = source.reconstructible_outer_loop_coordinates(edge_id) else {
                return Err(eyre::eyre!(
                    "contracted bridge edge {} retains an inner-loop energy coordinate and cannot be projected independently of the exact CFF source",
                    usize::from(edge_id),
                ));
            };
            if coordinates.iter().any(|coordinate| *coordinate != 0) {
                return Err(eyre::eyre!(
                    "contracted bridge edge {} has nonzero outer-loop coordinates {:?}; a CFF-external tree edge must have a purely external affine energy",
                    usize::from(edge_id),
                    coordinates,
                ));
            }
        }
        let energy_degree_bounds = self
            .automatic_numerator_energy_degree_bounds_in_atom_excluding_with_min_degree(
                analysis_numerator,
                initial_state_cut_edges.iter().chain(&bridge_edges).copied(),
                1,
            )
            .map_err(|error| {
                eyre::eyre!(
                    "could not analyze numerator in the contracted-source loop basis: {error}"
                )
            })?;
        if !energy_degree_bounds.is_empty() {
            return Err(eyre::eyre!(
                "exact 4D CFF generation for graph `{}` cannot soundly map nonconstant nonbridge numerator EMR energy-degree bounds {:?} to occurrence-local energies: an affine numerator-only EMR-to-occurrence contract is not implemented",
                self.name,
                energy_degree_bounds,
            ));
        }
        // The exact source has occurrence-local denominator IDs. Until it can
        // carry physical EMR numerator variables through an independent
        // affine map, only the explicit energy-constant nonbridge numerator
        // class is sound inside CFF. Bridge energies remain in the complete
        // factorized numerator and are projected independently as external
        // affine tree factors. This overrides production-wide bounds with the
        // class of the numerator actually owned by this 4D term.
        source_options.energy_degree_bounds = Some(energy_degree_bounds);
        let parsed = source.to_three_d_parsed_graph()?;
        let generated = three_dimensional_reps::generate_3d_expression(source, &source_options);
        generated.map_err(|error| {
            eyre::eyre!(
                "generalized CFF expression generation failed for exact 4D source in graph `{}` with source-edge numerator energy-degree bounds {:?}: {error}\n{}",
                self.name,
                source_options.energy_degree_bounds,
                three_d_source_summary(&parsed),
            )
        })
    }

    pub(crate) fn generate_3d_expression_for_integrand(
        &mut self,
        contract_edges: &[EdgeIndex],
        canonize_esurface: &Option<ShiftRewrite>,
        options: &Generate3DExpressionOptions,
        analysis_numerator: Option<&Atom>,
    ) -> Result<GeneratedThreeDExpression<Esurface, Hsurface>> {
        let initial_state_cut_edges = self
            .iter_edges_of(&self.initial_state_cut)
            .map(|(_, edge_id, _)| edge_id)
            .collect_vec();
        let source_contract_edges =
            source_contract_edges_for_3d_expression(contract_edges, &initial_state_cut_edges);
        if !source_contract_edges.is_empty() {
            debug!(
                "contracting {} internal edges before generalized 3D expression generation",
                source_contract_edges.len()
            );
        }
        let source = GraphThreeDSource::new(self, &source_contract_edges)?;
        let bridge_edges = self
            .iter_edges_of(&self.tree_edges)
            .map(|(_, edge_id, _)| edge_id)
            .filter(|edge_id| !initial_state_cut_edges.contains(edge_id))
            .collect::<HashSet<_>>();
        for edge_id in source_contract_edges
            .iter()
            .filter(|edge_id| bridge_edges.contains(edge_id))
        {
            let Some(coordinates) = source.reconstructible_outer_loop_coordinates(*edge_id) else {
                return Err(eyre::eyre!(
                    "contracted bridge edge {} retains an inner-loop energy coordinate and cannot be projected independently of the CFF source",
                    usize::from(*edge_id),
                ));
            };
            if coordinates.iter().any(|coordinate| *coordinate != 0) {
                return Err(eyre::eyre!(
                    "contracted bridge edge {} has nonzero outer-loop coordinates {:?}; a CFF-external tree edge must have a purely external affine energy",
                    usize::from(*edge_id),
                    coordinates,
                ));
            }
        }
        let mut source_options = options.clone();
        if let Some(numerator) = analysis_numerator {
            source_options.energy_degree_bounds = Some(
                self.automatic_numerator_energy_degree_bounds_in_atom_excluding_with_min_degree(
                    numerator,
                    initial_state_cut_edges.iter().chain(&bridge_edges).copied(),
                    1,
                )?,
            );
            debug!(
                graph = %self.name,
                bounds = ?source_options.energy_degree_bounds,
                "using source-edge numerator energy-degree bounds"
            );
        }
        if let Some(bounds) = &mut source_options.energy_degree_bounds {
            // Cut aliases and tree denominators are external to the CFF graph,
            // but their energies remain in the complete numerator. Removing
            // their EMR bounds here therefore changes only CFF capacity; in
            // particular, an explicitly bounded numerator class remains
            // `Some([])` when every bound belongs to this external sector.
            bounds.retain(|(edge_id, _)| {
                let edge_id = EdgeIndex(*edge_id);
                !initial_state_cut_edges.contains(&edge_id) && !bridge_edges.contains(&edge_id)
            });
        }
        for &(edge_id, degree) in source_options
            .energy_degree_bounds
            .as_deref()
            .unwrap_or(&[])
        {
            let edge_id = EdgeIndex(edge_id);
            if !source_contract_edges.contains(&edge_id) {
                continue;
            }
            if source
                .reconstructible_outer_loop_coordinates(edge_id)
                .is_some()
            {
                return Err(eyre::eyre!(
                    "numerator energy-degree bound {degree} for shrunken EMR edge {} cannot yet be supplied to generalized CFF generation: its exact outer affine energy is retained for numerator evaluation, but the standalone generator has no numerator-only/shrunken-edge bound channel",
                    usize::from(edge_id),
                ));
            }
            return Err(eyre::eyre!(
                "numerator energy-degree bound {degree} for shrunken EMR edge {} retains contracted inner-loop dependence and cannot be represented in the outer CFF source",
                usize::from(edge_id),
            ));
        }
        let mut generated = {
            let result = three_dimensional_reps::generate_3d_expression(&source, &source_options);
            result.map_err(|error| {
                    let source_summary = source
                        .to_three_d_parsed_graph()
                        .map(|parsed| three_d_source_summary(&parsed))
                        .unwrap_or_else(|source_error| {
                            format!("failed to rebuild 3D source summary: {source_error}")
                    });
                    eyre::eyre!(
                        "generalized CFF expression generation failed for graph `{}` with source-edge numerator energy-degree bounds {:?}: {error}\n{source_summary}",
                        self.name,
                        source_options.energy_degree_bounds,
                    )
                })
        }?;

        // Generic edge-index remapping leaves omitted parent edges at zero. A
        // contracted edge can still have an exact outer-energy map when its
        // source coordinates contain no inner-loop component; restore only
        // that source-owned information and leave genuinely inner-dependent
        // contracted edges unavailable.
        for edge_id in &source_contract_edges {
            let Some(coordinates) = source.reconstructible_outer_loop_coordinates(*edge_id) else {
                continue;
            };
            let signature = &self.loop_momentum_basis.edge_signatures[*edge_id];
            for orientation in generated.expression.orientations.iter_mut() {
                if coordinates.len() != orientation.loop_energy_map.len() {
                    return Err(eyre::eyre!(
                        "contracted edge {} has {} outer coordinates for {} generated loop-energy maps",
                        usize::from(*edge_id),
                        coordinates.len(),
                        orientation.loop_energy_map.len(),
                    ));
                }
                let mut edge_energy = coordinates.iter().zip(&orientation.loop_energy_map).fold(
                    LinearEnergyExpr::zero(),
                    |sum, (coefficient, loop_energy)| {
                        sum + loop_energy.clone().scale_rational(coefficient.clone())
                    },
                );
                for (external_edge, sign) in self
                    .loop_momentum_basis
                    .ext_edges
                    .iter()
                    .zip(&signature.external)
                {
                    let coefficient = match sign {
                        SignOrZero::Zero => 0,
                        SignOrZero::Plus => 1,
                        SignOrZero::Minus => -1,
                    };
                    edge_energy =
                        edge_energy + LinearEnergyExpr::external(*external_edge, coefficient);
                }
                let edge_energy_slot = orientation
                    .edge_energy_map
                    .get_mut(usize::from(*edge_id))
                    .ok_or_else(|| {
                    eyre::eyre!(
                        "contracted edge {} is outside the generated parent energy map",
                        usize::from(*edge_id),
                    )
                })?;
                *edge_energy_slot = edge_energy;
            }
        }

        self.convert_generated_expression_surfaces(
            generated,
            canonize_esurface,
            &initial_state_cut_edges,
        )
    }

    pub(crate) fn production_cff_3d_expression_options(
        &self,
        settings: &GenerationSettings,
    ) -> Result<Generate3DExpressionOptions> {
        self.cff_3d_expression_options(numerator_sampling_scale_mode(
            settings.uniform_numerator_sampling_scale,
        ))
    }

    pub fn cff_3d_expression_options(
        &self,
        numerator_sampling_scale: NumeratorSamplingScaleMode,
    ) -> Result<Generate3DExpressionOptions> {
        let numerator = self.production_numerator_atom_for_full_3d_expression();
        let cff_external_edges = self
            .iter_edges_of(&self.initial_state_cut)
            .chain(self.iter_edges_of(&self.tree_edges))
            .map(|(_, edge_id, _)| edge_id);
        // Analyze the factorized numerator compositionally while keeping every
        // active energy degree attached to its EMR source edge. Initial-cut
        // aliases and structural bridges stay in their separately projected
        // external sector. The standalone CFF generator remaps active source
        // IDs only after extracting its local graph, so a loop-momentum basis
        // never defines numerator ownership.
        let energy_degree_bounds = self
            .automatic_numerator_energy_degree_bounds_in_atom_excluding_with_min_degree(
                &numerator,
                cff_external_edges,
                1,
            )?;
        debug!(
            graph = %self.name,
            bounds = ?energy_degree_bounds,
            "using production CFF source-edge numerator energy-degree bounds"
        );
        Ok(Generate3DExpressionOptions {
            representation: RepresentationMode::Cff,
            energy_degree_bounds: Some(energy_degree_bounds),
            numerator_sampling_scale,
            preserve_internal_edges_as_four_d_denominators: Vec::new(),
        })
    }

    pub(crate) fn denominator_only_cff_3d_expression_options(&self) -> Generate3DExpressionOptions {
        Generate3DExpressionOptions {
            representation: RepresentationMode::Cff,
            energy_degree_bounds: Some(Vec::new()),
            numerator_sampling_scale: NumeratorSamplingScaleMode::None,
            preserve_internal_edges_as_four_d_denominators: Vec::new(),
        }
    }

    pub(crate) fn convert_generated_expression_surfaces(
        &mut self,
        generated: GeneratedThreeDExpression,
        canonize_esurface: &Option<ShiftRewrite>,
        initial_state_cut_edges: &[EdgeIndex],
    ) -> Result<GeneratedThreeDExpression<Esurface, Hsurface>> {
        self.convert_generated_expression_surfaces_impl(
            generated,
            canonize_esurface,
            initial_state_cut_edges,
            None,
        )
    }

    pub(crate) fn convert_4d_expression_surfaces(
        &mut self,
        generated: GeneratedThreeDExpression,
        physical_surfaces: &[Option<LinearSurface>],
    ) -> Result<GeneratedThreeDExpression<Esurface, Hsurface>> {
        let canonize_esurface = self.get_esurface_canonization(&self.loop_momentum_basis);
        let initial_state_cut_edges = self
            .iter_edges_of(&self.initial_state_cut)
            .map(|(_, edge_id, _)| edge_id)
            .collect_vec();
        self.convert_generated_expression_surfaces_impl(
            generated,
            &canonize_esurface,
            &initial_state_cut_edges,
            Some(physical_surfaces),
        )
    }

    fn convert_generated_expression_surfaces_impl(
        &mut self,
        generated: GeneratedThreeDExpression,
        canonize_esurface: &Option<ShiftRewrite>,
        initial_state_cut_edges: &[EdgeIndex],
        physical_surfaces: Option<&[Option<LinearSurface>]>,
    ) -> Result<GeneratedThreeDExpression<Esurface, Hsurface>> {
        let GeneratedThreeDExpression {
            mut expression,
            energy_factor_ownership,
        } = generated;
        if !expression.residual_denominators.is_empty() {
            return Err(eyre::eyre!(
                "GammaLoop production CFF conversion received residual four-dimensional denominator edges [{}]; production callers must contract tree edges and attach their projected denominators separately",
                expression
                    .residual_denominators
                    .iter()
                    .map(|denominator| denominator.edge_id.0)
                    .join(", "),
            ));
        }
        let mut linear_surface_map = BTreeMap::<LinearSurfaceID, SurfaceMapEntry>::new();
        let mut retained_linear_surfaces = Vec::new();
        for (linear_surface_id, surface) in
            expression.surfaces.linear_surface_cache.iter_enumerated()
        {
            let physical = physical_surfaces
                .and_then(|surfaces| surfaces.get(usize::from(linear_surface_id)))
                .and_then(Option::as_ref);
            if physical_surfaces.is_some() && physical.is_none() {
                retained_linear_surfaces.push((linear_surface_id, surface.clone()));
            } else {
                let converted = self.intern_generated_linear_surface(
                    physical.unwrap_or(surface),
                    canonize_esurface,
                    initial_state_cut_edges,
                )?;
                linear_surface_map.insert(linear_surface_id, converted);
            }
        }
        let mut surface_cache = self.surface_cache.clone();
        for (source_id, surface) in retained_linear_surfaces {
            let target_id = LinearSurfaceID(surface_cache.linear_surface_cache.len());
            surface_cache.linear_surface_cache.push(surface);
            linear_surface_map.insert(
                source_id,
                SurfaceMapEntry {
                    surface_id: HybridSurfaceID::Linear(target_id),
                    sign: 1,
                },
            );
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
                if energy_factor_ownership == CffEnergyFactorOwnership::GlobalSourceProduct {
                    // GammaLoop's CFF evaluator convention keeps the
                    // on-shell-energy factors as one global product
                    // 1/prod(-2E_i) for the ordinary CFF denominator sector. The
                    // generalized numerator and causal powered-pole algorithms
                    // instead attach those factors to individual variants.
                    variant.half_edges.clear();
                }
                for remapped_denominator in signed_denominators {
                    let mut signed_variant = variant.clone();
                    signed_variant.denominator = remapped_denominator.denominator;
                    signed_variant.denominator_surface_signs = remapped_denominator.surface_signs;
                    if remapped_denominator.prefactor_sign < 0 {
                        signed_variant.prefactor *= Atom::num(-1);
                    }
                    remapped_variants.push(signed_variant);
                }
            }
            orientation.variants = remapped_variants;
        }

        Ok(GeneratedThreeDExpression {
            expression: CFFExpression {
                orientations: expression.orientations,
                surfaces: surface_cache,
                residual_denominators: Vec::new(),
            },
            energy_factor_ownership,
        })
    }

    fn intern_generated_linear_surface(
        &mut self,
        surface: &LinearSurface,
        canonize_esurface: &Option<ShiftRewrite>,
        initial_state_cut_edges: &[EdgeIndex],
    ) -> Result<SurfaceMapEntry> {
        if surface.numerator_only
            && (!surface.expression.uniform_scale_coeff.is_zero()
                || !surface.expression.constant.is_zero())
        {
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
) -> Vec<EdgeIndex> {
    let non_contractible_edges = initial_state_cut_edges
        .iter()
        .copied()
        .collect::<HashSet<_>>();
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
) -> Vec<RemappedDenominatorTree> {
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

    let mut chains_by_signature =
        BTreeMap::<(i64, Vec<(HybridSurfaceID, i64)>), Vec<Vec<HybridSurfaceID>>>::new();
    for leaf in denominator.get_bottom_layer() {
        let mut sign = 1;
        let mut surface_signs = BTreeMap::<HybridSurfaceID, i64>::new();
        let mut chain = Vec::new();
        let mut current = Some(leaf);
        while let Some(node_id) = current {
            let mut surface_id = denominator.get_node(node_id).data;
            let node_sign = node_signs.get(&node_id).copied().unwrap_or(1);
            sign *= node_sign;
            remap_generated_surface_id(&mut surface_id, linear_surface_map);
            if surface_id != HybridSurfaceID::Unit {
                if node_sign < 0 {
                    let entry = surface_signs.entry(surface_id).or_insert(1);
                    *entry *= node_sign;
                    if *entry > 0 {
                        surface_signs.remove(&surface_id);
                    }
                }
                chain.push(surface_id);
            }
            current = denominator.get_node(node_id).parent;
        }
        chain.reverse();
        chains_by_signature
            .entry((sign, surface_signs.into_iter().collect()))
            .or_default()
            .push(chain);
    }

    if chains_by_signature.is_empty() {
        return vec![RemappedDenominatorTree {
            prefactor_sign: 1,
            surface_signs: BTreeMap::new(),
            denominator: Tree::from_root(HybridSurfaceID::Unit),
        }];
    }

    chains_by_signature
        .into_iter()
        .map(
            |((prefactor_sign, surface_signs), chains)| RemappedDenominatorTree {
                prefactor_sign,
                surface_signs: surface_signs.into_iter().collect(),
                denominator: denominator_tree_from_chains(&chains),
            },
        )
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        dot,
        graph::{
            FeynmanGraph,
            cuts::{CutSet, LuCutSelection},
            parse::from_dot::IntoGraph,
        },
        initialisation::test_initialise,
        settings::global::{GenerationSettings, OrientationPattern},
        utils::GS,
    };

    #[test]
    fn contracted_raised_generation_rejects_absent_emr_bound_without_aliasing() -> Result<()> {
        test_initialise()?;
        let graph: Graph = dot!(
            digraph factorized_source_bounds {
                edge [pdg=1000 num=1 mass=0]
                node [num=1]
                ext_in [style=invis]
                ext_out [style=invis]
                ext_in -> A [id=4]
                C -> ext_out [id=5]
                A -> B [id=0 num="(Q(0,spenso::cind(0))+1)*(Q(1,spenso::cind(0))+2)"]
                B -> C [id=1 num=1]
                C -> D [id=2]
                D -> A [id=3]
            },
            "scalars"
        )?;
        let options = graph.production_cff_3d_expression_options(&GenerationSettings::default())?;
        assert_eq!(options.energy_degree_bounds, Some(vec![(0, 1), (1, 1)]));

        let numerator = graph.production_numerator_atom_for_full_3d_expression();
        let mut contracted_graph = graph.clone();
        let error = contracted_graph
            .generate_3d_expression_for_integrand(
                &[EdgeIndex(0)],
                &None,
                &options,
                Some(&numerator),
            )
            .expect_err("a bound owned by a shrunken EMR edge must not be silently reassigned");
        let error = format!("{error:#}");
        assert!(
            error.contains("shrunken EMR edge 0")
                && error.contains("no numerator-only/shrunken-edge bound channel"),
            "unexpected error: {error}"
        );
        Ok(())
    }

    #[test]
    fn raised_lu_cff_uses_causal_powered_pole_maps() -> Result<()> {
        test_initialise()?;
        let mut graph: Graph = dot!(digraph raised_lu {
            edge [num=1 mass=1]
            node [num=1]
            incoming [style=invis]
            outgoing [style=invis]

            incoming -> v1 [id=0]
            v1 -> v2 [id=1 lmb_id=0]
            v2 -> v3 [id=2]
            v1 -> v3 [id=3]
            v3 -> outgoing [id=4]
        })?;
        let options = graph.denominator_only_cff_3d_expression_options();
        let numerator = GS.emr_mom(EdgeIndex(1), GS.cind(0)).pow(2);
        let canonization = graph.get_esurface_canonization(&graph.loop_momentum_basis);
        let generated = graph.generate_3d_expression_for_integrand(
            &[],
            &canonization,
            &options,
            Some(&numerator),
        )?;
        assert_eq!(
            generated.energy_factor_ownership,
            CffEnergyFactorOwnership::VariantLocal
        );
        let production = generated.expression;
        let lu_cut = graph
            .determine_raised_esurfaces_from_expression(&production)
            .raised_groups
            .into_iter()
            .find(|group| {
                group.max_occurence > 1
                    && group.esurface_ids.iter().any(|esurface_id| {
                        !production.surfaces.esurface_cache[*esurface_id]
                            .external_shift
                            .is_empty()
                    })
            })
            .expect("raised production CFF should contain a physical repeated LU surface");
        let mut cutset = CutSet::empty(graph.n_hedges());
        cutset.residue_selector.lu = Some(LuCutSelection {
            raised_group: lu_cut.clone(),
            cut_edge_alternatives: lu_cut
                .esurface_ids
                .iter()
                .map(|esurface_id| {
                    production.surfaces.esurface_cache[*esurface_id]
                        .energies
                        .clone()
                })
                .collect(),
        });
        let direct = graph.cff(
            &cutset.union,
            &cutset,
            &OrientationPattern::default(),
            &options,
            Some(&numerator),
        )?;

        assert_eq!(direct.terms.len(), lu_cut.max_occurence);
        assert!(direct.terms.keys().all(|index| {
            index.lu_cut_order.is_some()
                && index.left_threshold_order.is_none()
                && index.right_threshold_order.is_none()
        }));
        assert!(direct.terms.iter().any(|(index, term)| {
            index.lu_cut_order == Some(2) && !term.orientations.is_empty()
        }));
        let selected_variants = direct
            .terms
            .values()
            .flat_map(|term| &term.orientations)
            .flat_map(|orientation| &orientation.orientation.variants)
            .collect_vec();
        assert!(!selected_variants.is_empty());
        assert!(
            selected_variants
                .iter()
                .all(|variant| !variant.half_edges.is_empty()),
            "selected causal powered-pole variants must retain their local on-shell-energy factors"
        );
        for local in direct.terms.values().flat_map(|term| &term.orientations) {
            assert!(
                production.orientations.iter().any(|full| {
                    full.data.orientation == local.orientation.data.orientation
                        && full.loop_energy_map == local.orientation.loop_energy_map
                        && full.edge_energy_map == local.orientation.edge_energy_map
                }),
                "local LU map {:?} with loop energies {:?} and edge energies {:?} is absent from the raised production carrier",
                local.orientation.data.orientation,
                local.orientation.loop_energy_map,
                local.orientation.edge_energy_map,
            );
        }
        Ok(())
    }
}
