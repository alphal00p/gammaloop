use std::collections::BTreeMap;

use crate::{
    cff::{
        VertexSet,
        hsurface::{Hsurface, HsurfaceID},
        surface::{HybridSurfaceID, LinearSurface, LinearSurfaceID, LinearSurfaceKind},
    },
    graph::{Graph, GraphThreeDSource},
    settings::global::{GenerationSettings, UniformNumeratorSamplingScale},
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
                "contracting {} internal edges before generalized 3D CFF generation",
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
                        "generalized 3D CFF generation failed for graph `{}`: {error}\n{source_summary}",
                        self.name
                    )
                },
            )
        }?;

        self.convert_generated_expression_surfaces(
            expression,
            canonize_esurface,
            &initial_state_cut_edges,
        )
    }

    pub(crate) fn production_cff_3d_expression_options(
        &self,
        settings: &GenerationSettings,
    ) -> Result<Generate3DExpressionOptions> {
        self.three_d_expression_options(
            RepresentationMode::Cff,
            numerator_sampling_scale_mode(settings.uniform_numerator_sampling_scale),
        )
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
        representation: RepresentationMode,
    ) -> Vec<EdgeIndex> {
        if representation == RepresentationMode::Cff {
            self.external_tree_4d_denominator_edges()
        } else {
            Vec::new()
        }
    }

    fn normalized_preserved_4d_denominator_edges(
        &self,
        options: &Generate3DExpressionOptions,
    ) -> Vec<EdgeIndex> {
        let mut preserved_edges =
            self.preserved_4d_denominator_edges_for_3d_expression(options.representation);
        if options.representation == RepresentationMode::Cff {
            preserved_edges.extend(
                options
                    .preserve_internal_edges_as_four_d_denominators
                    .iter()
                    .copied()
                    .map(EdgeIndex),
            );
        }
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
            preserve_internal_edges_as_four_d_denominators: Vec::new(),
        }
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
            residual_denominators: expression.residual_denominators,
        })
    }

    fn intern_generated_linear_surface(
        &mut self,
        surface: &LinearSurface,
        canonize_esurface: &Option<ShiftRewrite>,
        initial_state_cut_edges: &[EdgeIndex],
    ) -> Result<HybridSurfaceID> {
        if surface.numerator_only {
            let surface_id =
                Into::<LinearSurfaceID>::into(self.surface_cache.linear_surface_cache.len());
            self.surface_cache
                .linear_surface_cache
                .push(surface.clone());
            return Ok(HybridSurfaceID::Linear(surface_id));
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

                Ok(HybridSurfaceID::Esurface(esurface_id))
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

                Ok(HybridSurfaceID::Hsurface(hsurface_id))
            }
        }
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
