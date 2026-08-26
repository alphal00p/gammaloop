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
    numerator::energy_degree::EnergyPowerAssignmentPlan,
    settings::global::{GenerationSettings, UniformNumeratorSamplingScale},
};
use ahash::HashSet;
use color_eyre::Result;
use itertools::Itertools;
use linnet::half_edge::involution::EdgeIndex;
use linnet::num_traits::SignOrZero;
use symbolica::atom::{Atom, AtomCore};
use three_dimensional_reps::{
    CffEnergyFactorOwnership, CffGenerationContext, EnergyEdgeIndexMap,
    Generate3DExpressionOptions, GeneratedThreeDExpression, NumeratorSamplingScaleMode,
    ParsedGraph, RepresentationMode, ThreeDGraphSource,
    tree::{NodeId, Tree},
};

use tracing::debug;

use super::{
    CffEnergyBoundSourceKind, CffEnergyDegreeBoundReport,
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
    prefactor: i64,
    surface_signs: BTreeMap<HybridSurfaceID, i64>,
    denominator: Tree<HybridSurfaceID>,
}

type ExactCffGenerationKey = (ParsedGraph, EnergyEdgeIndexMap, Generate3DExpressionOptions);

#[derive(Default)]
pub(crate) struct ExactCffGenerationCache {
    // Keep each occurrence's requested capacity, even when source provenance
    // identifies equal physical energies. Reuse requires canonical topology
    // and identical per-edge bounds; independent terms must not combine into
    // a larger Cartesian capacity or redistribute their numerator ownership.
    entries: BTreeMap<ExactCffGenerationKey, GeneratedThreeDExpression>,
    // Losing trials keep only their count, so repeated terms can compare the
    // same complete source capacity without retaining or rebuilding its trees.
    map_counts: BTreeMap<ExactCffGenerationKey, usize>,
}

impl ExactCffGenerationCache {
    pub(crate) fn len(&self) -> usize {
        self.entries.len()
    }

    fn generation_key(
        parsed: &ParsedGraph,
        energy_edges: &EnergyEdgeIndexMap,
        options: &Generate3DExpressionOptions,
    ) -> ExactCffGenerationKey {
        let mut topology = parsed.clone();
        topology.node_name_to_internal = topology
            .internal_edges
            .iter()
            .flat_map(|edge| [edge.tail, edge.head])
            .chain(
                topology
                    .external_edges
                    .iter()
                    .flat_map(|edge| edge.source.into_iter().chain(edge.destination)),
            )
            .collect::<std::collections::BTreeSet<_>>()
            .into_iter()
            .map(|node| (format!("__gammaloop_exact_node_{node}"), node))
            .collect();
        (topology, energy_edges.clone(), options.clone())
    }
}

struct ExactCffGenerationPreparation {
    parsed: ParsedGraph,
    energy_edges: EnergyEdgeIndexMap,
    source_options: Generate3DExpressionOptions,
    exact_source_energy_mapper: crate::graph::three_d_source::ExactSourceEnergyMapper,
    energy_assignment_plans: Vec<EnergyPowerAssignmentPlan>,
    physical_energy_degree_bounds: Vec<(usize, usize)>,
}

impl Graph {
    fn prepare_3d_expression_for_4d_term(
        &self,
        source: &GraphThreeDSource<'_>,
        options: &Generate3DExpressionOptions,
        analysis_numerator: &Atom,
    ) -> Result<ExactCffGenerationPreparation> {
        let source_options = options.clone();
        let initial_state_cut_edges = self
            .iter_edges_of(&self.initial_state_cut)
            .map(|(_, edge_id, _)| edge_id)
            .collect::<HashSet<_>>();
        let bridge_edges = self
            .iter_edges_of(&self.tree_edges)
            .map(|(_, edge_id, _)| edge_id)
            .filter(|edge_id| !initial_state_cut_edges.contains(edge_id))
            .collect::<HashSet<_>>();
        let excluded_numerator_edges = initial_state_cut_edges
            .iter()
            .chain(&bridge_edges)
            .copied()
            .chain(source.factorized_external_emr_edges())
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
        let physical_energy_degree_bounds = self
            .automatic_numerator_energy_degree_bounds_in_atoms_excluding_with_min_degree(
                [analysis_numerator],
                excluded_numerator_edges.iter().copied(),
                1,
            )
            .map_err(|error| {
                eyre::eyre!("could not analyze numerator in physical EMR energy variables: {error}")
            })?;
        // Parse first so a malformed exact rational source returns its
        // structural error instead of being hidden behind the mapper's
        // optional convenience API.
        let parsed = source.to_three_d_parsed_graph()?;
        let energy_edges = source
            .energy_edge_index_map(&parsed)
            .expect("exact 4D source has an occurrence-local energy map");
        let exact_source_energy_mapper = source
            .exact_source_energy_mapper()
            .expect("exact 4D source has an owned parent-energy mapper");
        // Exact sources have occurrence-local denominator IDs. Original
        // factors retain their base occurrence; only denominator-derived hard
        // factors may use serial copies of that owner. Analyze all physical
        // active edges first so unused, unrelated candidate groups cannot
        // reject a constant numerator.
        let candidates = exact_source_energy_mapper
            .equivalent_energy_candidates(
                physical_energy_degree_bounds
                    .iter()
                    .map(|(edge, _)| EdgeIndex(*edge)),
            )
            .map_err(|error| {
                eyre::eyre!(
                    "could not certify exact 4D CFF numerator energies for graph `{}`: {error}",
                    self.name,
                )
            })?;
        // Each immutable factor-local plan owns both its exact bounds and the
        // later numerator substitutions. This keeps the numerator factorized
        // and prevents generation from understating the expression actually
        // sampled in a residue or contact sector. Rank proposes a bounded set
        // of plans; the real source map count chooses between them below.
        let energy_assignment_plans = self
            .plan_numerator_energy_assignment_proposals_in_atom_excluding(
                analysis_numerator,
                excluded_numerator_edges.iter().copied(),
                &candidates,
            )
            .map_err(|error| {
                eyre::eyre!(
                    "could not plan exact 4D CFF numerator energies for graph `{}`: {error}",
                    self.name,
                )
            })?;
        debug!(
            graph = %self.name,
            physical_energy_degree_bounds = ?physical_energy_degree_bounds,
            equivalent_energy_candidates = ?candidates,
            candidate_bounds = ?energy_assignment_plans.iter().map(|plan| plan.energy_degree_bounds()).collect::<Vec<_>>(),
            "planned factorized exact-CFF numerator energy assignment proposals"
        );
        Ok(ExactCffGenerationPreparation {
            parsed,
            energy_edges,
            source_options,
            exact_source_energy_mapper,
            energy_assignment_plans,
            physical_energy_degree_bounds,
        })
    }

    pub(crate) fn generate_3d_expression_for_4d_term(
        &self,
        source: &GraphThreeDSource<'_>,
        options: &Generate3DExpressionOptions,
        analysis_numerator: &Atom,
        mut cache: Option<&mut ExactCffGenerationCache>,
    ) -> Result<(
        GeneratedThreeDExpression,
        crate::graph::three_d_source::ExactSourceEnergyMapper,
        EnergyPowerAssignmentPlan,
        CffEnergyDegreeBoundReport,
    )> {
        let ExactCffGenerationPreparation {
            parsed,
            energy_edges,
            mut source_options,
            exact_source_energy_mapper,
            energy_assignment_plans,
            physical_energy_degree_bounds,
        } = self.prepare_3d_expression_for_4d_term(source, options, analysis_numerator)?;
        let generate = |source_options: &Generate3DExpressionOptions| {
            crate::debug_tags!(#generation, #uv, #local, #four_d, #cff, #profile;
                graph = %self.name,
                term_local_bounds = ?source_options.energy_degree_bounds,
                file.parsed_source = ?parsed,
                "Generating exact CFF at its term-local capacity"
            );
            three_dimensional_reps::generate_3d_expression(source, source_options).map_err(
                |error| {
                    eyre::eyre!(
                        "generalized CFF expression generation failed for exact 4D source in graph `{}` with physical EMR bounds {:?} and term-local exact-occurrence bounds {:?}: {error}\n{}",
                        self.name,
                        physical_energy_degree_bounds,
                        source_options.energy_degree_bounds,
                        three_d_source_summary(&parsed),
                    )
                },
            )
        };
        let mut selected: Option<(
            usize,
            EnergyPowerAssignmentPlan,
            Option<GeneratedThreeDExpression>,
        )> = None;
        for (proposal, plan) in energy_assignment_plans.into_iter().enumerate() {
            source_options.energy_degree_bounds = Some(plan.energy_degree_bounds().to_vec());
            let key =
                ExactCffGenerationCache::generation_key(&parsed, &energy_edges, &source_options);
            let known_count = cache
                .as_deref()
                .and_then(|cache| cache.map_counts.get(&key).copied());
            let started = std::time::Instant::now();
            // A known contender needs no expression until it wins. Generated
            // expressions clone their tree containers; count-only loser records
            // avoid retaining or rebuilding those trees on repeated requests.
            let generated = if known_count.is_some() {
                None
            } else {
                Some(generate(&source_options)?)
            };
            let map_count = known_count.unwrap_or_else(|| {
                generated
                    .as_ref()
                    .expect("an unknown candidate was freshly generated")
                    .expression
                    .orientations
                    .len()
            });
            if let Some(cache) = cache.as_deref_mut() {
                cache.map_counts.insert(key, map_count);
            }
            crate::debug_tags!(#generation, #uv, #local, #four_d, #cff, #profile;
                graph = %self.name,
                proposal,
                native_source_maps = map_count,
                bounds = ?source_options.energy_degree_bounds,
                count_memo_hit = known_count.is_some(),
                elapsed_ms = started.elapsed().as_secs_f64() * 1000.0,
                "Scored bounded exact-CFF assignment proposal"
            );
            // Proposal order is the whole-envelope/deterministic tie break.
            // This finds the best of the bounded proposals, not a global optimum.
            if selected
                .as_ref()
                .is_none_or(|(best, _, _)| map_count < *best)
            {
                selected = Some((map_count, plan, generated));
            }
        }
        let (selected_count, energy_assignment_plan, generated) =
            selected.expect("rank planning always provides a baseline assignment");
        source_options.energy_degree_bounds =
            Some(energy_assignment_plan.energy_degree_bounds().to_vec());
        let key = ExactCffGenerationCache::generation_key(&parsed, &energy_edges, &source_options);
        let cached = cache.as_deref().and_then(|cache| cache.entries.get(&key));
        let cache_hit = cached.is_some();
        let generated = if let Some(generated) = generated {
            generated
        } else if let Some(cached) = cached {
            cached.clone()
        } else {
            // A previous losing key can win against a different proposal set.
            // Its count was enough to select it; obtain its payload only now.
            generate(&source_options)?
        };
        debug_assert_eq!(generated.expression.orientations.len(), selected_count);
        if !cache_hit && let Some(cache) = cache {
            // Reuse requires both canonical topology and identical occurrence
            // capacity; keep the term's assignment plan unchanged. Retain
            // only the winner: trial losers never enter the shared cache.
            cache.entries.insert(key, generated.clone());
        }
        let energy_degree_bound_report = CffEnergyDegreeBoundReport {
            source_kind: CffEnergyBoundSourceKind::ExactFourD,
            physical_parent_bounds: physical_energy_degree_bounds,
            assigned_cff_source_bounds: energy_assignment_plan.energy_degree_bounds().to_vec(),
        };
        Ok((
            generated,
            exact_source_energy_mapper,
            energy_assignment_plan,
            energy_degree_bound_report,
        ))
    }

    pub(crate) fn generate_3d_expression_for_integrand(
        &mut self,
        contract_edges: &[EdgeIndex],
        canonize_esurface: &Option<ShiftRewrite>,
        options: &Generate3DExpressionOptions,
        analysis_numerator: Option<&Atom>,
    ) -> Result<GeneratedThreeDExpression<Esurface, Hsurface>> {
        let generated = self.generate_raw_3d_expression_for_integrand(
            contract_edges,
            options,
            analysis_numerator,
        )?;
        let initial_state_cut_edges = self
            .iter_edges_of(&self.initial_state_cut)
            .map(|(_, edge, _)| edge)
            .collect::<Vec<_>>();
        self.convert_generated_expression_surfaces(
            generated,
            canonize_esurface,
            &initial_state_cut_edges,
        )
    }

    /// Generate source-owned maps before interning any physical graph surface.
    /// Bounded dispatch can discard a raw proposal without leaving graph state.
    pub(crate) fn generate_raw_3d_expression_for_integrand(
        &self,
        contract_edges: &[EdgeIndex],
        options: &Generate3DExpressionOptions,
        analysis_numerator: Option<&Atom>,
    ) -> Result<GeneratedThreeDExpression> {
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
                self.automatic_numerator_energy_degree_bounds_in_atoms_excluding_with_min_degree(
                    [numerator],
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

        Ok(generated)
    }

    pub(crate) fn production_cff_3d_expression_options(
        &self,
        settings: &GenerationSettings,
    ) -> Result<Generate3DExpressionOptions> {
        let mut options = self.cff_3d_expression_options(numerator_sampling_scale_mode(
            settings.uniform_numerator_sampling_scale,
        ))?;
        options.medium_mode = settings.medium.mode;
        Ok(options)
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
            .automatic_numerator_energy_degree_bounds_in_atoms_excluding_with_min_degree(
                [&numerator],
                cff_external_edges,
                1,
            )?;
        debug!(
            graph = %self.name,
            bounds = ?energy_degree_bounds,
            "using production CFF source-edge numerator energy-degree bounds"
        );
        Ok(Generate3DExpressionOptions {
            medium_mode: Default::default(),
            representation: RepresentationMode::Cff,
            cff_generation_context: CffGenerationContext::Standalone,
            energy_degree_bounds: Some(energy_degree_bounds),
            numerator_sampling_scale,
            preserve_internal_edges_as_four_d_denominators: Vec::new(),
        })
    }

    #[cfg(test)]
    pub(crate) fn denominator_only_cff_3d_expression_options(&self) -> Generate3DExpressionOptions {
        Generate3DExpressionOptions {
            medium_mode: Default::default(),
            representation: RepresentationMode::Cff,
            cff_generation_context: CffGenerationContext::Standalone,
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
        Ok(self
            .convert_generated_expression_surfaces_impl(
                generated,
                canonize_esurface,
                initial_state_cut_edges,
                None,
            )?
            .0)
    }

    pub(crate) fn convert_4d_expression_surfaces(
        &mut self,
        generated: GeneratedThreeDExpression,
        physical_surfaces: &[Option<LinearSurface>],
    ) -> Result<(
        GeneratedThreeDExpression<Esurface, Hsurface>,
        BTreeMap<EsurfaceID, CffEnergyFactorOwnership>,
    )> {
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
    ) -> Result<(
        GeneratedThreeDExpression<Esurface, Hsurface>,
        BTreeMap<EsurfaceID, CffEnergyFactorOwnership>,
    )> {
        let GeneratedThreeDExpression {
            mut expression,
            energy_factor_ownership,
            energy_factor_components,
            source_energy_degree_bounds,
            denominator_only_global_prefactor_sign,
            core_global_prefactor_sign,
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
        let mut linear_surface_components = BTreeMap::new();
        for (surface_id, surface) in expression.surfaces.linear_surface_cache.iter_enumerated() {
            let mut surface_component = None;
            for (edge_id, _) in &surface.expression.internal_terms {
                for (component_index, component) in energy_factor_components.iter().enumerate() {
                    if !component.internal_edge_ids.contains(&usize::from(*edge_id)) {
                        continue;
                    }
                    if let Some(previous) = surface_component
                        && previous != component_index
                    {
                        return Err(eyre::eyre!(
                            "generated causal surface {surface_id:?} spans disconnected energy-factor components"
                        ));
                    }
                    surface_component = Some(component_index);
                }
            }
            if let Some(component_index) = surface_component {
                linear_surface_components.insert(surface_id, component_index);
            }
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
        let mut mapped_surface_components = BTreeMap::<EsurfaceID, usize>::new();
        for (source_id, component_index) in linear_surface_components {
            let Some(SurfaceMapEntry {
                surface_id: HybridSurfaceID::Esurface(target_id),
                ..
            }) = linear_surface_map.get(&source_id)
            else {
                continue;
            };
            if let Some(previous) = mapped_surface_components.insert(*target_id, component_index)
                && previous != component_index
            {
                return Err(eyre::eyre!(
                    "physical causal surface {target_id:?} merges disconnected energy-factor components"
                ));
            }
        }
        let mapped_surface_ownership = mapped_surface_components
            .into_iter()
            .map(|(surface_id, component_index)| {
                (
                    surface_id,
                    energy_factor_components[component_index].ownership,
                )
            })
            .collect();
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
                match energy_factor_ownership {
                    CffEnergyFactorOwnership::GlobalSourceProduct => {
                        // GammaLoop's CFF evaluator convention keeps the
                        // on-shell-energy factors as one global product
                        // 1/prod(-2E_i) for the ordinary CFF denominator sector.
                        variant.half_edges.clear();
                    }
                    // Generalized variants retain their local positive-energy
                    // factors. Their branch-relative denominator convention is
                    // already consumed while each powered channel is lowered.
                    CffEnergyFactorOwnership::VariantLocal => {}
                }
                for remapped_denominator in signed_denominators {
                    let mut signed_variant = variant.clone();
                    signed_variant.denominator = remapped_denominator.denominator;
                    signed_variant.denominator_surface_signs = remapped_denominator.surface_signs;
                    signed_variant.prefactor *= Atom::num(remapped_denominator.prefactor);
                    remapped_variants.push(signed_variant);
                }
            }
            orientation.variants = remapped_variants;
        }

        Ok((
            GeneratedThreeDExpression {
                expression: CFFExpression {
                    orientations: expression.orientations,
                    surfaces: surface_cache,
                    residual_denominators: Vec::new(),
                },
                energy_factor_ownership,
                energy_factor_components,
                source_energy_degree_bounds,
                denominator_only_global_prefactor_sign,
                core_global_prefactor_sign,
            },
            mapped_surface_ownership,
        ))
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
        BTreeMap::<(i64, Vec<(HybridSurfaceID, i64)>), BTreeMap<Vec<HybridSurfaceID>, usize>>::new(
        );
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
            .entry(chain)
            .and_modify(|multiplicity| *multiplicity += 1)
            .or_insert(1);
    }

    if chains_by_signature.is_empty() {
        return vec![RemappedDenominatorTree {
            prefactor: 1,
            surface_signs: BTreeMap::new(),
            denominator: Tree::from_root(HybridSurfaceID::Unit),
        }];
    }

    chains_by_signature
        .into_iter()
        .flat_map(|((prefactor_sign, surface_signs), chains)| {
            // Surface projection is not injective: distinct occurrence-local
            // chains can become the same physical chain. The tree shares
            // equal branches structurally, so retain their additive count in
            // the variant prefactor and still factor chains with equal counts.
            let mut chains_by_multiplicity = BTreeMap::<usize, Vec<Vec<HybridSurfaceID>>>::new();
            for (chain, multiplicity) in chains {
                chains_by_multiplicity
                    .entry(multiplicity)
                    .or_default()
                    .push(chain);
            }
            chains_by_multiplicity
                .into_iter()
                .map(move |(multiplicity, chains)| RemappedDenominatorTree {
                    prefactor: prefactor_sign
                        * i64::try_from(multiplicity)
                            .expect("denominator-chain multiplicity should fit in i64"),
                    surface_signs: surface_signs.iter().copied().collect(),
                    denominator: denominator_tree_from_chains(&chains),
                })
        })
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
            ExactUvSubLmbFrame, FeynmanGraph, FourDDenominator,
            cuts::{CutSet, LuCutSelection},
            parse::from_dot::IntoGraph,
        },
        initialisation::test_initialise,
        settings::global::{GenerationSettings, OrientationPattern},
        utils::GS,
        uv::uv_graph::UVE,
    };
    use linnet::half_edge::subgraph::SubSetOps;
    use symbolica::atom::FunctionBuilder;
    use three_dimensional_reps::{MomentumSignature, graph_io::ParsedGraphInternalEdge};

    #[test]
    fn projected_denominator_chains_preserve_occurrence_multiplicity() {
        let local_surfaces = [0, 1, 2, 3].map(LinearSurfaceID);
        let physical_surfaces = [
            HybridSurfaceID::Esurface(EsurfaceID(0)),
            HybridSurfaceID::Esurface(EsurfaceID(1)),
            HybridSurfaceID::Esurface(EsurfaceID(2)),
        ];
        let denominator = denominator_tree_from_chains(&[
            vec![HybridSurfaceID::Linear(local_surfaces[0])],
            vec![HybridSurfaceID::Linear(local_surfaces[1])],
            vec![HybridSurfaceID::Linear(local_surfaces[2])],
            vec![HybridSurfaceID::Linear(local_surfaces[3])],
        ]);
        let surface_map = BTreeMap::from([
            (
                local_surfaces[0],
                SurfaceMapEntry {
                    surface_id: physical_surfaces[0],
                    sign: 1,
                },
            ),
            (
                local_surfaces[1],
                SurfaceMapEntry {
                    surface_id: physical_surfaces[0],
                    sign: 1,
                },
            ),
            (
                local_surfaces[2],
                SurfaceMapEntry {
                    surface_id: physical_surfaces[1],
                    sign: 1,
                },
            ),
            (
                local_surfaces[3],
                SurfaceMapEntry {
                    surface_id: physical_surfaces[2],
                    sign: 1,
                },
            ),
        ]);

        let remapped = remap_denominator_tree_surface_ids(&denominator, &surface_map);
        let actual = remapped.iter().fold(Atom::new(), |sum, term| {
            sum + Atom::num(term.prefactor) * term.denominator.to_atom_inv()
        });
        let expected = Atom::num(2) / Atom::from(physical_surfaces[0])
            + Atom::num(1) / Atom::from(physical_surfaces[1])
            + Atom::num(1) / Atom::from(physical_surfaces[2]);

        assert!(
            (actual.clone() - expected).expand().is_zero(),
            "projecting distinct local chains onto one physical chain lost multiplicity: {actual}"
        );
        assert_eq!(remapped.len(), 2);
        assert_eq!(
            remapped
                .iter()
                .find(|term| term.prefactor == 1)
                .expect("unit-multiplicity chains should stay factorized")
                .denominator
                .get_bottom_layer()
                .len(),
            2,
        );
        assert!(remapped.iter().all(|term| term.surface_signs.is_empty()));

        let negative_map = BTreeMap::from([
            (
                local_surfaces[0],
                SurfaceMapEntry {
                    surface_id: physical_surfaces[0],
                    sign: -1,
                },
            ),
            (
                local_surfaces[1],
                SurfaceMapEntry {
                    surface_id: physical_surfaces[0],
                    sign: -1,
                },
            ),
        ]);
        let negative_denominator = denominator_tree_from_chains(&[
            vec![HybridSurfaceID::Linear(local_surfaces[0])],
            vec![HybridSurfaceID::Linear(local_surfaces[1])],
        ]);
        let negative = remap_denominator_tree_surface_ids(&negative_denominator, &negative_map);
        assert_eq!(negative.len(), 1);
        assert_eq!(negative[0].prefactor, -2);
        assert_eq!(
            negative[0].surface_signs,
            BTreeMap::from([(physical_surfaces[0], -1)])
        );
    }

    #[test]
    fn exact_cff_cache_preserves_occurrence_capacities() {
        let parsed = ParsedGraph {
            internal_edges: (0..3)
                .map(|edge_id| ParsedGraphInternalEdge {
                    edge_id,
                    tail: edge_id,
                    head: edge_id + 1,
                    label: String::new(),
                    mass_key: Some("muv".to_owned()),
                    signature: MomentumSignature {
                        loop_signature: vec![1],
                        external_signature: vec![],
                    },
                    had_pow: false,
                })
                .collect(),
            external_edges: vec![],
            initial_state_cut_edges: vec![],
            loop_names: vec!["k".to_owned()],
            external_names: vec![],
            node_name_to_internal: BTreeMap::new(),
        };

        let energy_edges = EnergyEdgeIndexMap::identity(3, 0);
        let requests = [vec![(0, 4)], vec![(1, 4)], vec![(0, 2), (1, 2)]].map(|bounds| {
            Generate3DExpressionOptions {
                energy_degree_bounds: Some(bounds),
                ..Default::default()
            }
        });
        let keys = requests
            .iter()
            .map(|request| ExactCffGenerationCache::generation_key(&parsed, &energy_edges, request))
            .collect::<std::collections::BTreeSet<_>>();
        assert_eq!(
            keys.len(),
            requests.len(),
            "equal channel totals must not merge different occurrence capacities",
        );
    }

    #[test]
    fn exact_cff_batch_reuses_owner_relabelled_sub_lmb_topology() -> Result<()> {
        test_initialise()?;
        let graph: Graph = dot!(digraph exact_owner_relabelled_cache {
            edge [num=1 mass=1]
            node [num=1]

            a -> b [id=0 lmb_id=0]
            a -> b [id=1]
        })?;
        let denominators = [EdgeIndex(0), EdgeIndex(1)].map(|source_edge| FourDDenominator {
            source_edge,
            momentum: FunctionBuilder::new(GS.emr_mom)
                .add_arg(usize::from(source_edge))
                .finish(),
            mass_squared: graph.underlying[source_edge].mass_atom().pow(2),
            full_expr: Atom::one(),
        });
        let mut relabelled = denominators.clone();
        relabelled[1].source_edge = EdgeIndex(0);
        let source = GraphThreeDSource::from_exact_denominators_in_uv_sub_lmb(
            &graph,
            &denominators,
            [EdgeIndex(0), EdgeIndex(1)],
            [],
            &graph.loop_momentum_basis,
            ExactUvSubLmbFrame::TaylorVacuum,
        )?;
        let relabelled_source = GraphThreeDSource::from_exact_denominators_in_uv_sub_lmb(
            &graph,
            &relabelled,
            [EdgeIndex(0), EdgeIndex(1)],
            [],
            &graph.loop_momentum_basis,
            ExactUvSubLmbFrame::TaylorVacuum,
        )?;
        let parsed = source.to_three_d_parsed_graph()?;
        let relabelled_parsed = relabelled_source.to_three_d_parsed_graph()?;
        let energy_edges = source
            .energy_edge_index_map(&parsed)
            .expect("exact source has an occurrence-local energy map");
        let relabelled_energy_edges = relabelled_source
            .energy_edge_index_map(&relabelled_parsed)
            .expect("relabelled exact source has an occurrence-local energy map");
        let options = graph.denominator_only_cff_3d_expression_options();

        assert_eq!(energy_edges, relabelled_energy_edges);
        assert_eq!(
            ExactCffGenerationCache::generation_key(&parsed, &energy_edges, &options),
            ExactCffGenerationCache::generation_key(
                &relabelled_parsed,
                &relabelled_energy_edges,
                &options,
            ),
            "physical owner provenance must not enter the canonical exact-CFF cache key",
        );

        let mut cache = ExactCffGenerationCache::default();
        graph.generate_3d_expression_for_4d_term(
            &source,
            &options,
            &Atom::one(),
            Some(&mut cache),
        )?;
        graph.generate_3d_expression_for_4d_term(
            &relabelled_source,
            &options,
            &Atom::one(),
            Some(&mut cache),
        )?;
        assert_eq!(cache.len(), 1);
        Ok(())
    }

    #[test]
    fn bounded_exact_cff_dispatch_selects_actual_count_and_preserves_owner_contour() -> Result<()> {
        use crate::cff::{
            expression::GammaLoopOrientationExpression, surface::GammaLoopSurfaceCache,
        };
        use std::collections::BTreeSet;
        use three_dimensional_reps::CffGlobalPrefactorSign;

        test_initialise()?;
        let graph: Graph = dot!(digraph bounded_exact_dispatch {
            edge [num=1 mass=2]
            node [num=1]

            a -> b [id=0 lmb_id=0]
            b -> a [id=1]
            a -> a [id=2 lmb_id=1]
        })?;
        let owners = [EdgeIndex(0), EdgeIndex(1), EdgeIndex(2)];
        let base = owners.map(|source_edge| FourDDenominator {
            source_edge,
            momentum: graph.loop_momentum_basis.loop_atom::<Atom>(
                source_edge,
                GS.emr_mom,
                &[],
                true,
            ),
            mass_squared: graph.underlying[source_edge].mass_atom().pow(2),
            full_expr: Atom::one(),
        });
        let options = Generate3DExpressionOptions {
            cff_generation_context: CffGenerationContext::EmbeddedCffFactor,
            numerator_sampling_scale: NumeratorSamplingScaleMode::None,
            ..graph.denominator_only_cff_3d_expression_options()
        };
        for dotted_owner in [0, 1] {
            let owner = owners[dotted_owner];
            let mut denominators = base.to_vec();
            denominators.extend([base[dotted_owner].clone(), base[dotted_owner].clone()]);
            let source = GraphThreeDSource::from_exact_denominators_in_uv_sub_lmb(
                &graph,
                &denominators,
                owners,
                [],
                &graph.loop_momentum_basis,
                ExactUvSubLmbFrame::TaylorVacuum,
            )?;
            let tagged = |derived| {
                FunctionBuilder::new(GS.emr_mom)
                    .add_arg(GS.uv_momentum_provenance_tag(
                        Atom::num(dotted_owner as i64).as_view(),
                        derived,
                        base[dotted_owner].momentum.as_view(),
                    ))
                    .add_arg(GS.cind(0))
                    .finish()
            };
            let fixed = tagged(false);
            let derived = tagged(true);
            // Keep the original quartic factor on its retained occurrence; only
            // the two new denominator-derived powers may use the serial copies.
            let numerator = fixed.clone().pow(4) * derived.clone().pow(2);
            let preparation =
                graph.prepare_3d_expression_for_4d_term(&source, &options, &numerator)?;
            assert_eq!(preparation.energy_assignment_plans.len(), 3);
            let fixed_preparation =
                graph.prepare_3d_expression_for_4d_term(&source, &options, &fixed)?;
            let fixed_occurrence =
                fixed_preparation.energy_assignment_plans[0].energy_degree_bounds()[0].0;
            let physical_edges = source.physical_energy_edge_index_map().unwrap();
            let owner_occurrences = physical_edges
                .internal
                .iter()
                .filter_map(|(occurrence, physical)| {
                    (*physical == dotted_owner).then_some(*occurrence)
                })
                .collect::<BTreeSet<_>>();

            // Generate the complete reference catalogues independently of
            // the production selector. Do not replace this with mock costs or
            // an envelope formula: topology and full-map collisions own the count.
            let mut references = Vec::new();
            for plan in &preparation.energy_assignment_plans {
                let mut fixed_uses = 0;
                let mut derived_uses = 0;
                let reconstructed = plan
                    .map_factors(|factor, assignments| {
                        if factor == &fixed {
                            assert_eq!(assignments[&owner], fixed_occurrence);
                            fixed_uses += 1;
                        } else if factor == &derived {
                            assert_ne!(assignments[&owner], fixed_occurrence);
                            assert!(owner_occurrences.contains(&assignments[&owner]));
                            derived_uses += 1;
                        }
                        Ok::<_, ()>(factor.clone())
                    })
                    .unwrap();
                assert_eq!(reconstructed, numerator);
                assert_eq!((fixed_uses, derived_uses), (4, 2));
                let mut request = preparation.source_options.clone();
                request.energy_degree_bounds = Some(plan.energy_degree_bounds().to_vec());
                let generated = three_dimensional_reps::generate_3d_expression(&source, &request)?;
                references.push((
                    plan.energy_degree_bounds().to_vec(),
                    generated.expression.orientations.len(),
                    generated
                        .expression
                        .orientations
                        .iter()
                        .map(|orientation| orientation.residue_map_key())
                        .collect::<BTreeSet<_>>(),
                ));
            }
            let winner = (0..references.len())
                .min_by_key(|index| (references[*index].1, *index))
                .unwrap();
            // Copy assignments can tie when neither copy owns the loop basis.
            // Test the exact selector, report and cache invariants in both
            // source topologies without claiming a speedup from that tie.
            let mut cache = ExactCffGenerationCache::default();
            for phase in ["fresh", "full cache hit", "count memo without payload"] {
                if phase == "count memo without payload" {
                    cache.entries.clear();
                }
                let (generated, mapper, plan, report) = graph.generate_3d_expression_for_4d_term(
                    &source,
                    &options,
                    &numerator,
                    Some(&mut cache),
                )?;
                assert_eq!(plan.energy_degree_bounds(), references[winner].0.as_slice());
                assert_eq!(generated.source_energy_degree_bounds, references[winner].0);
                assert_eq!(report.assigned_cff_source_bounds, references[winner].0);
                assert_eq!(
                    report.physical_parent_bounds,
                    preparation.physical_energy_degree_bounds
                );
                assert_eq!(
                    generated.expression.orientations.len(),
                    references[winner].1
                );
                assert_eq!(
                    generated
                        .expression
                        .orientations
                        .iter()
                        .map(|orientation| orientation.residue_map_key())
                        .collect::<BTreeSet<_>>(),
                    references[winner].2,
                );
                assert_eq!(cache.entries.len(), 1);
                assert_eq!(cache.map_counts.len(), 3);

                let replacements = generated.expression.surfaces.get_all_replacements_gs(&[]);
                let mut contour = Atom::Zero;
                for orientation in &generated.expression.orientations {
                    let mapped = mapper.map_planned_numerator(
                        &orientation.loop_energy_map,
                        &orientation.edge_energy_map,
                        &plan,
                    )?;
                    contour += orientation
                        .to_atom_gs()
                        .replace_multiple(&replacements)
                        .replace_multiple(mapper.exact_ose_replacements())
                        * mapped;
                }
                let source_sign =
                    generated
                        .energy_factor_components
                        .iter()
                        .fold(1, |sign, component| {
                            let frame = match component.ownership {
                                CffEnergyFactorOwnership::GlobalSourceProduct => {
                                    component.core_global_prefactor_sign
                                }
                                CffEnergyFactorOwnership::VariantLocal => {
                                    component.denominator_only_global_prefactor_sign
                                }
                            };
                            sign * CffGlobalPrefactorSign::from_exponent(
                                component.internal_edge_ids.len(),
                            )
                            .product(frame)
                            .factor()
                        });
                contour *= Atom::num(source_sign);
                let energy = |physical| {
                    let occurrence = physical_edges
                        .internal
                        .iter()
                        .find_map(|(occurrence, owner)| (*owner == physical).then_some(*occurrence))
                        .unwrap();
                    crate::utils::ose_atom_from_index(EdgeIndex(occurrence))
                        .replace_multiple(mapper.exact_ose_replacements())
                };
                // Independent clockwise Below contours are -5/(32E_b)
                // for q^6/D_b^4 and -1/(2E_t) for the attached tadpole.
                // The original quartic factor stays on its original occurrence.
                let expected = Atom::num(5) / (Atom::num(64) * energy(2) * energy(dotted_owner));
                // The exact mapper keeps an unchanged propagator as E(owner),
                // while a rewritten equal-momentum denominator uses its literal
                // square root. The repeated-pole oracle requires the same positive
                // physical mass shell for both dialects.
                let physical_shell = |mut expression: Atom| {
                    for physical in owners {
                        let on_shell_energy = (1..=3)
                            .fold(
                                graph.underlying[physical].mass_atom().pow(2),
                                |square, spatial_index| {
                                    square
                                        + graph
                                            .loop_momentum_basis
                                            .loop_atom::<Atom>(
                                                physical,
                                                GS.emr_mom,
                                                &[GS.cind(spatial_index)],
                                                true,
                                            )
                                            .pow(2)
                                },
                            )
                            .sqrt();
                        expression = expression.replace(GS.ose(physical)).with(on_shell_energy);
                    }
                    expression
                };
                let difference = (physical_shell(contour) - physical_shell(expected)).together();
                assert!(
                    difference.is_zero(),
                    "dotted owner {dotted_owner}, {phase}: dispatch changed the analytic contour: {difference}"
                );
            }
        }
        Ok(())
    }

    #[test]
    fn direct_cff_carries_core_global_sign_as_production_metadata() -> Result<()> {
        test_initialise()?;
        let mut graph: Graph = dot!(digraph direct_duplicate_sign {
            edge [num=1 mass=1]
            node [num=1]

            a -> b [id=0 lmb_id=0]
            a -> b [id=1]
        })?;
        let options = graph.denominator_only_cff_3d_expression_options();
        let source = GraphThreeDSource::new(&graph, &[])?;
        let shared_core = three_dimensional_reps::generate_3d_expression(&source, &options)?;

        assert_eq!(shared_core.core_global_prefactor_sign.factor(), -1);
        assert!(
            shared_core
                .expression
                .orientations
                .iter()
                .flat_map(|orientation| &orientation.variants)
                .all(|variant| (&variant.prefactor + Atom::one()).is_zero()),
            "the shared-core global sign must remain in the generated expression"
        );

        let cutset = CutSet::empty(graph.n_hedges());
        let contract: linnet::half_edge::subgraph::SuBitGraph = graph.empty_subgraph();
        let cff = graph.cff(
            &contract,
            &cutset,
            &OrientationPattern::default(),
            &options,
            None,
        )?;
        assert_eq!(cff.production_prefactor_factor(), -1);

        let contracted_loop = graph
            .get_edge_subgraph(EdgeIndex(0))
            .union(&graph.get_edge_subgraph(EdgeIndex(1)));
        assert_eq!(graph.cyclotomatic_number(&contracted_loop), 1);
        let reduced = graph.cff(
            &contracted_loop,
            &cutset,
            &OrientationPattern::default(),
            &options,
            None,
        )?;
        assert_eq!(
            reduced.production_prefactor_factor(),
            1,
            "contracting one loop must regenerate the active one-loop source convention"
        );
        Ok(())
    }

    #[test]
    fn two_loop_cff_carries_core_loop_parity_as_production_metadata() -> Result<()> {
        test_initialise()?;
        let mut graph: Graph = dot!(digraph direct_two_loop_sign {
            edge [num=1 mass=1]
            node [num=1]

            a -> b [id=0 lmb_id=0]
            a -> b [id=1 lmb_id=1]
            b -> a [id=2]
        })?;
        let options = graph.denominator_only_cff_3d_expression_options();
        let source = GraphThreeDSource::new(&graph, &[])?;
        let shared_core = three_dimensional_reps::generate_3d_expression(&source, &options)?;

        assert_eq!(shared_core.core_global_prefactor_sign.factor(), -1);
        assert!(
            shared_core
                .expression
                .orientations
                .iter()
                .flat_map(|orientation| &orientation.variants)
                .all(|variant| (&variant.prefactor + Atom::one()).is_zero()),
            "the connected two-loop core parity must remain in the generated expression"
        );

        let cutset = CutSet::empty(graph.n_hedges());
        let contract: linnet::half_edge::subgraph::SuBitGraph = graph.empty_subgraph();
        let cff = graph.cff(
            &contract,
            &cutset,
            &OrientationPattern::default(),
            &options,
            None,
        )?;
        assert_eq!(cff.production_prefactor_factor(), -1);
        Ok(())
    }

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
                A -> B [id=0 num="Q(0,spenso::cind(0))+1"]
                B -> C [id=1 num="Q(1,spenso::cind(0))+2"]
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
    fn selected_raw_cff_matches_normal_generation_after_discarding_a_trial() -> Result<()> {
        test_initialise()?;
        let mut graph: Graph = dot!(digraph selected_raw_cff {
            edge [num=1 mass=1]
            node [num=1]
            a -> b [id=0 lmb_id=0]
            a -> b [id=1]
            b -> b [id=2 lmb_id=1]
        })?;
        let mut ordinary_graph = graph.clone();
        let contract = graph.get_edge_subgraph(EdgeIndex(2));
        let options = Generate3DExpressionOptions {
            cff_generation_context: CffGenerationContext::EmbeddedCffFactor,
            ..graph.denominator_only_cff_3d_expression_options()
        };
        let numerator = GS.emr_mom(EdgeIndex(0), GS.cind(0)).pow(2);
        let untouched = bincode::encode_to_vec(&graph.surface_cache, bincode::config::standard())?;
        let discarded = graph.generate_raw_3d_expression_for_integrand(
            &[EdgeIndex(2)],
            &options,
            Some(&GS.emr_mom(EdgeIndex(0), GS.cind(0)).pow(3)),
        )?;
        drop(discarded);
        let selected = graph.generate_raw_3d_expression_for_integrand(
            &[EdgeIndex(2)],
            &options,
            Some(&numerator),
        )?;
        assert_eq!(
            bincode::encode_to_vec(&graph.surface_cache, bincode::config::standard())?,
            untouched,
            "discarded and selected raw proposals must not intern E, H or linear surfaces",
        );
        let canonization = graph.get_esurface_canonization(&graph.loop_momentum_basis);
        let selected = graph.convert_generated_expression_surfaces(selected, &canonization, &[])?;
        let cutset = CutSet::empty(graph.n_hedges());
        let pattern = OrientationPattern::default();
        let selected = graph
            .cff_from_generated_expression(selected, &contract, &cutset, &pattern, &options)?;
        let ordinary =
            ordinary_graph.cff(&contract, &cutset, &pattern, &options, Some(&numerator))?;
        assert_eq!(
            selected.energy_degree_bound_report,
            ordinary.energy_degree_bound_report
        );
        assert_eq!(
            selected.production_prefactor_factor(),
            ordinary.production_prefactor_factor()
        );
        assert_eq!(
            selected.terms.keys().collect::<Vec<_>>(),
            ordinary.terms.keys().collect::<Vec<_>>()
        );
        for (index, selected_term) in selected.terms {
            let ordinary_term = &ordinary.terms[&index];
            assert_eq!(
                selected_term.orientations.len(),
                ordinary_term.orientations.len()
            );
            for (selected, ordinary) in selected_term
                .orientations
                .iter()
                .zip(&ordinary_term.orientations)
            {
                assert_eq!(
                    selected.orientation.residue_map_key(),
                    ordinary.orientation.residue_map_key()
                );
                assert_eq!(selected.expression, ordinary.expression);
                assert_eq!(
                    selected.production_orientation_id,
                    ordinary.production_orientation_id
                );
            }
        }
        assert_eq!(
            bincode::encode_to_vec(&graph.surface_cache, bincode::config::standard())?,
            bincode::encode_to_vec(&ordinary_graph.surface_cache, bincode::config::standard())?,
            "only the chosen raw source may determine the persistent surface IDs",
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

    #[test]
    fn exact_massless_raised_lu_cff_retains_selected_residue() -> Result<()> {
        test_initialise()?;
        let mut graph: Graph = dot!(digraph exact_massless_raised_lu {
            edge [num=1 mass="UFO::ZERO"]
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
        let production = graph
            .generate_3d_expression_for_integrand(&[], &canonization, &options, Some(&numerator))?
            .expression;
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
        let denominators = [1, 2, 3].map(|edge| {
            let edge = EdgeIndex(edge);
            FourDDenominator {
                source_edge: edge,
                momentum: FunctionBuilder::new(GS.emr_mom)
                    .add_arg(usize::from(edge))
                    .finish(),
                mass_squared: graph.underlying[edge].mass_atom().pow(2),
                full_expr: Atom::one(),
            }
        });
        assert!(
            denominators
                .iter()
                .all(|denominator| denominator.mass_squared.is_zero()),
            "UFO::ZERO masses must be normalized before exact surface projection"
        );
        let source = GraphThreeDSource::from_exact_denominators(&graph, &denominators)?;
        let mapper = source
            .exact_source_energy_mapper()
            .expect("exact source has a physical EMR mapper");
        let candidates = mapper.equivalent_energy_candidates([EdgeIndex(1)])?;
        let plan = graph.plan_numerator_energy_assignment_in_atom_excluding(
            &numerator,
            std::iter::empty(),
            &candidates,
        )?;
        let assigned_occurrence = source
            .physical_energy_edge_index_map()
            .expect("exact source has a physical occurrence map")
            .internal
            .into_iter()
            .find_map(|(occurrence, owner)| (owner == 1).then_some(occurrence))
            .expect("the physical numerator carrier has one exact occurrence");
        assert_eq!(
            plan.energy_degree_bounds(),
            &[(assigned_occurrence, 2)],
            "canonical exact-edge IDs may change, but the quadratic bound must remain on the occurrence owned by physical edge 1",
        );

        let (exact, _) =
            graph.cff_from_4d_denominators(&denominators, &cutset, &options, &numerator)?;
        assert_eq!(
            exact.energy_degree_bound_report,
            CffEnergyDegreeBoundReport {
                source_kind: CffEnergyBoundSourceKind::ExactFourD,
                physical_parent_bounds: vec![(1, 2)],
                assigned_cff_source_bounds: plan.energy_degree_bounds().to_vec(),
            },
            "exact CFF diagnostics must keep physical-parent and assigned occurrence bounds distinct",
        );
        assert_eq!(exact.terms.len(), lu_cut.max_occurence);
        let selected = exact
            .terms
            .iter()
            .find(|(index, _)| index.lu_cut_order == Some(2))
            .expect("the exact raised CFF should retain its second LU residue");
        assert!(!selected.1.orientations.is_empty());
        assert!(
            selected
                .1
                .orientations
                .iter()
                .flat_map(|orientation| &orientation.orientation.variants)
                .next()
                .is_some(),
            "the selected exact LU residue must retain a causal variant"
        );
        Ok(())
    }
}
