use std::{
    collections::BTreeMap,
    ops::Neg,
    sync::Arc,
    time::{Duration, Instant},
};

use color_eyre::Result;
use eyre::eyre;
use linnet::half_edge::subgraph::SubSetLike;
use symbolica::atom::Atom;
use three_dimensional_reps::CffGenerationContext;

use crate::{
    cff::{
        CutCFFIndex,
        generation::{PreparationKey, PreparationValue},
    },
    debug_tags,
    graph::{ExactUvSubLmbFrame, Graph, cuts::CutSet},
    numerator::energy_degree::{EnergyPowerAnalyzer, EnergyPowerCapMap},
    utils::GS,
    uv::{
        UVgenerationSettings,
        approx::{
            local_3d::Localizer,
            local_4d::{CanonicalUvSector, FourDSector, Local4dCts},
        },
    },
};

/// Projected local-4D Taylor coefficients. They omit the untouched outer CFF,
/// which is attached only during final assembly. Each child contour has already
/// summed its own source maps, so coefficients carry no production host or cut key.
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub(crate) struct Projected4dCts(Vec<Projected4dSector>);

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub(crate) struct Projected4dSector {
    pub(crate) coefficient: Atom,
    /// Completed-loop localizers are multiplied after outer numerator mapping.
    pub(crate) frozen_factor: Atom,
}

impl Projected4dCts {
    pub(super) fn new(sectors: Vec<Projected4dSector>) -> Self {
        Self(sectors)
    }

    pub(crate) fn sectors(&self) -> &[Projected4dSector] {
        &self.0
    }

    #[cfg(test)]
    pub(crate) fn map<F: FnMut(&Atom) -> Result<Atom>>(&self, mut f: F) -> Result<Self> {
        Ok(Self(
            self.0
                .iter()
                .map(|sector| {
                    Ok(Projected4dSector {
                        coefficient: f(&sector.coefficient)?,
                        frozen_factor: sector.frozen_factor.clone(),
                    })
                })
                .collect::<Result<_>>()?,
        ))
    }
}

impl Neg for Projected4dCts {
    type Output = Self;

    fn neg(self) -> Self::Output {
        Self(
            self.0
                .into_iter()
                .map(|sector| Projected4dSector {
                    coefficient: -sector.coefficient,
                    frozen_factor: sector.frozen_factor,
                })
                .collect(),
        )
    }
}

/// Nonserialized preparation owned by one graph's complete UV computation.
/// Independent contours share reusable source payloads, never coefficients.
pub(crate) struct Local4dProjectionContext {
    pub(crate) generation_cache: crate::cff::generation::ExactCffGenerationCache,
    pub(crate) preparations:
        crate::cff::generation::GenerationCache<PreparationKey, PreparationValue>,
    pub(crate) numerator_rows:
        crate::cff::generation::GenerationCache<crate::cff::ExactNumeratorRowKey, Atom>,
    pub(crate) numerator_template_builds: usize,
    pub(crate) canonical_preparation_builds: usize,
    pub(crate) source_preparation_builds: usize,
    pub(crate) projected_component_requests: usize,
}

impl Default for Local4dProjectionContext {
    fn default() -> Self {
        Self {
            generation_cache: Default::default(),
            preparations: crate::cff::generation::GenerationCache::new(48 * 1024 * 1024, 4096),
            numerator_rows: crate::cff::generation::GenerationCache::new(16 * 1024 * 1024, 16384),
            numerator_template_builds: 0,
            canonical_preparation_builds: 0,
            source_preparation_builds: 0,
            projected_component_requests: 0,
        }
    }
}

impl Local4dProjectionContext {
    fn canonical_projection(
        &mut self,
        sector: FourDSector,
        graph: &Graph,
    ) -> Result<Arc<CanonicalUvSector>> {
        let started = Instant::now();
        let key = PreparationKey::Canonical(Arc::new(sector));
        if let Some(PreparationValue::Canonical(sector)) = self.preparations.get(&key) {
            let sector = Arc::clone(sector);
            drop(key);
            self.preparations.cache_time += started.elapsed();
            return Ok(sector);
        }
        let PreparationKey::Canonical(raw) = &key else {
            unreachable!()
        };
        let projection_started = Instant::now();
        let sector = Arc::new(raw.canonical_projection(graph)?);
        let projection_time = projection_started.elapsed();
        let bytes = raw.accounted_bytes()
            + sector.accounted_bytes()
            + 2 * std::mem::size_of::<PreparationKey>()
            + std::mem::size_of::<PreparationValue>()
            + 128;
        self.canonical_preparation_builds += 1;
        self.preparations
            .insert(key, PreparationValue::Canonical(Arc::clone(&sector)), bytes);
        self.preparations.cache_time += started.elapsed().saturating_sub(projection_time);
        Ok(sector)
    }
}

pub(crate) struct Projected4dApproximation<'a> {
    localizer: Localizer<'a>,
    graph: &'a mut Graph,
    settings: &'a UVgenerationSettings,
}

impl<'a> Projected4dApproximation<'a> {
    pub(crate) fn new(
        localizer: Localizer<'a>,
        graph: &'a mut Graph,
        settings: &'a UVgenerationSettings,
    ) -> Self {
        Self {
            localizer,
            graph,
            settings,
        }
    }
}

impl Localizer<'_> {
    /// Project one typed Taylor sector in the independent sub-LMB retained by
    /// each active component. Denominator ownership selects the component;
    /// original graph incidence and the stored coordinate LMB remain the sole
    /// topology and momentum authorities throughout the exact residue chain.
    pub(super) fn project_factorized_taylor_sector(
        self,
        graph: &mut Graph,
        sector: &CanonicalUvSector,
        context: &mut Local4dProjectionContext,
    ) -> Result<Atom> {
        if sector.active_components.is_empty() {
            return Err(eyre!(
                "typed local-4D projection requires at least one active Taylor component"
            ));
        }

        let child_cutset = CutSet::empty(graph.n_hedges());
        // Each Taylor component is an independently closed contour which is
        // multiplied into the outer CFF afterwards. EmbeddedCffFactor keeps
        // every causal term when finite-pole denominators remain. Only for a
        // connected, nonrepeated, full-rank terminal residue basis does it avoid
        // reopening equivalent closures: one deterministic Below representative
        // is then the complete integral.
        let mut options = self.orientation.cff_options()?;
        options.cff_generation_context = CffGenerationContext::EmbeddedCffFactor;
        let excluded_report_owners = graph
            .iter_edges_of(&graph.initial_state_cut)
            .chain(graph.iter_edges_of(&graph.tree_edges))
            .map(|(_, edge, _)| edge)
            .collect::<std::collections::BTreeSet<_>>();
        let mut terms = Vec::new();
        for term in &sector.terms {
            // Raw physical ownership is diagnostic provenance, separate from
            // canonical class algebra and the selected native capacities.
            // Compute it once, before residue states multiply across waves.
            let raw_degree_started = Instant::now();
            let raw_reports = sector
                .active_components
                .iter()
                .map(|(owners, _, _)| {
                    let owners = graph
                        .iter_edges_of(owners)
                        .filter_map(|(pair, edge, data)| {
                            (pair.is_paired()
                                && !data.data.is_dummy
                                && !excluded_report_owners.contains(&edge))
                            .then_some(edge)
                        });
                    EnergyPowerAnalyzer::for_physical_emr_edges(owners)
                        .analyze_atom(&term.source_numerator)
                })
                .collect::<Result<Vec<_>, _>>()?;
            debug_tags!(#generation, #profile, #uv, #local, #four_d;
                stage = "raw_physical_degree_report",
                components = raw_reports.len(),
                elapsed_ms = raw_degree_started.elapsed().as_secs_f64() * 1000.0,
                "Analyzed raw physical ownership for degree reporting"
            );
            let mut component_denominators = vec![Vec::new(); sector.active_components.len()];
            let mut residual_factor = Atom::one();
            for (denominator, class) in term.source_witness.iter().zip(&term.source_classes) {
                if let Some((class, _)) = class {
                    component_denominators[sector.class(*class).component]
                        .push(denominator.clone());
                } else {
                    // A denominator can be independent of its component's
                    // quotient loop after a pinch. Its full polynomial remains
                    // an ordinary factor, with the original physical binding.
                    residual_factor /= &denominator.full_expr;
                }
            }
            let classes = sector
                .classes
                .iter()
                .filter(|class| term.powers.contains_key(&class.id))
                .cloned()
                .collect::<Vec<_>>();
            terms.push((
                component_denominators,
                classes,
                term.powers.clone(),
                raw_reports,
                vec![(residual_factor, term.numerator.clone())],
            ));
        }
        // Each Taylor component owns an independent energy contour. Every
        // state entering one component retains its own occurrence capacities.
        // Reuse requires equal topology and capacity, so no registration pass
        // or combined bound envelope is needed before generation.
        // The graph-owned cache retains immutable source payloads across
        // waves; each request still carries its own numerator assignment.
        for (component, (_, source_scope, coordinate_lmb)) in
            sector.active_components.iter().enumerate()
        {
            let mut composition_time = Duration::ZERO;
            let uv_edges = graph
                .iter_edges_of(source_scope)
                .filter_map(|(pair, edge_id, edge)| {
                    (pair.is_paired() && !edge.data.is_dummy).then_some(edge_id)
                })
                .collect::<Vec<_>>();

            // Use the complete enclosing source scope for incidence:
            // exact-source reconstruction contracts the omitted prefix into
            // the quotient topology. Only the reduced owner set above selects
            // this component's active denominators.
            for (component_denominators, classes, _, raw_reports, states) in &mut terms {
                let denominators = &component_denominators[component];
                let active_classes = classes
                    .iter()
                    .filter(|class| class.component == component)
                    .cloned()
                    .collect::<Vec<_>>();
                if denominators.is_empty() {
                    return Err(eyre!(
                        "active Taylor component has no energy denominator in one 4D term"
                    ));
                }
                debug_tags!(#generation, #uv, #local, #four_d, #cff, #trace;
                    component,
                    source_scope = %source_scope.string_label(),
                    coordinate_lmb = ?coordinate_lmb,
                    denominators = ?denominators,
                    "Projecting one factorized local-4D Taylor component"
                );
                let mut next_states = Vec::new();
                for (carrier, numerator) in std::mem::take(states) {
                    debug_tags!(#generation, #uv, #local, #four_d, #cff, #trace;
                        component,
                        source_scope = %source_scope.string_label(),
                        log.numerator = &numerator,
                        "Preparing factorized local-4D child numerator for component {component} in {}: {}",
                        source_scope.string_label(),
                        numerator,
                    );
                    context.projected_component_requests += 1;
                    let (mut cff, _) = graph.cff_from_4d_denominators_in_uv_sub_lmb(
                        denominators,
                        uv_edges.iter().copied(),
                        [],
                        coordinate_lmb,
                        ExactUvSubLmbFrame::TaylorVacuum,
                        &child_cutset,
                        &options,
                        &numerator,
                        Some(&mut *context),
                        &active_classes,
                    )?;
                    cff.energy_degree_bound_report.physical_parent_bounds =
                        raw_reports[component].clone().into_generation_bounds();
                    self.orientation
                        .record_energy_degree_bound_report(&cff.energy_degree_bound_report);
                    // The sector already carries its forest subtraction sign.
                    // This bridge only converts this independent component's
                    // CFF energy-factor convention to the production convention.
                    let production_prefactor = Atom::num(cff.production_prefactor_factor());
                    for (index, cff_term) in cff.terms {
                        if index != CutCFFIndex::new_all_none() {
                            return Err(eyre!(
                                "an uncut UV-child CFF unexpectedly produced residue index {index}"
                            ));
                        }
                        for orientation in &cff_term.orientations {
                            let mapped_numerator = cff_term
                                .map_exact_source_numerator(
                                    &orientation.orientation,
                                    Some(&mut *context),
                                )
                                .map_err(|error| {
                                    eyre!(
                                        "{error}; exact UV-child component denominators are {:?}",
                                        denominators
                                    )
                                })?;
                            debug_tags!(#generation, #uv, #local, #four_d, #cff, #trace;
                                component,
                                source_scope = %source_scope.string_label(),
                                log.mapped_numerator = mapped_numerator,
                                "Mapped factorized local-4D child numerator for component {component} in {}: {}",
                                source_scope.string_label(),
                                mapped_numerator,
                            );
                            let composition_started = Instant::now();
                            next_states.push((
                                &carrier * &orientation.expression * &production_prefactor,
                                mapped_numerator,
                            ));
                            composition_time += composition_started.elapsed();
                        }
                    }
                }
                *states = next_states;
            }

            // Once a component is integrated, its old powers no longer
            // distinguish future requests. Canonical class IDs retain the
            // exact routing, mass and domain inside this common sector frame.
            // Keep the first deterministic witness, never concatenate source
            // denominators; the next source certifies the combined numerator.
            let composition_started = Instant::now();
            let mut remaining_requests = BTreeMap::new();
            for (mut denominators, mut classes, mut powers, mut raw_reports, states) in terms {
                denominators[component].clear();
                classes.retain(|class| class.component > component);
                powers.retain(|id, _| sector.class(*id).component > component);
                raw_reports[component] = Default::default();
                let (_, _, _, reports, combined) =
                    remaining_requests.entry(powers.clone()).or_insert_with(|| {
                        (
                            denominators,
                            classes,
                            powers,
                            vec![EnergyPowerCapMap::default(); raw_reports.len()],
                            Vec::new(),
                        )
                    });
                for (report, incoming) in reports.iter_mut().zip(raw_reports) {
                    report.max_assign(incoming);
                }
                combined.extend(states);
            }
            terms = remaining_requests.into_values().collect();
            for (_, _, _, _, states) in &mut terms {
                // Equal numerators share their summed carrier. Conversely,
                // equal carriers share one factorized numerator. All states here
                // have the same remaining component requests and contour frame,
                // including states from different original rational terms.
                // Combine numerator keys first: carrier-first grouping would
                // turn repeated identical N into different multiples of N and
                // conceal reuse at the next independent component.
                loop {
                    let previous_count = states.len();
                    let mut by_numerator = BTreeMap::<Atom, Atom>::new();
                    for (carrier, numerator) in std::mem::take(states) {
                        *by_numerator.entry(numerator).or_insert(Atom::Zero) += carrier;
                    }
                    let mut by_carrier = BTreeMap::<Atom, Atom>::new();
                    for (numerator, carrier) in by_numerator {
                        if !carrier.is_zero() {
                            *by_carrier.entry(carrier).or_insert(Atom::Zero) += numerator;
                        }
                    }
                    *states = by_carrier
                        .into_iter()
                        .filter(|(_, numerator)| !numerator.is_zero())
                        .collect();
                    // Summing numerators can create an equality which the
                    // first pass could not see. Every further useful pass
                    // strictly decreases the number of states.
                    if states.len() == previous_count {
                        break;
                    }
                }
            }

            composition_time += composition_started.elapsed();
            debug_tags!(#generation, #uv, #local, #four_d, #profile;
                stage = "component_composition",
                component,
                elapsed_ms = composition_time.as_secs_f64() * 1000.0,
                remaining_requests = terms.len(),
                remaining_states = terms.iter().map(|(_, _, _, _, states)| states.len()).sum::<usize>(),
                mapped_subtree_hits = context.numerator_rows.hits,
                mapped_subtree_misses = context.numerator_rows.misses,
                mapped_subtree_evictions = context.numerator_rows.evictions,
                retained_subtree_bytes = context.numerator_rows.retained_bytes(),
                "Composed and grouped independent UV component states"
            );
            context.numerator_rows.clear();
            debug_tags!(#generation, #uv, #local, #four_d, #cff, #summary;
                component,
                cached_exact_cff_expressions = context.generation_cache.len(),
                cached_preparations = context.preparations.len(),
                numerator_template_builds = context.numerator_template_builds,
                numerator_row_cache_hits = context.numerator_rows.hits,
                preparation_bytes = context.preparations.retained_bytes(),
                "Cached exact CFF expressions by topology and capacity in one local-4D component wave"
            );
        }

        let composition_started = Instant::now();
        let active = terms
            .into_iter()
            .flat_map(|(_, _, _, _, states)| states)
            .fold(Atom::Zero, |sum, (carrier, numerator)| {
                sum + carrier * numerator
            });
        debug_tags!(#generation, #uv, #local, #four_d, #profile;
            stage = "sector_composition",
            elapsed_ms = composition_started.elapsed().as_secs_f64() * 1000.0,
            "Composed factorized UV sector coefficient"
        );
        Ok(active)
    }
}

impl Projected4dApproximation<'_> {
    /// Project an already Taylor-expanded local 4D counterterm in the UV
    /// child's own energy frame. The contracted cograph is deliberately absent
    /// here: final assembly attaches it after the independent 4D Taylor
    /// operation has been converted to CFF.
    pub(crate) fn project_local_4d(
        &mut self,
        local: &Local4dCts,
        context: &mut Local4dProjectionContext,
    ) -> Result<Projected4dCts> {
        let projection_started = Instant::now();
        let preparation_cache_before = context.preparations.cache_time;
        let row_cache_before = context.numerator_rows.cache_time;
        if !self.settings.local_uv_cts_from_expanded_4d_integrands {
            return Err(eyre!(
                "the typed local-4D child projection is reserved for local counterterms requested from expanded 4D integrands"
            ));
        }
        if !self.localizer.orientation.explicit_orientation_sum_only {
            return Err(eyre!(
                "factorized local-4D-derived UV counterterms currently require `explicit_orientation_sum_only = true`; source-local child residues cannot be broadcast into ordinary orientation-local sectors"
            ));
        }

        // Select the typed local sectors directly. Recursive completions are
        // owned by the separately integrated branch; no cograph is attached here.
        if local.active_sectors().is_empty() {
            // Local4dCts prunes algebraically zero sectors and can also contain
            // only completed recursive factors. Its local contribution is then
            // a typed zero, which needs no child or production energy map.
            return Ok(Projected4dCts::new(vec![Projected4dSector {
                coefficient: Atom::Zero,
                frozen_factor: Atom::one(),
            }]));
        }
        if self.localizer.orientation.orientation_ids()?.is_empty() {
            return Err(eyre!(
                "orientation pattern selects no production energy maps"
            ));
        }

        let mut active_sectors = Vec::new();
        let normalization_started = Instant::now();
        let sectors = local.projection_sectors();
        debug_tags!(#generation, #uv, #local, #four_d, #profile;
            stage = "sector_grouping",
            raw_sectors = local.active_sectors().len(),
            grouped_sectors = sectors.len(),
            elapsed_ms = normalization_started.elapsed().as_secs_f64() * 1000.0,
            "Grouped compatible UV projection sectors"
        );
        for sector in sectors {
            let normalization_started = Instant::now();
            let cache_time_before = context.preparations.cache_time;
            let sector = context.canonical_projection(sector, self.graph)?;
            let cache_time = context.preparations.cache_time - cache_time_before;
            debug_tags!(#generation, #uv, #local, #four_d, #profile;
                stage = "canonical_normalization",
                buckets = sector.terms.len(),
                denominator_classes = sector.classes.len(),
                components = sector.active_components.len(),
                elapsed_ms = normalization_started.elapsed().as_secs_f64() * 1000.0,
                cache_work_ms = cache_time.as_secs_f64() * 1000.0,
                "Certified canonical factorized UV algebra"
            );
            let frozen_localizer = sector.frozen_lmbs.iter().fold(Atom::one(), |product, lmb| {
                product * GS.localizing_integrand(lmb)
            });
            let active = self
                .localizer
                .project_factorized_taylor_sector(self.graph, &sector, context)?;

            active_sectors.push(Projected4dSector {
                coefficient: active,
                frozen_factor: frozen_localizer,
            });
        }

        let (cff_hits, cff_misses, cff_evictions, retained_cff_bytes) =
            context.generation_cache.statistics();
        debug_tags!(#generation, #uv, #local, #four_d, #profile;
            stage = "local_projection_total",
            graph = %self.graph.name,
            elapsed_ms = projection_started.elapsed().as_secs_f64() * 1000.0,
            component_requests = context.projected_component_requests,
            template_builds = context.numerator_template_builds,
            canonical_preparation_builds = context.canonical_preparation_builds,
            source_preparation_builds = context.source_preparation_builds,
            preparation_hits = context.preparations.hits,
            preparation_misses = context.preparations.misses,
            preparation_evictions = context.preparations.evictions,
            retained_preparation_bytes = context.preparations.retained_bytes(),
            preparation_cache_ms = context.preparations.cache_time.as_secs_f64() * 1000.0,
            preparation_cache_work_ms = (context.preparations.cache_time - preparation_cache_before).as_secs_f64() * 1000.0,
            cff_hits,
            cff_misses,
            cff_evictions,
            retained_cff_bytes,
            native_generations = context.generation_cache.native_generations,
            row_hits = context.numerator_rows.hits,
            row_misses = context.numerator_rows.misses,
            row_evictions = context.numerator_rows.evictions,
            retained_row_bytes = context.numerator_rows.retained_bytes(),
            row_cache_ms = context.numerator_rows.cache_time.as_secs_f64() * 1000.0,
            row_cache_work_ms = (context.numerator_rows.cache_time - row_cache_before).as_secs_f64() * 1000.0,
            "Completed local four-dimensional UV projection"
        );
        Ok(Projected4dCts::new(active_sectors))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::uv::approx::local_4d::FourDSector;
    use crate::{
        dot,
        graph::{LMBext, parse::IntoGraph},
        initialisation::test_initialise,
        settings::global::OrientationPattern,
        utils::W_,
        uv::{Spinney, UltravioletGraph, approx::OrientationProjection},
    };
    use linnet::half_edge::{
        involution::EdgeIndex,
        subgraph::{InternalSubGraph, SubSetOps},
    };
    use symbolica::atom::{AtomCore, FunctionBuilder};

    #[test]
    fn canonical_source_and_template_preparation_share_retention_invariant_budget() -> Result<()> {
        test_initialise()?;
        let mut graph: Graph = dot!(digraph preparation_retention {
            edge [num=1 mass=1]
            node [num=1]
            a -> b [id=0 lmb_id=0]
            a -> b [id=1]
        })?;
        let mut coefficient = GS.emr_mom(EdgeIndex(0), GS.cind(0)).pow(2) + Atom::one();
        for edge in [EdgeIndex(0), EdgeIndex(1)] {
            let momentum = FunctionBuilder::new(GS.emr_mom)
                .add_arg(usize::from(edge))
                .finish();
            let polynomial = (1..=3).fold(
                GS.emr_mom(edge, GS.cind(0)).pow(2) - Atom::one(),
                |polynomial, index| polynomial - GS.emr_mom(edge, GS.cind(index)).pow(2),
            );
            coefficient *= GS
                .den(usize::from(edge), momentum, Atom::one(), polynomial)
                .pow(-2);
        }
        let full = graph.full_filter();
        let raw = FourDSector::new(
            coefficient,
            vec![(full.clone(), full, graph.loop_momentum_basis.clone())],
            Vec::new(),
        );
        let cutset = CutSet::empty(graph.n_hedges());
        let pattern = OrientationPattern::default();
        let production = Default::default();
        let options = graph.denominator_only_cff_3d_expression_options();
        let localizer = Localizer::new(
            &cutset,
            OrientationProjection::exact(&production, &options, &pattern, false),
        );
        let mut context = Local4dProjectionContext::default();
        let mut expected = None;
        for phase in ["cold", "warm", "disabled", "evicted", "evicted_again"] {
            if phase == "disabled" {
                context.preparations = crate::cff::generation::GenerationCache::new(0, 0);
            } else if phase == "evicted" {
                context.preparations =
                    crate::cff::generation::GenerationCache::new(48 * 1024 * 1024, 1);
            }
            let before = (
                context.canonical_preparation_builds,
                context.source_preparation_builds,
                context.numerator_template_builds,
                context.generation_cache.native_generations,
            );
            let sector = context.canonical_projection(raw.clone(), &graph)?;
            let coefficient = localizer
                .project_factorized_taylor_sector(&mut graph, &sector, &mut context)?
                .collect_factors();
            assert_eq!(
                expected.get_or_insert_with(|| coefficient.clone()),
                &coefficient,
                "{phase}"
            );
            assert!(context.preparations.retained_bytes() <= 48 * 1024 * 1024);
            if phase == "warm" {
                assert_eq!(
                    before,
                    (
                        context.canonical_preparation_builds,
                        context.source_preparation_builds,
                        context.numerator_template_builds,
                        context.generation_cache.native_generations,
                    )
                );
                assert_eq!(context.preparations.len(), 3);
            } else {
                assert_eq!(context.canonical_preparation_builds, before.0 + 1);
                assert_eq!(context.source_preparation_builds, before.1 + 1);
                assert_eq!(context.numerator_template_builds, before.2 + 1);
            }
            if phase == "disabled" {
                assert_eq!(context.preparations.len(), 0);
                assert_eq!(context.preparations.retained_bytes(), 0);
            } else if phase.starts_with("evicted") {
                assert_eq!(context.preparations.len(), 1);
                assert!(context.preparations.evictions >= 2);
            }
        }
        Ok(())
    }

    #[test]
    fn nested_banana_quotient_powered_component_has_the_analytic_one_energy_sign() -> Result<()> {
        test_initialise()?;
        let mut graph: Graph = dot!(digraph nested_banana_quotient_sign {
            edge [num=1 mass=1]
            node [num=1]
            incoming [style=invis]
            outgoing [style=invis]

            incoming -> a [id=0]
            a -> b [id=1 lmb_id=0]
            a -> b [id=2 lmb_id=1]
            a -> b [id=3]
            b -> outgoing [id=4]
        })?;
        let inner_filter = graph
            .get_edge_subgraph(EdgeIndex(1))
            .union(&graph.get_edge_subgraph(EdgeIndex(2)));
        let outer_filter = inner_filter.union(&graph.get_edge_subgraph(EdgeIndex(3)));
        let inner_subgraph =
            InternalSubGraph::cleaned_filter_optimist(inner_filter, graph.as_ref());
        let outer_subgraph =
            InternalSubGraph::cleaned_filter_optimist(outer_filter.clone(), graph.as_ref());
        let quotient_lmb = graph.shrunken_sub_lmb(
            &outer_filter,
            &inner_subgraph,
            graph.dummy_stripped_external_flows_of(&outer_subgraph),
            None,
        )?;
        assert_eq!(
            quotient_lmb.loop_edges.iter().copied().collect::<Vec<_>>(),
            vec![EdgeIndex(3)]
        );

        let edge = EdgeIndex(3);
        let momentum = FunctionBuilder::new(GS.emr_mom)
            .add_arg(usize::from(edge))
            .finish();
        let mass_squared = Atom::var(GS.m_uv_expansion).pow(2);
        let energy_squared = (1..=3).fold(mass_squared.clone(), |sum, spatial_index| {
            sum + GS.emr_mom(edge, GS.cind(spatial_index)).pow(2)
        });
        let full_denominator = GS.emr_mom(edge, GS.cind(0)).pow(2) - &energy_squared;
        let denominator = GS.den(
            usize::from(edge),
            momentum.clone(),
            mass_squared.clone(),
            full_denominator,
        );
        let provenance = GS.uv_momentum_provenance_tag(
            Atom::num(usize::from(edge) as i64).as_view(),
            true,
            momentum.as_view(),
        );
        let tagged_momentum = FunctionBuilder::new(GS.emr_mom)
            .add_arg(provenance.clone())
            .finish();
        let tagged_energy_squared = (1..=3).fold(mass_squared.clone(), |sum, spatial_index| {
            sum + FunctionBuilder::new(GS.emr_mom)
                .add_arg(provenance.clone())
                .add_arg(GS.cind(spatial_index))
                .finish()
                .pow(2)
        });
        let tagged_full_denominator = FunctionBuilder::new(GS.emr_mom)
            .add_arg(provenance)
            .add_arg(GS.cind(0))
            .finish()
            .pow(2)
            - tagged_energy_squared;
        let tagged_denominator = GS.den(
            usize::from(edge),
            tagged_momentum,
            mass_squared.clone(),
            tagged_full_denominator,
        );
        let constant = Atom::one() - mass_squared;
        let component = vec![(
            graph.get_edge_subgraph(edge),
            outer_filter,
            quotient_lmb.clone(),
        )];
        let powered_sector = FourDSector::new(
            (&tagged_denominator + &constant) * denominator.pow(-2),
            component.clone(),
            Vec::new(),
        );
        let cancelled_sector = FourDSector::new(
            tagged_denominator * denominator.pow(-2),
            component.clone(),
            Vec::new(),
        );
        let dotted_sector = FourDSector::new(denominator.pow(-2), component.clone(), Vec::new());
        let one_pole_sector = FourDSector::new(denominator.pow(-1), component, Vec::new());
        let pattern = OrientationPattern::default();
        let cutset = CutSet::empty(graph.n_hedges());
        // This isolated typed-sector projection consumes options, not a
        // production selector set; the sector carries its own exact sources.
        let production = Default::default();
        let projection_options = graph.denominator_only_cff_3d_expression_options();
        let localizer = Localizer::new(
            &cutset,
            OrientationProjection::exact(&production, &projection_options, &pattern, false),
        );
        let powered = {
            let canonical = powered_sector.canonical_projection(&graph)?;
            localizer.project_factorized_taylor_sector(
                &mut graph,
                &canonical,
                &mut Local4dProjectionContext::default(),
            )
        }?;
        let cancelled = {
            let canonical = cancelled_sector.canonical_projection(&graph)?;
            localizer.project_factorized_taylor_sector(
                &mut graph,
                &canonical,
                &mut Local4dProjectionContext::default(),
            )
        }?;
        let dotted = {
            let canonical = dotted_sector.canonical_projection(&graph)?;
            localizer.project_factorized_taylor_sector(
                &mut graph,
                &canonical,
                &mut Local4dProjectionContext::default(),
            )
        }?;
        let one_pole = {
            let canonical = one_pole_sector.canonical_projection(&graph)?;
            localizer.project_factorized_taylor_sector(
                &mut graph,
                &canonical,
                &mut Local4dProjectionContext::default(),
            )
        }?;

        // D/D^2 = 1/D, while the repeated-pole remainder obeys
        // CFF[1/D^2] = -CFF[1/D]/(2 E^2). Compare each term directly to the
        // absolute one-pole normalization.
        let unwrap_denominator = |atom: Atom| {
            GS.erase_uv_momentum_provenance(
                &atom.replace(GS.den(W_.a_, W_.b_, W_.c_, W_.d_)).with(W_.d_),
            )
        };
        let powered = unwrap_denominator(powered);
        let cancelled = unwrap_denominator(cancelled);
        let dotted = unwrap_denominator(dotted);
        let one_pole_projection = unwrap_denominator(one_pole);
        for (mass, spatial, energy) in [
            (
                Atom::num(2),
                [Atom::Zero, Atom::Zero, Atom::Zero],
                Atom::num(2),
            ),
            (
                Atom::num(symbolica::domains::rational::Rational::from((3, 5))),
                [
                    Atom::num(symbolica::domains::rational::Rational::from((4, 5))),
                    Atom::Zero,
                    Atom::Zero,
                ],
                Atom::one(),
            ),
        ] {
            let fixed_point = |mut atom: Atom| {
                atom = atom.replace(GS.m_uv_expansion).with(mass.clone());
                for (spatial_index, value) in spatial.iter().enumerate() {
                    atom = atom
                        .replace(GS.emr_mom(edge, GS.cind(spatial_index + 1)))
                        .with(value.clone());
                }
                atom.together()
            };
            let powered = fixed_point(powered.clone());
            let cancelled = fixed_point(cancelled.clone());
            let dotted = fixed_point(dotted.clone());
            let one_pole_projection = fixed_point(one_pole_projection.clone());
            let energy_squared = energy.pow(2);
            let constant = Atom::one() - mass.pow(2);
            let one_pole = -Atom::i() / (Atom::num(16) * Atom::var(GS.pi).pow(3) * energy);
            let expected_powered = &one_pole * (Atom::num(2) * &energy_squared - &constant)
                / (Atom::num(2) * &energy_squared);
            let difference = (&powered - &expected_powered).together();
            let opposite = (&powered + expected_powered).together();
            let cancellation = (&cancelled - &one_pole).together();
            let dotted_difference =
                (&dotted + &one_pole / (Atom::num(2) * energy_squared)).together();
            let linearity = (&cancelled + constant * &dotted - &powered).together();
            let cancellation_to_generated_one_pole = (&cancelled - &one_pole_projection).together();
            assert!(
                [
                    &difference,
                    &cancellation,
                    &dotted_difference,
                    &linearity,
                    &cancellation_to_generated_one_pole,
                ]
                .into_iter()
                .all(Atom::is_zero),
                "the isolated powered quotient violates its analytic one-energy form at mass={mass}, spatial={spatial:?}: powered difference={difference}, opposite-sign diagnostic={opposite}, D/D^2 difference={cancellation}, dotted difference={dotted_difference}, linearity difference={linearity}, generated one-pole={one_pole_projection}, D/D^2-to-generated-one-pole difference={cancellation_to_generated_one_pole}"
            );
        }
        Ok(())
    }

    #[test]
    fn typed_taylor_wave_batches_genuine_owner_relabelled_terms() -> Result<()> {
        test_initialise()?;
        let mut graph: Graph = dot!(
            digraph G {
                edge [particle="scalar_1"];
                node [num=1];
                a -> b [id=0];
                a -> b [id=1];
            },
            "scalars"
        )?;
        let full = graph.full_filter();
        let momentum = |edge: EdgeIndex| {
            FunctionBuilder::new(GS.emr_mom)
                .add_arg(usize::from(edge))
                .finish()
        };
        let denominator = |source_edge: EdgeIndex, momentum_edge: EdgeIndex| {
            GS.den(
                usize::from(source_edge),
                momentum(momentum_edge),
                graph.underlying[momentum_edge].particle.mass_atom().pow(2),
                Atom::one(),
            )
            .pow(-1)
        };
        let term_a = Atom::var(symbolica::symbol!("typed_taylor_cache::a"))
            * denominator(EdgeIndex(0), EdgeIndex(0))
            * denominator(EdgeIndex(1), EdgeIndex(1));
        let term_b = Atom::var(symbolica::symbol!("typed_taylor_cache::b"))
            * denominator(EdgeIndex(0), EdgeIndex(0))
            * denominator(EdgeIndex(0), EdgeIndex(1));
        let sector = FourDSector::new(
            &term_a + &term_b,
            vec![(
                full.clone(),
                full.clone(),
                graph.loop_momentum_basis.clone(),
            )],
            Vec::new(),
        );
        let pattern = OrientationPattern::default();
        let cutset = CutSet::empty(graph.n_hedges());
        // This isolated typed-sector projection consumes options, not a
        // production selector set; the sector carries its own exact sources.
        let production = Default::default();
        let projection_options = graph.denominator_only_cff_3d_expression_options();
        let localizer = Localizer::new(
            &cutset,
            OrientationProjection::exact(&production, &projection_options, &pattern, false),
        );
        let options = graph.denominator_only_cff_3d_expression_options();
        let uv_edges = graph
            .iter_edges_of(&full)
            .filter_map(|(pair, edge, data)| {
                (pair.is_paired() && !data.data.is_dummy).then_some(edge)
            })
            .collect::<Vec<_>>();
        let mut structural_cache = Local4dProjectionContext::default();
        for term in sector.physical_terms()? {
            graph.cff_from_4d_denominators_in_uv_sub_lmb(
                &term.denominators,
                uv_edges.iter().copied(),
                [],
                &graph.loop_momentum_basis.clone(),
                ExactUvSubLmbFrame::TaylorVacuum,
                &cutset,
                &options,
                &term.numerator,
                Some(&mut structural_cache),
                &[],
            )?;
        }

        let batched = {
            let canonical = sector.canonical_projection(&graph)?;
            localizer.project_factorized_taylor_sector(
                &mut graph,
                &canonical,
                &mut Local4dProjectionContext::default(),
            )
        }?;
        let mut sequential = Atom::Zero;
        for term in sector.physical_terms()? {
            let (cff, _) = graph.cff_from_4d_denominators_in_uv_sub_lmb(
                &term.denominators,
                uv_edges.iter().copied(),
                [],
                &graph.loop_momentum_basis.clone(),
                ExactUvSubLmbFrame::TaylorVacuum,
                &cutset,
                &options,
                &term.numerator,
                None,
                &[],
            )?;
            let production_prefactor = Atom::num(cff.production_prefactor_factor());
            for (index, cff_term) in cff.terms {
                assert_eq!(index, CutCFFIndex::new_all_none());
                for orientation in &cff_term.orientations {
                    sequential += &orientation.expression
                        * &production_prefactor
                        * cff_term.map_exact_source_numerator(&orientation.orientation, None)?;
                }
            }
        }
        assert!(
            (batched.collect_factors() - sequential.collect_factors()).is_zero(),
            "batched production projection must equal independently generated sequential CFFs"
        );

        let permuted = FourDSector::new(
            term_b + term_a,
            sector.active_components.clone(),
            Vec::new(),
        );
        let permuted = {
            let canonical = permuted.canonical_projection(&graph)?;
            localizer.project_factorized_taylor_sector(
                &mut graph,
                &canonical,
                &mut Local4dProjectionContext::default(),
            )
        }?;
        assert_eq!(
            permuted, batched,
            "component-wave registration and output must be invariant under term permutation"
        );
        Ok(())
    }

    #[test]
    fn typed_taylor_next_component_reuses_one_topology_for_prior_residue_states() -> Result<()> {
        test_initialise()?;
        let mut graph: Graph = dot!(digraph typed_taylor_component_waves {
            edge [particle="scalar_1"];
            node [num=1];
            a -> b [id=0 lmb_id=0];
            a -> b [id=1 lmb_id=1];
            a -> b [id=2];
        }, "scalars")?;
        let inner_filter = graph
            .get_edge_subgraph(EdgeIndex(0))
            .union(&graph.get_edge_subgraph(EdgeIndex(1)));
        let inner_subgraph =
            InternalSubGraph::cleaned_filter_optimist(inner_filter, graph.as_ref());
        let full = graph.full_filter();
        let outer_subgraph =
            InternalSubGraph::cleaned_filter_optimist(full.clone(), graph.as_ref());
        let inner_spinney =
            Spinney::new(inner_subgraph.clone(), &graph, &graph.loop_momentum_basis)
                .ok_or_else(|| eyre!("the inner component fixture has no compatible sub-LMB"))?;
        let shell_lmb = graph.shrunken_sub_lmb(
            &full,
            &inner_subgraph,
            graph.dummy_stripped_external_flows_of(&outer_subgraph),
            None,
        )?;
        let shell = full.subtract(&inner_subgraph.filter);
        assert_eq!(
            graph.n_loops(&shell),
            0,
            "deleting the inner bubble leaves a tree in the original incidence"
        );
        assert_eq!(
            shell_lmb.loop_edges.iter().copied().collect::<Vec<_>>(),
            vec![EdgeIndex(2)],
            "contracting the inner bubble turns that same shell into one quotient loop"
        );
        let denominator = |edge: EdgeIndex| {
            let momentum_edge = if edge == EdgeIndex(1) {
                EdgeIndex(0)
            } else {
                edge
            };
            let momentum = FunctionBuilder::new(GS.emr_mom)
                .add_arg(usize::from(momentum_edge))
                .finish();
            let momentum = if edge == EdgeIndex(1) {
                -momentum
            } else {
                momentum
            };
            GS.den(
                usize::from(edge),
                &momentum,
                graph.underlying[edge].particle.mass_atom().pow(2),
                Atom::one(),
            )
            .pow(-1)
        };
        let numerator = Atom::var(symbolica::symbol!("typed_taylor_cache::nested"));
        let coefficient = &numerator
            * denominator(EdgeIndex(0))
            * denominator(EdgeIndex(1))
            * denominator(EdgeIndex(2));
        let raised_coefficient = &coefficient * denominator(EdgeIndex(0));
        let sector = FourDSector::new(
            coefficient.clone(),
            vec![
                (
                    inner_subgraph.filter.clone(),
                    inner_subgraph.filter.clone(),
                    inner_spinney.lmb.clone(),
                ),
                (
                    graph.get_edge_subgraph(EdgeIndex(2)),
                    full.clone(),
                    shell_lmb.clone(),
                ),
            ],
            Vec::new(),
        );

        let pattern = OrientationPattern::default();
        let cutset = CutSet::empty(graph.n_hedges());
        // This isolated typed-sector projection consumes options, not a
        // production selector set; the sector carries its own exact sources.
        let production = Default::default();
        let projection_options = graph.denominator_only_cff_3d_expression_options();
        let localizer = Localizer::new(
            &cutset,
            OrientationProjection::exact(&production, &projection_options, &pattern, false),
        );
        let batched = {
            let canonical = sector.canonical_projection(&graph)?;
            localizer.project_factorized_taylor_sector(
                &mut graph,
                &canonical,
                &mut Local4dProjectionContext::default(),
            )
        }?;
        let options = graph.denominator_only_cff_3d_expression_options();
        let [term] = sector
            .physical_terms()?
            .try_into()
            .map_err(|terms: Vec<_>| {
                eyre!(
                    "the component-wave fixture produced {} physical terms",
                    terms.len()
                )
            })?;
        let (shell_denominators, inner_denominators): (Vec<_>, Vec<_>) = term
            .denominators
            .into_iter()
            .partition(|denominator| denominator.source_edge == EdgeIndex(2));
        let inner_edges = [EdgeIndex(0), EdgeIndex(1)];
        let mut states = Vec::new();
        let (inner_cff, _) = graph.cff_from_4d_denominators_in_uv_sub_lmb(
            &inner_denominators,
            inner_edges,
            [],
            &inner_spinney.lmb,
            ExactUvSubLmbFrame::TaylorVacuum,
            &cutset,
            &options,
            &term.numerator,
            None,
            &[],
        )?;
        let inner_prefactor = Atom::num(inner_cff.production_prefactor_factor());
        for (index, cff_term) in inner_cff.terms {
            assert_eq!(index, CutCFFIndex::new_all_none());
            for orientation in &cff_term.orientations {
                states.push((
                    &orientation.expression * &inner_prefactor,
                    cff_term.map_exact_source_numerator(&orientation.orientation, None)?,
                ));
            }
        }
        let uv_edges = [EdgeIndex(0), EdgeIndex(1), EdgeIndex(2)];
        let mut shell_cache = Local4dProjectionContext::default();

        let mut cached = Atom::Zero;
        let mut sequential = Atom::Zero;
        for (carrier, state_numerator) in states {
            for (cache, sum) in [
                (Some(&mut shell_cache), &mut cached),
                (None, &mut sequential),
            ] {
                let (cff, _) = graph.cff_from_4d_denominators_in_uv_sub_lmb(
                    &shell_denominators,
                    uv_edges,
                    [],
                    &shell_lmb,
                    ExactUvSubLmbFrame::TaylorVacuum,
                    &cutset,
                    &options,
                    &state_numerator,
                    cache,
                    &[],
                )?;
                let prefactor = Atom::num(cff.production_prefactor_factor());
                for (index, cff_term) in cff.terms {
                    assert_eq!(index, CutCFFIndex::new_all_none());
                    for orientation in &cff_term.orientations {
                        *sum += &carrier
                            * &orientation.expression
                            * &prefactor
                            * cff_term
                                .map_exact_source_numerator(&orientation.orientation, None)?;
                    }
                }
            }
        }
        assert!((cached.collect_factors() - sequential.collect_factors()).is_zero());
        assert!(
            (batched.collect_factors() - sequential.collect_factors()).is_zero(),
            "the production two-pass waves must equal fully uncached sequential component projection"
        );

        // Different first-component powers become the same remaining request.
        // Both terms retain one independent inner contour and one outer frame;
        // merging their completed states must avoid a duplicate outer request.
        let combined = FourDSector::new(
            &coefficient + &raised_coefficient,
            sector.active_components.clone(),
            Vec::new(),
        )
        .canonical_projection(&graph)?;
        assert_eq!(combined.terms.len(), 2);
        let mut combined_context = Local4dProjectionContext::default();
        let combined = localizer.project_factorized_taylor_sector(
            &mut graph,
            &combined,
            &mut combined_context,
        )?;
        let mut separate_context = Local4dProjectionContext::default();
        let mut separate = Atom::Zero;
        for coefficient in [coefficient, raised_coefficient] {
            let canonical =
                FourDSector::new(coefficient, sector.active_components.clone(), Vec::new())
                    .canonical_projection(&graph)?;
            separate += localizer.project_factorized_taylor_sector(
                &mut graph,
                &canonical,
                &mut separate_context,
            )?;
        }
        assert!(
            (combined.collect_factors() - separate.collect_factors())
                .together()
                .is_zero()
        );
        assert!(
            combined_context.projected_component_requests
                < separate_context.projected_component_requests
        );
        Ok(())
    }
}
