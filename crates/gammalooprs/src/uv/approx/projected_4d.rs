use std::{
    collections::BTreeMap,
    ops::Neg,
    sync::Arc,
    time::{Duration, Instant},
};

use color_eyre::Result;
use eyre::eyre;
use linnet::half_edge::subgraph::SubSetLike;
use spenso::structure::slot::IsAbstractSlot;
use symbolica::{
    atom::{Atom, AtomCore, Indeterminate},
    id::Replacement,
    symbol,
};
use three_dimensional_reps::CffGenerationContext;

use crate::{
    cff::{
        CutCFFIndex,
        generation::{PreparationKey, PreparationValue},
    },
    debug_tags,
    graph::{ExactUvSubLmbFrame, Graph, cuts::CutSet},
    integrands::process::param_builder::FnMapEntry,
    numerator::energy_degree::{EnergyPowerAnalyzer, EnergyPowerCapMap},
    utils::GS,
    uv::{
        Integrands, UVgenerationSettings,
        approx::{
            direct_3d::DirectResidueBranches,
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
    pub(crate) coefficient: Integrands,
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
                        coefficient: sector.coefficient.fallible_map(&mut f)?.map_numerators(
                            |tensor| Ok(tensor.with_checked_expression(f(tensor.expression())?)?),
                        )?,
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
    #[cfg(test)]
    pub(crate) certify_denominator_reconstruction: bool,
    #[cfg(test)]
    pub(crate) denominator_reconstruction_certificates: usize,
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
            #[cfg(test)]
            certify_denominator_reconstruction: false,
            #[cfg(test)]
            denominator_reconstruction_certificates: 0,
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
    ) -> Result<Integrands> {
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
        let mut options = self.orientation.cff_options()?.clone();
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
                        .analyze_atom(term.source_numerator.expression())
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
                term.numerator.clone(),
                term.on_shell_points.clone(),
                Vec::<Atom>::new(),
                vec![(residual_factor, Vec::<Atom>::new())],
            ));
        }
        // Each component still owns its exact source and occurrence assignment.
        // Previous residues bind scalar arguments; they never specialize the
        // ordinary numerator used to prepare the next source. A zero leading
        // coefficient in one row can therefore overestimate that row's rank,
        // but cannot understate any expression subsequently sampled.
        for (component, (_, source_scope, coordinate_lmb)) in
            sector.active_components.iter().enumerate()
        {
            let mut composition_time = Duration::ZERO;
            // A coefficient position has the same formal name throughout this
            // wave. Different affine meanings remain visible in the ordinary
            // body and in each ordered argument row before families can merge.
            let scope = DirectResidueBranches::numerator_scope().1;
            let uv_edges = graph
                .iter_edges_of(source_scope)
                .filter_map(|(pair, edge_id, edge)| {
                    (pair.is_paired() && !edge.data.is_dummy).then_some(edge_id)
                })
                .collect::<Vec<_>>();
            // Incidence uses the enclosing source scope, while denominator
            // owners retain this component's quotient topology and capacities.
            for (
                component_denominators,
                classes,
                _,
                raw_reports,
                numerator,
                on_shell_points,
                parameters,
                rows,
            ) in &mut terms
            {
                if rows.is_empty() {
                    continue;
                }
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
                    log.numerator = numerator.expression(),
                    "Preparing shared factorized local-4D child numerator"
                );
                context.projected_component_requests += 1;
                let fixed_external_momenta = on_shell_points
                    .iter()
                    .filter(|point| point.component == component)
                    .map(|point| point.binding())
                    .collect::<BTreeMap<_, _>>();
                let boundaries = if !fixed_external_momenta.is_empty() {
                    graph
                        .dummy_stripped_external_flows_of(source_scope)
                        .included_iter()
                        .collect::<Vec<_>>()
                } else {
                    Vec::new()
                };
                let frame = if !fixed_external_momenta.is_empty() {
                    ExactUvSubLmbFrame::RetainedPhysicalCrown
                } else {
                    ExactUvSubLmbFrame::TaylorVacuum
                };
                let (mut cff, _) = graph.cff_from_4d_denominators_in_uv_sub_lmb(
                    denominators,
                    uv_edges.iter().copied(),
                    boundaries,
                    coordinate_lmb,
                    frame,
                    &child_cutset,
                    &options,
                    numerator.expression(),
                    Some(&mut *context),
                    &active_classes,
                    &fixed_external_momenta,
                )?;
                cff.energy_degree_bound_report.physical_parent_bounds =
                    raw_reports[component].clone().into_generation_bounds();
                self.orientation
                    .record_energy_degree_bound_report(&cff.energy_degree_bound_report);
                // The sector already carries its forest subtraction sign. Only
                // the child CFF energy-factor convention changes at this bridge.
                let production_prefactor = Atom::num(cff.production_prefactor_factor());
                if cff.terms.is_empty()
                    || cff.terms.values().all(|term| term.orientations.is_empty())
                {
                    rows.clear();
                    continue;
                }
                if cff.terms.len() != 1 {
                    return Err(eyre!("an uncut UV-child CFF must have one residue index"));
                }
                let (index, cff_term) = cff.terms.into_iter().next().unwrap();
                if index != CutCFFIndex::new_all_none() {
                    return Err(eyre!(
                        "an uncut UV-child CFF unexpectedly produced residue index {index}"
                    ));
                }
                let (mapped, new_parameters, argument_rows) = cff_term
                    .prepare_exact_source_numerator(&scope)
                    .map_err(|error| {
                        eyre!(
                            "{error}; exact UV-child component denominators are {:?}",
                            denominators
                        )
                    })?;
                let composition_started = Instant::now();
                let mut next_rows = BTreeMap::<Vec<Atom>, Atom>::new();
                for (carrier, arguments) in std::mem::take(rows) {
                    for (orientation, incoming) in cff_term.orientations.iter().zip(&argument_rows)
                    {
                        let arguments = arguments
                            .iter()
                            .cloned()
                            .chain(incoming.iter().cloned())
                            .collect();
                        *next_rows.entry(arguments).or_insert(Atom::Zero) +=
                            &carrier * &orientation.expression * &production_prefactor;
                    }
                }
                *rows = next_rows
                    .into_iter()
                    .filter_map(|(arguments, carrier)| {
                        (!carrier.is_zero()).then_some((carrier, arguments))
                    })
                    .collect();
                *numerator = numerator.with_checked_expression(mapped)?;
                on_shell_points.retain(|point| point.component != component);
                parameters.extend(new_parameters);
                // Equal complete argument rows share one summed carrier. Other
                // Taylor families remain separate ordinary bodies; adding them
                // must not create another numerator copy per residue sample.
                composition_time += composition_started.elapsed();
            }
            let composition_started = Instant::now();
            let families_before = terms.len();
            let mut remaining = BTreeMap::new();
            for (
                mut denominators,
                mut classes,
                mut powers,
                mut reports,
                numerator,
                on_shell_points,
                parameters,
                rows,
            ) in terms
            {
                if rows.is_empty() {
                    continue;
                }
                // Completed powers no longer distinguish future source
                // requests. Equal remaining class IDs retain the same exact
                // frame and denominator witness; never concatenate witnesses.
                denominators[component].clear();
                classes.retain(|class| class.component > component);
                powers.retain(|class, _| sector.class(*class).component > component);
                reports[component] = Default::default();
                // Equal scalar spellings are insufficient for ranked zeros or
                // differently ordered open ports. Keep the complete interface.
                let interface = numerator
                    .structure()
                    .logical_slots()
                    .into_iter()
                    .map(|slot| slot.to_atom())
                    .collect::<Vec<_>>();
                let key = (
                    powers.clone(),
                    numerator.expression().clone(),
                    interface,
                    on_shell_points.clone(),
                    parameters.clone(),
                );
                let (_, _, _, combined_reports, _, _, _, combined_rows) =
                    remaining.entry(key).or_insert_with(|| {
                        (
                            denominators,
                            classes,
                            powers,
                            vec![EnergyPowerCapMap::default(); reports.len()],
                            numerator,
                            on_shell_points,
                            parameters,
                            BTreeMap::<Vec<Atom>, Atom>::new(),
                        )
                    });
                for (combined, report) in combined_reports.iter_mut().zip(reports) {
                    combined.max_assign(report);
                }
                // Merge only exact body/formal equality. Distinct argument
                // rows, including maps with the same physical host, survive.
                for (carrier, arguments) in rows {
                    *combined_rows.entry(arguments).or_insert(Atom::Zero) += carrier;
                }
            }
            terms = remaining
                .into_values()
                .map(
                    |(
                        denominators,
                        classes,
                        powers,
                        reports,
                        numerator,
                        on_shell_points,
                        parameters,
                        rows,
                    )| {
                        (
                            denominators,
                            classes,
                            powers,
                            reports,
                            numerator,
                            on_shell_points,
                            parameters,
                            rows.into_iter()
                                .filter_map(|(arguments, carrier)| {
                                    (!carrier.is_zero()).then_some((carrier, arguments))
                                })
                                .collect(),
                        )
                    },
                )
                .collect();
            composition_time += composition_started.elapsed();
            debug_tags!(#generation, #uv, #local, #four_d, #profile;
                stage = "component_composition",
                component,
                elapsed_ms = composition_time.as_secs_f64() * 1000.0,
                families_before,
                remaining_requests = terms.len(),
                remaining_states = terms.iter().map(|(_, _, _, _, _, _, _, rows)| rows.len()).sum::<usize>(),
                "Composed scalar argument rows for shared UV numerator families"
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
        let mut active = Atom::Zero;
        let mut definitions = Vec::new();
        for (_, _, _, _, numerator, _, parameters, rows) in terms {
            if rows.is_empty() || numerator.is_zero() {
                continue;
            }
            let mut used = vec![false; parameters.len()];
            let _ = numerator.expression().replace_map(|view, _, _| {
                for (index, parameter) in parameters.iter().enumerate() {
                    used[index] |= parameter.as_view() == view;
                }
            });
            let parameters = parameters
                .into_iter()
                .zip(&used)
                .filter_map(|(parameter, used)| used.then_some(parameter))
                .collect::<Vec<_>>();
            let scope = DirectResidueBranches::numerator_scope().1;
            let lhs = symbol!("gammalooprs::uv::numerator_family")
                .call_args(std::iter::once(scope.clone()).chain(parameters.iter().cloned()));
            let entry = Arc::new(FnMapEntry {
                lhs,
                rhs: numerator,
                args: parameters
                    .iter()
                    .cloned()
                    .map(Indeterminate::try_from)
                    .collect::<std::result::Result<Vec<_>, _>>()
                    .map_err(|error| eyre!(error))?,
                tags: vec![scope],
                inlining: symbolica::evaluate::InliningPolicy::Always,
                is_alias: false,
            });
            for (carrier, arguments) in rows {
                let arguments = arguments
                    .into_iter()
                    .zip(&used)
                    .filter_map(|(argument, used)| used.then_some(argument));
                let call = entry.lhs.replace_multiple(
                    parameters
                        .iter()
                        .zip(arguments)
                        .map(|(parameter, value)| Replacement::new(parameter.to_pattern(), value)),
                );
                active += carrier * call;
            }
            definitions.push(entry);
        }
        debug_tags!(#generation, #uv, #local, #four_d, #profile;
            stage = "sector_composition",
            numerator_families = definitions.len(),
            elapsed_ms = composition_started.elapsed().as_secs_f64() * 1000.0,
            "Composed factorized UV sector coefficient"
        );
        Integrands::from_iter([(CutCFFIndex::new_all_none(), active)]).with_numerators(definitions)
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
        let normalization_started = Instant::now();
        let sectors = local.projection_sectors()?;
        debug_tags!(#generation, #uv, #local, #four_d, #profile;
            stage = "sector_grouping",
            raw_sectors = local.active_sectors().len(),
            grouped_sectors = sectors.len(),
            elapsed_ms = normalization_started.elapsed().as_secs_f64() * 1000.0,
            "Grouped compatible UV projection sectors"
        );
        if sectors.is_empty() {
            // Local4dCts prunes algebraically zero sectors and can also contain
            // only completed recursive factors. Grouping can cancel nonzero
            // sectors with identical bindings. The local contribution is then
            // a typed zero, which needs no child or production energy map.
            return Ok(Projected4dCts::new(vec![Projected4dSector {
                coefficient: Integrands::from_iter([(CutCFFIndex::new_all_none(), Atom::Zero)]),
                frozen_factor: Atom::one(),
            }]));
        }
        if self.localizer.orientation.orientation_ids()?.is_empty() {
            return Err(eyre!(
                "orientation pattern selects no production energy maps"
            ));
        }

        let mut active_sectors = Vec::new();
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
    use crate::uv::approx::{
        Rooted, UVCtx,
        local_4d::{FourDSector, Full4dCts, on_shell_forest_tests, uv_limit},
    };
    use crate::{
        finalized_runtime_dot,
        graph::LMBext,
        initialisation::test_initialise,
        settings::global::OrientationPattern,
        utils::W_,
        uv::{
            Spinney, UltravioletGraph, approx::OrientationProjection, hedge_poset::OwnedForestNode,
        },
    };
    use linnet::half_edge::{
        involution::EdgeIndex,
        subgraph::{InternalSubGraph, SuBitGraph, SubSetOps},
    };
    use symbolica::{
        atom::{AtomCore, FunctionBuilder},
        domains::rational::Rational,
        function,
    };

    impl Local4dProjectionContext {
        fn assert_physical_family(
            &mut self,
            graph: &mut Graph,
            local: &Local4dCts,
            components: usize,
        ) -> Result<()> {
            let before = local.clone();
            let cutset = CutSet::empty(graph.n_hedges());
            let pattern = OrientationPattern::default();
            let production = Default::default();
            let options = graph.denominator_only_cff_3d_expression_options();
            let localizer = Localizer::new(
                &cutset,
                OrientationProjection::exact(&production, &options, &pattern, false),
            );
            let mut nonzero = false;
            let requests_before = self.projected_component_requests;
            for sector in local.projection_sectors()? {
                // Both canonical reconstruction and the immutable CFF plan's
                // signed diagonal certificate are production boundaries.
                let canonical = sector.canonical_projection(graph)?;
                assert_eq!(canonical.active_components.len(), components);
                let projected =
                    localizer.project_factorized_taylor_sector(graph, &canonical, self)?;
                let has_coefficient = projected
                    .iter()
                    .any(|(_, coefficient)| !coefficient.is_zero());
                nonzero |= has_coefficient;
                if has_coefficient {
                    assert!(
                        !projected.numerators().is_empty(),
                        "physical forest coefficients must retain shared numerator bodies"
                    );
                }
                // Finite tensor execution is deliberately not performed here.
                // Spinor factors must remain in their shared bodies while the
                // scalar CFF rows carry only the residue weights and calls.
                for (_, coefficient) in projected.iter() {
                    for tensor in [
                        idenso::dirac::AGS.gamma,
                        idenso::dirac::AGS.gamma0,
                        spenso::network::tags::SPENSO_TAG.trace,
                    ] {
                        assert!(
                            !coefficient.contains_symbol(tensor),
                            "a physical graph numerator was copied into the scalar CFF rows"
                        );
                    }
                }
            }
            assert!(
                nonzero,
                "the projected physical forest family must not vanish structurally"
            );
            assert!(self.projected_component_requests >= requests_before + components);
            assert_eq!(
                local, &before,
                "projection must preserve the raw completed child sectors"
            );
            Ok(())
        }

        fn assert_nested_physical_families(
            &mut self,
            fixture: (Graph, OwnedForestNode, OwnedForestNode, OwnedForestNode),
        ) -> Result<()> {
            let (mut graph, parent, child, root) = fixture;
            let settings = UVgenerationSettings {
                generate_integrated: false,
                ..Default::default()
            };
            let (child_local, parent_local, nested_local) = {
                let ctx = UVCtx::new(&graph, &graph.model, &settings);
                let seed = Full4dCts::from_factorized_local(&Local4dCts::root());
                let child_local = uv_limit(&seed, &ctx, &child, &root, &child, &root)?;
                let parent_local = uv_limit(&seed, &ctx, &parent, &root, &parent, &root)?;
                let completed = Full4dCts::from_factorized_local(&child_local);
                let nested_local = uv_limit(&completed, &ctx, &parent, &child, &parent, &child)?;
                (child_local, parent_local, nested_local)
            };
            for (family, local, components) in [
                ("child", &child_local, 1),
                ("parent", &parent_local, 1),
                ("completed child inside parent", &nested_local, 2),
            ] {
                self.assert_physical_family(&mut graph, local, components)
                    .map_err(|error| eyre!("{} {family}: {error:#}", graph.name))?;
            }
            Ok(())
        }
    }

    #[test]
    fn os_quark_in_ir_gluon_physical_families_project_to_cff() -> Result<()> {
        test_initialise()?;
        Local4dProjectionContext::default()
            .assert_nested_physical_families(on_shell_forest_tests::nested_gluon_fixture()?)
    }

    #[test]
    fn os_quark_in_os_quark_physical_families_project_to_cff() -> Result<()> {
        test_initialise()?;
        Local4dProjectionContext::default()
            .assert_nested_physical_families(on_shell_forest_tests::nested_quark_fixture(false)?)
    }

    #[test]
    fn ir_gluon_in_os_quark_physical_families_project_to_cff() -> Result<()> {
        test_initialise()?;
        Local4dProjectionContext::default()
            .assert_nested_physical_families(on_shell_forest_tests::nested_quark_fixture(true)?)
    }

    #[test]
    fn disjoint_os_quark_physical_families_project_to_cff() -> Result<()> {
        test_initialise()?;
        let (mut graph, left, right, root) = on_shell_forest_tests::disjoint_quark_fixture()?;
        let settings = UVgenerationSettings {
            generate_integrated: false,
            ..Default::default()
        };
        let (left_local, right_local, product) = {
            let ctx = UVCtx::new(&graph, &graph.model, &settings);
            let seed = Full4dCts::from_factorized_local(&Local4dCts::root());
            let left_local = uv_limit(&seed, &ctx, &left, &root, &left, &root)?;
            let right_local = uv_limit(&seed, &ctx, &right, &root, &right, &root)?;
            let product = Local4dCts::from_full_product([
                Full4dCts::from_factorized_local(&left_local),
                Full4dCts::from_factorized_local(&right_local),
            ])?;
            (left_local, right_local, product)
        };
        let mut context = Local4dProjectionContext::default();
        for (family, local, components) in [
            ("left child", &left_local, 1),
            ("right child", &right_local, 1),
            ("disjoint product", &product, 2),
        ] {
            context
                .assert_physical_family(&mut graph, local, components)
                .map_err(|error| eyre!("{} {family}: {error:#}", graph.name))?;
        }
        Ok(())
    }

    #[test]
    fn on_shell_point_keeps_a_live_crown_carrier_in_the_same_numerator() -> Result<()> {
        use crate::utils::symbols::UvMomentumProvenanceRole;
        test_initialise()?;
        let model = crate::utils::load_generic_model("sm");
        let mut graph: Graph = finalized_runtime_dot!(digraph os_live_crown {
            graph [projector=1]
            edge [num=1 mass=1]
            node [num=1]
            incoming [style=invis]
            outgoing [style=invis]
            incoming -> a [id=0 sink="{ufo_order:0}"]
            a -> b [id=1 lmb_id=0 source="{ufo_order:1}" sink="{ufo_order:0}"]
            a -> b [id=2 source="{ufo_order:2}" sink="{ufo_order:1}"]
            b -> outgoing [id=3 source="{ufo_order:2}"]
        }, &model)?;
        let owners = graph
            .full_filter()
            .subtract(&graph.external_filter::<SuBitGraph>());
        let subgraph = InternalSubGraph::cleaned_filter_optimist(owners.clone(), graph.as_ref());
        let lmb = graph.try_compatible_sub_lmb(
            &subgraph,
            graph.dummy_stripped_external_flows_of(&subgraph),
            &graph.loop_momentum_basis,
        )?;
        let external = *lmb.ext_edges.iter().next().unwrap();
        let mass = Atom::var(symbol!("os_live_crown::mass"));
        let loop_momentum = function!(GS.emr_mom, 1);
        let live_tag = GS.uv_momentum_provenance_tag(
            usize::from(external),
            UvMomentumProvenanceRole::PhysicalSourceFixed,
            function!(GS.emr_mom, usize::from(external)),
        );
        let quadratic = |momentum: &Atom, mass_squared: Atom| {
            (0..=3).fold(-mass_squared, |sum, index| {
                sum + Atom::num(if index == 0 { 1 } else { -1 })
                    * GS.indexed_momentum(momentum, &[GS.cind(index)]).pow(2)
            })
        };
        let cutset = CutSet::empty(graph.n_hedges());
        let pattern = OrientationPattern::default();
        let production = Default::default();
        let options = graph.denominator_only_cff_3d_expression_options();
        let localizer = Localizer::new(
            &cutset,
            OrientationProjection::exact(&production, &options, &pattern, false),
        );
        let mut context = Local4dProjectionContext::default();
        for velocity in [["1", "0", "0", "0"], ["5/4", "1/4", "1/2", "1/2"]] {
            let velocity = UVgenerationSettings {
                on_shell_subtraction_velocity: velocity.map(str::to_owned),
                ..Default::default()
            }
            .on_shell_velocity()?;
            for sign in [1, -1] {
                let point = GS.on_shell_momentum(external, &(Atom::num(sign) * &mass), &velocity);
                let shifted = &point - &loop_momentum;
                let quark = GS.den(
                    1,
                    &loop_momentum,
                    mass.pow(2),
                    quadratic(&loop_momentum, mass.pow(2)),
                );
                let denominator = &quark * GS.den(2, &shifted, 0, quadratic(&shifted, Atom::Zero));
                // The binding certificate preserves the complete pole polynomial,
                // including explicit prescription data. This is a reconstruction
                // check only; it does not enable additional contour prescriptions.
                let prescription = Atom::var(symbol!("os_live_crown::prescription"));
                let prescribed = [-1, 1]
                    .into_iter()
                    .map(|sign| {
                        Atom::one()
                            / (&quark
                                * GS.den(
                                    2,
                                    &shifted,
                                    0,
                                    quadratic(&shifted, Atom::Zero)
                                        + Atom::num(sign) * Atom::i() * &prescription,
                                ))
                    })
                    .sum::<Atom>();
                let certificate = FourDSector::new(
                    crate::graph::parse::NumeratorInterface::admit(prescribed, None)?
                        .try_into_indexed()?,
                    vec![(owners.clone(), owners.clone(), lmb.clone())],
                    Vec::new(),
                )
                .canonical_projection(&graph)?;
                assert_eq!(
                    certificate
                        .classes
                        .iter()
                        .filter(|class| class.members.contains(&EdgeIndex(2)))
                        .count(),
                    2,
                    "opposite explicit prescriptions must remain distinct denominator classes"
                );
                let invalid = &quark
                    * GS.den(
                        2,
                        &shifted,
                        0,
                        quadratic(&shifted, Atom::Zero) + GS.emr_mom(external, GS.cind(1)),
                    );
                let error = FourDSector::new(
                    crate::graph::parse::NumeratorInterface::admit(Atom::one() / invalid, None)?
                        .try_into_indexed()?,
                    vec![(owners.clone(), owners.clone(), lmb.clone())],
                    Vec::new(),
                )
                .canonical_projection(&graph)
                .expect_err(
                    "a live crown in the source polynomial cannot masquerade as the fixed point",
                );
                assert!(
                    error
                        .to_string()
                        .contains("OS denominator point-binding certificate failed")
                );
                for index in [0, 1] {
                    let live = function!(GS.emr_mom, &live_tag, GS.cind(index));
                    for use_shifted_carrier in [false, true] {
                        let source = if use_shifted_carrier {
                            function!(
                                GS.emr_mom,
                                GS.uv_momentum_provenance_tag(
                                    2,
                                    UvMomentumProvenanceRole::TaylorFixed,
                                    &shifted
                                ),
                                GS.cind(0)
                            )
                        } else {
                            GS.indexed_momentum(&point, &[GS.cind(0)])
                        };
                        let numerator = source * &live;
                        let sector = FourDSector::new(
                            crate::graph::parse::NumeratorInterface::admit(
                                &numerator / &denominator,
                                None,
                            )?
                            .try_into_indexed()?,
                            vec![(owners.clone(), owners.clone(), lmb.clone())],
                            Vec::new(),
                        );
                        let canonical = sector.canonical_projection(&graph)?;
                        assert_eq!(canonical.terms.len(), 1);
                        let term = &canonical.terms[0];
                        assert_eq!(
                            crate::graph::parse::NumeratorInterface::export(&term.source_numerator)
                                .0,
                            numerator,
                            "the raw factorized witness must retain both momentum domains"
                        );
                        // Factor sigma using sigma^2=1. This fixes numeric signs
                        // independently of the production routing/certificate.
                        let expected = Atom::num(sign)
                            * (&mass * &velocity[0]
                                - Atom::num(sign)
                                    * if use_shifted_carrier {
                                        GS.emr_mom(EdgeIndex(1), GS.cind(0))
                                    } else {
                                        Atom::Zero
                                    })
                            * GS.emr_mom(external, GS.cind(index));
                        for numerator in [&term.source_numerator, &term.numerator] {
                            let neutral = canonical.neutral_numerator(
                                numerator,
                                &graph,
                                &term.on_shell_points,
                            )?;
                            assert!(
                                crate::uv::overlap_numerators::NumeratorJetCache::default()
                                    .canonicalize(&(neutral.expression() - &expected))?
                                    .is_zero,
                                "the fixed point and live crown must agree exactly, allowing equivalent factored signs"
                            );
                        }
                        let projected = localizer.project_factorized_taylor_sector(
                            &mut graph,
                            &canonical,
                            &mut context,
                        )?;
                        assert!(!projected.numerators().is_empty());
                        for entry in projected.numerators() {
                            assert!(
                                entry
                                    .rhs
                                    .expression()
                                    .pattern_match(&live.to_pattern(), None, None)
                                    .next()
                                    .is_some(),
                                "source mapping must retain the typed live parent factor"
                            );
                            if !use_shifted_carrier {
                                assert_eq!(
                                    entry.rhs.expression(),
                                    &(Atom::num(sign) * &mass * &velocity[0] * &live),
                                    "source point binding must precede any lowering of the live parent carrier"
                                );
                            }
                        }
                    }
                }
            }
        }
        Ok(())
    }

    #[test]
    fn on_shell_point_projection_matches_independent_bubble_residues() -> Result<()> {
        test_initialise()?;
        let model = crate::utils::load_generic_model("sm");
        let mut graph: Graph = finalized_runtime_dot!(digraph on_shell_bubble {
            graph [projector=1]
            edge [num=1 mass=1]
            node [num=1]
            incoming [style=invis]
            outgoing [style=invis]
            incoming -> a [id=0 sink="{ufo_order:0}"]
            a -> b [id=1 lmb_id=0 source="{ufo_order:1}" sink="{ufo_order:0}"]
            a -> b [id=2 source="{ufo_order:2}" sink="{ufo_order:1}"]
            b -> outgoing [id=3 source="{ufo_order:2}"]
        }, &model)?;
        let owners = graph
            .full_filter()
            .subtract(&graph.external_filter::<SuBitGraph>());
        let subgraph = InternalSubGraph::cleaned_filter_optimist(owners.clone(), graph.as_ref());
        let lmb = graph.try_compatible_sub_lmb(
            &subgraph,
            graph.dummy_stripped_external_flows_of(&subgraph),
            &graph.loop_momentum_basis,
        )?;
        let external = *lmb.ext_edges.iter().next().unwrap();
        let mass = Atom::var(symbol!("os_projection_test::mass"));
        let scalar_parameter = symbol!("os_projection_test::a");
        let scalar = Atom::var(scalar_parameter);
        let spectator = Atom::var(symbol!("os_projection_test::b"))
            + Atom::var(symbol!("os_projection_test::c"));
        let loop_momentum = function!(GS.emr_mom, 1);
        let loop_energy = GS.emr_mom(EdgeIndex(1), GS.cind(0));
        let spatial_squared = (1..=3)
            .map(|index| GS.emr_mom(EdgeIndex(1), GS.cind(index)).pow(2))
            .sum::<Atom>();
        let gluon_energy = Atom::var(symbol!("os_projection_test::E_gluon"));
        let cutset = CutSet::empty(graph.n_hedges());
        let pattern = OrientationPattern::default();
        let production = Default::default();
        let options = graph.denominator_only_cff_3d_expression_options();
        let localizer = Localizer::new(
            &cutset,
            OrientationProjection::exact(&production, &options, &pattern, false),
        );
        let mut context = Local4dProjectionContext::default();
        for (request, sign) in [-1, 1, 0, -1].into_iter().enumerate() {
            let branch_mass = if sign == 0 { Atom::Zero } else { mass.clone() };
            let fermion_invariant = &spatial_squared + branch_mass.pow(2);
            let fermion_ose = fermion_invariant.sqrt();
            let fermion_energy = if sign == 0 {
                gluon_energy.clone()
            } else {
                Atom::var(symbol!("os_projection_test::E_fermion"))
            };
            let point_energy = Atom::num(sign) * &branch_mass;
            let point = GS.on_shell_momentum(
                external,
                &point_energy,
                &UVgenerationSettings::default().on_shell_velocity()?,
            );
            let gluon_momentum = &point - &loop_momentum;
            let quadratic = |momentum: &Atom, mass_squared: Atom| {
                (0..=3).fold(-mass_squared, |sum, index| {
                    sum + Atom::num(if index == 0 { 1 } else { -1 })
                        * GS.indexed_momentum(momentum, &[GS.cind(index)]).pow(2)
                })
            };
            let fermion_denominator = GS.den(
                1,
                &loop_momentum,
                branch_mass.pow(2),
                quadratic(&loop_momentum, branch_mass.pow(2)),
            );
            let gluon_denominator = GS.den(
                2,
                &gluon_momentum,
                0,
                quadratic(&gluon_momentum, Atom::Zero),
            );
            // The point remains indexed until projection. The independent
            // numerator factors are never combined with denominator algebra.
            let numerator = (&loop_energy + &scalar) * &spectator;
            let sector = FourDSector::new(
                crate::graph::parse::NumeratorInterface::admit(
                    numerator / fermion_denominator / gluon_denominator,
                    None,
                )?
                .try_into_indexed()?,
                vec![(owners.clone(), owners.clone(), lmb.clone())],
                Vec::new(),
            );
            let canonical = sector.canonical_projection(&graph)?;
            assert!(
                canonical
                    .terms
                    .iter()
                    .all(|term| !term.on_shell_points.is_empty())
            );
            for term in &canonical.terms {
                assert_eq!(term.source_witness.len(), 2);
                assert_eq!(term.powers.values().sum::<usize>(), 2);
                assert_eq!(
                    term.source_witness
                        .iter()
                        .map(|denominator| denominator.mass_squared.clone())
                        .collect::<std::collections::BTreeSet<_>>(),
                    std::collections::BTreeSet::from([Atom::Zero, branch_mass.pow(2)]),
                );
            }
            let projected = localizer
                .project_factorized_taylor_sector(&mut graph, &canonical, &mut context)?
                .resolved()?;
            assert_eq!(
                context.source_preparation_builds,
                (request + 1).min(3),
                "different fixed points need distinct prepared mappers; repeating a point reuses its exact binding"
            );
            let actual = projected
                .iter()
                .next()
                .unwrap()
                .1
                .clone()
                .replace(crate::utils::ose_atom_from_index(EdgeIndex(1)))
                .with(fermion_ose.clone())
                // Both exact positive source energies are compared in a
                // neutral formal frame. Replace every half-integer invariant
                // power, including inverse roots in the residue measure.
                .replace_map(|view, _, output| {
                    let symbolica::atom::AtomView::Pow(power) = view else {
                        return;
                    };
                    let energy = if power.get_base() == fermion_invariant.as_view() {
                        &fermion_energy
                    } else if power.get_base() == spatial_squared.as_view() {
                        &gluon_energy
                    } else {
                        return;
                    };
                    let twice_exponent = power.get_exp().to_owned() * Atom::num(2);
                    if let Ok(exponent) = i64::try_from(twice_exponent.as_view()) {
                        **output = energy.pow(exponent);
                    }
                });
            assert!(!actual.contains_symbol(GS.os_mom));
            assert!(!actual.contains_symbol(GS.emr_mom));
            // Clockwise Feynman contour: the poles below are k0=E_m and
            // k0=p0+E_0. The production bridge includes i/(2 pi)^3.
            let expected = if sign == 0 {
                // At m=pOS=0 there is one double pole. Differentiating
                // (z+a)/(z+E)^2 at z=E gives -a/(4*E^3).
                Atom::i() * &scalar * &spectator
                    / (Atom::num(32) * Atom::var(GS.pi).pow(3) * gluon_energy.pow(3))
            } else {
                let residue = (&fermion_energy + &scalar)
                    / (Atom::num(2)
                        * &fermion_energy
                        * ((&point_energy - &fermion_energy).pow(2) - gluon_energy.pow(2)))
                    + (&point_energy + &gluon_energy + &scalar)
                        / (Atom::num(2)
                            * &gluon_energy
                            * ((&point_energy + &gluon_energy).pow(2) - fermion_energy.pow(2)));
                -Atom::i() * residue * &spectator / (Atom::num(8) * Atom::var(GS.pi).pow(3))
            };
            // Certify the two exact coefficients of the degree-one numerator
            // separately. Rational cancellation then acts on scalar contour
            // weights only, never on the complete graph numerator sum.
            let actual = actual / &spectator;
            let expected = expected / &spectator;
            for value in [&actual, &expected] {
                assert!(
                    value
                        .derivative(scalar_parameter)
                        .derivative(scalar_parameter)
                        .is_zero()
                );
            }
            for coefficient in [1, 0] {
                let extract = |value: &Atom| {
                    if coefficient == 0 {
                        value.replace(scalar_parameter).with(Atom::Zero)
                    } else {
                        value.derivative(scalar_parameter)
                    }
                };
                let actual_coefficient = extract(&actual);
                let expected_coefficient = extract(&expected);
                let difference = (&actual_coefficient - &expected_coefficient).together();
                assert!(
                    difference.is_zero(),
                    "OS branch {sign}, numerator coefficient {coefficient}: actual={actual_coefficient}, expected={expected_coefficient}, difference={difference}"
                );
                if sign != 0 {
                    // On the positive soft ray k=rho*r, choose m=|r|=1.
                    // The independent residue has weights proportional to
                    // 1/(E*w^2) and sigma/(w^2*(E+w)). The radial measure
                    // rho^2 therefore leaves finite, nonzero limits 1 and
                    // sigma. Apply this check to the certified production
                    // coefficients, without expanding any graph numerator.
                    let rho = symbol!("os_projection_test::soft_rho");
                    let radial = Atom::var(rho);
                    let weighted = (&actual_coefficient * Atom::num(32) * Atom::var(GS.pi).pow(3)
                        / Atom::i()
                        * radial.pow(2))
                    .replace(fermion_energy.clone())
                    .with((Atom::one() + radial.pow(2)).sqrt())
                    .replace(gluon_energy.clone())
                    .with(radial)
                    .replace(mass.clone())
                    .with(Atom::one());
                    let soft_series = weighted.series(rho, Atom::Zero, 0)?;
                    assert_eq!(
                        soft_series.get_trailing_exponent(),
                        Rational::from(0),
                        "OS branch {sign}, coefficient {coefficient} must have exactly integrable rho^-2 soft scaling"
                    );
                    assert_eq!(
                        soft_series.coefficient(Rational::from(0)),
                        Some(Atom::num(if coefficient == 1 { 1 } else { sign })),
                        "the massive OS soft-limit check must retain its nonzero leading term"
                    );
                }
            }
        }
        Ok(())
    }

    #[test]
    fn boosted_doubled_bubble_projection_matches_independent_soft_residues() -> Result<()> {
        test_initialise()?;
        let model = crate::utils::load_generic_model("sm");
        let mut graph: Graph = finalized_runtime_dot!(digraph doubled_on_shell_bubble {
            graph [projector=1]
            edge [num=1 mass=1]
            node [num=1]
            incoming [style=invis]
            outgoing [style=invis]
            incoming -> a [id=0 sink="{ufo_order:0}"]
            a -> b [id=1 lmb_id=0 source="{ufo_order:1}" sink="{ufo_order:0}"]
            a -> b [id=2 source="{ufo_order:2}" sink="{ufo_order:1}"]
            b -> outgoing [id=3 source="{ufo_order:2}"]
        }, &model)?;
        let owners = graph
            .full_filter()
            .subtract(&graph.external_filter::<SuBitGraph>());
        let subgraph = InternalSubGraph::cleaned_filter_optimist(owners.clone(), graph.as_ref());
        let lmb = graph.try_compatible_sub_lmb(
            &subgraph,
            graph.dummy_stripped_external_flows_of(&subgraph),
            &graph.loop_momentum_basis,
        )?;
        let external = *lmb.ext_edges.iter().next().unwrap();
        let mass = Atom::num(12);
        let constant_parameter = symbol!("os_double_bubble_test::constant");
        let linear_parameter = symbol!("os_double_bubble_test::linear");
        let constant = Atom::var(constant_parameter);
        let linear = Atom::var(linear_parameter);
        let spectator_head = symbol!("os_double_bubble_test::spectator"; Scalar);
        let spectator = function!(
            spectator_head,
            (Atom::var(symbol!("os_double_bubble_test::u"))
                + Atom::var(symbol!("os_double_bubble_test::v")))
            .pow(11)
        );
        let loop_momentum = function!(GS.emr_mom, 1);
        let loop_energy = GS.emr_mom(EdgeIndex(1), GS.cind(0));
        let fermion_invariant = mass.pow(2)
            + (1..=3)
                .map(|index| GS.emr_mom(EdgeIndex(1), GS.cind(index)).pow(2))
                .sum::<Atom>();
        let fermion_energy = Atom::var(symbol!("os_double_bubble_test::E"));
        let gluon_energy = Atom::var(symbol!("os_double_bubble_test::W"));
        let energy_variable = symbol!("os_double_bubble_test::z");
        let energy = Atom::var(energy_variable);
        let radius_variable = symbol!("os_double_bubble_test::rho");
        let radius = Atom::var(radius_variable);
        let cutset = CutSet::empty(graph.n_hedges());
        let pattern = OrientationPattern::default();
        let production = Default::default();
        let options = graph.denominator_only_cff_3d_expression_options();
        let localizer = Localizer::new(
            &cutset,
            OrientationProjection::exact(&production, &options, &pattern, false),
        );
        let mut context = Local4dProjectionContext::default();
        for (velocity_components, positive_spatial_point) in [
            (["5/4", "3/4", "0", "0"], [9, 0, 0]),
            (["5/4", "1/4", "1/2", "1/2"], [3, 6, 6]),
        ] {
            let velocity = UVgenerationSettings {
                on_shell_subtraction_velocity: velocity_components.map(str::to_owned),
                ..Default::default()
            }
            .on_shell_velocity()?;
            for sign in [-1, 1] {
                let point_energy = Atom::num(sign * 15);
                let spatial_point = positive_spatial_point.map(|component| sign * component);
                let point = GS.on_shell_momentum(external, &Atom::num(sign * 12), &velocity);
                let gluon_momentum = &point - &loop_momentum;
                let quadratic = |momentum: &Atom, mass_squared: Atom| {
                    (0..=3).fold(-mass_squared, |sum, index| {
                        sum + Atom::num(if index == 0 { 1 } else { -1 })
                            * GS.indexed_momentum(momentum, &[GS.cind(index)]).pow(2)
                    })
                };
                let fermion_denominator = GS.den(
                    1,
                    &loop_momentum,
                    mass.pow(2),
                    quadratic(&loop_momentum, mass.pow(2)),
                );
                let gluon_denominator = GS.den(
                    2,
                    &gluon_momentum,
                    0,
                    quadratic(&gluon_momentum, Atom::Zero),
                );
                let gluon_invariant = (1..=3)
                    .map(|index| {
                        (GS.emr_mom(EdgeIndex(1), GS.cind(index))
                            - Atom::num(spatial_point[index - 1]))
                        .pow(2)
                    })
                    .sum::<Atom>();
                let numerator = (&constant + &linear * (&loop_energy - &point_energy)) * &spectator;
                let source = crate::graph::parse::NumeratorInterface::admit(
                    &numerator / fermion_denominator.pow(2) / gluon_denominator,
                    None,
                )?
                .try_into_indexed()?;
                assert!(source.is_scalar());
                let sector = FourDSector::new(
                    source,
                    vec![(owners.clone(), owners.clone(), lmb.clone())],
                    Vec::new(),
                );
                // The existing canonical owner certifies signed routing, complete
                // pole polynomials, powers and source incidence before CFF dispatch.
                let canonical = sector.canonical_projection(&graph)?;
                assert_eq!(canonical.terms.len(), 1);
                let term = &canonical.terms[0];
                assert_eq!(
                    crate::graph::parse::NumeratorInterface::export(&term.source_numerator).0,
                    numerator
                );
                assert!(!term.on_shell_points.is_empty());
                assert_eq!(term.powers.values().sum::<usize>(), 3);
                assert_eq!(term.powers.values().copied().max(), Some(2));
                assert_eq!(term.source_witness.len(), 3);
                for (edge, multiplicity, mass_squared) in [
                    (EdgeIndex(1), 2, mass.pow(2)),
                    (EdgeIndex(2), 1, Atom::Zero),
                ] {
                    let occurrences = term
                        .source_witness
                        .iter()
                        .filter(|denominator| denominator.source_edge == edge)
                        .collect::<Vec<_>>();
                    assert_eq!(occurrences.len(), multiplicity);
                    assert!(
                        occurrences
                            .iter()
                            .all(|denominator| denominator.mass_squared == mass_squared)
                    );
                }
                let projected = localizer.project_factorized_taylor_sector(
                    &mut graph,
                    &canonical,
                    &mut context,
                )?;
                assert!(
                    projected.numerators().iter().any(|entry| {
                        entry
                            .rhs
                            .expression()
                            .pattern_match(&spectator.to_pattern(), None, None)
                            .next()
                            .is_some()
                    }),
                    "the opaque powered spectator stays in a shared numerator body"
                );
                let resolved = projected.resolved()?;
                let actual = resolved
                    .iter()
                    .map(|(_, coefficient)| coefficient.clone())
                    .sum::<Atom>()
                    .replace(crate::utils::ose_atom_from_index(EdgeIndex(1)))
                    .with(fermion_invariant.sqrt())
                    .replace_map(|view, _, output| {
                        let symbolica::atom::AtomView::Pow(power) = view else {
                            return;
                        };
                        if !power.get_base().contains_symbol(GS.emr_mom) {
                            return;
                        }
                        // Expand only a scalar energy polynomial for equality;
                        // neither numerator nor spectator enters this comparison.
                        let base = power.get_base().to_owned();
                        let positive_energy = if (&base - &fermion_invariant).expand().is_zero() {
                            &fermion_energy
                        } else if (&base - &gluon_invariant).expand().is_zero() {
                            &gluon_energy
                        } else {
                            return;
                        };
                        let twice_exponent = power.get_exp().to_owned() * Atom::num(2);
                        if let Ok(exponent) = i64::try_from(twice_exponent.as_view()) {
                            **output = positive_energy.pow(exponent);
                        }
                    });
                assert!(!actual.contains_symbol(GS.os_mom));
                assert!(!actual.contains_symbol(GS.emr_mom));
                assert!(!actual.contains_symbol(GS.ose));
                assert!(actual.contains_symbol(spectator_head));

                // Independent clockwise Feynman residues: one double pole at E
                // and one simple pole at q0+W. Keep both until after their sum.
                let energy_numerator = &constant + &linear * (&energy - &point_energy);
                let double_pole = (&energy_numerator
                    / ((&energy + &fermion_energy).pow(2)
                        * ((&energy - &point_energy).pow(2) - gluon_energy.pow(2))))
                .derivative(energy_variable)
                .replace(energy_variable)
                .with(fermion_energy.clone());
                let gluon_pole = energy_numerator
                    .replace(energy_variable)
                    .with(&point_energy + &gluon_energy)
                    / (Atom::num(2)
                        * &gluon_energy
                        * ((&point_energy + &gluon_energy).pow(2) - fermion_energy.pow(2)).pow(2));
                let contour_factor = Atom::i() / (Atom::num(8) * Atom::var(GS.pi).pow(3));
                let expected = -&contour_factor * (double_pole + gluon_pole);
                let actual = actual / &spectator;
                for value in [&actual, &expected] {
                    assert!(
                        value
                            .derivative(constant_parameter)
                            .derivative(constant_parameter)
                            .is_zero()
                    );
                    assert!(
                        value
                            .derivative(linear_parameter)
                            .derivative(linear_parameter)
                            .is_zero()
                    );
                    assert!(
                        value
                            .derivative(constant_parameter)
                            .derivative(linear_parameter)
                            .is_zero()
                    );
                    assert!(
                        value
                            .replace(constant_parameter)
                            .with(0)
                            .replace(linear_parameter)
                            .with(0)
                            .is_zero()
                    );
                }
                for (parameter, radial_power, leading_sign) in
                    [(constant_parameter, 3, -1), (linear_parameter, 2, sign)]
                {
                    let actual_coefficient = actual.derivative(parameter);
                    let expected_coefficient = expected.derivative(parameter);
                    assert!(
                        (&actual_coefficient - &expected_coefficient)
                            .together()
                            .is_zero(),
                        "velocity={velocity_components:?}, sign={sign}, parameter={parameter}: exact scalar CFF coefficient"
                    );
                    for ray in [[1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0]] {
                        let scalar_product: i32 =
                            spatial_point.iter().zip(ray).map(|(p, n)| p * n).sum();
                        let energy_squared = Atom::num(225)
                            + Atom::num(2 * scalar_product) * &radius
                            + radius.pow(2);
                        // Positive source energies stay exact until after the
                        // complete residue sum. These unit rays have W=rho.
                        let weighted = (&actual_coefficient / &contour_factor
                            * radius.pow(radial_power))
                        .replace(fermion_energy.clone())
                        .with(energy_squared.sqrt())
                        .replace(gluon_energy.clone())
                        .with(radius.clone());
                        let series = weighted.series(radius_variable, Atom::Zero, 0)?;
                        assert_eq!(series.get_trailing_exponent(), Rational::from(0));
                        assert_eq!(
                            series.coefficient(Rational::from(0)),
                            Some(
                                Atom::num(leading_sign)
                                    / Atom::num(8 * (15 + scalar_product).pow(2))
                            ),
                            "velocity={velocity_components:?}, sign={sign}, ray={ray:?}, parameter={parameter}: independent nonzero radial coefficient"
                        );
                    }
                }
            }
        }
        Ok(())
    }

    #[test]
    fn projected_rows_retain_one_ordinary_numerator_and_scalar_arguments() -> Result<()> {
        test_initialise()?;
        let mut graph: Graph = finalized_runtime_dot!(digraph projected_parametric_numerator { graph [projector=1]
            edge [num=1 mass=1]
            node [num=1]
            a -> b [id=0 lmb_id=0 source="{ufo_order:0}" sink="{ufo_order:0}"]
            a -> b [id=1 source="{ufo_order:1}" sink="{ufo_order:1}"]
        })?;
        let common = symbol!("projected_parametric_test::common"; Scalar);
        let numerator = Atom::var(common)
            * (GS.emr_mom(EdgeIndex(0), GS.cind(0)).pow(2) + Atom::one())
            * (GS.emr_mom(EdgeIndex(1), GS.cind(0)) + Atom::num(3)).pow(2);
        let mut coefficient = numerator.clone();
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
                .pow(-1);
        }
        let full = graph.full_filter();
        let sector = FourDSector::new(
            crate::numerator::AppliedFeynmanRule::new(coefficient)?.expr,
            vec![(full.clone(), full, graph.loop_momentum_basis.clone())],
            Vec::new(),
        );
        let canonical = sector.canonical_projection(&graph)?;
        let cutset = CutSet::empty(graph.n_hedges());
        let pattern = OrientationPattern::default();
        let production = Default::default();
        let options = graph.denominator_only_cff_3d_expression_options();
        let localizer = Localizer::new(
            &cutset,
            OrientationProjection::exact(&production, &options, &pattern, false),
        );
        let mut context = Local4dProjectionContext::default();
        let projected =
            localizer.project_factorized_taylor_sector(&mut graph, &canonical, &mut context)?;
        assert_eq!(context.projected_component_requests, 1);
        assert_eq!(projected.numerators().len(), 1);
        let entry = &projected.numerators()[0];
        assert!(entry.rhs.expression().contains_symbol(common));
        assert!(
            !entry
                .rhs
                .expression()
                .contains_symbol(symbol!("gammalooprs::uv::numerator_family"))
        );
        assert!(!entry.args.is_empty());
        let root = projected.iter().next().unwrap().1;
        assert!(
            !root.contains_symbol(common),
            "complete graph factors must live only in the shared body"
        );
        let opaque_numerator = symbol!("projected_parametric_test::mapped_numerator");
        let replacements = [entry.replacement()];
        let mut calls = 0;
        let actual = root.replace_map(|view, _, out| {
            if let symbolica::atom::AtomView::Fun(call) = view
                && call.get_symbol() == symbol!("gammalooprs::uv::numerator_family")
            {
                calls += 1;
                assert!(
                    call.iter()
                        .skip(entry.tags.len())
                        .all(|argument| { matches!(argument, symbolica::atom::AtomView::Num(_)) })
                );
                **out = symbolica::function!(
                    opaque_numerator,
                    view.replace_multiple(&replacements).collect_factors()
                );
            }
        });
        assert!(
            calls > 1,
            "the fixture must share its numerator across distinct residue samples"
        );
        // The fixture's two unit-mass lines have opposite loop momenta,
        // so both physical OSE aliases denote the same spatial energy.
        let energy_squared = (1..=3).fold(Atom::one(), |square, index| {
            square + GS.emr_mom(EdgeIndex(0), GS.cind(index)).pow(2)
        });
        let on_shell = [EdgeIndex(0), EdgeIndex(1)].map(|edge| {
            Replacement::new(
                crate::utils::ose_atom_from_index(edge).to_pattern(),
                energy_squared.sqrt().to_pattern(),
            )
        });
        let mut expected = Atom::Zero;
        for term in sector.physical_terms()? {
            let (cff, _) = graph.cff_from_4d_denominators(
                &term.denominators,
                &cutset,
                &options,
                term.numerator.expression(),
            )?;
            let prefactor = Atom::num(cff.production_prefactor_factor());
            for term in cff.terms.values() {
                for orientation in &term.orientations {
                    expected += (&orientation.expression * &prefactor).replace_multiple(&on_shell)
                        * symbolica::function!(
                            opaque_numerator,
                            term.map_exact_source_numerator(&orientation.orientation, None)?
                                .replace_multiple(&on_shell)
                                .collect_factors()
                        );
                }
            }
        }
        let resolved = projected.resolved()?;
        // Collect each side's complete numerator independently inside an opaque
        // function: expansion leaves its factorized argument untouched and only
        // compares scalar CFF weights after clearing their common denominator.
        let denominator = Atom::var(GS.pi).pow(3) * energy_squared.pow(2);
        assert!(
            ((actual - expected) * denominator).expand().is_zero(),
            "shared projected rows must equal independently mapped exact source numerators"
        );

        // The final family has a new tag, but its inherited coefficient scopes
        // remain actual atoms in the persisted formals and ordinary body. They
        // must not depend on unrelated, unused State reservations surviving.
        for argument in &entry.args {
            let argument = Atom::from(argument.clone());
            let symbolica::atom::AtomView::Fun(parameter) = argument.as_view() else {
                panic!("an exact-map formal must be a scoped coefficient");
            };
            assert!(matches!(
                parameter.get(0),
                symbolica::atom::AtomView::Var(_)
            ));
            assert_ne!(parameter.get(0), entry.tags[0].as_view());
        }
        let encoded = bincode::encode_to_vec(&projected, bincode::config::standard())?;
        let mut state = Vec::new();
        symbolica::state::State::export(&mut state)?;
        let state_map = symbolica::state::State::import(&mut std::io::Cursor::new(state), None)?;
        let model = crate::model::Model::empty("test");
        let (decoded, consumed): (Integrands, _) = bincode::decode_from_slice_with_context(
            &encoded,
            bincode::config::standard(),
            crate::GammaLoopContextContainer {
                state_map: &state_map,
                model: &model,
            },
        )?;
        assert_eq!(consumed, encoded.len());
        assert_eq!(decoded, projected);
        assert_eq!(decoded.resolved()?, resolved);
        Ok(())
    }

    #[test]
    fn canonical_source_and_template_preparation_share_retention_invariant_budget() -> Result<()> {
        test_initialise()?;
        let mut graph: Graph = finalized_runtime_dot!(digraph preparation_retention {
            projector=1;
            edge [num=1 mass=1]
            node [num=1]
            a -> b [id=0 lmb_id=0 source="{ufo_order:0}" sink="{ufo_order:0}"]
            a -> b [id=1 source="{ufo_order:1}" sink="{ufo_order:1}"]
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
            crate::numerator::AppliedFeynmanRule::new(coefficient)?.expr,
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
            let sector = context.canonical_projection(raw.clone(), &graph)?;
            let coefficient = localizer
                .project_factorized_taylor_sector(&mut graph, &sector, &mut context)
                .and_then(|coefficient| {
                    Ok(coefficient.resolved()?.iter().next().unwrap().1.clone())
                })?
                .collect_factors();
            assert_eq!(
                expected.get_or_insert_with(|| coefficient.clone()),
                &coefficient,
                "{phase}"
            );
            let budget = if phase == "disabled" {
                0
            } else {
                48 * 1024 * 1024
            };
            assert!(context.preparations.retained_bytes() <= budget);
        }
        Ok(())
    }

    #[test]
    fn nested_banana_quotient_powered_component_has_the_analytic_one_energy_sign() -> Result<()> {
        test_initialise()?;
        let mut graph: Graph = finalized_runtime_dot!(digraph nested_banana_quotient_sign { graph [projector=1]
            edge [num=1 mass=1]
            node [num=1]
            incoming [style=invis]
            outgoing [style=invis]

            incoming -> a [id=0 sink="{ufo_order:3}"]
            a -> b [id=1 lmb_id=0 source="{ufo_order:0}" sink="{ufo_order:0}"]
            a -> b [id=2 lmb_id=1 source="{ufo_order:1}" sink="{ufo_order:1}"]
            a -> b [id=3 source="{ufo_order:2}" sink="{ufo_order:2}"]
            b -> outgoing [id=4 source="{ufo_order:3}"]
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
            crate::numerator::AppliedFeynmanRule::new(
                (&tagged_denominator + &constant) * denominator.pow(-2),
            )?
            .expr,
            component.clone(),
            Vec::new(),
        );
        let cancelled_sector = FourDSector::new(
            crate::numerator::AppliedFeynmanRule::new(tagged_denominator * denominator.pow(-2))?
                .expr,
            component.clone(),
            Vec::new(),
        );
        let dotted_sector = FourDSector::new(
            crate::numerator::AppliedFeynmanRule::new(denominator.pow(-2))?.expr,
            component.clone(),
            Vec::new(),
        );
        let one_pole_sector = FourDSector::new(
            crate::numerator::AppliedFeynmanRule::new(denominator.pow(-1))?.expr,
            component,
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
        let powered = {
            let canonical = powered_sector.canonical_projection(&graph)?;
            localizer
                .project_factorized_taylor_sector(
                    &mut graph,
                    &canonical,
                    &mut Local4dProjectionContext::default(),
                )
                .and_then(|coefficient| {
                    Ok(coefficient.resolved()?.iter().next().unwrap().1.clone())
                })
        }?;
        let cancelled = {
            let canonical = cancelled_sector.canonical_projection(&graph)?;
            localizer
                .project_factorized_taylor_sector(
                    &mut graph,
                    &canonical,
                    &mut Local4dProjectionContext::default(),
                )
                .and_then(|coefficient| {
                    Ok(coefficient.resolved()?.iter().next().unwrap().1.clone())
                })
        }?;
        let dotted = {
            let canonical = dotted_sector.canonical_projection(&graph)?;
            localizer
                .project_factorized_taylor_sector(
                    &mut graph,
                    &canonical,
                    &mut Local4dProjectionContext::default(),
                )
                .and_then(|coefficient| {
                    Ok(coefficient.resolved()?.iter().next().unwrap().1.clone())
                })
        }?;
        let one_pole = {
            let canonical = one_pole_sector.canonical_projection(&graph)?;
            localizer
                .project_factorized_taylor_sector(
                    &mut graph,
                    &canonical,
                    &mut Local4dProjectionContext::default(),
                )
                .and_then(|coefficient| {
                    Ok(coefficient.resolved()?.iter().next().unwrap().1.clone())
                })
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
        let mut graph: Graph = finalized_runtime_dot!(
            digraph G { graph [projector=1]
                edge [particle="scalar_1" num="1𝑖"];
                node [num=1];
                a -> b [id=0 source="{ufo_order:0}" sink="{ufo_order:0}" lmb_id=0];
                a -> b [id=1 source="{ufo_order:1}" sink="{ufo_order:1}"];
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
                graph.underlying[momentum_edge]
                    .particle
                    .mass_atom(&graph.model)
                    .pow(2),
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
            crate::numerator::AppliedFeynmanRule::new(&term_a + &term_b)?.expr,
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
                term.numerator.expression(),
                Some(&mut structural_cache),
                &[],
                &BTreeMap::new(),
            )?;
        }

        let batched = {
            let canonical = sector.canonical_projection(&graph)?;
            localizer
                .project_factorized_taylor_sector(
                    &mut graph,
                    &canonical,
                    &mut Local4dProjectionContext::default(),
                )
                .and_then(|coefficient| {
                    Ok(coefficient.resolved()?.iter().next().unwrap().1.clone())
                })
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
                term.numerator.expression(),
                None,
                &[],
                &BTreeMap::new(),
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
            crate::numerator::AppliedFeynmanRule::new(term_b + term_a)?.expr,
            sector.active_components.clone(),
            Vec::new(),
        );
        let permuted = {
            let canonical = permuted.canonical_projection(&graph)?;
            localizer
                .project_factorized_taylor_sector(
                    &mut graph,
                    &canonical,
                    &mut Local4dProjectionContext::default(),
                )
                .and_then(|coefficient| {
                    Ok(coefficient.resolved()?.iter().next().unwrap().1.clone())
                })
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
        let mut graph: Graph = finalized_runtime_dot!(digraph typed_taylor_component_waves { graph [projector=1]
            edge [particle="scalar_1" num="1𝑖"];
            node [num=1];
            a -> b [id=0 lmb_id=0 source="{ufo_order:0}" sink="{ufo_order:0}"];
            a -> b [id=1 lmb_id=1 source="{ufo_order:1}" sink="{ufo_order:1}"];
            a -> b [id=2 source="{ufo_order:2}" sink="{ufo_order:2}"];
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
            Spinney::new(inner_subgraph.clone(), &graph, &graph.loop_momentum_basis)?
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
                graph.underlying[edge]
                    .particle
                    .mass_atom(&graph.model)
                    .pow(2),
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
            crate::numerator::AppliedFeynmanRule::new(coefficient.clone())?.expr,
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
            localizer
                .project_factorized_taylor_sector(
                    &mut graph,
                    &canonical,
                    &mut Local4dProjectionContext::default(),
                )
                .and_then(|coefficient| {
                    Ok(coefficient.resolved()?.iter().next().unwrap().1.clone())
                })
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
            term.numerator.expression(),
            None,
            &[],
            &BTreeMap::new(),
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
                    &BTreeMap::new(),
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
        // merging their completed states must preserve the complete outer value.
        let combined = FourDSector::new(
            crate::numerator::AppliedFeynmanRule::new(&coefficient + &raised_coefficient)?.expr,
            sector.active_components.clone(),
            Vec::new(),
        )
        .canonical_projection(&graph)?;
        let mut combined_context = Local4dProjectionContext::default();
        let combined = localizer.project_factorized_taylor_sector(
            &mut graph,
            &combined,
            &mut combined_context,
        )?;
        assert_eq!(
            combined.numerators().len(),
            1,
            "equal ordinary bodies must retain one coefficient family"
        );
        assert_eq!(
            combined_context.projected_component_requests, 3,
            "two first-component requests must share one remaining-component preparation"
        );
        let combined = combined.resolved()?.iter().next().unwrap().1.clone();
        let mut separate_context = Local4dProjectionContext::default();
        let mut separate = Atom::Zero;
        for coefficient in [coefficient, raised_coefficient] {
            let canonical = FourDSector::new(
                crate::numerator::AppliedFeynmanRule::new(coefficient)?.expr,
                sector.active_components.clone(),
                Vec::new(),
            )
            .canonical_projection(&graph)?;
            separate += localizer
                .project_factorized_taylor_sector(&mut graph, &canonical, &mut separate_context)
                .and_then(|coefficient| {
                    Ok(coefficient.resolved()?.iter().next().unwrap().1.clone())
                })?;
        }
        assert_eq!(separate_context.projected_component_requests, 4);
        assert!(
            (combined.collect_factors() - separate.collect_factors())
                .together()
                .is_zero()
        );
        Ok(())
    }
}
