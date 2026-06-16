use std::ops::Neg;

use color_eyre::Result;
use eyre::eyre;
use linnet::half_edge::subgraph::{Inclusion, SubSetLike};
use symbolica::atom::Atom;
use three_dimensional_reps::CffGenerationContext;

use crate::{
    cff::CutCFFIndex,
    debug_tags,
    graph::{ExactUvSubLmbFrame, Graph, cuts::CutSet},
    utils::GS,
    uv::{
        UVgenerationSettings,
        approx::{
            local_3d::Localizer,
            local_4d::{FourDSector, Local4dCts},
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
        sector: &FourDSector,
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
        let mut terms = Vec::new();
        for term in sector.physical_terms()? {
            let mut component_denominators = vec![Vec::new(); sector.active_components.len()];
            let mut residual_factor = Atom::one();
            for denominator in &term.denominators {
                let memberships = sector
                    .active_components
                    .iter()
                    .enumerate()
                    .filter(|(_, (owners, _, _))| {
                        owners.includes(&graph[&denominator.source_edge].1)
                    })
                    .map(|(component, _)| component)
                    .collect::<Vec<_>>();
                match memberships.as_slice() {
                    [] => residual_factor /= &denominator.full_expr,
                    [component] => {
                        // A denominator owned by this Taylor component can
                        // nevertheless be independent of its quotient loop
                        // coordinates after a pinch. It then has no contour
                        // pole and remains an ordinary factor, exactly as in
                        // the former single-frame projection.
                        if denominator
                            .momentum_signature_in_lmb(
                                &sector.active_components[*component].2,
                                true,
                            )?
                            .loop_signature
                            .iter()
                            .any(|coefficient| *coefficient != 0)
                        {
                            component_denominators[*component].push(denominator.clone());
                        } else {
                            residual_factor /= &denominator.full_expr;
                        }
                    }
                    _ => {
                        return Err(eyre!(
                            "4D denominator owner {} belongs to overlapping active Taylor components {:?}",
                            usize::from(denominator.source_edge),
                            memberships
                        ));
                    }
                }
            }

            terms.push((
                component_denominators,
                vec![(residual_factor, term.numerator)],
            ));
        }
        // Each Taylor component owns an independent energy contour. Every
        // state entering one component retains its own occurrence capacities.
        // Reuse requires equal topology and capacity, so no registration pass
        // or combined bound envelope is needed before generation.
        // The next component starts a fresh cache because its
        // numerators are immutable outputs of this wave's exact residue maps.
        for (component, (_, source_scope, coordinate_lmb)) in
            sector.active_components.iter().enumerate()
        {
            let uv_edges = graph
                .iter_edges_of(source_scope)
                .filter_map(|(pair, edge_id, edge)| {
                    (pair.is_paired() && !edge.data.is_dummy).then_some(edge_id)
                })
                .collect::<Vec<_>>();
            let mut generation_cache = crate::cff::generation::ExactCffGenerationCache::default();

            // Use the complete enclosing source scope for incidence:
            // exact-source reconstruction contracts the omitted prefix into
            // the quotient topology. Only the reduced owner set above selects
            // this component's active denominators.
            for (component_denominators, states) in &mut terms {
                let denominators = &component_denominators[component];
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
                    let (cff, _) = graph.cff_from_4d_denominators_in_uv_sub_lmb(
                        denominators,
                        uv_edges.iter().copied(),
                        [],
                        coordinate_lmb,
                        ExactUvSubLmbFrame::TaylorVacuum,
                        &child_cutset,
                        &options,
                        &numerator,
                        Some(&mut generation_cache),
                    )?;
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
                                .map_exact_source_numerator(&orientation.orientation)
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
                            next_states.push((
                                &carrier * &orientation.expression * &production_prefactor,
                                mapped_numerator,
                            ));
                        }
                    }
                }
                *states = next_states;
            }

            debug_tags!(#generation, #uv, #local, #four_d, #cff, #summary;
                component,
                cached_exact_cff_expressions = generation_cache.len(),
                "Cached exact CFF expressions by topology and capacity in one local-4D component wave"
            );
        }

        let active = terms
            .into_iter()
            .flat_map(|(_, states)| states)
            .fold(Atom::Zero, |sum, (carrier, numerator)| {
                sum + carrier * numerator
            });
        Ok(active)
    }
}

impl Projected4dApproximation<'_> {
    /// Project an already Taylor-expanded local 4D counterterm in the UV
    /// child's own energy frame. The contracted cograph is deliberately absent
    /// here: final assembly attaches it after the independent 4D Taylor
    /// operation has been converted to CFF.
    pub(crate) fn project_local_4d(&mut self, local: &Local4dCts) -> Result<Projected4dCts> {
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
        for sector in local.active_sectors() {
            let frozen_localizer = sector
                .frozen_lmbs()
                .iter()
                .fold(Atom::one(), |product, lmb| {
                    product * GS.localizing_integrand(lmb)
                });
            let active = self
                .localizer
                .project_factorized_taylor_sector(self.graph, sector)?;

            active_sectors.push(Projected4dSector {
                coefficient: active,
                frozen_factor: frozen_localizer,
            });
        }

        Ok(Projected4dCts::new(active_sectors))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        cff::generation::ExactCffGenerationCache,
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
        let powered = localizer.project_factorized_taylor_sector(&mut graph, &powered_sector)?;
        let cancelled =
            localizer.project_factorized_taylor_sector(&mut graph, &cancelled_sector)?;
        let dotted = localizer.project_factorized_taylor_sector(&mut graph, &dotted_sector)?;
        let one_pole = localizer.project_factorized_taylor_sector(&mut graph, &one_pole_sector)?;

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
        let mut structural_cache = ExactCffGenerationCache::default();
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
            )?;
        }

        let batched = localizer.project_factorized_taylor_sector(&mut graph, &sector)?;
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
            )?;
            let production_prefactor = Atom::num(cff.production_prefactor_factor());
            for (index, cff_term) in cff.terms {
                assert_eq!(index, CutCFFIndex::new_all_none());
                for orientation in &cff_term.orientations {
                    sequential += &orientation.expression
                        * &production_prefactor
                        * cff_term.map_exact_source_numerator(&orientation.orientation)?;
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
        let permuted = localizer.project_factorized_taylor_sector(&mut graph, &permuted)?;
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
        let sector = FourDSector::new(
            &numerator
                * denominator(EdgeIndex(0))
                * denominator(EdgeIndex(1))
                * denominator(EdgeIndex(2)),
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
        let batched = localizer.project_factorized_taylor_sector(&mut graph, &sector)?;
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
        )?;
        let inner_prefactor = Atom::num(inner_cff.production_prefactor_factor());
        for (index, cff_term) in inner_cff.terms {
            assert_eq!(index, CutCFFIndex::new_all_none());
            for orientation in &cff_term.orientations {
                states.push((
                    &orientation.expression * &inner_prefactor,
                    cff_term.map_exact_source_numerator(&orientation.orientation)?,
                ));
            }
        }
        let uv_edges = [EdgeIndex(0), EdgeIndex(1), EdgeIndex(2)];
        let mut shell_cache = ExactCffGenerationCache::default();

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
                )?;
                let prefactor = Atom::num(cff.production_prefactor_factor());
                for (index, cff_term) in cff.terms {
                    assert_eq!(index, CutCFFIndex::new_all_none());
                    for orientation in &cff_term.orientations {
                        *sum += &carrier
                            * &orientation.expression
                            * &prefactor
                            * cff_term.map_exact_source_numerator(&orientation.orientation)?;
                    }
                }
            }
        }
        assert!((cached.collect_factors() - sequential.collect_factors()).is_zero());
        assert!(
            (batched.collect_factors() - sequential.collect_factors()).is_zero(),
            "the production two-pass waves must equal fully uncached sequential component projection"
        );
        Ok(())
    }
}
