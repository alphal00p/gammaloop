use std::collections::BTreeMap;
use std::sync::Arc;
use std::time::Instant;

use color_eyre::Result;
use eyre::{WrapErr, eyre};
use idenso::{CookMode, CookSettings, IndexTooling, tensor::SymbolicTensor};
use linnet::half_edge::subgraph::{Inclusion, SubSetLike, SubSetOps};
use spenso::structure::IndexedStructure;
use spenso::structure::partial::PartialStructure;
use spenso::structure::slot::IsAbstractSlot;
use symbolica::atom::{Atom, AtomCore, AtomView};
use symbolica::symbol;
use three_dimensional_reps::RepresentationMode;
use vakint::Vakint;

use crate::integrands::process::param_builder::NumeratorDefinition;

use crate::{
    graph::{Graph, cuts::CutSet, feynman_graph::FeynmanGraph},
    model::Model,
    numerator::aind::Aind,
    processes::GraphGenerationStats,
    settings::global::GenerationSettings,
    uv::{
        Integrands, RenormalizationPart, UVOrchestrator, UVgenerationSettings, UltravioletGraph,
        approx::{
            CutStructure, OrientationProjection, local_3d::Localizer,
            projected_4d::Local4dProjectionContext,
        },
        forest::ParametricIntegrands,
        hedge_poset::Wood as HedgePosetWood,
        marker::UvMarker,
        settings::FinalIntegrandDimension,
        wood::CutWoods,
    },
};

impl Graph {
    pub(crate) fn parametric_integrands(
        &mut self,
        model: &Model,
        cut_structure: CutStructure,
        vakint: &Vakint,
        orientations: &[OrientationProjection<'_>],
        settings: &GenerationSettings,
        stats: &mut GraphGenerationStats,
    ) -> Result<BTreeMap<RepresentationMode, Vec<ParametricIntegrands>>> {
        let graph = self;
        if !matches!(settings.uv.final_integrand, FinalIntegrandDimension::ThreeD) {
            return Err(eyre!(
                "4D parametric UV integrands are not supported yet; this mode is planned for a future implementation"
            ));
        }
        validate_local_counterterm_schemes(graph, &cut_structure, &settings.uv)?;

        settings.validate_medium_uv_compatibility()?;

        let result = match settings.uv.orchestrator {
            UVOrchestrator::LegacyDagForest => legacy_parametric_integrands(
                graph,
                model,
                cut_structure,
                vakint,
                orientations,
                settings,
                stats,
            ),
            UVOrchestrator::HedgePoset => hedge_poset_parametric_integrands(
                graph,
                model,
                cut_structure,
                vakint,
                orientations,
                settings,
                stats,
            ),
            UVOrchestrator::Compare => compare_parametric_integrands(
                graph,
                model,
                cut_structure,
                vakint,
                orientations,
                settings,
                stats,
            ),
        }?;
        let marker = UvMarker::new(&settings.uv);
        result
            .into_iter()
            .map(|(representation, integrands)| {
                let started = Instant::now();
                let integrands = integrands
                    .into_iter()
                    .map(|integrands| integrands.map_expressions(|atom| Ok(marker.finish(atom))))
                    .collect::<Result<_>>()?;
                stats.representation_mut(representation).total_time += started.elapsed();
                Ok((representation, integrands))
            })
            .collect()
    }
}

impl UVOrchestrator {
    pub(crate) fn renormalization_part(
        self,
        graph: &mut Graph,
        model: &Model,
        orientation: OrientationProjection<'_>,
        settings: &UVgenerationSettings,
    ) -> Result<RenormalizationPart> {
        let settings = UVgenerationSettings {
            final_integrand: FinalIntegrandDimension::FourD,
            ..settings.clone()
        };
        validate_local_counterterm_schemes(graph, &CutStructure::empty(graph), &settings)?;
        let mut result = match self {
            Self::LegacyDagForest => {
                legacy_renormalization_part(graph, model, orientation, &settings)
            }
            Self::HedgePoset => hedge_poset_renormalization_part(graph, model, &settings),
            Self::Compare => compare_renormalization_part(graph, model, orientation, &settings),
        }?;
        result.expression = UvMarker::new(&settings).finish(&result.expression);
        Ok(result)
    }
}

pub(crate) fn validate_local_counterterm_schemes(
    graph: &Graph,
    cuts: &CutStructure,
    settings: &UVgenerationSettings,
) -> Result<()> {
    if !settings.subtract_uv {
        return Ok(());
    }

    let mut uncut = graph.full_filter();
    uncut.subtract_with(&graph.initial_state_cut.left);
    for (cut_index, cut) in cuts.cuts.iter().enumerate() {
        let subgraph = uncut.subtract(&cut.union);
        let spinneys = graph
            .classified_spinneys(
                &subgraph,
                &graph.model,
                settings,
                &graph.loop_momentum_basis,
            )
            .wrap_err_with(|| {
                format!(
                    "graph '{}' cut {cut_index}: failed to classify UV counterterm components",
                    graph.name
                )
            })?;
        let ir_components = spinneys.iter().filter(|spinney| {
            spinney.n_components() == 1
                && spinney.renormalization_scheme == crate::uv::ApproximationType::IR
        });
        let os_components = spinneys.iter().filter(|spinney| {
            spinney.n_components() == 1
                && spinney.renormalization_scheme == crate::uv::ApproximationType::OS
        });
        let pole_components = spinneys.iter().filter(|spinney| {
            spinney.n_components() == 1
                && spinney.renormalization_scheme == crate::uv::ApproximationType::PolePart
        });
        if let Some((child, parent)) = ir_components.clone().find_map(|child| {
            (child.dod > 0).then(|| {
                spinneys
                    .iter()
                    .find(|parent| {
                        parent.n_components() == 1
                            && parent.dod > 0
                            && parent.renormalization_scheme == crate::uv::ApproximationType::MUV
                            && parent.filter().includes(child.filter())
                            && parent.filter() != child.filter()
                    })
                    .map(|parent| (child, parent))
            })?
        }) {
            return Err(eyre!(
                "graph '{}' cut {cut_index}: positive-degree MUV component {} (d={}) cannot contain local soft/IR component {} (d={}) in phase 1; its reduced-cograph subtraction is not local without the containing soft branch, so assign IR to the containing component or wait for the deferred finite scheme-change/integrated policy",
                graph.name,
                parent.filter().string_label(),
                parent.dod,
                child.filter().string_label(),
                child.dod,
            ));
        }
        let ir_pole_wood = ir_components.clone().find_map(|ir| {
            pole_components.clone().find_map(|pole| {
                let nested =
                    ir.filter().includes(pole.filter()) || pole.filter().includes(ir.filter());
                let disjoint_union = !ir.filter().intersects(pole.filter())
                    && spinneys.iter().any(|candidate| {
                        candidate.n_components() > 1
                            && candidate.filter() == &ir.filter().union(pole.filter())
                    });
                if nested {
                    Some((ir, pole, "nested"))
                } else if disjoint_union {
                    Some((ir, pole, "disconnected"))
                } else {
                    None
                }
            })
        });
        if let Some((ir, pole, relationship)) = ir_pole_wood {
            return Err(eyre!(
                "graph '{}' cut {cut_index}: local soft/IR and PolePart counterterms cannot be combined in one UV wood before the scheme-specific integrated policy is implemented; component {} and component {} form a {relationship} wood",
                graph.name,
                ir.filter().string_label(),
                pole.filter().string_label(),
            ));
        }
        if os_components.count() > 0 && settings.generate_integrated {
            return Err(eyre!(
                "local on-shell/OS counterterms are local-only in phase 1; set uv.generate_integrated=false"
            ));
        }
    }

    Ok(())
}

fn legacy_parametric_integrands(
    graph: &mut Graph,
    model: &Model,
    cut_structure: CutStructure,
    vakint: &Vakint,
    orientations: &[OrientationProjection<'_>],
    settings: &GenerationSettings,
    stats: &mut GraphGenerationStats,
) -> Result<BTreeMap<RepresentationMode, Vec<ParametricIntegrands>>> {
    let cut_woods = CutWoods::new(cut_structure, graph, model, settings)?;
    let mut cut_forests = cut_woods.unfold(graph);
    let mut projection_context = Local4dProjectionContext::default();
    orientations
        .iter()
        .enumerate()
        .map(|(index, &orientation)| {
            let started = Instant::now();
            let shared = cut_forests.compute(
                graph,
                vakint,
                orientation,
                &settings.uv,
                index == 0,
                &mut projection_context,
            )?;
            let integrands = cut_forests.orientation_parametric_exprs(graph, &settings.uv)?;
            let representation = orientation.cff_options()?.representation;
            stats.representation_mut(representation).total_time +=
                started.elapsed().saturating_sub(shared);
            Ok((representation, integrands))
        })
        .collect()
}

fn hedge_poset_parametric_integrands(
    graph: &mut Graph,
    model: &Model,
    cut_structure: CutStructure,
    vakint: &Vakint,
    orientations: &[OrientationProjection<'_>],
    settings: &GenerationSettings,
    stats: &mut GraphGenerationStats,
) -> Result<BTreeMap<RepresentationMode, Vec<ParametricIntegrands>>> {
    let wood = HedgePosetWood::new(cut_structure, graph, model, settings)?;
    let mut forests = wood.unfold();
    let mut projection_context = Local4dProjectionContext::default();
    orientations
        .iter()
        .enumerate()
        .map(|(index, &orientation)| {
            let started = Instant::now();
            let shared_before = projection_context.shared_preparation_time;
            if index == 0 {
                forests.prepare_four_d(
                    graph,
                    model,
                    vakint,
                    &settings.uv,
                    &mut projection_context,
                )?;
            }
            forests.compute(graph, orientation, &settings.uv, &mut projection_context)?;
            let shared = projection_context.shared_preparation_time - shared_before;
            let integrands = forests.orientation_parametric_exprs(graph, &settings.uv)?;
            let representation = orientation.cff_options()?.representation;
            stats.representation_mut(representation).total_time +=
                started.elapsed().saturating_sub(shared);
            Ok((representation, integrands))
        })
        .collect()
}

fn compare_parametric_integrands(
    graph: &mut Graph,
    model: &Model,
    cut_structure: CutStructure,
    vakint: &Vakint,
    orientations: &[OrientationProjection<'_>],
    settings: &GenerationSettings,
    stats: &mut GraphGenerationStats,
) -> Result<BTreeMap<RepresentationMode, Vec<ParametricIntegrands>>> {
    let mut hedge_graph = graph.clone();
    let legacy = legacy_parametric_integrands(
        graph,
        model,
        cut_structure.clone(),
        vakint,
        orientations,
        settings,
        stats,
    )?;
    let hedge = hedge_poset_parametric_integrands(
        &mut hedge_graph,
        model,
        cut_structure,
        vakint,
        orientations,
        settings,
        stats,
    )?;

    for (representation, legacy) in &legacy {
        let started = Instant::now();
        ParametricIntegrandsComparison {
            legacy,
            hedge: &hedge[representation],
        }
        .compare()?;
        stats.representation_mut(*representation).total_time += started.elapsed();
    }
    Ok(legacy)
}

fn legacy_renormalization_part(
    graph: &mut Graph,
    model: &Model,
    orientation: OrientationProjection<'_>,
    settings: &UVgenerationSettings,
) -> Result<RenormalizationPart> {
    let mut vk_settings = settings.vakint.true_settings();
    vk_settings.project_onto_tensor_integrals =
        settings.project_integrated_uv_cts_onto_tensor_integrals;
    let wood = graph
        .wood_with_settings(
            &graph.no_dummy(),
            model,
            settings,
            &graph.loop_momentum_basis,
        )
        .wrap_err_with(|| {
            format!(
                "graph '{}': failed to classify uncut UV counterterm components",
                graph.name
            )
        })?;
    // MUV renormalization extracts the finite term, so retain one term beyond
    // the maximal pole order, as in the other forest integration paths.
    vk_settings.number_of_terms_in_epsilon_expansion = wood.max_loops as i64 + 1;

    let mut forest = wood.unfold(graph, &graph.loop_momentum_basis);
    let vk = (crate::utils::vakint()?, &vk_settings);
    let cuts = CutSet::empty(graph.n_hedges());
    forest.compute(
        graph,
        vk,
        Localizer::new(&cuts, orientation),
        settings,
        &mut super::approx::projected_4d::Local4dProjectionContext::default(),
        true,
    )?;

    forest.renormalization_part_of_ends(graph, settings)
}

fn hedge_poset_renormalization_part(
    graph: &mut Graph,
    model: &Model,
    settings: &UVgenerationSettings,
) -> Result<RenormalizationPart> {
    let cuts = CutStructure::empty(graph);
    let generation_settings = GenerationSettings {
        uv: settings.clone(),
        ..GenerationSettings::default()
    };
    let wood = HedgePosetWood::new(cuts, graph, model, &generation_settings)?;
    let mut forest = wood.unfold();
    forest.integrate(graph, model, crate::utils::vakint()?, settings)?;
    forest.renormalization_part_of_ends(graph, settings)
}

fn compare_renormalization_part(
    graph: &mut Graph,
    model: &Model,
    orientation: OrientationProjection<'_>,
    settings: &UVgenerationSettings,
) -> Result<RenormalizationPart> {
    let mut hedge_graph = graph.clone();
    let legacy = legacy_renormalization_part(graph, model, orientation, settings)?;
    let hedge = hedge_poset_renormalization_part(&mut hedge_graph, model, settings)?;

    RenormalizationComparison {
        legacy: &legacy,
        hedge: &hedge,
    }
    .compare()?;
    Ok(legacy)
}

struct ParametricIntegrandsComparison<'a> {
    legacy: &'a [ParametricIntegrands],
    hedge: &'a [ParametricIntegrands],
}

impl ParametricIntegrandsComparison<'_> {
    fn compare(&self) -> Result<()> {
        if self.legacy.len() != self.hedge.len() {
            return Err(eyre!(
                "UV orchestrator compare mismatch: legacy produced {} cut integrands, hedge-poset produced {}",
                self.legacy.len(),
                self.hedge.len()
            ));
        }

        for (cut_index, (legacy, hedge)) in self.legacy.iter().zip(self.hedge.iter()).enumerate() {
            if legacy.cuts != hedge.cuts {
                return Err(eyre!(
                    "UV orchestrator compare mismatch at cut {}: cut structures differ",
                    cut_index
                ));
            }
            IntegrandMapComparison {
                cut_index,
                legacy: &legacy.integrands,
                hedge: &hedge.integrands,
            }
            .compare()?;
        }

        Ok(())
    }
}

struct IntegrandMapComparison<'a> {
    cut_index: usize,
    legacy: &'a Integrands,
    hedge: &'a Integrands,
}

impl IntegrandMapComparison<'_> {
    fn compare(&self) -> Result<()> {
        let legacy = self.legacy.resolved_scalars()?;
        let hedge = self.hedge.resolved_scalars()?;
        let mut legacy_cache = BTreeMap::new();
        let mut hedge_cache = BTreeMap::new();
        legacy
            .checked_zip(&hedge, |key, legacy_expr, hedge_expr| {
                let legacy_tensor = ComparableExpr::materialize(
                    legacy_expr.as_view(),
                    legacy.numerators(),
                    &mut legacy_cache,
                )?;
                let hedge_tensor = ComparableExpr::materialize(
                    hedge_expr.as_view(),
                    hedge.numerators(),
                    &mut hedge_cache,
                )?;
                if !ComparableExpr::new(&legacy_tensor)
                    .equivalent_to(&ComparableExpr::new(&hedge_tensor))?
                {
                    crate::debug_tags!(#uv, #compare, #mismatch;
                        cut_index = self.cut_index,
                        residue = ?key,
                        file.legacy = legacy_expr.to_canonical_string(),
                        file.hedge = hedge_expr.to_canonical_string(),
                        "UV orchestrator expressions differ at the shared residue boundary"
                    );
                    return Err(eyre!(
                        "UV orchestrator compare mismatch at cut {} residue {:?}",
                        self.cut_index,
                        key
                    ));
                }
                Ok(Atom::Zero)
            })
            .wrap_err_with(|| {
                format!(
                    "while comparing UV orchestrator residue structure at cut {}",
                    self.cut_index
                )
            })?;

        Ok(())
    }
}

struct RenormalizationComparison<'a> {
    legacy: &'a RenormalizationPart,
    hedge: &'a RenormalizationPart,
}

impl RenormalizationComparison<'_> {
    fn compare(&self) -> Result<()> {
        if !ComparableExpr::new(&self.legacy.cooked_tensor()?)
            .equivalent_to(&ComparableExpr::new(&self.hedge.cooked_tensor()?))?
        {
            return Err(eyre!(
                "UV orchestrator compare mismatch in integrated renormalization part"
            ));
        }
        Ok(())
    }
}

struct ComparableExpr<'a> {
    tensor: &'a SymbolicTensor<IndexedStructure>,
}

impl<'a> ComparableExpr<'a> {
    fn new(tensor: &'a SymbolicTensor<IndexedStructure>) -> Self {
        Self { tensor }
    }

    /// Compose retained families through the raw residue shell. Only syntax
    /// without a retained owner is admitted; expanded numerator bodies keep
    /// their established interfaces and never pass through root inference.
    fn materialize(
        atom: AtomView<'_>,
        numerators: &[Arc<NumeratorDefinition>],
        cache: &mut BTreeMap<Atom, SymbolicTensor<IndexedStructure>>,
    ) -> Result<SymbolicTensor<IndexedStructure>> {
        let original = atom.to_owned();
        if let Some(value) = cache.get(&original) {
            return Ok(value.clone());
        }
        let family = symbol!("gammalooprs::uv::numerator_family");
        let cooking = CookSettings::indices()
            .with_mode(CookMode::ReversibleEncoding)
            .with_representation_payloads(true, true);
        let value = if !atom.contains_symbol(family) {
            // This is the first typed boundary for the raw residue shell or
            // an actual argument supplied to a retained function.
            SymbolicTensor::infer_indexed(cooking.try_cook(atom)?)?
        } else {
            match atom {
                AtomView::Add(sum) => {
                    let mut terms = sum.iter();
                    let mut result = Self::materialize(
                        terms.next().expect("nonempty symbolic sum"),
                        numerators,
                        cache,
                    )?;
                    for term in terms {
                        result = result.try_add(&Self::materialize(term, numerators, cache)?)?;
                    }
                    result
                }
                AtomView::Mul(product) => {
                    let mut factors = product.iter();
                    let mut result = Self::materialize(
                        factors.next().expect("nonempty symbolic product"),
                        numerators,
                        cache,
                    )?;
                    for factor in factors {
                        result = result.multiply(Self::materialize(factor, numerators, cache)?)?;
                    }
                    result
                }
                AtomView::Pow(power) => {
                    let (base, exponent) = power.get_base_exp();
                    Self::materialize(base, numerators, cache)?
                        .pow(&Self::materialize(exponent, numerators, cache)?)?
                }
                AtomView::Fun(call) if call.get_symbol() == family => {
                    let mut entries = numerators.iter().filter(|entry| {
                        entry.tags.len() + entry.args.len() == call.get_nargs()
                            && call
                                .iter()
                                .zip(&entry.tags)
                                .all(|(arg, tag)| arg == tag.as_view())
                    });
                    let entry = entries
                        .next()
                        .ok_or_else(|| eyre!("unresolved retained numerator call {atom}"))?;
                    if entries.next().is_some() {
                        return Err(eyre!("ambiguous retained numerator call {atom}"));
                    }
                    let mut scalar_binding = true;
                    for argument in call.iter().skip(entry.tags.len()) {
                        let argument = Self::materialize(argument, numerators, cache)?;
                        let argument: SymbolicTensor<PartialStructure> = argument.into();
                        scalar_binding &=
                            argument.is_scalar() && argument.explicit_indices()?.next().is_none();
                    }
                    let preserves_ports =
                        entry.rhs.structure().logical_slots().iter().all(|slot| {
                            let slot = cooking.uncook(slot.to_atom().as_view());
                            atom.replace_multiple([entry.with_rhs(slot.clone()).replacement()])
                                == slot
                        });
                    let expression =
                        cooking.try_cook(atom.replace_multiple([entry.replacement()]).as_view())?;
                    if expression.contains_symbol(family) {
                        return Err(eyre!("unresolved or cyclic retained numerator family call"));
                    }
                    if scalar_binding && preserves_ports {
                        entry.rhs.with_materialized_expression(
                            "retained numerator binding",
                            expression,
                        )?
                    } else {
                        // Arbitrary actual arguments or index-valued formals
                        // are a user substitution boundary, not a built-in
                        // scalar binding with an established interface.
                        SymbolicTensor::infer_indexed(expression)?
                    }
                }
                _ => {
                    return Err(eyre!(
                        "retained numerator call occurs outside tensor arithmetic: {atom}"
                    ));
                }
            }
        };
        cache.insert(original, value.clone());
        Ok(value)
    }

    fn equivalent_to(&self, other: &Self) -> Result<bool> {
        let left = self.normalized()?;
        let right = other.normalized()?;

        // Backend-local topology labels may differ on contracted indices.
        if left.collect_factors() == right.collect_factors() {
            return Ok(true);
        }

        Ok(left.canonize(Aind::Dummy)?.collect_factors()
            == right.canonize(Aind::Dummy)?.collect_factors())
    }

    fn normalized(&self) -> Result<Atom> {
        let tensor = self.tensor.map_representation_dimensions(|_, dimension| {
            dimension.replace(crate::utils::GS.dim).with(4)
        })?;
        let cooking = CookSettings::indices()
            .with_mode(CookMode::ReversibleEncoding)
            .with_representation_payloads(true, true);
        let normalized = tensor
            .with_produced_expression(tensor.expression().replace(crate::utils::GS.dim).with(4))?
            .simplify_algebra(&idenso::tensor::AlgebraSettings {
                color: Some(Default::default()),
                ..Default::default()
            })?
            .contract(idenso::tensor::ContractSettings {
                collect_chains: false,
                collect_traces: false,
                ..Default::default()
            })?
            .to_dots()?;
        Ok(cooking.uncook(normalized.expression().as_view()))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        finalized_runtime_dot,
        initialisation::test_initialise,
        utils::load_generic_model,
        uv::{ApproximationType, CTRenormalizationRule, RenormalizationPrescriptionSettings},
    };
    use idenso::{bis, gamma};
    use linnet::half_edge::{involution::EdgeIndex, subgraph::SubSetLike};
    use spenso::{chain, mink};
    use std::collections::BTreeSet;
    use symbolica::evaluate::InliningPolicy;
    use symbolica::symbol;

    use crate::integrands::process::param_builder::FnMapEntry;

    #[test]
    fn compare_canonicalizes_contracted_uv_indices() {
        crate::initialisation::test_initialise().unwrap();
        let spectator = symbolica::parse!("2*(compare_a+compare_b)*(compare_c+compare_d)");
        let expression = |topology| {
            let contracted = mink!(4, Atom::from(Aind::UVTerm(topology, 2)));
            let fixed = mink!(4, Atom::from(Aind::Edge(2, 1)));
            let start = bis!(4, Atom::from(Aind::Hedge(0, 0)));
            let end = bis!(4, Atom::from(Aind::Hedge(1, 0)));
            let common = &spectator * Atom::var(symbol!("compare_common_factor"));
            let term = |index: Atom| {
                common.clone()
                    * chain!(start.clone(), end.clone(), gamma!(index.clone()))
                    * crate::utils::GS.emr_vec(EdgeIndex::from(0), index)
            };
            term(contracted) + term(fixed)
        };
        let legacy = expression(1);
        let hedge = expression(0);
        assert_ne!(legacy, hedge, "the source dummy-index charts must differ");
        let cooking = CookSettings::indices()
            .with_mode(CookMode::ReversibleEncoding)
            .with_representation_payloads(true, true);
        let legacy =
            SymbolicTensor::infer_indexed(cooking.try_cook(legacy.as_view()).unwrap()).unwrap();
        let hedge =
            SymbolicTensor::infer_indexed(cooking.try_cook(hedge.as_view()).unwrap()).unwrap();
        let legacy = ComparableExpr::new(&legacy);
        let hedge = ComparableExpr::new(&hedge);

        // Shared tensor reduction now contracts these indices during normalization.
        assert_eq!(legacy.normalized().unwrap(), hedge.normalized().unwrap());
        for expression in [&legacy, &hedge] {
            let normalized = expression.normalized().unwrap();
            for factor in [
                symbolica::parse!("compare_a+compare_b"),
                symbolica::parse!("compare_c+compare_d"),
            ] {
                assert!(normalized.contains(&factor), "lost scalar factor {factor}");
            }
        }
        assert!(
            legacy
                .equivalent_to(&hedge)
                .expect("test expressions should canonicalize")
        );
    }

    #[test]
    fn compare_retains_family_and_renormalization_open_ports() -> Result<()> {
        test_initialise()?;
        let cooking = CookSettings::indices()
            .with_mode(CookMode::ReversibleEncoding)
            .with_representation_payloads(true, true);
        let start = bis!(4, Atom::from(Aind::Hedge(0, 0)));
        let end = bis!(4, Atom::from(Aind::Hedge(1, 0)));
        let index = mink!(4, Atom::from(Aind::UVTerm(0, 2)));
        let raw = chain!(start, end, gamma!(index.clone()))
            * crate::utils::GS.emr_vec(EdgeIndex::from(0), index);
        let parameter = symbol!("compare_retained_parameter");
        let source = SymbolicTensor::infer_indexed(
            cooking.try_cook((raw * Atom::var(parameter)).as_view())?,
        )?;
        let family = symbol!("gammalooprs::uv::numerator_family");
        let call = family.call_args([Atom::num(17), Atom::num(3)]);
        for source in [source.clone(), source.with_produced_expression(Atom::Zero)?] {
            let definition = Arc::new(FnMapEntry {
                lhs: family.call_args([Atom::num(17), Atom::var(parameter)]),
                rhs: source,
                tags: vec![Atom::num(17)],
                args: vec![parameter.into()],
                inlining: InliningPolicy::Always,
                is_alias: false,
            });
            let root = &call * symbolica::parse!("(compare_e+compare_f)*(compare_g+compare_h)");
            let value = ComparableExpr::materialize(
                root.as_view(),
                &[Arc::clone(&definition)],
                &mut BTreeMap::new(),
            )?;
            assert_eq!(value.structure(), definition.rhs.structure());
            assert_eq!(value.rank(), 2);
            let part = RenormalizationPart::legacy(
                cooking.uncook(value.expression().as_view()),
                value.structure().clone(),
            );
            let recovered = part.cooked_tensor()?;
            assert_eq!(recovered.structure(), value.structure());
            assert!(ComparableExpr::new(&recovered).equivalent_to(&ComparableExpr::new(&value))?);
        }
        Ok(())
    }

    #[test]
    fn integrated_soft_scheme_is_allowed_before_forest_generation() {
        test_initialise().unwrap();
        let model = load_generic_model("scalars");
        let graph: Graph = Graph::import_dot(
            include_str!("../../../../tests/resources/graphs/scalar/dod2_bubble.dot"),
            &model,
        )
        .map(|mut graphs| {
            assert_eq!(graphs.len(), 1, "one fixture graph expected");
            graphs.remove(0)
        })
        .unwrap();
        let settings = UVgenerationSettings {
            generate_integrated: true,
            renormalization_prescription: RenormalizationPrescriptionSettings {
                log_divergent: ApproximationType::IR,
                massive_power_divergent: ApproximationType::IR,
                massless_power_divergent: ApproximationType::IR,
                ..Default::default()
            },
            ..Default::default()
        };

        validate_local_counterterm_schemes(&graph, &CutStructure::empty(&graph), &settings)
            .expect("the shared integration path supports soft counterterms");
    }

    #[test]
    fn incompatible_selected_spinney_route_is_reported_with_context() {
        test_initialise().unwrap();
        let model = load_generic_model("scalars");
        let mut graph: Graph = Graph::import_dot(
            include_str!("../../../../tests/resources/graphs/scalar/dod2_bubble.dot"),
            &model,
        )
        .map(|mut graphs| {
            assert_eq!(graphs.len(), 1, "one fixture graph expected");
            graphs.remove(0)
        })
        .unwrap();
        let settings = UVgenerationSettings {
            generate_integrated: false,
            renormalization_prescription: RenormalizationPrescriptionSettings {
                log_divergent: ApproximationType::IR,
                massive_power_divergent: ApproximationType::IR,
                massless_power_divergent: ApproximationType::IR,
                ..Default::default()
            },
            ..Default::default()
        };
        let cuts = CutStructure::empty(&graph);
        let mut uncut = graph.full_filter();
        uncut.subtract_with(&graph.initial_state_cut.left);
        let component = graph
            .spinneys(&uncut)
            .into_iter()
            .find(|spinney| {
                !spinney.is_empty()
                    && graph.compute_dod(&spinney.filter) >= 0
                    && graph.underlying.connected_components(spinney).len() == 1
            })
            .expect("the DOD-2 bubble must contain a divergent connected component");
        let component_label = component.string_label();
        let dod = graph.compute_dod(&component.filter);
        assert_eq!(
            graph.approximation_scheme(&component.filter, &graph.model, &settings, dod),
            ApproximationType::IR
        );

        graph.loop_momentum_basis.loop_edges.clear();
        let error = validate_local_counterterm_schemes(&graph, &cuts, &settings)
            .expect_err("an incompatible selected component route must fail closed");
        let message = format!("{error:#}");
        assert!(
            message.contains(&format!("graph '{}'", graph.name)),
            "{message}"
        );
        assert!(message.contains("cut 0"), "{message}");
        assert!(message.contains(&component_label), "{message}");
        assert!(message.contains("scheme IR"), "{message}");
        assert!(message.contains(&format!("DOD {dod}")), "{message}");
        assert!(
            message.contains("no loop-momentum basis compatible with the parent basis"),
            "{message}"
        );
    }

    #[test]
    fn integrated_on_shell_scheme_is_rejected_before_forest_generation() {
        test_initialise().unwrap();
        let model = load_generic_model("scalars");
        let graph: Graph = Graph::import_dot(
            include_str!("../../../../tests/resources/graphs/scalar/dod2_bubble.dot"),
            &model,
        )
        .map(|mut graphs| {
            assert_eq!(graphs.len(), 1, "one fixture graph expected");
            graphs.remove(0)
        })
        .unwrap();
        let settings = UVgenerationSettings {
            generate_integrated: true,
            renormalization_prescription: RenormalizationPrescriptionSettings {
                log_divergent: ApproximationType::OS,
                massive_power_divergent: ApproximationType::OS,
                massless_power_divergent: ApproximationType::OS,
                ..Default::default()
            },
            ..Default::default()
        };

        let error =
            validate_local_counterterm_schemes(&graph, &CutStructure::empty(&graph), &settings)
                .expect_err("integrated local on-shell generation must fail at the orchestrator");
        let message = error.to_string();
        assert!(message.contains("local on-shell/OS counterterms"));
        assert!(message.contains("generate_integrated=false"));
    }

    #[test]
    fn soft_child_rejects_positive_degree_muv_parent() {
        test_initialise().unwrap();
        let graph: Graph = finalized_runtime_dot!(digraph nested_positive_degree_hu { graph [projector=1]
                edge [particle=scalar_1 num=1]
                node [num=1]
                ext [style=invis]
                ext -> v1:0 [id=0 sink="{ufo_order:0}"]
                v4:1 -> ext [id=1 source="{ufo_order:0}"]
                v1 -> v2 [id=2 lmb_id=1 num="Q(2,spenso::cind(0))^2" source="{ufo_order:1}" sink="{ufo_order:0}"]
                v2 -> v3 [id=3 particle=scalar_2 lmb_id=0 num="Q(3,spenso::cind(0))" source="{ufo_order:1}" sink="{ufo_order:0}"]
                v2 -> v3 [id=4 particle=scalar_2 source="{ufo_order:2}" sink="{ufo_order:1}"]
                v3 -> v4 [id=5 source="{ufo_order:2}" sink="{ufo_order:1}"]
                v1 -> v4 [id=6 source="{ufo_order:2}" sink="{ufo_order:2}"]
            }, "scalars")
        .unwrap();
        let baseline = UVgenerationSettings {
            generate_integrated: false,
            ..Default::default()
        };
        let spinneys = graph
            .classified_spinneys(
                &graph.full_filter(),
                &graph.model,
                &baseline,
                &graph.loop_momentum_basis,
            )
            .unwrap();
        let (child, parent) = spinneys
            .iter()
            .find_map(|child| {
                if child.n_components() == 1 && child.dod > 0 {
                    spinneys
                        .iter()
                        .find(|parent| {
                            parent.n_components() == 1
                                && parent.dod > 0
                                && parent.filter() != child.filter()
                                && parent.filter().includes(child.filter())
                        })
                        .map(|parent| (child, parent))
                } else {
                    None
                }
            })
            .expect("the fixture must contain a positive-degree child and parent");
        assert_ne!(
            graph.ct_identifier(child.filter(), &graph.model),
            graph.ct_identifier(parent.filter(), &graph.model),
            "the fixture must allow the child and parent schemes to be selected independently"
        );
        let settings = UVgenerationSettings {
            generate_integrated: false,
            renormalization_prescription: RenormalizationPrescriptionSettings {
                overrides: vec![CTRenormalizationRule::new(
                    graph.ct_identifier(child.filter(), &graph.model),
                    ApproximationType::IR,
                )],
                ..Default::default()
            },
            ..Default::default()
        };
        let classified = graph
            .classified_spinneys(
                &graph.full_filter(),
                &graph.model,
                &settings,
                &graph.loop_momentum_basis,
            )
            .unwrap();
        assert!(classified.iter().any(|spinney| {
            spinney.filter() == child.filter()
                && spinney.renormalization_scheme == ApproximationType::IR
        }));
        assert!(classified.iter().any(|spinney| {
            spinney.filter() == parent.filter()
                && spinney.renormalization_scheme == ApproximationType::MUV
        }));
        // Appendix B.15--B.16 needs the parent's S-US branch to cancel the
        // reduced-cograph limit of the finite child refinement. A positive-
        // degree ordinary parent has no such branch in the local-only phase.
        let error =
            validate_local_counterterm_schemes(&graph, &CutStructure::empty(&graph), &settings)
                .expect_err("a positive-degree U parent over an H child must be rejected");
        let message = error.to_string();
        assert!(
            message.contains("positive-degree MUV component"),
            "{message}"
        );
        assert!(
            message.contains("deferred finite scheme-change/integrated policy"),
            "{message}"
        );
    }

    #[test]
    fn soft_and_pole_part_components_cannot_share_a_wood() {
        test_initialise().unwrap();
        let model = load_generic_model("sm");
        let graph: Graph = Graph::import_dot(
            include_str!("../../../../tests/resources/graphs/dgse.dot"),
            &model,
        )
        .map(|mut graphs| {
            assert_eq!(graphs.len(), 1, "one fixture graph expected");
            graphs.remove(0)
        })
        .unwrap();
        let identifiers = graph
            .spinneys(&graph.full_filter())
            .into_iter()
            .filter(|spinney| {
                graph.compute_dod(&spinney.filter) >= 0
                    && graph.underlying.connected_components(spinney).len() == 1
            })
            .map(|spinney| graph.ct_identifier(&spinney.filter, &graph.model))
            .collect::<BTreeSet<_>>();
        assert!(
            identifiers.len() >= 2,
            "DGSE must expose two distinct connected counterterm identifiers"
        );
        let mut identifiers = identifiers.into_iter();
        let soft = identifiers.next().unwrap();
        let pole = identifiers.next().unwrap();
        let settings = UVgenerationSettings {
            generate_integrated: false,
            renormalization_prescription: RenormalizationPrescriptionSettings {
                log_divergent: ApproximationType::Unsubtracted,
                massive_power_divergent: ApproximationType::Unsubtracted,
                massless_power_divergent: ApproximationType::Unsubtracted,
                overrides: vec![
                    CTRenormalizationRule::new(soft, ApproximationType::IR),
                    CTRenormalizationRule::new(pole, ApproximationType::PolePart),
                ],
            },
            ..Default::default()
        };

        let error =
            validate_local_counterterm_schemes(&graph, &CutStructure::empty(&graph), &settings)
                .expect_err("soft and pole-part components must not share a phase-1 wood");
        assert!(error.to_string().contains("cannot be combined"));
    }

    #[test]
    fn disjoint_soft_and_pole_part_components_cannot_share_a_wood() {
        test_initialise().unwrap();
        let graph: Graph = Graph::import_dot(
            stringify!(digraph asymmetric_scalar_spectacles {
                num = "-1𝑖/2"
                node [num = "1"]
                edge [particle = scalar_1]

                ext [style = invis]
                ext -> a:0 [id = 5]
                d:1 -> ext [id = 6]
                a -> b [
                    id = 0
                    num = "Q(0,spenso::cind(1))^2+Q(0,spenso::cind(2))^2+Q(0,spenso::cind(3))^2"
                ]
                a -> b [id = 1]
                b -> c [id = 2 particle = scalar_0]
                c -> d [
                    id = 3
                    particle = scalar_2
                    num = "Q(3,spenso::cind(1))^2+Q(3,spenso::cind(2))^2+Q(3,spenso::cind(3))^2"
                ]
                c -> d [id = 4 particle = scalar_2]
            }),
            &crate::utils::load_generic_model("scalars"),
        )
        .map(|mut graphs| {
            assert_eq!(graphs.len(), 1, "one fixture graph expected");
            graphs.remove(0)
        })
        .unwrap();
        let edge_ids = |filter: &linnet::half_edge::subgraph::SuBitGraph| {
            filter
                .included_iter()
                .map(|hedge| graph.underlying[&hedge].0)
                .collect::<BTreeSet<_>>()
        };
        let union = graph
            .spinneys(&graph.full_filter())
            .into_iter()
            .find(|spinney| graph.underlying.connected_components(spinney).len() == 2)
            .expect("the asymmetric spectacles fixture must contain one disconnected spinney");
        assert_eq!(edge_ids(&union.filter), BTreeSet::from([0, 1, 3, 4]));

        let mut components = graph
            .underlying
            .connected_components(&union)
            .into_iter()
            .map(|component| (edge_ids(&component), component))
            .collect::<Vec<_>>();
        components.sort_by_key(|(edges, _)| edges.clone());
        assert_eq!(components.len(), 2);
        assert_eq!(components[0].0, BTreeSet::from([0, 1]));
        assert_eq!(components[1].0, BTreeSet::from([3, 4]));
        let left = &components[0].1;
        let right = &components[1].1;
        assert!(
            !left.intersects(right),
            "the selected components must exercise the disjoint-union validator branch"
        );
        assert_eq!(union.filter, left.union(right));

        let left_identifier = graph.ct_identifier(left, &graph.model);
        let right_identifier = graph.ct_identifier(right, &graph.model);
        assert_ne!(
            left_identifier, right_identifier,
            "different internal scalar species must make the two component rules independent"
        );
        let settings = UVgenerationSettings {
            generate_integrated: false,
            renormalization_prescription: RenormalizationPrescriptionSettings {
                log_divergent: ApproximationType::Unsubtracted,
                massive_power_divergent: ApproximationType::Unsubtracted,
                massless_power_divergent: ApproximationType::Unsubtracted,
                overrides: vec![
                    CTRenormalizationRule::new(left_identifier, ApproximationType::IR),
                    CTRenormalizationRule::new(right_identifier, ApproximationType::PolePart),
                ],
            },
            ..Default::default()
        };
        let classified = graph
            .classified_spinneys(
                &graph.full_filter(),
                &graph.model,
                &settings,
                &graph.loop_momentum_basis,
            )
            .unwrap();
        assert_eq!(
            classified.len(),
            3,
            "both connected factors and their retained disconnected union must be classified"
        );
        let selected_left = classified
            .iter()
            .find(|spinney| spinney.filter() == left)
            .expect("the left spectacles component must remain selected");
        let selected_right = classified
            .iter()
            .find(|spinney| spinney.filter() == right)
            .expect("the right spectacles component must remain selected");
        let selected_union = classified
            .iter()
            .find(|spinney| spinney.filter() == &union.filter)
            .expect("the disconnected union must remain selected when both factors are selected");
        assert_eq!(selected_left.renormalization_scheme, ApproximationType::IR);
        assert_eq!(
            selected_right.renormalization_scheme,
            ApproximationType::PolePart
        );
        assert_eq!(
            selected_union.renormalization_scheme,
            ApproximationType::MUV,
            "the disconnected union must retain only its neutral aggregate scheme"
        );

        let error =
            validate_local_counterterm_schemes(&graph, &CutStructure::empty(&graph), &settings)
                .expect_err("disjoint soft and pole-part factors must not share a phase-1 wood");
        let message = format!("{error:#}");
        assert!(
            message.contains("local soft/IR and PolePart counterterms")
                && message.contains("cannot be combined in one UV wood")
                && message.contains("scheme-specific integrated policy"),
            "the disjoint-union rejection lacks scheme-policy context: {message}"
        );
    }
}
