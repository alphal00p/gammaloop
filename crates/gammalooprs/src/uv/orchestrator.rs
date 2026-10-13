use color_eyre::Result;
use eyre::{WrapErr, eyre};
use idenso::{
    IndexTooling,
    color::ColorSimplifier,
    shorthands::{metric::MetricSimplifier, schoonschip::Schoonschip},
};
use linnet::half_edge::subgraph::{Inclusion, InternalSubGraph, SubSetLike, SubSetOps};
use symbolica::atom::{Atom, AtomCore};
use vakint::Vakint;

use crate::{
    graph::{Graph, cuts::CutSet, feynman_graph::FeynmanGraph},
    numerator::aind::Aind,
    uv::{
        Integrands, RenormalizationPart, UVOrchestrator, UVgenerationSettings, UltravioletGraph,
        approx::{CutStructure, OrientationProjection, local_3d::Localizer},
        forest::ParametricIntegrands,
        hedge_poset::Wood as HedgePosetWood,
        marker::UvMarker,
        settings::FinalIntegrandDimension,
        wood::CutWoods,
    },
};

impl UVOrchestrator {
    pub(crate) fn parametric_integrands(
        self,
        graph: &mut Graph,
        cut_structure: CutStructure,
        vakint: &Vakint,
        orientation: OrientationProjection<'_>,
        settings: &UVgenerationSettings,
    ) -> Result<Vec<ParametricIntegrands>> {
        if !matches!(settings.final_integrand, FinalIntegrandDimension::ThreeD) {
            return Err(eyre!(
                "4D parametric UV integrands are not supported yet; this mode is planned for a future implementation"
            ));
        }
        validate_local_counterterm_schemes(graph, &cut_structure, settings)?;

        let result = match self {
            Self::LegacyDagForest => {
                legacy_parametric_integrands(graph, cut_structure, vakint, orientation, settings)
            }
            Self::HedgePoset => hedge_poset_parametric_integrands(
                graph,
                cut_structure,
                vakint,
                orientation,
                settings,
            ),
            Self::Compare => {
                compare_parametric_integrands(graph, cut_structure, vakint, orientation, settings)
            }
        }?;
        let marker = UvMarker::new(settings);
        result
            .into_iter()
            .map(|integrands| integrands.map_expressions(|atom| Ok(marker.finish(atom))))
            .collect()
    }

    pub(crate) fn renormalization_part(
        self,
        graph: &mut Graph,
        orientation: OrientationProjection<'_>,
        settings: &UVgenerationSettings,
    ) -> Result<RenormalizationPart> {
        let settings = UVgenerationSettings {
            final_integrand: FinalIntegrandDimension::FourD,
            ..settings.clone()
        };
        validate_local_counterterm_schemes(graph, &CutStructure::empty(graph), &settings)?;
        let mut result = match self {
            Self::LegacyDagForest => legacy_renormalization_part(graph, orientation, &settings),
            Self::HedgePoset => hedge_poset_renormalization_part(graph, &settings),
            Self::Compare => compare_renormalization_part(graph, orientation, &settings),
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
    let direct_three_d = matches!(settings.final_integrand, FinalIntegrandDimension::ThreeD)
        && !settings.local_uv_cts_from_expanded_4d_integrands;

    let mut uncut = graph.full_filter();
    uncut.subtract_with(&graph.initial_state_cut.left);
    for (cut_index, cut) in cuts.cuts.iter().enumerate() {
        let subgraph = uncut.subtract(&cut.union);
        let spinneys = graph
            .classified_spinneys(&subgraph, settings, &graph.loop_momentum_basis)
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
        if direct_three_d
            && os_components.clone().nth(2).is_some()
            && let Some(union) = spinneys.iter().find(|spinney| {
                spinney.n_components() > 2
                    && graph
                        .underlying
                        .connected_components(spinney.filter())
                        .into_iter()
                        .map(|component| {
                            InternalSubGraph::cleaned_filter_optimist(component, graph.as_ref())
                        })
                        .filter(|component| {
                            os_components
                                .clone()
                                .any(|os| os.filter() == &component.filter)
                        })
                        .take(3)
                        .count()
                        > 2
            })
        {
            return Err(eyre!(
                "graph '{}' cut {cut_index}: direct 3D on-shell/OS replay supports at most two disjoint OS components; union {} requires the projected-4D route",
                graph.name,
                union.filter().string_label(),
            ));
        }
        for component in os_components {
            if component.dod != 1 {
                return Err(eyre!(
                    "graph '{}' cut {cut_index}: on-shell/OS component {} requires a degree-one massive-fermion self-energy",
                    graph.name,
                    component.filter().string_label(),
                ));
            }
            let related = spinneys.iter().find(|other| {
                let other = *other;
                if other.n_components() != 1 || other.filter() == component.filter() {
                    return false;
                }
                let nested = component.filter().includes(other.filter())
                    || other.filter().includes(component.filter());
                let disjoint = !component.filter().intersects(other.filter())
                    && spinneys.iter().any(|candidate| {
                        candidate.n_components() > 1
                            && candidate.filter() == &component.filter().union(other.filter())
                    });
                if nested {
                    let (parent, child) = if component.filter().includes(other.filter()) {
                        (component, other)
                    } else {
                        (other, component)
                    };
                    let supported =
                        match (parent.renormalization_scheme, child.renormalization_scheme) {
                            (
                                crate::uv::ApproximationType::OS,
                                crate::uv::ApproximationType::OS,
                            )
                            | (
                                crate::uv::ApproximationType::OS,
                                crate::uv::ApproximationType::IR,
                            ) => true,
                            (
                                crate::uv::ApproximationType::IR,
                                crate::uv::ApproximationType::OS,
                            ) => {
                                parent.dod == 2
                                    && graph.has_qcd_self_energy_insertion_ward_identity(
                                        parent.filter(),
                                        child.filter(),
                                    )
                            }
                            _ => false,
                        };
                    !supported
                } else {
                    disjoint && other.renormalization_scheme != crate::uv::ApproximationType::OS
                }
            });
            if let Some(other) = related {
                return Err(eyre!(
                    "graph '{}' cut {cut_index}: on-shell/OS component {} and component {} require an unsupported nested or disconnected scheme policy; IR-over-OS requires the canonical-QCD self-energy-insertion Ward certificate",
                    graph.name,
                    component.filter().string_label(),
                    other.filter().string_label(),
                ));
            }
            let loops = graph.n_loops(component.filter());
            // A proper multiloop source can replay its contained prefix, but
            // disjoint replay currently reconstructs only two primitive
            // operations. Admitting the nested source does not certify its
            // product with an independent component.
            if direct_three_d
                && loops > 1
                && let Some(other) = spinneys.iter().find(|other| {
                    other.n_components() == 1
                        && other.renormalization_scheme == crate::uv::ApproximationType::OS
                        && !component.filter().intersects(other.filter())
                        && spinneys.iter().any(|union| {
                            union.n_components() > 1
                                && union.filter() == &component.filter().union(other.filter())
                        })
                })
            {
                return Err(eyre!(
                    "graph '{}' cut {cut_index}: direct 3D disjoint on-shell/OS replay requires primitive one-loop components; component {} has {loops} loops and is disjoint from component {}, so use the projected-4D route",
                    graph.name,
                    component.filter().string_label(),
                    other.filter().string_label(),
                ));
            }
            let mut proper_descendants = spinneys.iter().filter(|child| {
                child.n_components() == 1
                    && child.filter() != component.filter()
                    && component.filter().includes(child.filter())
            });
            // A proper source replays its original contour and complete
            // primitive OS prefix before attaching the live cograph. This
            // pure OS history needs no gluon Ward completion; all other
            // mixed-prescription checks above still apply independently.
            let certified_proper_history = loops == 2
                && proper_descendants.next().is_some_and(|child| {
                    child.renormalization_scheme == crate::uv::ApproximationType::OS
                        && child.dod == 1
                        && graph.n_loops(child.filter()) == 1
                })
                && proper_descendants.next().is_none();
            if direct_three_d
                && (loops == 0
                    || (loops > 1
                        && graph.paired_edges(component.filter())
                            != graph.paired_edges(&graph.full_filter())
                        && !certified_proper_history))
            {
                return Err(eyre!(
                    "graph '{}' cut {cut_index}: direct 3D on-shell/OS component {} requires a primitive one-loop component, the complete paired graph domain, or a proper two-loop component with one primitive OS child; other proper histories require the projected-4D route",
                    graph.name,
                    component.filter().string_label(),
                ));
            }
        }
    }

    Ok(())
}

fn legacy_parametric_integrands(
    graph: &mut Graph,
    cut_structure: CutStructure,
    vakint: &Vakint,
    orientation: OrientationProjection<'_>,
    settings: &UVgenerationSettings,
) -> Result<Vec<ParametricIntegrands>> {
    let cut_woods = CutWoods::new(cut_structure, graph, settings)?;
    let mut cut_forests = cut_woods.unfold(graph);
    cut_forests.compute(graph, vakint, orientation, settings)?;
    cut_forests.orientation_parametric_exprs(graph, settings)
}

fn hedge_poset_parametric_integrands(
    graph: &mut Graph,
    cut_structure: CutStructure,
    vakint: &Vakint,
    orientation: OrientationProjection<'_>,
    settings: &UVgenerationSettings,
) -> Result<Vec<ParametricIntegrands>> {
    let wood = HedgePosetWood::new(cut_structure, graph, settings)?;
    let mut forests = wood.unfold();
    forests.compute(graph, vakint, orientation, settings)?;
    forests.orientation_parametric_exprs(graph, settings)
}

fn compare_parametric_integrands(
    graph: &mut Graph,
    cut_structure: CutStructure,
    vakint: &Vakint,
    orientation: OrientationProjection<'_>,
    settings: &UVgenerationSettings,
) -> Result<Vec<ParametricIntegrands>> {
    let mut hedge_graph = graph.clone();
    let legacy =
        legacy_parametric_integrands(graph, cut_structure.clone(), vakint, orientation, settings)?;
    let hedge = hedge_poset_parametric_integrands(
        &mut hedge_graph,
        cut_structure,
        vakint,
        orientation,
        settings,
    )?;

    ParametricIntegrandsComparison {
        legacy: &legacy,
        hedge: &hedge,
    }
    .compare()?;
    Ok(legacy)
}

fn legacy_renormalization_part(
    graph: &mut Graph,
    orientation: OrientationProjection<'_>,
    settings: &UVgenerationSettings,
) -> Result<RenormalizationPart> {
    let mut vk_settings = settings.vakint.true_settings();
    vk_settings.project_onto_tensor_integrals =
        settings.project_integrated_uv_cts_onto_tensor_integrals;
    let wood = graph
        .wood_with_settings(&graph.no_dummy(), settings, &graph.loop_momentum_basis)
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
    )?;

    forest.renormalization_part_of_ends(graph, settings)
}

fn hedge_poset_renormalization_part(
    graph: &mut Graph,
    settings: &UVgenerationSettings,
) -> Result<RenormalizationPart> {
    let cuts = CutStructure::empty(graph);
    let wood = HedgePosetWood::new(cuts, graph, settings)?;
    let mut forest = wood.unfold();
    forest.integrate(graph, crate::utils::vakint()?, settings)?;
    forest.renormalization_part_of_ends(graph, settings)
}

fn compare_renormalization_part(
    graph: &mut Graph,
    orientation: OrientationProjection<'_>,
    settings: &UVgenerationSettings,
) -> Result<RenormalizationPart> {
    let mut hedge_graph = graph.clone();
    let legacy = legacy_renormalization_part(graph, orientation, settings)?;
    let hedge = hedge_poset_renormalization_part(&mut hedge_graph, settings)?;

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
        let legacy = self.legacy.resolved()?;
        let hedge = self.hedge.resolved()?;
        legacy
            .checked_zip(&hedge, |key, legacy_expr, hedge_expr| {
                if !ComparableExpr::new(legacy_expr)
                    .equivalent_to(&ComparableExpr::new(hedge_expr))?
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
        if !ComparableExpr::new(&self.legacy.expression)
            .equivalent_to(&ComparableExpr::new(&self.hedge.expression))?
        {
            return Err(eyre!(
                "UV orchestrator compare mismatch in integrated renormalization part"
            ));
        }
        Ok(())
    }
}

struct ComparableExpr<'a> {
    atom: &'a Atom,
}

impl<'a> ComparableExpr<'a> {
    fn new(atom: &'a Atom) -> Self {
        Self { atom }
    }

    fn equivalent_to(&self, other: &Self) -> Result<bool> {
        let left = self.normalized();
        let right = other.normalized();

        // Backend-local topology labels may differ on contracted indices.
        if left.collect_factors() == right.collect_factors() {
            return Ok(true);
        }

        Ok(left.canonize(Aind::Dummy)?.collect_factors()
            == right.canonize(Aind::Dummy)?.collect_factors())
    }

    fn normalized(&self) -> Atom {
        self.atom
            .replace(crate::utils::GS.dim)
            .with(4)
            .simplify_metrics()
            .to_dots()
            .simplify_color()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        dot,
        graph::parse::IntoGraph,
        initialisation::test_initialise,
        utils::load_generic_model,
        uv::{ApproximationType, CTRenormalizationRule, RenormalizationPrescriptionSettings},
    };
    use idenso::{bis, gamma};
    use linnet::half_edge::subgraph::SubSetLike;
    use spenso::{chain, mink, p};
    use std::collections::BTreeSet;
    use symbolica::symbol;

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
                    * p!(index)
            };
            term(contracted) + term(fixed)
        };
        let legacy = expression(1);
        let hedge = expression(0);
        let legacy = ComparableExpr::new(&legacy);
        let hedge = ComparableExpr::new(&hedge);

        assert_ne!(legacy.normalized(), hedge.normalized());
        for expression in [&legacy, &hedge] {
            assert!(
                expression
                    .normalized()
                    .pattern_match(&spectator.to_pattern(), None, None)
                    .next()
                    .is_some()
            );
        }
        assert!(
            legacy
                .equivalent_to(&hedge)
                .expect("test expressions should canonicalize")
        );
    }

    #[test]
    fn integrated_soft_scheme_is_allowed_before_forest_generation() {
        test_initialise().unwrap();
        let model = load_generic_model("scalars");
        let graph: Graph =
            include_str!("../../../../tests/resources/graphs/scalar/dod2_bubble.dot")
                .into_graph(&model)
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
        let mut graph: Graph =
            include_str!("../../../../tests/resources/graphs/scalar/dod2_bubble.dot")
                .into_graph(&model)
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
            graph.approximation_scheme(&component.filter, &settings, dod),
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
    fn on_shell_scheme_rejects_nonfermionic_degree_before_forest_generation() {
        test_initialise().unwrap();
        let model = load_generic_model("scalars");
        let graph: Graph =
            include_str!("../../../../tests/resources/graphs/scalar/dod2_bubble.dot")
                .into_graph(&model)
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
                .expect_err("a quadratic scalar component cannot use the fermionic OS operator");
        let message = error.to_string();
        assert!(message.contains("on-shell/OS component"));
        assert!(message.contains("degree-one massive-fermion self-energy"));
    }

    #[test]
    fn on_shell_forest_dispatch_requires_certified_mixed_prescriptions() {
        test_initialise().unwrap();
        let model = load_generic_model("sm");
        for source in [
            include_str!(
                "../../../../tests/resources/graphs/paper_appendix_b1_nested_gluon_self_energy.dot"
            ),
            include_str!(
                "../../../../tests/resources/graphs/local_nested_soft_top_self_energy.dot"
            ),
            include_str!(
                "../../../../tests/resources/graphs/local_soft_gluon_in_top_self_energy.dot"
            ),
            include_str!(
                "../../../../tests/resources/graphs/local_disjoint_os_top_self_energies.dot"
            ),
        ] {
            let graph: Graph = source.into_graph(&model).unwrap();
            let baseline = UVgenerationSettings::default();
            let components = graph
                .classified_spinneys(&graph.full_filter(), &baseline, &graph.loop_momentum_basis)
                .unwrap();
            let overrides = components
                .iter()
                .filter(|component| component.n_components() == 1 && component.dod > 0)
                .map(|component| {
                    CTRenormalizationRule::new(
                        graph.ct_identifier(component.filter()),
                        if component.dod == 1 {
                            ApproximationType::OS
                        } else {
                            ApproximationType::IR
                        },
                    )
                })
                .collect();
            let mut settings = UVgenerationSettings {
                generate_integrated: true,
                final_integrand: FinalIntegrandDimension::ThreeD,
                local_uv_cts_from_expanded_4d_integrands: true,
                renormalization_prescription: RenormalizationPrescriptionSettings {
                    overrides,
                    ..Default::default()
                },
                ..baseline
            };
            let selected = graph
                .classified_spinneys(&graph.full_filter(), &settings, &graph.loop_momentum_basis)
                .unwrap();
            assert!(
                selected
                    .iter()
                    .any(|component| component.renormalization_scheme == ApproximationType::OS)
            );
            for projected in [true, false] {
                settings.local_uv_cts_from_expanded_4d_integrands = projected;
                validate_local_counterterm_schemes(&graph, &CutStructure::empty(&graph), &settings)
                    .unwrap_or_else(|error| {
                        panic!("{}, projected={projected}: {error:#}", graph.name)
                    });
            }
        }
    }

    #[test]
    fn direct_on_shell_dispatch_limits_disjoint_os_history_only() {
        test_initialise().unwrap();
        let model = load_generic_model("sm");
        let graph: Graph = r#"digraph three_disjoint_os_components {
            num = "1/3*spenso::g(spenso::dind(spenso::cof(3,gammalooprs::hedge(0))),spenso::cof(3,gammalooprs::hedge(1)))"
            ext [style=invis]
            ext -> left_in:0 [id=0 particle=t]
            right_out:1 -> ext [id=1 particle=t]
            left_in:2 -> left_out:3 [id=2 particle=t lmb_id=0]
            left_in:4 -> left_out:5 [id=3 particle=g]
            left_out:6 -> middle_in:7 [id=4 particle=t]
            middle_in:8 -> middle_out:9 [id=5 particle=t lmb_id=1]
            middle_in:10 -> middle_out:11 [id=6 particle=g]
            middle_out:12 -> right_in:13 [id=7 particle=t]
            right_in:14 -> right_out:15 [id=8 particle=t lmb_id=2]
            right_in:16 -> right_out:17 [id=9 particle=g]
        }"#
        .into_graph(&model)
        .unwrap();
        let mut settings = UVgenerationSettings {
            final_integrand: FinalIntegrandDimension::ThreeD,
            local_uv_cts_from_expanded_4d_integrands: false,
            renormalization_prescription: RenormalizationPrescriptionSettings {
                overrides: vec![CTRenormalizationRule::new(
                    crate::uv::CTIdentifier::new(
                        [-6, 6].into_iter().collect(),
                        Some([6, 21].into_iter().collect()),
                    ),
                    ApproximationType::OS,
                )],
                ..Default::default()
            },
            ..Default::default()
        };
        let cuts = CutStructure::empty(&graph);
        let error = validate_local_counterterm_schemes(&graph, &cuts, &settings)
            .expect_err("three disjoint OS factors exceed the certified direct replay history");
        assert!(
            error
                .to_string()
                .contains("at most two disjoint OS components")
        );

        settings.local_uv_cts_from_expanded_4d_integrands = true;
        validate_local_counterterm_schemes(&graph, &cuts, &settings)
            .expect("the direct replay bound must not change projected dispatch");
        settings.local_uv_cts_from_expanded_4d_integrands = false;
        settings.renormalization_prescription.overrides.clear();
        validate_local_counterterm_schemes(&graph, &cuts, &settings)
            .expect("non-OS woods must retain their existing disjoint support");
    }

    #[test]
    fn integrated_on_shell_self_energy_accepts_all_three_d_routes() {
        test_initialise().unwrap();
        let model = load_generic_model("sm");
        let graph: Graph =
            include_str!("../../../../tests/resources/graphs/local_os_top_self_energy.dot")
                .into_graph(&model)
                .unwrap();
        let mut settings = UVgenerationSettings {
            generate_integrated: true,
            final_integrand: FinalIntegrandDimension::ThreeD,
            local_uv_cts_from_expanded_4d_integrands: true,
            renormalization_prescription: RenormalizationPrescriptionSettings {
                massive_power_divergent: ApproximationType::OS,
                ..Default::default()
            },
            ..Default::default()
        };
        let components = graph
            .classified_spinneys(&graph.full_filter(), &settings, &graph.loop_momentum_basis)
            .unwrap();
        assert!(components.iter().any(|component| {
            component.renormalization_scheme == ApproximationType::OS
                && component.n_components() == 1
        }));
        validate_local_counterterm_schemes(&graph, &CutStructure::empty(&graph), &settings)
            .expect("the primitive OS self-energy has an explicit vacuum integration policy");
        settings.local_uv_cts_from_expanded_4d_integrands = false;
        validate_local_counterterm_schemes(&graph, &CutStructure::empty(&graph), &settings)
            .expect("the complete primitive OS component has validated direct 3D operators");
        settings.final_integrand = FinalIntegrandDimension::FourD;
        validate_local_counterterm_schemes(&graph, &CutStructure::empty(&graph), &settings)
            .expect("symbolic 4D construction does not require a 3D projection route");
    }

    #[test]
    fn direct_on_shell_dispatch_rejects_disjoint_multiloop_component() {
        test_initialise().unwrap();
        // Insert a primitive quark self-energy into the left bubble of the
        // existing disjoint pair. Each component is supported separately,
        // but their product exceeds the primitive disjoint replay history.
        let graph: Graph = r#"digraph disjoint_nested_os_components {
            num = "1/3*spenso::g(spenso::dind(spenso::cof(3,gammalooprs::hedge(0))),spenso::cof(3,gammalooprs::hedge(1)))"
            ext [style=invis]
            ext -> left_in:0 [id=0 particle=t]
            right_out:1 -> ext [id=1 particle=t]
            left_in:2 -> child_in:3 [id=2 particle=t lmb_id=0]
            child_in:4 -> child_out:5 [id=3 particle=t lmb_id=1]
            child_in:6 -> child_out:7 [id=4 particle=g]
            child_out:8 -> left_out:9 [id=5 particle=t]
            left_in:10 -> left_out:11 [id=6 particle=g]
            left_out:12 -> right_in:13 [id=7 particle=t]
            right_in:14 -> right_out:15 [id=8 particle=t lmb_id=2]
            right_in:16 -> right_out:17 [id=9 particle=g]
        }"#
        .into_graph(&load_generic_model("sm"))
        .unwrap();
        let mut settings = UVgenerationSettings {
            generate_integrated: true,
            final_integrand: FinalIntegrandDimension::ThreeD,
            local_uv_cts_from_expanded_4d_integrands: false,
            renormalization_prescription: RenormalizationPrescriptionSettings {
                massive_power_divergent: ApproximationType::OS,
                ..Default::default()
            },
            ..Default::default()
        };
        let classified = graph
            .classified_spinneys(&graph.full_filter(), &settings, &graph.loop_momentum_basis)
            .unwrap();
        let parent = classified
            .iter()
            .find(|node| node.n_components() == 1 && graph.n_loops(node.filter()) == 2)
            .expect("the left bubble is a proper two-loop component");
        let children = classified
            .iter()
            .filter(|node| {
                node.n_components() == 1
                    && node.filter() != parent.filter()
                    && parent.filter().includes(node.filter())
            })
            .collect::<Vec<_>>();
        assert_eq!(children.len(), 1);
        assert_eq!(graph.n_loops(children[0].filter()), 1);
        assert_eq!(children[0].renormalization_scheme, ApproximationType::OS);
        let sibling = classified
            .iter()
            .find(|node| node.n_components() == 1 && !parent.filter().intersects(node.filter()))
            .expect("the right bubble is independent of the nested left component");
        assert_eq!(graph.n_loops(sibling.filter()), 1);
        assert_eq!(parent.renormalization_scheme, ApproximationType::OS);
        assert_eq!(sibling.renormalization_scheme, ApproximationType::OS);
        assert!(classified.iter().any(|node| {
            node.n_components() == 2 && node.filter() == &parent.filter().union(sibling.filter())
        }));
        let cuts = CutStructure::empty(&graph);
        let error = validate_local_counterterm_schemes(&graph, &cuts, &settings)
            .expect_err("a disjoint multiloop OS factor needs its own replay certificate");
        assert!(
            error
                .to_string()
                .contains("disjoint on-shell/OS replay requires primitive one-loop components")
        );
        settings.local_uv_cts_from_expanded_4d_integrands = true;
        validate_local_counterterm_schemes(&graph, &cuts, &settings)
            .expect("the direct replay restriction must not change projected dispatch");
    }

    #[test]
    fn direct_on_shell_dispatch_accepts_only_certified_proper_two_loop_history() {
        test_initialise().unwrap();
        let (graph, _, _, _) =
            crate::uv::approx::local_4d::on_shell_forest_tests::three_loop_quark_fixture().unwrap();
        let baseline = UVgenerationSettings {
            generate_integrated: true,
            final_integrand: FinalIntegrandDimension::ThreeD,
            local_uv_cts_from_expanded_4d_integrands: false,
            renormalization_prescription: RenormalizationPrescriptionSettings {
                massive_power_divergent: ApproximationType::OS,
                ..Default::default()
            },
            ..Default::default()
        };
        let cuts = CutStructure::empty(&graph);
        let classified = graph
            .classified_spinneys(&graph.full_filter(), &baseline, &graph.loop_momentum_basis)
            .unwrap();
        let parent = classified
            .iter()
            .find(|node| node.n_components() == 1 && graph.n_loops(node.filter()) == 2)
            .expect("the real three-loop graph has a proper two-loop component");
        let children = classified
            .iter()
            .filter(|node| {
                node.n_components() == 1
                    && node.filter() != parent.filter()
                    && parent.filter().includes(node.filter())
            })
            .collect::<Vec<_>>();
        assert_eq!(children.len(), 1);
        assert_eq!(graph.n_loops(children[0].filter()), 1);
        assert_eq!(children[0].renormalization_scheme, ApproximationType::OS);
        validate_local_counterterm_schemes(&graph, &cuts, &baseline).expect(
            "bare, completed and finite proper OS sources are certified in both direct routes",
        );

        // CT identifiers contain particle sets, not multiplicities. Give the
        // actual innermost loop a photon so its prescription can be selected
        // independently while the two quark parents still contain gluons.
        let negative_source = include_str!(
            "../../../../tests/resources/graphs/local_os_three_loop_top_self_energy.dot"
        )
        .replace("[id=5 particle=g]", "[id=5 particle=a]");
        let graph: Graph = negative_source
            .as_str()
            .into_graph(&load_generic_model("sm"))
            .unwrap();
        let cuts = CutStructure::empty(&graph);
        let classified = graph
            .classified_spinneys(&graph.full_filter(), &baseline, &graph.loop_momentum_basis)
            .unwrap();
        let child = classified
            .iter()
            .find(|node| node.n_components() == 1 && graph.n_loops(node.filter()) == 1)
            .expect("the photon child is independently classifiable");
        let parent = classified
            .iter()
            .find(|node| node.n_components() == 1 && graph.n_loops(node.filter()) == 2)
            .expect("the proper parent still has two loops");
        assert_ne!(
            graph.ct_identifier(child.filter()),
            graph.ct_identifier(parent.filter())
        );
        let mut settings = baseline.clone();
        settings.renormalization_prescription.overrides = vec![CTRenormalizationRule::new(
            graph.ct_identifier(child.filter()),
            ApproximationType::IR,
        )];
        let selected = graph
            .classified_spinneys(&graph.full_filter(), &settings, &graph.loop_momentum_basis)
            .unwrap();
        assert!(selected.iter().any(|node| node.filter() == child.filter()
            && node.renormalization_scheme == ApproximationType::IR));
        assert!(selected.iter().any(|node| node.filter() == parent.filter()
            && node.renormalization_scheme == ApproximationType::OS));
        let error = validate_local_counterterm_schemes(&graph, &cuts, &settings)
            .expect_err("a proper IR prefix is outside the certified direct history");
        assert!(error.to_string().contains("one primitive OS child"));
        settings.local_uv_cts_from_expanded_4d_integrands = true;
        validate_local_counterterm_schemes(&graph, &cuts, &settings)
            .expect("the direct history limit does not restrict the projected route");
        settings.local_uv_cts_from_expanded_4d_integrands = false;
        settings.renormalization_prescription.overrides = vec![CTRenormalizationRule::new(
            graph.ct_identifier(child.filter()),
            ApproximationType::MUV,
        )];
        let error = validate_local_counterterm_schemes(&graph, &cuts, &settings)
            .expect_err("the existing mixed MUV/OS scheme guard must still run first");
        assert!(
            error
                .to_string()
                .contains("unsupported nested or disconnected scheme policy")
        );
    }

    #[test]
    fn direct_on_shell_dispatch_rejects_multiple_proper_os_descendants() {
        test_initialise().unwrap();
        let graph: Graph = r#"digraph uncertified_proper_os_history {
            num = "1/3*spenso::g(spenso::dind(spenso::cof(3,gammalooprs::hedge(0))),spenso::cof(3,gammalooprs::hedge(1)))"
            ext [style=invis]
            ext -> v1:0 [id=0 particle=t]
            v8:1 -> ext [id=1 particle=t]
            v1:2 -> v2:3 [id=2 particle=t lmb_id=0]
            v2:4 -> v3:5 [id=3 particle=t lmb_id=1]
            v3:6 -> v4:7 [id=4 particle=t lmb_id=2]
            v4:8 -> v5:9 [id=5 particle=t lmb_id=3]
            v4:10 -> v5:11 [id=6 particle=g]
            v5:12 -> v6:13 [id=7 particle=t]
            v3:14 -> v6:15 [id=8 particle=g]
            v6:16 -> v7:17 [id=9 particle=t]
            v2:18 -> v7:19 [id=10 particle=g]
            v7:20 -> v8:21 [id=11 particle=t]
            v1:22 -> v8:23 [id=12 particle=g]
        }"#
        .into_graph(&load_generic_model("sm"))
        .unwrap();
        let mut settings = UVgenerationSettings {
            generate_integrated: false,
            final_integrand: FinalIntegrandDimension::ThreeD,
            local_uv_cts_from_expanded_4d_integrands: false,
            renormalization_prescription: RenormalizationPrescriptionSettings {
                massive_power_divergent: ApproximationType::OS,
                ..Default::default()
            },
            ..Default::default()
        };
        let classified = graph
            .classified_spinneys(&graph.full_filter(), &settings, &graph.loop_momentum_basis)
            .unwrap();
        let parent = classified
            .iter()
            .find(|node| node.n_components() == 1 && graph.n_loops(node.filter()) == 3)
            .expect("the graph has a proper component with two nested descendants");
        assert_eq!(
            classified
                .iter()
                .filter(|node| {
                    node.n_components() == 1
                        && node.filter() != parent.filter()
                        && parent.filter().includes(node.filter())
                })
                .count(),
            2
        );
        let cuts = CutStructure::empty(&graph);
        let error = validate_local_counterterm_schemes(&graph, &cuts, &settings)
            .expect_err("multiple proper prefix steps need their own direct source certificate");
        assert!(error.to_string().contains("one primitive OS child"));
        settings.local_uv_cts_from_expanded_4d_integrands = true;
        validate_local_counterterm_schemes(&graph, &cuts, &settings)
            .expect("the projected route retains its existing forest support");
    }

    #[test]
    fn soft_child_rejects_positive_degree_muv_parent() {
        test_initialise().unwrap();
        let graph: Graph = dot!(
            digraph nested_positive_degree_hu {
                edge [particle=scalar_1 num=1]
                node [num=1]
                ext [style=invis]
                ext -> v1:0 [id=0]
                v4:1 -> ext [id=1]
                v1 -> v2 [id=2 lmb_id=1 num="Q(2,spenso::cind(0))^2"]
                v2 -> v3 [id=3 particle=scalar_2 lmb_id=0 num="Q(3,spenso::cind(0))"]
                v2 -> v3 [id=4 particle=scalar_2]
                v3 -> v4 [id=5]
                v1 -> v4 [id=6]
            },
            "scalars"
        )
        .unwrap();
        let baseline = UVgenerationSettings {
            generate_integrated: false,
            ..Default::default()
        };
        let spinneys = graph
            .classified_spinneys(&graph.full_filter(), &baseline, &graph.loop_momentum_basis)
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
            graph.ct_identifier(child.filter()),
            graph.ct_identifier(parent.filter()),
            "the fixture must allow the child and parent schemes to be selected independently"
        );
        let settings = UVgenerationSettings {
            generate_integrated: false,
            renormalization_prescription: RenormalizationPrescriptionSettings {
                overrides: vec![CTRenormalizationRule::new(
                    graph.ct_identifier(child.filter()),
                    ApproximationType::IR,
                )],
                ..Default::default()
            },
            ..Default::default()
        };
        let classified = graph
            .classified_spinneys(&graph.full_filter(), &settings, &graph.loop_momentum_basis)
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
        let graph: Graph = include_str!("../../../../tests/resources/graphs/dgse.dot")
            .into_graph(&model)
            .unwrap();
        let identifiers = graph
            .spinneys(&graph.full_filter())
            .into_iter()
            .filter(|spinney| {
                graph.compute_dod(&spinney.filter) >= 0
                    && graph.underlying.connected_components(spinney).len() == 1
            })
            .map(|spinney| graph.ct_identifier(&spinney.filter))
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
        let graph: Graph = dot!(
            digraph asymmetric_scalar_spectacles {
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
            },
            "scalars"
        )
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

        let left_identifier = graph.ct_identifier(left);
        let right_identifier = graph.ct_identifier(right);
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
            .classified_spinneys(&graph.full_filter(), &settings, &graph.loop_momentum_basis)
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
