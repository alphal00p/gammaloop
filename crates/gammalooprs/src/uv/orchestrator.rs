use color_eyre::Result;
use eyre::{WrapErr, eyre};
use idenso::{
    IndexTooling,
    color::ColorSimplifier,
    shorthands::{metric::MetricSimplifier, schoonschip::Schoonschip},
};
use symbolica::atom::{Atom, AtomCore};
use vakint::Vakint;

use crate::{
    graph::{Graph, cuts::CutSet, feynman_graph::FeynmanGraph},
    model::Model,
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
        model: &Model,
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

        let result = match self {
            Self::LegacyDagForest => legacy_parametric_integrands(
                graph,
                model,
                cut_structure,
                vakint,
                orientation,
                settings,
            ),
            Self::HedgePoset => hedge_poset_parametric_integrands(
                graph,
                model,
                cut_structure,
                vakint,
                orientation,
                settings,
            ),
            Self::Compare => compare_parametric_integrands(
                graph,
                model,
                cut_structure,
                vakint,
                orientation,
                settings,
            ),
        }?;
        let marker = UvMarker::new(settings);
        Ok(result
            .into_iter()
            .map(|integrands| integrands.map(|atom| marker.finish(&atom)))
            .collect())
    }

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

fn legacy_parametric_integrands(
    graph: &mut Graph,
    model: &Model,
    cut_structure: CutStructure,
    vakint: &Vakint,
    orientation: OrientationProjection<'_>,
    settings: &UVgenerationSettings,
) -> Result<Vec<ParametricIntegrands>> {
    let cut_woods = CutWoods::new(cut_structure, graph, model, settings);
    let mut cut_forests = cut_woods.unfold(graph);
    cut_forests.compute(graph, model, vakint, orientation, settings)?;
    cut_forests.orientation_parametric_exprs(graph, settings)
}

fn hedge_poset_parametric_integrands(
    graph: &mut Graph,
    model: &Model,
    cut_structure: CutStructure,
    vakint: &Vakint,
    orientation: OrientationProjection<'_>,
    settings: &UVgenerationSettings,
) -> Result<Vec<ParametricIntegrands>> {
    let wood = HedgePosetWood::new(cut_structure, graph, model, settings);
    let mut forests = wood.unfold();
    forests.compute(graph, model, vakint, orientation, settings)?;
    forests.orientation_parametric_exprs(graph, settings)
}

fn compare_parametric_integrands(
    graph: &mut Graph,
    model: &Model,
    cut_structure: CutStructure,
    vakint: &Vakint,
    orientation: OrientationProjection<'_>,
    settings: &UVgenerationSettings,
) -> Result<Vec<ParametricIntegrands>> {
    let mut hedge_graph = graph.clone();
    let legacy = legacy_parametric_integrands(
        graph,
        model,
        cut_structure.clone(),
        vakint,
        orientation,
        settings,
    )?;
    let hedge = hedge_poset_parametric_integrands(
        &mut hedge_graph,
        model,
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
    model: &Model,
    orientation: OrientationProjection<'_>,
    settings: &UVgenerationSettings,
) -> Result<RenormalizationPart> {
    let mut vk_settings = settings.vakint.true_settings();
    let wood = graph.wood_with_settings(
        &graph.no_dummy(),
        model,
        settings,
        &graph.loop_momentum_basis,
    );
    vk_settings.number_of_terms_in_epsilon_expansion = wood.max_loops as i64;

    let mut forest = wood.unfold(graph, &graph.loop_momentum_basis);
    let vk = (crate::utils::vakint()?, &vk_settings);
    let cuts = CutSet::empty(graph.n_hedges());
    forest.compute(
        graph,
        vk,
        Localizer::new(&cuts, orientation, model),
        settings,
    )?;

    forest.renormalization_part_of_ends(graph, settings)
}

fn hedge_poset_renormalization_part(
    graph: &mut Graph,
    model: &Model,
    settings: &UVgenerationSettings,
) -> Result<RenormalizationPart> {
    let cuts = CutStructure::empty(graph);
    let wood = HedgePosetWood::new(cuts, graph, model, settings);
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
        self.legacy
            .checked_zip(self.hedge, |key, legacy_expr, hedge_expr| {
                if !ComparableExpr::new(legacy_expr).equivalent_to(&ComparableExpr::new(hedge_expr))
                {
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
            .equivalent_to(&ComparableExpr::new(&self.hedge.expression))
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

    fn equivalent_to(&self, other: &Self) -> bool {
        let left = self.normalized();
        let right = other.normalized();

        // Backend-local topology labels may differ on contracted indices.
        left.collect_factors() == right.collect_factors()
            || left.canonize(Aind::Dummy).collect_factors()
                == right.canonize(Aind::Dummy).collect_factors()
    }

    fn normalized(&self) -> Atom {
        self.atom
            .replace(crate::utils::GS.dim)
            .with(4)
            .simplify_metrics()
            .to_dots()
            .simplify_color()
            .expand_num()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use idenso::{bis, gamma};
    use linnet::half_edge::involution::EdgeIndex;
    use spenso::{chain, mink};
    use symbolica::symbol;

    #[test]
    fn compare_canonicalizes_contracted_uv_indices() {
        crate::initialisation::test_initialise().unwrap();
        let expression = |topology| {
            let contracted = mink!(4, Atom::from(Aind::UVTerm(topology, 2)));
            let fixed = mink!(4, Atom::from(Aind::Edge(2, 1)));
            let start = bis!(4, Atom::from(Aind::Hedge(0, 0)));
            let end = bis!(4, Atom::from(Aind::Hedge(1, 0)));
            let common = Atom::var(symbol!("compare_common_factor"));
            let term = |index: Atom| {
                common.clone()
                    * chain!(start.clone(), end.clone(), gamma!(index.clone()))
                    * crate::utils::GS.emr_vec_index(EdgeIndex::from(0), index)
            };
            term(contracted) + term(fixed)
        };
        let legacy = expression(1);
        let hedge = expression(0);
        let legacy = ComparableExpr::new(&legacy);
        let hedge = ComparableExpr::new(&hedge);

        assert_ne!(legacy.normalized(), hedge.normalized());
        assert!(legacy.equivalent_to(&hedge));
    }
}
