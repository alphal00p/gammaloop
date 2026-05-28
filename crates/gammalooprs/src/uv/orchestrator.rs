use color_eyre::Result;
use eyre::{WrapErr, eyre};
use idenso::{
    color::ColorSimplifier,
    shorthands::{metric::MetricSimplifier, schoonschip::Schoonschip},
};
use symbolica::atom::{Atom, AtomCore};
use vakint::Vakint;

use crate::{
    cff::ResidueSelectedTerms,
    graph::{Graph, cuts::CutSet, feynman_graph::FeynmanGraph},
    uv::{
        RenormalizationPart, UVOrchestrator, UVgenerationSettings, UltravioletGraph,
        approx::{CutStructure, OrientationProjection},
        forest::ParametricIntegrands,
        hedge_poset::Wood as HedgePosetWood,
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
        match self {
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
        }
    }

    pub(crate) fn renormalization_part(
        self,
        graph: &mut Graph,
        orientation: OrientationProjection<'_>,
        settings: &UVgenerationSettings,
    ) -> Result<RenormalizationPart> {
        match self {
            Self::LegacyDagForest => legacy_renormalization_part(graph, orientation, settings),
            Self::HedgePoset => hedge_poset_renormalization_part(graph, settings),
            Self::Compare => compare_renormalization_part(graph, orientation, settings),
        }
    }
}

fn legacy_parametric_integrands(
    graph: &mut Graph,
    cut_structure: CutStructure,
    vakint: &Vakint,
    orientation: OrientationProjection<'_>,
    settings: &UVgenerationSettings,
) -> Result<Vec<ParametricIntegrands>> {
    let cut_woods = CutWoods::new(cut_structure, graph, settings);
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
    let wood = HedgePosetWood::new(cut_structure, graph, settings);
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
    let wood = graph.wood_with_settings(&graph.no_dummy(), settings, &graph.loop_momentum_basis);
    vk_settings.number_of_terms_in_epsilon_expansion = wood.max_loops as i64;

    let mut forest = wood.unfold(graph, &graph.loop_momentum_basis);
    let vk = (crate::utils::vakint()?, &vk_settings);
    let cuts = CutSet::empty(graph.n_hedges());
    forest.compute(graph, vk, &cuts, orientation, settings)?;

    forest.pole_part_of_ends(graph, settings.pole_part)
}

fn hedge_poset_renormalization_part(
    graph: &mut Graph,
    settings: &UVgenerationSettings,
) -> Result<RenormalizationPart> {
    let cuts = CutStructure::empty(graph);
    let wood = HedgePosetWood::new(cuts, graph, settings);
    let mut forest = wood.unfold();
    forest.integrate(graph, crate::utils::vakint()?, settings)?;
    forest.pole_part_of_ends(graph)
}

fn compare_renormalization_part(
    graph: &mut Graph,
    orientation: OrientationProjection<'_>,
    settings: &UVgenerationSettings,
) -> Result<RenormalizationPart> {
    let mut legacy_graph = graph.clone();
    let legacy = legacy_renormalization_part(&mut legacy_graph, orientation, settings)?;
    let hedge = hedge_poset_renormalization_part(graph, settings)?;

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
    legacy: &'a ResidueSelectedTerms,
    hedge: &'a ResidueSelectedTerms,
}

impl IntegrandMapComparison<'_> {
    fn compare(&self) -> Result<()> {
        if self.legacy.keys().ne(self.hedge.keys()) {
            return Err(eyre!(
                "UV orchestrator compare mismatch at cut {}: residue keys differ",
                self.cut_index
            ));
        }

        for (key, legacy_expr) in self.legacy {
            let hedge_expr = self.hedge.get(*key).wrap_err_with(|| {
                format!(
                    "while comparing UV orchestrator outputs at cut {}",
                    self.cut_index
                )
            })?;
            if ComparableExpr::new(legacy_expr).normalized()
                != ComparableExpr::new(hedge_expr).normalized()
            {
                return Err(eyre!(
                    "UV orchestrator compare mismatch at cut {} residue {:?}",
                    self.cut_index,
                    key
                ));
            }
        }

        Ok(())
    }
}

struct RenormalizationComparison<'a> {
    legacy: &'a RenormalizationPart,
    hedge: &'a RenormalizationPart,
}

impl RenormalizationComparison<'_> {
    fn compare(&self) -> Result<()> {
        if ComparableExpr::new(&self.legacy.expression).normalized()
            != ComparableExpr::new(&self.hedge.expression).normalized()
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

    fn normalized(&self) -> Atom {
        self.atom
            .replace(crate::utils::GS.dim)
            .with(4)
            .simplify_metrics()
            .to_dots()
            .simplify_color()
            .expand_num()
            .collect_factors()
    }
}
