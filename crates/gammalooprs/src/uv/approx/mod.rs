use crate::{
    cff::{
        CffEnergyDegreeBoundReport,
        esurface::Esurface,
        expression::{
            GammaLoopOrientationExpression, OrientationExpression, OrientationID,
            OrientationSelector, energy_map_replacements_gs,
        },
        hsurface::Hsurface,
        surface::LinearEnergyExpr,
    },
    debug_tags,
    graph::{Graph, LoopMomentumBasis, cuts::CutSet},
    momentum::Sign,
    settings::global::OrientationPattern,
    utils::GS,
    uv::{
        ApproximationType, Spinney, UVgenerationSettings,
        approx::{
            direct_3d::Direct3dApproximation,
            final_integrand::{FinalIntegrandBuilder, FinalIntegrands},
            integrated::{Integrated, IntegratedCts},
            local_3d::{Local3DApproximation, Local3DCts, Localizer},
            local_4d::{Full4dCts, Local4dCts},
        },
        marker::UvMarker,
        settings::FinalIntegrandDimension,
    },
};
use color_eyre::Result;
use eyre::eyre;
use gammaloop_tracing_filter::{LogMessage, debug_instrument};

use std::{hash::Hash, sync::Mutex};

use symbolica::{
    atom::{Atom, AtomCore, AtomOrView},
    function,
};

use linnet::half_edge::involution::{EdgeIndex, EdgeVec, Orientation};
use linnet::half_edge::subgraph::{InternalSubGraph, SuBitGraph, SubSetLike, SubSetOps};
#[cfg(test)]
use three_dimensional_reps::CffGenerationContext;
use three_dimensional_reps::{Generate3DExpressionOptions, GeneratedThreeDExpression};
use typed_index_collections::TiVec;

use super::IntegrandExpr;
use vakint::Vakint;

#[cfg(test)]
use self::direct_3d::DirectResidueBranches;

pub mod direct_3d;
pub mod final_integrand;
pub mod integrated;
pub mod local_3d;
pub mod local_4d;
mod projected_4d;

pub trait Rooted {
    fn root() -> Self;
}

pub trait ForestNodeLike: LogMessage {
    fn subgraph(&self) -> &SuBitGraph;
    fn lmb(&self) -> &LoopMomentumBasis;
    fn lmb_id(&self) -> EdgeIndex {
        *self.lmb().loop_edges.first().unwrap()
    }
    // fn lmb_given(&self, subgraph: &SuBitGraph) -> &LoopMomentumBasis;
    fn dod(&self) -> i32;
    fn renormalization_scheme(&self) -> ApproximationType;
    fn topo_order(&self) -> usize;
    fn reduced_subgraph(&self, given: &Self) -> SuBitGraph;
}

pub trait ApproximationKernel<C> {
    fn kernel<S: ForestNodeLike>(
        &self,
        ctx: &C,
        current: &S,
        given: &S,
        atom: &Atom,
    ) -> Result<Atom>;
}

pub struct UVCtx<'a> {
    pub graph: &'a Graph,
    pub settings: &'a UVgenerationSettings,
}

impl<'a> UVCtx<'a> {
    pub fn new(graph: &'a Graph, settings: &'a UVgenerationSettings) -> Self {
        Self { graph, settings }
    }
}

pub trait ApproxKernel {
    fn apply<'a, A: Into<AtomOrView<'a>>>(&self, atom: A) -> Atom;
}

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub enum ApproxOp {
    NotComputed,
    // Union operations must be computed before use.
    Union {
        t_args: Vec<IntegrandExpr>,
        subgraphs: Vec<InternalSubGraph>,
    },
    Dependent {
        t_arg: IntegrandExpr,
        subgraph: InternalSubGraph,
    },
    Root,
}

#[derive(Clone)]
pub struct SimpleApprox {
    t_args: Vec<Atom>,
    pub sign: Sign,
    graph: InternalSubGraph,
}

impl SimpleApprox {
    pub(crate) fn expr(&self, bigger_graph: &SuBitGraph) -> Atom {
        let reduced = UvMarker::subgraph(bigger_graph, &self.graph.filter);
        let mut mul = Atom::num(1);
        for i in &self.t_args {
            mul *= i;
        }
        reduced * mul
    }

    pub(crate) fn t_op(&self, bigger_graph: &SuBitGraph) -> Atom {
        function!(GS.uv_approx, self.expr(bigger_graph))
    }

    pub(crate) fn root(subgraph: InternalSubGraph) -> Self {
        if !subgraph.is_empty() {
            panic!(
                "Root approximation must be empty {} {:?}",
                subgraph.string_label(),
                subgraph
            )
        }
        SimpleApprox {
            t_args: vec![],
            sign: Sign::Positive,
            graph: subgraph,
        }
    }

    pub(crate) fn dependent(&self, bigger_graph: InternalSubGraph) -> Self {
        Self {
            t_args: vec![self.t_op(&bigger_graph.filter)],
            sign: -self.sign,
            graph: bigger_graph,
        }
    }
}

#[derive(Clone)]
pub struct Approximation {
    pub spinney: Spinney,
    local_3d: Option<Local3DCts>,
    local: Option<Local4dCts>,
    integrated: Option<IntegratedCts>,
    final_integrand: Option<FinalIntegrands>,
    pub topo_order: usize,
    pub simple_approx: Option<SimpleApprox>,
}

impl Approximation {
    pub(crate) fn integrated(&self, graph: &Graph) -> Result<&IntegratedCts> {
        self.integrated
            .as_ref()
            .ok_or_else(|| eyre!("No integrated CT for {}", self.simple_display(graph)))
    }

    pub(crate) fn local(&self, graph: &Graph) -> Result<&Local4dCts> {
        self.local
            .as_ref()
            .ok_or_else(|| eyre!("No local CT for {}", self.simple_display(graph)))
    }

    pub(crate) fn recursion_input_4d(&self, graph: &Graph) -> Result<Full4dCts> {
        Full4dCts::recursion_input(
            self.local(graph)?,
            self.integrated(graph)?,
            self.renormalization_scheme(),
            self.spinney.subgraph.is_empty(),
            self.lmb(),
        )
    }

    pub(crate) fn local_3d(&self, graph: &Graph) -> Result<&Local3DCts> {
        self.local_3d
            .as_ref()
            .ok_or_else(|| eyre!("No local 3D CT for {}", self.simple_display(graph)))
    }

    pub(crate) fn final_integrand(&self, graph: &Graph) -> Result<&FinalIntegrands> {
        self.final_integrand
            .as_ref()
            .ok_or_else(|| eyre!("No final integrand for {}", self.simple_display(graph)))
    }

    pub fn simple_display(&self, graph: &Graph) -> String {
        format!(
            "{} of {}",
            self.simple_approx
                .as_ref()
                .unwrap()
                .expr(&graph.full_filter()),
            graph.name
        )
    }
}

impl LogMessage for Approximation {
    fn log_display(&self) -> String {
        format!(
            "subgraph={}, topo_order={}, dod={}",
            self.spinney.filter().string_label(),
            self.topo_order,
            self.spinney.dod
        )
    }
}

impl ForestNodeLike for Approximation {
    fn dod(&self) -> i32 {
        self.spinney.dod
    }

    fn renormalization_scheme(&self) -> ApproximationType {
        self.spinney.renormalization_scheme
    }

    fn lmb(&self) -> &LoopMomentumBasis {
        &self.spinney.lmb
    }

    fn reduced_subgraph(&self, given: &Self) -> SuBitGraph {
        self.spinney
            .subgraph
            .subtract(&given.spinney.subgraph)
            .filter
    }

    fn subgraph(&self) -> &SuBitGraph {
        self.spinney.filter()
    }

    fn topo_order(&self) -> usize {
        self.topo_order
    }
}

#[derive(Clone)]
pub struct CutStructure {
    pub cuts: Vec<CutSet>,
}

#[derive(Clone, Copy, Debug)]
enum OrientationProjectionSource<'a> {
    Exact {
        orientations: &'a TiVec<OrientationID, OrientationExpression>,
        root_expression: Option<&'a GeneratedThreeDExpression<Esurface, Hsurface>>,
    },
    Coarse(&'a [EdgeVec<Orientation>]),
}

#[derive(Clone, Copy, Debug)]
pub(crate) struct OrientationProjection<'a> {
    source: OrientationProjectionSource<'a>,
    options: Option<&'a Generate3DExpressionOptions>,
    energy_degree_bound_reports: Option<&'a Mutex<Vec<CffEnergyDegreeBoundReport>>>,
    pub(crate) orientation_pattern: &'a OrientationPattern,
    pub(crate) explicit_orientation_sum_only: bool,
}

impl<'a> OrientationProjection<'a> {
    /// Construct the ordinary coarse-orientation projector used by legacy
    /// exports and isolated UV tests. Production 3D UV uses [`Self::exact`].
    pub(crate) fn new(
        valid_orientations: &'a [EdgeVec<Orientation>],
        orientation_pattern: &'a OrientationPattern,
    ) -> Self {
        Self {
            source: OrientationProjectionSource::Coarse(valid_orientations),
            options: None,
            energy_degree_bound_reports: None,
            orientation_pattern,
            explicit_orientation_sum_only: false,
        }
    }

    #[cfg(test)]
    pub(crate) fn exact(
        orientations: &'a TiVec<OrientationID, OrientationExpression>,
        options: &'a Generate3DExpressionOptions,
        orientation_pattern: &'a OrientationPattern,
        explicit_orientation_sum_only: bool,
    ) -> Self {
        Self {
            source: OrientationProjectionSource::Exact {
                orientations,
                root_expression: None,
            },
            options: Some(options),
            energy_degree_bound_reports: None,
            orientation_pattern,
            explicit_orientation_sum_only,
        }
    }

    pub(crate) fn exact_expression(
        expression: &'a GeneratedThreeDExpression<Esurface, Hsurface>,
        options: &'a Generate3DExpressionOptions,
        orientation_pattern: &'a OrientationPattern,
        explicit_orientation_sum_only: bool,
    ) -> Self {
        Self {
            source: OrientationProjectionSource::Exact {
                orientations: &expression.expression.orientations,
                root_expression: Some(expression),
            },
            options: Some(options),
            energy_degree_bound_reports: None,
            orientation_pattern,
            explicit_orientation_sum_only,
        }
    }

    pub(crate) fn with_energy_degree_bound_reports(
        mut self,
        reports: &'a Mutex<Vec<CffEnergyDegreeBoundReport>>,
    ) -> Self {
        self.energy_degree_bound_reports = Some(reports);
        self
    }

    pub(crate) fn record_energy_degree_bound_report(self, report: &CffEnergyDegreeBoundReport) {
        let Some(reports) = self.energy_degree_bound_reports else {
            return;
        };
        let mut report = report.clone();
        report.physical_parent_bounds.sort_unstable();
        report.assigned_cff_source_bounds.sort_unstable();
        let mut reports = reports
            .lock()
            .unwrap_or_else(|poisoned| poisoned.into_inner());
        if !reports.contains(&report) {
            reports.push(report);
        }
    }

    pub(crate) fn cff_options(self, graph: &Graph) -> Generate3DExpressionOptions {
        self.options
            .cloned()
            .unwrap_or_else(|| graph.denominator_only_cff_3d_expression_options())
    }

    pub(crate) fn orientation_ids(self) -> Vec<OrientationID> {
        match self.source {
            OrientationProjectionSource::Exact { orientations, .. } => orientations
                .iter_enumerated()
                .filter_map(|(id, orientation)| {
                    self.orientation_pattern
                        .filter_orientation(&orientation.data.orientation)
                        .then_some(id)
                })
                .collect(),
            OrientationProjectionSource::Coarse([]) => {
                vec![OrientationID(0)]
            }
            OrientationProjectionSource::Coarse(orientations) => orientations
                .iter()
                .enumerate()
                .filter_map(|(id, orientation)| {
                    self.orientation_pattern
                        .filter_orientation(orientation)
                        .then_some(OrientationID(id))
                })
                .collect(),
        }
    }

    pub(crate) fn exact_orientations(
        self,
    ) -> Option<&'a TiVec<OrientationID, OrientationExpression>> {
        match self.source {
            OrientationProjectionSource::Exact { orientations, .. } => Some(orientations),
            OrientationProjectionSource::Coarse(_) => None,
        }
    }

    pub(crate) fn root_expression(
        self,
    ) -> Option<&'a GeneratedThreeDExpression<Esurface, Hsurface>> {
        match self.source {
            OrientationProjectionSource::Exact {
                root_expression, ..
            } => root_expression,
            OrientationProjectionSource::Coarse(_) => None,
        }
    }

    pub(crate) fn orientation(self, id: OrientationID) -> Option<&'a EdgeVec<Orientation>> {
        match self.source {
            OrientationProjectionSource::Exact { orientations, .. } => orientations
                .get(id)
                .map(|orientation| &orientation.data.orientation),
            OrientationProjectionSource::Coarse(orientations) => orientations.get(id.0),
        }
    }

    /// Apply one branch-owned exact energy map only to the newly owned
    /// numerator fragment. The source map, when present, is authoritative;
    /// otherwise the production selector also owns the numerator map.
    pub(crate) fn map_numerator(
        self,
        graph: &Graph,
        selector_id: OrientationID,
        source_edge_energy_map: Option<&[LinearEnergyExpr]>,
        numerator: &Atom,
    ) -> Result<Atom> {
        match self.source {
            OrientationProjectionSource::Exact { orientations, .. } => {
                let orientation = orientations.get(selector_id).ok_or_else(|| {
                    eyre!(
                        "missing production energy map for orientation {}",
                        selector_id.0
                    )
                })?;
                let replacements = source_edge_energy_map.map_or_else(
                    || orientation.energy_replacements_gs(graph),
                    |edge_energy_map| energy_map_replacements_gs(edge_energy_map, graph),
                );
                Ok(numerator.replace_multiple(replacements))
            }
            OrientationProjectionSource::Coarse(_) if source_edge_energy_map.is_some() => {
                Err(eyre!(
                    "a source-local exact energy map cannot be carried by a coarse orientation projector"
                ))
            }
            OrientationProjectionSource::Coarse(_) => Ok(numerator.clone()),
        }
    }
}

impl CutStructure {
    pub(crate) fn empty(graph: &Graph) -> Self {
        Self {
            cuts: vec![CutSet::empty(graph.n_hedges())],
        }
    }
}

impl Approximation {
    pub(crate) fn root(
        &mut self,
        graph: &mut Graph,
        localizer: Localizer<'_>,
        settings: &UVgenerationSettings,
    ) -> Result<()> {
        self.simple_approx = Some(SimpleApprox::root(self.spinney.subgraph.clone()));
        self.local = Some(Local4dCts::root());
        let integrated = IntegratedCts::root();
        if let FinalIntegrandDimension::ThreeD = settings.final_integrand {
            // The root is the original factorized 3D integrand, not a local UV
            // approximation. Expanded-4D projection applies only to proper UV
            // nodes and therefore must not change this identity element.
            let local_3d = Local3DCts::root(graph, localizer)?;
            self.final_integrand = Some(FinalIntegrandBuilder::new(localizer, settings).build_3d(
                graph,
                self,
                &local_3d,
                &integrated,
            )?);
            self.local_3d = Some(local_3d);
        }
        self.integrated = Some(integrated);

        Ok(())
    }

    pub(crate) fn new(spinney: Spinney) -> Approximation {
        Approximation {
            spinney,
            topo_order: 0,
            final_integrand: None,
            simple_approx: None,
            local: None,
            local_3d: None,
            integrated: None,
        }
    }

    #[debug_instrument(
        graph = %graph.log_display(),
        current = %self.log_display(),
        given = %dependent.log_display(),
        reduced = ?self.reduced_subgraph(dependent),
    )]
    pub(crate) fn compute_4d(
        &mut self,
        graph: &Graph,
        vakint: (&Vakint, &vakint::VakintSettings),
        dependent: &Self,
        settings: &UVgenerationSettings,
    ) -> Result<()> {
        let ctx = UVCtx { graph, settings };
        debug_tags!(#generation,#uv,#fourd;
            simple = %self.simple_display(graph),
            "Computing 4D",
        );

        let old_full = dependent.recursion_input_4d(graph)?;
        let local = local_4d::uv_limit(&old_full, &ctx, self, dependent, self, dependent)?;
        let integrated = if settings.generate_integrated {
            Integrated::new(vakint.0, vakint.1)
                .run(&local, &ctx, self, dependent, self, dependent)?
        } else {
            IntegratedCts::root()
        };

        self.local = Some(local);
        self.integrated = Some(integrated);

        Ok(())
    }

    /// Computes the 3d approximation of the UV
    #[allow(clippy::too_many_arguments)]
    #[debug_instrument(
        graph = %graph.log_display(),
        current = %self.log_display(),
        given = %dependent.log_display(),
    )]
    pub(crate) fn compute_3d(
        &mut self,
        dependent: &Self,
        graph: &mut Graph,
        localizer: Localizer<'_>,
        settings: &UVgenerationSettings,
    ) -> Result<()> {
        let local_3d = if settings.local_uv_cts_from_expanded_4d_integrands {
            let local_4d = self.local(graph)?;
            Local3DApproximation::new(localizer, graph, settings)
                .project_local_4d(local_4d, self)?
        } else {
            let parent_local = dependent.local_3d(graph)?;
            let parent_integrated = dependent.integrated(graph)?;
            Local3DCts::Direct(Direct3dApproximation::new(localizer, graph, settings).run(
                parent_local.direct()?,
                parent_integrated,
                self,
                dependent,
                self,
                dependent,
            )?)
        };

        let integrated = self.integrated(graph)?;

        self.final_integrand = Some(
            FinalIntegrandBuilder::new(localizer, settings)
                .build_3d(graph, self, &local_3d, integrated)?,
        );
        self.local_3d = Some(local_3d);
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use std::collections::{BTreeMap, BTreeSet};

    use super::*;
    use crate::{
        cff::CutCFFIndex,
        dot,
        graph::{FeynmanGraph, cuts::LuCutSelection, parse::IntoGraph},
        initialisation::test_initialise,
        integrands::{
            evaluation::EvaluationMetaData,
            process::{
                GenericEvaluatorFloat,
                evaluators::{EvaluatorStack, SingleOrAllOrientations},
            },
        },
        momentum::{
            ThreeMomentum,
            sample::{BareMomentumSample, ExternalFourMomenta, LoopMomenta, MomentumSample},
        },
        processes::{EvaluatorSettings, cross_section::build_derivative_structure},
        settings::RuntimeSettings,
        utils::{
            ArbPrec, F, FloatLike, W_, cut_energy,
            hyperdual_utils::{DualOrNot, extract_t_derivatives_complex, simple_n_deriv_shape},
        },
        uv::UltravioletGraph,
    };
    use idenso::shorthands::{metric::MetricSimplifier, schoonschip::Schoonschip};
    use linnet::half_edge::subgraph::subset::SubSet;
    use spenso::{
        algebra::{algebraic_traits::IsZero, complex::Complex},
        structure::representation::{Minkowski, RepName},
    };
    use symbolica::domains::dual::HyperDual;
    use symbolica::{
        atom::AtomView,
        domains::{
            float::{Complex as SymComplex, Real},
            integer::IntegerRing,
            rational::{Fraction, Rational},
        },
        evaluate::ExpressionEvaluator,
        parse_lit,
    };
    use three_dimensional_reps::CffEnergyFactorOwnership;

    #[test]
    fn expanded_4d_setting_does_not_change_the_empty_forest_root() -> Result<()> {
        test_initialise()?;
        let graph: Graph = dot!(digraph root_identity {
            edge [num=1 mass=1]
            node [num=1]

            a -> b [id=0 lmb_id=0]
            a -> b [id=1]
        })?;
        let mut direct_graph = graph.clone();
        let mut projected_graph = graph;
        let cutset = CutSet::empty(direct_graph.n_hedges());
        let orientation_pattern = OrientationPattern::default();
        let localizer = Localizer::new(
            &cutset,
            OrientationProjection::new(&[], &orientation_pattern),
        );
        let mut direct = Approximation::new(Spinney::empty(&direct_graph));
        let mut projected = Approximation::new(Spinney::empty(&projected_graph));
        let direct_settings = UVgenerationSettings::default();
        let projected_settings = UVgenerationSettings {
            local_uv_cts_from_expanded_4d_integrands: true,
            ..direct_settings.clone()
        };

        direct.root(&mut direct_graph, localizer, &direct_settings)?;
        projected.root(&mut projected_graph, localizer, &projected_settings)?;

        assert_eq!(projected.local_3d, direct.local_3d);
        assert_eq!(projected.final_integrand, direct.final_integrand);
        Ok(())
    }

    #[test]
    fn selected_raised_lu_projection_preserves_quadratic_contact_family() -> Result<()> {
        test_initialise()?;
        let mut graph: Graph = dot!(digraph raised_lu_contact_family {
            num = 1
            edge [particle="scalar_1" num=1]
            node [num=1]

            ext0 [style=invis is_cut=0]
            v2 -> ext0 [id=0]
            ext0 -> v3
            v0 -> v1 [id=1 lmb_id=0]
            v0 -> v1 [id=2]
            v3 -> v0 [id=3 lmb_id=1]
            v1 -> v2 [id=4]
            v2 -> v3 [id=5]
        }, "scalars")?;
        let options = graph.denominator_only_cff_3d_expression_options();
        let production_numerator =
            GS.emr_mom(EdgeIndex(1), GS.cind(0)) * GS.emr_mom(EdgeIndex(3), GS.cind(0));
        let canonization = graph.get_esurface_canonization(&graph.loop_momentum_basis);
        let production = graph.generate_3d_expression_for_integrand(
            &[],
            &canonization,
            &options,
            Some(&production_numerator),
        )?;
        let lu_cut = graph
            .determine_raised_esurfaces_from_expression(&production.expression)
            .raised_groups
            .into_iter()
            .find(|group| {
                group.max_occurence == 2
                    && group.esurface_ids.iter().any(|surface_id| {
                        !production.expression.surfaces.esurface_cache[*surface_id]
                            .external_shift
                            .is_empty()
                    })
            })
            .expect("the repeated outer channel supplies a raised LU surface");
        let mut cutset = CutSet::empty(graph.n_hedges());
        cutset.residue_selector.lu = Some(LuCutSelection {
            cut_edge_alternatives: lu_cut
                .esurface_ids
                .iter()
                .map(|surface_id| {
                    production.expression.surfaces.esurface_cache[*surface_id]
                        .energies
                        .clone()
                })
                .collect(),
            raised_group: lu_cut,
        });
        let contracted = graph
            .get_edge_subgraph(EdgeIndex(1))
            .union(&graph.get_edge_subgraph(EdgeIndex(2)));
        let analysis_numerator = GS.emr_mom(EdgeIndex(3), GS.cind(0))
            * GS.emr_mom(EdgeIndex(4), GS.cind(0))
            + Atom::num(Rational::from((7, 11)));
        let contract_subgraph = contracted
            .union(&graph.tree_edges)
            .subtract(&graph.initial_state_cut);
        let contract_edges = graph
            .iter_edges_of(&contract_subgraph)
            .filter_map(|(pair, edge, _)| pair.is_paired().then_some(edge))
            .collect::<Vec<_>>();
        let contracted_generated = graph.generate_3d_expression_for_integrand(
            &contract_edges,
            &canonization,
            &options,
            Some(&analysis_numerator),
        )?;
        assert_eq!(
            contracted_generated.energy_factor_ownership,
            CffEnergyFactorOwnership::VariantLocal,
            "the quadratic contact family owns its on-shell factors per generalized residue-map variant",
        );
        let mut raw_graph = graph.clone();
        let raw = raw_graph.cff(
            &contract_subgraph,
            &cutset,
            &OrientationPattern::default(),
            &options,
            Some(&analysis_numerator),
        )?;
        let lu1 = raw
            .terms
            .iter()
            .find(|(index, _)| index.lu_cut_order == Some(1))
            .expect("the raw CFF has an LU-order-one term")
            .1;
        let contact_samples = lu1
            .orientations
            .iter()
            .flat_map(|orientation| {
                orientation
                    .orientation
                    .variants
                    .iter()
                    .filter(|variant| variant.numerator_surfaces.is_empty())
                    .map(|variant| {
                        (
                            orientation.orientation.edge_energy_map[3].clone(),
                            variant.prefactor.clone(),
                        )
                    })
            })
            .collect::<Vec<_>>();
        let expected_contact_samples = [
            (LinearEnergyExpr::ose(EdgeIndex(3), -1), Atom::num(2)),
            (LinearEnergyExpr::zero(), Atom::num(-4)),
            (LinearEnergyExpr::ose(EdgeIndex(3), 1), Atom::num(2)),
        ];
        assert!(
            contact_samples.len() == expected_contact_samples.len()
                && expected_contact_samples
                    .iter()
                    .all(|sample| contact_samples.contains(sample)),
            "the selected LU1 residue needs its complete minus/zero/plus quadratic sampling family"
        );
        assert!(
            raw.terms
                .iter()
                .find(|(index, _)| index.lu_cut_order == Some(2))
                .expect("the raw CFF has an LU-order-two term")
                .1
                .orientations
                .iter()
                .flat_map(|orientation| &orientation.orientation.variants)
                .all(|variant| !variant.numerator_surfaces.is_empty()),
            "LU2 is the ordinary physical residue and has no quadratic contact family"
        );

        let orientation_pattern = OrientationPattern::default();
        let localizer = Localizer::new(
            &cutset,
            OrientationProjection::exact_expression(
                &production,
                &options,
                &orientation_pattern,
                true,
            ),
        );
        let projected = localizer.projected_cff(
            &mut graph,
            &contracted,
            &analysis_numerator,
            CffGenerationContext::Standalone,
        )?;
        let projected_lu1_samples = projected
            .iter_orientations()
            .filter(|(_, _, integrands)| {
                integrands
                    .iter()
                    .any(|(index, atom)| index.lu_cut_order == Some(1) && !atom.is_zero())
            })
            .filter_map(|(_, map, _)| map.map(|map| map[3].clone()))
            .fold(Vec::new(), |mut distinct, sample| {
                if !distinct.contains(&sample) {
                    distinct.push(sample);
                }
                distinct
            });
        let expected_projected_samples = [
            LinearEnergyExpr::ose(EdgeIndex(3), -1),
            LinearEnergyExpr::zero(),
            LinearEnergyExpr::ose(EdgeIndex(3), 1),
        ];
        assert!(
            projected_lu1_samples.len() == expected_projected_samples.len()
                && expected_projected_samples
                    .iter()
                    .all(|sample| projected_lu1_samples.contains(sample)),
            "production hosting must preserve the three distinct minus/zero/plus members of the selected LU1 sampling family; projected {projected_lu1_samples:?}"
        );

        let raw_expressions = raw
            .terms
            .iter()
            .map(|(index, term)| {
                (
                    *index,
                    term.orientations
                        .iter()
                        .map(|orientation| {
                            &orientation.expression
                                * analysis_numerator.replace_multiple(
                                    orientation.orientation.energy_replacements_gs(&raw_graph),
                                )
                        })
                        .fold(Atom::Zero, |sum, term| sum + term)
                        * Atom::num(raw.production_prefactor_factor()),
                )
            })
            .collect::<BTreeMap<_, _>>();
        let mut denominator_only_graph = graph.clone();
        let denominator_only = denominator_only_graph.cff(
            &contract_subgraph,
            &cutset,
            &OrientationPattern::default(),
            &options,
            Some(&Atom::one()),
        )?;
        let denominator_only_expressions = denominator_only
            .terms
            .iter()
            .map(|(index, term)| {
                (
                    *index,
                    term.orientations
                        .iter()
                        .map(|orientation| {
                            &orientation.expression
                                * analysis_numerator.replace_multiple(
                                    orientation
                                        .orientation
                                        .energy_replacements_gs(&denominator_only_graph),
                                )
                        })
                        .fold(Atom::Zero, |sum, term| sum + term)
                        * Atom::num(denominator_only.production_prefactor_factor()),
                )
            })
            .collect::<BTreeMap<_, _>>();
        let mut projected_expressions = BTreeMap::<CutCFFIndex, Atom>::new();
        for (orientation_id, source_edge_energy_map, integrands) in projected.iter_orientations() {
            let mapped_numerator = localizer.map_numerator(
                &graph,
                orientation_id,
                source_edge_energy_map,
                &analysis_numerator,
            )?;
            for (index, expression) in integrands.iter() {
                *projected_expressions.entry(*index).or_default() += expression * &mapped_numerator;
            }
        }
        let fixed_point = |mut expression: Atom| -> Result<Atom> {
            expression = expression
                .replace(GS.dim)
                .with(4)
                .replace(GS.dim_epsilon)
                .with(0)
                .replace(GS.numerator_sampling_scale)
                .with(Atom::num(Rational::from((13, 10))))
                .replace(function!(GS.tree_denom_wrapper, W_.x_))
                .with(W_.x_)
                .replace(GS.den(W_.a_, W_.b_, W_.c_, W_.d_))
                .with(W_.d_);
            for edge in 0..graph.underlying.n_edges() {
                let edge = EdgeIndex(edge);
                expression = expression
                    .replace(graph.underlying[edge].particle.mass_atom())
                    .with(Atom::one())
                    .replace(GS.ose(edge))
                    .with(Atom::num(Rational::from((5, 4))))
                    .replace(cut_energy(edge))
                    .with(Atom::num(Rational::from((5, 4))))
                    .replace(GS.emr_mom(edge, GS.cind(0)))
                    .with(Atom::num(Rational::from((17, 7))));
            }
            Ok(expression)
        };
        let evaluate_arb = |expression: Atom| -> Result<Complex<F<ArbPrec>>> {
            let parameters = [Atom::var(GS.pi)];
            let rational: ExpressionEvaluator<SymComplex<Fraction<IntegerRing>>> = expression
                .evaluator(&parameters)
                .build()
                .map_err(|error| eyre!("failed to build raised-contact evaluator: {error}"))?;
            let mut arb: ExpressionEvaluator<Complex<F<ArbPrec>>> =
                rational.map_coeff(&|coefficient| {
                    Complex::new(F::from(&coefficient.re), F::from(&coefficient.im))
                });
            let zero = F(ArbPrec::default());
            Ok(arb.evaluate_single(&[Complex::new(zero.clone().pi(), zero)]))
        };
        let tolerance = F(ArbPrec::default().epsilon()).sqrt().sqrt().sqrt();
        for order in [1, 2] {
            let index = raw_expressions
                .keys()
                .find(|index| index.lu_cut_order == Some(order))
                .copied()
                .expect("both raised orders are present");
            let raw_value = evaluate_arb(fixed_point(raw_expressions[&index].clone())?)?;
            let denominator_only_value =
                evaluate_arb(fixed_point(denominator_only_expressions[&index].clone())?)?;
            let projected_value =
                evaluate_arb(fixed_point(projected_expressions[&index].clone())?)?;
            let scale = raw_value
                .clone()
                .norm()
                .re
                .max(denominator_only_value.clone().norm().re)
                .max(projected_value.clone().norm().re);
            assert!(
                (raw_value.clone() - &denominator_only_value).norm().re <= &tolerance * &scale,
                "the generalized and denominator-only selected LU{order} residues differ: generalized={raw_value}, denominator-only={denominator_only_value}"
            );
            assert!(
                (raw_value - projected_value).norm().re <= &tolerance * scale,
                "production hosting changed the complete selected LU{order} residue"
            );
        }
        Ok(())
    }

    #[test]
    fn gl24_direct_3d_modes_preserve_orientation_selector_contracts() -> Result<()> {
        test_initialise()?;
        let mut graph: Graph = dot!(digraph gl24_selector_contract {
            edge [num=1 mass=1]
            node [num=1]

            v0 -> v1 [id=0 is_cut=0]
            v0 -> v3 [id=1 lmb_id=0 num="(Q(1,spenso::cind(0))^2-Q(1,spenso::cind(1))^2-Q(1,spenso::cind(2))^2-Q(1,spenso::cind(3))^2)^2"]
            v0 -> v5 [id=2]
            v1 -> v2 [id=3 lmb_id=1]
            v1 -> v3 [id=4]
            v2 -> v4 [id=5 lmb_id=2]
            v2 -> v4 [id=6]
            v3 -> v5 [id=7]
            v4 -> v5 [id=8 mass=2]
        })?;
        let options = graph.denominator_only_cff_3d_expression_options();
        let numerator = graph.production_numerator_atom_for_full_3d_expression();
        let canonization = graph.get_esurface_canonization(&graph.loop_momentum_basis);
        let production = graph.generate_3d_expression_for_integrand(
            &[],
            &canonization,
            &options,
            Some(&numerator),
        )?;

        let uv_filter = graph
            .get_edge_subgraph(EdgeIndex(5))
            .union(&graph.get_edge_subgraph(EdgeIndex(6)));
        let uv_subgraph = InternalSubGraph::cleaned_filter_optimist(uv_filter, graph.as_ref());
        let child_spinney = Spinney::new(uv_subgraph, &graph, &graph.loop_momentum_basis)
            .expect("the GL24 e5/e6 bubble has a compatible nonempty UV sub-LMB");
        assert!(!child_spinney.subgraph.is_empty());

        let target_support = [EdgeIndex(3), EdgeIndex(4)]
            .into_iter()
            .collect::<BTreeSet<_>>();
        let lu_surface = graph
            .determine_raised_esurfaces_from_expression(&production.expression)
            .raised_groups
            .into_iter()
            .find(|group| {
                group.esurface_ids.iter().any(|surface_id| {
                    production.expression.surfaces.esurface_cache[*surface_id]
                        .energies
                        .iter()
                        .copied()
                        .collect::<BTreeSet<_>>()
                        == target_support
                })
            })
            .ok_or_else(|| eyre!("GL24 q3/q4 LU surface was not generated"))?;
        let mut cutset = CutSet::empty(graph.n_hedges());
        cutset.residue_selector.lu = Some(LuCutSelection {
            raised_group: lu_surface,
            cut_edge_alternatives: vec![vec![EdgeIndex(3), EdgeIndex(4)]],
        });
        let expected_keys = cutset
            .residue_selector
            .generate_allowed_keys()
            .into_iter()
            .collect::<BTreeSet<_>>();
        assert_eq!(
            expected_keys
                .iter()
                .filter_map(|index| index.lu_cut_order)
                .collect::<BTreeSet<_>>(),
            BTreeSet::from([1]),
        );

        let orientation_pattern = OrientationPattern::default();
        let settings = UVgenerationSettings {
            generate_integrated: false,
            local_uv_cts_from_expanded_4d_integrands: false,
            add_marker: false,
            ..Default::default()
        };
        let check = |explicit_orientation_sum_only| -> Result<(usize, Atom)> {
            let mut route_graph = graph.clone();
            let localizer = Localizer::new(
                &cutset,
                OrientationProjection::exact_expression(
                    &production,
                    &options,
                    &orientation_pattern,
                    explicit_orientation_sum_only,
                ),
            );
            let mut root = Approximation::new(Spinney::empty(&route_graph));
            root.root(&mut route_graph, localizer, &settings)?;
            let mut child = Approximation::new(child_spinney.clone());
            child.simple_approx = Some(
                root.simple_approx
                    .as_ref()
                    .expect("the root approximation is initialized")
                    .dependent(child.spinney.subgraph.clone()),
            );
            let vakint_settings = vakint::VakintSettings::default();
            child.compute_4d(
                &route_graph,
                (crate::utils::vakint()?, &vakint_settings),
                &root,
                &settings,
            )?;
            child.compute_3d(&root, &mut route_graph, localizer, &settings)?;
            let sectors = child.local_3d(&route_graph)?.direct()?.sectors()?;
            let [sector] = sectors else {
                return Err(eyre!(
                    "the isolated GL24 bubble must produce exactly one complete-CFF sector"
                ));
            };
            // Both direct modes Taylor-expand the same complete CFF. Exact
            // residue-map keys stay in the sparse branch sidecar while T acts;
            // neither branch body contains a physical-theta selector.
            let integrands = sector.combine()?;

            let mut nonzero_bodies = 0;
            for (key, branch) in integrands.iter_keys() {
                assert_eq!(
                    branch
                        .iter()
                        .map(|(index, _)| *index)
                        .collect::<BTreeSet<_>>(),
                    expected_keys
                );
                for (_, body) in branch.iter().filter(|(_, body)| !body.is_zero()) {
                    nonzero_bodies += 1;
                    assert!(
                        !body.contains_symbol(GS.theta)
                            && !body.contains_symbol(OrientationID::symbol()),
                        "a direct Taylor body must retain its complete map key in the sparse sidecar"
                    );
                }
                assert!(
                    usize::from(key.selector_host.0) < production.expression.orientations.len(),
                    "the sparse branch selector host must be a production residue-map key"
                );
            }
            let final_integrand = child
                .final_integrand(&route_graph)?
                .iter()
                .find(|(index, _)| expected_keys.contains(index))
                .ok_or_else(|| eyre!("the selected GL24 residue has no final integrand"))?
                .1;
            Ok((nonzero_bodies, final_integrand))
        };
        let (ordinary_count, ordinary_integrand) = check(false)?;
        let (explicit_count, explicit_integrand) = check(true)?;
        assert!(ordinary_count > 0);
        assert!(explicit_count > 0);

        let production_map_ids = production
            .expression
            .orientations
            .iter_enumerated()
            .map(|(id, _)| id)
            .collect::<Vec<_>>();
        let ordinary_sum = production_map_ids
            .iter()
            .map(|id| id.select(ordinary_integrand.as_view()))
            .fold(Atom::Zero, |sum, atom| sum + atom);
        assert!(
            (ordinary_sum - explicit_integrand).expand().is_zero(),
            "summing the orientation-local final integrand must recover the explicit final integrand"
        );

        Ok(())
    }

    #[test]
    fn nested_scalar_bubble_direct_3d_modes_match_without_a_cut() -> Result<()> {
        test_initialise()?;
        let mut graph: Graph = dot!(
            digraph nested_scalar_bubble_selector_contract {
                edge [particle="scalar_1"];
                node [num=1];

                v1 -> v2 [id=0 lmb_id=0];
                v1 -> v2 [id=1];
                v2 -> v3 [id=2 lmb_id=1];
                v3 -> v1 [id=3];
            },
            "scalars"
        )?;
        let options = graph.denominator_only_cff_3d_expression_options();
        let numerator = graph.production_numerator_atom_for_full_3d_expression();
        let canonization = graph.get_esurface_canonization(&graph.loop_momentum_basis);
        let production = graph.generate_3d_expression_for_integrand(
            &[],
            &canonization,
            &options,
            Some(&numerator),
        )?;
        let uv_filter = graph
            .get_edge_subgraph(EdgeIndex(0))
            .union(&graph.get_edge_subgraph(EdgeIndex(1)));
        let uv_subgraph = InternalSubGraph::cleaned_filter_optimist(uv_filter, graph.as_ref());
        let child_spinney = Spinney::new(uv_subgraph, &graph, &graph.loop_momentum_basis)
            .expect("the doubled edge is a logarithmic one-loop UV subgraph");
        let orientation_pattern = OrientationPattern::default();
        let settings = UVgenerationSettings {
            generate_integrated: false,
            add_marker: false,
            ..Default::default()
        };
        let build = |mut route_graph: Graph,
                     explicit_orientation_sum_only|
         -> Result<(
            Vec<OrientationID>,
            Vec<(OrientationID, Atom)>,
            Vec<(OrientationID, Atom)>,
            Atom,
        )> {
            let cutset = CutSet::empty(route_graph.n_hedges());
            let localizer = Localizer::new(
                &cutset,
                OrientationProjection::exact_expression(
                    &production,
                    &options,
                    &orientation_pattern,
                    explicit_orientation_sum_only,
                ),
            );
            let mut root = Approximation::new(Spinney::empty(&route_graph));
            root.root(&mut route_graph, localizer, &settings)?;
            let mut child = Approximation::new(child_spinney.clone());
            child.simple_approx = Some(
                root.simple_approx
                    .as_ref()
                    .expect("the root approximation is initialized")
                    .dependent(child.spinney.subgraph.clone()),
            );
            child.compute_4d(
                &route_graph,
                (crate::utils::vakint()?, &vakint::VakintSettings::default()),
                &root,
                &settings,
            )?;
            child.compute_3d(&root, &mut route_graph, localizer, &settings)?;
            let sectors = child.local_3d(&route_graph)?.direct()?.sectors()?;
            let mut hosts = sectors
                .iter()
                .flat_map(|sector| sector.active.iter_keys())
                .filter(|(_, integrands)| integrands.iter().any(|(_, atom)| !atom.is_zero()))
                .map(|(key, _)| key.selector_host)
                .collect::<Vec<_>>();
            hosts.sort_unstable();
            hosts.dedup();
            let summarize_direct = |integrands: &DirectResidueBranches| {
                let mut by_host = BTreeMap::<OrientationID, Atom>::new();
                for (key, branch) in integrands.iter_keys() {
                    *by_host.entry(key.selector_host).or_default() +=
                        branch.iter().fold(Atom::Zero, |sum, (_, atom)| sum + atom);
                }
                by_host.into_iter().collect::<Vec<_>>()
            };
            let summarize_projected = |integrands: &local_3d::OrientationIntegrands| {
                let mut by_host = BTreeMap::<OrientationID, Atom>::new();
                for (host, _, branch) in integrands.iter_orientations() {
                    *by_host.entry(host).or_default() +=
                        branch.iter().fold(Atom::Zero, |sum, (_, atom)| sum + atom);
                }
                by_host.into_iter().collect::<Vec<_>>()
            };
            let child_branches = sectors
                .iter()
                .flat_map(|sector| summarize_direct(&sector.active))
                .collect::<Vec<_>>();
            let mut outer_branches = Vec::new();
            if !explicit_orientation_sum_only {
                for sector in sectors {
                    let reduced = route_graph
                        .full_filter()
                        .subtract(child.subgraph())
                        .subtract(&route_graph.initial_state_cut);
                    let resnum = route_graph
                        .numerator(&reduced, child.subgraph())
                        .get_single_atom()
                        .expect("the scalar outer numerator is available")
                        * route_graph.global_atom();
                    let analysis_numerator = sector.active.factorized_capacity_envelope() * &resnum;
                    let outer = localizer.projected_cff(
                        &mut route_graph,
                        child.subgraph(),
                        &analysis_numerator,
                        CffGenerationContext::Standalone,
                    )?;
                    outer_branches.extend(summarize_projected(&outer));
                }
            }
            let final_integrand = child
                .final_integrand(&route_graph)?
                .iter()
                .find(|(index, _)| *index == CutCFFIndex::new_all_none())
                .ok_or_else(|| eyre!("the scalar-bubble child has no root residue"))?
                .1
                .clone();
            Ok((hosts, child_branches, outer_branches, final_integrand))
        };
        let build_full_cff_oracle = |mut route_graph: Graph| -> Result<Atom> {
            let cutset = CutSet::empty(route_graph.n_hedges());
            let localizer = Localizer::new(
                &cutset,
                OrientationProjection::exact_expression(
                    &production,
                    &options,
                    &orientation_pattern,
                    false,
                ),
            );
            let mut root = Approximation::new(Spinney::empty(&route_graph));
            root.root(&mut route_graph, localizer, &settings)?;
            let child = Approximation::new(child_spinney.clone());
            let root_local = root.local_3d(&route_graph)?.direct()?.clone();
            let local = Local3DCts::Direct(
                Direct3dApproximation::new(localizer, &mut route_graph, &settings).run(
                    &root_local,
                    &IntegratedCts::root(),
                    &child,
                    &root,
                    &child,
                    &root,
                )?,
            );
            Ok(FinalIntegrandBuilder::new(localizer, &settings)
                .build_3d(&mut route_graph, &child, &local, &IntegratedCts::root())?
                .iter()
                .find(|(index, _)| *index == CutCFFIndex::new_all_none())
                .ok_or_else(|| eyre!("the full-CFF selector oracle has no root residue"))?
                .1)
        };
        let (ordinary_hosts, child_branches, outer_branches, ordinary) =
            build(graph.clone(), false)?;
        let (_, _, _, explicit) = build(graph.clone(), true)?;
        let full_cff_oracle = build_full_cff_oracle(graph.clone())?;

        let production_map_ids = production
            .expression
            .orientations
            .iter_enumerated()
            .map(|(id, _)| id)
            .collect::<Vec<_>>();
        let ordinary_sum = production_map_ids
            .iter()
            .map(|id| id.select(ordinary.as_view()))
            .fold(Atom::Zero, |sum, atom| sum + atom);
        let full_cff_sum = production_map_ids
            .iter()
            .map(|id| id.select(full_cff_oracle.as_view()))
            .fold(Atom::Zero, |sum, atom| sum + atom);
        assert!(
            (&full_cff_sum - &explicit).expand().is_zero(),
            "the established full-CFF selector oracle must recover the explicit child sum"
        );
        let mismatched_map_keys = production_map_ids
            .iter()
            .filter_map(|id| {
                let ordinary = id.select(ordinary.as_view());
                let oracle = id.select(full_cff_oracle.as_view());
                let difference = (&ordinary - &oracle).expand();
                (!difference.is_zero()).then_some((
                    *id,
                    production.expression.orientations[*id]
                        .data
                        .orientation
                        .clone(),
                    ordinary,
                    oracle,
                    difference,
                ))
            })
            .collect::<Vec<_>>();
        assert!(
            (ordinary_sum - explicit).expand().is_zero(),
            "the isolated scalar-bubble child differs between direct local-3D modes; ordinary hosts {ordinary_hosts:?}, child branches {child_branches:?}, outer branches {outer_branches:?}, mismatched production map keys {mismatched_map_keys:?}"
        );

        Ok(())
    }

    #[test]
    fn factorized_owned_dot_child_cff_matches_direct_3d_for_uncut_self_energy() -> Result<()> {
        test_initialise()?;
        let mut graph: Graph = dot!(digraph factorized_child_cff {
            edge [num=1 mass=1]
            node [num=1]
            incoming [style=invis]
            outgoing [style=invis]

            incoming -> v3 [id=0]
            v0 -> v1 [id=1 lmb_id=0 num="Q(1,spenso::mink(4,1))"]
            v0 -> v1 [id=2 num="Q(2,spenso::mink(4,1))"]
            v3 -> v0 [id=3 lmb_id=1]
            v1 -> v2 [id=4]
            v2 -> v3 [id=5]
            v2 -> outgoing [id=6]
        })?;
        let options = graph.denominator_only_cff_3d_expression_options();
        let numerator = graph.production_numerator_atom_for_full_3d_expression();
        let canonization = graph.get_esurface_canonization(&graph.loop_momentum_basis);
        let production = graph.generate_3d_expression_for_integrand(
            &[],
            &canonization,
            &options,
            Some(&numerator),
        )?;
        let uv_filter = graph
            .get_edge_subgraph(EdgeIndex(1))
            .union(&graph.get_edge_subgraph(EdgeIndex(2)));
        let uv_subgraph = InternalSubGraph::cleaned_filter_optimist(uv_filter, graph.as_ref());
        let child_spinney = Spinney::new(uv_subgraph, &graph, &graph.loop_momentum_basis)
            .expect("the self-energy bubble has a compatible sub-LMB");
        assert_eq!(child_spinney.dod, 2);
        let cutset = CutSet::empty(graph.n_hedges());
        let orientation_pattern = OrientationPattern::default();
        let build = |mut route_graph: Graph,
                     from_expanded_4d|
         -> Result<(Atom, Vec<CffEnergyDegreeBoundReport>, Vec<Vec<usize>>)> {
            let settings = UVgenerationSettings {
                generate_integrated: false,
                local_uv_cts_from_expanded_4d_integrands: from_expanded_4d,
                ..Default::default()
            };
            let bound_reports = Mutex::new(Vec::new());
            let localizer = Localizer::new(
                &cutset,
                OrientationProjection::exact_expression(
                    &production,
                    &options,
                    &orientation_pattern,
                    true,
                )
                .with_energy_degree_bound_reports(&bound_reports),
            );
            let mut root = Approximation::new(Spinney::empty(&route_graph));
            root.root(&mut route_graph, localizer, &settings)?;
            let mut child = Approximation::new(child_spinney.clone());
            child.simple_approx = Some(
                root.simple_approx
                    .as_ref()
                    .expect("the root approximation is initialized")
                    .dependent(child.spinney.subgraph.clone()),
            );
            let vakint_settings = vakint::VakintSettings::default();
            child.compute_4d(
                &route_graph,
                (crate::utils::vakint()?, &vakint_settings),
                &root,
                &settings,
            )?;
            let projected_denominator_owners = if from_expanded_4d {
                child
                    .local(&route_graph)?
                    .active_sectors()
                    .iter()
                    .map(|sector| sector.physical_terms())
                    .collect::<Result<Vec<_>>>()?
                    .into_iter()
                    .flatten()
                    .map(|term| {
                        term.denominators
                            .into_iter()
                            .map(|denominator| usize::from(denominator.source_edge))
                            .collect()
                    })
                    .collect()
            } else {
                Vec::new()
            };
            child.compute_3d(&root, &mut route_graph, localizer, &settings)?;
            match (from_expanded_4d, child.local_3d(&route_graph)?) {
                (false, Local3DCts::Direct(_)) | (true, Local3DCts::Projected4d(_)) => {}
                (false, _) => {
                    return Err(eyre!(
                        "the direct self-energy must retain complete-CFF sectors"
                    ));
                }
                (true, _) => {
                    return Err(eyre!(
                        "the projected self-energy must retain factorized local-4D coefficients"
                    ));
                }
            }
            let integrand = child
                .final_integrand(&route_graph)?
                .iter()
                .find(|(index, _)| *index == CutCFFIndex::new_all_none())
                .expect("the uncut self-energy has one residue sector")
                .1;
            let reports = bound_reports
                .lock()
                .unwrap_or_else(|poisoned| poisoned.into_inner())
                .clone();
            Ok((integrand, reports, projected_denominator_owners))
        };
        let (direct_expression, direct_reports, _) = build(graph.clone(), false)?;
        let (projected_expression, projected_reports, projected_denominator_owners) =
            build(graph.clone(), true)?;
        assert_eq!(projected_denominator_owners, vec![vec![1, 1, 2, 2, 2]]);

        let direct_bound = direct_reports
            .iter()
            .find(|report| {
                report.source_kind == crate::cff::CffEnergyBoundSourceKind::PhysicalGraph
                    && report.physical_parent_bounds == vec![(1, 1), (2, 1)]
            })
            .expect("the direct child CFF must retain its owned-dot energy bounds");
        assert_eq!(
            direct_bound.assigned_cff_source_bounds,
            vec![(1, 1), (2, 1)]
        );
        let projected_bound = projected_reports
            .iter()
            .find(|report| {
                report.source_kind == crate::cff::CffEnergyBoundSourceKind::ExactFourD
                    && report.physical_parent_bounds == vec![(1, 3), (2, 5)]
            })
            .expect("the projected child CFF must report its rank-eight exact source");
        let mut occurrence_degrees = projected_bound
            .assigned_cff_source_bounds
            .iter()
            .map(|(_, degree)| *degree)
            .collect::<Vec<_>>();
        occurrence_degrees.sort_unstable();
        assert_eq!(occurrence_degrees, vec![1, 1, 2, 2, 2]);

        // Component expansion is confined to these test copies. Production
        // keeps the owned-dot numerator factorized throughout construction.
        let normalize = |expression: Atom| {
            expression
                .replace(GS.dim)
                .with(4)
                .replace(GS.dim_epsilon)
                .with(0)
                .replace(GS.m_uv_expansion)
                .with(Atom::num(Rational::from((7, 5))))
                .replace(GS.m_uv_vacuum)
                .with(Atom::num(Rational::from((7, 5))))
                .replace(function!(GS.tree_denom_wrapper, W_.x_))
                .with(W_.x_)
                .replace(GS.den(W_.a_, W_.b_, W_.c_, W_.d_))
                .with(W_.d_)
                .replace(function!(GS.ose, W_.mass_, W_.prop_))
                .with(W_.prop_)
                .expand_dots()
                .expect("test-only component expansion must succeed")
        };
        let expressions = [
            normalize(direct_expression),
            normalize(projected_expression),
        ];
        assert_eq!(graph.loop_momentum_basis.loop_edges.len(), 2);
        assert_eq!(graph.loop_momentum_basis.ext_edges.len(), 2);
        let rational =
            |numerator, denominator| F::<ArbPrec>::from(&Rational::from((numerator, denominator)));
        let loop_moms: LoopMomenta<F<ArbPrec>> = [
            ThreeMomentum::new(rational(31, 100), rational(-47, 100), rational(83, 100)),
            ThreeMomentum::new(rational(-19, 100), rational(37, 100), rational(61, 100)),
        ]
        .into_iter()
        .collect();
        let external_moms: ExternalFourMomenta<F<ArbPrec>> = [
            [
                rational(12, 5),
                rational(0, 1),
                rational(0, 1),
                rational(0, 1),
            ]
            .into(),
            [
                rational(-12, 5),
                rational(0, 1),
                rational(0, 1),
                rational(0, 1),
            ]
            .into(),
        ]
        .into_iter()
        .collect();
        let sample = MomentumSample {
            sample: BareMomentumSample {
                loop_moms,
                dual_loop_moms: None,
                loop_mom_cache_id: 0,
                loop_mom_base_cache_id: 0,
                external_moms,
                external_mom_cache_id: 0,
                external_mom_base_cache_id: 0,
                jacobian: rational(1, 1),
                orientation: None,
                parameterization_branch: None,
            },
        };
        let orientations = TiVec::<OrientationID, EdgeVec<Orientation>>::new();
        let orientation_filter = SubSet::full(orientations.len());
        let mut param_builder = graph.param_builder.clone();
        let (mut evaluator, _) = EvaluatorStack::new_explicit_sum_with_timings(
            &expressions,
            &param_builder,
            None,
            &EvaluatorSettings::default(),
        )?;
        let input = <ArbPrec as GenericEvaluatorFloat>::get_parameters(
            &mut param_builder,
            (false, false),
            &graph,
            &sample,
            &[],
            &[],
            None,
            None,
            None,
        );
        let values = evaluator
            .evaluate(
                input,
                SingleOrAllOrientations::All {
                    all: &orientations,
                    filter: &orientation_filter,
                },
                &RuntimeSettings::default(),
                &mut EvaluationMetaData::new_empty(),
                false,
            )?
            .into_iter()
            .map(|value| value.unwrap_real())
            .collect::<Vec<_>>();
        let distance = (values[0].clone() - values[1].clone()).norm().re;
        let direct_norm = values[0].clone().norm().re;
        let projected_norm = values[1].clone().norm().re;
        let scale = if direct_norm > projected_norm {
            direct_norm
        } else {
            projected_norm
        };
        let relative_distance = if scale.is_zero() {
            distance
        } else {
            distance / scale
        };
        let tolerance = F(ArbPrec::default().epsilon()).sqrt().sqrt().sqrt();
        assert!(
            relative_distance <= tolerance,
            "kinematically consistent owned-dot self-energy mismatch: direct={:e}, projected={:e}, relative delta={relative_distance:e}, tolerance={tolerance:e}",
            values[0],
            values[1],
        );
        Ok(())
    }

    #[test]
    fn nested_scalar_banana_direct_3d_matches_typed_local_4d_with_arb() -> Result<()> {
        test_initialise()?;
        let mut graph: Graph = dot!(digraph nested_scalar_banana {
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
        let global_marker = symbolica::symbol!("nested_scalar_banana::global");
        graph.global_prefactor.num = Atom::var(global_marker);

        let inner_filter = graph
            .get_edge_subgraph(EdgeIndex(1))
            .union(&graph.get_edge_subgraph(EdgeIndex(2)));
        let outer_filter = inner_filter.union(&graph.get_edge_subgraph(EdgeIndex(3)));
        let inner_subgraph =
            InternalSubGraph::cleaned_filter_optimist(inner_filter, graph.as_ref());
        let outer_subgraph =
            InternalSubGraph::cleaned_filter_optimist(outer_filter, graph.as_ref());
        let inner_spinney = Spinney::new(inner_subgraph, &graph, &graph.loop_momentum_basis)
            .expect("the two-edge scalar bubble has a compatible sub-LMB");
        let outer_spinney = Spinney::new(outer_subgraph, &graph, &graph.loop_momentum_basis)
            .expect("the three-edge scalar banana has a compatible sub-LMB");

        let options = graph.denominator_only_cff_3d_expression_options();
        let numerator = graph.production_numerator_atom_for_full_3d_expression();
        let canonization = graph.get_esurface_canonization(&graph.loop_momentum_basis);
        let production = graph.generate_3d_expression_for_integrand(
            &[],
            &canonization,
            &options,
            Some(&numerator),
        )?;
        let cutset = CutSet::empty(graph.n_hedges());
        let orientation_pattern = OrientationPattern::default();

        let build = |mut route_graph: Graph, from_expanded_4d| -> Result<(Atom, Atom)> {
            let settings = UVgenerationSettings {
                generate_integrated: false,
                local_uv_cts_from_expanded_4d_integrands: from_expanded_4d,
                ..Default::default()
            };
            let localizer = Localizer::new(
                &cutset,
                OrientationProjection::exact_expression(
                    &production,
                    &options,
                    &orientation_pattern,
                    true,
                ),
            );
            let vakint_settings = vakint::VakintSettings::default();
            let mut root = Approximation::new(Spinney::empty(&route_graph));
            root.root(&mut route_graph, localizer, &settings)?;

            let mut inner = Approximation::new(inner_spinney.clone());
            inner.simple_approx = Some(
                root.simple_approx
                    .as_ref()
                    .expect("the root approximation is initialized")
                    .dependent(inner.spinney.subgraph.clone()),
            );
            inner.compute_4d(
                &route_graph,
                (crate::utils::vakint()?, &vakint_settings),
                &root,
                &settings,
            )?;
            inner.compute_3d(&root, &mut route_graph, localizer, &settings)?;
            match (from_expanded_4d, inner.local_3d(&route_graph)?) {
                (false, Local3DCts::Direct(direct)) => {
                    let [sector] = direct.sectors()? else {
                        return Err(eyre!(
                            "the direct inner bubble must retain one complete-CFF sector"
                        ));
                    };
                    assert_eq!(sector.active_subgraph, *inner.subgraph());
                    assert_eq!(route_graph.n_loops(&sector.active_subgraph), 1);
                    assert!(
                        sector
                            .frozen_integrands
                            .iter()
                            .all(|(_, atom)| atom.is_one()),
                        "an unintegrated direct CFF sector has no frozen localization factor"
                    );
                    assert!(
                        sector.active.iter_keys().any(|(_, integrands)| integrands
                            .iter()
                            .any(|(_, atom)| !atom.is_zero())),
                        "the complete inner CFF Taylor coefficient must be nonzero"
                    );
                }
                (true, Local3DCts::Projected4d(sectors)) => {
                    assert!(
                        !sectors.is_empty(),
                        "the projected inner bubble must retain its typed local-4D sectors"
                    );
                }
                (false, _) => {
                    return Err(eyre!(
                        "the direct inner bubble must use complete-CFF Taylor sectors"
                    ));
                }
                (true, _) => {
                    return Err(eyre!(
                        "the local-4D route must retain projected typed sectors"
                    ));
                }
            }
            let inner_final = inner
                .final_integrand(&route_graph)?
                .iter()
                .find(|(index, _)| *index == CutCFFIndex::new_all_none())
                .ok_or_else(|| eyre!("the inner bubble has no uncut final integrand"))?
                .1
                .clone();

            let mut outer = Approximation::new(outer_spinney.clone());
            outer.simple_approx = Some(
                inner
                    .simple_approx
                    .as_ref()
                    .expect("the inner approximation is initialized")
                    .dependent(outer.spinney.subgraph.clone()),
            );
            outer.compute_4d(
                &route_graph,
                (crate::utils::vakint()?, &vakint_settings),
                &inner,
                &settings,
            )?;
            outer.compute_3d(&inner, &mut route_graph, localizer, &settings)?;
            match (from_expanded_4d, outer.local_3d(&route_graph)?) {
                (false, Local3DCts::Direct(direct)) => {
                    let [sector] = direct.sectors()? else {
                        return Err(eyre!(
                            "the nested direct banana must retain one complete-CFF sector"
                        ));
                    };
                    assert_eq!(sector.active_subgraph, *outer.subgraph());
                    assert_eq!(route_graph.n_loops(&sector.active_subgraph), 2);
                    assert!(
                        sector
                            .frozen_integrands
                            .iter()
                            .all(|(_, atom)| atom.is_one())
                    );
                    assert!(
                        !sector
                            .active
                            .factorized_sum()
                            .contains_symbol(global_marker),
                        "the untouched global numerator must remain outside the Taylor-expanded CFF"
                    );
                }
                (true, Local3DCts::Projected4d(sectors)) => {
                    assert!(
                        !sectors.is_empty(),
                        "the nested projected banana must retain typed local-4D sectors"
                    );
                }
                (false, _) => {
                    return Err(eyre!(
                        "the nested direct banana must use complete-CFF Taylor sectors"
                    ));
                }
                (true, _) => {
                    return Err(eyre!(
                        "the nested local-4D route must retain projected typed sectors"
                    ));
                }
            }
            let outer_final = outer
                .final_integrand(&route_graph)?
                .iter()
                .find(|(index, _)| *index == CutCFFIndex::new_all_none())
                .ok_or_else(|| eyre!("the nested banana has no uncut final integrand"))?
                .1
                .clone();
            Ok((inner_final, outer_final))
        };

        let (direct_inner, direct_outer) = build(graph.clone(), false)?;
        let (projected_inner, projected_outer) = build(graph.clone(), true)?;
        // Test-only component expansion makes the independently constructed
        // expressions evaluable by one Arb stack. Production keeps all
        // numerator factors intact throughout CFF and UV construction.
        let normalize = |expression: Atom| {
            expression
                .replace(global_marker)
                .with(1)
                .replace(GS.dim)
                .with(4)
                .replace(GS.dim_epsilon)
                .with(0)
                .replace(GS.m_uv_expansion)
                .with(Atom::num(Rational::from((7, 5))))
                .replace(GS.m_uv_vacuum)
                .with(Atom::num(Rational::from((7, 5))))
                .replace(function!(GS.tree_denom_wrapper, W_.x_))
                .with(W_.x_)
                .replace(GS.den(W_.a_, W_.b_, W_.c_, W_.d_))
                .with(W_.d_)
                .replace(function!(GS.ose, W_.mass_, W_.prop_))
                .with(W_.prop_)
                .expand_dots()
                .expect("test-only component expansion must succeed")
        };
        let expressions = [
            normalize(direct_inner),
            normalize(projected_inner),
            normalize(direct_outer),
            normalize(projected_outer),
        ];

        let rational =
            |numerator, denominator| F::<ArbPrec>::from(&Rational::from((numerator, denominator)));
        let loop_moms: LoopMomenta<F<ArbPrec>> = [
            ThreeMomentum::new(rational(31, 100), rational(-47, 100), rational(83, 100)),
            ThreeMomentum::new(rational(-19, 100), rational(37, 100), rational(61, 100)),
        ]
        .into_iter()
        .collect();
        let external_moms: ExternalFourMomenta<F<ArbPrec>> = [
            [
                rational(12, 5),
                rational(11, 100),
                rational(-8, 100),
                rational(5, 100),
            ]
            .into(),
            [
                rational(-12, 5),
                rational(-11, 100),
                rational(8, 100),
                rational(-5, 100),
            ]
            .into(),
        ]
        .into_iter()
        .collect();
        let sample = MomentumSample {
            sample: BareMomentumSample {
                loop_moms,
                dual_loop_moms: None,
                loop_mom_cache_id: 0,
                loop_mom_base_cache_id: 0,
                external_moms,
                external_mom_cache_id: 0,
                external_mom_base_cache_id: 0,
                jacobian: rational(1, 1),
                orientation: None,
                parameterization_branch: None,
            },
        };
        let orientations = TiVec::<OrientationID, EdgeVec<Orientation>>::new();
        let orientation_filter = SubSet::full(orientations.len());
        let mut param_builder = graph.param_builder.clone();
        let (mut evaluator, _) = EvaluatorStack::new_explicit_sum_with_timings(
            &expressions,
            &param_builder,
            None,
            &EvaluatorSettings::default(),
        )?;
        let input = <ArbPrec as GenericEvaluatorFloat>::get_parameters(
            &mut param_builder,
            (false, false),
            &graph,
            &sample,
            &[],
            &[],
            None,
            None,
            None,
        );
        let values = evaluator
            .evaluate(
                input,
                SingleOrAllOrientations::All {
                    all: &orientations,
                    filter: &orientation_filter,
                },
                &RuntimeSettings::default(),
                &mut EvaluationMetaData::new_empty(),
                false,
            )?
            .into_iter()
            .map(|value| value.unwrap_real())
            .collect::<Vec<_>>();
        let tolerance = F(ArbPrec::default().epsilon()).sqrt().sqrt().sqrt();
        for (label, direct, projected) in [
            ("inner bubble", &values[0], &values[1]),
            ("nested banana", &values[2], &values[3]),
        ] {
            let distance = (direct.clone() - projected.clone()).norm().re;
            let direct_norm = direct.clone().norm().re;
            let projected_norm = projected.clone().norm().re;
            let scale = if direct_norm > projected_norm {
                direct_norm
            } else {
                projected_norm
            };
            let relative_distance = if scale.is_zero() {
                distance
            } else {
                distance / scale
            };
            assert!(
                relative_distance <= tolerance,
                "{label} direct complete-CFF vs projected local-4D mismatch: direct={direct:e}, projected={projected:e}, relative delta={relative_distance:e}, tolerance={tolerance:e}"
            );
        }
        Ok(())
    }

    #[test]
    fn direct_nested_scalar_banana_replays_complete_cff_at_depth_three() -> Result<()> {
        test_initialise()?;
        let mut graph: Graph = dot!(digraph depth_three_scalar_banana {
            edge [num=1 mass=1]
            node [num=1]
            incoming [style=invis]
            outgoing [style=invis]

            incoming -> a [id=0]
            a -> b [id=1 lmb_id=0]
            a -> b [id=2 lmb_id=1]
            a -> b [id=3 lmb_id=2]
            a -> b [id=4]
            b -> outgoing [id=5]
        })?;
        let global_marker = symbolica::symbol!("depth_three_scalar_banana::global");
        graph.global_prefactor.num = Atom::var(global_marker);

        let inner_filter = graph
            .get_edge_subgraph(EdgeIndex(1))
            .union(&graph.get_edge_subgraph(EdgeIndex(2)));
        let middle_filter = inner_filter.union(&graph.get_edge_subgraph(EdgeIndex(3)));
        let outer_filter = middle_filter.union(&graph.get_edge_subgraph(EdgeIndex(4)));
        let inner_subgraph =
            InternalSubGraph::cleaned_filter_optimist(inner_filter, graph.as_ref());
        let middle_subgraph =
            InternalSubGraph::cleaned_filter_optimist(middle_filter, graph.as_ref());
        let outer_subgraph =
            InternalSubGraph::cleaned_filter_optimist(outer_filter, graph.as_ref());
        let inner_spinney = Spinney::new(inner_subgraph, &graph, &graph.loop_momentum_basis)
            .expect("the inner bubble has a compatible sub-LMB");
        let middle_spinney = Spinney::new(middle_subgraph, &graph, &graph.loop_momentum_basis)
            .expect("the middle banana has a compatible sub-LMB");
        let outer_spinney = Spinney::new(outer_subgraph, &graph, &graph.loop_momentum_basis)
            .expect("the outer banana has a compatible sub-LMB");

        let options = graph.denominator_only_cff_3d_expression_options();
        let numerator = graph.production_numerator_atom_for_full_3d_expression();
        let canonization = graph.get_esurface_canonization(&graph.loop_momentum_basis);
        let production = graph.generate_3d_expression_for_integrand(
            &[],
            &canonization,
            &options,
            Some(&numerator),
        )?;
        let cutset = CutSet::empty(graph.n_hedges());
        let orientation_pattern = OrientationPattern::default();
        let localizer = Localizer::new(
            &cutset,
            OrientationProjection::exact_expression(
                &production,
                &options,
                &orientation_pattern,
                true,
            ),
        );
        let settings = UVgenerationSettings {
            generate_integrated: false,
            local_uv_cts_from_expanded_4d_integrands: false,
            ..Default::default()
        };
        let vakint_settings = vakint::VakintSettings::default();
        let mut root = Approximation::new(Spinney::empty(&graph));
        root.root(&mut graph, localizer, &settings)?;

        let assert_complete_cff = |graph: &Graph,
                                   approximation: &Approximation,
                                   expected_loops: usize,
                                   label: &str|
         -> Result<()> {
            let sectors = approximation.local_3d(graph)?.direct()?.sectors()?;
            let [sector] = sectors else {
                return Err(eyre!(
                    "{label} must have one complete-CFF sector when integrated generation is disabled"
                ));
            };
            assert_eq!(
                sector.active_subgraph,
                *approximation.subgraph(),
                "{label} must carry the complete nested UV subgraph"
            );
            assert_eq!(
                graph.n_loops(&sector.active_subgraph),
                expected_loops,
                "{label} must retain every loop already reached by the nested replay"
            );
            assert!(
                sector
                    .frozen_integrands
                    .iter()
                    .all(|(_, atom)| atom.is_one()),
                "{label} has no integrated prefix to freeze"
            );
            assert!(
                sector
                    .active
                    .iter_keys()
                    .any(|(_, integrands)| integrands.iter().any(|(_, atom)| !atom.is_zero())),
                "{label} must contain a nonzero Taylor coefficient of the complete CFF"
            );
            assert!(
                sector.active.iter_keys().all(|(_, integrands)| integrands
                    .iter()
                    .all(|(index, _)| *index == CutCFFIndex::new_all_none())),
                "{label} must preserve the uncut CFF residue key through every Taylor replay"
            );
            let combined = sector.combine()?.factorized_sum();
            assert_eq!(
                approximation
                    .local_3d(graph)?
                    .direct()?
                    .branches()?
                    .factorized_sum(),
                combined,
                "{label} materialization must be exactly its stored complete-CFF sector"
            );
            assert!(
                !combined.contains_symbol(global_marker),
                "{label} must leave the untouched global numerator outside the UV Taylor operation"
            );
            Ok(())
        };

        let mut inner = Approximation::new(inner_spinney);
        inner.simple_approx = Some(
            root.simple_approx
                .as_ref()
                .expect("the root approximation is initialized")
                .dependent(inner.spinney.subgraph.clone()),
        );
        inner.compute_4d(
            &graph,
            (crate::utils::vakint()?, &vakint_settings),
            &root,
            &settings,
        )?;
        inner.compute_3d(&root, &mut graph, localizer, &settings)?;
        assert_complete_cff(&graph, &inner, 1, "inner bubble")?;

        let mut middle = Approximation::new(middle_spinney);
        middle.simple_approx = Some(
            inner
                .simple_approx
                .as_ref()
                .expect("the inner approximation is initialized")
                .dependent(middle.spinney.subgraph.clone()),
        );
        middle.compute_4d(
            &graph,
            (crate::utils::vakint()?, &vakint_settings),
            &inner,
            &settings,
        )?;
        middle.compute_3d(&inner, &mut graph, localizer, &settings)?;
        assert_complete_cff(&graph, &middle, 2, "middle banana")?;

        let mut outer = Approximation::new(outer_spinney);
        outer.simple_approx = Some(
            middle
                .simple_approx
                .as_ref()
                .expect("the middle approximation is initialized")
                .dependent(outer.spinney.subgraph.clone()),
        );
        outer.compute_4d(
            &graph,
            (crate::utils::vakint()?, &vakint_settings),
            &middle,
            &settings,
        )?;
        outer.compute_3d(&middle, &mut graph, localizer, &settings)?;
        assert_complete_cff(&graph, &outer, 3, "outer banana")?;

        let original = outer.local_3d(&graph)?;
        let original_sum = original.direct()?.branches()?.factorized_sum();
        let transformed = (-original.clone()).map(|atom| Ok(atom.clone()))?;
        let Local3DCts::Direct(transformed_direct) = &transformed else {
            return Err(eyre!(
                "negation and atom mapping must preserve direct complete-CFF sectors"
            ));
        };
        let transformed_sectors = transformed_direct.sectors()?;
        let original_sectors = original.direct()?.sectors()?;
        assert_eq!(
            transformed_sectors
                .iter()
                .map(|sector| &sector.active_subgraph)
                .collect::<Vec<_>>(),
            original_sectors
                .iter()
                .map(|sector| &sector.active_subgraph)
                .collect::<Vec<_>>(),
            "algebraic maps must preserve the nested complete-CFF active subgraphs"
        );
        assert!(
            (transformed.direct()?.branches()?.factorized_sum() + original_sum)
                .expand()
                .is_zero(),
            "mapping a negated complete-CFF replay must negate its full materialized expression"
        );

        let final_integrand = outer
            .final_integrand(&graph)?
            .iter()
            .find(|(index, _)| *index == CutCFFIndex::new_all_none())
            .ok_or_else(|| eyre!("the outer banana has no uncut final integrand"))?
            .1;
        assert!(
            final_integrand.contains_symbol(global_marker),
            "final assembly must restore the global numerator left outside every local Taylor operation"
        );
        Ok(())
    }

    #[test]
    fn complete_self_energy_taylor_sum_matches_direct_3d_for_raised_lu_jets() -> Result<()> {
        test_initialise()?;
        let base_graph: Graph = dot!(
            digraph complete_self_energy_taylor_sum {
                num = 1
                edge [particle="scalar_1" num=1]
                node [num=1]

                ext0 [style=invis is_cut=0]
                v2 -> ext0 [id=0]
                ext0 -> v3
                v0 -> v1 [id=1 lmb_id=0 num="Q(1,spenso::mink(4,1))"]
                v0 -> v1 [id=2]
                v3 -> v0 [id=3 lmb_id=1 num="Q(3,spenso::mink(4,1))"]
                v1 -> v2 [id=4]
                v2 -> v3 [id=5]
            },
            "scalars"
        )?;
        assert_eq!(
            base_graph.loop_momentum_basis.edge_signatures[EdgeIndex(3)],
            base_graph.loop_momentum_basis.edge_signatures[EdgeIndex(4)],
            "the GL0-like cograph must expose a genuinely repeated outer channel"
        );

        let outer_factors = [
            "unit outer numerator",
            "factorized cubic outer numerator",
            "factorized affine cubic outer numerator",
        ];
        let mut failures = Vec::new();
        for outer_label in outer_factors {
            let mut graph = base_graph.clone();
            graph.underlying = graph.underlying.map_data_ref(
                |_, _, vertex| {
                    let mut vertex = vertex.clone();
                    if outer_label == "factorized affine cubic outer numerator"
                        && vertex.name.value == "v2"
                    {
                        vertex.num.value *= GS.emr_mom(EdgeIndex(5), GS.cind(0))
                            + GS.emr_mom(EdgeIndex(0), GS.cind(0));
                    }
                    vertex
                },
                |_, edge, _, data| {
                    data.map(|data| {
                        let mut data = data.clone();
                        let factor = match (outer_label, edge) {
                            ("factorized cubic outer numerator", EdgeIndex(3..=5)) => {
                                GS.emr_mom(edge, GS.cind(0))
                            }
                            ("factorized affine cubic outer numerator", EdgeIndex(3 | 5)) => {
                                GS.emr_mom(edge, GS.cind(0))
                            }
                            _ => Atom::one(),
                        };
                        data.num.value *= factor;
                        data
                    })
                },
                |_, data| data.clone(),
            );
            let options = graph.denominator_only_cff_3d_expression_options();
            let numerator = graph.production_numerator_atom_for_full_3d_expression();
            let canonization = graph.get_esurface_canonization(&graph.loop_momentum_basis);
            let production = graph.generate_3d_expression_for_integrand(
                &[],
                &canonization,
                &options,
                Some(&numerator),
            )?;
            let lu_cut = graph
                .determine_raised_esurfaces_from_expression(&production.expression)
                .raised_groups
                .into_iter()
                .find(|group| {
                    group.max_occurence > 1
                        && group.esurface_ids.iter().any(|surface_id| {
                            !production.expression.surfaces.esurface_cache[*surface_id]
                                .external_shift
                                .is_empty()
                        })
                })
                .expect("the repeated outer channel must supply a physical raised LU surface");
            assert_eq!(
                lu_cut.max_occurence, 2,
                "the GL0-like cograph must expose LU orders one and two"
            );
            let mut cutset = CutSet::empty(graph.n_hedges());
            cutset.residue_selector.lu = Some(LuCutSelection {
                raised_group: lu_cut.clone(),
                cut_edge_alternatives: lu_cut
                    .esurface_ids
                    .iter()
                    .map(|surface_id| {
                        production.expression.surfaces.esurface_cache[*surface_id]
                            .energies
                            .clone()
                    })
                    .collect(),
            });
            let uv_filter = graph
                .get_edge_subgraph(EdgeIndex(1))
                .union(&graph.get_edge_subgraph(EdgeIndex(2)));
            let uv_subgraph = InternalSubGraph::cleaned_filter_optimist(uv_filter, graph.as_ref());
            let child_spinney = Spinney::with_scheme(
                uv_subgraph,
                &graph,
                &graph.loop_momentum_basis,
                ApproximationType::MUV,
                1,
            )
            .expect("the self-energy bubble has a compatible sub-LMB");
            let orientation_pattern = OrientationPattern::default();

            let build_child = |mut route_graph: Graph,
                               route_cutset: &CutSet,
                               local_uv_cts_from_expanded_4d_integrands|
             -> Result<BTreeMap<CutCFFIndex, Atom>> {
                let settings = UVgenerationSettings {
                    generate_integrated: false,
                    local_uv_cts_from_expanded_4d_integrands,
                    ..Default::default()
                };
                let localizer = Localizer::new(
                    route_cutset,
                    OrientationProjection::exact_expression(
                        &production,
                        &options,
                        &orientation_pattern,
                        true,
                    ),
                );
                let mut root = Approximation::new(Spinney::empty(&route_graph));
                root.root(&mut route_graph, localizer, &settings)?;
                let mut child = Approximation::new(child_spinney.clone());
                child.simple_approx = Some(
                    root.simple_approx
                        .as_ref()
                        .expect("the root approximation is initialized")
                        .dependent(child.spinney.subgraph.clone()),
                );
                let vakint_settings = vakint::VakintSettings::default();
                child.compute_4d(
                    &route_graph,
                    (crate::utils::vakint()?, &vakint_settings),
                    &root,
                    &settings,
                )?;
                let local_terms = child.local(&route_graph)?.terms()?;
                assert_eq!(
                    local_terms.len(),
                    1,
                    "the DOD-one self-energy Taylor sum must retain one factorized outer term"
                );
                assert_eq!(
                    local_terms[0].denominators.len(),
                    3,
                    "the collected coefficient must retain one simple and one raised UV denominator"
                );
                assert_eq!(
                    [EdgeIndex(1), EdgeIndex(2)].map(|edge| local_terms[0]
                        .denominators
                        .iter()
                        .filter(|denominator| denominator.source_edge == edge)
                        .count()),
                    [1, 2]
                );
                let reduced = child.reduced_subgraph(&root);
                let t_arg = route_graph
                    .numerator(&reduced, root.subgraph())
                    .to_d_dim(GS.dim)
                    .get_single_atom()
                    .expect("the scalar child numerator is available")
                    / route_graph.denominator(&reduced, |_| 1);
                let explicit_taylor_sum = route_graph
                    .uv_rescaled(
                        &reduced,
                        route_graph.n_loops(child.subgraph()),
                        child.lmb(),
                        &t_arg,
                    )
                    .series(GS.rescale, Atom::Zero, 0)
                    .expect("the scalar child Taylor series exists")
                    .to_atom()
                    .replace(GS.rescale)
                    .with(Atom::one())
                    .simplify_metrics()
                    .to_dots()
                    .normalize_dots();
                let local_atom = child.local(&route_graph)?.atom();
                let explicit_taylor_sum = explicit_taylor_sum.expand();
                let explicit_taylor_leaves = match explicit_taylor_sum.as_view() {
                    AtomView::Add(add) => add.iter().count(),
                    _ => 1,
                };
                assert_eq!(
                    explicit_taylor_leaves, 2,
                    "the test-only Taylor expansion must retain its simple and raised contributions"
                );
                assert!(
                    (local_atom + explicit_taylor_sum).together().is_zero(),
                    "the collected local counterterm must equal the negative independently reconstructed Taylor sum"
                );
                child.compute_3d(&root, &mut route_graph, localizer, &settings)?;
                Ok(child.final_integrand(&route_graph)?.iter().collect())
            };
            let direct = build_child(graph.clone(), &cutset, false)?;
            let expanded = build_child(graph.clone(), &cutset, true)?;
            let expected_orders = BTreeSet::from([1, 2]);
            assert_eq!(
                direct
                    .keys()
                    .filter_map(|index| index.lu_cut_order)
                    .collect::<BTreeSet<_>>(),
                expected_orders,
                "the direct child must retain both raised LU orders"
            );
            assert_eq!(
                direct.keys().collect::<Vec<_>>(),
                expanded.keys().collect::<Vec<_>>(),
                "the two routes must finalize the same residue sectors"
            );

            let mass_squared = Atom::num(Rational::from((4, 9)));
            let m_uv = Atom::num(Rational::from((7, 5)));
            let base_spatial_momentum = |edge: EdgeIndex| match usize::from(edge) {
                1 => Atom::num(Rational::from((4, 3))),
                2 => Atom::num(Rational::from((-7, 12))),
                3..=5 => Atom::num(Rational::from((3, 4))),
                _ => Atom::Zero,
            };
            let scaled_on_shell_energy = |edge: EdgeIndex, rescale: Atom| {
                (mass_squared.clone() + (base_spatial_momentum(edge) * rescale).pow(2)).sqrt()
            };
            let representative =
                &production.expression.surfaces.esurface_cache[lu_cut.esurface_ids[0]];
            let root_energy_sum = representative
                .energies
                .iter()
                .map(|edge| scaled_on_shell_energy(*edge, Atom::one()))
                .fold(Atom::Zero, |sum, energy| sum + energy);
            let external_shift_coefficient = representative
                .external_shift
                .iter()
                .map(|(_, coefficient)| *coefficient)
                .sum::<i64>();
            assert_ne!(
                external_shift_coefficient, 0,
                "the selected raised LU surface must have a nonzero external shift"
            );
            let external_energy = -root_energy_sum / Atom::num(external_shift_coefficient);
            let external_edges = graph
                .external_momentum_edge_order()
                .into_iter()
                .collect::<BTreeSet<_>>();
            let rescale_expression = |mut expression: Atom| -> Result<Atom> {
                expression = expression
                    .replace(GS.dim)
                    .with(4)
                    .expand()
                    .simplify_metrics()
                    .to_dots()
                    .normalize_dots()
                    .expand_dots()?;
                let rescale = Atom::var(GS.rescale);
                let compact_minkowski = parse_lit!(spenso::mink(4));
                let time_direction =
                    function!(GS.delta_vec, GS.cind(0), compact_minkowski.as_view());
                expression = expression
                    .replace(function!(
                        GS.dot,
                        time_direction.as_view(),
                        time_direction.as_view()
                    ))
                    .with(Atom::one());
                for left in 0..graph.underlying.n_edges() {
                    let left = EdgeIndex(left);
                    let left_spatial = if !external_edges.contains(&left) {
                        base_spatial_momentum(left) * &rescale
                    } else {
                        Atom::Zero
                    };
                    for right in left.0..graph.underlying.n_edges() {
                        let right = EdgeIndex(right);
                        let right_spatial = if !external_edges.contains(&right) {
                            base_spatial_momentum(right) * &rescale
                        } else {
                            Atom::Zero
                        };
                        expression = expression
                            .replace(
                                Minkowski {}
                                    .new_rep(4)
                                    .inner_product(GS.emr_vec(left), GS.emr_vec(right)),
                            )
                            .with(-&left_spatial * right_spatial);
                    }
                    expression = expression
                        .replace(function!(
                            GS.dot,
                            time_direction.as_view(),
                            GS.emr_vec_index(left, compact_minkowski.as_view())
                        ))
                        .with(Atom::Zero);
                }
                expression = expression
                    .replace(Atom::var(GS.m_uv_expansion))
                    .with(m_uv.clone())
                    .replace(Atom::var(GS.m_uv_vacuum))
                    .with(m_uv.clone())
                    .replace(Atom::var(GS.integrated_loop_scale))
                    .with(Atom::one())
                    .replace(function!(GS.tree_denom_wrapper, W_.x_))
                    .with(W_.x_)
                    .replace(GS.den(W_.a_, W_.b_, W_.c_, W_.d_))
                    .with(W_.d_);
                for edge in 0..graph.underlying.n_edges() {
                    let edge = EdgeIndex(edge);
                    expression = expression
                        .replace(graph.underlying[edge].particle.mass_atom())
                        .with(mass_squared.clone().sqrt());
                    if external_edges.contains(&edge) {
                        expression = expression
                            .replace(GS.emr_mom(edge, GS.cind(0)))
                            .with(external_energy.clone());
                        for spatial_index in 1..=3 {
                            expression = expression
                                .replace(GS.emr_mom(edge, GS.cind(spatial_index)))
                                .with(Atom::Zero)
                                .replace(GS.emr_vec_index(edge, GS.cind(spatial_index)))
                                .with(Atom::Zero);
                        }
                    } else {
                        let on_shell_energy = scaled_on_shell_energy(edge, rescale.clone());
                        expression = expression
                            .replace(GS.emr_mom(edge, GS.cind(1)))
                            .with(base_spatial_momentum(edge) * &rescale)
                            .replace(GS.emr_vec_index(edge, GS.cind(1)))
                            .with(base_spatial_momentum(edge) * &rescale)
                            .replace(GS.ose(edge))
                            .with(on_shell_energy.clone())
                            .replace(cut_energy(edge))
                            .with(on_shell_energy);
                        for spatial_index in 2..=3 {
                            expression = expression
                                .replace(GS.emr_mom(edge, GS.cind(spatial_index)))
                                .with(Atom::Zero)
                                .replace(GS.emr_vec_index(edge, GS.cind(spatial_index)))
                                .with(Atom::Zero);
                        }
                    }
                }
                Ok(expression)
            };
            let eta = rescale_expression(representative.to_atom(&[]))?;
            assert!(
                eta.replace(GS.rescale).with(Atom::one()).expand().is_zero(),
                "the comparison point must lie on the selected raised LU surface"
            );
            let value_and_t_derivative = |expression: Atom| -> Result<[Atom; 2]> {
                let series = rescale_expression(expression)?
                    .series(GS.rescale, Atom::one(), 1)
                    .map_err(|error| eyre!("failed to build raised-LU child t jet: {error}"))?;
                Ok([
                    series.coefficient(Rational::from(0)),
                    series.coefficient(Rational::from(1)),
                ])
            };
            let evaluate_arb = |expression: Atom| -> Result<Complex<F<ArbPrec>>> {
                let parameters = [Atom::var(GS.pi)];
                let rational: ExpressionEvaluator<SymComplex<Fraction<IntegerRing>>> = expression
                    .evaluator(&parameters)
                    .build()
                    .map_err(|error| eyre!("failed to build complete-child evaluator: {error}"))?;
                let mut arb: ExpressionEvaluator<Complex<F<ArbPrec>>> =
                    rational.map_coeff(&|coefficient| {
                        Complex::new(F::from(&coefficient.re), F::from(&coefficient.im))
                    });
                Ok(arb.evaluate_single(&[Complex::new(
                    F::<ArbPrec>::from_f64(0.0).pi(),
                    F::<ArbPrec>::from_f64(0.0),
                )]))
            };
            // One eighth of ArbPrec's requested precision permits substantial,
            // route-dependent loss while preserving a precision-scaled oracle.
            let tolerance = F(ArbPrec::default().epsilon()).sqrt().sqrt().sqrt();

            // Exercise the production boundary which the symbolic comparison below bypasses:
            // forest output is tensor-preprocessed, lowered into an EvaluatorStack, and fed the
            // same first-order HyperDual loop-momentum jet as raised-LU pass one.  Construct the
            // independent Arb oracle directly from the generation LMB, so this remains a
            // comparison of two evaluation routes rather than two uses of the compiled evaluator.
            let production_rational =
                |numerator, denominator| Atom::num(Rational::from((numerator, denominator)));
            let base_loop_atoms = [
                [
                    production_rational(31, 100),
                    production_rational(-47, 100),
                    production_rational(83, 100),
                ],
                [
                    production_rational(-19, 100),
                    production_rational(37, 100),
                    production_rational(61, 100),
                ],
            ];
            let external_atoms = [[
                production_rational(12, 5),
                production_rational(11, 100),
                production_rational(-8, 100),
                production_rational(5, 100),
            ]];
            assert_eq!(
                graph.loop_momentum_basis.loop_edges.len(),
                base_loop_atoms.len()
            );
            assert_eq!(
                graph.loop_momentum_basis.ext_edges.len(),
                external_atoms.len()
            );

            let production_rescale = |mut expression: Atom| -> Result<Atom> {
                expression = expression
                    .replace(GS.dim)
                    .with(4)
                    // Scalarize only the independent test oracle. The generated production
                    // numerator above and the EvaluatorStack input remain factorized.
                    .expand()
                    .simplify_metrics()
                    .to_dots()
                    .normalize_dots()
                    .replace(GS.dim_epsilon)
                    .with(0)
                    .replace(GS.m_uv_expansion)
                    .with(production_rational(7, 5))
                    .replace(GS.m_uv_vacuum)
                    .with(production_rational(7, 5))
                    .replace(function!(GS.tree_denom_wrapper, W_.x_))
                    .with(W_.x_)
                    .replace(GS.den(W_.a_, W_.b_, W_.c_, W_.d_))
                    .with(W_.d_)
                    .replace(function!(GS.ose, W_.mass_, W_.prop_))
                    .with(W_.prop_)
                    .expand_dots()?;
                let t = Atom::var(GS.rescale);
                let loop_components = (0..3)
                    .map(|component| {
                        base_loop_atoms
                            .iter()
                            .map(|momentum| &momentum[component] * &t)
                            .collect::<Vec<_>>()
                    })
                    .collect::<Vec<_>>();
                let external_components = (0..3)
                    .map(|component| {
                        external_atoms
                            .iter()
                            .map(|momentum| momentum[component + 1].clone())
                            .collect::<Vec<_>>()
                    })
                    .collect::<Vec<_>>();
                let external_temporal = external_atoms
                    .iter()
                    .map(|momentum| momentum[0].clone())
                    .collect::<Vec<_>>();
                for edge in 0..graph.underlying.n_edges() {
                    let edge = EdgeIndex(edge);
                    let signature = &graph.loop_momentum_basis.edge_signatures[edge];
                    let spatial = (0..3)
                        .map(|component| {
                            signature
                                .try_compute_momentum(
                                    &loop_components[component],
                                    &external_components[component],
                                )
                                .unwrap_or(Atom::Zero)
                        })
                        .collect::<Vec<_>>();
                    let on_shell_energy = spatial
                        .iter()
                        .fold(Atom::one(), |norm, component| norm + component.pow(2))
                        .sqrt();
                    let temporal = if signature.internal.iter().any(|sign| sign.is_sign()) {
                        on_shell_energy.clone()
                    } else {
                        signature
                            .external
                            .try_apply(&external_temporal)
                            .unwrap_or(Atom::Zero)
                    };
                    expression = expression
                        .replace(GS.ose(edge))
                        .with(on_shell_energy.clone())
                        .replace(cut_energy(edge))
                        .with(on_shell_energy)
                        .replace(GS.emr_mom(edge, GS.cind(0)))
                        .with(temporal);
                    for (component, spatial_component) in spatial.iter().enumerate() {
                        expression = expression
                            .replace(GS.emr_mom(edge, GS.cind(component + 1)))
                            .with(spatial_component.clone())
                            .replace(GS.emr_vec_index(edge, GS.cind(component + 1)))
                            .with(spatial_component.clone());
                    }
                }
                Ok(expression)
            };
            let production_jet = |expression: Atom| -> Result<[Complex<F<ArbPrec>>; 2]> {
                let series = production_rescale(expression)?
                    .series(GS.rescale, Atom::one(), 1)
                    .map_err(|error| {
                        eyre!("failed to build production-boundary Arb jet: {error}")
                    })?;
                Ok([
                    evaluate_arb(series.coefficient(Rational::from(0)))?,
                    evaluate_arb(series.coefficient(Rational::from(1)))?,
                ])
            };
            let production_expressions = [&direct, &expanded]
                .into_iter()
                .map(|route| {
                    route
                        .iter()
                        .find(|(index, _)| index.lu_cut_order == Some(2))
                        .expect("the production-boundary route has an order-two LU sector")
                        .1
                        .clone()
                        .replace(GS.dim)
                        .with(4)
                        .replace(GS.dim_epsilon)
                        .with(0)
                        .replace(GS.m_uv_expansion)
                        .with(production_rational(7, 5))
                        .replace(GS.m_uv_vacuum)
                        .with(production_rational(7, 5))
                        .replace(function!(GS.tree_denom_wrapper, W_.x_))
                        .with(W_.x_)
                        .replace(GS.den(W_.a_, W_.b_, W_.c_, W_.d_))
                        .with(W_.d_)
                        .replace(function!(GS.ose, W_.mass_, W_.prop_))
                        .with(W_.prop_)
                })
                .collect::<Vec<_>>();
            let expected_production_jets = production_expressions
                .iter()
                .cloned()
                .map(production_jet)
                .collect::<Result<Vec<_>>>()?;
            let arb_rational = |numerator, denominator| {
                F::<ArbPrec>::from(&Rational::from((numerator, denominator)))
            };
            let loop_moms: LoopMomenta<F<ArbPrec>> = [
                ThreeMomentum::new(
                    arb_rational(31, 100),
                    arb_rational(-47, 100),
                    arb_rational(83, 100),
                ),
                ThreeMomentum::new(
                    arb_rational(-19, 100),
                    arb_rational(37, 100),
                    arb_rational(61, 100),
                ),
            ]
            .into_iter()
            .collect();
            let external_moms: ExternalFourMomenta<F<ArbPrec>> = [[
                arb_rational(12, 5),
                arb_rational(11, 100),
                arb_rational(-8, 100),
                arb_rational(5, 100),
            ]
            .into()]
            .into_iter()
            .collect();
            let dual_t =
                HyperDual::new(simple_n_deriv_shape(1)).variable(0, F::<ArbPrec>::from_f64(1.0));
            let dual_loop_moms = loop_moms.rescale_with_hyper_dual(&dual_t, None);
            let sample = MomentumSample {
                sample: BareMomentumSample {
                    loop_moms,
                    dual_loop_moms: Some(dual_loop_moms),
                    loop_mom_cache_id: 0,
                    loop_mom_base_cache_id: 0,
                    external_moms,
                    external_mom_cache_id: 0,
                    external_mom_base_cache_id: 0,
                    jacobian: arb_rational(1, 1),
                    orientation: None,
                    parameterization_branch: None,
                },
            };
            // The complete factorized child is deliberately deeper than a normal unit-test
            // thread stack. Use the same 32 MiB worker size as UV profiling rather than flattening
            // the numerator before it reaches the production evaluator.
            let actual_production_jets = std::thread::scope(|scope| -> Result<_> {
                let worker = std::thread::Builder::new()
                    .name("raised-lu-production-boundary".to_owned())
                    .stack_size(32 * 1024 * 1024)
                    .spawn_scoped(scope, || -> Result<_> {
                        let orientations = TiVec::<OrientationID, EdgeVec<Orientation>>::new();
                        let orientation_filter = SubSet::full(orientations.len());
                        let mut param_builder = graph.param_builder.clone();
                        param_builder.initialize_duals(2);
                        let (mut production_evaluator, _) =
                            EvaluatorStack::new_explicit_sum_with_timings(
                                &production_expressions,
                                &param_builder,
                                Some(simple_n_deriv_shape(1)),
                                &EvaluatorSettings::default(),
                            )?;
                        let input = <ArbPrec as GenericEvaluatorFloat>::get_parameters(
                            &mut param_builder,
                            (false, false),
                            &graph,
                            &sample,
                            &[],
                            &[],
                            None,
                            None,
                            None,
                        );
                        production_evaluator.evaluate(
                            input,
                            SingleOrAllOrientations::All {
                                all: &orientations,
                                filter: &orientation_filter,
                            },
                            &RuntimeSettings::default(),
                            &mut EvaluationMetaData::new_empty(),
                            false,
                        )
                    })?;
                worker
                    .join()
                    .map_err(|_| eyre!("raised-LU production-boundary worker panicked"))?
            })?;
            for (route, actual, expected) in ["direct 3D", "projected local 4D"]
                .into_iter()
                .zip(actual_production_jets)
                .zip(&expected_production_jets)
                .map(|((route, actual), expected)| (route, actual, expected))
            {
                let DualOrNot::Dual(actual) = actual else {
                    panic!("{route} order-two production evaluator did not return a HyperDual jet")
                };
                let actual = extract_t_derivatives_complex(actual);
                for (component, actual, expected) in ["f2", "f2 prime"]
                    .into_iter()
                    .zip(actual)
                    .zip(expected)
                    .map(|((component, actual), expected)| (component, actual, expected))
                {
                    let distance = (actual.clone() - expected.clone()).norm().re;
                    let scale = actual.clone().norm().re.max(expected.clone().norm().re);
                    let relative_distance = if scale.is_zero() {
                        distance
                    } else {
                        distance / scale
                    };
                    assert!(
                        relative_distance <= tolerance,
                        "{route} production-boundary {component} mismatch: actual={actual:e}, expected={expected:e}, relative delta={relative_distance:e}, tolerance={tolerance:e}"
                    );
                }
            }

            let consistent_eta = production_rescale(representative.to_atom(&[]))?;
            let consistent_eta_series = consistent_eta
                .series(GS.rescale, Atom::one(), 2)
                .map_err(|error| eyre!("failed to build production-boundary eta jet: {error}"))?;
            let consistent_eta_prime =
                evaluate_arb(consistent_eta_series.coefficient(Rational::from(1)))?;
            let consistent_eta_second =
                evaluate_arb(consistent_eta_series.coefficient(Rational::from(2)) * Atom::num(2))?;
            let mut pass_two_evaluator =
                build_derivative_structure(2, -1, &EvaluatorSettings::default());
            for (route, jet) in ["direct 3D", "projected local 4D"]
                .into_iter()
                .zip(&expected_production_jets)
            {
                let actual = <ArbPrec as GenericEvaluatorFloat>::get_evaluator_single(
                    &mut pass_two_evaluator,
                )(&[
                    jet[0].clone(),
                    jet[1].clone(),
                    consistent_eta_prime.clone(),
                    consistent_eta_second.clone(),
                ]);
                let expected = jet[1].clone() / consistent_eta_prime.clone().pow(2)
                    - jet[0].clone() * &consistent_eta_second / consistent_eta_prime.clone().pow(3);
                let distance = (actual.clone() - expected.clone()).norm().re;
                let scale = actual.clone().norm().re.max(expected.clone().norm().re);
                let relative_distance = if scale.is_zero() {
                    distance
                } else {
                    distance / scale
                };
                assert!(
                    relative_distance <= tolerance,
                    "{route} production pass-two mismatch: actual={actual:e}, expected={expected:e}, relative delta={relative_distance:e}, tolerance={tolerance:e}"
                );
            }
            for (index, direct_expression) in &direct {
                let expanded_expression = expanded
                    .get(index)
                    .expect("the expanded route has every direct residue sector");
                let direct_jet = value_and_t_derivative(direct_expression.clone())?;
                let expanded_jet = value_and_t_derivative(expanded_expression.clone())?;
                for (component, direct_expression, expanded_expression) in
                    ["value", "first t derivative"]
                        .into_iter()
                        .zip(direct_jet)
                        .zip(expanded_jet)
                        .map(|((component, direct), expanded)| (component, direct, expanded))
                {
                    let direct_value = evaluate_arb(direct_expression)?;
                    let expanded_value = evaluate_arb(expanded_expression)?;
                    let distance = (direct_value.clone() - expanded_value.clone()).norm().re;
                    let direct_norm = direct_value.clone().norm().re;
                    let expanded_norm = expanded_value.clone().norm().re;
                    let scale = if direct_norm > expanded_norm {
                        direct_norm
                    } else {
                        expanded_norm
                    };
                    let relative_distance = if scale.is_zero() {
                        distance
                    } else {
                        distance / scale
                    };
                    if relative_distance > tolerance {
                        failures.push(format!(
                            "{outer_label}, {index} {component}: direct={direct_value:e}, expanded={expanded_value:e}, relative delta={relative_distance:e}, tolerance={tolerance:e}"
                        ));
                    }
                }
            }
        }
        assert!(
            failures.is_empty(),
            "complete child Taylor-sum projection failures:\n{}",
            failures.join("\n")
        );
        Ok(())
    }
}
