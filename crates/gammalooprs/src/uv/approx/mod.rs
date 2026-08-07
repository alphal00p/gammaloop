use crate::{
    cff::expression::{
        GammaLoopOrientationExpression, OrientationExpression, OrientationID, OrientationSelector,
    },
    debug_tags,
    graph::{Graph, LoopMomentumBasis, cuts::CutSet},
    momentum::Sign,
    settings::global::OrientationPattern,
    utils::GS,
    uv::{
        ApproximationType, Spinney, UVgenerationSettings,
        approx::{
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

use std::hash::Hash;

use symbolica::{
    atom::{Atom, AtomCore, AtomOrView},
    function,
};

use linnet::half_edge::involution::{EdgeIndex, EdgeVec, Orientation};
use linnet::half_edge::subgraph::{InternalSubGraph, SuBitGraph, SubSetLike, SubSetOps};
use three_dimensional_reps::Generate3DExpressionOptions;
use typed_index_collections::TiVec;

use super::IntegrandExpr;
use vakint::Vakint;

pub mod final_integrand;
pub mod integrated;
pub mod local_3d;
pub mod local_4d;

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
    Exact(&'a TiVec<OrientationID, OrientationExpression>),
    Coarse(&'a [EdgeVec<Orientation>]),
}

#[derive(Clone, Copy, Debug)]
pub(crate) struct OrientationProjection<'a> {
    source: OrientationProjectionSource<'a>,
    options: Option<&'a Generate3DExpressionOptions>,
    pub(crate) orientation_pattern: &'a OrientationPattern,
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
            orientation_pattern,
        }
    }

    pub(crate) fn exact(
        orientations: &'a TiVec<OrientationID, OrientationExpression>,
        options: &'a Generate3DExpressionOptions,
        orientation_pattern: &'a OrientationPattern,
    ) -> Self {
        Self {
            source: OrientationProjectionSource::Exact(orientations),
            options: Some(options),
            orientation_pattern,
        }
    }

    pub(crate) fn cff_options(self, graph: &Graph) -> Generate3DExpressionOptions {
        self.options
            .cloned()
            .unwrap_or_else(|| graph.denominator_only_cff_3d_expression_options())
    }

    pub(crate) fn orientation_ids(self) -> Vec<OrientationID> {
        match self.source {
            OrientationProjectionSource::Exact(orientations) => orientations
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
            OrientationProjectionSource::Exact(orientations) => Some(orientations),
            OrientationProjectionSource::Coarse(_) => None,
        }
    }

    pub(crate) fn orientation(self, id: OrientationID) -> Option<&'a EdgeVec<Orientation>> {
        match self.source {
            OrientationProjectionSource::Exact(orientations) => orientations
                .get(id)
                .map(|orientation| &orientation.data.orientation),
            OrientationProjectionSource::Coarse(orientations) => orientations.get(id.0),
        }
    }

    /// Apply one production energy map only to the newly owned numerator
    /// fragment. Coarse projectors retain the ordinary UV behavior.
    pub(crate) fn map_numerator(
        self,
        graph: &Graph,
        id: OrientationID,
        numerator: &Atom,
    ) -> Result<Atom> {
        match self.source {
            OrientationProjectionSource::Exact(orientations) => {
                let orientation = orientations.get(id).ok_or_else(|| {
                    eyre!("missing production energy map for orientation {}", id.0)
                })?;
                Ok(numerator.replace_multiple(orientation.energy_replacements_gs(graph)))
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
            let local_3d = if settings.local_uv_cts_from_expanded_4d_integrands {
                let cograph = graph
                    .full_filter()
                    .subtract(self.subgraph())
                    .subtract(&graph.initial_state_cut);
                let source = Full4dCts::with_cograph(self.local(graph)?, graph, &cograph);
                localizer.project_4d(&source, graph, false)?
            } else {
                Local3DCts::root(graph, localizer)?
            };
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
            let cograph = graph
                .full_filter()
                .subtract(self.subgraph())
                .subtract(&graph.initial_state_cut);
            let source = Full4dCts::with_cograph(self.local(graph)?, graph, &cograph);
            localizer.project_4d(&source, graph, true)?
        } else {
            let parent_local = dependent.local_3d(graph)?;
            let parent_integrated = dependent.integrated(graph)?;
            Local3DApproximation::new(localizer, graph, settings).run(
                parent_local,
                parent_integrated,
                self,
                dependent,
                self,
                dependent,
            )?
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
