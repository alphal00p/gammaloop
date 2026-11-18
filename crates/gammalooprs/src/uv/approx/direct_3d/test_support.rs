//! Fixtures shared by the kernel tests of the direct local-3D route: thermal
//! CFFs and forest nodes routed like the spinneys of a wood.

use color_eyre::Result;
use eyre::eyre;
use linnet::half_edge::subgraph::InternalSubGraph;
use symbolica::atom::Atom;
use three_dimensional_reps::{
    Generate3DExpressionOptions, GeneratedThreeDExpression, MediumMode, NumeratorSamplingScaleMode,
    RepresentationMode,
};

use crate::cff::esurface::Esurface;
use crate::cff::hsurface::Hsurface;
use crate::graph::cuts::CutSet;
use crate::graph::{FeynmanGraph, Graph};
use crate::settings::global::{GenerationSettings, MediumSettings, OrientationPattern};
use crate::uv::approx::OrientationProjection;
use crate::uv::approx::direct_3d::{Direct3dCts, DirectResidueBranches};
use crate::uv::approx::local_3d::Localizer;
use crate::uv::hedge_poset::OwnedForestNode;
use crate::uv::thermal_zero_mode::ThermalZeroModeRegion;
use crate::uv::{ApproximationType, Spinney};

/// A CFF of a graph in a medium, which owns the residue maps of its
/// branches.
pub(super) struct ThermalCff {
    production: GeneratedThreeDExpression<Esurface, Hsurface>,
    options: Generate3DExpressionOptions,
    pattern: OrientationPattern,
}

impl ThermalCff {
    /// The denominator-only thermal CFF of `graph`.
    pub(super) fn new(graph: &mut Graph) -> Result<Self> {
        Self::in_medium(graph, MediumMode::ThermodynamicEquilibrium)
    }

    /// The denominator-only CFF of `graph` in `medium_mode`, which also
    /// becomes the medium of `graph`.
    pub(super) fn in_medium(graph: &mut Graph, medium_mode: MediumMode) -> Result<Self> {
        graph.set_medium_mode(medium_mode)?;
        let mut options = graph.denominator_only_cff_3d_expression_options();
        options.medium_mode = medium_mode;
        Self::generate(graph, options, Some(&Atom::one()))
    }

    /// The production thermal CFF of `graph`: its numerator's energy-degree
    /// bounds add the contact terms of every numerator it can carry.
    pub(super) fn production(graph: &mut Graph, vacuum_subtraction: bool) -> Result<Self> {
        let settings = GenerationSettings {
            medium: MediumSettings {
                mode: MediumMode::ThermodynamicEquilibrium,
                vacuum_subtraction,
            },
            ..Default::default()
        };
        let options = graph.production_3d_expression_options(&settings, RepresentationMode::Cff)?;
        Self::generate(graph, options, None)
    }

    /// The thermal CFF of `graph` for the energy-degree bounds of
    /// `numerator`, which samples every numerator energy at multiples of the
    /// uniform scale `M`.
    pub(super) fn sampled(graph: &mut Graph, numerator: &Atom) -> Result<Self> {
        let mut options = graph.denominator_only_cff_3d_expression_options();
        options.medium_mode = MediumMode::ThermodynamicEquilibrium;
        options.numerator_sampling_scale = NumeratorSamplingScaleMode::All;
        Self::generate(graph, options, Some(numerator))
    }

    /// Generate the CFF of `graph` with `options`; an `analysis_numerator`
    /// replaces their energy-degree bounds by its own.
    fn generate(
        graph: &mut Graph,
        options: Generate3DExpressionOptions,
        analysis_numerator: Option<&Atom>,
    ) -> Result<Self> {
        let canonization = graph.get_esurface_canonization(&graph.loop_momentum_basis);
        let production = graph.generate_3d_expression_for_integrand(
            &[],
            &canonization,
            &options,
            analysis_numerator,
        )?;
        Ok(Self {
            production,
            options,
            pattern: OrientationPattern::default(),
        })
    }

    pub(super) fn orientation(&self, vacuum_subtraction: bool) -> OrientationProjection<'_> {
        OrientationProjection::exact_expression(
            &self.production,
            &self.options,
            &self.pattern,
            false,
        )
        .with_vacuum_subtraction(vacuum_subtraction)
    }

    /// The production residue branches of `graph` under
    /// [`Self::orientation`].
    pub(super) fn branches(
        &self,
        graph: &Graph,
        vacuum_subtraction: bool,
    ) -> Result<DirectResidueBranches> {
        let cutset = CutSet::empty(graph.n_hedges());
        Direct3dCts::root(
            graph,
            Localizer::new(&cutset, self.orientation(vacuum_subtraction)),
        )?
        .branches()
    }
}

/// A forest node on `edges`, routed compatibly with the graph basis as the
/// spinneys of a wood are.
pub(super) fn forest_node(
    graph: &Graph,
    edges: &[usize],
    scheme: ApproximationType,
    dod: i32,
) -> Result<OwnedForestNode> {
    let subgraph = InternalSubGraph::try_new(graph.edge_subgraph(edges), graph.as_ref())
        .ok_or_else(|| eyre!("edges {edges:?} do not form an internal subgraph"))?;
    let spinney = Spinney::with_scheme(subgraph, graph, &graph.loop_momentum_basis, scheme, dod)?
        .ok_or_else(|| eyre!("edges {edges:?} do not form a divergent spinney"))?;
    Ok(OwnedForestNode {
        spinney,
        topo_order: 0,
    })
}

/// The zero-mode host on `edges` of degree `dod` that hosts `region`.
pub(super) fn host_node(
    graph: &Graph,
    edges: &[usize],
    dod: i32,
    region: ThermalZeroModeRegion,
) -> Result<OwnedForestNode> {
    let OwnedForestNode {
        spinney,
        topo_order,
    } = forest_node(graph, edges, ApproximationType::IR, dod)?;
    Ok(OwnedForestNode {
        spinney: spinney.with_thermal_zero_mode(region),
        topo_order,
    })
}

/// The empty forest node, the root of every forest.
pub(super) fn empty_forest_node(graph: &Graph) -> OwnedForestNode {
    OwnedForestNode {
        spinney: Spinney::empty(graph),
        topo_order: 0,
    }
}
