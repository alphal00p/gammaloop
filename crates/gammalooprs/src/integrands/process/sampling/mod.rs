//! Graph-aware sampling maps, their canonical channel catalogue, and validation.
//!
//! Map kernels, conditional geometry, symbolic evaluators and multichannel
//! partitions share this boundary. Process-specific bindings remain with the
//! amplitude and cross-section implementations.

pub mod context;
pub mod evaluator;
pub mod joint;
pub mod maps;
pub mod partition;
pub mod reference;
pub mod selection;

pub use context::{PreparedSurfaceStatus, SamplingCutSide};
pub use evaluator::{SamplingDualValue, SamplingExpressionEvaluator};
pub use joint::{SharedEnergyJointGeometry, SharedEnergyJointMap};
pub use maps::{
    ImplicitSurfaceContextPreparer, ImplicitSurfaceRadialContextEvaluator,
    ImplicitSurfaceRadialEvaluator, ImplicitSurfaceRadialMap, SamplingJacobian,
    SamplingMapAcceptanceReport, SamplingMapAffine, SamplingMapComponent, SamplingMapComposition,
    SamplingMapContextTransform, SamplingMapContract, SamplingMapDefinition, SamplingMapEmbedding,
    SamplingMapEvaluation, SamplingMapKernel, SamplingMapPoint, SamplingSupport, SurfaceRadialMap,
    SurfaceRadialPoint,
};
pub use partition::{
    SamplingChannelScore, SamplingPartition, SamplingPartitionMode, SamplingScoreFunction,
};
pub use reference::{
    GaussianReferenceFunction, ReferenceMoments, ReferenceSampleEvaluation, ReferenceSamplingReport,
};
pub use selection::{
    CompiledSamplingChannel, CompiledSamplingMap, ResolvedNamedSamplingChannel,
    ResolvedSamplingBlock, ResolvedSamplingChannelSelection, SamplingCatalogueEntry,
    SamplingChannelBridge, SamplingChannelBridgeAcceptanceReport, SamplingChannelBridgeError,
    SamplingChannelBridgeEvaluation, SamplingChannelCatalogue, SamplingChannelCompileContext,
    SamplingChannelCompileError, SamplingChannelId, SamplingChannelInspection,
    SamplingChannelPreset, SamplingChannelRuntimeContexts, SamplingChannelSelector,
    SamplingCoverageReport, SamplingGeometryKey, SamplingMomentumSampleContext,
    SamplingSelectionError, build_sampling_channel_catalogue,
    build_sampling_channel_catalogue_with_surfaces,
    build_sampling_channel_catalogue_with_surfaces_and_coverage, explicitly_selected_graphs,
    graph_channel_definitions, resolve_sampling_channel_selection,
    resolve_sampling_channel_selection_replacing_default,
};
