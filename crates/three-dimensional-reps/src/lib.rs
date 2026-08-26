pub mod cff_recursion;
mod cut_structure;
#[cfg(feature = "display")]
pub mod display;
pub mod energy_bounds;
#[cfg(feature = "eval")]
pub mod eval;
pub mod expression;
pub mod generation;
pub mod graph_io;
mod graph_signatures;
pub mod surface;
pub mod symbols;
pub mod tree;
pub mod utils;
pub mod validator;

#[cfg(feature = "display")]
pub use display::{DisplayOptions, NumeratorDisplay, render_expression_summary};
pub use energy_bounds::{
    EnergyDirectionReport, EnergyDivergenceReport, auto_numerator_expr_for_bounds,
    energy_divergence_report, normalize_energy_degree_bounds,
};
pub use expression::{
    CFFExpression, GraphOrientation, OrientationData, OrientationExpression, OrientationID,
    OrientationSelector, ResidualDenominator, ThreeDExpression,
};
pub use generation::{
    CffEnergyFactorOwnership, CffGlobalPrefactorSign, Generate3DExpressionOptions,
    GeneratedThreeDExpression, GenerationError, NumeratorSamplingScaleMode, RepresentationMode,
    generate_3d_expression,
};
pub use graph_io::{
    EnergyEdgeIndexMap, GraphInfo, ParsedGraph, RepeatedGroup, ThreeDGraphSource, graph_info,
    repeated_groups,
};
pub use graph_signatures::MomentumSignature;
pub use surface::{
    EsurfaceCollection, EsurfaceID, HsurfaceCollection, HsurfaceID, HybridSurface, HybridSurfaceID,
    LinearEnergyExpr, LinearSurface, LinearSurfaceID, LinearSurfaceKind, SurfaceAtom, SurfaceCache,
};
pub use utils::StringSerializedAtom;
pub use validator::{GraphValidation, validate_parsed_graph};
