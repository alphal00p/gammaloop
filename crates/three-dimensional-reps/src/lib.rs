pub mod cff_recursion;
pub mod cut_structure;
pub mod diagnostics;
#[cfg(feature = "display")]
pub mod display;
pub mod energy_bounds;
#[cfg(feature = "eval")]
pub mod eval;
pub mod expression;
pub mod generation;
pub mod graph_io;
pub mod graph_signatures;
pub mod surface;
pub mod symbols;
pub mod tree;
pub mod utils;
pub mod validator;

pub use cut_structure::{ContourClosure, Residue, ltd_residues};
pub use diagnostics::{ProfileWarningRow, profile_warnings};
#[cfg(feature = "display")]
pub use display::{DisplayOptions, NumeratorDisplay, render_expression_summary};
pub use energy_bounds::{
    EnergyDirectionReport, EnergyDivergenceReport, auto_numerator_expr_for_bounds,
    energy_divergence_report, normalize_energy_degree_bounds,
};
pub use expression::{
    CFFExpression, GraphOrientation, OrientationData, OrientationExpression, OrientationID,
    OrientationSelector, ThreeDExpression,
};
pub use generation::{
    Generate3DExpressionOptions, GenerationDiagnostic, GenerationWarning,
    NumeratorSamplingScaleMode, RepresentationMode, generate_3d_expression,
    generate_3d_expression_from_parsed,
};
pub use graph_io::{
    EnergyEdgeIndexMap, GraphInfo, ParsedGraph, RepeatedGroup, ThreeDGraphSource, graph_info,
    repeated_groups,
};
pub use graph_signatures::{
    ExtractedSignatureExpression, MomentumSignature, ReconstructDotFormat, ReconstructDotOptions,
    extract_signatures_and_masses_from_symbolica_expression, reconstruct_dot,
    reconstruct_dot_from_expression,
};
pub use surface::{
    EsurfaceCollection, EsurfaceID, HsurfaceCollection, HsurfaceID, HybridSurface, HybridSurfaceID,
    LinearEnergyExpr, LinearSurface, LinearSurfaceID, LinearSurfaceKind, SurfaceAtom, SurfaceCache,
};
pub use utils::StringSerializedAtom;
pub use validator::{GraphValidation, validate_parsed_graph};
