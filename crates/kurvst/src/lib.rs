mod curve_api;
pub mod marks;

pub use curve_api::{
    CubicBezierSpec, CurvePathOutput, CurvePoint, FittedCoilInput, HobbySplineSpec,
    ParallelPathSpec, PathFrame, PathFramesSpec, PathLengthSpec, PathTrimmer, PatternInput,
    PatternPathOutput, PatternPathSpec, RegionPart, RegionSamplesSpec, TrimPathSpec,
    curve_fitted_coil_points_bytes, curve_hobby_spline_bytes, curve_hobby_through_bytes,
    curve_parallel_path_bytes, curve_path_frames_bytes, curve_path_intersections_bytes,
    curve_path_length_bytes, curve_pattern_cetz_bytes, curve_pattern_path_bytes,
    curve_region_samples_bytes, curve_stroke_outline_bytes, curve_trim_path_bytes,
    curve_trim_paths_bytes,
};
pub use marks::{mark_geometry_bytes, mark_geometry_packed_bytes};

#[cfg(all(target_arch = "wasm32", feature = "typst-plugin"))]
use wasm_minimal_protocol::*;

#[cfg(all(target_arch = "wasm32", feature = "typst-plugin"))]
initiate_protocol!();

#[cfg(all(target_arch = "wasm32", feature = "typst-plugin"))]
#[wasm_func]
pub fn curve_path_frames(arg: &[u8]) -> Result<Vec<u8>, String> {
    curve_path_frames_bytes(arg)
}

#[cfg(all(target_arch = "wasm32", feature = "typst-plugin"))]
#[wasm_func]
pub fn curve_trim_path(arg: &[u8]) -> Result<Vec<u8>, String> {
    curve_trim_path_bytes(arg)
}

#[cfg(all(target_arch = "wasm32", feature = "typst-plugin"))]
#[wasm_func]
pub fn curve_trim_paths(arg: &[u8]) -> Result<Vec<u8>, String> {
    curve_trim_paths_bytes(arg)
}

#[cfg(all(target_arch = "wasm32", feature = "typst-plugin"))]
#[wasm_func]
pub fn curve_hobby_through(arg: &[u8]) -> Result<Vec<u8>, String> {
    curve_hobby_through_bytes(arg)
}

#[cfg(all(target_arch = "wasm32", feature = "typst-plugin"))]
#[wasm_func]
pub fn curve_hobby_spline(arg: &[u8]) -> Result<Vec<u8>, String> {
    curve_hobby_spline_bytes(arg)
}

#[cfg(all(target_arch = "wasm32", feature = "typst-plugin"))]
#[wasm_func]
pub fn curve_fitted_coil_points(arg: &[u8]) -> Result<Vec<u8>, String> {
    curve_fitted_coil_points_bytes(arg)
}

#[cfg(all(target_arch = "wasm32", feature = "typst-plugin"))]
#[wasm_func]
pub fn curve_pattern_path(arg: &[u8]) -> Result<Vec<u8>, String> {
    curve_pattern_path_bytes(arg)
}

#[cfg(all(target_arch = "wasm32", feature = "typst-plugin"))]
#[wasm_func]
pub fn curve_pattern_cetz(arg: &[u8]) -> Result<Vec<u8>, String> {
    curve_pattern_cetz_bytes(arg)
}

#[cfg(all(target_arch = "wasm32", feature = "typst-plugin"))]
#[wasm_func]
pub fn curve_parallel_path(arg: &[u8]) -> Result<Vec<u8>, String> {
    curve_parallel_path_bytes(arg)
}

#[cfg(all(target_arch = "wasm32", feature = "typst-plugin"))]
#[wasm_func]
pub fn curve_stroke_outline(arg: &[u8]) -> Result<Vec<u8>, String> {
    curve_stroke_outline_bytes(arg)
}

#[cfg(all(target_arch = "wasm32", feature = "typst-plugin"))]
#[wasm_func]
pub fn curve_path_length(arg: &[u8]) -> Result<Vec<u8>, String> {
    curve_path_length_bytes(arg)
}

#[cfg(all(target_arch = "wasm32", feature = "typst-plugin"))]
#[wasm_func]
pub fn curve_path_intersections(arg: &[u8]) -> Result<Vec<u8>, String> {
    curve_path_intersections_bytes(arg)
}

#[cfg(all(target_arch = "wasm32", feature = "typst-plugin"))]
#[wasm_func]
pub fn curve_region_samples(arg: &[u8]) -> Result<Vec<u8>, String> {
    curve_region_samples_bytes(arg)
}

#[cfg(all(target_arch = "wasm32", feature = "typst-plugin"))]
#[wasm_func]
pub fn mark_geometry(arg: &[u8]) -> Result<Vec<u8>, String> {
    mark_geometry_bytes(arg)
}

#[cfg(all(target_arch = "wasm32", feature = "typst-plugin"))]
#[wasm_func]
pub fn mark_geometry_packed(arg: &[u8]) -> Result<Vec<u8>, String> {
    mark_geometry_packed_bytes(arg)
}
