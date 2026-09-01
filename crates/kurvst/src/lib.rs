mod curve_api;

pub use curve_api::{
    curve_hobby_spline_bytes, curve_hobby_through_bytes, curve_parallel_path_bytes,
    curve_path_length_bytes, curve_pattern_path_bytes, curve_trim_path_bytes,
};

#[cfg(all(target_arch = "wasm32", feature = "typst-plugin"))]
use wasm_minimal_protocol::*;

#[cfg(all(target_arch = "wasm32", feature = "typst-plugin"))]
initiate_protocol!();

#[cfg(all(target_arch = "wasm32", feature = "typst-plugin"))]
#[wasm_func]
pub fn curve_trim_path(arg: &[u8]) -> Result<Vec<u8>, String> {
    curve_trim_path_bytes(arg)
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
pub fn curve_pattern_path(arg: &[u8]) -> Result<Vec<u8>, String> {
    curve_pattern_path_bytes(arg)
}

#[cfg(all(target_arch = "wasm32", feature = "typst-plugin"))]
#[wasm_func]
pub fn curve_parallel_path(arg: &[u8]) -> Result<Vec<u8>, String> {
    curve_parallel_path_bytes(arg)
}

#[cfg(all(target_arch = "wasm32", feature = "typst-plugin"))]
#[wasm_func]
pub fn curve_path_length(arg: &[u8]) -> Result<Vec<u8>, String> {
    curve_path_length_bytes(arg)
}
