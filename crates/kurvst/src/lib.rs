mod curve_api;

pub use curve_api::{
    curve_hobby_through_bytes, curve_parallel_cubic_bytes, curve_pattern_cubic_bytes,
    curve_split_cubic_bytes, curve_split_quadratic_through_bytes, curve_trim_cubic_bytes,
};

#[cfg(target_arch = "wasm32")]
use wasm_minimal_protocol::*;

#[cfg(target_arch = "wasm32")]
initiate_protocol!();

#[cfg(target_arch = "wasm32")]
#[wasm_func]
pub fn curve_split_cubic(arg: &[u8]) -> Result<Vec<u8>, String> {
    curve_split_cubic_bytes(arg)
}

#[cfg(target_arch = "wasm32")]
#[wasm_func]
pub fn curve_trim_cubic(arg: &[u8]) -> Result<Vec<u8>, String> {
    curve_trim_cubic_bytes(arg)
}

#[cfg(target_arch = "wasm32")]
#[wasm_func]
pub fn curve_split_quadratic_through(arg: &[u8]) -> Result<Vec<u8>, String> {
    curve_split_quadratic_through_bytes(arg)
}

#[cfg(target_arch = "wasm32")]
#[wasm_func]
pub fn curve_hobby_through(arg: &[u8]) -> Result<Vec<u8>, String> {
    curve_hobby_through_bytes(arg)
}

#[cfg(target_arch = "wasm32")]
#[wasm_func]
pub fn curve_pattern_cubic(arg: &[u8]) -> Result<Vec<u8>, String> {
    curve_pattern_cubic_bytes(arg)
}

#[cfg(target_arch = "wasm32")]
#[wasm_func]
pub fn curve_parallel_cubic(arg: &[u8]) -> Result<Vec<u8>, String> {
    curve_parallel_cubic_bytes(arg)
}
