use miniz_oxide::inflate::decompress_to_vec_zlib_with_limit;
use wasm_minimal_protocol::*;

initiate_protocol!();

const MAX_DECOMPRESSED_SIZE: usize = 64 * 1024 * 1024;

/// Expand a zlib-wrapped DEFLATE stream into WebAssembly module bytes.
#[wasm_func]
pub fn decompress(input: &[u8]) -> Result<Vec<u8>, String> {
    decompress_to_vec_zlib_with_limit(input, MAX_DECOMPRESSED_SIZE)
        .map_err(|error| format!("could not decompress the plugin: {error}"))
}
