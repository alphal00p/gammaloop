use std::collections::BTreeMap;
use std::env;
use std::path::PathBuf;

use symbolica_data::Asset;

fn main() {
    let manifest =
        PathBuf::from(env::var_os("CARGO_MANIFEST_DIR").expect("Cargo sets CARGO_MANIFEST_DIR"));
    let output = PathBuf::from(env::var_os("OUT_DIR").expect("Cargo sets OUT_DIR"));
    Asset::generate(&manifest, &output.join("data_assets.rs"), &BTreeMap::new())
        .expect("bundled graph font assets must generate successfully");
}
