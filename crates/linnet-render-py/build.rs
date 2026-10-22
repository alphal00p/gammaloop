use std::{collections::BTreeMap, env, path::PathBuf};

fn main() {
    symbolica_data::Asset::generate(
        &PathBuf::from(env::var_os("CARGO_MANIFEST_DIR").unwrap()),
        &PathBuf::from(env::var_os("OUT_DIR").unwrap()).join("data_assets.rs"),
        &BTreeMap::new(),
    )
    .unwrap();
}
