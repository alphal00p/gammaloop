use std::{collections::BTreeMap, env, path::PathBuf};

fn main() {
    let root = PathBuf::from(env::var_os("CARGO_MANIFEST_DIR").unwrap());
    let manifest: serde_json::Value =
        serde_json::from_slice(&std::fs::read(root.join("data.json")).unwrap()).unwrap();
    let external = manifest["assets"]
        .as_array()
        .unwrap()
        .iter()
        .zip(typst_assets::fonts())
        .map(|(entry, bytes)| (entry["name"].as_str().unwrap().to_owned(), bytes))
        .collect::<BTreeMap<_, _>>();
    symbolica_data::Asset::generate(
        &root,
        &PathBuf::from(env::var_os("OUT_DIR").unwrap()).join("data_assets.rs"),
        &external,
    )
    .unwrap();
}
