use std::{env, fs::File, path::PathBuf};

#[path = "../linnest/build_assets.rs"]
mod build_assets;

fn main() -> std::io::Result<()> {
    let manifest = PathBuf::from(env::var_os("CARGO_MANIFEST_DIR").expect("manifest directory"));
    println!("cargo:rerun-if-changed=../linnest/build_assets.rs");
    let mut output = File::create(
        PathBuf::from(env::var_os("OUT_DIR").expect("output directory")).join("embedded_assets.rs"),
    )?;
    build_assets::embed_dir(
        &mut output,
        "TEMPLATE_ASSETS",
        &manifest.join("templates"),
        |_| true,
    )?;
    build_assets::embed_typst_packages(
        &mut output,
        &manifest.join(".."),
        &["ec-layout.wasm", "LICENSE.ec-layout", "LICENSE.clarabel"],
    )?;
    Ok(())
}
