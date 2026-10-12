use std::{env, path::PathBuf};

#[path = "../linnest/build_assets.rs"]
mod build_assets;

fn main() -> std::io::Result<()> {
    let manifest = PathBuf::from(env::var("CARGO_MANIFEST_DIR").expect("manifest directory"));
    println!("cargo:rerun-if-changed=../linnest/build_assets.rs");
    let mut output = std::fs::File::create(
        PathBuf::from(env::var_os("OUT_DIR").expect("output directory")).join("embedded_assets.rs"),
    )?;
    build_assets::embed_typst_packages(&mut output, &manifest.join(".."), &["LICENSE"])?;
    build_assets::embed_dir(
        &mut output,
        "TYPST_PACKAGE_ASSETS",
        &manifest.join("vendor/typst-packages"),
        |_| true,
    )
}
