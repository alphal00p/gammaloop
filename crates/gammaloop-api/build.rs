use std::{env, path::PathBuf};
use walkdir::WalkDir;

#[path = "../linnest/build_assets.rs"]
mod build_assets;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("cargo:rerun-if-env-changed=EXTRA_MACOS_LIBS_FOR_GNU_GCC");

    #[cfg(feature = "python")]
    pyo3_build_config::add_extension_module_link_args();

    #[cfg(target_os = "macos")]
    if std::env::var_os("EXTRA_MACOS_LIBS_FOR_GNU_GCC").is_some() {
        println!("cargo:rustc-link-lib=gcc_s");
        println!("cargo:rustc-link-lib=gcc");
        println!("cargo:rustc-link-lib=gfortran");
    }

    // models at: $CARGO_MANIFEST_DIR/../../assets/models
    let manifest = PathBuf::from(env::var("CARGO_MANIFEST_DIR").unwrap());
    let models_dir = manifest.join("../../assets/models");

    println!("cargo:rerun-if-changed=../linnest/build_assets.rs");
    let mut embedded_assets = std::fs::File::create(
        PathBuf::from(env::var_os("OUT_DIR").unwrap()).join("embedded_assets.rs"),
    )?;
    build_assets::embed_dir(
        &mut embedded_assets,
        "GAMMALOOP_TEMPLATE_ASSETS",
        &manifest.join("../../assets/embedded/drawing/templates"),
        |path| path.extension().is_some_and(|extension| extension == "typ"),
    )?;
    build_assets::embed_typst_packages(&mut embedded_assets, &manifest.join(".."), &[])?;

    // Re-run if directory metadata changes
    println!("cargo:rerun-if-changed={}", models_dir.display());

    // Re-run if any file inside changes (recursive)
    for entry in WalkDir::new(&models_dir).into_iter().filter_map(Result::ok) {
        if entry.file_type().is_file() {
            println!("cargo:rerun-if-changed={}", entry.path().display());
        }
    }

    Ok(())
}
