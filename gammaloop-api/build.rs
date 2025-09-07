use std::{env, path::PathBuf};
use walkdir::WalkDir;

fn main() {
    if cfg!(feature = "python_api") {
        pyo3_build_config::add_extension_module_link_args();
    }
    // if cfg!(feature = "fjcore") {
    //     println!("cargo:rustc-link-search=./python/gammaloop/dependencies/fjcore");
    //     println!("cargo:rustc-link-lib=stdc++");
    // }
    #[cfg(target_os = "macos")]
    if std::env::var_os("EXTRA_MACOS_LIBS_FOR_GNU_GCC").is_some() {
        println!("cargo:rustc-link-lib=gcc_s");
        println!("cargo:rustc-link-lib=gcc");
        println!("cargo:rustc-link-lib=gfortran");
    }

    // models at: $CARGO_MANIFEST_DIR/../models
    let manifest = PathBuf::from(env::var("CARGO_MANIFEST_DIR").unwrap());
    let models_dir = manifest.join("../models");

    // Re-run if directory metadata changes
    println!("cargo:rerun-if-changed={}", models_dir.display());

    // Re-run if any file inside changes (recursive)
    for entry in WalkDir::new(&models_dir).into_iter().filter_map(Result::ok) {
        if entry.file_type().is_file() {
            println!("cargo:rerun-if-changed={}", entry.path().display());
        }
    }
}
