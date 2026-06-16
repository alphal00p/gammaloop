use std::{env, path::PathBuf};
use walkdir::WalkDir;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    emit_macos_gcc_runtime_link();

    #[cfg(feature = "python_api")]
    pyo3_build_config::add_extension_module_link_args();

    // if cfg!(feature = "fjcore") {
    //     println!("cargo:rustc-link-search=./python/gammaloop/dependencies/fjcore");
    //     println!("cargo:rustc-link-lib=stdc++");
    // }
    // models at: $CARGO_MANIFEST_DIR/../../assets/models
    let manifest = PathBuf::from(env::var("CARGO_MANIFEST_DIR").unwrap());
    let models_dir = manifest.join("../../assets/models");

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

fn emit_macos_gcc_runtime_link() {
    println!("cargo:rerun-if-env-changed=EXTRA_MACOS_LIBS_FOR_GNU_GCC");

    if std::env::var_os("EXTRA_MACOS_LIBS_FOR_GNU_GCC").is_none() {
        return;
    }
    if std::env::var("CARGO_CFG_TARGET_OS").as_deref() != Ok("macos") {
        return;
    }

    let lib_dir = "/opt/local/lib/libgcc";
    println!("cargo:rustc-link-search=native={lib_dir}");
    println!("cargo:rustc-link-lib=dylib=gcc_s.1.1");
    println!("cargo:rustc-link-arg=-lgcc_s.1.1");
    println!("cargo:rustc-link-arg=-Wl,-rpath,{lib_dir}");
}
