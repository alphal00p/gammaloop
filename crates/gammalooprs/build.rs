use std::{env, path::PathBuf};
use walkdir::WalkDir;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("cargo:rerun-if-env-changed=EXTRA_MACOS_LIBS_FOR_GNU_GCC");

    #[cfg(target_os = "macos")]
    if std::env::var_os("EXTRA_MACOS_LIBS_FOR_GNU_GCC").is_some() {
        println!("cargo:rustc-link-lib=gcc_s");
        println!("cargo:rustc-link-lib=gcc");
        println!("cargo:rustc-link-lib=gfortran");
    }

    // models at: $CARGO_MANIFEST_DIR/assets/models
    let manifest = PathBuf::from(env::var("CARGO_MANIFEST_DIR").unwrap());
    let models_dir = manifest.join("assets/models");

    // Re-run if directory metadata changes
    println!("cargo:rerun-if-changed={}", models_dir.display());

    // Re-run if any file inside changes (recursive)
    for entry in WalkDir::new(&models_dir)
        .follow_links(true)
        .into_iter()
        .filter_map(Result::ok)
    {
        if entry.file_type().is_file() {
            println!("cargo:rerun-if-changed={}", entry.path().display());
        }
    }

    Ok(())
}
