use std::{env, path::PathBuf};

use walkdir::WalkDir;

fn main() {
    let manifest = PathBuf::from(env::var("CARGO_MANIFEST_DIR").expect("manifest directory"));
    for directory in [
        manifest.join("../linnest/typst"),
        manifest.join("../kurvst/typst"),
        manifest.join("vendor/typst-packages"),
    ] {
        println!("cargo:rerun-if-changed={}", directory.display());
        for entry in WalkDir::new(directory).into_iter().filter_map(Result::ok) {
            if entry.file_type().is_file() {
                println!("cargo:rerun-if-changed={}", entry.path().display());
            }
        }
    }
}
