use std::path::Path;

fn main() {
    embed_runtime_data();
    println!("cargo:rerun-if-env-changed=EXTRA_MACOS_LIBS_FOR_GNU_GCC");

    if cfg!(target_os = "macos") && std::env::var_os("EXTRA_MACOS_LIBS_FOR_GNU_GCC").is_some() {
        let lib_dir = "/opt/local/lib/libgcc";
        let lib_name = "gcc_s.1.1";
        let lib_path = format!("{}/lib{}.dylib", lib_dir, lib_name);
        if Path::new(&lib_path).exists() {
            println!("cargo:rustc-link-search=native={}", lib_dir);
            println!("cargo:rustc-link-lib={}", lib_name);
            println!("cargo:rustc-link-arg=-l{}", lib_name);
        }
    }
}

// Keep ordinary Cargo builds self-contained. Only the explicitly vendored
// runtime files are embedded; downloading and preparing caches are user actions.
fn embed_runtime_data() {
    fn visit(root: &Path, directory: &Path, files: &mut Vec<std::path::PathBuf>) {
        let mut entries = std::fs::read_dir(directory)
            .expect("vendored runtime data directory exists")
            .map(|entry| entry.expect("read vendored data entry").path())
            .collect::<Vec<_>>();
        entries.sort();
        for path in entries {
            if path.is_dir() {
                visit(root, &path, files);
            } else {
                files.push(path.strip_prefix(root).unwrap().to_owned());
            }
        }
    }
    let root =
        std::path::PathBuf::from(std::env::var_os("CARGO_MANIFEST_DIR").unwrap()).join("data");
    println!("cargo:rerun-if-changed={}", root.display());
    let mut files = Vec::new();
    visit(&root, &root, &mut files);
    let mut source = String::from("pub(super) static FILES: &[(&str, &[u8])] = &[\n");
    for path in files {
        source.push_str(&format!(
            "({:?}, include_bytes!({:?})),\n",
            path.to_str().unwrap().replace('\\', "/"),
            root.join(path),
        ));
    }
    source.push_str("];\n");
    let output = std::path::PathBuf::from(std::env::var_os("OUT_DIR").unwrap());
    std::fs::write(output.join("vakint_vendored_data.rs"), source).unwrap();
}
