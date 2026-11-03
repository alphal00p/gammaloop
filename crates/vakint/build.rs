use std::path::{Path, PathBuf};

fn main() {
    embed_runtime_data();
    embed_forcer_sources();
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

fn vendored_files(root: &Path) -> Vec<PathBuf> {
    fn visit(root: &Path, directory: &Path, files: &mut Vec<PathBuf>) {
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
    let mut files = Vec::new();
    visit(root, root, &mut files);
    files
}

fn manifest_directory() -> PathBuf {
    PathBuf::from(std::env::var_os("CARGO_MANIFEST_DIR").unwrap())
}

fn out_directory() -> PathBuf {
    PathBuf::from(std::env::var_os("OUT_DIR").unwrap())
}

// Keep ordinary Cargo builds self-contained. Only the explicitly vendored
// runtime files are embedded; downloading and preparing caches are user actions.
fn embed_runtime_data() {
    let root = manifest_directory().join("data");
    println!("cargo:rerun-if-changed={}", root.display());
    let mut source = String::from("pub(super) static FILES: &[(&str, &[u8])] = &[\n");
    for path in vendored_files(&root) {
        source.push_str(&format!(
            "({:?}, include_bytes!({:?})),\n",
            path.to_str().unwrap().replace('\\', "/"),
            root.join(path),
        ));
    }
    source.push_str("];\n");
    std::fs::write(out_directory().join("vakint_vendored_data.rs"), source).unwrap();
}

// The vendored Forcer sources (about 10 MB of FORM code) are embedded
// gzip-compressed and written out once per process when the Forcer backend
// first runs (`src/p_integral/forcer.rs`).
fn embed_forcer_sources() {
    use std::io::Write;
    let root = manifest_directory().join("form_src").join("forcer");
    println!("cargo:rerun-if-changed={}", root.display());
    let compressed = out_directory().join("forcer");
    std::fs::create_dir_all(&compressed).unwrap();
    let mut source = String::from("pub(super) static FORCER_FILES: &[(&str, &[u8])] = &[\n");
    for path in vendored_files(&root) {
        let name = path.to_str().unwrap().replace('\\', "/");
        let target = compressed.join(format!("{}.gz", name.replace('/', "__")));
        let mut encoder = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
        encoder
            .write_all(&std::fs::read(root.join(&path)).unwrap())
            .unwrap();
        std::fs::write(&target, encoder.finish().unwrap()).unwrap();
        source.push_str(&format!("({name:?}, include_bytes!({target:?})),\n"));
    }
    source.push_str("];\n");
    std::fs::write(out_directory().join("vakint_forcer_sources.rs"), source).unwrap();
}
