fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("cargo:rerun-if-env-changed=EXTRA_MACOS_LIBS_FOR_GNU_GCC");

    // Standalone programs must use the same locally developed upstream crates
    // as their generator. Cargo does not propagate workspace patches into the
    // independent manifest embedded in an exported rust-script.
    let workspace = std::path::PathBuf::from(std::env::var("CARGO_MANIFEST_DIR")?)
        .join("../..")
        .canonicalize()?;
    let manifest_path = workspace.join("Cargo.toml");
    println!("cargo:rerun-if-changed={}", manifest_path.display());
    let manifest = match std::fs::read_to_string(&manifest_path) {
        Ok(contents) => contents.parse::<toml::Table>()?,
        Err(error) if error.kind() == std::io::ErrorKind::NotFound => toml::Table::new(),
        Err(error) => return Err(error.into()),
    };
    let mut patches = toml::Table::new();
    if let Some(overrides) = manifest
        .get("patch")
        .and_then(|patch| patch.get("crates-io"))
        .and_then(toml::Value::as_table)
    {
        for name in ["symbolica", "numerica", "graphica", "symjit"] {
            if let Some(mut source) = overrides.get(name).cloned() {
                if let Some(path) = source.get("path").and_then(toml::Value::as_str) {
                    source["path"] = workspace
                        .join(path)
                        .canonicalize()?
                        .to_string_lossy()
                        .into_owned()
                        .into();
                }
                patches.insert(name.into(), source);
            }
        }
    }
    let dependencies = if patches.is_empty() {
        String::new()
    } else {
        let patch = toml::Table::from_iter([("crates-io".into(), patches.into())]);
        let manifest = toml::Table::from_iter([("patch".into(), patch.into())]);
        toml::to_string(&manifest)?
    };
    std::fs::write(
        std::path::PathBuf::from(std::env::var("OUT_DIR")?).join("standalone_dependencies.toml"),
        dependencies,
    )?;

    #[cfg(target_os = "macos")]
    if std::env::var_os("EXTRA_MACOS_LIBS_FOR_GNU_GCC").is_some() {
        println!("cargo:rustc-link-lib=gcc_s");
        println!("cargo:rustc-link-lib=gcc");
        println!("cargo:rustc-link-lib=gfortran");
    }
    Ok(())
}
