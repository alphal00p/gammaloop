use std::fs;

const STUB_PATH: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/python/symbolica/community/hepkit/__init__.pyi"
);

fn main() -> pyo3_stub_gen::Result<()> {
    let info = feynkit_py::stub_info()?;
    let module = info
        .modules
        .get("symbolica.community.hepkit")
        .ok_or_else(|| {
            std::io::Error::new(
                std::io::ErrorKind::NotFound,
                "FeynKit did not contribute a symbolica.community.hepkit stub module",
            )
        })?;
    fs::write(STUB_PATH, feynkit_py::stub_source(module))?;
    Ok(())
}
