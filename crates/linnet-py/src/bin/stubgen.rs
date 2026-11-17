use std::env;
use std::fs;
use std::path::PathBuf;

fn main() -> pyo3_stub_gen::Result<()> {
    let root = env::var_os("GAMMALOOP_WORKSPACE_ROOT")
        .or_else(|| env::var_os("NEXTEST_WORKSPACE_ROOT"))
        .map(|root| PathBuf::from(root).join("crates/linnet-py"))
        .unwrap_or_else(|| PathBuf::from(env!("CARGO_MANIFEST_DIR")));
    let info = linnet_py::stub_info()?;
    let module = info
        .modules
        .get("symbolica.community.graph")
        .expect("linnet StubInfo must contain the linnet module");
    fs::write(
        root.join("python/symbolica/community/graph/__init__.pyi"),
        linnet_py::canonical_stub(module)?,
    )?;
    Ok(())
}
