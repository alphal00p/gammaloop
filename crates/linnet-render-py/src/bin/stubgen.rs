use std::env;
use std::fs;
use std::path::PathBuf;

fn main() -> pyo3_stub_gen::Result<()> {
    let path = env::var_os("GAMMALOOP_WORKSPACE_ROOT")
        .or_else(|| env::var_os("NEXTEST_WORKSPACE_ROOT"))
        .map(|root| PathBuf::from(root).join("crates/linnet-render-py"))
        .unwrap_or_else(|| PathBuf::from(env!("CARGO_MANIFEST_DIR")))
        .join("python/symbolica/community/graph/__init__.pyi");
    fs::write(path, linnet_render_py::canonical_stub()?)?;
    Ok(())
}
