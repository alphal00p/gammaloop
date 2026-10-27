//! Regenerate the Community graph package and documentation stubs.

use std::fs;
use std::path::PathBuf;

use linnet_py as graph;
use pyo3_stub_gen::Result;

fn main() -> Result<()> {
    let root = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let info = graph::stub_info()?;
    let module = info
        .modules
        .get("symbolica.community.graph")
        .expect("Linnet StubInfo must contain the graph module");
    let stub = graph::canonical_stub(module)?;
    for path in [
        root.join("python/symbolica/community/graph/__init__.pyi"),
        root.join("../../docs/api/python/linnet-graph.pyi"),
        root.join(
            "../../examples/notebooks/symbolica-host/python/symbolica/community/graph/__init__.pyi",
        ),
    ] {
        fs::write(path, &stub)?;
    }
    Ok(())
}
