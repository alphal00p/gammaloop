//! Writes the standalone and Community Linnet stubs; `alphal00p-docs-python-exporter --features
//! linnet` writes and checks the same files. The final extension selects Symbolica's arithmetic
//! backend, so name one when running the generator:
//!
//! ```text
//! cargo run -p linnet-py --bin stubgen \
//!     --features community,python_stubgen,symbolica/integer-gmp,symbolica/float-mpfr
//! ```

use linnet_py::{stub, stub_info, PythonModule};

fn main() -> pyo3_stub_gen::Result<()> {
    let root = env!("CARGO_MANIFEST_DIR");
    let info = stub_info()?;
    let module = info
        .modules
        .get(PythonModule::COMPILED.name())
        .expect("linnet StubInfo must contain the linnet module");
    std::fs::write(
        format!("{root}/linnet_py.pyi"),
        stub(module, PythonModule::Standalone)?,
    )?;
    std::fs::write(
        format!("{root}/python/symbolica/community/graph/__init__.pyi"),
        stub(module, PythonModule::Community)?,
    )?;
    Ok(())
}
