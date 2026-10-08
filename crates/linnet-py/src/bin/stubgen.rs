use linnet_py::PythonModule;

fn main() -> pyo3_stub_gen::Result<()> {
    let root = env!("CARGO_MANIFEST_DIR");
    std::fs::write(
        format!("{root}/linnet_py.pyi"),
        PythonModule::Standalone.stub()?,
    )?;
    std::fs::write(
        format!("{root}/python/symbolica/community/graph/__init__.pyi"),
        PythonModule::Community.stub()?,
    )?;
    Ok(())
}
