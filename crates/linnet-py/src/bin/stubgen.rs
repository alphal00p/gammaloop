fn main() -> pyo3_stub_gen::Result<()> {
    std::fs::write(
        concat!(env!("CARGO_MANIFEST_DIR"), "/linnet_py.pyi"),
        linnet_py::canonical_stub()?,
    )?;
    Ok(())
}
