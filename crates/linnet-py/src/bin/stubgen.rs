fn main() -> pyo3_stub_gen::Result<()> {
    std::fs::write(
        concat!(env!("CARGO_MANIFEST_DIR"), "/linnet.pyi"),
        linnet_py::canonical_stub()?,
    )?;
    Ok(())
}
