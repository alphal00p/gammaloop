fn main() -> pyo3_stub_gen::Result<()> {
    let info = linnet_py::stub_info()?;
    let module = info
        .modules
        .get("symbolica.community.graph")
        .expect("linnet StubInfo must contain the linnet module");
    std::fs::write(
        concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/python/symbolica/community/graph/__init__.pyi"
        ),
        linnet_py::canonical_stub(module)?,
    )?;
    Ok(())
}
