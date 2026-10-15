fn main() -> pyo3_stub_gen::Result<()> {
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("python/symbolica/community/graph/__init__.pyi");
    std::fs::write(path, linnet_render_py::canonical_stub()?)?;
    Ok(())
}
