use pyo3::{
    Bound, PyResult, Python, pyfunction, pymodule,
    types::{PyAnyMethods, PyModule, PyModuleMethods},
    wrap_pyfunction,
};
use symbolica::api::python::{Citation, SymbolicaCommunityModule, create_symbolica_module};

macro_rules! register_module {
    ($m:expr, $module_type:ty) => {{
        let native_name = format!("{}_native", <$module_type>::get_name());

        #[pyfunction]
        fn initialize_module(py: Python) -> PyResult<()> {
            <$module_type>::initialize(py)
        }

        let child_module = PyModule::new($m.py(), &native_name)?;
        child_module.add_function(wrap_pyfunction!(initialize_module, &child_module)?)?;

        <$module_type>::register_module(&child_module)?;
        $m.add_submodule(&child_module)?;

        $m.py().import("sys")?.getattr("modules")?.set_item(
            format!("symbolica.community.{}", native_name),
            &child_module,
        )?;
    }};
}

#[pymodule]
fn core(m: &Bound<'_, PyModule>) -> PyResult<()> {
    create_symbolica_module(m)?;
    m.add_function(pyo3::wrap_pyfunction!(get_citations, m)?)?;
    register_module!(m, feynkit_py::FeynkitModule);
    register_module!(m, spynso3::SpensoModule);
    Ok(())
}

#[pyo3::pyfunction]
fn get_citations() -> Vec<Citation> {
    let mut citations = vec![Citation {
        id: "doi:10.5281/zenodo.17054381".into(),
        reference: "Ben Ruijl. Symbolica (2025). doi:10.5281/zenodo.17054381.".into(),
        bibtex: r#"@software{ruijl_symbolica_2025,
  author = {Ruijl, Ben},
  title = {Symbolica},
  year = {2025},
  version = {0.18.0},
  doi = {10.5281/zenodo.17054381},
  url = {https://zenodo.org/records/17054381}
}"#
        .into(),
        reasons: vec!["Symbolic and numerical computation with Symbolica.".into()],
        description: String::new(),
        relevance: None,
    }];
    citations.extend(feynkit_py::FeynkitModule::get_citations());
    citations.extend(spynso3::SpensoModule::get_citations());
    citations
}
