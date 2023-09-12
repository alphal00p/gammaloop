use crate::{cli_functions::cli, model::Model};
use git_version::git_version;
use symbolica;
const GIT_VERSION: &str = git_version!();

#[allow(unused)]
use pyo3::{
    exceptions,
    prelude::*,
    pyclass,
    pyclass::CompareOp,
    pyfunction, pymethods, pymodule,
    types::{PyModule, PyTuple, PyType},
    wrap_pyfunction, FromPyObject, IntoPy, PyObject, PyRef, PyResult, Python,
};

#[pyfunction]
#[pyo3(name = "rust_cli")]
fn cli_wrapper(py: Python) -> PyResult<()> {
    /*
    // This one includes python and the name of the wrapper script itself, e.g.
    // `["/home/ferris/.venv/bin/python", "/home/ferris/.venv/bin/print_cli_args", "a", "b", "c"]`
    println!("{:?}", env::args().collect::<Vec<_>>());
    // This one includes only the name of the wrapper script itself, e.g.
    // `["/home/ferris/.venv/bin/print_cli_args", "a", "b", "c"])`
    println!(
        "{:?}",
        py.import("sys")?
            .getattr("argv")?
            .extract::<Vec<String>>()?
    );
    Ok(())
    */

    cli(&py
        .import("sys")?
        .getattr("argv")?
        .extract::<Vec<String>>()?)
    .map_err(|e| exceptions::PyException::new_err(e.to_string()))
}

#[pymodule]
#[pyo3(name = "_gammaloop")]
fn gammalooprs(_py: Python, m: &PyModule) -> PyResult<()> {
    // TODO: Verify that indeed Python logger level is used in that case.
    pyo3_log::init();
    m.add_class::<PythonWorker>()?;
    m.add("git_version", GIT_VERSION)?;
    m.add_wrapped(wrap_pyfunction!(cli_wrapper))?;
    Ok(())
}

#[pyclass(name = "Worker")]
pub struct PythonWorker {
    pub model: Model,
    sb_state: symbolica::state::State,
    sb_workspace: symbolica::state::Workspace,
}
// TODO: Improve error broadcasting to Python so as to show rust backtrace

#[pymethods]
impl PythonWorker {
    #[classmethod]
    pub fn new(_cls: &PyType) -> PyResult<PythonWorker> {
        Ok(PythonWorker {
            model: Model::default(),
            sb_state: symbolica::state::State::new(),
            sb_workspace: symbolica::state::Workspace::new(),
        })
    }

    pub fn load_model(&mut self, file_path: &str) -> PyResult<()> {
        Model::from_file(
            String::from(file_path),
            &mut self.sb_state,
            &self.sb_workspace,
        )
        .map_err(|e| exceptions::PyException::new_err(e.to_string()))
        .map(|m| {
            self.model = m;
            ()
        })
    }

    pub fn load_model_from_yaml_str(&mut self, yaml_str: &str) -> PyResult<()> {
        Model::from_yaml_str(
            String::from(yaml_str),
            &mut self.sb_state,
            &self.sb_workspace,
        )
        .map_err(|e| exceptions::PyException::new_err(e.root_cause().to_string()))
        .map(|m| {
            self.model = m;
            ()
        })
    }
}
