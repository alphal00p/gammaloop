use pyo3::{
    Bound, PyResult, pyfunction,
    types::{PyModule, PyModuleMethods},
    wrap_pyfunction,
};
#[cfg(feature = "python_stubgen")]
use pyo3_stub_gen::derive::gen_stub_pyfunction;
use symbolica::{
    api::python::PythonExpression,
    atom::{Atom, AtomCore},
};

use crate::selective_expand::SelectiveExpand;

pub type PythonTerm = (PythonExpression, PythonExpression);

pub(super) fn python_terms(terms: Vec<(Atom, Atom)>) -> Vec<PythonTerm> {
    terms
        .into_iter()
        .map(|(structure, coefficient)| (structure.into(), coefficient.into()))
        .collect()
}

#[cfg_attr(
    feature = "python_stubgen",
    gen_stub_pyfunction(module = "symbolica.community.idenso")
)]
#[pyfunction]
/// Selectively expand around the supplied expression patterns.
///
/// Results retain the Rust API's `(structure, coefficient)` factorization.
pub fn expand_in_patterns(
    expression: &PythonExpression,
    patterns: Vec<PythonExpression>,
) -> Vec<PythonTerm> {
    let patterns = patterns
        .iter()
        .map(|pattern| pattern.expr.to_pattern())
        .collect::<Vec<_>>();
    python_terms(expression.expr.expand_in_patterns(&patterns))
}

pub(super) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_function(wrap_pyfunction!(expand_in_patterns, module)?)?;
    Ok(())
}
