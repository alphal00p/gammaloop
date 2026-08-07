use pyo3::{Bound, PyResult, pyfunction, types::PyModule, types::PyModuleMethods, wrap_pyfunction};
#[cfg(feature = "python_stubgen")]
use pyo3_stub_gen::derive::gen_stub_pyfunction;
use symbolica::{api::python::PythonExpression, atom::AtomCore};

use crate::{color::ColorSimplifier, selective_expand::SelectiveExpand};

type PythonTerm = (PythonExpression, PythonExpression);

fn python_terms(terms: Vec<(symbolica::atom::Atom, symbolica::atom::Atom)>) -> Vec<PythonTerm> {
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
/// Selectively expands around the supplied expression patterns.
///
/// Each result is returned as `(structure, coefficient)` so callers can retain
/// the same factorization boundary as the Rust API.
fn expand_in_patterns(
    expression: &PythonExpression,
    patterns: Vec<PythonExpression>,
) -> Vec<PythonTerm> {
    let patterns = patterns
        .iter()
        .map(|pattern| pattern.expr.to_pattern())
        .collect::<Vec<_>>();
    python_terms(expression.expr.expand_in_patterns(&patterns))
}

macro_rules! term_expansion {
    ($name:ident, $method:ident) => {
        #[cfg_attr(
            feature = "python_stubgen",
            gen_stub_pyfunction(module = "symbolica.community.idenso")
        )]
        #[pyfunction]
        fn $name(expression: &PythonExpression) -> Vec<PythonTerm> {
            python_terms(SelectiveExpand::$method(&expression.expr))
        }
    };
}

term_expansion!(expand_mink_terms, expand_mink);
term_expansion!(expand_bis_terms, expand_bis);
term_expansion!(expand_mink_bis_terms, expand_mink_bis);
term_expansion!(expand_metrics_terms, expand_metrics);

#[cfg_attr(
    feature = "python_stubgen",
    gen_stub_pyfunction(module = "symbolica.community.idenso")
)]
#[pyfunction]
fn expand_color_terms(expression: &PythonExpression) -> Vec<PythonTerm> {
    python_terms(ColorSimplifier::expand_color(&expression.expr))
}

pub(super) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_function(wrap_pyfunction!(expand_in_patterns, module)?)?;
    module.add_function(wrap_pyfunction!(expand_mink_terms, module)?)?;
    module.add_function(wrap_pyfunction!(expand_bis_terms, module)?)?;
    module.add_function(wrap_pyfunction!(expand_mink_bis_terms, module)?)?;
    module.add_function(wrap_pyfunction!(expand_metrics_terms, module)?)?;
    module.add_function(wrap_pyfunction!(expand_color_terms, module)?)?;
    Ok(())
}
