use symbolica::{api::python::PythonExpression, atom::Atom};

pub(crate) type PythonTerm = (PythonExpression, PythonExpression);

pub(crate) fn python_terms(terms: Vec<(Atom, Atom)>) -> Vec<PythonTerm> {
    terms
        .into_iter()
        .map(|(structure, coefficient)| (structure.into(), coefficient.into()))
        .collect()
}
