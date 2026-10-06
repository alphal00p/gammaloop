use pyo3::{prelude::*, types::PyDict};
use spynso3::SpensoModule;
use symbolica::api::python::SymbolicaCommunityModule;

#[test]
fn citations_follow_successful_tensor_and_algebra_operations() {
    Python::initialize();
    Python::attach(|py| -> PyResult<()> {
        SpensoModule::initialize(py)?;
        let module = PyModule::new(py, "spenso")?;
        SpensoModule::register_module(&module)?;
        let globals = PyDict::new(py);
        globals.set_item("sp", module)?;
        assert!(SpensoModule::get_citations().is_empty());

        py.run(
            c"
settings = sp.DisplaySettings()
repr(settings)
sp.set_symbolica_rayon_enabled(sp.SymbolicParallelism.Serial)
",
            Some(&globals),
            None,
        )?;
        assert!(SpensoModule::get_citations().is_empty());

        py.run(
            c"
space = sp.Representation.euc(2)
A = sp.TensorName('citation_A')(space, space)
tensor = sp.Tensor.dense(A, [1.0, 2.0, 3.0, 4.0])
assert tensor.shape == (2, 2)
repr(tensor)
",
            Some(&globals),
            None,
        )?;
        let citations = SpensoModule::get_citations();
        assert_eq!(citations.len(), 1);
        assert_eq!(citations[0].id, "10.5281/zenodo.18248388");

        py.run(
            c"
try:
    A.simplify_algebra(gamma_output='invalid')
except ValueError:
    pass
else:
    raise AssertionError('invalid algebra settings succeeded')
",
            Some(&globals),
            None,
        )?;
        assert_eq!(SpensoModule::get_citations().len(), 1);

        py.run(
            c"
network = sp.TensorNetwork(tensor)
network.execute()
assert network.result_tensor()[1, 0] == 3.0
evaluator = tensor.evaluator([], jit_compile=False)
assert evaluator.evaluate([[]])[0][1, 0] == 3.0
",
            Some(&globals),
            None,
        )?;
        let citations = SpensoModule::get_citations();
        assert_eq!(citations.len(), 1);
        assert!(
            citations[0]
                .reasons
                .iter()
                .any(|reason| reason.contains("Tensor-network contractions"))
        );
        assert!(
            citations[0]
                .reasons
                .iter()
                .any(|reason| reason.contains("Numerical evaluation of tensor components"))
        );

        py.run(c"A.contract()", Some(&globals), None)?;
        let citations = SpensoModule::get_citations();
        assert_eq!(citations.len(), 2);
        assert_eq!(citations[1].id, "10.5281/zenodo.18248409");
        assert_eq!(citations[1].reasons.len(), 1);
        assert!(citations[1].reasons[0].contains("Symbolic tensor contractions"));

        py.run(
            c"
spin = sp.Representation.bis(4)
lorentz = sp.Representation.mink(4)
gamma = sp.TensorExpression.dirac_gamma(4)
word = gamma(spin('i'), spin('j'), lorentz('mu')) * gamma(spin('j'), spin('i'), lorentz('nu'))
word.simplify_algebra(gamma=True, color=False)
word.simplify_algebra(gamma=True, color=False)
",
            Some(&globals),
            None,
        )?;
        let citations = SpensoModule::get_citations();
        assert_eq!(citations.len(), 4);
        assert_eq!(citations[2].id, "arXiv:2601.19982");
        assert_eq!(citations[3].id, "arXiv:1707.06453");
        let reasons = &citations[1].reasons;
        assert_eq!(reasons.len(), 3);
        assert!(
            reasons
                .iter()
                .any(|reason| reason.contains("Dirac gamma algebra"))
        );
        assert!(
            !reasons
                .iter()
                .any(|reason| reason.contains("Color algebra"))
        );
        assert!(
            !reasons
                .iter()
                .any(|reason| reason.contains("Levi-Civita identities"))
        );

        py.run(c"A.canonize()", Some(&globals), None)?;
        let citations = SpensoModule::get_citations();
        assert_eq!(citations[1].reasons.len(), 4);
        assert_eq!(
            SpensoModule::get_citations()[1].reasons,
            citations[1].reasons
        );

        py.run(
            c"
f = sp.TensorExpression.color_f(8)
color_word = f('a', 'c', 'd') * f('b', 'c', 'd')
color_word.simplify_algebra(gamma=False, color=True)
color_word.simplify_algebra(gamma=False, color=True)
",
            Some(&globals),
            None,
        )?;
        let citations = SpensoModule::get_citations();
        assert_eq!(citations.len(), 5);
        assert_eq!(citations[2].id, "arXiv:2601.19982");
        assert_eq!(citations[3].id, "arXiv:1707.06453");
        assert_eq!(citations[4].id, "arXiv:hep-ph/9802376");
        assert_eq!(citations[4].reasons, ["Color algebra."]);
        Ok(())
    })
    .unwrap();
}
