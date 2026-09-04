use std::ffi::CString;

use pyo3::{
    exceptions::{PyTypeError, PyValueError},
    prelude::*,
    types::{PyDict, PyList, PyModule},
};
use spenso::structure::{
    representation::{LibraryRep, Representation},
    slot::IsAbstractSlot,
};
use spynso3::{
    SpensoModule,
    structure::{SpensoRepresentation, SpensoSlot},
};
use symbolica::{
    api::python::{PythonExpression, SymbolicaCommunityModule},
    atom::{Atom, AtomView},
    symbol,
};

fn spenso_module<'py>(py: Python<'py>) -> PyResult<Bound<'py, PyModule>> {
    SpensoModule::initialize(py)?;
    let module = PyModule::new(py, "spenso")?;
    SpensoModule::register_module(&module)?;
    Ok(module)
}

fn logical_representations(
    expression: &Bound<'_, PyAny>,
) -> PyResult<Vec<Representation<LibraryRep>>> {
    expression
        .getattr("interface")?
        .extract::<Vec<SpensoRepresentation>>()
        .map(|representations| {
            representations
                .into_iter()
                .map(|representation| representation.representation)
                .collect()
        })
}

fn logical_slot_representations(
    expression: &Bound<'_, PyAny>,
) -> PyResult<Vec<Representation<LibraryRep>>> {
    expression
        .getattr("interface")?
        .extract::<Vec<SpensoSlot>>()
        .map(|slots| slots.into_iter().map(|slot| slot.slot.rep()).collect())
}

fn symbolic_argument_count(expression: &Bound<'_, PyAny>) -> PyResult<usize> {
    let expression = expression
        .call_method0("to_expression")?
        .extract::<PythonExpression>()?;
    let AtomView::Fun(function) = expression.expr.as_view() else {
        panic!("a predefined tensor factory did not produce a function")
    };
    Ok(function.iter().count())
}

#[test]
fn typed_factories_expose_logical_interfaces_without_key_arguments() {
    Python::initialize();
    Python::attach(|py| -> PyResult<()> {
        let module = spenso_module(py)?;
        let tensor_expression = module.getattr("TensorExpression")?;
        let representation = module.getattr("Representation")?;

        let coad8 = representation
            .call_method1("coad", (8,))?
            .extract::<SpensoRepresentation>()?
            .representation;
        let cof_object = representation.call_method1("cof", (3,))?;
        let cof3 = cof_object.extract::<SpensoRepresentation>()?.representation;
        let coaf3 = cof_object
            .call_method0("dual")?
            .extract::<SpensoRepresentation>()?
            .representation;

        let f = tensor_expression.call_method1("f", (8,))?;
        assert_eq!(logical_representations(&f)?, vec![coad8, coad8, coad8]);
        assert_eq!(symbolic_argument_count(&f)?, 3);

        let t = tensor_expression.call_method1("t", (8, 3))?;
        let expected_t = vec![coad8, cof3, coaf3];
        assert_eq!(logical_representations(&t)?, expected_t);
        assert_eq!(symbolic_argument_count(&t)?, 3);

        let indexed = t.call1(("a", "i", "j"))?;
        assert_eq!(logical_slot_representations(&indexed)?, expected_t);
        let error = t
            .call1(("a", "i"))
            .expect_err("the typed tensor must retain its fixed arity");
        assert!(error.is_instance_of::<PyValueError>(py));
        Ok(())
    })
    .unwrap();
}

#[test]
fn predefined_names_are_heads_and_user_names_remain_constructible() {
    Python::initialize();
    Python::attach(|py| -> PyResult<()> {
        let module = spenso_module(py)?;
        let tensor_name = module.getattr("TensorName")?;
        let representation = module
            .getattr("Representation")?
            .call_method1("euc", (2,))?;

        let predefined = tensor_name.call_method0("t")?;
        let error = predefined
            .call1(())
            .expect_err("predefined heads must use their typed factory");
        assert!(error.is_instance_of::<PyTypeError>(py));
        assert!(error.to_string().contains("TensorExpression.t"));

        let globals = PyDict::new(py);
        globals.set_item("spenso", &module)?;
        let source = CString::new("user_name = spenso.TensorName('A')").unwrap();
        py.run(source.as_c_str(), Some(&globals), None)?;
        let user_name = globals
            .get_item("user_name")?
            .expect("the Python snippet must define user_name");
        let expression = user_name.call1((&representation,))?;
        assert_eq!(expression.getattr("rank")?.extract::<usize>()?, 1);
        Ok(())
    })
    .unwrap();
}

#[test]
fn public_patterns_match_canonicalized_builtin_tensors() {
    Python::initialize();
    Python::attach(|py| -> PyResult<()> {
        let module = spenso_module(py)?;
        let tensor_expression = module.getattr("TensorExpression")?;
        let tensor_pattern = module.getattr("TensorPattern")?;

        let pattern = tensor_pattern.call_method1(
            "t",
            (
                PythonExpression::from(Atom::var(symbol!("factory_adjoint_dimension_"))),
                PythonExpression::from(Atom::var(symbol!("factory_fundamental_dimension_"))),
                PythonExpression::from(Atom::var(symbol!("factory_adjoint_index_"))),
                PythonExpression::from(Atom::var(symbol!("factory_fundamental_index_"))),
                PythonExpression::from(Atom::var(symbol!("factory_antifundamental_index_"))),
            ),
        )?;
        let target = tensor_expression
            .call_method1("t", (8, 3))?
            .call1(("a", "i", "j"))?
            .call_method0("to_expression")?;
        let replaced = target
            .call_method1("replace", (pattern, 37))?
            .extract::<PythonExpression>()?;
        assert_eq!(replaced.expr, Atom::num(37));

        let ports = PyList::new(
            py,
            [PythonExpression::from(Atom::var(symbol!(
                "factory_pattern_ports___"
            )))],
        )?;
        let kwargs = PyDict::new(py);
        kwargs.set_item("ports", ports)?;
        let any_tensor =
            tensor_pattern.call_method("any", ("factory_pattern_head_",), Some(&kwargs))?;
        let replaced = target
            .call_method1("replace", (any_tensor, 41))?
            .extract::<PythonExpression>()?;
        assert_eq!(replaced.expr, Atom::num(41));
        Ok(())
    })
    .unwrap();
}

#[test]
fn unresolved_factory_ports_survive_composition() {
    Python::initialize();
    Python::attach(|py| -> PyResult<()> {
        let module = spenso_module(py)?;
        let tensor_expression = module.getattr("TensorExpression")?;

        let generator = tensor_expression.call_method1("t", (8, 3))?;
        let gamma = tensor_expression.call_method1("gamma", (4,))?;
        let outer = generator.call_method1("outer", (&gamma,))?;
        let indexed = outer.call1(("a", "i", "j", "mu", "r", "s"))?;
        assert_eq!(indexed.getattr("rank")?.extract::<usize>()?, 6);

        let minkowski = module
            .getattr("Representation")?
            .call_method1("mink", (4,))?;
        let globals = PyDict::new(py);
        globals.set_item("spenso", &module)?;
        let source =
            CString::new("factory_vector = spenso.TensorName('factory_contraction_vector')")
                .unwrap();
        py.run(source.as_c_str(), Some(&globals), None)?;
        let vector = globals
            .get_item("factory_vector")?
            .expect("the Python snippet must define factory_vector")
            .call1((&minkowski,))?;
        let kwargs = PyDict::new(py);
        kwargs.set_item("left", 0)?;
        kwargs.set_item("right", 0)?;
        let contracted = gamma.call_method("contract", (vector,), Some(&kwargs))?;
        let indexed = contracted.call1(("r", "s"))?;
        assert_eq!(indexed.getattr("rank")?.extract::<usize>()?, 2);

        let euclidean = module
            .getattr("Representation")?
            .call_method1("euc", (4,))?;
        let flat = tensor_expression.call_method1("flat", (&euclidean,))?;
        let scaled = flat.call_method1("__rmul__", (2,))?;
        assert_eq!(
            scaled
                .call1(("i", "j"))?
                .getattr("rank")?
                .extract::<usize>()?,
            2
        );
        assert!(
            flat.call_method0("trace")?
                .getattr("is_scalar")?
                .is_truthy()?
        );
        let outer = gamma.call_method1("outer", (&flat,))?;
        assert_eq!(
            outer
                .call1(("mu", "r", "s", "i", "j"))?
                .getattr("rank")?
                .extract::<usize>()?,
            5
        );
        Ok(())
    })
    .unwrap();
}
