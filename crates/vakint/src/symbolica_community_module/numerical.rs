//! Preserve input precision while adapting Python numbers to the shared evaluator.

use std::collections::BTreeSet;

use pyo3::{
    exceptions::{PyTypeError, PyValueError},
    prelude::*,
    types::{
        PyBool, PyComplex, PyComplexMethods, PyFloat, PyFloatMethods, PyInt, PyMapping,
        PyMappingMethods,
    },
};
use symbolica::{
    api::python::PythonExpression,
    atom::{Atom, AtomCore},
    domains::{
        float::{Complex, Float, PythonComplexFloat, PythonFloat},
        integer::Integer,
    },
    id::Replacement,
};

use crate::{Momentum, NumericalEvaluationResult, Vakint, VakintSettings, precision};

use super::vakint_to_python_error;

enum Input {
    Expression(Atom),
    Floating(Complex<Float>),
}

impl Input {
    fn extract(value: &Bound<'_, PyAny>) -> PyResult<Self> {
        if let Ok(value) = value.extract::<PyRef<PythonComplexFloat>>() {
            return Ok(Self::Floating(value.0.clone()));
        }
        if let Ok(value) = value.extract::<PyRef<PythonFloat>>() {
            return Ok(Self::Floating(Complex::new(
                value.0.clone(),
                Float::new(value.0.prec()),
            )));
        }
        if value.is_instance_of::<PyBool>() {
            return Err(PyTypeError::new_err(
                "boolean values are not numerical parameters",
            ));
        }
        if let Ok(value) = value.cast::<PyFloat>() {
            return Ok(Self::Floating(Complex::new(
                Float::with_val(53, value.value()),
                Float::new(53),
            )));
        }
        if let Ok(value) = value.cast::<PyComplex>() {
            return Ok(Self::Floating(Complex::new(
                Float::with_val(53, value.real()),
                Float::with_val(53, value.imag()),
            )));
        }
        if value.is_instance_of::<PyInt>() {
            return Ok(Self::Expression(Atom::num(value.extract::<Integer>()?)));
        }
        if let Ok(value) = value.extract::<PythonExpression>() {
            return Ok(Self::Expression(value.expr));
        }
        Err(PyTypeError::new_err(
            "expected Symbolica Float, ComplexFloat or Expression, or a Python int, float or complex",
        ))
    }
}

/// Only input adaptation lives here; contraction and tracked arithmetic stay in Rust's core.
pub(super) fn evaluate(
    engine: &Vakint,
    settings: &VakintSettings,
    expression: &Atom,
    params: &Bound<'_, PyMapping>,
    externals: Option<&Bound<'_, PyMapping>>,
) -> PyResult<(NumericalEvaluationResult, Option<NumericalEvaluationResult>)> {
    let data = engine
        .rustred_data()
        .map_err(|error| vakint_to_python_error(error.into()))?;
    let mut floating = ahash::HashMap::default();
    let mut exact = Vec::new();
    let mut identities = BTreeSet::new();
    let epsilon = crate::vakint_symbol!(&settings.epsilon_symbol);
    for item in params.items()?.try_iter()? {
        let (name, value): (String, Bound<'_, PyAny>) = item?.extract()?;
        let identity = crate::vakint_parse!(&name).map_err(PyValueError::new_err)?;
        crate::inspection::validate_parameter(identity.as_view(), Some(&data.masters))
            .map_err(vakint_to_python_error)?;
        if !identities.insert(identity.clone()) {
            return Err(PyValueError::new_err(format!(
                "duplicate canonical parameter {name}"
            )));
        }
        match Input::extract(&value)? {
            Input::Expression(value) => {
                let value = precision::legacy::stage_table(
                    identity.as_view(),
                    value.as_view(),
                    epsilon,
                    "Input",
                )
                .map_err(vakint_to_python_error)?;
                exact.push((identity, value));
            }
            Input::Floating(value) => {
                floating.insert(name, value);
            }
        }
    }
    let mut momenta = ahash::HashMap::default();
    if let Some(externals) = externals {
        for item in externals.items()?.try_iter()? {
            let (id, values): (usize, Bound<'_, PyAny>) = item?.extract()?;
            i64::try_from(id).map_err(|_| {
                PyValueError::new_err("external momentum ID exceeds the supported integer range")
            })?;
            let values = values.try_iter()?.take(5).collect::<PyResult<Vec<_>>>()?;
            if values.len() != 4 {
                return Err(PyValueError::new_err(format!(
                    "external momentum {id} requires four components (time, x, y, z)"
                )));
            }
            let mut components = Vec::with_capacity(4);
            for (component, value) in values.iter().enumerate() {
                components.push(match Input::extract(value)? {
                    Input::Floating(value) => value,
                    Input::Expression(value) => {
                        let key = precision::external_component(id, component);
                        let value = precision::legacy::stage_table(
                            key.as_view(),
                            value.as_view(),
                            epsilon,
                            "Input",
                        )
                        .map_err(vakint_to_python_error)?;
                        exact.push((key, value));
                        Complex::new(Float::new(53), Float::new(53))
                    }
                });
            }
            let [time, x, y, z]: [Complex<Float>; 4] = components
                .try_into()
                .map_err(|_| PyValueError::new_err("external momentum requires four components"))?;
            momenta.insert(id, Momentum::Complex((time, x, y, z)));
        }
    }
    let (expression, mut parameters) = precision::evaluation_inputs(
        expression.as_view(),
        epsilon,
        &ahash::HashMap::default(),
        &floating,
        externals.map(|_| &momenta),
    )
    .map_err(vakint_to_python_error)?;
    let replacements = exact
        .into_iter()
        .map(|(source, value)| {
            parameters.remove(&source);
            Replacement::new(source, value)
        })
        .collect::<Vec<_>>();
    let expression = expression.replace_multiple(&replacements);
    precision::evaluate_with_errors(
        expression.as_view(),
        &parameters,
        settings,
        Some(&data.masters),
    )
    .map_err(vakint_to_python_error)
}
