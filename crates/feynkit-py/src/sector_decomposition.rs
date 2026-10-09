//! Optional sector-decomposition entry points on the existing HEPKit owners.
//!
//! Resolve the backend at call time: FeynKit does not depend on FastSecDec.

use pyo3::{
    exceptions::{PyImportError, PyModuleNotFoundError},
    prelude::*,
    types::PyDict,
};

use crate::{graph::PyFeynmanDiagram, integrals::PyIntegralFamily};

const MODULE: &str = "symbolica.community.hepkit.sector_decomposition";

fn forward(input: Bound<'_, PyAny>, kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<Py<PyAny>> {
    let py = input.py();
    let backend = py.import(MODULE).map_err(|error| {
        let missing = error.is_instance_of::<PyModuleNotFoundError>(py)
            && error
                .value(py)
                .getattr("name")
                .ok()
                .and_then(|name| name.extract::<String>().ok())
                .is_some_and(|name| name == MODULE || MODULE.starts_with(&(name + ".")));
        if missing {
            let unavailable = PyImportError::new_err(
                "Sector decomposition requires a community wheel with the FastSecDec backend",
            );
            unavailable.set_cause(py, Some(error));
            unavailable
        } else {
            error
        }
    })?;
    Ok(backend
        .getattr("sector_decompose")?
        .call((input,), kwargs)?
        .unbind())
}

#[spenso_macros::track_usage(crate::record_usage)]
#[pymethods]
impl PyFeynmanDiagram {
    /// Generate Laurent integrands through the optional sector-decomposition backend.
    ///
    /// The complete native diagram, all keyword values and observer exceptions are
    /// forwarded unchanged to ``hepkit.sector_decomposition.sector_decompose``.
    /// Supply numerical kinematics and external-state bindings explicitly. The
    /// backend uses this diagram's numerator and weights; compilation and numerical
    /// integration remain separate operations on the returned GeneratedIntegral.
    ///
    /// Examples
    /// --------
    /// With a complete diagram, an admitted numerical point and regulator ``eps``:
    ///
    /// >>> generated = diagram.sector_decompose(regulator=eps, kinematics=kinematics)
    /// >>> kernels = generated.compile()
    ///
    /// Parameters
    /// ----------
    /// kwargs : keyword arguments
    ///     Arguments of ``sector_decomposition.sector_decompose`` after its input.
    #[pyo3(signature = (**kwargs), text_signature = "($self, *, regulator, kinematics=None, dimension=None, powers=None, numerator=None, scalar_values=None, auxiliary_momenta=None, measure_multiplier=None, runtime_parameters=None, model_parameters='runtime', max_order=0, coefficient_expansion='full_expression', mode='symbolic', subtraction='taylor', contour=False, observer=None, progress='auto')")]
    fn sector_decompose(
        slf: PyRef<'_, Self>,
        kwargs: Option<&Bound<'_, PyDict>>,
    ) -> PyResult<Py<PyAny>> {
        let py = slf.py();
        forward(slf.into_pyobject(py)?.into_any(), kwargs)
    }
}

#[spenso_macros::track_usage(crate::record_usage)]
#[pymethods]
impl PyIntegralFamily {
    /// Generate Laurent integrands for an explicitly specified member of this family.
    ///
    /// Forward to ``hepkit.sector_decomposition.sector_decompose`` without changing
    /// any Python values. Supply signed powers in native denominator order and an
    /// already weighted scalar numerator. A completed family's auxiliary slots do
    /// not implicitly acquire power one. The backend owns validation and algebra.
    ///
    /// Examples
    /// --------
    /// For a family with numerical external kinematics and regulator ``eps``:
    ///
    /// >>> generated = family.sector_decompose(regulator=eps, powers=[1, 1], numerator=E("1"))
    /// >>> kernels = generated.compile()
    ///
    /// Parameters
    /// ----------
    /// kwargs : keyword arguments
    ///     Arguments of ``sector_decomposition.sector_decompose`` after its input.
    #[pyo3(signature = (**kwargs), text_signature = "($self, *, regulator, kinematics=None, dimension=None, powers=None, numerator=None, scalar_values=None, auxiliary_momenta=None, measure_multiplier=None, runtime_parameters=None, model_parameters='runtime', max_order=0, coefficient_expansion='full_expression', mode='symbolic', subtraction='taylor', contour=False, observer=None, progress='auto')")]
    fn sector_decompose(
        slf: PyRef<'_, Self>,
        kwargs: Option<&Bound<'_, PyDict>>,
    ) -> PyResult<Py<PyAny>> {
        let py = slf.py();
        forward(slf.into_pyobject(py)?.into_any(), kwargs)
    }
}

// Match the backend's discoverable keyword interface while keeping the runtime
// forwarders opaque: conversion and scientific validation belong to the backend.
#[cfg(feature = "python_stubgen")]
pyo3_stub_gen::inventory::submit! {
    pyo3_stub_gen::derive::gen_methods_from_python! {
        r#"
        import typing
        import symbolica
        import symbolica.community.hepkit.sector_decomposition

        class PyFeynmanDiagram:
            def sector_decompose(self, *, regulator: symbolica.Expression,
                kinematics: pyo3_stub_gen.RustType["Option<crate::kinematics::PyKinematics>"] = None,
                dimension: typing.Optional[symbolica.Expression] = None,
                powers: typing.Optional[typing.Dict[int, int]] = None,
                numerator: None = None,
                scalar_values: typing.Optional[typing.Dict[symbolica.Expression, symbolica.Expression]] = None,
                auxiliary_momenta: typing.Optional[typing.Sequence[symbolica.Expression]] = None,
                measure_multiplier: typing.Optional[symbolica.Expression] = None,
                runtime_parameters: typing.Optional[typing.List[symbolica.Expression]] = None,
                model_parameters: str = "runtime",
                max_order: int = 0, coefficient_expansion: str = "full_expression",
                mode: str = "symbolic", subtraction: str = "taylor", contour: bool = False,
                observer: typing.Optional[typing.Callable[[symbolica.community.hepkit.sector_decomposition.GenerationSnapshot], typing.Optional[bool]]] = None,
                progress: typing.Union[typing.Literal["auto"], typing.Callable[[symbolica.community.hepkit.sector_decomposition.GenerationSnapshot], typing.Optional[bool]], None] = "auto",
            ) -> symbolica.community.hepkit.sector_decomposition.GeneratedIntegral:
                """Generate Laurent integrands using the complete native diagram.

                Examples
                --------
                With a complete diagram and admitted numerical kinematics:

                >>> generated = diagram.sector_decompose(regulator=eps, kinematics=kinematics)
                >>> kernels = generated.compile()

                Parameters
                ----------
                regulator : Expression
                    Dimensional regulator symbol.
                kinematics : Kinematics
                    Numerical external point retaining its symbolic tensor dimension.
                dimension : Expression or None
                    Integration dimension; defaults to 4 - 2*regulator.
                powers : mapping[int, int] or None
                    Positive propagator powers by stable diagram edge ID.
                numerator : None
                    Use the diagram's numerator; explicit overrides are rejected.
                scalar_values : mapping[Expression, Expression] or None
                    Explicit masses, couplings and invariant substitutions.
                auxiliary_momenta : sequence[Expression] or None
                    Additional external vector heads, such as polarizations.
                measure_multiplier : Expression or None
                    Explicit multiplicative measure convention, applied once.
                runtime_parameters : list[Expression] or None
                    Scalar inputs retained for binding after generation.
                model_parameters : str
                    Retain model parameters at runtime or substitute their supplied values.
                max_order : int
                    Largest signed epsilon power retained.
                coefficient_expansion : str
                    Native physical or package coefficient convention.
                mode : str
                    Symbolic or numerical-dual coefficient construction.
                subtraction : str
                    Native Taylor or IBP endpoint subtraction.
                contour : bool
                    Generate optional contour-deformation capability for later binding.
                observer : callable or None
                    Native generation events; False cancels at an event boundary.
                progress : "auto", callable or None
                    Automatic marimo display unless observer is supplied; None disables it.
                    A callable receives every native event after observer and may cancel.
                """

        "#
    }
}

#[cfg(feature = "python_stubgen")]
pyo3_stub_gen::inventory::submit! {
    pyo3_stub_gen::derive::gen_methods_from_python! {
        r#"
        import typing
        import symbolica
        import symbolica.community.hepkit.sector_decomposition

        class PyIntegralFamily:
            def sector_decompose(self, *, regulator: symbolica.Expression,
                kinematics: pyo3_stub_gen.RustType["Option<crate::kinematics::PyKinematics>"] = None,
                dimension: typing.Optional[symbolica.Expression] = None,
                powers: typing.Optional[typing.Sequence[int]] = None,
                numerator: typing.Optional[symbolica.Expression] = None,
                scalar_values: typing.Optional[typing.Dict[symbolica.Expression, symbolica.Expression]] = None,
                auxiliary_momenta: typing.Optional[typing.Sequence[symbolica.Expression]] = None,
                measure_multiplier: typing.Optional[symbolica.Expression] = None,
                runtime_parameters: typing.Optional[typing.List[symbolica.Expression]] = None,
                model_parameters: str = "runtime",
                max_order: int = 0, coefficient_expansion: str = "full_expression",
                mode: str = "symbolic", subtraction: str = "taylor", contour: bool = False,
                observer: typing.Optional[typing.Callable[[symbolica.community.hepkit.sector_decomposition.GenerationSnapshot], typing.Optional[bool]]] = None,
                progress: typing.Union[typing.Literal["auto"], typing.Callable[[symbolica.community.hepkit.sector_decomposition.GenerationSnapshot], typing.Optional[bool]], None] = "auto",
            ) -> symbolica.community.hepkit.sector_decomposition.GeneratedIntegral:
                """Generate Laurent integrands for an explicit family member.

                Examples
                --------
                For a numerical family, with powers in native denominator order:

                >>> generated = family.sector_decompose(regulator=eps, powers=[1, 1], numerator=E("1"))
                >>> kernels = generated.compile()

                Parameters
                ----------
                regulator : Expression
                    Dimensional regulator symbol.
                kinematics : Kinematics or None
                    Explicit point; None retains the family's scoped kinematics.
                dimension : Expression or None
                    Integration dimension; defaults to 4 - 2*regulator.
                powers : sequence[int]
                    Required signed powers; zero omits a slot and negative moves it to the numerator.
                numerator : Expression
                    Required scalar numerator including explicit physical weights once.
                scalar_values : mapping[Expression, Expression] or None
                    Explicit scalar substitutions applied consistently to the family.
                auxiliary_momenta : sequence[Expression] or None
                    Additional external vector heads, such as polarizations.
                measure_multiplier : Expression or None
                    Explicit multiplicative measure convention, applied once.
                runtime_parameters : list[Expression] or None
                    Scalar inputs retained for binding after generation.
                model_parameters : str
                    Retain model parameters at runtime or substitute their supplied values.
                max_order : int
                    Largest signed epsilon power retained.
                coefficient_expansion : str
                    Native physical or package coefficient convention.
                mode : str
                    Symbolic or numerical-dual coefficient construction.
                subtraction : str
                    Native Taylor or IBP endpoint subtraction.
                contour : bool
                    Generate optional contour-deformation capability for later binding.
                observer : callable or None
                    Native generation events; False cancels at an event boundary.
                progress : "auto", callable or None
                    Automatic marimo display unless observer is supplied; None disables it.
                    A callable receives every native event after observer and may cancel.
                """
        "#
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use feynkit_graph::FeynmanDiagram;
    use feynkit_kinematics::Kinematics;
    use feynkit_model::Model;
    use std::ffi::CString;

    #[test]
    fn forwards_existing_owner_and_keyword_identities_and_exceptions() {
        Python::initialize();
        Python::attach(|py| {
            let diagram = FeynmanDiagram::from_dot(
                Model::phi4(),
                r#"digraph { a -> a [particle="phi"]; a -> a [particle="phi"]; }"#,
            )
            .unwrap();
            let family = diagram.propagator_family(&Kinematics::new()).unwrap();
            let locals = PyDict::new(py);
            locals
                .set_item(
                    "diagram",
                    Py::new(py, PyFeynmanDiagram::from(diagram)).unwrap(),
                )
                .unwrap();
            locals
                .set_item(
                    "family",
                    Py::new(py, PyIntegralFamily { inner: family }).unwrap(),
                )
                .unwrap();
            locals.set_item("backend_name", MODULE).unwrap();
            let script = CString::new(r#"
import inspect
import sys
import types

absent = object()
previous = sys.modules.get(backend_name, absent)
# The native owner unit test must not require an installed community wheel.
created_parents = []
for name in ('symbolica', 'symbolica.community', 'symbolica.community.hepkit'):
    if name not in sys.modules:
        parent = types.ModuleType(name)
        parent.__path__ = []
        sys.modules[name] = parent
        created_parents.append(name)
backend = types.ModuleType(backend_name)
calls = []
answer, regulator, observer, powers, progress, runtime_parameters = (object() for _ in range(6))
def capture(input, **kwargs):
    calls.append((input, kwargs))
    return answer
backend.sector_decompose = capture
sys.modules[backend_name] = backend
try:
    for input in (diagram, family):
        assert input.sector_decompose(regulator=regulator, powers=powers, observer=observer, progress=progress, runtime_parameters=runtime_parameters, model_parameters="runtime", mode="numerical_dual", subtraction="ibp", contour=True) is answer
        actual, kwargs = calls.pop()
        assert actual is input
        assert kwargs == dict(regulator=regulator, powers=powers, observer=observer, progress=progress, runtime_parameters=runtime_parameters, model_parameters="runtime", mode="numerical_dual", subtraction="ibp", contour=True)
        signature = inspect.signature(input.sector_decompose)
        assert signature.parameters['regulator'].kind == inspect.Parameter.KEYWORD_ONLY
        assert signature.parameters['max_order'].default == 0
        assert signature.parameters['coefficient_expansion'].default == 'full_expression'
        assert signature.parameters['runtime_parameters'].default is None
        assert signature.parameters['model_parameters'].default == 'runtime'
        assert signature.parameters['mode'].default == 'symbolic'
        assert signature.parameters['subtraction'].default == 'taylor'
        assert signature.parameters['contour'].default is False
        assert signature.parameters['progress'].default == 'auto'
    failure = ValueError('backend input or observer failure')
    def fail(*args, **kwargs):
        raise failure
    backend.sector_decompose = fail
    for input in (diagram, family):
        try:
            input.sector_decompose(regulator=regulator)
        except ValueError as caught:
            assert caught is failure
        else:
            raise AssertionError('backend failure was swallowed')
    sys.modules[backend_name] = None
    try:
        family.sector_decompose(regulator=regulator)
    except ImportError as caught:
        assert 'FastSecDec backend' in str(caught)
        assert isinstance(caught.__cause__, ModuleNotFoundError)
    else:
        raise AssertionError('missing optional backend was accepted')
finally:
    if previous is absent:
        del sys.modules[backend_name]
    else:
        sys.modules[backend_name] = previous
    for name in reversed(created_parents):
        del sys.modules[name]
"#).unwrap();
            py.run(&script, Some(&locals), Some(&locals)).unwrap();
        });
    }
}
