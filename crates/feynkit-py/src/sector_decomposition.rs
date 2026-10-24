//! Optional sector-decomposition entry points on the existing HEPKit owners.
//!
//! The containing extension registers its linked backend without a dependency
//! from FeynKit to FastSecDec or any Python module namespace.

use pyo3::types::PyDict;
use pyo3::{Bound, Py, PyAny, PyRef, PyResult, pymethods};
use symbolica::api::python::{ExtensionPoint, PythonExtensionFunction};

use crate::{graph::PyFeynmanDiagram, integrals::PyIntegralFamily};

/// Optional sector decomposition supplied by a backend linked into the same host.
///
/// Register a static `ExtensionImplementation<PythonExtensionFunction>` during
/// host initialization. The backend retains ownership of argument validation,
/// defaults, and scientific computation; no Python module lookup is performed.
///
/// ```no_run
/// use feynkit_py::SECTOR_DECOMPOSE;
/// use pyo3::types::{PyDict, PyTuple, PyTupleMethods};
/// use pyo3::{Bound, Py, PyAny, PyErr, PyResult};
/// use symbolica::api::python::{ExtensionImplementation, PythonExtensionFunction};
///
/// fn backend(
///     args: &Bound<'_, PyTuple>,
///     _kwargs: Option<&Bound<'_, PyDict>>,
/// ) -> PyResult<Py<PyAny>> {
///     Ok(args.get_item(0)?.unbind())
/// }
/// static BACKEND: ExtensionImplementation<PythonExtensionFunction> =
///     ExtensionImplementation::new("linked-sector-backend", backend);
/// SECTOR_DECOMPOSE.register(&BACKEND)?;
/// # Ok::<(), PyErr>(())
/// ```
pub static SECTOR_DECOMPOSE: ExtensionPoint<PythonExtensionFunction> =
    ExtensionPoint::new("hepkit.sector_decompose");

#[spenso_macros::track_usage(crate::record_usage)]
#[pymethods]
impl PyFeynmanDiagram {
    /// Generate Laurent integrands through the optional sector-decomposition backend.
    ///
    /// The complete native diagram, all keyword values and observer exceptions are
    /// forwarded unchanged to the host's registered ``hepkit.sector_decompose`` backend.
    /// Supply numerical kinematics and external-state bindings explicitly. The
    /// backend uses this diagram's numerator and weights; compilation and numerical
    /// integration remain separate operations on the returned GeneratedIntegral.
    ///
    /// Examples
    /// --------
    /// With a complete diagram, an admitted numerical point and regulator ``eps``:
    ///
    /// >>> from symbolica import E, S
    /// >>> from symbolica.community import hepkit as hep
    /// >>> model = hep.Model.phi3()
    /// >>> diagram = model.process(["phi"], ["phi"]).generate_diagrams(loops=1, max_vertices=2, maximum_bridges=0, progress=None).diagrams[0]
    /// >>> from symbolica.community.hepkit.sector_decomposition import with_diagram_expressions
    /// >>> diagram = with_diagram_expressions(diagram, overall_factor=diagram.overall_factor_expression(evaluate=True))
    /// >>> eps, D = S("eps", "D")
    /// >>> P = hep.Kinematics.external_momentum()
    /// >>> point = hep.Kinematics(D, momenta=[P(0)]).with_scalar_product(P(0), P(0), E("-1"))
    /// >>> values = {model.particle("phi").mass: E("1"), model.parameter("g").symbol: E("1")}
    /// >>> generated = diagram.sector_decompose(regulator=eps, kinematics=point, scalar_values=values, progress=None)
    /// >>> kernels = generated.compile()
    ///
    /// Parameters
    /// ----------
    /// kwargs : keyword arguments
    ///     Arguments of ``sector_decomposition.sector_decompose`` after its input.
    #[pyo3(signature = (**kwargs), text_signature = "($self, *, regulator, kinematics=None, dimension=None, powers=None, numerator=None, scalar_values=None, auxiliary_momenta=None, measure_multiplier=None, max_order=0, coefficient_expansion='physical', observer=None, progress='auto')")]
    fn sector_decompose(
        slf: PyRef<'_, Self>,
        kwargs: Option<&Bound<'_, PyDict>>,
    ) -> PyResult<Py<PyAny>> {
        let py = slf.py();
        SECTOR_DECOMPOSE.call(py, (slf,), kwargs)
    }
}

#[spenso_macros::track_usage(crate::record_usage)]
#[pymethods]
impl PyIntegralFamily {
    /// Generate Laurent integrands for an explicitly specified member of this family.
    ///
    /// Forward to the host's registered ``hepkit.sector_decompose`` backend without changing
    /// any Python values. Supply signed powers in native denominator order and an
    /// already weighted scalar numerator. A completed family's auxiliary slots do
    /// not implicitly acquire power one. The backend owns validation and algebra.
    ///
    /// Examples
    /// --------
    /// For a family with numerical external kinematics and regulator ``eps``:
    ///
    /// >>> from symbolica import E, S
    /// >>> from symbolica.community import hepkit as hep
    /// >>> eps, D, k = S("eps", "D", "k")
    /// >>> point = hep.Kinematics(D, momenta=[k])
    /// >>> family = hep.IntegralFamily([k], [], [point.scalar_product(k, k) - 1], kinematics=point)
    /// >>> generated = family.sector_decompose(regulator=eps, powers=[2], numerator=E("1"), progress=None)
    /// >>> kernels = generated.compile()
    ///
    /// Parameters
    /// ----------
    /// kwargs : keyword arguments
    ///     Arguments of ``sector_decomposition.sector_decompose`` after its input.
    #[pyo3(signature = (**kwargs), text_signature = "($self, *, regulator, kinematics=None, dimension=None, powers=None, numerator=None, scalar_values=None, auxiliary_momenta=None, measure_multiplier=None, max_order=0, coefficient_expansion='physical', observer=None, progress='auto')")]
    fn sector_decompose(
        slf: PyRef<'_, Self>,
        kwargs: Option<&Bound<'_, PyDict>>,
    ) -> PyResult<Py<PyAny>> {
        let py = slf.py();
        SECTOR_DECOMPOSE.call(py, (slf,), kwargs)
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
                max_order: int = 0, coefficient_expansion: str = "physical",
                observer: typing.Optional[typing.Callable[[symbolica.community.hepkit.sector_decomposition.GenerationSnapshot], typing.Optional[bool]]] = None,
                progress: typing.Union[typing.Literal["auto"], typing.Callable[[symbolica.community.hepkit.sector_decomposition.GenerationSnapshot], typing.Optional[bool]], None] = "auto",
            ) -> symbolica.community.hepkit.sector_decomposition.GeneratedIntegral:
                """Generate Laurent integrands using the complete native diagram.

                Examples
                --------
                With a complete diagram and admitted numerical kinematics:

                >>> from symbolica import E, S
                >>> from symbolica.community import hepkit as hep
                >>> model = hep.Model.phi3()
                >>> diagram = model.process(["phi"], ["phi"]).generate_diagrams(loops=1, max_vertices=2, maximum_bridges=0, progress=None).diagrams[0]
                >>> from symbolica.community.hepkit.sector_decomposition import with_diagram_expressions
                >>> diagram = with_diagram_expressions(diagram, overall_factor=diagram.overall_factor_expression(evaluate=True))
                >>> eps, D = S("eps", "D")
                >>> P = hep.Kinematics.external_momentum()
                >>> point = hep.Kinematics(D, momenta=[P(0)]).with_scalar_product(P(0), P(0), E("-1"))
                >>> values = {model.particle("phi").mass: E("1"), model.parameter("g").symbol: E("1")}
                >>> generated = diagram.sector_decompose(regulator=eps, kinematics=point, scalar_values=values, progress=None)
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
                max_order : int
                    Largest signed epsilon power retained.
                coefficient_expansion : str
                    Native physical or package coefficient convention.
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
                max_order: int = 0, coefficient_expansion: str = "physical",
                observer: typing.Optional[typing.Callable[[symbolica.community.hepkit.sector_decomposition.GenerationSnapshot], typing.Optional[bool]]] = None,
                progress: typing.Union[typing.Literal["auto"], typing.Callable[[symbolica.community.hepkit.sector_decomposition.GenerationSnapshot], typing.Optional[bool]], None] = "auto",
            ) -> symbolica.community.hepkit.sector_decomposition.GeneratedIntegral:
                """Generate Laurent integrands for an explicit family member.

                Examples
                --------
                For a numerical family, with powers in native denominator order:

                >>> from symbolica import E, S
                >>> from symbolica.community import hepkit as hep
                >>> eps, D, k = S("eps", "D", "k")
                >>> point = hep.Kinematics(D, momenta=[k])
                >>> family = hep.IntegralFamily([k], [], [point.scalar_product(k, k) - 1], kinematics=point)
                >>> generated = family.sector_decompose(regulator=eps, powers=[2], numerator=E("1"), progress=None)
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
                max_order : int
                    Largest signed epsilon power retained.
                coefficient_expansion : str
                    Native physical or package coefficient convention.
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
    use std::ffi::CString;

    use super::*;
    use feynkit_graph::FeynmanDiagram;
    use feynkit_kinematics::Kinematics;
    use feynkit_model::Model;
    use pyo3::exceptions::PyRuntimeError;
    use pyo3::types::{PyAnyMethods, PyCFunction, PyDictMethods, PyTuple};
    use pyo3::{PyErr, Python};
    use symbolica::api::python::ExtensionImplementation;

    fn observe(
        args: &Bound<'_, PyTuple>,
        kwargs: Option<&Bound<'_, PyDict>>,
    ) -> PyResult<Py<PyAny>> {
        let kwargs = kwargs.expect("the test supplies keyword arguments");
        let observer = kwargs
            .get_item("observer")?
            .expect("the test supplies an observer");
        Ok(observer.call(args, Some(kwargs))?.unbind())
    }

    fn conflicting_backend(
        _args: &Bound<'_, PyTuple>,
        _kwargs: Option<&Bound<'_, PyDict>>,
    ) -> PyResult<Py<PyAny>> {
        Err(PyRuntimeError::new_err(
            "this backend must not be installed",
        ))
    }

    static OBSERVER: ExtensionImplementation<PythonExtensionFunction> =
        ExtensionImplementation::new("observer", observe);
    static CONFLICTING: ExtensionImplementation<PythonExtensionFunction> =
        ExtensionImplementation::new("conflicting", conflicting_backend);

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
            locals
                .set_item(
                    "install_backend",
                    PyCFunction::new_closure(py, None, None, |_, _| {
                        SECTOR_DECOMPOSE.register(&OBSERVER)?;
                        assert!(SECTOR_DECOMPOSE.register(&CONFLICTING).is_err());
                        Ok::<_, PyErr>(())
                    })
                    .unwrap(),
                )
                .unwrap();
            let script = CString::new(r#"
import builtins
import inspect

original_import = builtins.__import__
def reject_symbolica_import(name, *args, **kwargs):
    if name == 'symbolica' or name.startswith('symbolica.'):
        raise AssertionError('canonical Symbolica namespace was imported: ' + name)
    return original_import(name, *args, **kwargs)

calls = []
answer, regulator, powers, progress = (object() for _ in range(4))
def observer(input, **kwargs):
    calls.append((input, kwargs))
    return answer
builtins.__import__ = reject_symbolica_import
try:
    for input in (diagram, family):
        try:
            input.sector_decompose(regulator=regulator)
        except ImportError as caught:
            assert 'hepkit.sector_decompose' in str(caught)
            assert caught.__cause__ is None
        else:
            raise AssertionError('missing optional backend was accepted')
    install_backend()
    install_backend()  # Repeated initialization may register the same function.
    for input in (diagram, family):
        assert input.sector_decompose(regulator=regulator, powers=powers, observer=observer, progress=progress) is answer
        actual, kwargs = calls.pop()
        assert actual is input
        assert kwargs == dict(regulator=regulator, powers=powers, observer=observer, progress=progress)
        signature = inspect.signature(input.sector_decompose)
        assert signature.parameters['regulator'].kind == inspect.Parameter.KEYWORD_ONLY
        assert signature.parameters['max_order'].default == 0
        assert signature.parameters['coefficient_expansion'].default == 'physical'
        assert signature.parameters['progress'].default == 'auto'
    failure = ValueError('backend input or observer failure')
    def fail(*args, **kwargs):
        raise failure
    for input in (diagram, family):
        try:
            input.sector_decompose(regulator=regulator, observer=fail)
        except ValueError as caught:
            assert caught is failure
        else:
            raise AssertionError('backend failure was swallowed')
finally:
    builtins.__import__ = original_import
"#).unwrap();
            py.run(&script, Some(&locals), Some(&locals)).unwrap();
        });
    }
}
