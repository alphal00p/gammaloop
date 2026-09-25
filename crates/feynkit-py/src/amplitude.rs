use std::{collections::BTreeMap, sync::Arc};

use feynkit_amplitude::{Amplitude, AmplitudeLeg, AmplitudeOptions, SquaredAmplitude};
use feynkit_graph::ExternalState;
use pyo3::{prelude::*, types::PyModule};
use spynso3::{expression::TensorExpression, structure::SpensoSlot};
use symbolica::{
    api::python::{ConvertibleToExpression, PythonExpression},
    atom::Atom,
};

#[cfg(feature = "python_stubgen")]
use pyo3_stub_gen::derive::{gen_stub_pyclass, gen_stub_pymethods};

use crate::{error, graph::PyFeynmanDiagram, model::PyParticle};

/// One physical external state of a symbolic amplitude.
///
/// Examples
/// --------
/// >>> leg = amplitude.legs[0]
/// >>> leg.particle.spin_sum(leg.momentum, leg.tensor_index, S("conjugate_index"))
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    name = "AmplitudeLeg",
    module = "symbolica.community.feynkit",
    frozen,
    from_py_object
)]
#[derive(Clone)]
pub struct PyAmplitudeLeg {
    inner: AmplitudeLeg,
    model: Arc<feynkit_model::Model>,
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[pymethods]
impl PyAmplitudeLeg {
    /// Stable external-leg label shared by every diagram.
    #[getter]
    fn index(&self) -> usize {
        self.inner.index
    }
    /// Model particle, including particle/antiparticle identity.
    #[getter]
    fn particle(&self) -> PyParticle {
        PyParticle::new(self.inner.particle, self.model.clone())
    }
    /// Whether the state is incoming or outgoing.
    #[getter]
    fn state(&self) -> &'static str {
        match self.inner.state {
            ExternalState::Incoming => "incoming",
            ExternalState::Outgoing => "outgoing",
        }
    }
    /// Unindexed physical momentum P(index).
    #[getter]
    fn momentum(&self) -> PythonExpression {
        PythonExpression {
            expr: self.inner.momentum(),
        }
    }
    /// Bare label shared by this leg's spin and color slots.
    #[getter]
    fn tensor_index(&self) -> PythonExpression {
        PythonExpression {
            expr: self.inner.index_atom(),
        }
    }
    /// Typed open spin/color slots attached to this physical state.
    #[getter]
    fn slots(&self) -> Vec<SpensoSlot> {
        self.inner
            .slots
            .iter()
            .map(|&slot| SpensoSlot { slot })
            .collect()
    }
    /// Describe the state and its physical leg identity.
    ///
    /// Examples
    /// --------
    /// >>> print(amplitude.legs[0])
    fn __repr__(&self) -> String {
        format!(
            "AmplitudeLeg(index={}, particle={:?}, state={:?})",
            self.inner.index,
            self.model.particle_by_id(self.inner.particle).unwrap().name,
            self.state()
        )
    }
}

/// A coherent sum of amputated, unintegrated Feynman-diagram operators.
///
/// Diagrams must describe the same external states in the same model. Named
/// couplings are expanded; model-declared real parameters and physical momenta
/// are real under conjugation. Complex parameters remain complex. Quadratic
/// denominators use the graph convention without widths or an i0 prescription.
/// Already-sewn forward diagrams are rejected.
///
/// Examples
/// --------
/// >>> amplitude = Amplitude(generated.diagrams)
/// >>> operator = amplitude.expression()
/// >>> conjugate = amplitude.conjugate().expression()
/// >>> squared = amplitude.squared().sum_spins(average_initial=True).sum_colors()
/// >>> scalar = squared.expression().simplify_gamma().simplify_color()
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    name = "Amplitude",
    module = "symbolica.community.feynkit",
    frozen,
    from_py_object
)]
#[derive(Clone)]
pub struct PyAmplitude {
    pub(crate) inner: Amplitude,
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[pymethods]
impl PyAmplitude {
    /// Align external ports and sum the weighted diagram operators.
    ///
    /// Examples
    /// --------
    /// >>> amplitude = Amplitude(generated.diagrams)
    ///
    /// Parameters
    /// ----------
    /// diagrams : list[FeynmanDiagram]
    ///     Nonempty collection of complete, unsewn diagrams.
    /// dimension : int or Expression, optional
    ///     Lorentz dimension, default four. Bispinor spaces retain dimension four.
    /// real : list[Expression] or None, optional
    ///     Additional scalar expressions assumed real under conjugation.
    #[new]
    #[pyo3(signature = (diagrams, *, dimension=None, real=None))]
    pub(crate) fn new(
        diagrams: Vec<PyFeynmanDiagram>,
        dimension: Option<ConvertibleToExpression>,
        real: Option<Vec<PythonExpression>>,
    ) -> PyResult<Self> {
        let mut options = AmplitudeOptions::default();
        if let Some(dimension) = dimension {
            options.dimension = dimension
                .to_expression()
                .expr
                .as_view()
                .try_into()
                .map_err(|_| {
                    error::AmplitudeError::new_err(
                        "dimension must be a positive integer or a symbol",
                    )
                })?;
        }
        options.real = real
            .unwrap_or_default()
            .into_iter()
            .map(|r| r.expr)
            .collect();
        if diagrams.iter().any(|d| !d.is_whole_diagram()) {
            return Err(error::AmplitudeError::new_err(
                "an amplitude requires complete diagrams, not selected regions",
            ));
        }
        Ok(Self {
            inner: Amplitude::new(diagrams.into_iter().map(|d| d.inner), options)
                .map_err(error::amplitude)?,
        })
    }

    /// Construct an amplitude from one complete diagram.
    ///
    /// Examples
    /// --------
    /// >>> amplitude = Amplitude.from_diagram(diagram)
    ///
    /// Parameters
    /// ----------
    /// diagram : FeynmanDiagram
    ///     Complete, unsewn source diagram.
    /// dimension : int or Expression, optional
    ///     Lorentz dimension, default four.
    /// real : list[Expression] or None, optional
    ///     Additional scalar reality assumptions.
    #[staticmethod]
    #[pyo3(signature = (diagram, *, dimension=None, real=None))]
    fn from_diagram(
        diagram: PyFeynmanDiagram,
        dimension: Option<ConvertibleToExpression>,
        real: Option<Vec<PythonExpression>>,
    ) -> PyResult<Self> {
        Self::new(vec![diagram], dimension, real)
    }

    /// Source diagrams, retaining weights, routing, and graph provenance.
    #[getter]
    fn diagrams(&self) -> Vec<PyFeynmanDiagram> {
        self.inner
            .diagrams()
            .iter()
            .map(|d| (**d).clone().into())
            .collect()
    }
    /// Physical external states in increasing external-label order.
    #[getter]
    fn legs(&self) -> Vec<PyAmplitudeLeg> {
        self.inner
            .legs()
            .iter()
            .map(|l| PyAmplitudeLeg {
                inner: l.clone(),
                model: self.inner.diagrams()[0].model_arc(),
            })
            .collect()
    }
    /// Individual weighted operators, with aligned external tensor ports.
    #[getter]
    fn terms(&self, py: Python<'_>) -> PyResult<Vec<Py<TensorExpression>>> {
        self.inner
            .terms()
            .iter()
            .map(|t| TensorExpression::from_atom_interface(py, t.clone(), None))
            .collect()
    }
    /// Whether this amplitude is the physical adjoint of its source diagrams.
    #[getter]
    fn is_conjugated(&self) -> bool {
        self.inner.is_conjugated()
    }

    /// Return the complete operator as a Spenso tensor expression.
    ///
    /// Examples
    /// --------
    /// >>> operator = amplitude.expression().factor()
    fn expression(&self, py: Python<'_>) -> PyResult<Py<TensorExpression>> {
        TensorExpression::from_atom_interface(py, self.inner.expression(), None)
    }
    /// Conjugate scalar coefficients, color tensors, and Dirac chains.
    ///
    /// Physical external-leg labels are preserved across different fermion pairings.
    ///
    /// Examples
    /// --------
    /// >>> adjoint = amplitude.conjugate()
    fn conjugate(&self) -> PyResult<Self> {
        Ok(Self {
            inner: self.inner.conjugate().map_err(error::amplitude)?,
        })
    }
    /// Form the coherent square, including all interferences, with distinct ports.
    ///
    /// Spin/color sums and initial-state averages remain explicit operations.
    ///
    /// Examples
    /// --------
    /// >>> squared = amplitude.squared().sum_spins().sum_colors()
    fn squared(&self) -> PyResult<PySquaredAmplitude> {
        Ok(PySquaredAmplitude {
            inner: self.inner.squared().map_err(error::amplitude)?,
        })
    }
    /// Describe the retained diagrams, external states, and conjugation state.
    ///
    /// Examples
    /// --------
    /// >>> print(amplitude)
    fn __repr__(&self) -> String {
        format!(
            "Amplitude(diagrams={}, external_legs={}, conjugated={})",
            self.inner.diagrams().len(),
            self.inner.legs().len(),
            self.inner.is_conjugated()
        )
    }
    /// Render the operator using Spenso's existing tensor printer.
    ///
    /// Examples
    /// --------
    /// >>> amplitude  # notebook output
    fn _repr_html_(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        Ok(self
            .expression(py)?
            .bind(py)
            .call_method0("_repr_html_")?
            .unbind())
    }
}

/// A coherent amplitude square with independent ket and bra tensor indices.
///
/// State sums return new objects; an already-summed leg raises AmplitudeError.
/// No phase-space integration, symmetry factor, flux, or state average is implicit.
///
/// Examples
/// --------
/// >>> squared = Amplitude(generated.diagrams).squared()
/// >>> unpolarized = squared.sum_spins(average_initial=True).sum_colors(average_initial=True)
/// >>> tensor = unpolarized.expression()
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    name = "SquaredAmplitude",
    module = "symbolica.community.feynkit",
    frozen,
    from_py_object
)]
#[derive(Clone)]
pub struct PySquaredAmplitude {
    inner: SquaredAmplitude,
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[pymethods]
impl PySquaredAmplitude {
    /// Construct the coherent square of a single unsewn diagram.
    ///
    /// Examples
    /// --------
    /// >>> squared = SquaredAmplitude.from_diagram(diagram)
    ///
    /// Parameters
    /// ----------
    /// diagram : FeynmanDiagram
    ///     Complete, unsewn source diagram.
    /// dimension : int or Expression, optional
    ///     Lorentz dimension, default four.
    /// real : list[Expression] or None, optional
    ///     Additional scalar reality assumptions.
    #[staticmethod]
    #[pyo3(signature = (diagram, *, dimension=None, real=None))]
    fn from_diagram(
        diagram: PyFeynmanDiagram,
        dimension: Option<ConvertibleToExpression>,
        real: Option<Vec<PythonExpression>>,
    ) -> PyResult<Self> {
        PyAmplitude::from_diagram(diagram, dimension, real)?.squared()
    }
    /// The original coherent amplitude and its source diagrams.
    #[getter]
    fn amplitude(&self) -> PyAmplitude {
        PyAmplitude {
            inner: self.inner.amplitude().clone(),
        }
    }
    /// External labels whose spin states have been summed.
    #[getter]
    fn spin_summed(&self) -> Vec<usize> {
        self.inner.spin_summed().iter().copied().collect()
    }
    /// External labels whose color states have been summed.
    #[getter]
    fn color_summed(&self) -> Vec<usize> {
        self.inner.color_summed().iter().copied().collect()
    }
    /// Return the current tensor expression for further Spenso simplification.
    ///
    /// Examples
    /// --------
    /// >>> tensor = squared.expression().simplify_gamma().simplify_color()
    fn expression(&self, py: Python<'_>) -> PyResult<Py<TensorExpression>> {
        TensorExpression::from_atom_interface(py, self.inner.expression().clone(), None)
    }
    /// Sum selected physical spin states, optionally averaging incoming states.
    ///
    /// Omitted legs selects all unsummed legs. Omitting a massless vector reference
    /// uses the covariant sum and assumes a gauge-invariant amplitude.
    ///
    /// Examples
    /// --------
    /// >>> unpolarized = squared.sum_spins(average_initial=True)
    /// >>> photons = squared.sum_spins([2, 3], references={2: P(3), 3: P(2)})
    ///
    /// Parameters
    /// ----------
    /// legs : list[int] or None, optional
    ///     External labels to sum; default all remaining labels.
    /// average_initial : bool, optional
    ///     Divide each selected incoming completeness tensor by its state count.
    /// references : dict[int, Expression] or None, optional
    ///     Unindexed axial reference momentum per selected vector leg.
    /// spin_vectors : dict[int, Expression] or None, optional
    ///     Physical spin vector per selected massive Dirac leg.
    #[pyo3(signature = (legs=None, *, average_initial=false, references=None, spin_vectors=None))]
    fn sum_spins(
        &self,
        legs: Option<Vec<usize>>,
        average_initial: bool,
        references: Option<BTreeMap<usize, PythonExpression>>,
        spin_vectors: Option<BTreeMap<usize, PythonExpression>>,
    ) -> PyResult<Self> {
        let legs = legs.unwrap_or_else(|| {
            self.inner
                .amplitude()
                .legs()
                .iter()
                .map(|l| l.index)
                .filter(|i| !self.inner.spin_summed().contains(i))
                .collect()
        });
        let convert = |values: Option<BTreeMap<usize, PythonExpression>>| -> BTreeMap<usize, Atom> {
            values
                .unwrap_or_default()
                .into_iter()
                .map(|(i, e)| (i, e.expr))
                .collect()
        };
        Ok(Self {
            inner: self
                .inner
                .sum_spins(
                    &legs,
                    average_initial,
                    &convert(references),
                    &convert(spin_vectors),
                )
                .map_err(error::amplitude)?,
        })
    }
    /// Sum selected color states, optionally averaging incoming states.
    ///
    /// Examples
    /// --------
    /// >>> color_averaged = squared.sum_colors(average_initial=True)
    ///
    /// Parameters
    /// ----------
    /// legs : list[int] or None, optional
    ///     External labels to sum; default all remaining labels.
    /// average_initial : bool, optional
    ///     Divide by the selected incoming color-space dimensions.
    #[pyo3(signature = (legs=None, *, average_initial=false))]
    fn sum_colors(&self, legs: Option<Vec<usize>>, average_initial: bool) -> PyResult<Self> {
        let legs = legs.unwrap_or_else(|| {
            self.inner
                .amplitude()
                .legs()
                .iter()
                .map(|l| l.index)
                .filter(|i| !self.inner.color_summed().contains(i))
                .collect()
        });
        Ok(Self {
            inner: self
                .inner
                .sum_colors(&legs, average_initial)
                .map_err(error::amplitude)?,
        })
    }
    /// Describe which external states have already been summed.
    ///
    /// Examples
    /// --------
    /// >>> print(squared)
    fn __repr__(&self) -> String {
        format!(
            "SquaredAmplitude(diagrams={}, spin_summed={:?}, color_summed={:?})",
            self.inner.amplitude().diagrams().len(),
            self.spin_summed(),
            self.color_summed()
        )
    }
    /// Render the current tensor expression with Spenso's printer.
    ///
    /// Examples
    /// --------
    /// >>> squared  # notebook output
    fn _repr_html_(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        Ok(self
            .expression(py)?
            .bind(py)
            .call_method0("_repr_html_")?
            .unbind())
    }
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PyAmplitude>()?;
    module.add_class::<PySquaredAmplitude>()?;
    module.add_class::<PyAmplitudeLeg>()?;
    Ok(())
}
