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
/// >>> from symbolica import S, E
/// >>> from symbolica.community import hepkit as hep
/// >>> model = hep.Model.phi4()
/// >>> process = model.process(["phi", "phi"], ["phi", "phi"])
/// >>> generated = process.generate_diagrams()
/// >>> diagram = generated.diagrams[0]
/// >>> amplitude = hep.Amplitude(generated.diagrams)
/// >>> leg = amplitude.legs[0]
/// >>> state = (leg.index, leg.particle.name, leg.state)
/// >>> spin_sum = leg.particle.spin_sum(leg.momentum, leg.tensor_index, S("conjugate_index"))
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    name = "AmplitudeLeg",
    module = "symbolica.community.hepkit",
    frozen,
    from_py_object
)]
#[derive(Clone)]
pub struct PyAmplitudeLeg {
    inner: AmplitudeLeg,
    model: Arc<feynkit_model::Model>,
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[spenso_macros::track_usage(crate::record_usage)]
#[pymethods]
impl PyAmplitudeLeg {
    /// Stable external-leg label shared by every diagram.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``AmplitudeLeg`` class example:
    ///
    /// >>> label = leg.index
    /// >>> assert label in [state.index for state in amplitude.legs]
    #[getter]
    fn index(&self) -> usize {
        self.inner.index
    }
    /// Model particle, including particle/antiparticle identity.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``AmplitudeLeg`` class example:
    ///
    /// >>> particle_name = leg.particle.name
    #[getter]
    fn particle(&self) -> PyParticle {
        PyParticle::new(self.inner.particle, self.model.clone())
    }
    /// Whether the state is incoming or outgoing.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``AmplitudeLeg`` class example:
    ///
    /// >>> incoming = [state for state in amplitude.legs if state.state == "incoming"]
    #[getter]
    fn state(&self) -> &'static str {
        match self.inner.state {
            ExternalState::Incoming => "incoming",
            ExternalState::Outgoing => "outgoing",
        }
    }
    /// Unindexed physical momentum P(index).
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``AmplitudeLeg`` class example:
    ///
    /// >>> external_momentum = leg.momentum
    #[getter]
    fn momentum(&self) -> PythonExpression {
        PythonExpression {
            expr: self.inner.momentum(),
        }
    }
    /// Bare label shared by this leg's spin and color slots.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``AmplitudeLeg`` class example:
    ///
    /// >>> spin_sum = leg.particle.spin_sum(leg.momentum, leg.tensor_index, S("bra"))
    #[getter]
    fn tensor_index(&self) -> PythonExpression {
        PythonExpression {
            expr: self.inner.index_atom(),
        }
    }
    /// Typed open spin/color slots attached to this physical state.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``AmplitudeLeg`` class example:
    ///
    /// >>> open_slots = leg.slots
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
    /// Using the setup in the ``AmplitudeLeg`` class example:
    ///
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
/// >>> from symbolica import S, E
/// >>> from symbolica.community import hepkit as hep
/// >>> model = hep.Model.phi4()
/// >>> process = model.process(["phi", "phi"], ["phi", "phi"])
/// >>> generated = process.generate_diagrams()
/// >>> diagram = generated.diagrams[0]
/// >>> amplitude = hep.Amplitude(generated.diagrams)
/// >>> operator = amplitude.expression()
/// >>> adjoint = amplitude.conjugate().expression()
/// >>> squared = amplitude.squared().sum_spins(average_initial=True).sum_colors()
/// >>> scalar = squared.expression()
/// >>> scalar = scalar.contract(collect_chains=False, collect_traces=False)
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    name = "Amplitude",
    module = "symbolica.community.hepkit",
    frozen,
    from_py_object
)]
#[derive(Clone)]
pub struct PyAmplitude {
    pub(crate) inner: Amplitude,
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[spenso_macros::track_usage(crate::record_usage)]
#[pymethods]
impl PyAmplitude {
    /// Sum weighted diagram operators with matching graph external indices.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Amplitude`` class example:
    ///
    /// >>> amplitude = hep.Amplitude(generated.diagrams)
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
    /// Using the setup in the ``Amplitude`` class example:
    ///
    /// >>> amplitude = hep.Amplitude.from_diagram(diagram)
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
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Amplitude`` class example:
    ///
    /// >>> assert len(amplitude.diagrams) == len(generated.diagrams)
    #[getter]
    fn diagrams(&self) -> Vec<PyFeynmanDiagram> {
        self.inner
            .diagrams()
            .iter()
            .map(|d| (**d).clone().into())
            .collect()
    }
    /// Physical external states in increasing external-label order.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Amplitude`` class example:
    ///
    /// >>> external_states = [(leg.index, leg.particle.name) for leg in amplitude.legs]
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
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Amplitude`` class example:
    ///
    /// >>> operators = amplitude.terms
    /// >>> assert len(operators) == len(amplitude.diagrams)
    #[getter]
    fn terms(&self, py: Python<'_>) -> PyResult<Vec<Py<TensorExpression>>> {
        self.inner
            .terms()
            .iter()
            .map(|t| TensorExpression::from_atom_interface(py, t.clone(), None))
            .collect()
    }
    /// Whether this amplitude is the physical adjoint of its source diagrams.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Amplitude`` class example:
    ///
    /// >>> assert not amplitude.is_conjugated
    /// >>> assert amplitude.conjugate().is_conjugated
    #[getter]
    fn is_conjugated(&self) -> bool {
        self.inner.is_conjugated()
    }

    /// Return the complete operator as a Spenso tensor expression.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Amplitude`` class example:
    ///
    /// >>> operator = amplitude.expression()
    fn expression(&self, py: Python<'_>) -> PyResult<Py<TensorExpression>> {
        TensorExpression::from_atom_interface(py, self.inner.expression(), self.inner.structure())
    }
    /// Conjugate scalar coefficients, color tensors, and Dirac chains.
    ///
    /// Physical external-leg labels are preserved across different fermion pairings.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Amplitude`` class example:
    ///
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
    /// Using the setup in the ``Amplitude`` class example:
    ///
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
    /// Using the setup in the ``Amplitude`` class example:
    ///
    /// >>> print(amplitude)
    fn __repr__(&self) -> String {
        format!(
            "Amplitude(diagrams={}, external_legs={}, conjugated={})",
            self.inner.diagrams().len(),
            self.inner.legs().len(),
            self.inner.is_conjugated()
        )
    }
    /// Render a configurable snapshot of the amplitude's diagrams and weighted terms.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Amplitude`` class example:
    ///
    /// >>> drawing = amplitude.render(config=hep.RenderSettings(node_radius=5), max_diagrams=2)
    /// >>> html = drawing.to_html()
    ///
    /// Parameters
    /// ----------
    /// config : RenderSettings, optional
    ///     Layout, labels, and stroke settings shared by all diagrams.
    /// max_diagrams : int or None, optional
    ///     Maximum displayed contributions (default 6); None includes all, 0 none.
    /// term_settings : DisplaySettings, optional
    ///     Tensor notation settings for each weighted contribution.
    ///
    #[pyo3(signature = (*, config=None, max_diagrams=Some(6), term_settings=None),
        text_signature = "($self, *, config=None, max_diagrams=6, term_settings=None)")]
    fn render(
        &self,
        py: Python<'_>,
        config: Option<&crate::PyRenderSettings>,
        max_diagrams: Option<usize>,
        term_settings: Option<&spynso3::display::DisplaySettings>,
    ) -> PyResult<PyAmplitudeRender> {
        let limit = max_diagrams.unwrap_or(self.inner.diagrams().len());
        let kwargs = pyo3::types::PyDict::new(py);
        if let Some(settings) = term_settings {
            kwargs.set_item("settings", settings.clone())?;
        }
        let terms = self
            .inner
            .terms()
            .iter()
            .take(limit)
            .map(|term| {
                TensorExpression::from_atom_interface(py, term.clone(), None)?
                    .bind(py)
                    .call_method("to_html", (), Some(&kwargs))?
                    .extract::<String>()
            })
            .collect::<PyResult<Vec<_>>>()?;
        let (html, diagrams) = crate::display::collection_html(
            py,
            if self.inner.is_conjugated() {
                "Conjugate amplitude"
            } else {
                "Amplitude"
            },
            &format!(
                "{} diagrams · {} external legs",
                self.inner.diagrams().len(),
                self.inner.legs().len()
            ),
            self.inner
                .diagrams()
                .iter()
                .map(|diagram| PyFeynmanDiagram::from(diagram.as_ref().clone())),
            Some(&terms),
            config,
            limit,
        )?;
        Ok(PyAmplitudeRender { html, diagrams })
    }
    /// Render compact diagram rows with weighted operators and expandable graphs.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Amplitude`` class example:
    ///
    /// >>> from IPython.display import display
    /// >>> display(amplitude)
    fn _repr_html_(&self, py: Python<'_>) -> PyResult<String> {
        Ok(self
            .render(py, None, Some(crate::display::PREVIEW_LIMIT), None)?
            .html)
    }
}

/// A coherent amplitude square with independent ket and bra tensor indices.
///
/// State sums return new objects; an already-summed leg raises AmplitudeError.
/// No phase-space integration, symmetry factor, flux, or state average is implicit.
///
/// Examples
/// --------
/// >>> from symbolica import S, E
/// >>> from symbolica.community import hepkit as hep
/// >>> model = hep.Model.phi4()
/// >>> process = model.process(["phi", "phi"], ["phi", "phi"])
/// >>> generated = process.generate_diagrams()
/// >>> diagram = generated.diagrams[0]
/// >>> amplitude = hep.Amplitude(generated.diagrams)
/// >>> squared = amplitude.squared()
/// >>> unpolarized = squared.sum_spins(average_initial=True).sum_colors(average_initial=True)
/// >>> tensor = unpolarized.expression()
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    name = "SquaredAmplitude",
    module = "symbolica.community.hepkit",
    frozen,
    from_py_object
)]
#[derive(Clone)]
pub struct PySquaredAmplitude {
    inner: SquaredAmplitude,
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[spenso_macros::track_usage(crate::record_usage)]
#[pymethods]
impl PySquaredAmplitude {
    /// Construct the coherent square of a single unsewn diagram.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``SquaredAmplitude`` class example:
    ///
    /// >>> squared = hep.SquaredAmplitude.from_diagram(diagram)
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
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``SquaredAmplitude`` class example:
    ///
    /// >>> source_diagrams = squared.amplitude.diagrams
    #[getter]
    fn amplitude(&self) -> PyAmplitude {
        PyAmplitude {
            inner: self.inner.amplitude().clone(),
        }
    }
    /// External labels whose spin states have been summed.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``SquaredAmplitude`` class example:
    ///
    /// >>> summed = squared.sum_spins()
    /// >>> assert set(summed.spin_summed) == {leg.index for leg in amplitude.legs}
    #[getter]
    fn spin_summed(&self) -> Vec<usize> {
        self.inner.spin_summed().iter().copied().collect()
    }
    /// External labels whose color states have been summed.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``SquaredAmplitude`` class example:
    ///
    /// >>> summed = squared.sum_colors()
    /// >>> assert set(summed.color_summed) == {leg.index for leg in amplitude.legs}
    #[getter]
    fn color_summed(&self) -> Vec<usize> {
        self.inner.color_summed().iter().copied().collect()
    }
    /// Return the current tensor expression for further Spenso simplification.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``SquaredAmplitude`` class example:
    ///
    /// >>> tensor = squared.expression()
    /// >>> tensor = tensor.contract(collect_chains=False, collect_traces=False)
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
    /// Using the setup in the ``SquaredAmplitude`` class example:
    ///
    /// >>> unpolarized = squared.sum_spins(average_initial=True)
    /// >>> assert set(unpolarized.spin_summed) == {leg.index for leg in amplitude.legs}
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
    /// Using the setup in the ``SquaredAmplitude`` class example:
    ///
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
    /// Using the setup in the ``SquaredAmplitude`` class example:
    ///
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
    /// Using the setup in the ``SquaredAmplitude`` class example:
    ///
    /// >>> from IPython.display import display
    /// >>> display(squared)
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
    module.add_class::<PyAmplitudeRender>()?;
    module.add_class::<PySquaredAmplitude>()?;
    module.add_class::<PyAmplitudeLeg>()?;
    Ok(())
}

/// A rendered amplitude collection retaining its configured diagram snapshots.
///
/// Examples
/// --------
/// >>> from symbolica.community import hepkit as hep
/// >>> process = hep.Model.phi4().process(["phi", "phi"], ["phi", "phi"])
/// >>> amplitude = hep.Amplitude(process.generate_diagrams().diagrams)
/// >>> drawing = amplitude.render(max_diagrams=None)
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    name = "AmplitudeRender",
    module = "symbolica.community.hepkit",
    frozen
)]
pub struct PyAmplitudeRender {
    html: String,
    /// Rendered diagram snapshots in contribution order, bounded by max_diagrams.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``AmplitudeRender`` class example:
    ///
    /// >>> svg = drawing.diagrams[0].to_svg()
    #[pyo3(get)]
    diagrams: Vec<crate::PyDiagramRender>,
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[pymethods]
impl PyAmplitudeRender {
    /// Export the configured collection as interactive notebook HTML.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``AmplitudeRender`` class example:
    ///
    /// >>> html = drawing.to_html()
    fn to_html(&self) -> &str {
        &self.html
    }

    /// Display the configured collection in IPython and Jupyter.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``AmplitudeRender`` class example:
    ///
    /// >>> from IPython.display import display
    /// >>> display(drawing)
    fn _repr_html_(&self) -> &str {
        self.to_html()
    }

    /// Display the configured collection in Marimo.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``AmplitudeRender`` class example:
    ///
    /// >>> import marimo as mo
    /// >>> mo.as_html(drawing)
    fn _mime_(&self) -> (&str, &str) {
        ("text/html", self.to_html())
    }

    /// Summarize the number of rendered contributions in text-only frontends.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``AmplitudeRender`` class example:
    ///
    /// >>> text = repr(drawing)
    fn __repr__(&self) -> String {
        format!("AmplitudeRender(diagrams={})", self.diagrams.len())
    }
}
