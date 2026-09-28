use std::{collections::BTreeMap, path::PathBuf, sync::Arc};

use feynkit_model::{
    ComplexValue, Coupling, CouplingId, EvaluatedValues, EvaluationRequest, LorentzStructure,
    LorentzStructureId, Model, ModelEvaluator, ModelExpression, ModelFormFactor, ModelFormFactorId,
    ModelFunction, ModelFunctionId, Parameter, ParameterCard, ParameterId, ParameterNature,
    ParameterType, Particle, ParticleId, Propagator, PropagatorId, VertexRule, VertexRuleId,
};
use pyo3::{
    exceptions::PyValueError,
    prelude::*,
    types::{PyAny, PyComplex, PyModule},
};

#[cfg(feature = "python_stubgen")]
use pyo3_stub_gen::{
    derive::{
        gen_methods_from_python, gen_stub_pyclass, gen_stub_pyclass_enum, gen_stub_pymethods,
    },
    inventory::submit,
};

use crate::{
    display::{escape_html, model_expression_html, model_record_html},
    error,
    generation::{PyProcess, SelectorInput, VertexInput},
};
use spynso3::{display::DisplaySettings, expression::TensorExpression};
use symbolica::{
    api::python::{ConvertibleToExpression, PythonExpression},
    atom::{Atom, AtomCore},
    symbol,
};

fn complex_value<'py>(py: Python<'py>, value: ComplexValue) -> Bound<'py, PyComplex> {
    PyComplex::from_doubles(py, value.re, value.im)
}

fn display_value(value: Option<ComplexValue>) -> String {
    match value {
        Some(v) if v.im == 0.0 => v.re.to_string(),
        Some(v) => format!("{} {:+}i", v.re, v.im),
        None => "not evaluated".to_owned(),
    }
}

/// A particle species in a loaded interaction model.
///
/// Particle records expose the signed PDG code, spin and color
/// representations, electric charge, and the parameters used for mass and
/// width.
///
/// Examples
/// --------
/// >>> from symbolica import S, E
/// >>> from symbolica.community import hep
/// >>> model = hep.Model.standard_model()
/// >>> electron = model.particle_by_pdg(11)
/// >>> assert electron.name == "e-"
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    name = "Particle",
    module = "symbolica.community.feynkit",
    frozen,
    from_py_object
)]
#[derive(Clone)]
pub struct PyParticle {
    id: ParticleId,
    model: Arc<Model>,
}

impl PyParticle {
    pub(crate) fn new(id: ParticleId, model: Arc<Model>) -> Self {
        Self { id, model }
    }

    fn inner(&self) -> &Particle {
        self.model.particle_by_id(self.id).unwrap()
    }

    pub(crate) fn signed_pdg(&self) -> i64 {
        self.inner().pdg_code
    }
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[pymethods]
impl PyParticle {
    /// Construct the identity on this particle's color space.
    ///
    /// The bare ``left`` index carries the particle representation; ``right``
    /// carries its dual. Antiquarks and antisextets reverse the dual orientation.
    /// Supports UFO singlet, fundamental, sextet and adjoint representations.
    /// ``average=True`` divides by the number of color states (1, 3, 6 or 8).
    /// This sums color only; spin sums and color-algebra simplification remain
    /// separate operations. To close an existing tensor, use indices whose
    /// slots are dual to the open color slots of that tensor.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Particle`` class example:
    ///
    /// >>> from symbolica import S
    /// >>> i, j = S("i", "j")
    /// >>> color_projector = model.particle_by_pdg(5).color_sum(i, j, average=True)
    /// >>> model.particle_by_pdg(11).color_sum(i, j) == 1
    /// True
    ///
    /// Parameters
    /// ----------
    /// left : Expression
    ///     Bare index in the particle's color representation.
    /// right : Expression
    ///     Bare index in its dual representation.
    /// average : bool, optional
    ///     Divide by the number of color states. Defaults to ``False``.
    ///
    /// Returns
    /// -------
    /// Expression
    ///     Spenso color identity, or one for a color singlet.
    ///
    /// Raises
    /// ------
    /// ValueError
    ///     If the particle's UFO color representation is unsupported.
    #[pyo3(signature = (left, right, *, average=false))]
    fn color_sum(
        &self,
        left: &PythonExpression,
        right: &PythonExpression,
        average: bool,
    ) -> PyResult<PythonExpression> {
        let sum = feynkit_amplitude::ColorSum::new(self.inner())
            .map_err(|error| PyValueError::new_err(error.to_string()))?
            .averaged(average);
        Ok(PythonExpression {
            expr: sum.expression([left.expr.clone(), right.expr.clone()]),
        })
    }

    /// Construct this particle's external-state spin or polarization sum.
    ///
    /// Return an ordinary Symbolica expression using Spenso gamma matrices
    /// and metrics. Indices are bare symbols; momentum and reference are
    /// unindexed symbols or labeled calls such as ``Q(1)``. The calculation
    /// defaults to four-dimensional external states. ``dimension`` changes the
    /// Lorentz dimension while Dirac spinor slots retain dimension four.
    /// Massive vectors use the Proca
    /// projector. For massless vectors, supply a reference for a physical
    /// axial sum, or omit it for the covariant sum of a gauge-invariant
    /// amplitude. Subsequent kinematic substitutions must enforce on-shell
    /// conditions and a nonzero momentum-reference scalar product.
    ///
    /// For a massive Dirac particle, ``spin_vector`` selects one physical spin
    /// state instead of summing states. Supply a dimensionless unindexed vector
    /// satisfying ``p.s = 0`` and ``s.s = -1``: the rest-frame spin direction,
    /// boosted with the particle. The same projector sign applies to fermions
    /// and antifermions; do not reverse this vector for an antiparticle.
    /// This option requires ``average=False`` and four Lorentz dimensions.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Particle`` class example:
    ///
    /// >>> from symbolica import S
    /// >>> p, i, j = S("p", "i", "j")
    /// >>> projector = model.particle("e-").spin_sum(p, i, j, average=True)
    /// >>> polarized = model.particle("ta-").spin_sum(p, i, j, spin_vector=S("s"))
    /// >>> dimensional = model.particle("g").spin_sum(p, i, j, dimension=S("D"))
    /// >>> six_dimensional = model.particle("g").spin_sum(p, i, j, dimension=6)
    ///
    /// Parameters
    /// ----------
    /// momentum : Expression
    ///     Unindexed external momentum.
    /// left : Expression
    ///     Open index on the amplitude.
    /// right : Expression
    ///     Open index on the conjugate amplitude.
    /// average : bool
    ///     Divide by two for Dirac fermions, D-2 for massless vectors, or D-1
    ///     for massive vectors. Scalars have one state.
    /// dimension : Expression | int | None
    ///     Integer or symbolic Lorentz dimension; defaults to four. Dirac
    ///     spinor slots and their trace dimension stay four. To use a fixed
    ///     two-state vector average at symbolic D, leave average=False and
    ///     divide the result by two.
    /// reference : Expression | None
    ///     Axial reference momentum for a massless vector; need not be null.
    /// covariant : bool
    ///     Use the Feynman-gauge vector numerator even for a massive vector.
    /// spin_vector : Expression | None
    ///     Physical spin vector of a massive Dirac state, with ``p.s = 0`` and
    ///     ``s.s = -1``. Requires ``average=False``; ``None`` sums both states.
    ///
    /// Raises
    /// ------
    /// ValueError
    ///     If ``spin_vector`` is used with averaging, a massless particle, or
    ///     a particle other than a Dirac fermion, non-four-dimensional Lorentz
    ///     slots, or is not an unindexed name. Also raised for an invalid
    ///     dimension or a concrete dimension with no physical vector states.
    #[pyo3(signature = (momentum, left, right, *, average=false, reference=None, covariant=false, spin_vector=None, dimension=None))]
    #[allow(clippy::too_many_arguments)]
    fn spin_sum(
        &self,
        momentum: &PythonExpression,
        left: &PythonExpression,
        right: &PythonExpression,
        average: bool,
        reference: Option<&PythonExpression>,
        covariant: bool,
        spin_vector: Option<&PythonExpression>,
        dimension: Option<ConvertibleToExpression>,
    ) -> PyResult<PythonExpression> {
        let sum = feynkit_amplitude::SpinSum::new(self.inner(), &self.model)
            .map_err(|error| PyValueError::new_err(error.to_string()))?
            .with_dimension(
                &dimension
                    .map_or_else(|| symbolica::atom::Atom::num(4), |d| d.to_expression().expr),
            )
            .map_err(|error| PyValueError::new_err(error.to_string()))?
            .averaged(average)
            .covariant(covariant);
        sum.expression(
            &momentum.expr,
            [left.expr.clone(), right.expr.clone()],
            reference.map(|reference| &reference.expr),
            spin_vector.map(|spin_vector| &spin_vector.expr),
        )
        .map(|expr| PythonExpression { expr })
        .map_err(|error| PyValueError::new_err(error.to_string()))
    }

    /// Sum paired generated external wavefunctions for one edge.
    ///
    /// Replace this particle's wavefunction and its adjoint using the same
    /// completeness relation as ``spin_sum``. The expression may be a sewn
    /// diagram's ``projector_expression()`` or a squared amplitude. Only pairs
    /// with the supplied edge label are replaced; unpaired wavefunctions stay
    /// unchanged. Scalar particles have no external wavefunction factors.
    /// Vector wavefunction slots must use the requested Lorentz dimension.
    /// External states default to four dimensions; reference and gauge conventions
    /// are those of ``spin_sum``. For a massive Dirac particle, ``spin_vector``
    /// selects the same physical spin state as in ``spin_sum``, including for
    /// antiparticles. This does not sum color, conjugate amplitudes, or apply
    /// graph symmetry factors.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Particle`` class example:
    ///
    /// >>> from symbolica import S, E
    /// >>> from symbolica.community import hep
    /// >>> model = hep.Model.standard_model()
    /// >>> diagram = model.process(["e-", "e+"], ["mu-", "mu+"]).generate_diagrams().diagrams[0]
    /// >>> electron = model.particle("e-")
    /// >>> projector = diagram.projector_expression()
    /// >>> spin_summed = electron.sum_spins(projector, S("gammalooprs::P")(1), edge=1)
    ///
    /// Parameters
    /// ----------
    /// expression : Expression
    ///     Projector or squared expression containing paired wavefunctions.
    /// momentum : Expression
    ///     Unindexed external momentum in the physical particle direction.
    /// edge : int
    ///     Generated edge label of the pair to replace.
    /// average : bool
    ///     Divide by two for Dirac fermions, D-2 for massless vectors, or D-1
    ///     for massive vectors. Scalars have one state.
    /// dimension : Expression | int | None
    ///     Integer or symbolic Lorentz dimension; defaults to four. Dirac
    ///     spinor slots and their trace dimension stay four. To use a fixed
    ///     two-state vector average at symbolic D, leave average=False and
    ///     divide the result by two.
    /// reference : Expression | None
    ///     Axial reference for a massless vector; need not be null.
    /// covariant : bool
    ///     Use the Feynman-gauge vector numerator even for a massive vector.
    /// spin_vector : Expression | None
    ///     Physical spin vector of a massive Dirac state, with ``p.s = 0`` and
    ///     ``s.s = -1``. Requires ``average=False``; ``None`` sums both states.
    ///
    /// Raises
    /// ------
    /// ValueError
    ///     If ``spin_vector`` is used with averaging, a massless particle, or
    ///     a particle other than a Dirac fermion, non-four-dimensional Lorentz
    ///     slots, or is not an unindexed name. Also raised for an invalid
    ///     dimension or a concrete dimension with no physical vector states.
    #[pyo3(signature = (expression, momentum, *, edge, average=false, reference=None, covariant=false, spin_vector=None, dimension=None))]
    #[allow(clippy::too_many_arguments)]
    fn sum_spins(
        &self,
        expression: &PythonExpression,
        momentum: &PythonExpression,
        edge: usize,
        average: bool,
        reference: Option<&PythonExpression>,
        covariant: bool,
        spin_vector: Option<&PythonExpression>,
        dimension: Option<ConvertibleToExpression>,
    ) -> PyResult<PythonExpression> {
        feynkit_amplitude::SpinSum::new(self.inner(), &self.model)
            .map_err(|error| PyValueError::new_err(error.to_string()))?
            .with_dimension(
                &dimension
                    .map_or_else(|| symbolica::atom::Atom::num(4), |d| d.to_expression().expr),
            )
            .map_err(|error| PyValueError::new_err(error.to_string()))?
            .averaged(average)
            .covariant(covariant)
            .apply(
                &expression.expr,
                &momentum.expr,
                edge,
                reference.map(|reference| &reference.expr),
                spin_vector.map(|spin_vector| &spin_vector.expr),
            )
            .map(|expr| PythonExpression { expr })
            .map_err(|error| PyValueError::new_err(error.to_string()))
    }

    /// Return the particle name used by the model.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Particle`` class example:
    ///
    /// >>> assert electron.name == "e-"
    #[getter]
    fn name(&self) -> &str {
        &self.inner().name
    }

    /// Return the name of the corresponding antiparticle.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Particle`` class example:
    ///
    /// >>> assert electron.antiname == "e+"
    #[getter]
    fn antiname(&self) -> &str {
        &self
            .model
            .particle_by_id(self.inner().antiparticle)
            .unwrap()
            .name
    }

    /// Return the corresponding antiparticle from the same model.
    ///
    /// Self-conjugate particles map to themselves.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Particle`` class example:
    ///
    /// >>> electron = model.particle_by_pdg(11)
    /// >>> electron.antiparticle.name
    /// 'e+'
    /// >>> model.particle_by_pdg(22).antiparticle.name  # the photon is self-conjugate
    /// 'a'
    #[getter]
    fn antiparticle(&self) -> PyResult<PyParticle> {
        Ok(Self::new(
            self.inner().antiparticle,
            Arc::clone(&self.model),
        ))
    }

    /// Return the signed Particle Data Group code.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Particle`` class example:
    ///
    /// >>> model.particle("e-").pdg_code
    /// 11
    #[getter]
    fn pdg_code(&self) -> i64 {
        self.inner().pdg_code
    }

    /// Return the UFO spin code ``2S + 1`` (negative for ghost fields).
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Particle`` class example:
    ///
    /// >>> model.particle_by_pdg(11).spin  # spin-1/2 electron
    /// 2
    #[getter]
    fn spin(&self) -> i64 {
        self.inner().spin
    }

    /// Return the signed UFO SU(3) color representation code.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Particle`` class example:
    ///
    /// >>> model.particle_by_pdg(11).color  # color-singlet electron
    /// 1
    #[getter]
    fn color(&self) -> i64 {
        self.inner().color
    }

    /// Return the name of the particle's mass parameter.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Particle`` class example:
    ///
    /// >>> particle = model.particle_by_pdg(13)
    /// >>> mass = model.parameter(particle.mass_parameter)
    #[getter]
    fn mass_parameter(&self) -> &str {
        &self.model.parameter_by_id(self.inner().mass).unwrap().name
    }

    /// Return the name of the particle's width parameter.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Particle`` class example:
    ///
    /// >>> particle = model.particle_by_pdg(23)
    /// >>> width = model.parameter(particle.width_parameter)
    #[getter]
    fn width_parameter(&self) -> &str {
        &self.model.parameter_by_id(self.inner().width).unwrap().name
    }

    /// Exact symbolic mass, with the UFO ZERO parameter represented as zero.
    /// Other parameters remain symbolic, even when their current value is zero.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Particle`` class example:
    ///
    /// >>> model.particle("c").mass_expression
    #[getter]
    fn mass_expression(&self) -> PythonExpression {
        self.inner().symbolic_mass(&self.model).into()
    }

    /// Electric charge in units of e, as an exact Symbolica expression.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Particle`` class example:
    ///
    /// >>> model.particle("c").charge
    #[getter]
    fn charge(&self) -> PythonExpression {
        Atom::num(self.inner().charge.clone()).into()
    }

    /// Hypercharge in Q = T3 + Y/2; left-handed for fermions.
    /// None means absent, undefined, or unspecified, rather than zero.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Particle`` class example:
    ///
    /// >>> model.particle("c").y_charge
    #[getter]
    fn y_charge(&self) -> Option<PythonExpression> {
        self.inner().y_charge.clone().map(|y| Atom::num(y).into())
    }

    /// Right-handed fermion hypercharge in Q = T3 + Y/2, if specified.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Particle`` class example:
    ///
    /// >>> model.particle("c").y_charge_right
    #[getter]
    fn y_charge_right(&self) -> Option<PythonExpression> {
        self.inner()
            .y_charge_right
            .clone()
            .map(|y| Atom::num(y).into())
    }

    /// Third weak-isospin component Q - Y/2; left-handed for fermions.
    /// Antiparticle chiralities are exchanged by charge conjugation.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Particle`` class example:
    ///
    /// >>> model.particle("c").weak_isospin
    #[getter]
    fn weak_isospin(&self) -> Option<PythonExpression> {
        self.inner().weak_isospin().map(|t| Atom::num(t).into())
    }

    /// Right-handed fermion third weak-isospin component, if specified.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Particle`` class example:
    ///
    /// >>> model.particle("c").weak_isospin_right
    #[getter]
    fn weak_isospin_right(&self) -> Option<PythonExpression> {
        self.inner()
            .weak_isospin_right()
            .map(|t| Atom::num(t).into())
    }

    /// Report whether this object represents an antiparticle.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Particle`` class example:
    ///
    /// >>> model.particle_by_pdg(11).is_antiparticle
    /// False
    #[getter]
    fn is_antiparticle(&self) -> bool {
        self.inner().is_antiparticle()
    }

    /// Report whether the particle is its own antiparticle.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Particle`` class example:
    ///
    /// >>> model.particle_by_pdg(22).is_self_antiparticle
    /// True
    #[getter]
    fn is_self_antiparticle(&self) -> bool {
        self.model.particle_is_self_conjugate(self.id)
    }

    /// Report whether the particle has fermionic spin.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Particle`` class example:
    ///
    /// >>> model.particle_by_pdg(11).is_fermion
    /// True
    #[getter]
    fn is_fermion(&self) -> bool {
        self.inner().is_fermion()
    }

    /// Report whether the particle's mass parameter is zero.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Particle`` class example:
    ///
    /// >>> model.particle_by_pdg(22).is_massless
    /// True
    #[getter]
    fn is_massless(&self) -> bool {
        self.model.particle_is_massless(self.id)
    }

    /// Summarize the defining data as well as the name and PDG code.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Particle`` class example:
    ///
    /// >>> print(model.particle_by_pdg(11))
    fn __repr__(&self) -> String {
        format!(
            "Particle({:?}, pdg={}, antiparticle={:?}, spin={}, color={}, charge={}, mass={}, width={})",
            self.name(),
            self.pdg_code(),
            self.model
                .particle_by_id(self.inner().antiparticle)
                .unwrap()
                .name,
            self.spin(),
            self.color(),
            self.inner().charge,
            self.mass_parameter(),
            self.width_parameter()
        )
    }

    /// Display the model member's physical data and defining expressions.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Particle`` class example:
    ///
    /// >>> from IPython.display import display
    /// >>> display(electron)
    fn _repr_html_(&self, py: Python<'_>) -> PyResult<String> {
        let p = self.inner();
        let spin = if p.is_ghost() {
            "ghost".to_owned()
        } else {
            format!(
                "{} (UFO {})",
                symbolica::domains::rational::Rational::from((p.spin - 1, 2)),
                p.spin
            )
        };
        let mut rows = vec![
            ("PDG", p.pdg_code.to_string()),
            (
                "Antiparticle",
                escape_html(&self.model.particle_by_id(p.antiparticle).unwrap().name),
            ),
            ("Spin", spin),
            ("Color", p.color.to_string()),
            ("Charge", model_expression_html(py, self.charge())?),
            ("Mass", model_expression_html(py, self.mass_expression())?),
            ("Width", escape_html(self.width_parameter())),
            (
                "Ghost / lepton number",
                format!("{} / {}", p.ghost_number, p.lepton_number),
            ),
            (
                "Propagating / Goldstone",
                format!("{} / {}", p.propagating, p.goldstone),
            ),
        ];
        if let Some(y) = self.y_charge() {
            rows.push(("Hypercharge (left)", model_expression_html(py, y)?));
        }
        if let Some(y) = self.y_charge_right() {
            rows.push(("Hypercharge (right)", model_expression_html(py, y)?));
        }
        Ok(model_record_html("Particle", self.name(), &rows))
    }

    /// Write the complete text summary to an IPython pretty printer.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Particle`` class example:
    ///
    /// >>> from IPython.lib.pretty import pretty
    /// >>> text = pretty(electron)
    ///
    /// Parameters
    /// ----------
    /// pretty : object
    ///     IPython pretty printer receiving the text.
    /// cycle : bool
    ///     Whether the object occurs recursively in the current display.
    fn _repr_pretty_(&self, pretty: &Bound<'_, PyAny>, cycle: bool) -> PyResult<()> {
        pretty.call_method1(
            "text",
            (if cycle {
                "Particle(...)".to_owned()
            } else {
                self.__repr__()
            },),
        )?;
        Ok(())
    }
}

/// Whether a model parameter is supplied externally or derived internally.
///
/// Examples
/// --------
/// >>> from symbolica import S, E
/// >>> from symbolica.community import hep
/// >>> model = hep.Model.standard_model()
/// >>> nature = hep.ParameterNature.EXTERNAL
/// >>> external = [p for p in model.parameters if p.nature == nature]
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass_enum)]
#[pyclass(
    name = "ParameterNature",
    module = "symbolica.community.feynkit",
    rename_all = "SCREAMING_SNAKE_CASE",
    frozen,
    eq,
    eq_int,
    from_py_object
)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PyParameterNature {
    External,
    Internal,
}

impl From<ParameterNature> for PyParameterNature {
    fn from(value: ParameterNature) -> Self {
        match value {
            ParameterNature::External => Self::External,
            ParameterNature::Internal => Self::Internal,
        }
    }
}

/// Whether a model parameter is real-valued or complex-valued.
///
/// Examples
/// --------
/// >>> from symbolica import S, E
/// >>> from symbolica.community import hep
/// >>> model = hep.Model.standard_model()
/// >>> real_parameters = [p for p in model.parameters
/// ...                    if p.parameter_type == hep.ParameterType.REAL]
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass_enum)]
#[pyclass(
    name = "ParameterType",
    module = "symbolica.community.feynkit",
    rename_all = "SCREAMING_SNAKE_CASE",
    frozen,
    eq,
    eq_int,
    from_py_object
)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PyParameterType {
    Real,
    Complex,
}

impl From<ParameterType> for PyParameterType {
    fn from(value: ParameterType) -> Self {
        match value {
            ParameterType::Real => Self::Real,
            ParameterType::Complex => Self::Complex,
        }
    }
}

/// A numerical or symbolic parameter in a particle-physics model.
///
/// External parameters carry parameter-card coordinates and values; internal
/// parameters carry expressions derived from the external inputs.
///
/// Examples
/// --------
/// >>> from symbolica import S, E
/// >>> from symbolica.community import hep
/// >>> model = hep.Model.standard_model()
/// >>> mass = model.parameter("MM")
/// >>> assert mass.nature == hep.ParameterNature.EXTERNAL
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    name = "Parameter",
    module = "symbolica.community.feynkit",
    frozen,
    from_py_object
)]
#[derive(Clone)]
pub struct PyParameter {
    id: ParameterId,
    model: Arc<Model>,
}

impl PyParameter {
    fn new(id: ParameterId, model: Arc<Model>) -> Self {
        Self { id, model }
    }

    fn inner(&self) -> &Parameter {
        self.model.parameter_by_id(self.id).unwrap()
    }
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[pymethods]
impl PyParameter {
    /// Return the symbolic reference used by this model's expressions.
    /// This preserves the parameter or function identity without substituting
    /// its defining expression or numerical value.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import hep
    /// >>> model = hep.Model.standard_model()
    /// >>> reference = model.parameter("ee").symbol
    #[getter]
    fn symbol(&self) -> PythonExpression {
        Atom::var(symbol!(&format!("UFO::{}", self.inner().name))).into()
    }

    /// Return the parameter name.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Parameter`` class example:
    ///
    /// >>> assert mass.name == "MM"
    #[getter]
    fn name(&self) -> &str {
        &self.inner().name
    }

    /// The model's LaTeX display label, or None when no label was supplied.
    /// MiTeX renders this label in Typst and notebook math output.
    ///
    /// >>> hep.Model.standard_model().parameter("ee").texname
    /// 'e'
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Parameter`` class example:
    ///
    /// >>> label = model.parameter("ee").texname
    /// >>> assert label == "e"
    #[getter]
    fn texname(&self) -> Option<&str> {
        self.inner().texname.as_deref()
    }

    /// Return the Les Houches block name, when defined.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Parameter`` class example:
    ///
    /// >>> mass = model.parameter(model.particle_by_pdg(13).mass_parameter)
    /// >>> mass.lhablock
    /// 'MASS'
    #[getter]
    fn lhablock(&self) -> Option<String> {
        self.inner().lhablock.clone()
    }

    /// Return the Les Houches entry indices, when defined.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Parameter`` class example:
    ///
    /// >>> mass = model.parameter(model.particle_by_pdg(13).mass_parameter)
    /// >>> mass.lhacode
    /// [13]
    #[getter]
    fn lhacode(&self) -> Option<Vec<usize>> {
        self.inner().lhacode.clone()
    }

    /// Return whether the parameter is external or internal.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Parameter`` class example:
    ///
    /// >>> mass = model.parameter(model.particle_by_pdg(13).mass_parameter)
    /// >>> mass.nature == hep.ParameterNature.EXTERNAL
    /// True
    #[getter]
    fn nature(&self) -> PyParameterNature {
        self.inner().nature.clone().into()
    }

    /// Return whether the parameter is real or complex.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Parameter`` class example:
    ///
    /// >>> mass = model.parameter(model.particle_by_pdg(13).mass_parameter)
    /// >>> mass.parameter_type == hep.ParameterType.REAL
    /// True
    #[getter]
    fn parameter_type(&self) -> PyParameterType {
        self.inner().parameter_type.clone().into()
    }

    /// Return the evaluated value as a native Python complex number.
    ///
    /// Raises :class:`ModelError` when this parameter has not been evaluated.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Parameter`` class example:
    ///
    /// >>> mass_value = model.parameter("MM").value
    /// >>> print("mass:", mass_value.real, "width component:", mass_value.imag)
    #[getter]
    fn value<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyComplex>> {
        self.inner()
            .value
            .map(|value| complex_value(py, value))
            .ok_or_else(|| {
                error::ModelError::new_err(format!(
                    "parameter {:?} has not been evaluated",
                    self.inner().name
                ))
            })
    }

    /// Return the defining expression for an internal parameter.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Parameter`` class example:
    ///
    /// >>> internal = next(p for p in model.parameters if p.expression is not None)
    /// >>> formula = internal.expression
    #[getter]
    fn expression(&self) -> Option<PythonExpression> {
        self.inner()
            .expression
            .clone()
            .map(|expr| PythonExpression { expr })
    }

    /// Summarize the defining data as well as the parameter name.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Parameter`` class example:
    ///
    /// >>> print(model.parameter(model.particle_by_pdg(13).mass_parameter))
    fn __repr__(&self) -> String {
        format!(
            "Parameter({:?}, {:?}, {:?}, expression={}, value={})",
            self.name(),
            self.inner().nature,
            self.inner().parameter_type,
            self.inner().expression.as_ref().map_or_else(
                || "external input".to_owned(),
                AtomCore::to_canonical_string
            ),
            display_value(self.inner().value)
        )
    }

    /// Display the model member's physical data and defining expressions.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Parameter`` class example:
    ///
    /// >>> from IPython.display import display
    /// >>> display(mass)
    fn _repr_html_(&self, py: Python<'_>) -> PyResult<String> {
        let p = self.inner();
        let mut rows = vec![
            (
                "Symbol",
                model_expression_html(
                    py,
                    PythonExpression {
                        expr: Atom::var(symbol!(&format!("UFO::{}", p.name))),
                    },
                )?,
            ),
            (
                "Nature / type",
                format!("{:?} / {:?}", p.nature, p.parameter_type),
            ),
            ("Value", escape_html(&display_value(p.value))),
        ];
        if let Some(expression) = self.expression() {
            rows.push(("Definition", model_expression_html(py, expression)?));
        }
        if let Some(block) = &p.lhablock {
            rows.push((
                "LHA entry",
                escape_html(&format!(
                    "{} {:?}",
                    block,
                    p.lhacode.as_deref().unwrap_or_default()
                )),
            ));
        }
        Ok(model_record_html("Parameter", self.name(), &rows))
    }

    /// Write the complete text summary to an IPython pretty printer.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Parameter`` class example:
    ///
    /// >>> from IPython.lib.pretty import pretty
    /// >>> text = pretty(mass)
    ///
    /// Parameters
    /// ----------
    /// pretty : object
    ///     IPython pretty printer receiving the text.
    /// cycle : bool
    ///     Whether the object occurs recursively in the current display.
    fn _repr_pretty_(&self, pretty: &Bound<'_, PyAny>, cycle: bool) -> PyResult<()> {
        pretty.call_method1(
            "text",
            (if cycle {
                "Parameter(...)".to_owned()
            } else {
                self.__repr__()
            },),
        )?;
        Ok(())
    }
}

/// A coupling coefficient associated with interaction vertices.
///
/// Couplings retain their symbolic expression, perturbative coupling orders,
/// and evaluated complex value when the model has been numerically resolved.
///
/// Examples
/// --------
/// >>> from symbolica import S, E
/// >>> from symbolica.community import hep
/// >>> model = hep.Model.standard_model()
/// >>> coupling = model.couplings[0]
/// >>> orders = coupling.orders
/// >>> formula = coupling.expression
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    name = "Coupling",
    module = "symbolica.community.feynkit",
    frozen,
    from_py_object
)]
#[derive(Clone)]
pub struct PyCoupling {
    id: CouplingId,
    model: Arc<Model>,
}

impl PyCoupling {
    fn new(id: CouplingId, model: Arc<Model>) -> Self {
        Self { id, model }
    }

    fn inner(&self) -> &Coupling {
        self.model.coupling_by_id(self.id).unwrap()
    }
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[pymethods]
impl PyCoupling {
    /// Return the symbolic reference used by this model's expressions.
    /// This preserves the parameter or function identity without substituting
    /// its defining expression or numerical value.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import hep
    /// >>> model = hep.Model.standard_model()
    /// >>> reference = model.couplings[0].symbol
    #[getter]
    fn symbol(&self) -> PythonExpression {
        Atom::var(symbol!(&format!("UFO::{}", self.inner().name))).into()
    }

    /// Return the coupling name.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Coupling`` class example:
    ///
    /// >>> coupling_by_name = {item.name: item.expression for item in model.couplings}
    #[getter]
    fn name(&self) -> &str {
        &self.inner().name
    }

    /// Return the expression defining the coupling.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Coupling`` class example:
    ///
    /// >>> formula = coupling.expression
    #[getter]
    fn expression(&self) -> PythonExpression {
        PythonExpression {
            expr: self.inner().expression.clone(),
        }
    }

    /// Return the coupling-order powers keyed by order name.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Coupling`` class example:
    ///
    /// >>> qed_couplings = [c for c in model.couplings if c.orders.get("QED", 0) > 0]
    #[getter]
    fn orders(&self) -> BTreeMap<String, usize> {
        self.inner().orders.clone()
    }

    /// Return the evaluated coupling as a native Python complex number.
    ///
    /// Raises :class:`ModelError` when this coupling has not been evaluated.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Coupling`` class example:
    ///
    /// >>> coupling_value = model.couplings[0].value
    /// >>> print("coupling:", coupling_value.real, coupling_value.imag)
    #[getter]
    fn value<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyComplex>> {
        self.inner()
            .value
            .map(|value| complex_value(py, value))
            .ok_or_else(|| {
                error::ModelError::new_err(format!(
                    "coupling {:?} has not been evaluated",
                    self.inner().name
                ))
            })
    }

    /// Summarize the defining data as well as the coupling name.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Coupling`` class example:
    ///
    /// >>> print(next(c for c in model.couplings if c.orders.get("QED", 0) > 0))
    fn __repr__(&self) -> String {
        format!(
            "Coupling({:?}, expression={}, orders={:?}, value={})",
            self.name(),
            self.inner().expression.to_canonical_string(),
            self.inner().orders,
            display_value(self.inner().value)
        )
    }

    /// Display the model member's physical data and defining expressions.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Coupling`` class example:
    ///
    /// >>> from IPython.display import display
    /// >>> display(coupling)
    fn _repr_html_(&self, py: Python<'_>) -> PyResult<String> {
        Ok(model_record_html(
            "Coupling",
            self.name(),
            &[
                ("Definition", model_expression_html(py, self.expression())?),
                ("Orders", escape_html(&format!("{:?}", self.inner().orders))),
                ("Value", escape_html(&display_value(self.inner().value))),
            ],
        ))
    }

    /// Write the complete text summary to an IPython pretty printer.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Coupling`` class example:
    ///
    /// >>> from IPython.lib.pretty import pretty
    /// >>> text = pretty(coupling)
    ///
    /// Parameters
    /// ----------
    /// pretty : object
    ///     IPython pretty printer receiving the text.
    /// cycle : bool
    ///     Whether the object occurs recursively in the current display.
    fn _repr_pretty_(&self, pretty: &Bound<'_, PyAny>, cycle: bool) -> PyResult<()> {
        pretty.call_method1(
            "text",
            (if cycle {
                "Coupling(...)".to_owned()
            } else {
                self.__repr__()
            },),
        )?;
        Ok(())
    }
}

/// An interaction vertex rule from a particle model.
///
/// A vertex rule links its participating particles to color tensors, Lorentz
/// structures, and the corresponding matrix of coupling names.
///
/// Examples
/// --------
/// >>> from symbolica import S, E
/// >>> from symbolica.community import hep
/// >>> model = hep.Model.standard_model()
/// >>> vertex = model.vertex_rules[0]
/// >>> particles = [model.particle(name) for name in vertex.particles]
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    name = "VertexRule",
    module = "symbolica.community.feynkit",
    frozen,
    from_py_object
)]
#[derive(Clone)]
pub struct PyVertexRule {
    id: VertexRuleId,
    model: Arc<Model>,
}

impl PyVertexRule {
    pub(crate) fn new(id: VertexRuleId, model: Arc<Model>) -> Self {
        Self { id, model }
    }

    fn inner(&self) -> &VertexRule {
        self.model.vertex_rule_by_id(self.id).unwrap()
    }

    pub(crate) fn vertex_rule_name(&self) -> &str {
        &self.inner().name
    }
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[pymethods]
impl PyVertexRule {
    /// Return the vertex-rule name.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``VertexRule`` class example:
    ///
    /// >>> rules_by_name = {rule.name: rule for rule in model.vertex_rules}
    #[getter]
    fn name(&self) -> &str {
        &self.inner().name
    }

    /// Return the ordered particle names attached to the vertex.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``VertexRule`` class example:
    ///
    /// >>> vertex = model.vertex_rules[0]
    /// >>> particles = [model.particle(name) for name in vertex.particles]
    #[getter]
    fn particles(&self) -> Vec<String> {
        self.inner()
            .particles
            .iter()
            .map(|id| self.model.particle_by_id(*id).unwrap().name.clone())
            .collect()
    }

    /// Return the color structures used by the vertex.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``VertexRule`` class example:
    ///
    /// >>> color_basis = vertex.color_structures
    #[getter]
    fn color_structures(&self) -> Vec<PythonExpression> {
        self.inner()
            .color_structures
            .iter()
            .cloned()
            .map(|expr| PythonExpression { expr })
            .collect()
    }

    /// Return the Lorentz-structure names used by the vertex.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``VertexRule`` class example:
    ///
    /// >>> vertex = model.vertex_rules[0]
    /// >>> tensors = [model.lorentz_structure(name) for name in vertex.lorentz_structures]
    #[getter]
    fn lorentz_structures(&self) -> Vec<String> {
        self.inner()
            .lorentz_structures
            .iter()
            .map(|id| {
                self.model
                    .lorentz_structure_by_id(*id)
                    .unwrap()
                    .name
                    .clone()
            })
            .collect()
    }

    /// Return the coupling-name matrix indexed by color and Lorentz structure.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``VertexRule`` class example:
    ///
    /// >>> vertex = model.vertex_rules[0]
    /// >>> names = [name for row in vertex.couplings for name in row if name is not None]
    /// >>> couplings = [model.coupling(name) for name in names]
    #[getter]
    fn couplings(&self) -> Vec<Vec<Option<String>>> {
        self.inner()
            .couplings
            .iter()
            .map(|row| {
                row.iter()
                    .map(|id| id.map(|id| self.model.coupling_by_id(id).unwrap().name.clone()))
                    .collect()
            })
            .collect()
    }

    /// Combine the coupling-order powers referenced by this vertex.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``VertexRule`` class example:
    ///
    /// >>> vertex.coupling_orders()
    /// {'QED': 2}
    fn coupling_orders(&self) -> BTreeMap<String, usize> {
        self.inner().coupling_orders(&self.model)
    }

    /// Summarize the defining data as well as the vertex-rule name.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``VertexRule`` class example:
    ///
    /// >>> print(next(v for v in model.vertex_rules if "e-" in v.particles))
    fn __repr__(&self) -> String {
        let colors = self
            .inner()
            .color_structures
            .iter()
            .map(AtomCore::to_canonical_string)
            .collect::<Vec<_>>();
        let lorentz = self
            .inner()
            .lorentz_structures
            .iter()
            .map(|id| {
                let l = self.model.lorentz_structure_by_id(*id).unwrap();
                format!("{}={}", l.name, l.structure.to_canonical_string())
            })
            .collect::<Vec<_>>();
        let couplings = self
            .inner()
            .couplings
            .iter()
            .enumerate()
            .flat_map(|(c, row)| {
                row.iter().enumerate().filter_map(move |(l, id)| {
                    id.map(|id| {
                        let g = self.model.coupling_by_id(id).unwrap();
                        format!(
                            "({c},{l}): {}={}",
                            g.name,
                            g.expression.to_canonical_string()
                        )
                    })
                })
            })
            .collect::<Vec<_>>();
        format!(
            "VertexRule({:?}, particles={:?}, color=[{}], lorentz=[{}], couplings=[{}])",
            self.name(),
            self.particles(),
            colors.join(", "),
            lorentz.join(", "),
            couplings.join(", ")
        )
    }

    /// Display the model member's physical data and defining expressions.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``VertexRule`` class example:
    ///
    /// >>> from IPython.display import display
    /// >>> display(vertex)
    fn _repr_html_(&self, py: Python<'_>) -> PyResult<String> {
        let particles = self
            .particles()
            .iter()
            .enumerate()
            .map(|(i, p)| format!("{}: {}", i + 1, escape_html(p)))
            .collect::<Vec<_>>()
            .join(" &nbsp;·&nbsp; ");
        let colors = self
            .color_structures()
            .iter()
            .enumerate()
            .map(|(i, c)| {
                Ok(format!(
                    "<div>C<sub>{i}</sub> = {}</div>",
                    model_expression_html(py, c.clone())?
                ))
            })
            .collect::<PyResult<Vec<_>>>()?
            .join("");
        let lorentz = self
            .inner()
            .lorentz_structures
            .iter()
            .enumerate()
            .map(|(i, id)| {
                let l = self.model.lorentz_structure_by_id(*id).unwrap();
                let expr = model_expression_html(
                    py,
                    PythonExpression {
                        expr: l.structure.clone(),
                    },
                )?;
                Ok(format!(
                    "<div>L<sub>{i}</sub> ({}) = {expr}</div>",
                    escape_html(&l.name)
                ))
            })
            .collect::<PyResult<Vec<_>>>()?
            .join("");
        let mut couplings = String::new();
        for (c, row) in self.inner().couplings.iter().enumerate() {
            for (l, id) in row.iter().enumerate() {
                if let Some(id) = id {
                    let g = self.model.coupling_by_id(*id).unwrap();
                    let expr = model_expression_html(
                        py,
                        PythonExpression {
                            expr: g.expression.clone(),
                        },
                    )?;
                    couplings.push_str(&format!("<div><code>{}</code> C<sub>{c}</sub> L<sub>{l}</sub>; &nbsp; <code>{}</code> = {expr}</div>", escape_html(&g.name), escape_html(&g.name)));
                }
            }
        }
        if couplings.is_empty() {
            couplings.push_str("No nonzero couplings");
        }
        Ok(model_record_html(
            "Vertex rule",
            self.name(),
            &[
                ("Particles (leg order)", particles),
                ("Color structures", colors),
                ("Lorentz structures", lorentz),
                ("Coupling terms", couplings),
                (
                    "Orders",
                    escape_html(&format!("{:?}", self.coupling_orders())),
                ),
            ],
        ))
    }

    /// Write the complete text summary to an IPython pretty printer.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``VertexRule`` class example:
    ///
    /// >>> from IPython.lib.pretty import pretty
    /// >>> text = pretty(vertex)
    ///
    /// Parameters
    /// ----------
    /// pretty : object
    ///     IPython pretty printer receiving the text.
    /// cycle : bool
    ///     Whether the object occurs recursively in the current display.
    fn _repr_pretty_(&self, pretty: &Bound<'_, PyAny>, cycle: bool) -> PyResult<()> {
        pretty.call_method1(
            "text",
            (if cycle {
                "VertexRule(...)".to_owned()
            } else {
                self.__repr__()
            },),
        )?;
        Ok(())
    }
}

/// A reusable spin-dependent tensor structure in an interaction rule.
///
/// The structure records the spin representation of each leg and the UFO
/// expression used when constructing a diagram numerator.
///
/// Examples
/// --------
/// >>> from symbolica import S, E
/// >>> from symbolica.community import hep
/// >>> model = hep.Model.standard_model()
/// >>> lorentz = model.lorentz_structures[0]
/// >>> spins = lorentz.spins
/// >>> formula = lorentz.structure
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    name = "LorentzStructure",
    module = "symbolica.community.feynkit",
    frozen,
    from_py_object
)]
#[derive(Clone)]
pub struct PyLorentzStructure {
    id: LorentzStructureId,
    model: Arc<Model>,
}

impl PyLorentzStructure {
    fn new(id: LorentzStructureId, model: Arc<Model>) -> Self {
        Self { id, model }
    }

    fn inner(&self) -> &LorentzStructure {
        self.model.lorentz_structure_by_id(self.id).unwrap()
    }
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[pymethods]
impl PyLorentzStructure {
    /// Return the Lorentz-structure name.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``LorentzStructure`` class example:
    ///
    /// >>> structures = {item.name: item.structure for item in model.lorentz_structures}
    #[getter]
    fn name(&self) -> &str {
        &self.inner().name
    }

    /// Return the UFO ``2S + 1`` spin codes of the structure's particle slots.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``LorentzStructure`` class example:
    ///
    /// >>> vertex = model.vertex_rules[0]
    /// >>> lorentz = model.lorentz_structure(vertex.lorentz_structures[0])
    /// >>> len(lorentz.spins) == len(vertex.particles)
    /// True
    #[getter]
    fn spins(&self) -> Vec<i64> {
        self.inner().spins.clone()
    }

    /// Return the symbolic Lorentz expression.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``LorentzStructure`` class example:
    ///
    /// >>> formula = lorentz.structure
    #[getter]
    fn structure(&self) -> PythonExpression {
        PythonExpression {
            expr: self.inner().structure.clone(),
        }
    }

    /// Summarize the defining data as well as the structure name.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``LorentzStructure`` class example:
    ///
    /// >>> vertex = model.vertex_rules[0]
    /// >>> print(model.lorentz_structure(vertex.lorentz_structures[0]))
    fn __repr__(&self) -> String {
        format!(
            "LorentzStructure({:?}, spins={:?}, structure={})",
            self.name(),
            self.inner().spins,
            self.inner().structure.to_canonical_string()
        )
    }

    /// Display the model member's physical data and defining expressions.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``LorentzStructure`` class example:
    ///
    /// >>> from IPython.display import display
    /// >>> display(lorentz)
    fn _repr_html_(&self, py: Python<'_>) -> PyResult<String> {
        Ok(model_record_html(
            "Lorentz structure",
            self.name(),
            &[
                ("Spins (UFO)", format!("{:?}", self.inner().spins)),
                ("Structure", model_expression_html(py, self.structure())?),
            ],
        ))
    }

    /// Write the complete text summary to an IPython pretty printer.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``LorentzStructure`` class example:
    ///
    /// >>> from IPython.lib.pretty import pretty
    /// >>> text = pretty(lorentz)
    ///
    /// Parameters
    /// ----------
    /// pretty : object
    ///     IPython pretty printer receiving the text.
    /// cycle : bool
    ///     Whether the object occurs recursively in the current display.
    fn _repr_pretty_(&self, pretty: &Bound<'_, PyAny>, cycle: bool) -> PyResult<()> {
        pretty.call_method1(
            "text",
            (if cycle {
                "LorentzStructure(...)".to_owned()
            } else {
                self.__repr__()
            },),
        )?;
        Ok(())
    }
}

/// A particle propagator with symbolic numerator and denominator.
///
/// Propagator formulas are kept in their model representation so that diagram
/// construction can combine them with interaction Lorentz structures.
///
/// Examples
/// --------
/// >>> from symbolica import S, E
/// >>> from symbolica.community import hep
/// >>> model = hep.Model.standard_model()
/// >>> propagator = model.propagators[0]
/// >>> formula = propagator.numerator / propagator.denominator
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    name = "Propagator",
    module = "symbolica.community.feynkit",
    frozen,
    from_py_object
)]
#[derive(Clone)]
pub struct PyPropagator {
    id: PropagatorId,
    model: Arc<Model>,
}

impl PyPropagator {
    pub(crate) fn new(id: PropagatorId, model: Arc<Model>) -> Self {
        Self { id, model }
    }

    fn inner(&self) -> &Propagator {
        self.model.propagator_by_id(self.id).unwrap()
    }
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[pymethods]
impl PyPropagator {
    /// Return the propagator name.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Propagator`` class example:
    ///
    /// >>> propagators_by_name = {item.name: item for item in model.propagators}
    #[getter]
    fn name(&self) -> &str {
        &self.inner().name
    }

    /// Return the name of the particle using this propagator.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Propagator`` class example:
    ///
    /// >>> propagator = model.propagators[0]
    /// >>> particle = model.particle(propagator.particle)
    #[getter]
    fn particle(&self) -> &str {
        &self
            .model
            .particle_by_id(self.inner().particle)
            .unwrap()
            .name
    }

    /// Return the symbolic propagator numerator.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Propagator`` class example:
    ///
    /// >>> formula = propagator.numerator / propagator.denominator
    #[getter]
    fn numerator(&self) -> PythonExpression {
        PythonExpression {
            expr: self.inner().numerator.clone(),
        }
    }

    /// Return the symbolic propagator denominator.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Propagator`` class example:
    ///
    /// >>> inverse_propagator = propagator.denominator
    #[getter]
    fn denominator(&self) -> PythonExpression {
        PythonExpression {
            expr: self.inner().denominator.clone(),
        }
    }

    /// Summarize the defining data as well as the propagator name.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Propagator`` class example:
    ///
    /// >>> print(model.propagators[0])
    fn __repr__(&self) -> String {
        format!(
            "Propagator({:?}, particle={:?}, numerator={}, denominator={})",
            self.name(),
            self.particle(),
            self.inner().numerator.to_canonical_string(),
            self.inner().denominator.to_canonical_string()
        )
    }

    /// Display the model member's physical data and defining expressions.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Propagator`` class example:
    ///
    /// >>> from IPython.display import display
    /// >>> display(propagator)
    fn _repr_html_(&self, py: Python<'_>) -> PyResult<String> {
        Ok(model_record_html(
            "Propagator",
            self.name(),
            &[
                ("Particle", escape_html(self.particle())),
                ("Numerator", model_expression_html(py, self.numerator())?),
                (
                    "Denominator",
                    model_expression_html(py, self.denominator())?,
                ),
            ],
        ))
    }

    /// Write the complete text summary to an IPython pretty printer.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Propagator`` class example:
    ///
    /// >>> from IPython.lib.pretty import pretty
    /// >>> text = pretty(propagator)
    ///
    /// Parameters
    /// ----------
    /// pretty : object
    ///     IPython pretty printer receiving the text.
    /// cycle : bool
    ///     Whether the object occurs recursively in the current display.
    fn _repr_pretty_(&self, pretty: &Bound<'_, PyAny>, cycle: bool) -> PyResult<()> {
        pretty.call_method1(
            "text",
            (if cycle {
                "Propagator(...)".to_owned()
            } else {
                self.__repr__()
            },),
        )?;
        Ok(())
    }
}

/// A named helper function used by model expressions.
///
/// Functions expose their formal arguments and, when available, a symbolic
/// implementation that an evaluator can use during parameter recomputation.
///
/// Examples
/// --------
/// >>> from symbolica import S, E
/// >>> from symbolica.community import hep
/// >>> model = hep.Model.standard_model()
/// >>> functions = {function.name: function.arguments for function in model.functions}
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    name = "ModelFunction",
    module = "symbolica.community.feynkit",
    frozen,
    from_py_object
)]
#[derive(Clone)]
pub struct PyModelFunction {
    id: ModelFunctionId,
    model: Arc<Model>,
}

impl PyModelFunction {
    fn new(id: ModelFunctionId, model: Arc<Model>) -> Self {
        Self { id, model }
    }

    fn inner(&self) -> &ModelFunction {
        self.model.function_by_id(self.id).unwrap()
    }
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[pymethods]
impl PyModelFunction {
    /// Return the symbolic reference used by this model's expressions.
    /// This preserves the parameter or function identity without substituting
    /// its defining expression or numerical value.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import hep
    /// >>> model = hep.Model.standard_model()
    /// >>> import json
    /// >>> definition = json.loads(model.to_json())
    /// >>> definition["functions"] = [{"name": "square", "arguments": ["x"], "expression": "x^2"}]
    /// >>> custom = hep.Model.from_json(json.dumps(definition))
    /// >>> reference = custom.function("square").symbol
    #[getter]
    fn symbol(&self) -> PythonExpression {
        Atom::var(symbol!(&format!("UFO::{}", self.inner().name))).into()
    }

    /// Return the function name.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``ModelFunction`` class example:
    ///
    /// >>> names = [function.name for function in model.functions]
    #[getter]
    fn name(&self) -> &str {
        &self.inner().name
    }

    /// Return the ordered argument names.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``ModelFunction`` class example:
    ///
    /// >>> arguments_by_name = {function.name: function.arguments for function in model.functions}
    #[getter]
    fn arguments(&self) -> Vec<String> {
        self.inner().arguments.clone()
    }

    /// Return the function body, when one is defined by the model.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``ModelFunction`` class example:
    ///
    /// >>> implementations = {function.name: function.expression for function in model.functions}
    #[getter]
    fn expression(&self) -> Option<PythonExpression> {
        self.inner()
            .expression
            .clone()
            .map(|expr| PythonExpression { expr })
    }

    /// Summarize the defining data as well as the function name.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``ModelFunction`` class example:
    ///
    /// >>> summaries = [repr(function) for function in model.functions]
    fn __repr__(&self) -> String {
        format!(
            "ModelFunction({:?}, arguments={:?}, expression={})",
            self.name(),
            self.inner().arguments,
            self.inner()
                .expression
                .as_ref()
                .map_or_else(|| "not defined".to_owned(), AtomCore::to_canonical_string)
        )
    }

    /// Display the model member's physical data and defining expressions.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``ModelFunction`` class example:
    ///
    /// >>> from IPython.display import display
    /// >>> display(model.functions)
    fn _repr_html_(&self, py: Python<'_>) -> PyResult<String> {
        let expression = self
            .expression()
            .map(|e| model_expression_html(py, e))
            .transpose()?
            .unwrap_or_else(|| "Not defined in the model".to_owned());
        Ok(model_record_html(
            "Model function",
            self.name(),
            &[
                ("Arguments", escape_html(&self.inner().arguments.join(", "))),
                ("Definition", expression),
            ],
        ))
    }

    /// Write the complete text summary to an IPython pretty printer.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``ModelFunction`` class example:
    ///
    /// >>> from IPython.lib.pretty import pretty
    /// >>> text = pretty(model.functions)
    ///
    /// Parameters
    /// ----------
    /// pretty : object
    ///     IPython pretty printer receiving the text.
    /// cycle : bool
    ///     Whether the object occurs recursively in the current display.
    fn _repr_pretty_(&self, pretty: &Bound<'_, PyAny>, cycle: bool) -> PyResult<()> {
        pretty.call_method1(
            "text",
            (if cycle {
                "ModelFunction(...)".to_owned()
            } else {
                self.__repr__()
            },),
        )?;
        Ok(())
    }
}

/// A momentum-dependent form factor referenced by interaction rules.
///
/// Form factors keep their model type and symbolic value so specialized
/// amplitude backends can evaluate them at the relevant kinematic point.
///
/// Examples
/// --------
/// >>> from symbolica import S, E
/// >>> from symbolica.community import hep
/// >>> model = hep.Model.standard_model()
/// >>> form_factors = {ff.name: ff.value for ff in model.form_factors}
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    name = "FormFactor",
    module = "symbolica.community.feynkit",
    frozen,
    from_py_object
)]
#[derive(Clone)]
pub struct PyFormFactor {
    id: ModelFormFactorId,
    model: Arc<Model>,
}

impl PyFormFactor {
    fn new(id: ModelFormFactorId, model: Arc<Model>) -> Self {
        Self { id, model }
    }

    fn inner(&self) -> &ModelFormFactor {
        self.model.form_factor_by_id(self.id).unwrap()
    }
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[pymethods]
impl PyFormFactor {
    /// Return the form-factor name.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``FormFactor`` class example:
    ///
    /// >>> names = [ff.name for ff in model.form_factors]
    #[getter]
    fn name(&self) -> &str {
        &self.inner().name
    }

    /// Return the model-defined form-factor type, when present.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``FormFactor`` class example:
    ///
    /// >>> types = {ff.name: ff.type_name for ff in model.form_factors}
    #[getter]
    fn type_name(&self) -> Option<String> {
        self.inner().type_name.clone()
    }

    /// Return the symbolic form-factor value, when present.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``FormFactor`` class example:
    ///
    /// >>> formulas = {ff.name: ff.value for ff in model.form_factors}
    #[getter]
    fn value(&self) -> Option<PythonExpression> {
        self.inner()
            .value
            .clone()
            .map(|expr| PythonExpression { expr })
    }

    /// Summarize the defining data as well as the form-factor name.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``FormFactor`` class example:
    ///
    /// >>> summaries = [repr(ff) for ff in model.form_factors]
    fn __repr__(&self) -> String {
        format!(
            "FormFactor({:?}, type={:?}, value={})",
            self.name(),
            self.inner().type_name,
            self.inner()
                .value
                .as_ref()
                .map_or_else(|| "not defined".to_owned(), AtomCore::to_canonical_string)
        )
    }

    /// Display the model member's physical data and defining expressions.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``FormFactor`` class example:
    ///
    /// >>> from IPython.display import display
    /// >>> display(model.form_factors)
    fn _repr_html_(&self, py: Python<'_>) -> PyResult<String> {
        let expression = self
            .value()
            .map(|e| model_expression_html(py, e))
            .transpose()?
            .unwrap_or_else(|| "Not defined in the model".to_owned());
        Ok(model_record_html(
            "Form factor",
            self.name(),
            &[
                (
                    "Type",
                    escape_html(self.inner().type_name.as_deref().unwrap_or("unspecified")),
                ),
                ("Value", expression),
            ],
        ))
    }

    /// Write the complete text summary to an IPython pretty printer.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``FormFactor`` class example:
    ///
    /// >>> from IPython.lib.pretty import pretty
    /// >>> text = pretty(model.form_factors)
    ///
    /// Parameters
    /// ----------
    /// pretty : object
    ///     IPython pretty printer receiving the text.
    /// cycle : bool
    ///     Whether the object occurs recursively in the current display.
    fn _repr_pretty_(&self, pretty: &Bound<'_, PyAny>, cycle: bool) -> PyResult<()> {
        pretty.call_method1(
            "text",
            (if cycle {
                "FormFactor(...)".to_owned()
            } else {
                self.__repr__()
            },),
        )?;
        Ok(())
    }
}

/// A named symbolic expression awaiting numerical model evaluation.
///
/// These records describe internal parameters or couplings passed to a custom
/// ``Model.recompute_with`` callback.
///
/// Examples
/// --------
/// >>> from symbolica import S, E
/// >>> from symbolica.community import hep
/// >>> model = hep.Model.phi4()
/// >>> process = model.process(["phi", "phi"], ["phi", "phi"])
/// >>> requests = []
/// >>> def evaluate(request):
/// ...     requests.append(request)
/// ...     return hep.EvaluatedValues(
/// ...         couplings={"SCALAR_COUPLING": (0.0, -1.0)},
/// ...     )
/// >>> updated_model = model.recompute_with(evaluate)
/// >>> request = requests[0]
/// >>> item = request.couplings[0]
/// >>> formula = item.expression
/// >>> name = item.name
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    name = "ModelExpression",
    module = "symbolica.community.feynkit",
    frozen,
    from_py_object
)]
#[derive(Clone)]
pub struct PyModelExpression {
    inner: ModelExpression,
}

impl From<ModelExpression> for PyModelExpression {
    fn from(inner: ModelExpression) -> Self {
        Self { inner }
    }
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[pymethods]
impl PyModelExpression {
    /// Return the name assigned to the expression.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``ModelExpression`` class example:
    ///
    /// >>> formulas = {item.name: item.expression for item in request.couplings}
    #[getter]
    fn name(&self) -> &str {
        &self.inner.name
    }

    /// Return the symbolic expression text.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``ModelExpression`` class example:
    ///
    /// >>> formula = item.expression
    #[getter]
    fn expression(&self) -> PythonExpression {
        crate::graph::parse_symbolic_annotation(&self.inner.expression)
            .expect("evaluation requests contain canonical Symbolica expressions")
    }
}

/// The complete input supplied to a custom model evaluator.
///
/// An evaluation request contains already-known numerical parameters and the
/// internal parameters, couplings, functions, and form factors still needed.
///
/// Examples
/// --------
/// This evaluator implements the built-in scalar model at its default ``lam=1``.
/// A callback must return every requested internal parameter and coupling.
///
/// >>> from symbolica import S, E
/// >>> from symbolica.community import hep
/// >>> model = hep.Model.phi4()
/// >>> process = model.process(["phi", "phi"], ["phi", "phi"])
/// >>> requests = []
/// >>> def evaluate(request):
/// ...     requests.append(request)
/// ...     return hep.EvaluatedValues(
/// ...         couplings={"SCALAR_COUPLING": (0.0, -1.0)},
/// ...     )
/// >>> updated_model = model.recompute_with(evaluate)
/// >>> request = requests[0]
/// >>> formulas = {item.name: item.expression for item in request.couplings}
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    name = "EvaluationRequest",
    module = "symbolica.community.feynkit",
    frozen,
    from_py_object
)]
#[derive(Clone)]
pub struct PyEvaluationRequest {
    inner: EvaluationRequest,
    model: Arc<Model>,
}

impl PyEvaluationRequest {
    fn new(inner: EvaluationRequest, model: Arc<Model>) -> Self {
        Self { inner, model }
    }
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[pymethods]
impl PyEvaluationRequest {
    /// Return the known parameter values keyed by name.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``EvaluationRequest`` class example:
    ///
    /// >>> def evaluate(request):
    /// ...     alpha_s = request.known_parameters["aS"]
    /// ...     return hep.EvaluatedValues()
    #[getter]
    fn known_parameters(&self) -> BTreeMap<String, (f64, f64)> {
        self.inner
            .known_parameters
            .iter()
            .map(|(name, value)| (name.clone(), (value.re, value.im)))
            .collect()
    }

    /// Return the derived parameter expressions that must be evaluated.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``EvaluationRequest`` class example:
    ///
    /// >>> formulas = {item.name: item.expression for item in request.internal_parameters}
    #[getter]
    fn internal_parameters(&self) -> Vec<PyModelExpression> {
        self.inner
            .internal_parameters
            .iter()
            .cloned()
            .map(Into::into)
            .collect()
    }

    /// Return the coupling expressions that must be evaluated.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``EvaluationRequest`` class example:
    ///
    /// >>> coupling_formulas = {item.name: item.expression for item in request.couplings}
    #[getter]
    fn couplings(&self) -> Vec<PyModelExpression> {
        self.inner
            .couplings
            .iter()
            .cloned()
            .map(Into::into)
            .collect()
    }

    /// Return the function definitions available during evaluation.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``EvaluationRequest`` class example:
    ///
    /// >>> helper_functions = {function.name: function for function in request.functions}
    #[getter]
    fn functions(&self) -> Vec<PyModelFunction> {
        self.inner
            .functions
            .iter()
            .map(|function| {
                PyModelFunction::new(
                    self.model.function_id(&function.name).unwrap(),
                    Arc::clone(&self.model),
                )
            })
            .collect()
    }

    /// Return the form-factor definitions available during evaluation.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``EvaluationRequest`` class example:
    ///
    /// >>> form_factors = {factor.name: factor for factor in request.form_factors}
    #[getter]
    fn form_factors(&self) -> Vec<PyFormFactor> {
        self.inner
            .form_factors
            .iter()
            .map(|form_factor| {
                PyFormFactor::new(
                    self.model.form_factor_id(&form_factor.name).unwrap(),
                    Arc::clone(&self.model),
                )
            })
            .collect()
    }
}

/// Numerical values returned by a custom model evaluator.
///
/// Values are complex numbers represented as ``(real, imaginary)`` pairs and
/// are applied atomically to the model after the callback completes.
///
/// Examples
/// --------
/// >>> from symbolica.community import hep
/// >>> values = hep.EvaluatedValues(couplings={"GC_1": (0.3, 0.0)})
///
/// Parameters
/// ----------
/// internal_parameters : dict[str, tuple[float, float]], optional
///     Evaluated internal parameter values keyed by name.
/// couplings : dict[str, tuple[float, float]], optional
///     Evaluated coupling values keyed by name.
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    name = "EvaluatedValues",
    module = "symbolica.community.feynkit",
    frozen,
    from_py_object
)]
#[derive(Clone, Default)]
pub struct PyEvaluatedValues {
    inner: EvaluatedValues,
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[pymethods]
impl PyEvaluatedValues {
    /// Create values returned by a custom model evaluator.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``EvaluatedValues`` class example:
    ///
    /// >>> values = hep.EvaluatedValues(couplings={"GC_1": (1.0, 0.0)})
    ///
    /// Parameters
    /// ----------
    /// internal_parameters : dict[str, tuple[float, float]] or None
    ///     Evaluated internal parameters keyed by name.
    /// couplings : dict[str, tuple[float, float]] or None
    ///     Evaluated couplings keyed by name.
    #[new]
    #[pyo3(signature = (*, internal_parameters=None, couplings=None))]
    fn new(
        internal_parameters: Option<BTreeMap<String, (f64, f64)>>,
        couplings: Option<BTreeMap<String, (f64, f64)>>,
    ) -> Self {
        Self {
            inner: EvaluatedValues {
                internal_parameters: internal_parameters
                    .unwrap_or_default()
                    .into_iter()
                    .map(|(name, (re, im))| (name, ComplexValue::new(re, im)))
                    .collect(),
                couplings: couplings
                    .unwrap_or_default()
                    .into_iter()
                    .map(|(name, (re, im))| (name, ComplexValue::new(re, im)))
                    .collect(),
            },
        }
    }

    /// Return the evaluated internal parameters keyed by name.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``EvaluatedValues`` class example:
    ///
    /// >>> values = hep.EvaluatedValues(internal_parameters={"alpha": (0.1, 0.0)})
    /// >>> assert values.internal_parameters["alpha"] == (0.1, 0.0)
    #[getter]
    fn internal_parameters(&self) -> BTreeMap<String, (f64, f64)> {
        self.inner
            .internal_parameters
            .iter()
            .map(|(name, value)| (name.clone(), (value.re, value.im)))
            .collect()
    }

    /// Return the evaluated couplings keyed by name.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``EvaluatedValues`` class example:
    ///
    /// >>> values = hep.EvaluatedValues(couplings={"GC_1": (0.3, 0.0)})
    /// >>> assert values.couplings["GC_1"] == (0.3, 0.0)
    #[getter]
    fn couplings(&self) -> BTreeMap<String, (f64, f64)> {
        self.inner
            .couplings
            .iter()
            .map(|(name, value)| (name.clone(), (value.re, value.im)))
            .collect()
    }
}

struct PythonModelEvaluator<'a, 'py> {
    callable: &'a Bound<'py, PyAny>,
    model: Arc<Model>,
}

impl ModelEvaluator for PythonModelEvaluator<'_, '_> {
    type Error = PyErr;

    fn evaluate(&mut self, request: EvaluationRequest) -> Result<EvaluatedValues, Self::Error> {
        Ok(self
            .callable
            .call1((PyEvaluationRequest::new(request, Arc::clone(&self.model)),))?
            .extract::<PyEvaluatedValues>()
            .map(|values| values.inner)?)
    }
}

/// Mutable external parameter values for a particle model.
///
/// A parameter card can be initialized from a model, edited by parameter name,
/// serialized as JSON, and atomically applied to a model copy.
///
/// Examples
/// --------
/// >>> from symbolica import S, E
/// >>> from symbolica.community import hep
/// >>> model = hep.Model.standard_model()
/// >>> card = model.default_parameter_card()
/// >>> card.set("MM", 0.105658, 0.0)
/// >>> shifted_model = model.with_parameter_card(card)
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    name = "ParameterCard",
    module = "symbolica.community.feynkit",
    from_py_object
)]
#[derive(Clone, Default)]
pub struct PyParameterCard {
    pub(crate) inner: ParameterCard,
}

impl From<ParameterCard> for PyParameterCard {
    fn from(inner: ParameterCard) -> Self {
        Self { inner }
    }
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[pymethods]
impl PyParameterCard {
    /// Create an empty parameter card.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``ParameterCard`` class example:
    ///
    /// >>> card = hep.ParameterCard()
    #[new]
    fn new() -> Self {
        Self::default()
    }

    /// Parse a parameter card from its JSON representation.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``ParameterCard`` class example:
    ///
    /// >>> card = hep.ParameterCard.from_json('{"mass": [1.0, 0.0]}')
    ///
    /// Parameters
    /// ----------
    /// json : str
    ///     Serialized parameter-card object.
    #[staticmethod]
    fn from_json(json: &str) -> PyResult<Self> {
        ParameterCard::from_json(json)
            .map(Self::from)
            .map_err(error::model)
    }

    /// Load a parameter card from a JSON file.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``ParameterCard`` class example:
    ///
    /// >>> from pathlib import Path
    /// >>> from tempfile import TemporaryDirectory
    /// >>> with TemporaryDirectory() as directory:
    /// ...     path = Path(directory) / "parameters.json"
    /// ...     path.write_text(card.to_json())
    /// ...     restored = hep.ParameterCard.from_path(path)
    ///
    /// Parameters
    /// ----------
    /// path : str or os.PathLike
    ///     Path to the JSON parameter card.
    #[staticmethod]
    fn from_path(py: Python<'_>, path: PathBuf) -> PyResult<Self> {
        py.detach(move || ParameterCard::from_path(path))
            .map(Self::from)
            .map_err(error::model)
    }

    /// Return a parameter's complex value, or ``None`` when it is absent.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``ParameterCard`` class example:
    ///
    /// >>> card.get("mass")
    ///
    /// Parameters
    /// ----------
    /// name : str
    ///     Parameter name to look up.
    fn get(&self, name: &str) -> Option<(f64, f64)> {
        self.inner.get(name).map(|value| (value.re, value.im))
    }

    /// Store a complex parameter value.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``ParameterCard`` class example:
    ///
    /// >>> card.set("mass", 1.0)
    ///
    /// Parameters
    /// ----------
    /// name : str
    ///     Parameter name to store.
    /// real : float
    ///     Real component of the value.
    /// imaginary : float
    ///     Imaginary component of the value.
    #[pyo3(signature = (name, real, imaginary=0.0))]
    fn set(&mut self, name: String, real: f64, imaginary: f64) {
        self.inner.insert(name, ComplexValue::new(real, imaginary));
    }

    /// Remove and return a parameter value, if present.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``ParameterCard`` class example:
    ///
    /// >>> card.remove("mass")
    ///
    /// Parameters
    /// ----------
    /// name : str
    ///     Parameter name to remove.
    fn remove(&mut self, name: &str) -> Option<(f64, f64)> {
        self.inner.remove(name).map(|value| (value.re, value.im))
    }

    /// Return all parameter-card entries.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``ParameterCard`` class example:
    ///
    /// >>> card.items()
    fn items(&self) -> Vec<(String, (f64, f64))> {
        self.inner
            .iter()
            .map(|(name, value)| (name.clone(), (value.re, value.im)))
            .collect()
    }

    /// Serialize the parameter card as JSON.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``ParameterCard`` class example:
    ///
    /// >>> card.to_json(pretty=False)
    ///
    /// Parameters
    /// ----------
    /// pretty : bool
    ///     Indent the output when true.
    #[pyo3(signature = (pretty=true))]
    fn to_json(&self, pretty: bool) -> PyResult<String> {
        if pretty {
            self.inner.to_json_pretty()
        } else {
            self.inner.to_json()
        }
        .map_err(error::model)
    }

    /// Write the parameter card to a JSON file.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``ParameterCard`` class example:
    ///
    /// >>> card.write_json("parameters.json")
    ///
    /// Parameters
    /// ----------
    /// path : str or os.PathLike
    ///     Destination for the JSON parameter card.
    fn write_json(&self, py: Python<'_>, path: PathBuf) -> PyResult<()> {
        let card = self.inner.clone();
        py.detach(move || card.write_json(path))
            .map_err(error::model)
    }

    /// Return the number of entries in the parameter card.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``ParameterCard`` class example:
    ///
    /// >>> number_of_external_inputs = len(model.default_parameter_card())
    fn __len__(&self) -> usize {
        self.inner.len()
    }
}

/// A loaded particle model ready for diagram generation.
///
/// Models provide typed access to particles, parameters, couplings, interaction
/// rules, propagators, and numerical parameter-card updates.
///
/// Examples
/// --------
/// Built-in models need no external files. To import your own UFO directory,
/// see ``UfoLoader``; to restore a normalized JSON model, use ``Model(path)``
/// or ``Model.from_json``.
///
/// >>> from symbolica import S, E
/// >>> from symbolica.community import hep
/// >>> model = hep.Model.standard_model()
/// >>> photon = model.particle("a")
/// >>> process = model.process(["e-", "e+"], ["mu-", "mu+"])
/// >>> result = process.generate_diagrams()
/// >>> assert result.report.completed
/// >>> scalar_model = hep.Model.phi4()
/// >>> assert scalar_model.particle("phi").spin == 1
///
/// Parameters
/// ----------
/// path : str or os.PathLike
///     Path to a normalized HEP JSON model.
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    name = "Model",
    module = "symbolica.community.feynkit",
    frozen,
    from_py_object
)]
#[derive(Clone)]
pub struct PyModel {
    pub(crate) inner: Arc<Model>,
}

impl From<Model> for PyModel {
    fn from(inner: Model) -> Self {
        Arc::new(inner).into()
    }
}

impl From<Arc<Model>> for PyModel {
    fn from(inner: Arc<Model>) -> Self {
        for parameter in inner.parameters() {
            DisplaySettings::register_latex_name(
                symbol!(&format!("UFO::{}", parameter.name)),
                parameter.texname.as_deref().unwrap_or_default(),
            );
        }
        Self { inner }
    }
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "python_stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[pymethods]
impl PyModel {
    /// Load the complete embedded Standard Model with default parameters.
    /// No model files or UFO installation are required.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Model`` class example:
    ///
    /// >>> model = hep.Model.standard_model()
    #[staticmethod]
    fn standard_model() -> Self {
        Model::standard_model().into()
    }

    /// Load QCD with six quark flavors, gluons, and gluon ghosts, filtered from
    /// the Standard Model. Preserves its parameters and their default values.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Model`` class example:
    ///
    /// >>> model = hep.Model.qcd()
    #[staticmethod]
    fn qcd() -> Self {
        Model::qcd().into()
    }

    /// Load QED with photons and all charged fermions (including quarks), filtered
    /// from the Standard Model. Preserves its parameters and their default values.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Model`` class example:
    ///
    /// >>> model = hep.Model.qed()
    #[staticmethod]
    fn qed() -> Self {
        Model::qed().into()
    }

    /// Load the electroweak Standard Model, retaining quarks, Higgs, Goldstones,
    /// and electroweak ghosts, and removing gluons and gluon ghosts.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Model`` class example:
    ///
    /// >>> model = hep.Model.electroweak()
    #[staticmethod]
    fn electroweak() -> Self {
        Model::electroweak().into()
    }

    /// Load strong and electromagnetic interactions of all quarks and charged
    /// leptons, including gluon ghosts, with Standard Model parameters.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Model`` class example:
    ///
    /// >>> model = hep.Model.qcd_qed()
    #[staticmethod]
    fn qcd_qed() -> Self {
        Model::qcd_qed().into()
    }

    /// Load pure SU(3) Yang-Mills theory: gluons and gluon ghosts without quarks.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Model`` class example:
    ///
    /// >>> model = hep.Model.yang_mills()
    #[staticmethod]
    fn yang_mills() -> Self {
        Model::yang_mills().into()
    }

    /// Load a real scalar phi with L_int = -g phi^3 / 3!.
    /// The external parameters mass and g default to one; the width is zero.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Model`` class example:
    ///
    /// >>> model = hep.Model.phi3()
    #[staticmethod]
    fn phi3() -> Self {
        Model::phi3().into()
    }

    /// Load a real scalar phi with L_int = -lam phi^4 / 4!.
    /// The external parameters mass and lam default to one; the width is zero.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Model`` class example:
    ///
    /// >>> model = hep.Model.phi4()
    #[staticmethod]
    fn phi4() -> Self {
        Model::phi4().into()
    }

    /// Load a real scalar phi with L_int = -g phi^3 / 3! - lam phi^4 / 4!.
    /// Independent parameters mass, g, and lam default to one; the width is zero.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Model`` class example:
    ///
    /// >>> model = hep.Model.phi_3_4()
    #[staticmethod]
    fn phi_3_4() -> Self {
        Model::phi_3_4().into()
    }

    /// Load scalar QED in Feynman gauge with a, phi+, and phi-.
    /// Parameters mass=e=1 and lam=0; the scalar potential is
    /// mass^2 |phi|^2 + lam |phi|^4 / 4, and all widths vanish.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Model`` class example:
    ///
    /// >>> model = hep.Model.scalar_qed()
    #[staticmethod]
    fn scalar_qed() -> Self {
        Model::scalar_qed().into()
    }

    /// Load a model from a normalized HEP JSON file.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Model`` class example:
    ///
    /// >>> from pathlib import Path
    /// >>> from tempfile import TemporaryDirectory
    /// >>> with TemporaryDirectory() as directory:
    /// ...     path = Path(directory) / "model.json"
    /// ...     model.write_json(path)
    /// ...     restored = hep.Model(path)
    /// ...     assert restored.name == model.name
    ///
    /// Parameters
    /// ----------
    /// path : str or os.PathLike
    ///     Path to the JSON model.
    #[new]
    fn new(py: Python<'_>, path: PathBuf) -> PyResult<Self> {
        py.detach(move || Model::from_path(path))
            .map(Self::from)
            .map_err(error::model)
    }

    /// Parse a model from its JSON representation.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Model`` class example:
    ///
    /// >>> model_json = model.to_json()
    /// >>> restored = hep.Model.from_json(model_json)
    /// >>> assert restored.name == model.name
    ///
    /// Parameters
    /// ----------
    /// json : str
    ///     Serialized model object.
    #[staticmethod]
    fn from_json(json: &str) -> PyResult<Self> {
        Model::from_json(json).map(Self::from).map_err(error::model)
    }

    /// Define a process with model-validated external states and sector restrictions.
    /// The returned process is immutable; loops and other calculation choices are
    /// arguments to its generate_diagrams, generate_amplitude and generate_cross_section methods.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Model`` class example:
    ///
    /// >>> process = model.process(["e-", "e+"], ["a", "a"])
    /// >>> amplitude = process.generate_amplitude(loops=0)
    ///
    /// Parameters
    /// ----------
    /// incoming : sequence[Particle | ParticleSelector | str | int]
    ///     Ordered incoming external states.
    /// outgoing : sequence[Particle | ParticleSelector | str | int]
    ///     Ordered outgoing external states.
    /// particle_veto : sequence[Particle | ParticleSelector | str | int] or None, optional
    ///     Excluded species, including their antiparticles.
    /// vertex_allow : sequence[VertexRule | str] or None, optional
    ///     Allowed interactions. None allows all; an empty list allows none.
    /// vertex_veto : sequence[VertexRule | str] or None, optional
    ///     Excluded interactions.
    #[pyo3(signature = (incoming, outgoing, *, particle_veto=None, vertex_allow=None, vertex_veto=None))]
    fn process(
        &self,
        incoming: Vec<SelectorInput>,
        outgoing: Vec<SelectorInput>,
        particle_veto: Option<Vec<SelectorInput>>,
        vertex_allow: Option<Vec<VertexInput>>,
        vertex_veto: Option<Vec<VertexInput>>,
    ) -> PyResult<PyProcess> {
        PyProcess::from_model(
            self,
            incoming,
            outgoing,
            particle_veto,
            vertex_allow,
            vertex_veto,
        )
    }

    /// Return the model name.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Model`` class example:
    ///
    /// >>> model = hep.Model.standard_model()
    /// >>> model_name = model.name
    #[getter]
    fn name(&self) -> &str {
        self.inner.name()
    }

    /// Return the applied restriction name, when present.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Model`` class example:
    ///
    /// >>> model = hep.Model.standard_model()
    /// >>> restriction = model.restriction
    #[getter]
    fn restriction(&self) -> Option<String> {
        self.inner.restriction().map(str::to_owned)
    }

    /// Return all particle species in deterministic model order.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Model`` class example:
    ///
    /// >>> fermions = [particle for particle in model.particles if particle.spin == 2]
    #[getter]
    fn particles(&self) -> Vec<PyParticle> {
        self.inner
            .particles()
            .iter()
            .enumerate()
            .map(|(index, _)| {
                PyParticle::new(ParticleId::from_index(index), Arc::clone(&self.inner))
            })
            .collect()
    }

    /// Look up a particle by name.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Model`` class example:
    ///
    /// >>> electron = model.particle("e-")
    /// >>> assert electron.pdg_code == 11
    ///
    /// Parameters
    /// ----------
    /// name : str
    ///     Particle name or antiname.
    fn particle(&self, name: &str) -> PyResult<PyParticle> {
        self.inner
            .particle_id(name)
            .map(|id| PyParticle::new(id, Arc::clone(&self.inner)))
            .map_err(error::model)
    }

    /// Look up a particle by PDG code.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Model`` class example:
    ///
    /// >>> particle = model.particle_by_pdg(11)
    ///
    /// Parameters
    /// ----------
    /// pdg : int
    ///     Signed PDG particle code.
    fn particle_by_pdg(&self, pdg: i64) -> PyResult<PyParticle> {
        self.inner
            .particle_id_by_pdg(pdg)
            .map(|id| PyParticle::new(id, Arc::clone(&self.inner)))
            .map_err(error::model)
    }

    /// Return all external and internal parameters in model order.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Model`` class example:
    ///
    /// >>> external = [p for p in model.parameters if p.nature == hep.ParameterNature.EXTERNAL]
    #[getter]
    fn parameters(&self) -> Vec<PyParameter> {
        self.inner
            .parameters()
            .iter()
            .enumerate()
            .map(|(index, _)| {
                PyParameter::new(ParameterId::from_index(index), Arc::clone(&self.inner))
            })
            .collect()
    }

    /// Look up a parameter by name.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Model`` class example:
    ///
    /// >>> mass = model.parameter("MM")
    /// >>> assert mass.nature == hep.ParameterNature.EXTERNAL
    ///
    /// Parameters
    /// ----------
    /// name : str
    ///     Parameter name.
    fn parameter(&self, name: &str) -> PyResult<PyParameter> {
        self.inner
            .parameter_id(name)
            .map(|id| PyParameter::new(id, Arc::clone(&self.inner)))
            .map_err(error::model)
    }

    /// Return all interaction couplings in model order.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Model`` class example:
    ///
    /// >>> qed = [c for c in model.couplings if c.orders.get("QED", 0) > 0]
    #[getter]
    fn couplings(&self) -> Vec<PyCoupling> {
        self.inner
            .couplings()
            .iter()
            .enumerate()
            .map(|(index, _)| {
                PyCoupling::new(CouplingId::from_index(index), Arc::clone(&self.inner))
            })
            .collect()
    }

    /// Look up a coupling by name.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Model`` class example:
    ///
    /// >>> coupling = model.coupling("GC_1")
    ///
    /// Parameters
    /// ----------
    /// name : str
    ///     Coupling name.
    fn coupling(&self, name: &str) -> PyResult<PyCoupling> {
        self.inner
            .coupling_id(name)
            .map(|id| PyCoupling::new(id, Arc::clone(&self.inner)))
            .map_err(error::model)
    }

    /// Replace named UFO coefficients by their analytic model expressions.
    ///
    /// The input and the stored model are unchanged.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Model`` class example:
    ///
    /// >>> analytic = model.expand_couplings(S("UFO::GC_11"))
    ///
    /// Parameters
    /// ----------
    /// expression : Expression
    ///     Symbolica expression containing named couplings from this model.
    #[gen_stub(skip)]
    fn expand_couplings(
        &self,
        py: Python<'_>,
        expression: &Bound<'_, PyAny>,
    ) -> PyResult<Py<PyAny>> {
        if let Ok(tensor) = expression.extract::<PyRef<'_, TensorExpression>>() {
            let result = self
                .inner
                .expand_couplings(tensor.structured().expression());
            return TensorExpression::preserving_interface(&tensor, py, result).map(Py::into_any);
        }
        let input = expression.extract::<PyRef<'_, PythonExpression>>()?;
        let result = self.inner.expand_couplings(&input.expr);
        Py::new(py, PythonExpression { expr: result }).map(Py::into_any)
    }

    /// Return all interaction vertex rules in model order.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Model`` class example:
    ///
    /// >>> electron_vertices = [v for v in model.vertex_rules if "e-" in v.particles]
    #[getter]
    fn vertex_rules(&self) -> Vec<PyVertexRule> {
        self.inner
            .vertex_rules()
            .iter()
            .enumerate()
            .map(|(index, _)| {
                PyVertexRule::new(VertexRuleId::from_index(index), Arc::clone(&self.inner))
            })
            .collect()
    }

    /// Look up a vertex rule by name.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Model`` class example:
    ///
    /// >>> vertex = model.vertex_rule("V_1")
    ///
    /// Parameters
    /// ----------
    /// name : str
    ///     Vertex-rule name.
    fn vertex_rule(&self, name: &str) -> PyResult<PyVertexRule> {
        self.inner
            .vertex_rule_id(name)
            .map(|id| PyVertexRule::new(id, Arc::clone(&self.inner)))
            .map_err(error::model)
    }

    /// Return all reusable Lorentz structures in model order.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Model`` class example:
    ///
    /// >>> lorentz_by_name = {item.name: item for item in model.lorentz_structures}
    #[getter]
    fn lorentz_structures(&self) -> Vec<PyLorentzStructure> {
        self.inner
            .lorentz_structures()
            .iter()
            .enumerate()
            .map(|(index, _)| {
                PyLorentzStructure::new(
                    LorentzStructureId::from_index(index),
                    Arc::clone(&self.inner),
                )
            })
            .collect()
    }

    /// Look up a Lorentz structure by name.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Model`` class example:
    ///
    /// >>> lorentz = model.lorentz_structure("FFV1")
    ///
    /// Parameters
    /// ----------
    /// name : str
    ///     Lorentz-structure name.
    fn lorentz_structure(&self, name: &str) -> PyResult<PyLorentzStructure> {
        self.inner
            .lorentz_structure_id(name)
            .map(|id| PyLorentzStructure::new(id, Arc::clone(&self.inner)))
            .map_err(error::model)
    }

    /// Return all model-defined propagators in model order.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Model`` class example:
    ///
    /// >>> propagators_by_name = {item.name: item for item in model.propagators}
    #[getter]
    fn propagators(&self) -> Vec<PyPropagator> {
        self.inner
            .propagators()
            .iter()
            .enumerate()
            .map(|(index, _)| {
                PyPropagator::new(PropagatorId::from_index(index), Arc::clone(&self.inner))
            })
            .collect()
    }

    /// Look up a propagator by name.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Model`` class example:
    ///
    /// >>> name = model.propagators[0].name
    /// >>> propagator = model.propagator(name)
    ///
    /// Parameters
    /// ----------
    /// name : str
    ///     Propagator name.
    fn propagator(&self, name: &str) -> PyResult<PyPropagator> {
        self.inner
            .propagator_id(name)
            .map(|id| PyPropagator::new(id, Arc::clone(&self.inner)))
            .map_err(error::model)
    }

    /// Return all helper functions available to model expressions.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Model`` class example:
    ///
    /// >>> functions = {function.name: function for function in model.functions}
    #[getter]
    fn functions(&self) -> Vec<PyModelFunction> {
        self.inner
            .functions()
            .iter()
            .enumerate()
            .map(|(index, _)| {
                PyModelFunction::new(ModelFunctionId::from_index(index), Arc::clone(&self.inner))
            })
            .collect()
    }

    /// Look up a model function by name.
    ///
    /// Examples
    /// --------
    /// Using ``model`` from the class example. A model may have no helper functions:
    ///
    /// >>> functions = [model.function(item.name) for item in model.functions]
    ///
    /// Parameters
    /// ----------
    /// name : str
    ///     Function name.
    fn function(&self, name: &str) -> PyResult<PyModelFunction> {
        self.inner
            .function_id(name)
            .map(|id| PyModelFunction::new(id, Arc::clone(&self.inner)))
            .map_err(error::model)
    }

    /// Return all momentum-dependent model form factors in model order.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Model`` class example:
    ///
    /// >>> form_factors = {factor.name: factor for factor in model.form_factors}
    #[getter]
    fn form_factors(&self) -> Vec<PyFormFactor> {
        self.inner
            .form_factors()
            .iter()
            .enumerate()
            .map(|(index, _)| {
                PyFormFactor::new(
                    ModelFormFactorId::from_index(index),
                    Arc::clone(&self.inner),
                )
            })
            .collect()
    }

    /// Look up a form factor by name.
    ///
    /// Examples
    /// --------
    /// Using ``model`` from the class example. A model may have no form factors:
    ///
    /// >>> form_factors = [model.form_factor(item.name) for item in model.form_factors]
    ///
    /// Parameters
    /// ----------
    /// name : str
    ///     Form-factor name.
    fn form_factor(&self, name: &str) -> PyResult<PyFormFactor> {
        self.inner
            .form_factor_id(name)
            .map(|id| PyFormFactor::new(id, Arc::clone(&self.inner)))
            .map_err(error::model)
    }

    /// Build a parameter card from the model's current external values.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Model`` class example:
    ///
    /// >>> card = model.default_parameter_card()
    fn default_parameter_card(&self) -> PyResult<PyParameterCard> {
        self.inner
            .default_parameter_card()
            .map(Into::into)
            .map_err(error::model)
    }

    /// Return a copy with a parameter card applied atomically.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Model`` class example:
    ///
    /// >>> card = model.default_parameter_card()
    /// >>> card.set("MM", 0.105658, 0.0)
    /// >>> updated = model.with_parameter_card(card)
    ///
    /// Parameters
    /// ----------
    /// card : ParameterCard
    ///     External parameter values to apply.
    /// evaluator : Callable[[EvaluationRequest], EvaluatedValues] or None
    ///     Optional evaluator used to recompute dependent values.
    #[pyo3(signature = (card, evaluator=None))]
    fn with_parameter_card(
        &self,
        card: &PyParameterCard,
        #[gen_stub(override_type(
            type_repr = "collections.abc.Callable[[EvaluationRequest], EvaluatedValues] | None",
            imports = ("collections.abc")
        ))]
        evaluator: Option<&Bound<'_, PyAny>>,
    ) -> PyResult<Self> {
        let mut model = self.inner.as_ref().clone();
        if let Some(callable) = evaluator {
            let mut evaluator = PythonModelEvaluator {
                callable,
                model: Arc::clone(&self.inner),
            };
            model
                .apply_parameter_card_with(&card.inner, &mut evaluator)
                .map_err(error::recompute)?;
        } else {
            model
                .apply_parameter_card(&card.inner)
                .map_err(error::model)?;
        }
        Ok(model.into())
    }

    /// Return a copy with all dependent values recomputed by a callback.
    ///
    /// Examples
    /// --------
    /// This callback evaluates the built-in scalar model at its default coupling
    /// ``lam=1``. See ``EvaluationRequest`` for callback inputs.
    ///
    /// >>> from symbolica import S, E
    /// >>> from symbolica.community import hep
    /// >>> model = hep.Model.phi4()
    /// >>> process = model.process(["phi", "phi"], ["phi", "phi"])
    /// >>> requests = []
    /// >>> def evaluate(request):
    /// ...     requests.append(request)
    /// ...     return hep.EvaluatedValues(
    /// ...         couplings={"SCALAR_COUPLING": (0.0, -1.0)},
    /// ...     )
    /// >>> updated_model = model.recompute_with(evaluate)
    /// >>> request = requests[0]
    /// >>> formulas = {item.name: item.expression for item in request.couplings}
    ///
    /// Parameters
    /// ----------
    /// evaluator : Callable[[EvaluationRequest], EvaluatedValues]
    ///     Callback that evaluates every expression in a request.
    fn recompute_with(
        &self,
        #[gen_stub(override_type(
            type_repr = "collections.abc.Callable[[EvaluationRequest], EvaluatedValues]",
            imports = ("collections.abc")
        ))]
        evaluator: &Bound<'_, PyAny>,
    ) -> PyResult<Self> {
        let mut model = self.inner.as_ref().clone();
        let mut evaluator = PythonModelEvaluator {
            callable: evaluator,
            model: Arc::clone(&self.inner),
        };
        model
            .recompute_with(&mut evaluator)
            .map_err(error::recompute)?;
        Ok(model.into())
    }

    /// Serialize the model as JSON.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Model`` class example:
    ///
    /// >>> model.to_json(pretty=False)
    ///
    /// Parameters
    /// ----------
    /// pretty : bool
    ///     Indent the output when true.
    #[pyo3(signature = (pretty=true))]
    fn to_json(&self, pretty: bool) -> PyResult<String> {
        if pretty {
            self.inner.to_json_pretty()
        } else {
            self.inner.to_json()
        }
        .map_err(error::model)
    }

    /// Write the model to a JSON file.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Model`` class example:
    ///
    /// >>> model.write_json("model.json")
    ///
    /// Parameters
    /// ----------
    /// path : str or os.PathLike
    ///     Destination for the JSON model.
    fn write_json(&self, py: Python<'_>, path: PathBuf) -> PyResult<()> {
        let model = self.inner.clone();
        py.detach(move || model.write_json(path))
            .map_err(error::model)
    }

    /// Summarize the defining data as well as the name and particle count.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Model`` class example:
    ///
    /// >>> print(model)
    fn __repr__(&self) -> String {
        format!(
            "Model(name='{}', particles={})",
            self.inner.name(),
            self.inner.particles().len()
        )
    }

    /// Render a compact inventory of the model in notebook frontends.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Model`` class example:
    ///
    /// >>> from IPython.display import display
    /// >>> display(model)
    fn _repr_html_(&self) -> String {
        let restriction = self
            .inner
            .restriction()
            .map_or_else(|| "none".to_owned(), escape_html);
        format!(
            "<div class=\"feynkit-model\" style=\"display:inline-block;max-width:100%;\
             overflow-x:auto\"><strong>{}</strong><span style=\"margin-left:.45rem;\
             opacity:.7\">restriction: {restriction}</span><table style=\"border-collapse:\
             collapse;margin-top:.35rem\"><thead><tr><th style=\"padding:.2rem .65rem;\
             text-align:right\">particles</th><th style=\"padding:.2rem .65rem;text-align:\
             right\">parameters</th><th style=\"padding:.2rem .65rem;text-align:right\">\
             couplings</th><th style=\"padding:.2rem .65rem;text-align:right\">vertex rules\
             </th><th style=\"padding:.2rem .65rem;text-align:right\">Lorentz structures\
             </th></tr></thead><tbody><tr><td style=\"padding:.2rem .65rem;text-align:right\">\
             {}</td><td style=\"padding:.2rem .65rem;text-align:right\">{}</td><td style=\"\
             padding:.2rem .65rem;text-align:right\">{}</td><td style=\"padding:.2rem .65rem;\
             text-align:right\">{}</td><td style=\"padding:.2rem .65rem;text-align:right\">\
             {}</td></tr></tbody></table></div>",
            escape_html(self.inner.name()),
            self.inner.particles().len(),
            self.inner.parameters().len(),
            self.inner.couplings().len(),
            self.inner.vertex_rules().len(),
            self.inner.lorentz_structures().len(),
        )
    }

    /// Write the concise model summary to an IPython pretty printer.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``Model`` class example:
    ///
    /// >>> from IPython.lib.pretty import pretty
    /// >>> text = pretty(model)
    ///
    /// Parameters
    /// ----------
    /// pretty : object
    ///     The IPython pretty-printer object.
    /// cycle : bool
    ///     Whether this object is part of a recursive formatting cycle.
    fn _repr_pretty_(&self, pretty: &Bound<'_, PyAny>, cycle: bool) -> PyResult<()> {
        pretty.call_method1(
            "text",
            (if cycle {
                "...".to_owned()
            } else {
                self.__repr__()
            },),
        )?;
        Ok(())
    }
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PyParticle>()?;
    module.add_class::<PyParameterNature>()?;
    module.add_class::<PyParameterType>()?;
    module.add_class::<PyParameter>()?;
    module.add_class::<PyCoupling>()?;
    module.add_class::<PyVertexRule>()?;
    module.add_class::<PyLorentzStructure>()?;
    module.add_class::<PyPropagator>()?;
    module.add_class::<PyModelFunction>()?;
    module.add_class::<PyFormFactor>()?;
    module.add_class::<PyModelExpression>()?;
    module.add_class::<PyEvaluationRequest>()?;
    module.add_class::<PyEvaluatedValues>()?;
    module.add_class::<PyParameterCard>()?;
    module.add_class::<PyModel>()?;
    Ok(())
}

#[cfg(feature = "python_stubgen")]
submit! {
    gen_methods_from_python! {
        r#"
        import typing

        class PyModel:
            @typing.overload
            def expand_couplings(
                self,
                expression: pyo3_stub_gen.RustType["TensorExpression"],
            ) -> pyo3_stub_gen.RustType["TensorExpression"]:
                """
                Replace named UFO coefficients by their analytic model expressions.

                The input and the stored model are unchanged.

                Examples
                --------
                Using the setup in the ``Model`` class example:

                >>> analytic = model.expand_couplings(S("UFO::GC_11"))

                Parameters
                ----------
                expression : Expression
                    Symbolica expression containing named couplings from this model.
                """

            @typing.overload
            def expand_couplings(
                self,
                expression: pyo3_stub_gen.RustType["PythonExpression"],
            ) -> pyo3_stub_gen.RustType["PythonExpression"]:
                """
                Replace named UFO coefficients by their analytic model expressions.

                The input and the stored model are unchanged.

                Examples
                --------
                Using the setup in the ``Model`` class example:

                >>> analytic = model.expand_couplings(S("UFO::GC_11"))

                Parameters
                ----------
                expression : Expression
                    Symbolica expression containing named couplings from this model.
                """
        "#
    }
}

#[cfg(test)]
mod tests {
    use std::ffi::CString;

    use pyo3::types::PyDict;

    use super::*;

    const MODEL_JSON: &str = r#"{
        "name": "scalar",
        "restriction": null,
        "orders": [{"name":"QED","expansion_order":99,"hierarchy":1}],
        "parameters": [{
            "name":"ZERO","lhablock":null,"lhacode":null,"nature":"internal",
            "parameter_type":"real","value":[0.0,0.0],"expression":null
        },{
            "name":"mass","lhablock":"MASS","lhacode":[1],"nature":"external",
            "parameter_type":"real","value":[1.0,0.0],"expression":null
        },{
            "name":"double_mass","lhablock":null,"lhacode":null,"nature":"internal",
            "parameter_type":"complex","value":[2.0,0.0],"expression":"2*mass"
        }],
        "particles": [{
            "pdg_code":1,"name":"s","antiname":"s","spin":1,"color":1,
            "mass":"mass","width":"ZERO","texname":"s","antitexname":"s",
            "charge":0.0,"ghost_number":0,"lepton_number":0,"y_charge":0,
            "propagator":"s_prop"
        },{
            "pdg_code":11,"name":"e-","antiname":"e+","spin":2,"color":1,
            "mass":"mass","width":"ZERO","texname":"e^-","antitexname":"e^+",
            "charge":-1.0,"ghost_number":0,"lepton_number":1,"y_charge":-1
        },{
            "pdg_code":-11,"name":"e+","antiname":"e-","spin":2,"color":1,
            "mass":"mass","width":"ZERO","texname":"e^+","antitexname":"e^-",
            "charge":1.0,"ghost_number":0,"lepton_number":-1,"y_charge":1
        }],
        "propagators": [{
            "name":"s_prop","particle":"s","numerator":"1","denominator":"P^2-mass^2"
        }],
        "lorentz_structures": [{"name":"L1","spins":[1,1,1],"structure":"1"}],
        "couplings": [{
            "name":"GC1","expression":"double_mass","orders":[["QED",1]],"value":[2.0,0.0]
        }],
        "vertex_rules": [{
            "name":"V1","particles":["s","s","s"],"color_structures":["1"],
            "lorentz_structures":["L1"],"couplings":[["GC1"]]
        }],
        "functions":[{"name":"twice","arguments":["x"],"expression":"2*x"}],
        "form_factors":[{"name":"FF1","type":"complex","value":"twice(P(1)^2)"}]
    }"#;

    fn registered_module<'py>(py: Python<'py>) -> Bound<'py, PyModule> {
        let module = PyModule::new(py, "symbolica.community.feynkit").unwrap();
        crate::initialize_feynkit(&module).unwrap();
        module
    }

    #[test]
    fn model_constructor_loads_paths_and_reports_read_errors() {
        let directory = tempfile::tempdir().unwrap();
        let model_path = directory.path().join("model.json");
        let missing_path = directory.path().join("missing.json");
        std::fs::write(&model_path, MODEL_JSON).unwrap();

        Python::initialize();
        Python::attach(|py| {
            let module = registered_module(py);
            let locals = PyDict::new(py);
            locals.set_item("fk", &module).unwrap();
            locals
                .set_item("MODEL_PATH", model_path.to_string_lossy().as_ref())
                .unwrap();
            locals
                .set_item("MISSING_PATH", missing_path.to_string_lossy().as_ref())
                .unwrap();
            let code = CString::new(
                r#"
from pathlib import Path

model = fk.Model(Path(MODEL_PATH))
assert model.name == "scalar"
assert model.particle("e-").antiparticle.name == "e+"
assert not hasattr(fk.Model, "from_path")

try:
    fk.Model(MISSING_PATH)
except fk.ModelError as error:
    assert "missing.json" in str(error)
else:
    raise AssertionError("loading a missing model path must raise ModelError")
"#,
            )
            .unwrap();
            py.run(&code, Some(&locals), Some(&locals)).unwrap();
        });
    }

    #[test]
    fn exposes_typed_model_entities_and_atomic_python_evaluation() {
        Python::initialize();
        Python::attach(|py| {
            let module = registered_module(py);
            let locals = PyDict::new(py);
            locals.set_item("fk", &module).unwrap();
            locals.set_item("MODEL_JSON", MODEL_JSON).unwrap();
            let code = CString::new(
                r#"
model = fk.Model.from_json(MODEL_JSON)

def plain(expression):
    return expression.format(
        max_line_length=None,
        color_top_level_sum=False,
        color_builtin_symbols=False,
        bracket_level_colors=None,
        multiplication_operator="*",
        num_exp_as_superscript=False,
    )

electron = model.particle_by_pdg(11)
positron = electron.antiparticle
assert isinstance(positron, fk.Particle)
assert positron.name == "e+"
assert positron.pdg_code == -11
assert str(positron.charge) == "1"
assert positron.antiparticle.name == "e-"

scalar = model.particle("s")
assert scalar.is_self_antiparticle
assert scalar.antiparticle.name == scalar.name
assert scalar.antiparticle.pdg_code == scalar.pdg_code

# A particle keeps enough native model context to resolve its antiparticle
# after the Python Model object has gone out of scope.
detached_positron = fk.Model.from_json(MODEL_JSON).particle("e-").antiparticle
assert detached_positron.name == "e+"
assert detached_positron.antiparticle.name == "e-"

assert len(model.parameters) == 3
mass = model.parameter("mass")
assert isinstance(mass, fk.Parameter)
assert mass.name == "mass"
assert mass.lhablock == "MASS"
assert mass.lhacode == [1]
assert mass.nature == fk.ParameterNature.EXTERNAL
assert mass.parameter_type == fk.ParameterType.REAL
assert mass.value == 1.0 + 0.0j
try:
    mass.name = "changed"
except AttributeError:
    pass
else:
    raise AssertionError("model entity wrappers must be immutable")

internal = model.parameter("double_mass")
assert internal.nature == fk.ParameterNature.INTERNAL
assert internal.parameter_type == fk.ParameterType.COMPLEX
assert plain(internal.expression) == "2*mass"

coupling = model.coupling("GC1")
assert isinstance(coupling, fk.Coupling)
assert coupling.orders == {"QED": 1}
assert coupling.value == 2.0 + 0.0j

vertex = model.vertex_rule("V1")
assert isinstance(vertex, fk.VertexRule)
assert vertex.particles == ["s", "s", "s"]
assert [str(structure) for structure in vertex.color_structures] == ["1"]
assert vertex.lorentz_structures == ["L1"]
assert vertex.couplings == [["GC1"]]
assert vertex.coupling_orders() == {"QED": 1}

lorentz = model.lorentz_structure("L1")
assert isinstance(lorentz, fk.LorentzStructure)
assert lorentz.spins == [1, 1, 1]
assert plain(lorentz.structure) == "1"

propagator = model.propagator("s_prop")
assert isinstance(propagator, fk.Propagator)
assert propagator.particle == "s"
assert plain(propagator.numerator) == "1"
denominator = str(propagator.denominator)
assert "P^2" in denominator and "mass^2" in denominator

function = model.function("twice")
assert isinstance(function, fk.ModelFunction)
assert function.arguments == ["x"]
assert plain(function.expression) == "2*x"

form_factor = model.form_factor("FF1")
assert isinstance(form_factor, fk.FormFactor)
assert form_factor.type_name == "complex"
assert plain(form_factor.value) == "twice(P(1)^2)"

assert all(isinstance(value, fk.Parameter) for value in model.parameters)
assert all(isinstance(value, fk.Coupling) for value in model.couplings)
assert all(isinstance(value, fk.VertexRule) for value in model.vertex_rules)
assert all(isinstance(value, fk.LorentzStructure) for value in model.lorentz_structures)
assert all(isinstance(value, fk.Propagator) for value in model.propagators)
assert all(isinstance(value, fk.ModelFunction) for value in model.functions)
assert all(isinstance(value, fk.FormFactor) for value in model.form_factors)

requests = []
def evaluate(request):
    assert isinstance(request, fk.EvaluationRequest)
    assert all(isinstance(value, fk.ModelExpression) for value in request.internal_parameters)
    assert all(isinstance(value, fk.ModelExpression) for value in request.couplings)
    assert request.internal_parameters[0].name == "double_mass"
    assert plain(request.internal_parameters[0].expression) == "2*mass"
    assert request.couplings[0].name == "GC1"
    assert request.functions[0].name == "twice"
    assert request.form_factors[0].name == "FF1"
    requests.append(request)
    real, imaginary = request.known_parameters["mass"]
    value = (2.0 * real, 2.0 * imaginary)
    return fk.EvaluatedValues(
        internal_parameters={"double_mass": value},
        couplings={"GC1": value},
    )

recomputed = model.recompute_with(evaluate)
assert isinstance(recomputed, fk.Model)
assert recomputed.parameter("double_mass").value == 2.0 + 0.0j
assert recomputed.coupling("GC1").value == 2.0 + 0.0j
assert model.parameter("double_mass").value == 2.0 + 0.0j

card = fk.ParameterCard()
card.set("mass", 3.0, 0.5)
updated = model.with_parameter_card(card, evaluator=evaluate)
assert updated.parameter("mass").value == 3.0 + 0.5j
assert updated.parameter("double_mass").value == 6.0 + 1.0j
assert updated.coupling("GC1").value == 6.0 + 1.0j
assert requests[-1].known_parameters["mass"] == (3.0, 0.5)

invalidated = model.with_parameter_card(card)
assert invalidated.parameter("mass").value == 3.0 + 0.5j
try:
    invalidated.parameter("double_mass").value
except fk.ModelError:
    pass
else:
    raise AssertionError("invalidated parameters must reject value access")
try:
    invalidated.coupling("GC1").value
except fk.ModelError:
    pass
else:
    raise AssertionError("invalidated couplings must reject value access")

before = model.to_json(pretty=False)
def fail(_request):
    raise ValueError("evaluator failed")

try:
    model.with_parameter_card(card, evaluator=fail)
except ValueError as error:
    assert str(error) == "evaluator failed"
else:
    raise AssertionError("the evaluator exception must be preserved")
assert model.to_json(pretty=False) == before

try:
    model.recompute_with(lambda _request: fk.EvaluatedValues())
except fk.ModelError as error:
    assert "did not return a value" in str(error)
else:
    raise AssertionError("incomplete evaluated values must be rejected")
assert model.to_json(pretty=False) == before

try:
    model.recompute_with(lambda _request: {})
except TypeError:
    pass
else:
    raise AssertionError("the evaluator must return EvaluatedValues")
assert model.to_json(pretty=False) == before
"#,
            )
            .unwrap();
            py.run(&code, Some(&locals), Some(&locals)).unwrap();
        });
    }
}
