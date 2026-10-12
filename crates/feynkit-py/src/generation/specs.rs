//! Native-backed specifications for each explicit generation stage.

use std::collections::BTreeMap;

use feynkit_amplitude::AmplitudeOptions;
use feynkit_generator::{
    CutAmplitude, Cuts, GraphGroupingOptions, Grouping, Loops, Naming, NumeratorOptions, Orders,
    Range, Restrict, Selection, Symmetrize, TopologyRules, Veto,
};
use feynkit_graph::EdgeId;
use linnet_py::PyBridges;
use pyo3::{exceptions::PyValueError, prelude::*, types::PyModule};
use symbolica::api::python::{ConvertibleToExpression, PythonExpression};

#[cfg(feature = "python_stubgen")]
use pyo3_stub_gen::derive::{gen_stub_pyclass, gen_stub_pymethods};

use super::{DiagramSelectionInput, OrderRangeInput, SelectorInput, VertexInput};
use crate::error;

/// An inclusive nonnegative count interval. Unrestricted ranges have no upper bound; conflicting bounds retain an empty interval.
///
/// Examples
/// --------
/// >>> from symbolica.community import hepkit as hep
/// >>> bounds = hep.Range.between(1, 3)
/// >>> assert 2 in bounds and 4 not in bounds
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    name = "Range",
    module = "symbolica.community.hepkit",
    frozen,
    from_py_object
)]
#[derive(Clone, Copy, Default)]
pub struct PyRange {
    pub(super) inner: Range,
}

impl PyRange {
    fn integer(value: &Bound<'_, PyAny>) -> PyResult<usize> {
        if value.is_instance_of::<pyo3::types::PyBool>() {
            return Err(pyo3::exceptions::PyTypeError::new_err(
                "a non-negative integer is required, not bool",
            ));
        }
        value.extract()
    }
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[spenso_macros::track_usage(crate::record_usage)]
#[pymethods]
impl PyRange {
    /// Allow every nonnegative count.
    ///
    /// Examples
    /// --------
    /// >>> bounds = hep.Range.unrestricted()
    /// >>> assert 100 in bounds
    #[staticmethod]
    fn unrestricted() -> Self {
        Self::default()
    }

    /// Accept exactly one count.
    ///
    /// Examples
    /// --------
    /// >>> bounds = hep.Range.exact(2)
    /// >>> assert 2 in bounds and 3 not in bounds
    ///
    /// Parameters
    /// ----------
    /// value : int
    ///     Required nonnegative count.
    #[staticmethod]
    fn exact(#[pyo3(from_py_with=PyRange::integer)] value: usize) -> Self {
        Self {
            inner: Range::exact(value),
        }
    }

    /// Bound the count from above, inclusively.
    ///
    /// Examples
    /// --------
    /// >>> bounds = hep.Range.at_most(2)
    /// >>> assert 0 in bounds and 3 not in bounds
    ///
    /// Parameters
    /// ----------
    /// maximum : int
    ///     Largest accepted count.
    #[staticmethod]
    fn at_most(#[pyo3(from_py_with=PyRange::integer)] maximum: usize) -> Self {
        Self {
            inner: Range::at_most(maximum),
        }
    }

    /// Bound the count from below, inclusively.
    ///
    /// Examples
    /// --------
    /// >>> bounds = hep.Range.at_least(2)
    /// >>> assert 2 in bounds and 1 not in bounds
    ///
    /// Parameters
    /// ----------
    /// minimum : int
    ///     Smallest accepted count.
    #[staticmethod]
    fn at_least(#[pyo3(from_py_with=PyRange::integer)] minimum: usize) -> Self {
        Self {
            inner: Range::at_least(minimum),
        }
    }

    /// Set both inclusive bounds, preserving contradictory bounds as empty.
    ///
    /// Examples
    /// --------
    /// >>> bounds = hep.Range.between(3, 1)
    /// >>> assert bounds.is_empty()
    ///
    /// Parameters
    /// ----------
    /// minimum : int
    ///     Smallest accepted count.
    /// maximum : int
    ///     Largest accepted count.
    #[staticmethod]
    fn between(
        #[pyo3(from_py_with=PyRange::integer)] minimum: usize,
        #[pyo3(from_py_with=PyRange::integer)] maximum: usize,
    ) -> Self {
        Self {
            inner: Range::Inclusive {
                min: minimum,
                max: Some(maximum),
            },
        }
    }

    /// Smallest count allowed by the lower bound.
    ///
    /// Examples
    /// --------
    /// >>> assert hep.Range.at_least(2).minimum == 2
    #[getter]
    fn minimum(&self) -> usize {
        self.inner.bounds().0
    }

    /// Inclusive upper bound, or None for an unbounded interval.
    ///
    /// Examples
    /// --------
    /// >>> assert hep.Range.unrestricted().maximum is None
    /// >>> assert hep.Range.exact(2).maximum == 2
    #[getter]
    fn maximum(&self) -> Option<usize> {
        self.inner.bounds().1
    }

    /// Report whether the bounds exclude every count.
    ///
    /// Examples
    /// --------
    /// >>> assert hep.Range.between(3, 1).is_empty()
    fn is_empty(&self) -> bool {
        self.inner.is_empty()
    }

    /// Return the common interval without changing either operand.
    ///
    /// Examples
    /// --------
    /// >>> bounds = hep.Range.at_least(2).intersect(hep.Range.at_most(4))
    /// >>> assert bounds.minimum == 2 and bounds.maximum == 4
    ///
    /// Parameters
    /// ----------
    /// other : Range
    ///     Interval to intersect.
    fn intersect(&self, other: &Self) -> Self {
        Self {
            inner: self.inner.intersect(other.inner),
        }
    }

    /// Test membership in the inclusive interval.
    ///
    /// Examples
    /// --------
    /// >>> assert 2 in hep.Range.between(1, 3)
    ///
    /// Parameters
    /// ----------
    /// value : int
    ///     Nonnegative count to test.
    fn __contains__(&self, #[pyo3(from_py_with=PyRange::integer)] value: usize) -> bool {
        self.inner.contains(value)
    }

    /// Return a compact representation of the configured options.
    ///
    /// Examples
    /// --------
    /// >>> assert "Range" in repr(hep.Range.exact(2))
    fn __repr__(&self) -> String {
        format!("Range({:?})", self.inner)
    }
}

/// A finite loop interval with an explicit amplitude or forward-graph domain.
///
/// Examples
/// --------
/// >>> from symbolica.community import hepkit as hep
/// >>> loops = hep.Loops.total((0, 2))
/// >>> assert loops.minimum == 0 and loops.maximum == 2
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    name = "Loops",
    module = "symbolica.community.hepkit",
    frozen,
    from_py_object
)]
#[derive(Clone)]
pub struct PyLoops {
    pub(super) inner: Loops,
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[spenso_macros::track_usage(crate::record_usage)]
#[pymethods]
impl PyLoops {
    /// Specify total amplitude loops, including vacuum amplitudes.
    ///
    /// Examples
    /// --------
    /// >>> loops = hep.Loops.total((0, 2))
    /// >>> assert loops.minimum == 0 and loops.maximum == 2
    ///
    /// Parameters
    /// ----------
    /// range : int | tuple[int, int]
    ///     Exact loop count or finite inclusive bounds.
    #[staticmethod]
    fn total(range: OrderRangeInput) -> Self {
        Self {
            inner: Loops::total(range.minimum..=range.maximum.expect("finite loop input")),
        }
    }

    /// Specify loops of the sewn forward graph for a cross section.
    ///
    /// Examples
    /// --------
    /// >>> loops = hep.Loops.forward((0, 2))
    /// >>> assert loops.minimum == 0 and loops.maximum == 2
    ///
    /// Parameters
    /// ----------
    /// range : int | tuple[int, int]
    ///     Exact loop count or finite inclusive bounds.
    #[staticmethod]
    fn forward(range: OrderRangeInput) -> Self {
        Self {
            inner: Loops::forward(range.minimum..=range.maximum.expect("finite loop input")),
        }
    }

    /// Inclusive minimum loop count.
    ///
    /// Examples
    /// --------
    /// >>> assert hep.Loops.total((0, 2)).minimum == 0
    #[getter]
    fn minimum(&self) -> usize {
        *self.inner.range().start()
    }

    /// Inclusive maximum loop count.
    ///
    /// Examples
    /// --------
    /// >>> assert hep.Loops.total((0, 2)).maximum == 2
    #[getter]
    fn maximum(&self) -> usize {
        *self.inner.range().end()
    }

    /// Return a compact representation of the configured options.
    ///
    /// Examples
    /// --------
    /// >>> assert "Loops" in repr(hep.Loops.total(0))
    fn __repr__(&self) -> String {
        format!("Loops({:?})", self.inner)
    }
}

/// Immutable coupling-order constraints. Plain and squared bounds intersect; contradictory intervals remain empty.
///
/// Examples
/// --------
/// >>> from symbolica.community import hepkit as hep
/// >>> orders = hep.Orders().min("QED", 3).max("QED", 2)
/// >>> assert orders.is_empty()
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    name = "Orders",
    module = "symbolica.community.hepkit",
    frozen,
    from_py_object
)]
#[derive(Clone, Default)]
pub struct PyOrders {
    pub(super) inner: Orders,
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[spenso_macros::track_usage(crate::record_usage)]
#[pymethods]
impl PyOrders {
    /// Construct coupling-order bounds for discovery and resolved interference.
    ///
    /// Examples
    /// --------
    /// >>> orders = hep.Orders(plain={"QED": hep.Range.exact(2)})
    /// >>> assert not orders.is_empty()
    ///
    /// Parameters
    /// ----------
    /// plain : dict[str, Range] | None
    ///     Bounds on coupling powers of each amplitude or generated diagram.
    /// squared : dict[str, Range] | None
    ///     Bounds on summed coupling powers of interfering amplitude pairs.
    /// perturbative : dict[str, int] | None
    ///     Perturbative expansion powers for cross-section cut-content resolution.
    #[new]
    #[pyo3(signature = (*, plain=None, squared=None, perturbative=None))]
    fn new(
        plain: Option<BTreeMap<String, PyRange>>,
        squared: Option<BTreeMap<String, PyRange>>,
        perturbative: Option<BTreeMap<String, usize>>,
    ) -> Self {
        Self {
            inner: Orders {
                plain: plain
                    .unwrap_or_default()
                    .into_iter()
                    .map(|(name, range)| (name, range.inner))
                    .collect(),
                squared: squared
                    .unwrap_or_default()
                    .into_iter()
                    .map(|(name, range)| (name, range.inner))
                    .collect(),
                perturbative: perturbative.unwrap_or_default(),
            },
        }
    }

    /// Return a specification with an inclusive upper bound intersected into the named coupling constraint.
    ///
    /// Examples
    /// --------
    /// >>> orders = hep.Orders().max("QED", 2)
    /// >>> assert orders.min("QED", 3).max("QED", 2).is_empty()
    ///
    /// Parameters
    /// ----------
    /// name : str
    ///     Coupling-order label from the model.
    /// maximum : int
    ///     Nonnegative coupling power.
    fn max(&self, name: String, #[pyo3(from_py_with=PyRange::integer)] maximum: usize) -> Self {
        Self {
            inner: self.inner.clone().max(name, maximum),
        }
    }

    /// Return a specification with an inclusive lower bound intersected into the named coupling constraint.
    ///
    /// Examples
    /// --------
    /// >>> orders = hep.Orders().min("QED", 2)
    /// >>> assert orders.min("QED", 3).max("QED", 2).is_empty()
    ///
    /// Parameters
    /// ----------
    /// name : str
    ///     Coupling-order label from the model.
    /// minimum : int
    ///     Nonnegative coupling power.
    fn min(&self, name: String, #[pyo3(from_py_with=PyRange::integer)] minimum: usize) -> Self {
        Self {
            inner: self.inner.clone().min(name, minimum),
        }
    }

    /// Return a specification with an exact power intersected into the named coupling constraint.
    ///
    /// Examples
    /// --------
    /// >>> orders = hep.Orders().exact("QED", 2)
    /// >>> assert orders.min("QED", 3).max("QED", 2).is_empty()
    ///
    /// Parameters
    /// ----------
    /// name : str
    ///     Coupling-order label from the model.
    /// value : int
    ///     Nonnegative coupling power.
    fn exact(&self, name: String, #[pyo3(from_py_with=PyRange::integer)] value: usize) -> Self {
        Self {
            inner: self.inner.clone().exact(name, value),
        }
    }

    /// Return a specification with an inclusive upper bound on pair powers intersected into the named coupling constraint.
    ///
    /// Examples
    /// --------
    /// >>> orders = hep.Orders().squared_max("QED", 2)
    /// >>> assert orders.squared_min("QED", 3).squared_max("QED", 2).is_empty()
    ///
    /// Parameters
    /// ----------
    /// name : str
    ///     Coupling-order label from the model.
    /// maximum : int
    ///     Nonnegative coupling power.
    fn squared_max(
        &self,
        name: String,
        #[pyo3(from_py_with=PyRange::integer)] maximum: usize,
    ) -> Self {
        Self {
            inner: self.inner.clone().squared_max(name, maximum),
        }
    }

    /// Return a specification with an inclusive lower bound on pair powers intersected into the named coupling constraint.
    ///
    /// Examples
    /// --------
    /// >>> orders = hep.Orders().squared_min("QED", 2)
    /// >>> assert orders.squared_min("QED", 3).squared_max("QED", 2).is_empty()
    ///
    /// Parameters
    /// ----------
    /// name : str
    ///     Coupling-order label from the model.
    /// minimum : int
    ///     Nonnegative coupling power.
    fn squared_min(
        &self,
        name: String,
        #[pyo3(from_py_with=PyRange::integer)] minimum: usize,
    ) -> Self {
        Self {
            inner: self.inner.clone().squared_min(name, minimum),
        }
    }

    /// Return a specification with an exact pair power intersected into the named coupling constraint.
    ///
    /// Examples
    /// --------
    /// >>> orders = hep.Orders().squared_exact("QED", 2)
    /// >>> assert orders.squared_min("QED", 3).squared_max("QED", 2).is_empty()
    ///
    /// Parameters
    /// ----------
    /// name : str
    ///     Coupling-order label from the model.
    /// value : int
    ///     Nonnegative coupling power.
    fn squared_exact(
        &self,
        name: String,
        #[pyo3(from_py_with=PyRange::integer)] value: usize,
    ) -> Self {
        Self {
            inner: self.inner.clone().squared_exact(name, value),
        }
    }

    /// Set a perturbative expansion power for cross-section cut-content resolution. Amplitude requests reject this constraint.
    ///
    /// Examples
    /// --------
    /// >>> orders = hep.Orders().perturbative("QCD", 1)
    /// >>> model = hep.Model.standard_model()
    /// >>> request = model.process([21, 21], [21, 21]).cross_section(loops=hep.Loops.forward(1), orders=orders)
    ///
    /// Parameters
    /// ----------
    /// name : str
    ///     Coupling-order label from the model.
    /// power : int
    ///     Nonnegative perturbative expansion power.
    fn perturbative(
        &self,
        name: String,
        #[pyo3(from_py_with=PyRange::integer)] power: usize,
    ) -> Self {
        Self {
            inner: self.inner.clone().perturbative(name, power),
        }
    }

    /// Intersect ordinary bounds and keep the stricter perturbative power for each shared label.
    ///
    /// Examples
    /// --------
    /// >>> orders = hep.Orders().min("QED", 3).merge(hep.Orders().max("QED", 2))
    /// >>> assert orders.is_empty()
    ///
    /// Parameters
    /// ----------
    /// other : Orders
    ///     Additional order constraints.
    fn merge(&self, other: &Self) -> Self {
        Self {
            inner: self.inner.clone().merge(other.inner.clone()),
        }
    }

    /// Report contradictory plain or squared coupling bounds.
    ///
    /// Examples
    /// --------
    /// >>> assert hep.Orders().min("QED", 3).max("QED", 2).is_empty()
    fn is_empty(&self) -> bool {
        self.inner.is_empty()
    }

    /// Return a compact representation of the configured options.
    ///
    /// Examples
    /// --------
    /// >>> assert "Orders" in repr(hep.Orders())
    fn __repr__(&self) -> String {
        format!("Orders({:?})", self.inner)
    }
}

/// Particle-species and vertex restrictions resolved in the request model. Allow clauses intersect; exclusions unite. An explicit empty allow selects nothing.
///
/// Examples
/// --------
/// >>> from symbolica.community import hepkit as hep
/// >>> model = hep.Model.standard_model()
/// >>> restriction = hep.Restrict().only_particles([11, 22])
/// >>> request = model.process([11, -11], [22, 22]).amplitude(loops=hep.Loops.total(0), restrict=restriction)
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    name = "Restrict",
    module = "symbolica.community.hepkit",
    frozen,
    from_py_object
)]
#[derive(Clone, Default)]
pub struct PyRestrict {
    pub(super) inner: Restrict,
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[spenso_macros::track_usage(crate::record_usage)]
#[pymethods]
impl PyRestrict {
    /// Start with unrestricted particles and vertices.
    ///
    /// Examples
    /// --------
    /// >>> restriction = hep.Restrict()
    /// >>> request = model.process([11, -11], [22, 22]).amplitude(loops=hep.Loops.total(0), restrict=restriction)
    #[new]
    fn new() -> Self {
        Self::default()
    }

    /// Intersect the allowed particle species with this clause. Return a new restriction.
    ///
    /// Examples
    /// --------
    /// >>> restriction = hep.Restrict().only_particles([11, 22])
    /// >>> request = model.process([11, -11], [22, 22]).amplitude(loops=hep.Loops.total(0), restrict=restriction)
    ///
    /// Parameters
    /// ----------
    /// particles : sequence[Particle | ParticleSelector | str | int]
    ///     Model-bound records or selectors to resolve in the request model.
    fn only_particles(&self, particles: Vec<SelectorInput>) -> Self {
        Self {
            inner: self.inner.clone().only_particles(particles),
        }
    }

    /// Exclude the selected particle species. Return a new restriction.
    ///
    /// Examples
    /// --------
    /// >>> restriction = hep.Restrict().exclude_particles([11, 22])
    /// >>> request = model.process([11, -11], [22, 22]).amplitude(loops=hep.Loops.total(0), restrict=restriction)
    ///
    /// Parameters
    /// ----------
    /// particles : sequence[Particle | ParticleSelector | str | int]
    ///     Model-bound records or selectors to resolve in the request model.
    fn exclude_particles(&self, particles: Vec<SelectorInput>) -> Self {
        Self {
            inner: self.inner.clone().exclude_particles(particles),
        }
    }

    /// Intersect the allowed vertex rules with this clause. Return a new restriction.
    ///
    /// Examples
    /// --------
    /// >>> restriction = hep.Restrict().only_vertices([])
    /// >>> request = model.process([11, -11], [22, 22]).amplitude(loops=hep.Loops.total(0), restrict=restriction)
    ///
    /// Parameters
    /// ----------
    /// vertices : sequence[VertexRule | str]
    ///     Model-bound records or selectors to resolve in the request model.
    fn only_vertices(&self, vertices: Vec<VertexInput>) -> Self {
        Self {
            inner: self.inner.clone().only_vertices(vertices),
        }
    }

    /// Exclude the selected vertex rules. Return a new restriction.
    ///
    /// Examples
    /// --------
    /// >>> restriction = hep.Restrict().exclude_vertices([])
    /// >>> request = model.process([11, -11], [22, 22]).amplitude(loops=hep.Loops.total(0), restrict=restriction)
    ///
    /// Parameters
    /// ----------
    /// vertices : sequence[VertexRule | str]
    ///     Model-bound records or selectors to resolve in the request model.
    fn exclude_vertices(&self, vertices: Vec<VertexInput>) -> Self {
        Self {
            inner: self.inner.clone().exclude_vertices(vertices),
        }
    }

    /// Combine restrictions by intersecting allow clauses and uniting exclusions.
    ///
    /// Examples
    /// --------
    /// >>> restriction = hep.Restrict().only_particles([11, 22]).merge(hep.Restrict().exclude_particles([13]))
    /// >>> request = model.process([11, -11], [22, 22]).amplitude(loops=hep.Loops.total(0), restrict=restriction)
    ///
    /// Parameters
    /// ----------
    /// other : Restrict
    ///     Additional restrictions.
    fn merge(&self, other: &Self) -> Self {
        Self {
            inner: self.inner.clone().merge(other.inner.clone()),
        }
    }

    /// Return a compact representation of the configured options.
    ///
    /// Examples
    /// --------
    /// >>> assert "Restrict" in repr(hep.Restrict())
    fn __repr__(&self) -> String {
        format!("Restrict({:?})", self.inner)
    }
}

/// Mass qualifiers for an opt-in topology veto. Both qualifiers are disabled by default.
///
/// Examples
/// --------
/// >>> from symbolica.community import hepkit as hep
/// >>> assert not hep.Veto().enabled()
/// >>> assert hep.Veto.massless().enabled()
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    name = "Veto",
    module = "symbolica.community.hepkit",
    frozen,
    from_py_object
)]
#[derive(Clone, Copy, Default)]
pub struct PyVeto {
    pub(super) inner: Veto,
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[spenso_macros::track_usage(crate::record_usage)]
#[pymethods]
impl PyVeto {
    /// Choose the mass categories to reject.
    ///
    /// Examples
    /// --------
    /// >>> assert hep.Veto(massive=True).enabled()
    ///
    /// Parameters
    /// ----------
    /// massive : bool
    ///     Reject massive instances of this topology.
    /// massless : bool
    ///     Reject massless instances of this topology.
    #[new]
    #[pyo3(signature = (*, massive=false, massless=false))]
    fn new(massive: bool, massless: bool) -> Self {
        Self {
            inner: Veto { massive, massless },
        }
    }

    /// Reject both massive and massless instances.
    ///
    /// Examples
    /// --------
    /// >>> assert hep.Veto.all().enabled()
    #[staticmethod]
    fn all() -> Self {
        Self { inner: Veto::all() }
    }

    /// Reject massive instances only.
    ///
    /// Examples
    /// --------
    /// >>> assert hep.Veto.massive().enabled()
    #[staticmethod]
    fn massive() -> Self {
        Self {
            inner: Veto::massive(),
        }
    }

    /// Reject massless instances only.
    ///
    /// Examples
    /// --------
    /// >>> assert hep.Veto.massless().enabled()
    #[staticmethod]
    fn massless() -> Self {
        Self {
            inner: Veto::massless(),
        }
    }

    /// Report whether either mass category is vetoed.
    ///
    /// Examples
    /// --------
    /// >>> assert not hep.Veto().enabled()
    /// >>> assert hep.Veto.all().enabled()
    fn enabled(&self) -> bool {
        self.inner.enabled()
    }

    /// Return a compact representation of the configured options.
    ///
    /// Examples
    /// --------
    /// >>> assert "Veto" in repr(hep.Veto())
    fn __repr__(&self) -> String {
        format!("Veto({:?})", self.inner)
    }
}

/// Topology constraints with shared native defaults and no process-dependent presets.
///
/// Examples
/// --------
/// >>> from symbolica.community import hepkit as hep
/// >>> rules = hep.TopologyRules(self_loops=False, tadpoles=hep.Veto.all())
/// >>> assert not rules.self_loops and rules.tadpoles.enabled()
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    name = "TopologyRules",
    module = "symbolica.community.hepkit",
    frozen,
    from_py_object
)]
#[derive(Clone, Default)]
pub struct PyTopologyRules {
    pub(super) inner: TopologyRules,
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[spenso_macros::track_usage(crate::record_usage)]
#[pymethods]
impl PyTopologyRules {
    /// Set independent topology constraints for discovery.
    ///
    /// Examples
    /// --------
    /// >>> rules = hep.TopologyRules(max_vertices=4, fermion_loops=hep.Range.exact(0))
    /// >>> assert rules.max_vertices == 4 and 1 not in rules.fermion_loops
    ///
    /// Parameters
    /// ----------
    /// bridges : Bridges | None
    ///     Bridge constraints from symbolica.community.graph; None uses native unrestricted defaults.
    /// self_loops : bool
    ///     Allow graph edges joining a vertex to itself.
    /// zero_flow_edges : bool
    ///     Allow edges carrying identically zero momentum.
    /// max_vertices : int | None
    ///     Optional upper bound on the number of internal vertices.
    /// self_energies : Veto | None
    ///     Mass-qualified self-energy veto; None disables it.
    /// tadpoles : Veto | None
    ///     Mass-qualified tadpole veto; None disables it.
    /// snails : Veto | None
    ///     Mass-qualified snail veto; None disables it.
    /// fermion_loops : Range | None
    ///     Allowed number of closed fermion loops; None is unrestricted.
    /// factorized_loops : Range | None
    ///     Allowed number of factorized loops; None is unrestricted.
    #[new]
    #[pyo3(signature = (*, bridges=None, self_loops=true, zero_flow_edges=false, max_vertices=None, self_energies=None, tadpoles=None, snails=None, fermion_loops=None, factorized_loops=None))]
    #[allow(clippy::too_many_arguments)]
    fn new(
        bridges: Option<PyBridges>,
        self_loops: bool,
        zero_flow_edges: bool,
        max_vertices: Option<usize>,
        self_energies: Option<PyVeto>,
        tadpoles: Option<PyVeto>,
        snails: Option<PyVeto>,
        fermion_loops: Option<PyRange>,
        factorized_loops: Option<PyRange>,
    ) -> Self {
        let mut inner = TopologyRules::default();
        if let Some(bridges) = bridges {
            inner.bridges = bridges.into();
        }
        inner.self_loops = self_loops;
        inner.zero_flow_edges = zero_flow_edges;
        inner.max_vertices = max_vertices;
        inner.self_energies = self_energies.unwrap_or_default().inner;
        inner.tadpoles = tadpoles.unwrap_or_default().inner;
        inner.snails = snails.unwrap_or_default().inner;
        inner.fermion_loops = fermion_loops.unwrap_or_default().inner;
        inner.factorized_loops = factorized_loops.unwrap_or_default().inner;
        Self { inner }
    }

    /// Bridge constraints from symbolica.community.graph; None uses native unrestricted defaults.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import graph
    /// >>> rules = hep.TopologyRules(bridges=graph.Bridges.any())
    /// >>> assert repr(rules.bridges) == repr(graph.Bridges.any())
    #[getter]
    fn bridges(&self) -> PyBridges {
        self.inner.bridges.into()
    }

    /// Allow graph edges joining a vertex to itself.
    ///
    /// Examples
    /// --------
    /// >>> assert not hep.TopologyRules(self_loops=False).self_loops
    #[getter]
    fn self_loops(&self) -> bool {
        self.inner.self_loops
    }

    /// Allow edges carrying identically zero momentum.
    ///
    /// Examples
    /// --------
    /// >>> assert not hep.TopologyRules().zero_flow_edges
    #[getter]
    fn zero_flow_edges(&self) -> bool {
        self.inner.zero_flow_edges
    }

    /// Optional upper bound on the number of internal vertices.
    ///
    /// Examples
    /// --------
    /// >>> assert hep.TopologyRules(max_vertices=4).max_vertices == 4
    #[getter]
    fn max_vertices(&self) -> Option<usize> {
        self.inner.max_vertices
    }

    /// Mass-qualified self-energy veto; None disables it.
    ///
    /// Examples
    /// --------
    /// >>> assert not hep.TopologyRules().self_energies.enabled()
    #[getter]
    fn self_energies(&self) -> PyVeto {
        PyVeto {
            inner: self.inner.self_energies,
        }
    }

    /// Mass-qualified tadpole veto; None disables it.
    ///
    /// Examples
    /// --------
    /// >>> assert hep.TopologyRules(tadpoles=hep.Veto.all()).tadpoles.enabled()
    #[getter]
    fn tadpoles(&self) -> PyVeto {
        PyVeto {
            inner: self.inner.tadpoles,
        }
    }

    /// Mass-qualified snail veto; None disables it.
    ///
    /// Examples
    /// --------
    /// >>> assert not hep.TopologyRules().snails.enabled()
    #[getter]
    fn snails(&self) -> PyVeto {
        PyVeto {
            inner: self.inner.snails,
        }
    }

    /// Allowed number of closed fermion loops; None is unrestricted.
    ///
    /// Examples
    /// --------
    /// >>> assert 2 in hep.TopologyRules().fermion_loops
    #[getter]
    fn fermion_loops(&self) -> PyRange {
        PyRange {
            inner: self.inner.fermion_loops,
        }
    }

    /// Allowed number of factorized loops; None is unrestricted.
    ///
    /// Examples
    /// --------
    /// >>> assert 1 not in hep.TopologyRules(factorized_loops=hep.Range.exact(0)).factorized_loops
    #[getter]
    fn factorized_loops(&self) -> PyRange {
        PyRange {
            inner: self.inner.factorized_loops,
        }
    }

    /// Return a compact representation of the configured options.
    ///
    /// Examples
    /// --------
    /// >>> assert "TopologyRules" in repr(hep.TopologyRules())
    fn __repr__(&self) -> String {
        format!("TopologyRules({:?})", self.inner)
    }
}

/// Order and loop constraints applied to each amplitude on a physical cut.
///
/// Examples
/// --------
/// >>> from symbolica.community import hepkit as hep
/// >>> side = hep.CutAmplitude(loops=hep.Range.exact(0))
/// >>> assert 0 in side.loops and 1 not in side.loops
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    name = "CutAmplitude",
    module = "symbolica.community.hepkit",
    frozen,
    from_py_object
)]
#[derive(Clone, Default)]
pub struct PyCutAmplitude {
    pub(super) inner: CutAmplitude,
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[spenso_macros::track_usage(crate::record_usage)]
#[pymethods]
impl PyCutAmplitude {
    /// Constrain each amplitude separated by a physical cut.
    ///
    /// Examples
    /// --------
    /// >>> side = hep.CutAmplitude(orders=hep.Orders().max("QED", 2), loops=hep.Range.exact(0))
    /// >>> assert not side.orders.is_empty()
    ///
    /// Parameters
    /// ----------
    /// orders : Orders | None
    ///     Coupling-order constraints for each cut amplitude.
    /// loops : Range | None
    ///     Loop interval for each cut amplitude; None is unrestricted.
    #[new]
    #[pyo3(signature = (*, orders=None, loops=None))]
    fn new(orders: Option<PyOrders>, loops: Option<PyRange>) -> Self {
        Self {
            inner: CutAmplitude {
                orders: orders.unwrap_or_default().inner,
                loops: loops.unwrap_or_default().inner,
            },
        }
    }

    /// Coupling-order constraints for each cut amplitude.
    ///
    /// Examples
    /// --------
    /// >>> side = hep.CutAmplitude(orders=hep.Orders().min("QED", 3).max("QED", 2))
    /// >>> assert side.orders.is_empty()
    #[getter]
    fn orders(&self) -> PyOrders {
        PyOrders {
            inner: self.inner.orders.clone(),
        }
    }

    /// Allowed loops in each cut amplitude.
    ///
    /// Examples
    /// --------
    /// >>> assert 1 not in hep.CutAmplitude(loops=hep.Range.exact(0)).loops
    #[getter]
    fn loops(&self) -> PyRange {
        PyRange {
            inner: self.inner.loops,
        }
    }

    /// Return a compact representation of the configured options.
    ///
    /// Examples
    /// --------
    /// >>> assert "CutAmplitude" in repr(hep.CutAmplitude())
    fn __repr__(&self) -> String {
        format!("CutAmplitude({:?})", self.inner)
    }
}

/// Physical cut constraints. Use Selection for forced-cut selection after building.
///
/// Examples
/// --------
/// >>> from symbolica.community import hepkit as hep
/// >>> cuts = hep.Cuts(amplitude=hep.CutAmplitude(loops=hep.Range.exact(0)))
/// >>> assert 1 not in cuts.amplitude.loops
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    name = "Cuts",
    module = "symbolica.community.hepkit",
    frozen,
    from_py_object
)]
#[derive(Clone, Default)]
pub struct PyCuts {
    pub(super) inner: Cuts,
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[spenso_macros::track_usage(crate::record_usage)]
#[pymethods]
impl PyCuts {
    /// Configure physical cut discovery.
    ///
    /// Examples
    /// --------
    /// >>> cuts = hep.Cuts(spectators=hep.Range.exact(0), symmetrize_left_right=True)
    /// >>> assert 1 not in cuts.spectators and cuts.symmetrize_left_right
    ///
    /// Parameters
    /// ----------
    /// amplitude : CutAmplitude | None
    ///     Constraints on each amplitude separated by the cut.
    /// blobs : Range | None
    ///     Allowed number of cut blobs; None is unrestricted.
    /// spectators : Range | None
    ///     Allowed spectator count; None is unrestricted.
    /// sewn_tadpoles : bool
    ///     Allow sewn tadpole contributions.
    /// symmetrize_left_right : bool
    ///     Identify equivalent left-right cut assignments.
    #[new]
    #[pyo3(signature = (*, amplitude=None, blobs=None, spectators=None, sewn_tadpoles=false, symmetrize_left_right=false))]
    fn new(
        amplitude: Option<PyCutAmplitude>,
        blobs: Option<PyRange>,
        spectators: Option<PyRange>,
        sewn_tadpoles: bool,
        symmetrize_left_right: bool,
    ) -> Self {
        Self {
            inner: Cuts {
                amplitude: amplitude.unwrap_or_default().inner,
                blobs: blobs.unwrap_or_default().inner,
                spectators: spectators.unwrap_or_default().inner,
                sewn_tadpoles,
                symmetrize_left_right,
            },
        }
    }

    /// Constraints on each amplitude separated by the cut.
    ///
    /// Examples
    /// --------
    /// >>> assert 1 not in hep.Cuts(amplitude=hep.CutAmplitude(loops=hep.Range.exact(0))).amplitude.loops
    #[getter]
    fn amplitude(&self) -> PyCutAmplitude {
        PyCutAmplitude {
            inner: self.inner.amplitude.clone(),
        }
    }

    /// Allowed number of cut blobs; None is unrestricted.
    ///
    /// Examples
    /// --------
    /// >>> assert 3 in hep.Cuts().blobs
    #[getter]
    fn blobs(&self) -> PyRange {
        PyRange {
            inner: self.inner.blobs,
        }
    }

    /// Allowed spectator count; None is unrestricted.
    ///
    /// Examples
    /// --------
    /// >>> assert 1 not in hep.Cuts(spectators=hep.Range.exact(0)).spectators
    #[getter]
    fn spectators(&self) -> PyRange {
        PyRange {
            inner: self.inner.spectators,
        }
    }

    /// Allow sewn tadpole contributions.
    ///
    /// Examples
    /// --------
    /// >>> assert not hep.Cuts().sewn_tadpoles
    #[getter]
    fn sewn_tadpoles(&self) -> bool {
        self.inner.sewn_tadpoles
    }

    /// Identify equivalent left-right cut assignments.
    ///
    /// Examples
    /// --------
    /// >>> assert hep.Cuts(symmetrize_left_right=True).symmetrize_left_right
    #[getter]
    fn symmetrize_left_right(&self) -> bool {
        self.inner.symmetrize_left_right
    }

    /// Return a compact representation of the configured options.
    ///
    /// Examples
    /// --------
    /// >>> assert "Cuts" in repr(hep.Cuts())
    fn __repr__(&self) -> String {
        format!("Cuts({:?})", self.inner)
    }
}

/// Independent external-state symmetrizations, combined with the | operator.
///
/// Examples
/// --------
/// >>> from symbolica.community import hepkit as hep
/// >>> symmetry = hep.Symmetrize.INITIAL | hep.Symmetrize.FINAL
/// >>> assert symmetry.initial and symmetry.final_state
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    name = "Symmetrize",
    module = "symbolica.community.hepkit",
    frozen,
    from_py_object
)]
#[derive(Clone, Copy, Default)]
pub struct PySymmetrize {
    pub(super) inner: Symmetrize,
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[spenso_macros::track_usage(crate::record_usage)]
#[pymethods]
impl PySymmetrize {
    /// Choose independent external-state symmetries.
    ///
    /// Examples
    /// --------
    /// >>> symmetry = hep.Symmetrize(initial=True)
    /// >>> assert symmetry.initial and not symmetry.final_state
    ///
    /// Parameters
    /// ----------
    /// initial : bool
    ///     Symmetrize incoming external states.
    /// final_state : bool
    ///     Symmetrize outgoing external states.
    /// external_fermions : bool
    ///     Apply external-fermion symmetrization.
    #[new]
    #[pyo3(signature = (*, initial=false, final_state=false, external_fermions=false))]
    fn new(initial: bool, final_state: bool, external_fermions: bool) -> Self {
        Self {
            inner: Symmetrize {
                initial,
                final_state,
                external_fermions,
            },
        }
    }

    /// Enable only the initial symmetry.
    ///
    /// Examples
    /// --------
    /// >>> assert hep.Symmetrize.INITIAL.initial
    #[classattr]
    #[pyo3(name = "INITIAL")]
    fn initial_symmetry() -> PySymmetrize {
        Self {
            inner: Symmetrize::INITIAL,
        }
    }

    /// Enable only the final state symmetry.
    ///
    /// Examples
    /// --------
    /// >>> assert hep.Symmetrize.FINAL.final_state
    #[classattr]
    #[pyo3(name = "FINAL")]
    fn final_symmetry() -> PySymmetrize {
        Self {
            inner: Symmetrize::FINAL,
        }
    }

    /// Enable only the external fermions symmetry.
    ///
    /// Examples
    /// --------
    /// >>> assert hep.Symmetrize.EXTERNAL_FERMIONS.external_fermions
    #[classattr]
    #[pyo3(name = "EXTERNAL_FERMIONS")]
    fn external_fermion_symmetry() -> PySymmetrize {
        Self {
            inner: Symmetrize::EXTERNAL_FERMIONS,
        }
    }

    /// Whether initial symmetrization is enabled.
    ///
    /// Examples
    /// --------
    /// >>> assert hep.Symmetrize(initial=True).initial
    #[getter]
    fn initial(&self) -> bool {
        self.inner.initial
    }

    /// Whether final state symmetrization is enabled.
    ///
    /// Examples
    /// --------
    /// >>> assert hep.Symmetrize(final_state=True).final_state
    #[getter]
    fn final_state(&self) -> bool {
        self.inner.final_state
    }

    /// Whether external fermions symmetrization is enabled.
    ///
    /// Examples
    /// --------
    /// >>> assert hep.Symmetrize(external_fermions=True).external_fermions
    #[getter]
    fn external_fermions(&self) -> bool {
        self.inner.external_fermions
    }

    /// Combine enabled symmetry flags.
    ///
    /// Examples
    /// --------
    /// >>> symmetry = hep.Symmetrize.INITIAL | hep.Symmetrize.FINAL
    /// >>> assert symmetry.initial and symmetry.final_state
    ///
    /// Parameters
    /// ----------
    /// other : Symmetrize
    ///     Additional flags to enable.
    fn __or__(&self, other: &Self) -> Self {
        Self {
            inner: self.inner | other.inner,
        }
    }

    /// Return a compact representation of the configured options.
    ///
    /// Examples
    /// --------
    /// >>> assert "Symmetrize" in repr(hep.Symmetrize())
    fn __repr__(&self) -> String {
        format!("Symmetrize({:?})", self.inner)
    }
}

/// An explicit naming demand for built diagrams. Discovery and default builds do not compute names.
///
/// Examples
/// --------
/// >>> from symbolica.community import hepkit as hep
/// >>> assert repr(hep.Naming.none()) == "Naming.none()"
/// >>> assert "sample" in repr(hep.Naming.named("sample"))
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    name = "Naming",
    module = "symbolica.community.hepkit",
    frozen,
    from_py_object
)]
#[derive(Clone, Default)]
pub struct PyNaming {
    pub(super) inner: Naming,
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[spenso_macros::track_usage(crate::record_usage)]
#[pymethods]
impl PyNaming {
    /// Leave built graph names unrequested.
    ///
    /// Examples
    /// --------
    /// >>> assert repr(hep.Naming.none()) == "Naming.none()"
    #[staticmethod]
    fn none() -> Self {
        Self::default()
    }

    /// Request built diagram names with the supplied prefix.
    ///
    /// Examples
    /// --------
    /// >>> assert "sample" in repr(hep.Naming.named("sample"))
    ///
    /// Parameters
    /// ----------
    /// prefix : str
    ///     Prefix for deterministic generated names.
    #[staticmethod]
    #[pyo3(signature = (prefix="GL"))]
    fn named(prefix: &str) -> Self {
        Self {
            inner: Naming::prefix(prefix),
        }
    }

    /// Return a compact representation of the configured options.
    ///
    /// Examples
    /// --------
    /// >>> assert "Naming" in repr(hep.Naming.named("sample"))
    fn __repr__(&self) -> String {
        match &self.inner {
            Naming::None => "Naming.none()".to_owned(),
            Naming::Named { prefix } => format!("Naming.named({prefix:?})"),
        }
    }
}

/// Options for the explicit numerator-build stage. These do not construct graphs during discovery.
///
/// Examples
/// --------
/// >>> from symbolica.community import hepkit as hep
/// >>> options = hep.NumeratorOptions(naming=hep.Naming.named("sample"))
/// >>> assert options.prefactor == 1 and "sample" in repr(options.naming)
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    name = "NumeratorOptions",
    module = "symbolica.community.hepkit",
    frozen,
    from_py_object
)]
#[derive(Clone, Default)]
pub struct PyNumeratorOptions {
    pub(super) inner: NumeratorOptions,
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[spenso_macros::track_usage(crate::record_usage)]
#[pymethods]
impl PyNumeratorOptions {
    /// Set numerator, naming, color filtering, and loop-momentum-basis demands.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica import S
    /// >>> factor = S("factor")
    /// >>> options = hep.NumeratorOptions(prefactor=factor)
    /// >>> assert options.prefactor == factor
    ///
    /// Parameters
    /// ----------
    /// prefactor : Expression | None
    ///     Global numerator prefactor; None uses one.
    /// projector : Expression | None
    ///     Optional numerator projector.
    /// naming : Naming | None
    ///     Explicit graph naming policy; None requests no names.
    /// filter_zero_color : bool
    ///     Remove diagrams whose constructed color numerator vanishes.
    /// loop_momentum_bases : sequence[tuple[FeynmanDiagram | str, sequence[int]]] | None
    ///     Loop-momentum edge indices keyed by a diagram, stable ID string, or explicit name.
    #[new]
    #[pyo3(signature = (*, prefactor=None, projector=None, naming=None, filter_zero_color=false, loop_momentum_bases=None))]
    fn new(
        prefactor: Option<PythonExpression>,
        projector: Option<PythonExpression>,
        naming: Option<PyNaming>,
        filter_zero_color: bool,
        loop_momentum_bases: Option<Vec<(DiagramSelectionInput, Vec<usize>)>>,
    ) -> PyResult<Self> {
        let mut inner = NumeratorOptions::default();
        if let Some(prefactor) = prefactor {
            inner.prefactor = prefactor.expr;
        }
        inner.projector = projector.map(|value| value.expr);
        inner.naming = naming.unwrap_or_default().inner;
        inner.filter_zero_color = filter_zero_color;
        for (diagram, edges) in loop_momentum_bases.unwrap_or_default() {
            let (mut ids, mut names) = (Vec::new(), Vec::new());
            diagram.split(&mut ids, &mut names)?;
            let edges = edges.into_iter().map(EdgeId).collect::<Vec<_>>();
            inner = if let Some(id) = ids.into_iter().next() {
                inner.with_loop_momentum_basis(id, edges)
            } else {
                inner.with_named_loop_momentum_basis(names.pop().expect("one selector"), edges)
            };
        }
        Ok(Self { inner })
    }

    /// Global numerator prefactor; None uses one.
    ///
    /// Examples
    /// --------
    /// >>> assert hep.NumeratorOptions().prefactor == 1
    #[getter]
    fn prefactor(&self) -> PythonExpression {
        PythonExpression {
            expr: self.inner.prefactor.clone(),
        }
    }

    /// Optional numerator projector.
    ///
    /// Examples
    /// --------
    /// >>> assert hep.NumeratorOptions().projector is None
    #[getter]
    fn projector(&self) -> Option<PythonExpression> {
        self.inner
            .projector
            .clone()
            .map(|expr| PythonExpression { expr })
    }

    /// Explicit graph naming policy; None requests no names.
    ///
    /// Examples
    /// --------
    /// >>> assert repr(hep.NumeratorOptions().naming) == "Naming.none()"
    #[getter]
    fn naming(&self) -> PyNaming {
        PyNaming {
            inner: self.inner.naming.clone(),
        }
    }

    /// Remove diagrams whose constructed color numerator vanishes.
    ///
    /// Examples
    /// --------
    /// >>> assert hep.NumeratorOptions(filter_zero_color=True).filter_zero_color
    #[getter]
    fn filter_zero_color(&self) -> bool {
        self.inner.filter_zero_color
    }

    /// Return a compact representation of the configured options.
    ///
    /// Examples
    /// --------
    /// >>> assert "NumeratorOptions" in repr(hep.NumeratorOptions())
    fn __repr__(&self) -> String {
        format!("NumeratorOptions({:?})", self.inner)
    }
}

/// Selection after numerator construction. None leaves selection unrestricted; an explicit empty selection keeps no diagrams. Mixed stable-ID and name selectors form a union.
///
/// Examples
/// --------
/// >>> from symbolica.community import hepkit as hep
/// >>> selection = hep.Selection(select=[])
/// >>> assert "Selection" in repr(selection)
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    name = "Selection",
    module = "symbolica.community.hepkit",
    frozen,
    from_py_object
)]
#[derive(Clone, Default)]
pub struct PySelection {
    pub(super) inner: Selection,
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[spenso_macros::track_usage(crate::record_usage)]
#[pymethods]
impl PySelection {
    /// Specify diagram inclusion, vetoes, and optional forced cuts.
    ///
    /// Examples
    /// --------
    /// >>> selection = hep.Selection(select=["sample_0"], veto=["sample_1"])
    /// >>> assert "sample_0" in repr(selection)
    ///
    /// Parameters
    /// ----------
    /// select : sequence[FeynmanDiagram | str] | None
    ///     Diagrams, stable IDs, or explicit names to retain; None retains all and [] retains none.
    /// veto : sequence[FeynmanDiagram | str] | None
    ///     Diagrams, stable IDs, or names to exclude.
    /// forced_cuts : sequence[sequence[int]] | None
    ///     Accepted cuts represented by edge indices.
    #[new]
    #[pyo3(signature = (*, select=None, veto=None, forced_cuts=None))]
    fn new(
        select: Option<Vec<DiagramSelectionInput>>,
        veto: Option<Vec<DiagramSelectionInput>>,
        forced_cuts: Option<Vec<Vec<usize>>>,
    ) -> PyResult<Self> {
        let mut inner = Selection::default();
        if let Some(diagrams) = select {
            let (mut ids, mut names) = (Vec::new(), Vec::new());
            for diagram in diagrams {
                diagram.split(&mut ids, &mut names)?;
            }
            // Preserve None versus an explicit empty selection and the native ID/name union.
            if !ids.is_empty() || names.is_empty() {
                inner = inner.select_ids(ids);
            }
            if !names.is_empty() {
                inner = inner.select_names(names);
            }
        }
        if let Some(diagrams) = veto {
            let (mut ids, mut names) = (Vec::new(), Vec::new());
            for diagram in diagrams {
                diagram.split(&mut ids, &mut names)?;
            }
            inner = inner.veto_ids(ids).veto_names(names);
        }
        if let Some(cuts) = forced_cuts {
            inner = inner.with_forced_cuts(cuts.into_iter().map(|cut| cut.into_iter().map(EdgeId)));
        }
        Ok(Self { inner })
    }

    /// Return a compact representation of the configured options.
    ///
    /// Examples
    /// --------
    /// >>> assert "Selection" in repr(hep.Selection())
    fn __repr__(&self) -> String {
        format!("Selection({:?})", self.inner)
    }
}

/// Compare completed numerators exactly, up to a sign, or up to a scalar.
///
/// Examples
/// --------
/// >>> from symbolica.community import hepkit as hep
/// >>> grouping = hep.Grouping("up_to_sign", number_of_numerical_samples=7)
/// >>> assert "UpToSign" in repr(grouping)
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    name = "Grouping",
    module = "symbolica.community.hepkit",
    frozen,
    from_py_object
)]
#[derive(Clone, Default)]
pub struct PyGrouping {
    pub(super) inner: Grouping,
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[spenso_macros::track_usage(crate::record_usage)]
#[pymethods]
impl PyGrouping {
    /// Select the numerator equivalence relation and comparison settings.
    ///
    /// Examples
    /// --------
    /// >>> assert "Identical" in repr(hep.Grouping("identical"))
    ///
    /// Parameters
    /// ----------
    /// mode : str
    ///     One of none, identical, up_to_sign, or up_to_scalar.
    /// numerical_sample_seed : int
    ///     Seed for deterministic numerical comparison samples.
    /// number_of_numerical_samples : int
    ///     Number of comparison samples.
    /// differentiate_particle_masses_only : bool
    ///     Distinguish particle masses while sharing other substitutions.
    /// fully_numerical_substitution : bool
    ///     Substitute all supported numerator variables numerically.
    /// check_canonical_numerator : bool
    ///     Also compare canonicalized numerators.
    /// symmetric_polarizations : bool
    ///     Use symmetric external polarization treatment during comparison.
    #[new]
    #[pyo3(signature = (mode="none", *, numerical_sample_seed=3, number_of_numerical_samples=5, differentiate_particle_masses_only=true, fully_numerical_substitution=false, check_canonical_numerator=false, symmetric_polarizations=false))]
    #[allow(clippy::too_many_arguments)]
    fn new(
        mode: &str,
        numerical_sample_seed: u16,
        number_of_numerical_samples: usize,
        differentiate_particle_masses_only: bool,
        fully_numerical_substitution: bool,
        check_canonical_numerator: bool,
        symmetric_polarizations: bool,
    ) -> PyResult<Self> {
        let options = GraphGroupingOptions {
            numerical_sample_seed,
            number_of_numerical_samples,
            differentiate_particle_masses_only,
            fully_numerical_substitution,
            check_canonical_numerator,
            symmetric_polarizations,
        };
        let inner = match mode {
            "none" => Grouping::None,
            "identical" => Grouping::Identical(options),
            "up_to_sign" => Grouping::UpToSign(options),
            "up_to_scalar" => Grouping::UpToScalar(options),
            _ => {
                return Err(PyValueError::new_err(
                    "grouping mode must be 'none', 'identical', 'up_to_sign', or 'up_to_scalar'",
                ));
            }
        };
        Ok(Self { inner })
    }

    /// Return a compact representation of the configured options.
    ///
    /// Examples
    /// --------
    /// >>> assert "Grouping" in repr(hep.Grouping())
    fn __repr__(&self) -> String {
        format!("Grouping({:?})", self.inner)
    }
}

/// Lorentz dimension and additional reality assumptions for amplitude construction.
///
/// Examples
/// --------
/// >>> from symbolica.community import hepkit as hep
/// >>> from symbolica import S
/// >>> mass = S("mass")
/// >>> options = hep.AmplitudeOptions(dimension=4, real=[mass])
/// >>> assert options.real == [mass]
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    name = "AmplitudeOptions",
    module = "symbolica.community.hepkit",
    frozen,
    from_py_object
)]
#[derive(Clone, Default)]
pub struct PyAmplitudeOptions {
    pub(super) inner: AmplitudeOptions,
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[spenso_macros::track_usage(crate::record_usage)]
#[pymethods]
impl PyAmplitudeOptions {
    /// Set amplitude dimension and symbols treated as real.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica import S
    /// >>> mass = S("mass")
    /// >>> options = hep.AmplitudeOptions(dimension=4, real=[mass])
    /// >>> assert options.real == [mass]
    ///
    /// Parameters
    /// ----------
    /// dimension : int | Expression | None
    ///     Positive integer or symbol used as the Lorentz dimension; None uses the native default.
    /// real : sequence[Expression] | None
    ///     Additional expressions treated as real under conjugation.
    #[new]
    #[pyo3(signature = (*, dimension=None, real=None))]
    fn new(
        dimension: Option<ConvertibleToExpression>,
        real: Option<Vec<PythonExpression>>,
    ) -> PyResult<Self> {
        let mut inner = AmplitudeOptions::default();
        if let Some(dimension) = dimension {
            inner.dimension = dimension
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
        inner.real = real
            .unwrap_or_default()
            .into_iter()
            .map(|value| value.expr)
            .collect();
        Ok(Self { inner })
    }

    /// Additional expressions treated as real.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica import S
    /// >>> mass = S("mass")
    /// >>> assert hep.AmplitudeOptions(real=[mass]).real == [mass]
    #[getter]
    fn real(&self) -> Vec<PythonExpression> {
        self.inner
            .real
            .iter()
            .cloned()
            .map(|expr| PythonExpression { expr })
            .collect()
    }

    /// Return a compact representation of the configured options.
    ///
    /// Examples
    /// --------
    /// >>> assert "AmplitudeOptions" in repr(hep.AmplitudeOptions())
    fn __repr__(&self) -> String {
        format!("AmplitudeOptions({:?})", self.inner)
    }
}

pub(super) fn register(_py: Python<'_>, module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PyBridges>()?;
    module.add_class::<PyRange>()?;
    module.add_class::<PyLoops>()?;
    module.add_class::<PyOrders>()?;
    module.add_class::<PyRestrict>()?;
    module.add_class::<PyVeto>()?;
    module.add_class::<PyTopologyRules>()?;
    module.add_class::<PyCutAmplitude>()?;
    module.add_class::<PyCuts>()?;
    module.add_class::<PySymmetrize>()?;
    module.add_class::<PyNaming>()?;
    module.add_class::<PyNumeratorOptions>()?;
    module.add_class::<PySelection>()?;
    module.add_class::<PyGrouping>()?;
    module.add_class::<PyAmplitudeOptions>()?;
    Ok(())
}
