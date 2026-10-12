//! Exact counting results and the diagnostics produced when compiling a plan.

use linnet::half_edge::generate::{
    Counts, LoopCounts, PlanCriterion, PlanReport, PlanRole, PlanStage,
};
use pyo3::{prelude::*, types::PyDict};

/// Exact totals for a generation run, including whether the search completed.
#[cfg_attr(feature = "python_stubgen", pyo3_stub_gen::derive::gen_stub_pyclass)]
#[pyclass(name = "Counts", frozen)]
#[cfg_attr(feature = "community", pyo3(module = "symbolica.community.graph"))]
#[cfg_attr(not(feature = "community"), pyo3(module = "linnet_py"))]
pub(super) struct PyCounts(pub(super) Counts);

impl PyCounts {
    // A power-of-two radix avoids Python's decimal digit limit for exact large totals.
    fn integer(py: Python<'_>, value: impl std::fmt::LowerHex) -> PyResult<Bound<'_, PyAny>> {
        py.import("builtins")?
            .getattr("int")?
            .call1((format!("{value:x}"), 16))
    }

    fn fraction(
        py: Python<'_>,
        numerator: impl std::fmt::LowerHex,
        denominator: impl std::fmt::LowerHex,
    ) -> PyResult<Bound<'_, PyAny>> {
        let numerator = Self::integer(py, numerator)?;
        let denominator = Self::integer(py, denominator)?;
        py.import("fractions")?
            .getattr("Fraction")?
            .call1((numerator, denominator))
    }
}

#[cfg_attr(feature = "python_stubgen", pyo3_stub_gen::derive::gen_stub_pymethods)]
#[cfg_attr(not(feature = "python_stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[pymethods]
impl PyCounts {
    #[getter]
    #[gen_stub(override_return_type(type_repr = "builtins.int", imports = ("builtins")))]
    fn skeletons<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        Self::integer(py, &self.0.skeletons)
    }

    #[getter]
    #[gen_stub(override_return_type(type_repr = "builtins.int", imports = ("builtins")))]
    fn diagrams<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        Self::integer(py, &self.0.diagrams)
    }

    /// The exact sum of inverse automorphism orders, without process-specific factors.
    #[getter]
    #[gen_stub(override_return_type(type_repr = "fractions.Fraction", imports = ("fractions")))]
    fn sum_inv_aut<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        Self::fraction(py, self.0.sum_inv_aut.numer(), self.0.sum_inv_aut.denom())
    }

    #[getter]
    #[gen_stub(override_return_type(type_repr = "dict[int, LoopCounts]"))]
    fn per_loop<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let result = PyDict::new(py);
        for (&loops, counts) in &self.0.per_loop {
            result.set_item(loops, PyLoopCounts(counts.clone()))?;
        }
        Ok(result)
    }

    #[getter]
    fn complete(&self) -> bool {
        self.0.complete
    }
}

/// Exact totals for one loop order.
#[cfg_attr(feature = "python_stubgen", pyo3_stub_gen::derive::gen_stub_pyclass)]
#[pyclass(name = "LoopCounts", frozen)]
#[cfg_attr(feature = "community", pyo3(module = "symbolica.community.graph"))]
#[cfg_attr(not(feature = "community"), pyo3(module = "linnet_py"))]
pub(super) struct PyLoopCounts(pub(super) LoopCounts);

#[cfg_attr(feature = "python_stubgen", pyo3_stub_gen::derive::gen_stub_pymethods)]
#[cfg_attr(not(feature = "python_stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[pymethods]
impl PyLoopCounts {
    #[getter]
    #[gen_stub(override_return_type(type_repr = "builtins.int", imports = ("builtins")))]
    fn skeletons<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        PyCounts::integer(py, &self.0.skeletons)
    }

    #[getter]
    #[gen_stub(override_return_type(type_repr = "builtins.int", imports = ("builtins")))]
    fn diagrams<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        PyCounts::integer(py, &self.0.diagrams)
    }

    #[getter]
    #[gen_stub(override_return_type(type_repr = "fractions.Fraction", imports = ("fractions")))]
    fn sum_inv_aut<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        PyCounts::fraction(py, self.0.sum_inv_aut.numer(), self.0.sum_inv_aut.denom())
    }
}

/// Effective restrictions and propagation decisions from compiling a query.
#[cfg_attr(feature = "python_stubgen", pyo3_stub_gen::derive::gen_stub_pyclass)]
#[pyclass(name = "PlanReport", frozen)]
#[cfg_attr(feature = "community", pyo3(module = "symbolica.community.graph"))]
#[cfg_attr(not(feature = "community"), pyo3(module = "linnet_py"))]
pub(super) struct PyPlanReport(pub(super) PlanReport);

#[cfg_attr(feature = "python_stubgen", pyo3_stub_gen::derive::gen_stub_pymethods)]
#[cfg_attr(not(feature = "python_stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[pymethods]
impl PyPlanReport {
    #[getter]
    fn criteria(&self) -> Vec<PyPlanCriterion> {
        self.0
            .criteria
            .iter()
            .cloned()
            .map(PyPlanCriterion)
            .collect()
    }

    #[getter]
    fn empty_partitions(&self) -> Vec<(usize, String)> {
        self.0.empty_partitions.clone()
    }

    #[getter]
    fn infeasible(&self) -> Option<String> {
        self.0.infeasible.clone()
    }

    fn __str__(&self) -> String {
        self.0.to_string()
    }
}

/// One effective constraint and where it is enforced.
#[cfg_attr(feature = "python_stubgen", pyo3_stub_gen::derive::gen_stub_pyclass)]
#[pyclass(name = "PlanCriterion", frozen)]
#[cfg_attr(feature = "community", pyo3(module = "symbolica.community.graph"))]
#[cfg_attr(not(feature = "community"), pyo3(module = "linnet_py"))]
pub(super) struct PyPlanCriterion(pub(super) PlanCriterion);

#[cfg_attr(feature = "python_stubgen", pyo3_stub_gen::derive::gen_stub_pymethods)]
#[cfg_attr(not(feature = "python_stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[pymethods]
impl PyPlanCriterion {
    #[getter]
    fn criterion(&self) -> &str {
        &self.0.criterion
    }

    #[getter]
    fn effective(&self) -> &str {
        &self.0.effective
    }

    #[getter]
    fn stage(&self) -> &'static str {
        match self.0.stage {
            PlanStage::Planning => "planning",
            PlanStage::Skeletons => "skeletons",
            PlanStage::Dressing => "dressing",
            PlanStage::Rules => "rules",
            PlanStage::Selection => "selection",
        }
    }

    #[getter]
    fn hereditary(&self) -> bool {
        self.0.hereditary
    }

    #[getter]
    fn role(&self) -> &'static str {
        match self.0.role {
            PlanRole::Restriction => "restriction",
            PlanRole::Propagator => "propagator",
            PlanRole::FinalCheck => "final_check",
        }
    }
}

/// Progress at an exact committed-output boundary; returning false stops before further work.
#[cfg_attr(feature = "python_stubgen", pyo3_stub_gen::derive::gen_stub_pyclass)]
#[pyclass(name = "Progress", frozen)]
#[cfg_attr(feature = "community", pyo3(module = "symbolica.community.graph"))]
#[cfg_attr(not(feature = "community"), pyo3(module = "linnet_py"))]
pub(super) struct PyProgress(pub(super) linnet::half_edge::generate::Progress);
#[cfg_attr(feature = "python_stubgen", pyo3_stub_gen::derive::gen_stub_pymethods)]
#[cfg_attr(not(feature = "python_stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[pymethods]
impl PyProgress {
    #[getter]
    fn loops(&self) -> usize {
        self.0.loops
    }
    #[getter]
    #[gen_stub(override_return_type(type_repr = "int"))]
    fn topologies<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        PyCounts::integer(py, &self.0.topologies)
    }
    #[getter]
    #[gen_stub(override_return_type(type_repr = "int"))]
    fn diagrams<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        PyCounts::integer(py, &self.0.diagrams)
    }
}

pub(super) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PyCounts>()?;
    module.add_class::<PyProgress>()?;
    module.add_class::<PyLoopCounts>()?;
    module.add_class::<PyPlanReport>()?;
    module.add_class::<PyPlanCriterion>()?;
    Ok(())
}
