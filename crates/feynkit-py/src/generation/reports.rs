use super::*;
use feynkit_generator::{Counts, GenerationProgress, GenerationReport, LoopCounts, Stage};
pub use linnet_py::PyPlanReport;
use pyo3::types::PyDict;

/// A stage in the physics discovery, construction, selection, or grouping pipeline.
///
/// Examples
/// --------
/// >>> from symbolica.community import hepkit as hep
/// >>> model = hep.Model.standard_model()
/// >>> request = model.process([11, -11], [22, 22]).amplitude(loops=hep.Loops.total(0))
/// >>> events = []
/// >>> run = hep.Run(progress=lambda event: events.append(event) or True)
/// >>> counts = request.count(run=run)
/// >>> assert counts.complete and events
/// >>> assert any(event.stage == hep.Stage.Skeletons for event in events)
#[cfg_attr(
    feature = "python_stubgen",
    pyo3_stub_gen::derive::gen_stub_pyclass_enum
)]
#[pyclass(
    name = "Stage",
    from_py_object,
    module = "symbolica.community.hepkit",
    frozen,
    eq,
    eq_int,
    hash
)]
#[derive(Clone, Copy, PartialEq, Eq, Hash)]
pub enum PyStage {
    Skeletons,
    Dressing,
    Prepare,
    Numerators,
    Selection,
    ZeroColor,
    GroupingPreparation,
    GroupingSamples,
    GroupingComparison,
    Grouping,
    Done,
}
impl From<Stage> for PyStage {
    fn from(stage: Stage) -> Self {
        match stage {
            Stage::Skeletons => Self::Skeletons,
            Stage::Dressing => Self::Dressing,
            Stage::Prepare => Self::Prepare,
            Stage::Numerators => Self::Numerators,
            Stage::Selection => Self::Selection,
            Stage::ZeroColor => Self::ZeroColor,
            Stage::GroupingPreparation => Self::GroupingPreparation,
            Stage::GroupingSamples => Self::GroupingSamples,
            Stage::GroupingComparison => Self::GroupingComparison,
            Stage::Grouping => Self::Grouping,
            Stage::Done => Self::Done,
        }
    }
}

/// Exact discovery totals returned by Request.count(), before numerator construction. Integers have arbitrary precision; sum_inv_aut and physical_weight are exact fractions. Interrupted prefixes have complete=False.
///
/// Examples
/// --------
/// >>> from symbolica.community import hepkit as hep
/// >>> model = hep.Model.standard_model()
/// >>> request = model.process([11, -11], [22, 22]).amplitude(loops=hep.Loops.total(0))
/// >>> counts = request.count()
/// >>> assert counts.complete and counts.diagrams > 0
/// >>> assert counts.sum_inv_aut > 0
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(name = "Counts", module = "symbolica.community.hepkit", frozen)]
pub struct PyCounts(pub(super) Counts);
impl PyCounts {
    // Power-of-two radix bypasses Python's decimal digit limit for exact large counts.
    pub(super) fn integer(
        py: Python<'_>,
        value: impl std::fmt::LowerHex,
    ) -> PyResult<Bound<'_, PyAny>> {
        py.import("builtins")?
            .getattr("int")?
            .call1((format!("{value:x}"), 16))
    }
    fn fraction(
        py: Python<'_>,
        numerator: impl std::fmt::LowerHex,
        denominator: impl std::fmt::LowerHex,
    ) -> PyResult<Bound<'_, PyAny>> {
        py.import("fractions")?.getattr("Fraction")?.call1((
            Self::integer(py, numerator)?,
            Self::integer(py, denominator)?,
        ))
    }
}
/// Exact discovery totals for a single loop order, obtained from Counts.per_loop.
///
/// Examples
/// --------
/// >>> from symbolica.community import hepkit as hep
/// >>> model = hep.Model.standard_model()
/// >>> request = model.process([11, -11], [22, 22]).amplitude(loops=hep.Loops.total(0))
/// >>> counts = request.count()
/// >>> loop_counts = counts.per_loop[0]
/// >>> assert loop_counts.diagrams == counts.diagrams
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(name = "LoopCounts", module = "symbolica.community.hepkit", frozen)]
pub struct PyLoopCounts(LoopCounts);

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "python_stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[pymethods]
impl PyCounts {
    /// Number of complete untyped skeletons accepted by enumeration.
    ///
    /// Examples
    /// --------
    /// Using the setup in this class example:
    ///
    /// >>> assert counts.skeletons > 0
    #[getter]
    #[gen_stub(override_return_type(type_repr="builtins.int",imports=("builtins")))]
    fn skeletons<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        PyCounts::integer(py, &self.0.skeletons)
    }
    /// Number of complete field and interaction assignments visited before physical preparation filters.
    ///
    /// Examples
    /// --------
    /// Using the setup in this class example:
    ///
    /// >>> assert counts.dressings > 0
    #[getter]
    #[gen_stub(override_return_type(type_repr="builtins.int",imports=("builtins")))]
    fn dressings<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        PyCounts::integer(py, &self.0.dressings)
    }
    /// Number of physically prepared assignments retained before numerator construction.
    ///
    /// Examples
    /// --------
    /// Using the setup in this class example:
    ///
    /// >>> assert counts.diagrams > 0
    #[getter]
    #[gen_stub(override_return_type(type_repr="builtins.int",imports=("builtins")))]
    fn diagrams<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        PyCounts::integer(py, &self.0.diagrams)
    }
    /// Exact sum of inverse automorphism orders of retained prepared diagrams.
    ///
    /// Examples
    /// --------
    /// Using the setup in this class example:
    ///
    /// >>> assert counts.sum_inv_aut > 0
    #[getter]
    #[gen_stub(override_return_type(type_repr="fractions.Fraction",imports=("fractions")))]
    fn sum_inv_aut<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        PyCounts::fraction(py, self.0.sum_inv_aut.numer(), self.0.sum_inv_aut.denom())
    }
    /// Exact physical weight including signed process factors; distinct from the unsigned inverse-automorphism sum.
    ///
    /// Examples
    /// --------
    /// Using the setup in this class example:
    ///
    /// >>> weight = counts.physical_weight
    /// >>> assert weight.denominator > 0
    #[getter]
    #[gen_stub(override_return_type(type_repr="fractions.Fraction",imports=("fractions")))]
    fn physical_weight<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        PyCounts::fraction(
            py,
            self.0.physical_weight.numer(),
            self.0.physical_weight.denom(),
        )
    }
    /// Whether enumeration finished without cancellation.
    ///
    /// Examples
    /// --------
    /// Using the setup in this class example:
    ///
    /// >>> assert counts.complete
    #[getter]
    fn complete(&self) -> bool {
        self.0.complete
    }
    /// Loop order to exact totals mapping for the observed search.
    ///
    /// Examples
    /// --------
    /// Using the setup in this class example:
    ///
    /// >>> assert sum(c.diagrams for c in counts.per_loop.values()) == counts.diagrams
    #[getter]
    #[gen_stub(override_return_type(type_repr = "dict[int, LoopCounts]"))]
    fn per_loop<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let result = PyDict::new(py);
        for (&loops, counts) in &self.0.per_loop {
            result.set_item(loops, PyLoopCounts(counts.clone()))?;
        }
        Ok(result)
    }
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "python_stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[pymethods]
impl PyLoopCounts {
    /// Number of complete untyped skeletons accepted by enumeration.
    ///
    /// Examples
    /// --------
    /// Using the setup in this class example:
    ///
    /// >>> assert loop_counts.skeletons > 0
    #[getter]
    #[gen_stub(override_return_type(type_repr="builtins.int",imports=("builtins")))]
    fn skeletons<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        PyCounts::integer(py, &self.0.skeletons)
    }
    /// Number of complete field and interaction assignments visited before physical preparation filters.
    ///
    /// Examples
    /// --------
    /// Using the setup in this class example:
    ///
    /// >>> assert loop_counts.dressings > 0
    #[getter]
    #[gen_stub(override_return_type(type_repr="builtins.int",imports=("builtins")))]
    fn dressings<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        PyCounts::integer(py, &self.0.dressings)
    }
    /// Number of physically prepared assignments retained before numerator construction.
    ///
    /// Examples
    /// --------
    /// Using the setup in this class example:
    ///
    /// >>> assert loop_counts.diagrams > 0
    #[getter]
    #[gen_stub(override_return_type(type_repr="builtins.int",imports=("builtins")))]
    fn diagrams<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        PyCounts::integer(py, &self.0.diagrams)
    }
    /// Exact sum of inverse automorphism orders of retained prepared diagrams.
    ///
    /// Examples
    /// --------
    /// Using the setup in this class example:
    ///
    /// >>> assert loop_counts.sum_inv_aut > 0
    #[getter]
    #[gen_stub(override_return_type(type_repr="fractions.Fraction",imports=("fractions")))]
    fn sum_inv_aut<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        PyCounts::fraction(py, self.0.sum_inv_aut.numer(), self.0.sum_inv_aut.denom())
    }
    /// Exact physical weight including signed process factors; distinct from the unsigned inverse-automorphism sum.
    ///
    /// Examples
    /// --------
    /// Using the setup in this class example:
    ///
    /// >>> weight = loop_counts.physical_weight
    /// >>> assert weight.denominator > 0
    #[getter]
    #[gen_stub(override_return_type(type_repr="fractions.Fraction",imports=("fractions")))]
    fn physical_weight<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        PyCounts::fraction(
            py,
            self.0.physical_weight.numer(),
            self.0.physical_weight.denom(),
        )
    }
}

/// An owned progress snapshot delivered to Run.progress. It remains readable after the callback returns.
///
/// Examples
/// --------
/// >>> from symbolica.community import hepkit as hep
/// >>> model = hep.Model.standard_model()
/// >>> request = model.process([11, -11], [22, 22]).amplitude(loops=hep.Loops.total(0))
/// >>> events = []
/// >>> run = hep.Run(progress=lambda event: events.append(event) or True)
/// >>> counts = request.count(run=run)
/// >>> assert counts.complete and events
/// >>> assert all(event.completed >= 0 for event in events)
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    name = "GenerationProgress",
    module = "symbolica.community.hepkit",
    frozen,
    from_py_object
)]
#[derive(Clone)]
pub struct PyGenerationProgress {
    pub(super) inner: GenerationProgress,
}
#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "python_stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[pymethods]
impl PyGenerationProgress {
    /// Exact number of completed items in this progress stage.
    ///
    /// Examples
    /// --------
    /// Using the setup in this class example:
    ///
    /// >>> assert all(event.completed >= 0 for event in events)
    #[getter]
    #[gen_stub(override_return_type(type_repr="builtins.int",imports=("builtins")))]
    fn completed<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        PyCounts::integer(py, &self.inner.completed)
    }
    /// Exact number of zero numerators reported removed so far.
    ///
    /// Examples
    /// --------
    /// Using the setup in this class example:
    ///
    /// >>> assert all(event.zero_numerator_count == 0 for event in events)
    #[getter]
    #[gen_stub(override_return_type(type_repr="builtins.int",imports=("builtins")))]
    fn zero_numerator_count<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        PyCounts::integer(py, &self.inner.zero_numerator_count)
    }
    /// Stage that produced this snapshot.
    ///
    /// Examples
    /// --------
    /// Using the setup in this class example:
    ///
    /// >>> assert any(event.stage == hep.Stage.Skeletons for event in events)
    #[getter]
    fn stage(&self) -> PyStage {
        self.inner.stage.into()
    }
    /// Expected stage total when known, otherwise None.
    ///
    /// Examples
    /// --------
    /// Using the setup in this class example:
    ///
    /// >>> known_totals = [event.total for event in events if event.total is not None]
    /// >>> assert all(total >= 0 for total in known_totals)
    #[getter]
    #[gen_stub(override_return_type(type_repr = "int | None"))]
    fn total<'py>(&self, py: Python<'py>) -> PyResult<Option<Bound<'py, PyAny>>> {
        self.inner
            .total
            .as_ref()
            .map(|v| PyCounts::integer(py, v))
            .transpose()
    }
}

/// An owned report from a discovered or built diagram set. A cancelled run returns completed=False.
///
/// Examples
/// --------
/// >>> from symbolica.community import hepkit as hep
/// >>> model = hep.Model.standard_model()
/// >>> request = model.process([11, -11], [22, 22]).amplitude(loops=hep.Loops.total(0))
/// >>> discovered = request.generate()
/// >>> report = discovered.report
/// >>> assert report.completed and report.retained_count == len(discovered)
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    name = "GenerationReport",
    module = "symbolica.community.hepkit",
    frozen,
    from_py_object
)]
#[derive(Clone)]
pub struct PyGenerationReport {
    pub(super) inner: GenerationReport,
}
#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "python_stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[pymethods]
impl PyGenerationReport {
    /// Complete untyped topologies accepted during enumeration; infeasible loop partitions contribute none.
    ///
    /// Examples
    /// --------
    /// Using the setup in this class example:
    ///
    /// >>> assert report.topology_count > 0
    #[getter]
    #[gen_stub(override_return_type(type_repr="builtins.int",imports=("builtins")))]
    fn topology_count<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        PyCounts::integer(py, &self.inner.topology_count)
    }
    /// Complete assignments before cut and charge-conjugation filters.
    ///
    /// Examples
    /// --------
    /// Using the setup in this class example:
    ///
    /// >>> assert report.interaction_assignment_count > 0
    #[getter]
    #[gen_stub(override_return_type(type_repr="builtins.int",imports=("builtins")))]
    fn interaction_assignment_count<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        PyCounts::integer(py, &self.inner.interaction_assignment_count)
    }
    /// Diagrams retained at the stage that produced this report.
    ///
    /// Examples
    /// --------
    /// Using the setup in this class example:
    ///
    /// >>> assert report.retained_count > 0
    #[getter]
    #[gen_stub(override_return_type(type_repr="builtins.int",imports=("builtins")))]
    fn retained_count<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        PyCounts::integer(py, &self.inner.retained_count)
    }
    /// Exact number of zero numerators removed by the reported operation.
    ///
    /// Examples
    /// --------
    /// Using the setup in this class example:
    ///
    /// >>> assert report.zero_numerator_count == 0
    #[getter]
    #[gen_stub(override_return_type(type_repr="builtins.int",imports=("builtins")))]
    fn zero_numerator_count<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        PyCounts::integer(py, &self.inner.zero_numerator_count)
    }
    /// Whether the reported operation completed without cancellation.
    ///
    /// Examples
    /// --------
    /// Using the setup in this class example:
    ///
    /// >>> assert report.completed
    #[getter]
    fn completed(&self) -> bool {
        self.inner.completed
    }
    /// Return the report counters as a compact diagnostic.
    ///
    /// Examples
    /// --------
    /// Using the setup in this class example:
    ///
    /// >>> assert "retained_count" in repr(report)
    fn __repr__(&self) -> String {
        format!("{:?}", self.inner)
    }
    /// Render the report counters as a notebook table.
    ///
    /// Examples
    /// --------
    /// Using the setup in this class example:
    ///
    /// >>> assert "<table>" in report._repr_html_()
    pub(super) fn _repr_html_(&self) -> String {
        let rows = [
            ("Topologies", self.inner.topology_count.to_string()),
            (
                "Interaction assignments",
                self.inner.interaction_assignment_count.to_string(),
            ),
            ("retained diagrams", self.inner.retained_count.to_string()),
            (
                "Zero numerators removed",
                self.inner.zero_numerator_count.to_string(),
            ),
        ]
        .into_iter()
        .map(|(k, v)| format!("<tr><th>{k}</th><td>{v}</td></tr>"))
        .collect::<String>();
        format!(
            "<div class=\"feynkit-generation-report\"><strong>Generation report ({})</strong><table>{rows}</table></div>",
            if self.inner.completed {
                "completed"
            } else {
                "cancelled"
            }
        )
    }
}

pub(super) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PyStage>()?;
    module.add_class::<PyCounts>()?;
    module.add_class::<PyLoopCounts>()?;
    module.add_class::<PyGenerationProgress>()?;
    module.add_class::<PyGenerationReport>()?;
    Ok(())
}
