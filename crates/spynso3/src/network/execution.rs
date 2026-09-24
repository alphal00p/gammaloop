//! Read-only progress inspection and execution on independent graph copies.
use super::{ExecutionMode, SpensoNet, Spensor, SpensorFunctionLibrary, SpensorLibrary};
use linnet::half_edge::involution::HedgePair;
use pyo3::{prelude::*, types::PyTuple};
#[cfg(feature = "python_stubgen")]
use pyo3_stub_gen::derive::{gen_stub_pyclass, gen_stub_pymethods};
#[cfg(not(feature = "python_stubgen"))]
use pyo3_stub_gen_derive::remove_gen_stub;

/// Immutable snapshot of the current graph, independent of its source expression.
/// Counts describe remaining graph structure, not estimated FLOPs or memory costs.
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    frozen,
    from_py_object,
    name = "ExecutionStatus",
    module = "symbolica.community.spenso"
)]
#[derive(Clone)]
pub(crate) struct ExecutionStatus {
    /// Number of remaining graph nodes, including leaves and operations.
    #[pyo3(get)]
    nodes: usize,
    /// Number of remaining operation nodes.
    #[pyo3(get)]
    operations: usize,
    /// Number of internal slot edges, including self traces.
    #[pyo3(get)]
    contractions: usize,
    /// True when one leaf remains, without pending operations or contractions.
    #[pyo3(get)]
    complete: bool,
    ready: Vec<String>,
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "python_stubgen"), remove_gen_stub)]
#[pymethods]
impl ExecutionStatus {
    /// Operations whose children are leaves, in the native scheduler's traversal order.
    /// Preprocessing (including self traces) can run before these operations.
    #[getter]
    #[gen_stub(override_return_type(type_repr = "tuple[str, ...]"))]
    fn ready_operations<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyTuple>> {
        PyTuple::new(py, &self.ready)
    }

    fn __repr__(&self) -> String {
        format!(
            "ExecutionStatus(complete={}, nodes={}, operations={}, contractions={}, ready_operations={:?})",
            if self.complete { "True" } else { "False" },
            self.nodes,
            self.operations,
            self.contractions,
            self.ready
        )
    }
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "python_stubgen"), remove_gen_stub)]
#[pymethods]
impl SpensoNet {
    /// Number of external axes in the semantic source interface.
    #[getter]
    fn rank(&self) -> usize {
        self.structure.rank()
    }

    /// Whether the network's result has no external axes.
    #[getter]
    fn is_scalar(&self) -> bool {
        self.structure.is_scalar()
    }

    /// Dimensions in logical axis order; symbolic dimensions remain symbolic.
    #[getter]
    #[gen_stub(override_return_type(type_repr = "tuple[int | Expression, ...]"))]
    fn shape(&self, py: Python<'_>) -> PyResult<Py<PyTuple>> {
        self.structure().shape(py)
    }

    /// Snapshot of remaining operations and contractions; does not execute the graph.
    #[getter]
    fn status(&self) -> ExecutionStatus {
        let mut graph = self.network.graph.clone();
        graph.cache_expr_tree_roots();
        let nodes = graph.n_nodes();
        let operations = graph
            .graph
            .iter_nodes()
            .filter(|(_, _, data)| data.is_op())
            .count();
        let contractions = graph
            .graph
            .iter_edges()
            .filter(|(pair, _, edge)| {
                edge.data.is_slot() && matches!(pair, HedgePair::Paired { .. })
            })
            .count();
        let ready = graph
            .ready_operation_refs()
            .iter()
            .map(|op| op.op().to_string())
            .collect();
        ExecutionStatus {
            nodes,
            operations,
            contractions,
            complete: nodes == 1 && operations == 0 && contractions == 0,
            ready,
        }
    }

    /// Copy the graph, its component values, and execution progress.
    fn copy(&self) -> Self {
        self.clone()
    }

    /// Return the next intermediate graph, leaving this network unchanged.
    /// Runs one native Single step, including the executor's usual preprocessing.
    #[pyo3(signature=(library=None, *, function_library=None))]
    fn step(
        &self,
        library: Option<&SpensorLibrary>,
        function_library: Option<&SpensorFunctionLibrary>,
    ) -> PyResult<Self> {
        let mut result = self.clone();
        result.execute(library, function_library, Some(1), ExecutionMode::Single)?;
        Ok(result)
    }

    /// Execute an independent copy to completion and return its component tensor.
    #[pyo3(signature=(library=None, *, function_library=None))]
    fn to_tensor(
        &self,
        library: Option<&SpensorLibrary>,
        function_library: Option<&SpensorFunctionLibrary>,
    ) -> PyResult<Spensor> {
        let mut result = self.clone();
        result.execute(library, function_library, None, ExecutionMode::All)?;
        result.result_tensor(library)
    }
}
