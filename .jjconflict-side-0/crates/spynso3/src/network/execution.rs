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
///
/// Examples
/// --------
/// >>> from symbolica.community import tensor as sp
/// >>> r = sp.Representation.euc(2)
/// >>> A = sp.TensorName("docs::A")(r, r)
/// >>> tensor = sp.Tensor.dense(A, [1.0, 2.0, 3.0, 4.0])
/// >>> network = tensor("i", "j") * tensor("j", "k")
/// >>> status = network.status
/// >>> result = status.complete
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    frozen,
    from_py_object,
    name = "ExecutionStatus",
    module = "symbolica.community.tensor"
)]
#[derive(Clone)]
pub(crate) struct ExecutionStatus {
    /// Number of remaining graph nodes, including leaves and operations.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> r = sp.Representation.euc(2)
    /// >>> A = sp.TensorName("docs::A")(r, r)
    /// >>> tensor = sp.Tensor.dense(A, [1.0, 2.0, 3.0, 4.0])
    /// >>> network = tensor("i", "j") * tensor("j", "k")
    /// >>> status = network.status
    /// >>> result = status.nodes
    #[pyo3(get)]
    pub(crate) nodes: usize,
    /// Number of remaining operation nodes.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> r = sp.Representation.euc(2)
    /// >>> A = sp.TensorName("docs::A")(r, r)
    /// >>> tensor = sp.Tensor.dense(A, [1.0, 2.0, 3.0, 4.0])
    /// >>> network = tensor("i", "j") * tensor("j", "k")
    /// >>> status = network.status
    /// >>> result = status.operations
    #[pyo3(get)]
    pub(crate) operations: usize,
    /// Number of internal slot edges, including self traces.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> r = sp.Representation.euc(2)
    /// >>> A = sp.TensorName("docs::A")(r, r)
    /// >>> tensor = sp.Tensor.dense(A, [1.0, 2.0, 3.0, 4.0])
    /// >>> network = tensor("i", "j") * tensor("j", "k")
    /// >>> status = network.status
    /// >>> result = status.contractions
    #[pyo3(get)]
    pub(crate) contractions: usize,
    /// True when one leaf remains, without pending operations or contractions.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> r = sp.Representation.euc(2)
    /// >>> A = sp.TensorName("docs::A")(r, r)
    /// >>> tensor = sp.Tensor.dense(A, [1.0, 2.0, 3.0, 4.0])
    /// >>> network = tensor("i", "j") * tensor("j", "k")
    /// >>> status = network.status
    /// >>> result = status.complete
    #[pyo3(get)]
    pub(crate) complete: bool,
    pub(crate) ready: Vec<String>,
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "python_stubgen"), remove_gen_stub)]
#[spenso_macros::track_usage(crate::record_usage)]
#[pymethods]
impl ExecutionStatus {
    /// Operations whose children are leaves, in the native scheduler's traversal order.
    /// Preprocessing (including self traces) can run before these operations.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> r = sp.Representation.euc(2)
    /// >>> A = sp.TensorName("docs::A")(r, r)
    /// >>> tensor = sp.Tensor.dense(A, [1.0, 2.0, 3.0, 4.0])
    /// >>> network = tensor("i", "j") * tensor("j", "k")
    /// >>> status = network.status
    /// >>> result = status.ready_operations
    #[getter]
    #[gen_stub(override_return_type(type_repr = "tuple[str, ...]"))]
    fn ready_operations<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyTuple>> {
        PyTuple::new(py, &self.ready)
    }

    /// Return a readable object description for inspection.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> r = sp.Representation.euc(2)
    /// >>> A = sp.TensorName("docs::A")(r, r)
    /// >>> tensor = sp.Tensor.dense(A, [1.0, 2.0, 3.0, 4.0])
    /// >>> network = tensor("i", "j") * tensor("j", "k")
    /// >>> status = network.status
    /// >>> text = repr(status)
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
#[spenso_macros::track_usage(crate::record_usage)]
#[pymethods]
impl SpensoNet {
    /// Number of external axes in the semantic source interface.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> r = sp.Representation.euc(2)
    /// >>> A = sp.TensorName("docs::A")(r, r)
    /// >>> tensor = sp.Tensor.dense(A, [1.0, 2.0, 3.0, 4.0])
    /// >>> network = tensor("i", "j") * tensor("j", "k")
    /// >>> result = network.rank
    #[getter]
    fn rank(&self) -> usize {
        self.structure.rank()
    }

    /// Whether the network's result has no external axes.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> r = sp.Representation.euc(2)
    /// >>> A = sp.TensorName("docs::A")(r, r)
    /// >>> tensor = sp.Tensor.dense(A, [1.0, 2.0, 3.0, 4.0])
    /// >>> network = tensor("i", "j") * tensor("j", "k")
    /// >>> result = network.is_scalar
    #[getter]
    fn is_scalar(&self) -> bool {
        self.structure.is_scalar()
    }

    /// Dimensions in logical axis order; symbolic dimensions remain symbolic.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> r = sp.Representation.euc(2)
    /// >>> A = sp.TensorName("docs::A")(r, r)
    /// >>> tensor = sp.Tensor.dense(A, [1.0, 2.0, 3.0, 4.0])
    /// >>> network = tensor("i", "j") * tensor("j", "k")
    /// >>> result = network.shape
    #[getter]
    #[gen_stub(override_return_type(type_repr = "tuple[int | Expression, ...]"))]
    fn shape(&self, py: Python<'_>) -> PyResult<Py<PyTuple>> {
        self.structure().shape(py)
    }

    /// Snapshot of remaining operations and contractions; does not execute the graph.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> r = sp.Representation.euc(2)
    /// >>> A = sp.TensorName("docs::A")(r, r)
    /// >>> tensor = sp.Tensor.dense(A, [1.0, 2.0, 3.0, 4.0])
    /// >>> network = tensor("i", "j") * tensor("j", "k")
    /// >>> result = network.status
    #[getter]
    pub(crate) fn status(&self) -> ExecutionStatus {
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
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> r = sp.Representation.euc(2)
    /// >>> A = sp.TensorName("docs::A")(r, r)
    /// >>> tensor = sp.Tensor.dense(A, [1.0, 2.0, 3.0, 4.0])
    /// >>> network = tensor("i", "j") * tensor("j", "k")
    /// >>> independent = network.copy()
    /// >>> independent.execute()
    /// >>> assert not network.status.complete
    fn copy(&self) -> Self {
        self.clone()
    }

    /// Return the next intermediate graph, leaving this network unchanged.
    /// Runs one native Single step, including the executor's usual preprocessing.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> r = sp.Representation.euc(2)
    /// >>> A = sp.TensorName("docs::A")(r, r)
    /// >>> tensor = sp.Tensor.dense(A, [1.0, 2.0, 3.0, 4.0])
    /// >>> network = tensor("i", "j") * tensor("j", "k")
    /// >>> next_step = network.step()
    /// >>> assert not network.status.complete
    ///
    /// Parameters
    /// ----------
    /// library : TensorLibrary, optional
    ///     Component definitions for named tensors. None uses the built-in four-dimensional HEP and SU(3) library; supply a TensorLibrary for your own tensors.
    /// function_library : TensorFunctionLibrary, optional
    ///     Callbacks for custom elementwise functions during component evaluation. None uses only the available built-in functions.
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
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> r = sp.Representation.euc(2)
    /// >>> A = sp.TensorName("docs::A")(r, r)
    /// >>> tensor = sp.Tensor.dense(A, [1.0, 2.0, 3.0, 4.0])
    /// >>> network = tensor("i", "j") * tensor("j", "k")
    /// >>> result = network.to_tensor()
    /// >>> assert not network.status.complete
    ///
    /// Parameters
    /// ----------
    /// library : TensorLibrary, optional
    ///     Component definitions for named tensors. None uses the built-in four-dimensional HEP and SU(3) library; supply a TensorLibrary for your own tensors.
    /// function_library : TensorFunctionLibrary, optional
    ///     Callbacks for custom elementwise functions during component evaluation. None uses only the available built-in functions.
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
