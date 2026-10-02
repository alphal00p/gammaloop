//! Read-only progress inspection and execution on independent graph copies.
use super::{ExecutionMode, SpensoNet, Spensor, SpensorFunctionLibrary, SpensorLibrary};
use linnet::half_edge::involution::HedgePair;
use pyo3::{prelude::*, types::PyTuple};
#[cfg(feature = "python_stubgen")]
use pyo3_stub_gen::derive::{gen_stub_pyclass, gen_stub_pymethods};
#[cfg(not(feature = "python_stubgen"))]
use pyo3_stub_gen_derive::remove_gen_stub;
use spenso::structure::partial::PartialStructureExt;

/// Read-only snapshot of a TensorNetwork's remaining work.
///
/// Obtain this object from ``network.status``. Counts refer to graph structure,
/// not estimated arithmetic cost or memory use; a previously obtained snapshot
/// does not change when the network executes.
///
/// Examples
/// --------
/// >>> from symbolica.community.tensor import Representation, TensorName, TensorExpression
/// >>> space = Representation.euc(2)
/// >>> A = TensorName("M")(space, space)
/// >>> from symbolica.community.tensor import Tensor
/// >>> tensor = Tensor.dense(A, [1.0, 2.0, 3.0, 4.0])
/// >>> from symbolica.community.tensor import TensorNetwork
/// >>> network = TensorNetwork(tensor)
/// >>> progress = network.status
/// >>> progress.complete
/// True
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    frozen,
    from_py_object,
    name = "ExecutionStatus",
    module = "symbolica.community.tensor"
)]
#[derive(Clone)]
pub(crate) struct ExecutionStatus {
    /// Number of remaining nodes, including values and operations.
    ///
    /// Returns
    /// -------
    /// int
    ///     Current graph node count.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import Representation, TensorName, TensorExpression
    /// >>> space = Representation.euc(2)
    /// >>> A = TensorName("M")(space, space)
    /// >>> from symbolica.community.tensor import Tensor
    /// >>> tensor = Tensor.dense(A, [1.0, 2.0, 3.0, 4.0])
    /// >>> from symbolica.community.tensor import TensorNetwork
    /// >>> network = TensorNetwork(tensor)
    /// >>> snapshot_value = network.status.nodes
    #[pyo3(get)]
    pub(crate) nodes: usize,
    /// Number of remaining operation nodes.
    ///
    /// Returns
    /// -------
    /// int
    ///     Pending sums, products, broadcasts, and other operation nodes.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import Representation, TensorName, TensorExpression
    /// >>> space = Representation.euc(2)
    /// >>> A = TensorName("M")(space, space)
    /// >>> from symbolica.community.tensor import Tensor
    /// >>> tensor = Tensor.dense(A, [1.0, 2.0, 3.0, 4.0])
    /// >>> from symbolica.community.tensor import TensorNetwork
    /// >>> network = TensorNetwork(tensor)
    /// >>> snapshot_value = network.status.operations
    #[pyo3(get)]
    pub(crate) operations: usize,
    /// Number of remaining internal index connections.
    ///
    /// Returns
    /// -------
    /// int
    ///     Includes self-traces.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import Representation, TensorName, TensorExpression
    /// >>> space = Representation.euc(2)
    /// >>> A = TensorName("M")(space, space)
    /// >>> from symbolica.community.tensor import Tensor
    /// >>> tensor = Tensor.dense(A, [1.0, 2.0, 3.0, 4.0])
    /// >>> from symbolica.community.tensor import TensorNetwork
    /// >>> network = TensorNetwork(tensor)
    /// >>> snapshot_value = network.status.contractions
    #[pyo3(get)]
    pub(crate) contractions: usize,
    /// Whether one result value remains with no pending work.
    ///
    /// Returns
    /// -------
    /// bool
    ///     True when no operations or index contractions remain.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import Representation, TensorName, TensorExpression
    /// >>> space = Representation.euc(2)
    /// >>> A = TensorName("M")(space, space)
    /// >>> from symbolica.community.tensor import Tensor
    /// >>> tensor = Tensor.dense(A, [1.0, 2.0, 3.0, 4.0])
    /// >>> from symbolica.community.tensor import TensorNetwork
    /// >>> network = TensorNetwork(tensor)
    /// >>> snapshot_value = network.status.complete
    #[pyo3(get)]
    pub(crate) complete: bool,
    pub(crate) ready: Vec<String>,
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "python_stubgen"), remove_gen_stub)]
#[spenso_macros::track_usage(crate::record_usage)]
#[pymethods]
impl ExecutionStatus {
    /// Operations whose inputs are currently available values.
    ///
    /// Returns
    /// -------
    /// tuple of str
    ///     Operation descriptions in the scheduler's current traversal order.
    ///     Preprocessing, including self-traces, can run before these operations.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import Representation, TensorName, TensorExpression
    /// >>> space = Representation.euc(2)
    /// >>> A = TensorName("M")(space, space)
    /// >>> from symbolica.community.tensor import Tensor
    /// >>> tensor = Tensor.dense(A, [1.0, 2.0, 3.0, 4.0])
    /// >>> from symbolica.community.tensor import TensorNetwork
    /// >>> network = TensorNetwork(tensor)
    /// >>> snapshot_value = network.status.ready_operations
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
    /// Number of external tensor axes.
    ///
    /// Returns
    /// -------
    /// int
    ///     Zero for a scalar.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import Representation, TensorName, TensorExpression
    /// >>> space = Representation.euc(2)
    /// >>> A = TensorName("M")(space, space)
    /// >>> from symbolica.community.tensor import Tensor
    /// >>> tensor = Tensor.dense(A, [1.0, 2.0, 3.0, 4.0])
    /// >>> from symbolica.community.tensor import TensorNetwork
    /// >>> network = TensorNetwork(tensor)
    /// >>> network.rank
    /// 2
    #[getter]
    fn rank(&self) -> usize {
        self.structure.rank()
    }

    /// Whether the tensor has no external axes.
    ///
    /// Returns
    /// -------
    /// bool
    ///     True exactly when rank is zero.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import Representation, TensorName, TensorExpression
    /// >>> space = Representation.euc(2)
    /// >>> A = TensorName("M")(space, space)
    /// >>> from symbolica.community.tensor import Tensor
    /// >>> tensor = Tensor.dense(A, [1.0, 2.0, 3.0, 4.0])
    /// >>> from symbolica.community.tensor import TensorNetwork
    /// >>> network = TensorNetwork(tensor)
    /// >>> network.is_scalar
    /// False
    #[getter]
    fn is_scalar(&self) -> bool {
        self.structure.is_scalar()
    }

    /// Dimensions in logical axis order.
    ///
    /// Returns
    /// -------
    /// tuple of int or Expression
    ///     Axis sizes in the order of axes, the current component view.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import Representation, TensorName, TensorExpression
    /// >>> space = Representation.euc(2)
    /// >>> A = TensorName("M")(space, space)
    /// >>> from symbolica.community.tensor import Tensor
    /// >>> tensor = Tensor.dense(A, [1.0, 2.0, 3.0, 4.0])
    /// >>> from symbolica.community.tensor import TensorNetwork
    /// >>> network = TensorNetwork(tensor)
    /// >>> network.shape
    /// (2, 2)
    #[getter]
    #[gen_stub(override_return_type(type_repr = "tuple[int | Expression, ...]"))]
    fn shape(&self, py: Python<'_>) -> PyResult<Py<PyTuple>> {
        crate::metadata::axis_shape(py, self.structure.structure().logical_slots())
    }

    /// Read a snapshot of the remaining graph work without executing it.
    ///
    /// Returns
    /// -------
    /// ExecutionStatus
    ///     Counts and ready operations at the moment of inspection.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import Representation, TensorName, TensorExpression
    /// >>> space = Representation.euc(2)
    /// >>> A = TensorName("M")(space, space)
    /// >>> from symbolica.community.tensor import Tensor
    /// >>> tensor = Tensor.dense(A, [1.0, 2.0, 3.0, 4.0])
    /// >>> from symbolica.community.tensor import TensorNetwork
    /// >>> network = TensorNetwork(tensor)
    /// >>> network.status.complete
    /// True
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

    /// Copy this network and its current execution progress.
    ///
    /// Returns
    /// -------
    /// TensorNetwork
    ///     An independent copy. Changes to it do not alter the original.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import Representation, TensorName, TensorExpression
    /// >>> space = Representation.euc(2)
    /// >>> A = TensorName("M")(space, space)
    /// >>> from symbolica.community.tensor import Tensor
    /// >>> tensor = Tensor.dense(A, [1.0, 2.0, 3.0, 4.0])
    /// >>> from symbolica.community.tensor import TensorNetwork
    /// >>> network = TensorNetwork(tensor)
    /// >>> copied = network.copy()
    /// >>> copied.shape == network.shape
    /// True
    fn copy(&self) -> Self {
        self.clone()
    }

    /// Run one execution round on a copy of this network.
    ///
    /// Parameters
    /// ----------
    /// library : TensorLibrary, optional
    ///     Component definitions. Defaults to the built-in four-dimensional Dirac
    ///     and SU(3) library. Unregistered tensors receive symbolic components.
    /// function_library : TensorFunctionLibrary, optional
    ///     Numerical implementations of broadcast functions. Defaults to the
    ///     built-in function library.
    ///
    /// Returns
    /// -------
    /// TensorNetwork
    ///     Next intermediate computation; this network remains unchanged.
    ///
    /// Notes
    /// -----
    /// The round uses ExecutionMode.Single and includes the usual preprocessing;
    /// it is not a promise to remove exactly one graph node.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import Representation, TensorName, TensorExpression
    /// >>> space = Representation.euc(2)
    /// >>> A = TensorName("M")(space, space)
    /// >>> from symbolica.community.tensor import Tensor
    /// >>> tensor = Tensor.dense(A, [1.0, 2.0, 3.0, 4.0])
    /// >>> from symbolica.community.tensor import TensorNetwork
    /// >>> network = TensorNetwork(tensor)
    /// >>> intermediate = network.step()
    /// >>> intermediate.shape
    /// (2, 2)
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

    /// Evaluate an independent copy of this network and return its components.
    ///
    /// Parameters
    /// ----------
    /// library : TensorLibrary, optional
    ///     Component definitions. Defaults to the built-in four-dimensional Dirac
    ///     and SU(3) library. Unregistered tensors receive symbolic components.
    /// function_library : TensorFunctionLibrary, optional
    ///     Numerical implementations of broadcast functions. Defaults to the
    ///     built-in function library.
    ///
    /// Returns
    /// -------
    /// Tensor
    ///     Resulting components in logical axis order. Dimensions must be concrete.
    ///
    /// Notes
    /// -----
    /// The original network and its execution progress and the supplied libraries are unchanged.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import Representation, TensorName, TensorExpression
    /// >>> space = Representation.euc(2)
    /// >>> A = TensorName("M")(space, space)
    /// >>> from symbolica.community.tensor import Tensor
    /// >>> tensor = Tensor.dense(A, [1.0, 2.0, 3.0, 4.0])
    /// >>> from symbolica.community.tensor import TensorNetwork
    /// >>> network = TensorNetwork(tensor)
    /// >>> network.to_tensor().shape
    /// (2, 2)
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
