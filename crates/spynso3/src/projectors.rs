//! Checked factor groups for the existing chain and trace constructors.
use idenso::tensor::SymbolicTensor;
use pyo3::{
    exceptions::PyValueError,
    prelude::*,
    types::{PyGenericAlias, PyTuple, PyType},
};
#[cfg(feature = "python_stubgen")]
use pyo3_stub_gen::derive::*;
#[cfg(not(feature = "python_stubgen"))]
use pyo3_stub_gen_derive::remove_gen_stub;
use spenso::{shadowing::FactorProjector, structure::partial::PartialStructure};

use crate::{
    ModuleInit,
    expression::{TensorExpression, TensorOperand},
    network::{ConvertibleToSpensoNet, SpensoNet},
    structure::ArithmeticStructure,
};

/// A normalized permutation group of matrix factors for a chain or trace.
///
/// Symmetric and antisymmetric groups average all permutations with weight
/// 1/n!, with signs for the antisymmetric case. Cyclic groups average n cyclic
/// rotations with weight 1/n. Groups may be nested and remain compact until
/// expanded. They are factor groups, not standalone tensor expressions.
///
/// Symbolic factors stay symbolic. A group containing Tensor or TensorNetwork
/// data produces a TensorNetwork when inserted in a chain or trace. Supply
/// individual factors rather than a precomposed multi-factor chain.
///
/// Examples
/// --------
/// >>> from symbolica.community.tensor import Representation, TensorName, TensorExpression
/// >>> space = Representation.euc(2)
/// >>> A = TensorName("M")(space, space)
/// >>> from symbolica.community.tensor import FactorProjector, trace
/// >>> group = FactorProjector.symmetric(A, A)
/// >>> trace(space, group).is_scalar
/// True
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    frozen,
    from_py_object,
    name = "FactorProjector",
    module = "symbolica.community.tensor"
)]
#[derive(Clone)]
pub struct PyFactorProjector {
    pub(crate) structure: SymbolicTensor<PartialStructure>,
    pub(crate) network: Option<SpensoNet>,
    kind: FactorProjector,
    count: usize,
}

impl ModuleInit for PyFactorProjector {}

impl PyFactorProjector {
    fn build(
        py: Python<'_>,
        factors: &Bound<'_, PyTuple>,
        kind: FactorProjector,
    ) -> PyResult<Self> {
        let mut symbolic = Vec::new();
        let mut concrete = false;
        for factor in factors.iter() {
            if let Ok(projector) = factor.extract::<Self>() {
                concrete |= projector.network.is_some();
                symbolic.push(projector.structure);
            } else if let Ok(network) = factor.extract::<SpensoNet>() {
                concrete = true;
                symbolic.push(network.structure);
            } else if let Ok(tensor) = factor.extract::<crate::Spensor>() {
                concrete = true;
                symbolic.push(tensor.descriptor);
            } else {
                symbolic.push(TensorOperand::extract(&factor)?.into_structured()?);
            }
        }
        let structure = SymbolicTensor::project_factors(&symbolic, kind)
            .map_err(|error| PyValueError::new_err(error.to_string()))?;
        let network = if concrete {
            let factors = factors
                .iter()
                .map(|factor| {
                    factor
                        .extract::<ConvertibleToSpensoNet>()
                        .map(ConvertibleToSpensoNet::to_net)
                })
                .collect::<PyResult<Vec<_>>>()?;
            Some(SpensoNet::project_factors(py, factors, kind)?)
        } else {
            None
        };
        Ok(Self {
            structure,
            network,
            kind,
            count: symbolic.len(),
        })
    }

    pub(crate) fn to_network(&self, py: Python<'_>) -> PyResult<SpensoNet> {
        if let Some(network) = &self.network {
            return Ok(network.clone());
        }
        let expression = TensorExpression::from_structured(py, self.structure.clone())?;
        SpensoNet::from_arithmetic(ArithmeticStructure::Tensor(expression), None)
    }
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "python_stubgen"), remove_gen_stub)]
#[spenso_macros::track_usage(crate::record_usage, on_success)]
#[pymethods]
impl PyFactorProjector {
    // Generic in stubs; support evaluated annotations such as FactorProjector[TensorExpression].
    #[gen_stub(skip)]
    #[classmethod]
    fn __class_getitem__<'py>(
        cls: &Bound<'py, PyType>,
        argument: &Bound<'py, PyAny>,
    ) -> PyResult<Bound<'py, PyGenericAlias>> {
        PyGenericAlias::new(cls.py(), cls, argument)
    }

    /// Average all permutations with weight 1/n!.
    ///
    /// Parameters
    /// ----------
    /// *factors : TensorExpression, Tensor, TensorNetwork, or FactorProjector
    ///     Individual matrix factors with compatible channels, or nested groups.
    ///     Their order defines the starting permutation.
    ///
    /// Returns
    /// -------
    /// FactorProjector
    ///     Compact normalized group, retaining any supplied component data.
    ///     Insert it with chain() or trace().
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import Representation, TensorName, TensorExpression
    /// >>> space = Representation.euc(2)
    /// >>> A = TensorName("M")(space, space)
    /// >>> from symbolica.community.tensor import FactorProjector, trace
    /// >>> group = FactorProjector.symmetric(A, A)
    /// >>> result = trace(space, group)
    #[staticmethod]
    #[pyo3(signature = (*factors))]
    fn symmetric(
        py: Python<'_>,
        #[gen_stub(override_type(
            type_repr = "TensorExpression | Tensor | TensorNetwork | FactorProjector"
        ))]
        factors: &Bound<'_, PyTuple>,
    ) -> PyResult<Self> {
        Self::build(py, factors, FactorProjector::Symmetric)
    }

    /// Average all signed permutations with weight 1/n!.
    ///
    /// Parameters
    /// ----------
    /// *factors : TensorExpression, Tensor, TensorNetwork, or FactorProjector
    ///     Individual matrix factors with compatible channels, or nested groups.
    ///     Their order defines the starting permutation.
    ///
    /// Returns
    /// -------
    /// FactorProjector
    ///     Compact normalized group, retaining any supplied component data.
    ///     Insert it with chain() or trace().
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import Representation, TensorName, TensorExpression
    /// >>> space = Representation.euc(2)
    /// >>> A = TensorName("M")(space, space)
    /// >>> from symbolica.community.tensor import FactorProjector, trace
    /// >>> group = FactorProjector.antisymmetric(A, A)
    /// >>> result = trace(space, group)
    #[staticmethod]
    #[pyo3(signature = (*factors))]
    fn antisymmetric(
        py: Python<'_>,
        #[gen_stub(override_type(
            type_repr = "TensorExpression | Tensor | TensorNetwork | FactorProjector"
        ))]
        factors: &Bound<'_, PyTuple>,
    ) -> PyResult<Self> {
        Self::build(py, factors, FactorProjector::Antisymmetric)
    }

    /// Average cyclic rotations with weight 1/n.
    ///
    /// Parameters
    /// ----------
    /// *factors : TensorExpression, Tensor, TensorNetwork, or FactorProjector
    ///     Individual matrix factors with compatible channels, or nested groups.
    ///     Their order defines the starting permutation.
    ///
    /// Returns
    /// -------
    /// FactorProjector
    ///     Compact normalized group, retaining any supplied component data.
    ///     Insert it with chain() or trace().
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import Representation, TensorName, TensorExpression
    /// >>> space = Representation.euc(2)
    /// >>> A = TensorName("M")(space, space)
    /// >>> from symbolica.community.tensor import FactorProjector, trace
    /// >>> group = FactorProjector.cyclic(A, A)
    /// >>> result = trace(space, group)
    #[staticmethod]
    #[pyo3(signature = (*factors))]
    fn cyclic(
        py: Python<'_>,
        #[gen_stub(override_type(
            type_repr = "TensorExpression | Tensor | TensorNetwork | FactorProjector"
        ))]
        factors: &Bound<'_, PyTuple>,
    ) -> PyResult<Self> {
        Self::build(py, factors, FactorProjector::Cyclic)
    }

    fn __repr__(&self) -> String {
        format!("FactorProjector({:?}, {} factors)", self.kind, self.count)
    }

    /// Show the checked factor group using the tensor index printer.
    fn _repr_html_(&self, py: Python<'_>) -> PyResult<String> {
        TensorExpression::from_structured(py, self.structure.clone())?
            .bind(py)
            .call_method0("_repr_html_")?
            .extract()
    }
}
