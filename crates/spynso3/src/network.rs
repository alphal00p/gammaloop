use std::{collections::HashMap, ops::Deref};

pub(crate) mod execution;

use pyo3::{
    Borrowed,
    exceptions::{self, PyRuntimeError, PyTypeError, PyValueError},
    prelude::*,
    types::{PyDict, PyTuple},
};

#[cfg(not(feature = "python_stubgen"))]
use pyo3_stub_gen_derive::remove_gen_stub;

use spenso::{
    algebra::complex::{Complex, RealOrComplex},
    iterators::IteratableTensor,
    network::{
        ContractScalars, ExecutionResult, Network, Sequential, SingleSmallestDegree,
        SmallestDegree, Steps,
        library::symbolic::{ExplicitKey, TensorLibrary},
        parsing::{ParseSettings, ShadowedStructure},
        store::{NetworkStore, TensorScalarStoreMapping},
        tags::SPENSO_TAG,
    },
    structure::{
        HasStructure, ScalarTensor, TensorStructure,
        abstract_index::AbstractIndex,
        partial::{PartialIndex, PartialStructure, PartialStructureExt},
        slot::IsAbstractSlot,
    },
    tensors::data::SparseTensor,
    tensors::parametric::{
        AtomViewOrConcrete, ConcreteOrParam, MixedTensor, ParamOrConcrete, atomcore::TensorAtomMaps,
    },
};
use spenso_hep_lib::{FUN_LIB, HEP_LIB};
use symbolica::{
    api::python::{
        ConvertibleToReplaceWith, PythonCondition, PythonExpression, PythonFormattedOutput,
        PythonPatternRestriction,
    },
    atom::FunctionBuilder,
    domains::float::Complex as SymComplex,
    id::{Condition, PatternRestriction},
    prelude::*,
};

use symbolica::api::python::ConvertibleToExpression;

use crate::{display, expression::TensorExpression, library::SpensorFunctionLibrary};
use idenso::tensor::{
    SymbolicTensor,
    composition::{self, ProductPlan},
    inference::InterfaceInference,
};

use super::{Spensor, library::SpensorLibrary, structure::ArithmeticStructure};

use super::ModuleInit;

#[cfg(feature = "python_stubgen")]
use pyo3_stub_gen::{PyStubType, derive::*};

/// A graph of tensor operations that can be simplified and executed.
///
/// Named tensor expressions are resolved through a `TensorLibrary`. Register tensors
/// to supply their component values. Unregistered tensors with concrete dimensions
/// acquire symbolic components automatically; `TensorExpression.components()` exposes
/// the same components without requiring a surrounding network.
///
/// A network retains the semantic source expression and its public tensor interface
/// separately from the executable graph and its stored values. Value specialization
/// and graph execution therefore do not rewrite the source expression returned by
/// `expression()` or used by the semantic display methods. The default rich
/// display draws the current executable graph through Linnest.
///
/// Examples
/// --------
/// >>> from symbolica.community.spenso import (
/// ...     ExecutionMode,
/// ...     Representation,
/// ...     Tensor,
/// ...     TensorLibrary,
/// ...     TensorName,
/// ...     TensorNetwork,
/// ... )
/// >>> rep = Representation.euc(2)
/// >>> A = TensorName("A")
/// >>> structure = A(rep, rep)
/// >>> library = TensorLibrary()
/// >>> library.register(
/// ...     Tensor.dense(structure, [1.0, 0.0, 0.0, 1.0])
/// ... )
/// >>> network = TensorNetwork(
/// ...     A(rep("i"), rep("j")),
/// ...     library=library,
/// ... )
/// >>> network.execute(library=library, mode=ExecutionMode.All)
/// >>> result = network.result_tensor(library=library)
/// >>> len(result)
/// 4
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    from_py_object,
    name = "TensorNetwork",
    module = "symbolica.community.spenso"
)]
#[derive(Clone)]
#[allow(clippy::type_complexity)]
pub struct SpensoNet {
    pub network: Network<
        NetworkStore<MixedTensor<f64, ShadowedStructure<AbstractIndex>>, Atom>,
        ExplicitKey<AbstractIndex>,
        Symbol,
    >,
    /// Semantic source expression and public interface, including unresolved ports.
    pub(crate) structure: SymbolicTensor<PartialStructure>,
    /// The source expression with all graph-facing ports made explicit.
    pub(crate) materialized: SymbolicTensor<PartialStructure>,
    /// Optional stored-data identity. Composite results deliberately remain unnamed.
    pub(crate) descriptor: Option<(Symbol, Vec<Atom>)>,
}

fn insert_function_values<'a>(
    expression: AtomView<'a>,
    constants: &mut ahash::AHashMap<Atom, f64>,
    functions: &HashMap<Symbol, Py<PyAny>>,
) -> PyResult<()> {
    match expression {
        AtomView::Fun(fun) => {
            for arg in fun.iter() {
                insert_function_values(arg, constants, functions)?;
            }

            if let Some(function) = functions.get(&fun.get_symbol()) {
                let args = fun
                    .iter()
                    .map(|arg| {
                        arg.evaluate(constants).map_err(|err| {
                            exceptions::PyValueError::new_err(format!(
                                "Could not evaluate function argument in `{expression}`: {err}"
                            ))
                        })
                    })
                    .collect::<PyResult<Vec<_>>>()?;
                let value =
                    Python::attach(|py| function.call(py, (args,), None)?.extract::<f64>(py))?;
                constants.insert(expression.to_owned(), value);
            }
        }
        AtomView::Add(add) => {
            for arg in add.iter() {
                insert_function_values(arg, constants, functions)?;
            }
        }
        AtomView::Mul(mul) => {
            for arg in mul.iter() {
                insert_function_values(arg, constants, functions)?;
            }
        }
        AtomView::Pow(pow) => {
            let (base, exponent) = pow.get_base_exp();
            insert_function_values(base, constants, functions)?;
            insert_function_values(exponent, constants, functions)?;
        }
        _ => {}
    }

    Ok(())
}

/// Execution modes for tensor network evaluation.
///
/// Controls how the tensor network execution engine processes the computational graph.
///
/// Variants
/// --------
/// Single : Select one smallest-degree rewrite per step; without `n_steps`, continue until no work remains
/// Scalar : Only contract scalar operations, leaving tensor structure intact
/// All : Execute all possible contractions for complete evaluation
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass_enum)]
#[pyclass(
    from_py_object,
    name = "ExecutionMode",
    module = "symbolica.community.spenso"
)]
#[derive(Clone)]
pub enum ExecutionMode {
    Single,
    Scalar,
    All,
}

impl ModuleInit for ExecutionMode {}

impl ModuleInit for SpensoNet {
    fn init(m: &Bound<'_, PyModule>) -> PyResult<()> {
        m.add_class::<SpensoNet>()
    }
}

pub type ParsingNet = Network<
    NetworkStore<MixedTensor<f64, ShadowedStructure<AbstractIndex>>, Atom>,
    ExplicitKey<AbstractIndex>,
    Symbol,
>;

/// A Symbolica pattern restriction accepted by tensor replacement.
pub struct ReplacementCondition(pub(crate) Condition<PatternRestriction>);

impl<'a, 'py> FromPyObject<'a, 'py> for ReplacementCondition {
    type Error = PyErr;

    fn extract(value: Borrowed<'a, 'py, PyAny>) -> Result<Self, Self::Error> {
        if let Ok(restriction) = value.extract::<PythonPatternRestriction>() {
            Ok(Self(restriction.condition))
        } else if let Ok(condition) = value.extract::<PythonCondition>() {
            condition
                .condition
                .try_into()
                .map(Self)
                .map_err(|error: &'static str| PyValueError::new_err(error))
        } else {
            Err(PyTypeError::new_err(
                "expected a Symbolica PatternRestriction or Condition",
            ))
        }
    }
}

#[cfg(feature = "python_stubgen")]
pyo3_stub_gen::impl_stub_type!(ReplacementCondition = PythonPatternRestriction | PythonCondition);

fn non_scalar_zero_network(value: &SymbolicTensor<PartialStructure>) -> Option<ParsingNet> {
    if !value.expression().as_view().is_zero() || value.is_scalar() {
        return None;
    }

    let dummies = std::cell::OnceCell::new();
    let replacements = value
        .structure()
        .logical_slots()
        .into_iter()
        .filter_map(|slot| match slot.aind {
            PartialIndex::Explicit(_) => None,
            PartialIndex::Open(id) => Some((
                id,
                dummies
                    .get_or_init(|| SymbolicTensor::reserved_dummies([value]))
                    .fresh_index(),
            )),
        })
        .collect();
    let structure: ShadowedStructure<AbstractIndex> = value
        .structure()
        .materialize_open_ports(&replacements)
        .into_canonical()
        .into();
    let zero = SparseTensor::empty(structure, 0.0);
    Some(Network::from_tensor(MixedTensor::from(zero)))
}

fn network_from_arithmetic(
    expression: ArithmeticStructure,
    library: &TensorLibrary<MixedTensor<f64, ExplicitKey<AbstractIndex>>, AbstractIndex>,
) -> eyre::Result<SpensoNet> {
    match expression {
        ArithmeticStructure::Tensor(expression) => Python::attach(|py| {
            let expression = expression.bind(py).borrow();
            let structure = TensorExpression::structured(&expression).clone();
            let materialized = structure.materialized().map_err(composition_error)?;
            let descriptor = TensorExpression::descriptor_name(&expression)
                .map(|name| (name, TensorExpression::descriptor_args(&expression)));
            if let Some(network) = non_scalar_zero_network(&materialized) {
                return Ok(SpensoNet {
                    network,
                    structure,
                    materialized,
                    descriptor,
                });
            }

            Ok(SpensoNet {
                network: ParsingNet::try_from_view(
                    materialized.expression().as_view(),
                    library,
                    &ParseSettings::default(),
                )?,
                structure,
                materialized,
                descriptor,
            })
        }),
        expression => {
            let atom = expression.to_expression()?.expr;
            let network =
                ParsingNet::try_from_view(atom.as_view(), library, &ParseSettings::default())?;
            let value = SymbolicTensor::new(
                atom.clone(),
                PartialStructure::from_logical_slots(
                    network
                        .graph
                        .dangling_indices()
                        .into_iter()
                        .map(|slot| slot.rep().slot(PartialIndex::Explicit(slot.aind()))),
                ),
            )
            .normalize_closed_root_chain()?;
            Ok(SpensoNet {
                network,
                structure: value.clone(),
                materialized: value,
                descriptor: None,
            })
        }
    }
}

pub struct ConvertibleToSpensoNet(pub(crate) SpensoNet);

impl ConvertibleToSpensoNet {
    pub fn to_net(self) -> SpensoNet {
        self.0
    }
}

impl SpensoNet {
    pub(crate) fn from_tensor(value: Spensor) -> PyResult<Self> {
        let canonical = value.tensor.external_structure();
        let logical = value
            .descriptor
            .structure()
            .layout()
            .canonical_to_logical(&canonical);
        let replacements = logical
            .into_iter()
            .enumerate()
            .map(|(position, slot)| (position, slot.aind()))
            .collect::<HashMap<_, _>>();
        let materialized = value
            .descriptor
            .reindex_interface_ports(&replacements)
            .map_err(composition_error)?;
        let network = Self {
            network: Network::from_tensor(value.tensor),
            materialized,
            structure: value.descriptor,
            descriptor: value
                .descriptor_name
                .map(|name| (name, value.descriptor_args)),
        };
        let indices = network
            .semantic_slots()
            .into_iter()
            .enumerate()
            .filter_map(|(position, slot)| match slot.aind {
                PartialIndex::Explicit(index) => Some((position, index)),
                PartialIndex::Open(_) => None,
            })
            .collect();
        let network = network.relabel_ports(&indices)?;
        network.validate_graph_interface()?;
        Ok(network)
    }

    pub(crate) fn from_tensor_reference(value: Spensor) -> PyResult<Self> {
        Self::from_tensor(value)
    }

    pub(crate) fn from_arithmetic(
        expr: ArithmeticStructure,
        library: Option<&SpensorLibrary>,
    ) -> PyResult<Self> {
        let library = library
            .map(|value| &value.library)
            .unwrap_or(HEP_LIB.deref());
        network_from_arithmetic(expr, library)
            .map_err(|error| PyRuntimeError::new_err(error.to_string()))
    }
}

impl<'a, 'py> FromPyObject<'a, 'py> for ConvertibleToSpensoNet {
    type Error = PyErr;

    fn extract(ob: pyo3::Borrowed<'a, 'py, pyo3::PyAny>) -> Result<Self, Self::Error> {
        if let Ok(a) = ob.extract::<SpensoNet>() {
            Ok(ConvertibleToSpensoNet(a))
        } else if let Ok(num) = ob.extract::<Spensor>() {
            Ok(ConvertibleToSpensoNet(SpensoNet::from_tensor(num)?))
        } else if let Ok(a) = ob.extract::<ArithmeticStructure>() {
            SpensoNet::from_arithmetic(a, None).map(ConvertibleToSpensoNet)
        } else {
            Err(exceptions::PyTypeError::new_err(
                "Cannot convert to expression",
            ))
        }
    }
}

#[cfg(feature = "python_stubgen")]
impl PyStubType for ConvertibleToSpensoNet {
    fn type_output() -> pyo3_stub_gen::TypeInfo {
        ArithmeticStructure::type_output() | SpensoNet::type_output() | Spensor::type_output()
    }
}

fn composition_error(error: composition::TensorCompositionError) -> PyErr {
    PyValueError::new_err(error.to_string())
}

impl SpensoNet {
    fn semantic_slots(&self) -> Vec<spenso::structure::partial::PartialSlot> {
        self.structure.structure().logical_slots()
    }

    fn relabel_ports(mut self, indices: &HashMap<usize, AbstractIndex>) -> PyResult<Self> {
        let slots = self.materialized.structure().logical_slots();
        let mut replacements = indices
            .iter()
            .map(|(&position, &index)| {
                let slot = slots.get(position).copied().ok_or_else(|| {
                    PyValueError::new_err(format!(
                        "port position {position} is outside an interface of rank {}",
                        slots.len()
                    ))
                })?;
                let PartialIndex::Explicit(current) = slot.aind else {
                    return Err(PyRuntimeError::new_err(
                        "tensor network graph interface contains an unresolved port",
                    ));
                };
                Ok((
                    position,
                    slot.rep().slot(current).to_lib(),
                    slot.rep().slot(index).to_lib(),
                ))
            })
            .collect::<PyResult<Vec<_>>>()?;
        replacements.sort_by_key(|(position, _, _)| *position);
        if replacements.iter().all(|(_, from, to)| from == to) {
            return Ok(self);
        }
        let materialized = self
            .materialized
            .reindex_interface_ports(indices)
            .map_err(composition_error)?;
        self.network = self
            .network
            .reindex_ports(
                &replacements
                    .into_iter()
                    .map(|(_, from, to)| (from, to))
                    .collect::<Vec<_>>(),
            )
            .map_err(|error| PyRuntimeError::new_err(error.to_string()))?;
        self.materialized = materialized;
        Ok(self)
    }

    fn apply_product(self, right: Self, plan: ProductPlan) -> PyResult<Self> {
        let structure = plan
            .apply(&self.structure, &right.structure)
            .map_err(composition_error)?;
        let (left_indices, right_indices) = plan
            .alignment_indices(&self.structure, &right.structure)
            .map_err(composition_error)?;
        let left = self.relabel_ports(&left_indices)?;
        let right = right.relabel_ports(&right_indices)?;
        let materialized = plan
            .apply(&left.materialized, &right.materialized)
            .map_err(composition_error)?;
        Self::finish(left.network * right.network, structure, materialized)
    }

    fn validate_graph_interface(&self) -> PyResult<()> {
        let mut expected = self
            .materialized
            .structure()
            .logical_slots()
            .into_iter()
            .map(|slot| {
                let PartialIndex::Explicit(index) = slot.aind else {
                    return Err(PyRuntimeError::new_err(
                        "materialized network interface contains an unresolved port",
                    ));
                };
                Ok(slot.rep().slot(index).to_lib())
            })
            .collect::<PyResult<Vec<_>>>()?;
        let mut actual = self.network.graph.dangling_indices();
        expected.sort();
        actual.sort();
        if expected != actual {
            return Err(PyRuntimeError::new_err(format!(
                "tensor-network topology diverged from the planned interface: expected {expected:?}, got {actual:?}"
            )));
        }
        Ok(())
    }

    fn finish(
        mut network: ParsingNet,
        structure: SymbolicTensor<PartialStructure>,
        materialized: SymbolicTensor<PartialStructure>,
    ) -> PyResult<Self> {
        let structure = structure
            .normalize_closed_root_chain()
            .map_err(composition_error)?;
        let materialized = materialized
            .normalize_closed_root_chain()
            .map_err(composition_error)?;
        network.state = network.graph.state();
        let value = Self {
            network,
            structure,
            materialized,
            descriptor: None,
        };
        value.validate_graph_interface()?;
        Ok(value)
    }

    pub(crate) fn add_network(mut self, mut right: Self, subtract: bool) -> PyResult<Self> {
        if !InterfaceInference::additive_interfaces_match(
            self.structure.structure(),
            right.structure.structure(),
        ) {
            return Err(PyValueError::new_err(if subtract {
                "subtraction requires compatible tensor interfaces"
            } else {
                "addition requires compatible tensor interfaces"
            }));
        }
        let pairs = self
            .structure
            .structure()
            .open_positions()
            .into_iter()
            .map(|position| composition::PortPair {
                left: position,
                right: position,
            })
            .collect::<Vec<_>>();
        let (left_indices, right_indices) = self
            .structure
            .aligned_indices(&right.structure, &pairs)
            .map_err(composition_error)?;
        self = self.relabel_ports(&left_indices)?;
        right = right.relabel_ports(&right_indices)?;
        let structure = SymbolicTensor::new(
            if subtract {
                self.structure.expression() - right.structure.expression()
            } else {
                self.structure.expression() + right.structure.expression()
            },
            self.structure.structure().clone(),
        );
        let materialized = SymbolicTensor::new(
            if subtract {
                self.materialized.expression() - right.materialized.expression()
            } else {
                self.materialized.expression() + right.materialized.expression()
            },
            self.materialized.structure().clone(),
        );
        let network = if subtract {
            self.network - right.network
        } else {
            self.network + right.network
        };
        Self::finish(network, structure, materialized)
    }

    pub(crate) fn multiply_network(self, right: Self) -> PyResult<Self> {
        let plan = self
            .structure
            .product_plan(&right.structure)
            .map_err(composition_error)?;
        self.apply_product(right, plan)
    }

    fn reciprocal(self) -> PyResult<Self> {
        if !self.structure.is_scalar() {
            return Err(PyValueError::new_err(
                "a non-scalar tensor cannot be used as a denominator",
            ));
        }
        let structure = SymbolicTensor::new(
            Atom::num(1) / self.structure.expression(),
            self.structure.structure().clone(),
        );
        let materialized = SymbolicTensor::new(
            Atom::num(1) / self.materialized.expression(),
            self.materialized.structure().clone(),
        );
        Self::finish(self.network.pow(-1), structure, materialized)
    }

    pub(crate) fn broadcast_function(self, function: Symbol) -> PyResult<Self> {
        let structure = SymbolicTensor::new(
            FunctionBuilder::new(function)
                .add_arg(self.structure.expression().clone())
                .finish(),
            self.structure.structure().clone(),
        );
        let materialized = SymbolicTensor::new(
            FunctionBuilder::new(function)
                .add_arg(self.materialized.expression().clone())
                .finish(),
            self.materialized.structure().clone(),
        );
        Self::finish(self.network.fun(function), structure, materialized)
    }

    pub(crate) fn index_network(
        &self,
        indices: &Bound<'_, PyTuple>,
        cook_indices: Option<&crate::simplification::PyCookSettings>,
    ) -> PyResult<Self> {
        let replacements = TensorExpression::index_replacements(
            self.structure.structure(),
            indices,
            cook_indices,
        )?;
        self.set_port_indices(&replacements)
    }

    pub(crate) fn set_port_indices(
        &self,
        replacements: &HashMap<usize, AbstractIndex>,
    ) -> PyResult<Self> {
        let structure = self
            .structure
            .reindex_interface_ports(replacements)
            .map_err(composition_error)?;
        let indices = self
            .structure
            .graph_port_indices(replacements)
            .map_err(composition_error)?;
        let mut network = self.clone().relabel_ports(&indices)?;
        network.network.graph.sew_dangling_slots();
        network.network.state = network.network.graph.state();

        let initial_rank = structure.rank();
        (network.structure, network.materialized) = structure
            .retain_materialized_ports(
                network.materialized,
                &network.network.graph.dangling_indices(),
            )
            .map_err(composition_error)?;
        let contracted = network.structure.rank() < initial_rank;
        if contracted {
            network.descriptor = None;
        }
        network.validate_graph_interface()?;
        Ok(network)
    }

    pub(crate) fn chain_form(self, channel: composition::MatrixChannel) -> PyResult<Self> {
        let structure = self
            .structure
            .chain_form(channel)
            .map_err(composition_error)?;
        let materialized = self
            .materialized
            .chain_form(channel)
            .map_err(composition_error)?;
        Self::finish(self.network, structure, materialized)
    }
}

// #[gen_stub_pymethods]

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "python_stubgen"), remove_gen_stub)]
#[pymethods]
impl SpensoNet {
    #[new]
    /// Create a tensor network by parsing an arithmetic expression.
    ///
    /// Parses symbolic expressions containing tensor operations and converts them
    /// into an optimizable computational graph representation.
    ///
    /// Parameters
    /// ----------
    /// expr : ArithmeticStructure
    ///     The arithmetic expression or tensor structure to parse
    /// library : TensorLibrary, optional
    ///     Tensor library for resolving named tensor references. Defaults to the built-in
    ///     four-dimensional HEP and SU(3) library returned by `TensorLibrary.hep_lib()`.
    ///
    /// Returns
    /// -------
    /// TensorNetwork
    ///     A new TensorNetwork representing the parsed expression
    #[pyo3(signature = (expr, library=None))]
    pub fn from_expression(
        expr: &Bound<'_, PyAny>,
        library: Option<&SpensorLibrary>,
    ) -> PyResult<SpensoNet> {
        if let Ok(network) = expr.extract::<SpensoNet>() {
            return Ok(network);
        }
        if let Ok(tensor) = expr.extract::<Spensor>() {
            return SpensoNet::from_tensor(tensor);
        }
        let expr = expr.extract::<ArithmeticStructure>().map_err(|_| {
            PyTypeError::new_err(
                "TensorNetwork() expects a Tensor, TensorNetwork, TensorExpression, or scalar expression",
            )
        })?;
        SpensoNet::from_arithmetic(expr, library)
    }

    #[staticmethod]
    /// Create a tensor network representing the scalar value 1.
    ///
    /// Returns
    /// -------
    /// TensorNetwork
    ///     A TensorNetwork containing only the scalar 1
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.spenso import TensorNetwork
    /// >>> one_net = TensorNetwork.one()
    /// >>> result = one_net.result_scalar()
    pub fn one() -> SpensoNet {
        let value = SymbolicTensor::new(
            Atom::num(1),
            PartialStructure::from_logical_slots(std::iter::empty()),
        );
        SpensoNet {
            network: Network::one(),
            structure: value.clone(),
            materialized: value,
            descriptor: None,
        }
    }

    #[staticmethod]
    /// Return the symbolic head used for structured product brackets.
    pub fn bracket() -> PythonExpression {
        PythonExpression {
            expr: Atom::var(SPENSO_TAG.bracket),
        }
    }

    #[staticmethod]
    /// Create a Symbolica function symbol tagged for elementwise tensor broadcasting.
    pub fn broadcast(str: &str) -> PythonExpression {
        PythonExpression {
            expr: Atom::var(symbol!(str, tag = SPENSO_TAG.broadcast)),
        }
    }
    #[staticmethod]
    /// Create a tensor network representing the scalar value 0.
    ///
    /// Returns
    /// -------
    /// TensorNetwork
    ///     A TensorNetwork containing only the scalar 0
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.spenso import TensorNetwork
    /// >>> zero_net = TensorNetwork.zero()
    /// >>> result = zero_net.result_scalar()
    pub fn zero() -> SpensoNet {
        let value = SymbolicTensor::new(
            Atom::Zero,
            PartialStructure::from_logical_slots(std::iter::empty()),
        );
        SpensoNet {
            network: Network::zero(),
            structure: value.clone(),
            materialized: value,
            descriptor: None,
        }
    }

    /// Replace patterns in stored symbolic network values.
    ///
    /// Rewrites scalar coefficients and symbolic elements of parametric tensors in
    /// the execution store. Tensor identities, graph topology, the public interface,
    /// and the semantic source expression returned by `structure()` are unchanged.
    /// Rewrite a `TensorExpression` before constructing the network when the source
    /// tensor expression itself should change.
    ///
    /// Parameters
    /// ----------
    /// pattern : Expression
    ///     The symbolic pattern to match within stored values
    /// rhs : Expression
    ///     The replacement expression or pattern
    /// cond : PatternRestriction or Condition, optional
    ///     Additional restriction that each match must satisfy
    /// non_greedy_wildcards : list of Expression, optional
    ///     List of wildcard symbols to match non-greedily
    /// level_range : tuple of int, optional
    ///     Tuple specifying depth range for pattern matching
    /// level_is_tree_depth : bool, optional
    ///     Whether level refers to tree depth or expression depth
    /// allow_new_wildcards_on_rhs : bool, optional
    ///     Allow new wildcards in replacement pattern
    /// rhs_cache_size : int, optional
    ///     Size of cache for replacement pattern compilation
    /// repeat : bool, optional
    ///     Whether to repeatedly apply the replacement until no more matches
    ///
    /// Returns
    /// -------
    /// TensorNetwork
    ///     A new TensorNetwork with matching stored values replaced and the same
    ///     semantic source structure
    #[pyo3(signature = (pattern, rhs, cond = None, non_greedy_wildcards = None, level_range = None, level_is_tree_depth = None, allow_new_wildcards_on_rhs = None, rhs_cache_size = None, repeat = None))]
    #[allow(clippy::too_many_arguments)]
    pub fn replace(
        &self,
        pattern: ConvertibleToExpression,
        rhs: ConvertibleToReplaceWith,
        cond: Option<ReplacementCondition>,
        non_greedy_wildcards: Option<Vec<PythonExpression>>,
        level_range: Option<(usize, Option<usize>)>,
        level_is_tree_depth: Option<bool>,
        allow_new_wildcards_on_rhs: Option<bool>,
        rhs_cache_size: Option<usize>,
        repeat: Option<bool>,
    ) -> PyResult<SpensoNet> {
        let pattern = pattern.to_expression().expr.to_pattern();
        let ReplaceWith::Pattern(rhs) = &rhs.to_replace_with()? else {
            return Err(exceptions::PyTypeError::new_err(
                "Only normal patterns supported",
            ));
        };

        let mut setting_non_greedy_wildcards = Vec::new();
        let mut setting_level_range = (0, None);
        let mut setting_level_is_tree_depth = false;
        let mut setting_allow_new_wildcards_on_rhs = false;
        let mut setting_rhs_cache_size = 100;

        if let Some(ngw) = non_greedy_wildcards {
            setting_non_greedy_wildcards = ngw
                .iter()
                .map(|x| match x.expr.as_view() {
                    AtomView::Var(v) => {
                        let name = v.get_symbol();
                        if v.get_wildcard_level() == 0 {
                            return Err(exceptions::PyTypeError::new_err(
                                "Only wildcards can be restricted.",
                            ));
                        }
                        Ok(name)
                    }
                    _ => Err(exceptions::PyTypeError::new_err(
                        "Only wildcards can be restricted.",
                    )),
                })
                .collect::<Result<_, _>>()?;
        }
        if let Some(level_range) = level_range {
            setting_level_range = level_range;
        }
        if let Some(level_is_tree_depth) = level_is_tree_depth {
            setting_level_is_tree_depth = level_is_tree_depth;
        }
        if let Some(allow_new_wildcards_on_rhs) = allow_new_wildcards_on_rhs {
            setting_allow_new_wildcards_on_rhs = allow_new_wildcards_on_rhs;
        }
        if let Some(rhs_cache_size) = rhs_cache_size {
            setting_rhs_cache_size = rhs_cache_size;
        }

        let cond = cond.map(|condition| condition.0);

        Ok(SpensoNet {
            network: self.network.map_ref(
                |s| {
                    let r = s.replace(&pattern);
                    let r = if let Some(cond) = cond.as_ref() {
                        r.when(cond)
                    } else {
                        r
                    }
                    .non_greedy_wildcards(setting_non_greedy_wildcards.clone())
                    .min_level(setting_level_range.0)
                    .max_level(setting_level_range.1)
                    .level_is_tree_depth(setting_level_is_tree_depth)
                    .allow_new_wildcards_on_rhs(setting_allow_new_wildcards_on_rhs)
                    .rhs_cache_size(setting_rhs_cache_size);

                    let mut r = if let Some(true) = repeat {
                        r.repeat()
                    } else {
                        r
                    };

                    r.with(rhs.borrow())
                },
                |t| match t {
                    ParamOrConcrete::Param(p) => {
                        let r = p.replace(&pattern);
                        let r = if let Some(cond) = cond.as_ref() {
                            r.when(cond)
                        } else {
                            r
                        }
                        .non_greedy_wildcards(setting_non_greedy_wildcards.clone())
                        .level_range(setting_level_range)
                        .level_is_tree_depth(setting_level_is_tree_depth)
                        .allow_new_wildcards_on_rhs(setting_allow_new_wildcards_on_rhs)
                        .rhs_cache_size(setting_rhs_cache_size);

                        let r = if let Some(true) = repeat {
                            r.repeat()
                        } else {
                            r
                        };

                        ParamOrConcrete::Param(r.with(rhs.borrow()))
                    }
                    _ => t.clone(),
                },
            ),
            structure: self.structure.clone(),
            materialized: self.materialized.clone(),
            descriptor: self.descriptor.clone(),
        })
    }

    /// Evaluate symbolic tensor values in the execution store.
    ///
    /// Substitutes symbolic constants and functions in stored parametric tensor
    /// elements, converting them to concrete numerical tensors. Tensor identities,
    /// graph topology, the public interface, and the semantic source expression are
    /// retained.
    ///
    /// Parameters
    /// ----------
    /// constants : dict
    ///     Dict mapping symbolic expressions to their numerical values
    /// functions : dict
    ///     Dict mapping function symbols to Python callable objects
    ///
    /// Returns
    /// -------
    /// TensorNetwork
    ///     A new TensorNetwork with stored tensor values evaluated and the same
    ///     semantic source structure
    pub fn evaluate(
        &self,
        constants: HashMap<PythonExpression, f64>,
        functions: HashMap<PolyVariable, Py<PyAny>>,
    ) -> PyResult<Self> {
        let mut constants: ahash::AHashMap<Atom, f64> = constants
            .iter()
            .map(|(k, v)| (k.expr.clone(), *v))
            .collect();

        let functions = functions
            .into_iter()
            .map(|(k, v)| {
                let id = if let PolyVariable::Symbol(v) = k {
                    v
                } else {
                    Err(exceptions::PyValueError::new_err(format!(
                        "Expected function name instead of {:?}",
                        k
                    )))?
                };

                Ok((id, v))
            })
            .collect::<PyResult<_>>()?;

        let mut network = self.network.clone();
        for tensor in network.iter_tensors() {
            for (_, value) in tensor.iter_flat() {
                if let AtomViewOrConcrete::Atom(atom) = value {
                    insert_function_values(atom, &mut constants, &functions)?;
                }
            }
        }
        network.evaluate_real(&constants);
        Ok(SpensoNet {
            network,
            structure: self.structure.clone(),
            materialized: self.materialized.clone(),
            descriptor: self.descriptor.clone(),
        })
    }

    /// Execute the tensor network to perform tensor contractions and simplifications.
    ///
    /// Processes the computational graph by executing tensor operations such as
    /// contractions, additions, and multiplications. The execution can be controlled
    /// by mode and step limits.
    ///
    /// Parameters
    /// ----------
    /// library : TensorLibrary, optional
    ///     Optional tensor library for resolving tensor operations
    /// function_library : TensorFunctionLibrary or None, optional
    ///     Tensor function callbacks; None uses the built-in function library
    /// n_steps : int, optional
    ///     Maximum number of execution steps (None for complete execution)
    /// mode : ExecutionMode, optional
    ///     Execution strategy. ExecutionMode.Single selects one smallest-degree rewrite per
    ///     step; use n_steps to bound how many steps run.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.spenso import TensorNetwork, ExecutionMode, TensorLibrary
    /// >>> network = TensorNetwork(some_expression)
    /// >>> network.execute()
    /// >>> network.execute(n_steps=5)
    /// >>> network.execute(mode=ExecutionMode.Scalar)
    /// >>> lib = TensorLibrary.hep_lib()
    /// >>> network.execute(library=lib)
    #[pyo3(signature = (library=None,function_library=None, n_steps=None, mode=ExecutionMode::All))]
    pub(crate) fn execute(
        &mut self,
        library: Option<&SpensorLibrary>,
        function_library: Option<&SpensorFunctionLibrary>,
        n_steps: Option<usize>,
        mode: ExecutionMode,
    ) -> PyResult<()> {
        let lib = library.map(|l| &l.library).unwrap_or(HEP_LIB.deref());
        let fn_lib = function_library
            .map(|l| &l.library)
            .unwrap_or(FUN_LIB.deref());

        if let Some(n) = n_steps {
            for _ in 0..n {
                match mode {
                    ExecutionMode::All => {
                        self.network
                            .execute::<Steps<1>, SmallestDegree, _, _, _>(lib, fn_lib)
                            .map_err(|a| PyRuntimeError::new_err(a.to_string()))?;
                    }
                    ExecutionMode::Scalar => {
                        self.network
                            .execute::<Steps<1>, ContractScalars, _, _, _>(lib, fn_lib)
                            .map_err(|a| PyRuntimeError::new_err(a.to_string()))?;
                    }
                    ExecutionMode::Single => {
                        self.network
                            .execute::<Steps<1>, SingleSmallestDegree<false>, _, _, _>(lib, fn_lib)
                            .map_err(|a| PyRuntimeError::new_err(a.to_string()))?;
                    }
                }
            }
        } else {
            match mode {
                ExecutionMode::All => {
                    self.network
                        .execute::<Sequential, SmallestDegree, _, _, _>(lib, fn_lib)
                        .map_err(|a| PyRuntimeError::new_err(a.to_string()))?;
                }
                ExecutionMode::Scalar => {
                    self.network
                        .execute::<Sequential, ContractScalars, _, _, _>(lib, fn_lib)
                        .map_err(|a| PyRuntimeError::new_err(a.to_string()))?;
                }
                ExecutionMode::Single => {
                    self.network
                        .execute::<Sequential, SingleSmallestDegree<false>, _, _, _>(lib, fn_lib)
                        .map_err(|a| PyRuntimeError::new_err(a.to_string()))?;
                }
            }
        }
        Ok(())
    }
    /// Extract the final tensor result from the executed network.
    ///
    /// After network execution, retrieves the computed tensor result. The network
    /// should be executed before calling this method.
    ///
    /// Parameters
    /// ----------
    /// library : TensorLibrary, optional
    ///     Optional tensor library for resolving tensor structures
    ///
    /// Returns
    /// -------
    /// Tensor
    ///     The computed tensor result
    ///
    /// Raises
    /// ------
    /// RuntimeError
    ///     If the network execution resulted in an error
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.spenso import TensorNetwork, TensorLibrary
    /// >>> network = TensorNetwork(tensor_expression)
    /// >>> network.execute()
    /// >>> result = network.result_tensor()
    /// >>> lib = TensorLibrary.hep_lib()
    /// >>> result_with_lib = network.result_tensor(library=lib)
    #[pyo3(signature = (library=None))]
    pub(crate) fn result_tensor(&self, library: Option<&SpensorLibrary>) -> PyResult<Spensor> {
        let lib = library.map(|l| &l.library).unwrap_or(HEP_LIB.deref());
        let descriptor = self.structure.clone();
        let (name, args) = self
            .descriptor
            .clone()
            .map(|(name, args)| (Some(name), args))
            .unwrap_or_default();

        Ok(
            match self
                .network
                .result_tensor(lib)
                .map_err(|s| PyRuntimeError::new_err(s.to_string()))?
            {
                ExecutionResult::One => Spensor::scalar_with_descriptor(1., descriptor, name, args),
                ExecutionResult::Zero => {
                    Spensor::scalar_with_descriptor(0., descriptor, name, args)
                }
                ExecutionResult::Val(value) => {
                    let value = value.into_owned();
                    let value = match value.clone().scalar() {
                        Some(ConcreteOrParam::Param(atom)) => {
                            match SymComplex::<f64>::try_from(&atom) {
                                Ok(number) if number.im == 0. => ParamOrConcrete::new_scalar(
                                    ConcreteOrParam::Concrete(RealOrComplex::Real(number.re)),
                                ),
                                Ok(number) => {
                                    ParamOrConcrete::new_scalar(ConcreteOrParam::Concrete(
                                        RealOrComplex::Complex(Complex::new(number.re, number.im)),
                                    ))
                                }
                                Err(_) => value,
                            }
                        }
                        _ => value,
                    };
                    // Network tensors are already stored in canonical Spenso axis
                    // order. The descriptor permutations only record how that
                    // storage maps back to the public logical interface.
                    Spensor::from_storage_with_descriptor(value, descriptor, name, args)
                }
            },
        )
    }

    /// Extract the final scalar result from the executed network.
    ///
    /// For networks that evaluate to scalar expressions, retrieves the computed
    /// scalar value. The network should be executed before calling this method.
    ///
    /// Returns
    /// -------
    /// Expression
    ///     The computed scalar expression
    ///
    /// Raises
    /// ------
    /// RuntimeError
    ///     If the network execution resulted in an error
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.spenso import TensorNetwork
    /// >>> network = TensorNetwork(scalar_expression)
    /// >>> network.execute()
    /// >>> scalar_result = network.result_scalar()
    fn result_scalar(&self) -> PyResult<PythonExpression> {
        Ok(
            match self
                .network
                .result_scalar()
                .map_err(|s| PyRuntimeError::new_err(s.to_string()))?
            {
                ExecutionResult::One => Atom::num(1).into(),
                ExecutionResult::Zero => Atom::Zero.into(),
                ExecutionResult::Val(v) => v.into_owned().into(),
            },
        )
    }

    fn __repr__(&self) -> String {
        format!(
            "TensorNetwork({})",
            display::format_structured(&self.structure, false)
        )
    }

    /// Format the semantic tensor expression; use `to_dot()` for the executable graph.
    fn __str__(&self) -> String {
        display::format_structured(&self.structure, false)
    }

    fn _repr_latex_(&self) -> String {
        display::structured_to_latex(&self.structure, false, None)
    }

    fn _repr_pretty_(&self, pretty: &Bound<'_, PyAny>, cycle: bool) -> PyResult<()> {
        pretty.call_method1("text", (if cycle { "...".into() } else { self.__str__() },))?;
        Ok(())
    }

    /// Return the computational graph in Graphviz DOT format.
    fn to_dot(&self) -> String {
        self.network.dot_pretty()
    }

    /// Return the exact Linnest/Typst entrypoint for the current executable graph.
    /// Uses the same ``linnet.RenderConfig`` and asset pipeline as Feynman diagrams.
    #[pyo3(signature = (*, config=None))]
    fn to_linnest(
        &self,
        py: Python<'_>,
        #[gen_stub(override_type(type_repr="linnet.RenderConfig | None", imports=("linnet")))]
        config: Option<&Bound<'_, PyAny>>,
    ) -> PyResult<String> {
        self.prepare_render(py, config)?
            .getattr("typst_source")?
            .extract()
    }

    /// Render the current graph as interactive SVG, using Linnest's network style.
    /// Operator nodes, stored tensors, library references, and scalars retain
    /// their identities. ``expression().to_svg()`` renders the source formula.
    #[pyo3(signature = (*, config=None))]
    fn render(
        &self,
        py: Python<'_>,
        #[gen_stub(override_type(type_repr="linnet.RenderConfig | None", imports=("linnet")))]
        config: Option<&Bound<'_, PyAny>>,
    ) -> PyResult<String> {
        let svg: String = self
            .prepare_render(py, config)?
            .call_method0("to_svg")?
            .extract()?;
        Ok(display::network::svg_theme(&svg))
    }

    /// Format the exact semantic structure using compact Spenso notation.
    #[pyo3(signature = (show_dimensions = None, *, settings = None))]
    fn format_tensor(
        &self,
        show_dimensions: Option<bool>,
        settings: Option<PyRef<'_, display::DisplaySettings>>,
    ) -> PyResult<String> {
        let settings = display::resolved_settings(show_dimensions, settings.as_deref());
        display::validate_plain_source_settings(&settings)?;
        Ok(display::format_structured_settings(
            &self.structure,
            &settings,
        ))
    }

    /// Format the semantic source structure as Typst math source.
    #[pyo3(signature = (show_dimensions = None, *, settings = None))]
    fn to_typst(
        &self,
        show_dimensions: Option<bool>,
        settings: Option<PyRef<'_, display::DisplaySettings>>,
    ) -> PyResult<String> {
        let settings = display::resolved_settings(show_dimensions, settings.as_deref());
        display::validate_typst_source_settings(&settings)?;
        Ok(display::structured_to_typst_with_settings(
            &self.structure,
            &settings,
        ))
    }

    /// Build Symbolica's rich display wrapper for the semantic source structure.
    #[pyo3(signature = (show_dimensions = None, *, settings = None, notation_source = None))]
    fn formatted(
        &self,
        py: Python<'_>,
        show_dimensions: Option<bool>,
        settings: Option<PyRef<'_, display::DisplaySettings>>,
        notation_source: Option<String>,
    ) -> PythonFormattedOutput {
        let settings = display::resolved_settings(show_dimensions, settings.as_deref());
        display::format_structured_output_rich(
            py,
            &self.structure,
            &settings,
            notation_source.as_deref(),
        )
    }

    /// Render the executable graph as a notebook figure. For the symbolic source
    /// formula, use ``expression().to_html(settings=...)``.
    #[pyo3(signature = (*, config=None))]
    fn to_html(
        &self,
        py: Python<'_>,
        #[gen_stub(override_type(type_repr="linnet.RenderConfig | None", imports=("linnet")))]
        config: Option<&Bound<'_, PyAny>>,
    ) -> PyResult<String> {
        Ok(display::network::html(
            &self.render(py, config)?,
            &self.status(),
        ))
    }

    #[pyo3(signature = (show_dimensions = None, *, settings = None, notation_source = None))]
    fn to_svg(
        &self,
        py: Python<'_>,
        show_dimensions: Option<bool>,
        settings: Option<PyRef<'_, display::DisplaySettings>>,
        notation_source: Option<String>,
    ) -> PyResult<String> {
        let settings = display::resolved_settings(show_dimensions, settings.as_deref());
        display::structured_to_svg(py, &self.structure, &settings, notation_source.as_deref())
    }

    fn _repr_html_(&self, py: Python<'_>) -> PyResult<String> {
        self.to_html(py, None)
    }

    /// Return the semantic source expression and its public tensor interface.
    ///
    /// This expression records the tensor-aware structure used for composition and
    /// provenance. It is not reconstructed from the current execution store, so
    /// `replace()`, `evaluate()`, and `execute()` leave it unchanged. Use
    /// `result_scalar()` or `result_tensor()` to inspect the current computed value.
    fn expression(&self, py: Python<'_>) -> PyResult<Py<TensorExpression>> {
        TensorExpression::from_known_parts(
            py,
            self.structure.expression().clone(),
            self.structure.structure().clone(),
            self.descriptor.as_ref().map(|(name, _)| *name),
            self.descriptor
                .as_ref()
                .map(|(_, args)| args.clone())
                .unwrap_or_default(),
        )
    }

    /// Source tensor identity and ordered external ports, independent of execution.
    #[getter]
    fn structure(&self) -> crate::metadata::SpensoTensorStructure {
        crate::metadata::SpensoTensorStructure {
            interface: self.structure.structure().clone(),
            name: self.descriptor.as_ref().map(|(name, _)| *name),
            arguments: self
                .descriptor
                .as_ref()
                .map(|(_, args)| args.clone())
                .unwrap_or_default(),
        }
    }

    /// Fill the unresolved external ports with `indices` in interface order.
    ///
    /// Pass `AUTO` to leave a port unresolved. Pass `cook_indices=CookSettings.indices()` to flatten nested
    /// symbolic index payloads before insertion.
    #[pyo3(signature = (*indices, cook_indices = None))]
    fn index(
        &self,
        indices: &Bound<'_, PyTuple>,
        cook_indices: Option<&crate::simplification::PyCookSettings>,
    ) -> PyResult<Self> {
        self.index_network(indices, cook_indices)
    }

    /// Assign all external ports, retaining data and execution progress. AUTO leaves a port unchanged.
    #[pyo3(signature = (*indices, cook_indices=None))]
    pub(crate) fn reindex(
        &self,
        indices: &Bound<'_, PyTuple>,
        cook_indices: Option<&crate::simplification::PyCookSettings>,
    ) -> PyResult<Self> {
        let positions = (0..self.structure.rank()).collect::<Vec<_>>();
        let replacements = TensorExpression::port_replacements(
            self.structure.structure(),
            &positions,
            indices,
            cook_indices,
        )?;
        self.set_port_indices(&replacements)
    }

    /// Rename external indices simultaneously without changing rank or capturing dummy indices.
    #[pyo3(signature = (mapping, *, cook_indices=None))]
    pub(crate) fn rename_indices(
        &self,
        mapping: &Bound<'_, PyDict>,
        cook_indices: Option<&crate::simplification::PyCookSettings>,
    ) -> PyResult<Self> {
        let replacements = TensorExpression::named_replacements(
            self.structure.structure(),
            mapping,
            cook_indices,
        )?;
        let result = self.set_port_indices(&replacements)?;
        if result.structure.rank() != self.structure.rank() {
            return Err(PyValueError::new_err(
                "renaming would contract external ports; use reindex() to request a contraction",
            ));
        }
        Ok(result)
    }

    /// Return a copy with reordered external axes and the same execution progress.
    /// Unresolved ports acquire fresh identities; reindex() can assign preferred labels.
    pub(crate) fn permute_axes(&self, axes: Vec<usize>) -> PyResult<Self> {
        let structure = self.structure.permuted(&axes).map_err(composition_error)?;
        let new_slots = structure.structure().logical_slots();
        let old_slots = self.structure.structure().logical_slots();
        let replacements = axes
            .iter()
            .enumerate()
            .filter_map(|(new, &old)| {
                if matches!(old_slots[old].aind, PartialIndex::Open(_)) {
                    let PartialIndex::Explicit(index) = new_slots[new].aind else {
                        unreachable!()
                    };
                    Some((old, index))
                } else {
                    None
                }
            })
            .collect();
        let mut result = self.set_port_indices(&replacements)?;
        result.structure = structure;
        result.materialized = result
            .materialized
            .permuted(&axes)
            .map_err(composition_error)?;
        result.validate_graph_interface()?;
        Ok(result)
    }

    /// Copy the graph, values, and execution progress independently.
    fn __copy__(&self) -> Self {
        self.clone()
    }

    /// Fill the unresolved external ports with `indices` in interface order.
    #[pyo3(signature = (*indices, cook_indices = None))]
    fn __call__(
        &self,
        indices: &Bound<'_, PyTuple>,
        cook_indices: Option<&crate::simplification::PyCookSettings>,
    ) -> PyResult<Self> {
        self.index_network(indices, cook_indices)
    }

    pub fn __neg__(&self) -> PyResult<Self> {
        let structure = SymbolicTensor::new(
            Atom::num(-1) * self.structure.expression(),
            self.structure.structure().clone(),
        );
        let materialized = SymbolicTensor::new(
            Atom::num(-1) * self.materialized.expression(),
            self.materialized.structure().clone(),
        );
        Self::finish(-self.network.clone(), structure, materialized)
    }

    /// Add two tensor networks element-wise.
    ///
    /// Parameters
    /// ----------
    /// rhs : TensorNetwork
    ///     The tensor network to add (right-hand side)
    ///
    /// Returns
    /// -------
    /// TensorNetwork
    ///     A new TensorNetwork representing the sum
    ///
    /// Examples
    /// --------
    /// >>> net1 = TensorNetwork(expr1)
    /// >>> net2 = TensorNetwork(expr2)
    /// >>> sum_net = net1 + net2
    pub fn __add__(&self, rhs: ConvertibleToSpensoNet) -> PyResult<SpensoNet> {
        self.clone().add_network(rhs.to_net(), false)
    }

    /// Add two tensor networks element-wise (right-hand addition).
    pub fn __radd__(&self, rhs: ConvertibleToSpensoNet) -> PyResult<SpensoNet> {
        self.__add__(rhs)
    }

    /// Subtract one tensor network from another element-wise.
    ///
    /// Parameters
    /// ----------
    /// rhs : TensorNetwork
    ///     The tensor network to subtract (right-hand side)
    ///
    /// Returns
    /// -------
    /// TensorNetwork
    ///     A new TensorNetwork representing the difference
    ///
    /// Examples
    /// --------
    /// >>> net1 = TensorNetwork(expr1)
    /// >>> net2 = TensorNetwork(expr2)
    /// >>> diff_net = net1 - net2
    pub fn __sub__(&self, rhs: ConvertibleToSpensoNet) -> PyResult<SpensoNet> {
        self.clone().add_network(rhs.to_net(), true)
    }

    /// Subtract one tensor network from another (right-hand subtraction).
    pub fn __rsub__(&self, rhs: ConvertibleToSpensoNet) -> PyResult<SpensoNet> {
        rhs.to_net().add_network(self.clone(), true)
    }

    /// Multiply two tensor networks.
    ///
    /// Parameters
    /// ----------
    /// rhs : TensorNetwork
    ///     The tensor network to multiply with (right-hand side)
    ///
    /// Returns
    /// -------
    /// TensorNetwork
    ///     A new TensorNetwork representing the product
    ///
    /// Examples
    /// --------
    /// >>> net1 = TensorNetwork(expr1)
    /// >>> net2 = TensorNetwork(expr2)
    /// >>> product_net = net1 * net2
    pub fn __mul__(&self, rhs: ConvertibleToSpensoNet) -> PyResult<SpensoNet> {
        self.clone().multiply_network(rhs.to_net())
    }

    /// Multiply two tensor networks (right-hand multiplication).
    pub fn __rmul__(&self, rhs: ConvertibleToSpensoNet) -> PyResult<SpensoNet> {
        rhs.to_net().multiply_network(self.clone())
    }

    pub fn __truediv__(&self, rhs: ConvertibleToSpensoNet) -> PyResult<SpensoNet> {
        self.clone().multiply_network(rhs.to_net().reciprocal()?)
    }

    pub fn __rtruediv__(&self, lhs: ConvertibleToSpensoNet) -> PyResult<SpensoNet> {
        if !lhs.0.structure.is_scalar() {
            return Err(PyTypeError::new_err(
                "tensor division requires a scalar numerator",
            ));
        }
        lhs.to_net().multiply_network(self.clone().reciprocal()?)
    }

    /// Form an outer tensor product without contracting compatible ports.
    pub fn outer(&self, rhs: ConvertibleToSpensoNet) -> PyResult<Self> {
        self.clone().apply_product(rhs.to_net(), ProductPlan::Outer)
    }

    /// Contract one selected pair of public interface positions.
    #[pyo3(signature = (rhs, *, left, right))]
    pub fn contract(
        &self,
        rhs: ConvertibleToSpensoNet,
        left: usize,
        right: usize,
    ) -> PyResult<Self> {
        self.clone().apply_product(
            rhs.to_net(),
            ProductPlan::Contract(vec![composition::PortPair { left, right }]),
        )
    }

    /// Compose two explicitly selected matrix channels.
    #[pyo3(signature = (rhs, *, left, right))]
    pub fn compose(
        &self,
        rhs: ConvertibleToSpensoNet,
        left: (usize, usize),
        right: (usize, usize),
    ) -> PyResult<Self> {
        self.clone().apply_product(
            rhs.to_net(),
            ProductPlan::Compose(
                composition::MatrixChannel {
                    input: left.0,
                    output: left.1,
                },
                composition::MatrixChannel {
                    input: right.0,
                    output: right.1,
                },
            ),
        )
    }

    /// Contract two rank-one operands into the canonical dot form.
    pub fn dot(&self, rhs: ConvertibleToSpensoNet) -> PyResult<Self> {
        if self.structure.rank() != 1 || rhs.0.structure.rank() != 1 {
            return Err(PyValueError::new_err(format!(
                "dot() requires rank-one operands, got ranks {} and {}",
                self.structure.rank(),
                rhs.0.structure.rank()
            )));
        }
        self.contract(rhs, 0, 0)
    }

    /// Close a selected or uniquely inferred propagation channel.
    #[pyo3(signature = (*, channel = None))]
    pub fn trace(&self, channel: Option<(usize, usize)>) -> PyResult<Self> {
        let selected = match channel {
            Some((input, output)) => composition::MatrixChannel { input, output },
            None => self
                .structure
                .matrix_channel()
                .ok_or_else(|| PyValueError::new_err("tensor has no unique matrix channel"))?,
        };
        let structure = self
            .structure
            .trace_ports(selected)
            .map_err(composition_error)?;
        let indices = self
            .structure
            .trace_indices(selected)
            .map_err(composition_error)?;
        let mut network = self.clone().relabel_ports(&indices)?;
        let materialized = network
            .materialized
            .trace_ports(selected)
            .map_err(composition_error)?;
        network.network.graph.sew_dangling_slots();
        network.network.state = network.network.graph.state();
        Self::finish(network.network, structure, materialized)
    }

    // pub fn __pow__(&self, rhs: usize, number: Option<i64>) -> PyResult<PythonExpression> {
    //     if number.is_some() {
    //         return Err(exceptions::PyValueError::new_err(
    //             "Optional number argument not supported",
    //         ));
    //     }

    //     // let rhs = rhs.to_net();
    //     Ok(self.network.pow(&rhs).into())
    // }
}

#[cfg(feature = "python_stubgen")]
pyo3_stub_gen::define_stub_info_gatherer!(stub_info);

#[cfg(test)]
mod tests {
    use idenso::{dirac::AGS, representations::initialize};
    use spenso::network::parsing::ParseSettings;
    use spenso::structure::{
        OrderedStructure, TensorStructure,
        dimension::Dimension,
        partial::{PartialIndex, PartialStructure},
        representation::{ExtendibleReps, RepName},
    };
    use spenso_hep_lib::HEP_LIB;
    use symbolica::parse_lit;

    use super::*;

    fn data_tensor(descriptor: SymbolicTensor<PartialStructure>, name: Symbol) -> Spensor {
        let owner = AbstractIndex::fresh_open_owner();
        let storage = OrderedStructure::new(
            descriptor
                .structure()
                .logical_slots()
                .into_iter()
                .enumerate()
                .map(|(axis, slot)| slot.rep().slot(AbstractIndex::Open { owner, axis }))
                .collect(),
        )
        .map_canonical(|structure| ShadowedStructure {
            structure,
            global_name: Some(name),
            additional_args: None,
        })
        .map_canonical(|structure| SparseTensor::<f64, _>::empty(structure, 0.0).into())
        .into_canonical();
        Spensor::from_storage_with_descriptor(storage, descriptor, Some(name), Vec::new())
    }

    fn tensor_descriptor(
        name: &str,
        slots: impl IntoIterator<Item = spenso::structure::partial::PartialSlot>,
    ) -> (SymbolicTensor<PartialStructure>, Symbol) {
        let name = SPENSO_TAG.tensor_symbol(name);
        let interface = PartialStructure::from_logical_slots(slots);
        let atom = FunctionBuilder::new(name)
            .add_args(
                interface
                    .logical_slots()
                    .into_iter()
                    .map(composition::port_atom),
            )
            .finish();
        (SymbolicTensor::new(atom, interface), name)
    }

    #[test]
    fn ordinary_tensor_expression_infers_interface_from_graph() {
        initialize();
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let representation = ExtendibleReps::EUCLIDEAN.new_rep(Dimension::Concrete(2));
            let slot: spenso::structure::representation::LibrarySlot<AbstractIndex> =
                representation.slot(AbstractIndex::Normal(71));
            let atom = FunctionBuilder::new(
                SPENSO_TAG.tensor_symbol("ordinary_network_interface_inference"),
            )
            .add_arg(slot.to_atom())
            .finish();
            let expression = PythonExpression { expr: atom }.into_pyobject(py)?;

            let network = SpensoNet::from_expression(expression.as_any(), None)?;

            assert_eq!(network.structure.rank(), 1);
            assert_eq!(
                network.structure.structure().logical_slots()[0].aind,
                PartialIndex::Explicit(slot.aind)
            );
            assert_eq!(network.network.graph.dangling_indices(), vec![slot]);
            Ok(())
        })
        .unwrap();
    }

    #[test]
    fn implicit_network_conversion_uses_default_hep_library() {
        initialize();
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let gamma5 = AGS.gamma5_strct::<AbstractIndex>(4);
            let expression = TensorExpression::from_structure(py, &gamma5)?;

            let mut network = expression
                .bind(py)
                .as_any()
                .extract::<ConvertibleToSpensoNet>()?
                .to_net();

            network.execute(None, None, None, ExecutionMode::All)?;
            let result = network.result_tensor(None)?;
            assert!(matches!(&result.tensor, ParamOrConcrete::Concrete(_)));
            Ok(())
        })
        .unwrap();
    }

    #[test]
    fn repeated_indices_contract_tensor_and_network_ports() {
        initialize();
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let representation = ExtendibleReps::EUCLIDEAN.new_rep(Dimension::Concrete(2));
            let (descriptor, name) = tensor_descriptor(
                "network_repeated_index_trace",
                [
                    representation.slot(PartialIndex::open(0)),
                    representation.slot(PartialIndex::open(1)),
                ],
            );
            let expression = TensorExpression::from_known_parts(
                py,
                descriptor.expression().clone(),
                descriptor.structure().clone(),
                Some(name),
                Vec::new(),
            )?;

            let indexed_expression = expression.bind(py).call1(("i", "i"))?;
            let indexed_expression = indexed_expression.extract::<PyRef<'_, TensorExpression>>()?;
            assert!(indexed_expression.interface().canonical().is_scalar());
            drop(indexed_expression);

            let tensor = Spensor::dense(
                expression.bind(py).as_any().extract()?,
                crate::AtomsOrFloats::Floats(vec![1.0, 2.0, 3.0, 4.0]),
            )?;
            let relabeled_tensor = Py::new(py, tensor.clone())?
                .bind(py)
                .call1(("i", "j"))?
                .extract::<SpensoNet>()?;
            assert_eq!(relabeled_tensor.structure.rank(), 2);
            assert_eq!(relabeled_tensor.materialized.rank(), 2);
            assert_eq!(relabeled_tensor.network.graph.dangling_indices().len(), 2);
            assert!(relabeled_tensor.descriptor.is_some());

            let indexed_tensor = Py::new(py, tensor.clone())?
                .bind(py)
                .call1(("i", "i"))?
                .extract::<SpensoNet>()?;
            let indexed_network = Py::new(py, SpensoNet::from_tensor(tensor)?)?
                .bind(py)
                .call1(("i", "i"))?
                .extract::<SpensoNet>()?;

            for mut network in [indexed_tensor, indexed_network] {
                assert!(network.structure.is_scalar());
                assert!(network.materialized.is_scalar());
                assert!(network.network.graph.dangling_indices().is_empty());
                assert!(network.network.state.is_scalar());
                assert!(network.descriptor.is_none());

                network.execute(None, None, None, ExecutionMode::All)?;
                let result = SymComplex::<f64>::try_from(&network.result_scalar()?.expr).unwrap();
                assert_eq!(result.re, 5.0);
                assert_eq!(result.im, 0.0);
            }
            Ok(())
        })
        .unwrap();
    }

    #[test]
    fn equal_network_chain_endpoints_normalize_to_a_root_trace() {
        initialize();
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let representation = ExtendibleReps::EUCLIDEAN.new_rep(Dimension::Concrete(2));
            let factor = FunctionBuilder::new(SPENSO_TAG.tensor_symbol("network_chain_factor"))
                .add_arg(Atom::var(SPENSO_TAG.chain_in))
                .add_arg(Atom::var(SPENSO_TAG.chain_out))
                .finish();
            let expression = TensorExpression::from_structured(
                py,
                SymbolicTensor::new(
                    SPENSO_TAG.chain(
                        representation.to_symbolic([]),
                        representation.to_symbolic([]),
                        [factor],
                    ),
                    PartialStructure::from_logical_slots([
                        representation.slot(PartialIndex::open(0)),
                        representation.slot(PartialIndex::open(1)),
                    ]),
                ),
            )?;
            let network = SpensoNet::from_arithmetic(
                ArithmeticStructure::Tensor(expression),
                None,
            )?;
            let index = AbstractIndex::Normal(73);
            let closed = network.set_port_indices(&HashMap::from([(0, index), (1, index)]))?;

            assert!(closed.structure.is_scalar());
            assert!(closed.materialized.is_scalar());
            assert!(
                matches!(closed.structure.expression().as_view(), AtomView::Fun(function) if function.get_symbol() == SPENSO_TAG.trace)
            );
            assert!(
                matches!(closed.materialized.expression().as_view(), AtomView::Fun(function) if function.get_symbol() == SPENSO_TAG.trace)
            );
            assert!(closed.network.graph.dangling_indices().is_empty());
            assert!(closed.network.state.is_scalar());
            Ok(())
        })
        .unwrap();
    }

    #[test]
    fn lazy_library_indices_follow_network_port_relabeling() {
        initialize();
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let representation = ExtendibleReps::EUCLIDEAN.new_rep(Dimension::Concrete(2));
            let (descriptor, name) = tensor_descriptor(
                "lazy_library_relabeling",
                [
                    representation.slot(PartialIndex::open(0)),
                    representation.slot(PartialIndex::open(1)),
                ],
            );
            let expression = TensorExpression::from_known_parts(
                py,
                descriptor.expression().clone(),
                descriptor.structure().clone(),
                Some(name),
                Vec::new(),
            )?;
            let tensor = Py::new(
                py,
                Spensor::dense(
                    expression.bind(py).as_any().extract()?,
                    crate::AtomsOrFloats::Floats(vec![1., 2., 3., 4.]),
                )?,
            )?;
            let mut library = SpensorLibrary::new();
            library.register(tensor.bind(py).borrow())?;

            let network = SpensoNet::from_expression(expression.bind(py).as_any(), Some(&library))?;
            assert!(network.network.store.tensors.is_empty());
            let mut network = Py::new(py, network)?
                .bind(py)
                .call1(("i", "i"))?
                .extract::<SpensoNet>()?;
            assert!(network.network.graph.dangling_indices().is_empty());
            assert!(network.network.state.is_scalar());

            network.execute(Some(&library), None, None, ExecutionMode::All)?;
            let result = SymComplex::<f64>::try_from(&network.result_scalar()?.expr).unwrap();
            assert_eq!(result.re, 5.);
            assert_eq!(result.im, 0.);
            Ok(())
        })
        .unwrap();
    }

    #[test]
    fn lazy_library_indices_follow_cross_node_alignment() {
        initialize();
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let representation = ExtendibleReps::EUCLIDEAN.new_rep(Dimension::Concrete(2));
            let (matrix_descriptor, matrix_name) = tensor_descriptor(
                "lazy_library_matrix",
                [
                    representation.slot(PartialIndex::open(0)),
                    representation.slot(PartialIndex::open(1)),
                ],
            );
            let (vector_descriptor, vector_name) = tensor_descriptor(
                "lazy_library_vector",
                [representation.slot(PartialIndex::open(0))],
            );
            let matrix = TensorExpression::from_known_parts(
                py,
                matrix_descriptor.expression().clone(),
                matrix_descriptor.structure().clone(),
                Some(matrix_name),
                Vec::new(),
            )?;
            let vector = TensorExpression::from_known_parts(
                py,
                vector_descriptor.expression().clone(),
                vector_descriptor.structure().clone(),
                Some(vector_name),
                Vec::new(),
            )?;
            let matrix_data = Py::new(
                py,
                Spensor::dense(
                    matrix.bind(py).as_any().extract()?,
                    crate::AtomsOrFloats::Floats(vec![1., 2., 3., 4.]),
                )?,
            )?;
            let vector_data = Py::new(
                py,
                Spensor::dense(
                    vector.bind(py).as_any().extract()?,
                    crate::AtomsOrFloats::Floats(vec![10., 100.]),
                )?,
            )?;
            let mut library = SpensorLibrary::new();
            library.register(matrix_data.bind(py).borrow())?;
            library.register(vector_data.bind(py).borrow())?;

            let matrix = SpensoNet::from_expression(matrix.bind(py).as_any(), Some(&library))?;
            let vector = SpensoNet::from_expression(vector.bind(py).as_any(), Some(&library))?;
            assert!(matrix.network.store.tensors.is_empty());
            assert!(vector.network.store.tensors.is_empty());
            let mut product = matrix.contract(ConvertibleToSpensoNet(vector), 1, 0)?;

            product.execute(Some(&library), None, None, ExecutionMode::All)?;
            let result = Py::new(py, product.result_tensor(Some(&library))?)?;
            assert_eq!(result.bind(py).get_item(0)?.extract::<f64>()?, 210.);
            assert_eq!(result.bind(py).get_item(1)?.extract::<f64>()?, 430.);
            Ok(())
        })
        .unwrap();
    }

    #[test]
    fn repeated_indices_preserve_remaining_port_order() {
        initialize();
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let representation = ExtendibleReps::EUCLIDEAN.new_rep(Dimension::Concrete(2));
            let (descriptor, name) = tensor_descriptor(
                "network_partial_repeated_index_trace",
                [
                    representation.slot(PartialIndex::open(0)),
                    representation.slot(PartialIndex::open(1)),
                    representation.slot(PartialIndex::open(2)),
                ],
            );
            let expression = TensorExpression::from_known_parts(
                py,
                descriptor.expression().clone(),
                descriptor.structure().clone(),
                Some(name),
                Vec::new(),
            )?;
            let tensor = Spensor::dense(
                expression.bind(py).as_any().extract()?,
                crate::AtomsOrFloats::Floats((1..=8).map(f64::from).collect()),
            )?;
            let mut network = Py::new(py, tensor)?
                .bind(py)
                .call1(("i", "i", "j"))?
                .extract::<SpensoNet>()?;

            assert_eq!(network.structure.rank(), 1);
            assert_eq!(network.materialized.rank(), 1);
            assert_eq!(network.network.graph.dangling_indices().len(), 1);

            network.execute(None, None, None, ExecutionMode::All)?;
            let result = Py::new(py, network.result_tensor(None)?)?;
            assert_eq!(result.bind(py).get_item(0)?.extract::<f64>()?, 8.0);
            assert_eq!(result.bind(py).get_item(1)?.extract::<f64>()?, 10.0);
            Ok(())
        })
        .unwrap();
    }

    #[test]
    fn additive_operations_materialize_permuted_library_occurrences_in_canonical_slot_order() {
        initialize();
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let representation = ExtendibleReps::EUCLIDEAN.new_rep(Dimension::Concrete(2));
            let (descriptor, name) = tensor_descriptor(
                "network_canonical_interface_addition",
                [
                    representation.slot(PartialIndex::open(0)),
                    representation.slot(PartialIndex::open(1)),
                ],
            );
            let expression = TensorExpression::from_known_parts(
                py,
                descriptor.expression().clone(),
                descriptor.structure().clone(),
                Some(name),
                Vec::new(),
            )?;
            let tensor = Py::new(
                py,
                Spensor::dense(
                    expression.bind(py).as_any().extract()?,
                    crate::AtomsOrFloats::Floats(vec![1., 2., 3., 4.]),
                )?,
            )?;
            let mut library = SpensorLibrary::new();
            library.register(tensor.bind(py).borrow())?;

            let left = expression
                .bind(py)
                .call1(("i", "j"))?
                .extract::<Py<TensorExpression>>()?;
            let right = expression
                .bind(py)
                .call1(("j", "i"))?
                .extract::<Py<TensorExpression>>()?;
            let expected_interface = left.bind(py).borrow().interface().logical_slots();

            let sum = left
                .bind(py)
                .call_method1("__add__", (right.clone_ref(py),))?;
            let difference = left
                .bind(py)
                .call_method1("__sub__", (right.clone_ref(py),))?;
            for value in [&sum, &difference] {
                let expression = value.extract::<PyRef<'_, TensorExpression>>()?;
                assert_eq!(expression.interface().logical_slots(), expected_interface);
            }

            let assert_result = |mut network: SpensoNet, expected: &[f64]| -> PyResult<()> {
                assert_eq!(
                    network.structure.structure().logical_slots(),
                    expected_interface
                );
                network.execute(Some(&library), None, None, ExecutionMode::All)?;
                let result = Py::new(py, network.result_tensor(Some(&library))?)?;
                let values = (0..4)
                    .map(|index| result.bind(py).get_item(index)?.extract::<f64>())
                    .collect::<PyResult<Vec<_>>>()?;
                assert_eq!(values, expected);
                Ok(())
            };

            assert_result(
                SpensoNet::from_expression(sum.as_any(), Some(&library))?,
                &[2., 5., 5., 8.],
            )?;
            assert_result(
                SpensoNet::from_expression(difference.as_any(), Some(&library))?,
                &[0., -1., 1., 0.],
            )?;
            assert_result(
                SpensoNet::from_expression(left.bind(py).as_any(), Some(&library))?.add_network(
                    SpensoNet::from_expression(right.bind(py).as_any(), Some(&library))?,
                    false,
                )?,
                &[2., 5., 5., 8.],
            )?;
            assert_result(
                SpensoNet::from_expression(left.bind(py).as_any(), Some(&library))?.add_network(
                    SpensoNet::from_expression(right.bind(py).as_any(), Some(&library))?,
                    true,
                )?,
                &[0., -1., 1., 0.],
            )?;
            Ok(())
        })
        .unwrap();
    }

    #[test]
    fn replace_honors_pattern_restrictions() {
        initialize();
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let function = symbol!("network_conditional_replace_f");
            let wildcard = symbol!("network_conditional_replace_x_");
            let retained = function!(function, Atom::num(1));
            let replaced = function!(function, Atom::num(2));
            let expression = PythonExpression {
                expr: &retained + replaced,
            }
            .into_pyobject(py)?;
            let network = Py::new(py, SpensoNet::from_expression(expression.as_any(), None)?)?;
            let pattern = PythonExpression {
                expr: function!(function, Atom::var(wildcard)),
            }
            .into_pyobject(py)?;
            let wildcard = PythonExpression {
                expr: Atom::var(wildcard),
            }
            .into_pyobject(py)?;
            let condition = wildcard.call_method1("req_gt", (1,))?;

            let replaced = network
                .bind(py)
                .call_method1("replace", (pattern.as_any(), 0, condition.as_any()))?;
            let mut replaced = replaced.extract::<SpensoNet>()?;
            replaced.execute(None, None, None, ExecutionMode::All)?;

            assert_eq!(replaced.result_scalar()?.expr, retained);
            Ok(())
        })
        .unwrap();
    }

    #[test]
    fn test_parse() {
        initialize();
        let expr = parse_lit!(
            (-1 * gammalooprs::mUV
                ^ 2 + gammalooprs::Q(6, spenso::mink(4, gammalooprs::uv_mink_1337))
                    * gammalooprs::Q(7, spenso::mink(4, gammalooprs::uv_mink_1337)))
                * 2
        );

        // Use the native Rust API directly to avoid Python linking issues
        let value = SymbolicTensor::new(
            expr.clone(),
            PartialStructure::from_logical_slots(std::iter::empty()),
        );
        let net = SpensoNet {
            network: ParsingNet::try_from_view(
                expr.as_view(),
                &*HEP_LIB,
                &ParseSettings::default(),
            )
            .unwrap(),
            structure: value.clone(),
            materialized: value,
            descriptor: None,
        };

        println!("{}", net.network.dot_pretty())
    }

    #[test]
    fn non_scalar_zero_preserves_and_materializes_external_interface() {
        initialize();
        let representation = ExtendibleReps::EUCLIDEAN.new_rep(Dimension::Concrete(4));
        let explicit = AbstractIndex::Normal(37);
        let value = SymbolicTensor::new(
            Atom::Zero,
            PartialStructure::from_logical_slots([
                representation.slot(PartialIndex::open(0)),
                representation.slot(PartialIndex::Explicit(explicit)),
            ]),
        );

        let network = non_scalar_zero_network(&value).unwrap();
        assert!(network.state.is_tensor());
        assert_eq!(network.store.tensors.len(), 1);
        let slots = network.store.tensors[0].external_structure();
        assert_eq!(slots.len(), 2);
        assert!(slots.iter().any(|slot| slot.aind == explicit));
        assert_eq!(network.graph.dangling_indices().len(), 2);
    }

    #[test]
    fn open_tensor_product_contracts_only_the_planned_pair() {
        initialize();
        let euclidean = ExtendibleReps::EUCLIDEAN.new_rep(Dimension::Concrete(2));
        let minkowski = ExtendibleReps::MINKOWSKI.new_rep(Dimension::Concrete(2));
        let (left, left_name) = tensor_descriptor(
            "network_open_mixed_left",
            [
                euclidean.slot(PartialIndex::open(0)),
                minkowski.slot(PartialIndex::Explicit(AbstractIndex::Normal(41))),
            ],
        );
        let (right, right_name) = tensor_descriptor(
            "network_open_mixed_right",
            [euclidean.slot(PartialIndex::open(0))],
        );

        let product = SpensoNet::from_tensor(data_tensor(left, left_name))
            .unwrap()
            .multiply_network(SpensoNet::from_tensor(data_tensor(right, right_name)).unwrap())
            .unwrap();

        assert_eq!(product.structure.rank(), 1);
        assert_eq!(product.materialized.rank(), 1);
        assert_eq!(product.network.graph.dangling_indices().len(), 1);
        assert_eq!(
            product.structure.structure().logical_slots()[0].aind,
            PartialIndex::Explicit(AbstractIndex::Normal(41))
        );
    }

    #[test]
    fn tensor_product_contracts_permuted_explicit_indices() {
        initialize();
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let representation = ExtendibleReps::EUCLIDEAN.new_rep(Dimension::Concrete(2));
            let mut operands = Vec::new();
            for (name, indices) in [
                ("network_explicit_product_left", [71, 73, 79]),
                ("network_explicit_product_right", [73, 79, 71]),
            ] {
                let (descriptor, name) = tensor_descriptor(
                    name,
                    indices.map(|index| {
                        representation.slot(PartialIndex::Explicit(AbstractIndex::Normal(index)))
                    }),
                );
                let expression = TensorExpression::from_known_parts(
                    py,
                    descriptor.expression().clone(),
                    descriptor.structure().clone(),
                    Some(name),
                    Vec::new(),
                )?;
                let tensor = Spensor::dense(
                    expression.bind(py).as_any().extract()?,
                    crate::AtomsOrFloats::Floats((1..=8).map(f64::from).collect()),
                )?;
                operands.push(SpensoNet::from_tensor(tensor)?);
            }
            let mut product = operands.remove(0).multiply_network(operands.remove(0))?;

            assert!(product.structure.is_scalar());
            assert!(product.materialized.is_scalar());
            assert!(product.network.graph.dangling_indices().is_empty());
            product.execute(None, None, None, ExecutionMode::All)?;
            let result = SymComplex::<f64>::try_from(&product.result_scalar()?.expr).unwrap();
            // Sum T[i,j,k] U[j,k,i] with both arrays containing 1 through 8.
            assert_eq!(result.re, 190.0);
            assert_eq!(result.im, 0.0);
            Ok(())
        })
        .unwrap();
    }

    #[test]
    fn tensor_product_aligns_each_unique_open_representation_pair() {
        initialize();
        let euclidean = ExtendibleReps::EUCLIDEAN.new_rep(Dimension::Concrete(2));
        let minkowski = ExtendibleReps::MINKOWSKI.new_rep(Dimension::Concrete(2));
        let (left, left_name) = tensor_descriptor(
            "network_disjoint_open_left",
            [
                euclidean.slot(PartialIndex::open(0)),
                minkowski.slot(PartialIndex::open(1)),
            ],
        );
        let (right, right_name) = tensor_descriptor(
            "network_disjoint_open_right",
            [
                minkowski.slot(PartialIndex::open(0)),
                euclidean.slot(PartialIndex::open(1)),
            ],
        );

        let product = SpensoNet::from_tensor(data_tensor(left, left_name))
            .unwrap()
            .multiply_network(SpensoNet::from_tensor(data_tensor(right, right_name)).unwrap())
            .unwrap();

        assert!(product.structure.is_scalar());
        assert!(product.materialized.is_scalar());
        assert!(product.network.graph.dangling_indices().is_empty());
    }

    #[test]
    fn sequential_explicit_relabels_do_not_alias_storage_ports() {
        initialize();
        let representation = ExtendibleReps::EUCLIDEAN.new_rep(Dimension::Concrete(2));
        let first = AbstractIndex::Normal(1);
        let second = AbstractIndex::Normal(2);
        let (descriptor, name) = tensor_descriptor(
            "network_sequential_explicit_relabels",
            [
                representation.slot(PartialIndex::Explicit(first)),
                representation.slot(PartialIndex::Explicit(second)),
            ],
        );

        let network = SpensoNet::from_tensor(data_tensor(descriptor, name)).unwrap();
        let mut dangling = network.network.graph.dangling_indices();
        dangling.sort();

        assert_eq!(
            dangling
                .into_iter()
                .map(|slot| slot.aind())
                .collect::<Vec<_>>(),
            [first, second]
        );
    }

    #[test]
    fn explicit_port_swaps_use_collision_free_relabeling() {
        initialize();
        let representation = ExtendibleReps::EUCLIDEAN.new_rep(Dimension::Concrete(2));
        let first = AbstractIndex::Normal(61);
        let second = AbstractIndex::Normal(62);
        let (descriptor, name) = tensor_descriptor(
            "network_explicit_port_swap",
            [
                representation.slot(PartialIndex::Explicit(first)),
                representation.slot(PartialIndex::Explicit(second)),
            ],
        );
        let network = SpensoNet::from_tensor(data_tensor(descriptor, name)).unwrap();

        let swapped = network
            .set_port_indices(&HashMap::from([(0, second), (1, first)]))
            .unwrap();
        assert_eq!(
            swapped
                .structure
                .structure()
                .logical_slots()
                .into_iter()
                .map(|slot| slot.aind)
                .collect::<Vec<_>>(),
            [
                PartialIndex::Explicit(second),
                PartialIndex::Explicit(first),
            ]
        );
        assert_eq!(swapped.network.graph.dangling_indices().len(), 2);
        assert!(swapped.descriptor.is_some());
    }

    #[test]
    fn open_outer_product_does_not_alias_explicit_dummy_symbol() {
        initialize();
        let representation = ExtendibleReps::EUCLIDEAN.new_rep(Dimension::Concrete(2));
        let explicit = AbstractIndex::Symbol(symbol!("d_0").into());
        let (left, left_name) = tensor_descriptor(
            "network_open_outer_left",
            [representation.slot(PartialIndex::open(0))],
        );
        let (right, right_name) = tensor_descriptor(
            "network_explicit_dummy_outer_right",
            [representation.slot(PartialIndex::Explicit(explicit))],
        );
        let left = SpensoNet::from_tensor(data_tensor(left, left_name)).unwrap();
        let right = SpensoNet::from_tensor(data_tensor(right, right_name)).unwrap();

        let product = left.outer(ConvertibleToSpensoNet(right)).unwrap();
        let dangling = product.network.graph.dangling_indices();

        assert_eq!(product.structure.rank(), 2);
        assert_eq!(dangling.len(), 2);
        assert!(dangling.iter().any(|slot| slot.aind() == explicit));
        assert!(
            dangling
                .iter()
                .any(|slot| matches!(slot.aind(), AbstractIndex::Open { .. }))
        );
    }

    #[test]
    fn unequal_explicit_tensor_ports_remain_an_outer_product() {
        initialize();
        let representation = ExtendibleReps::EUCLIDEAN.new_rep(Dimension::Concrete(2));
        let (left, left_name) = tensor_descriptor(
            "network_explicit_outer_left",
            [representation.slot(PartialIndex::Explicit(AbstractIndex::Normal(51)))],
        );
        let (right, right_name) = tensor_descriptor(
            "network_explicit_outer_right",
            [representation.slot(PartialIndex::Explicit(AbstractIndex::Normal(53)))],
        );

        let product = SpensoNet::from_tensor(data_tensor(left, left_name))
            .unwrap()
            .multiply_network(SpensoNet::from_tensor(data_tensor(right, right_name)).unwrap())
            .unwrap();

        assert_eq!(product.structure.rank(), 2);
        assert_eq!(product.materialized.rank(), 2);
        assert_eq!(product.network.graph.dangling_indices().len(), 2);
        assert_eq!(
            product
                .structure
                .structure()
                .logical_slots()
                .iter()
                .map(|slot| slot.aind)
                .collect::<Vec<_>>(),
            [
                PartialIndex::Explicit(AbstractIndex::Normal(51)),
                PartialIndex::Explicit(AbstractIndex::Normal(53)),
            ]
        );
        let mut stored_indices = product
            .network
            .store
            .tensors
            .iter()
            .flat_map(TensorStructure::external_indices_iter)
            .collect::<Vec<_>>();
        stored_indices.sort();
        assert_eq!(
            stored_indices,
            [AbstractIndex::Normal(51), AbstractIndex::Normal(53)]
        );
    }
}
