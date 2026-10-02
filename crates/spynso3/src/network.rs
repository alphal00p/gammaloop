use std::{collections::HashMap, ops::Deref, sync::Arc};

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
    iterators::IteratableTensor,
    network::{
        ContractScalars, ExecutionResult, MinIntermediateCost, Network, Sequential,
        SingleSmallestDegree, Steps,
        library::symbolic::{ExplicitKey, TensorLibrary},
        parsing::{ParseSettings, ShadowedStructure},
        store::{NetworkStore, TensorScalarStoreMapping},
        tags::SPENSO_TAG,
    },
    structure::{
        TensorStructure,
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
    id::{Condition, PatternRestriction},
    prelude::*,
};

use symbolica::api::python::ConvertibleToExpression;
#[cfg(test)]
use symbolica::domains::float::Complex as SymComplex;

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

/// An executable tensor calculation that retains component data.
///
/// Networks combine tensor contractions, sums, products, and elementwise
/// functions. Arithmetic creates new networks. ``execute()`` advances a network
/// in place; ``step()`` and ``to_tensor()`` work on copies. ``expression()``
/// returns the symbolic source, while ``status`` reports remaining work.
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
/// >>> result = network.to_tensor()
/// >>> result[0, 1]
/// 2.0
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    from_py_object,
    name = "TensorNetwork",
    module = "symbolica.community.tensor"
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

/// Choose how network execution schedules contractions.
///
/// ``All`` uses the default contraction schedule. ``Single`` selects one
/// contraction at a time. ``Scalar`` restricts contraction to scalar work.
/// Use ``execute(n_steps=...)`` to bound the number of execution rounds;
/// the mode by itself does not bound the total number of rounds.
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
/// >>> from symbolica.community.tensor import ExecutionMode
/// >>> network.execute(n_steps=1, mode=ExecutionMode.Single)
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass_enum)]
#[pyclass(
    from_py_object,
    name = "ExecutionMode",
    module = "symbolica.community.tensor"
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
    pub(crate) fn project_factors(
        py: Python<'_>,
        mut factors: Vec<Self>,
        kind: spenso::shadowing::FactorProjector,
    ) -> PyResult<Self> {
        let sources = factors
            .iter()
            .map(|factor| factor.structure.clone())
            .collect::<Vec<_>>();
        let structure =
            SymbolicTensor::project_factors(&sources, kind).map_err(composition_error)?;
        // Give spectator ports stable identities across permutations. Aligning
        // anonymous axes by position would exchange the spectator data as well.
        let indices = SymbolicTensor::reserved_dummies(
            factors
                .iter()
                .flat_map(|factor| [&factor.structure, &factor.materialized]),
        );
        for factor in &mut factors {
            let channel = factor.structure.matrix_channel().ok_or_else(|| {
                composition_error(composition::TensorCompositionError::NoMatrixChannel)
            })?;
            let replacements = factor
                .structure
                .structure()
                .logical_slots()
                .iter()
                .enumerate()
                .filter(|(position, slot)| {
                    *position != channel.input
                        && *position != channel.output
                        && matches!(slot.aind, PartialIndex::Open(_))
                })
                .map(|(position, _)| (position, indices.fresh_index()))
                .collect();
            *factor = factor.set_port_indices(&replacements)?;
        }
        let terms = kind
            .permutations(factors.len())
            .ok_or_else(|| PyValueError::new_err("factor projector normalization is too large"))?;
        let endpoints = HashMap::from([(0, indices.fresh_index()), (1, indices.fresh_index())]);
        let mut result: Option<Self> = None;
        for (coefficient, order) in terms {
            let first = factors[order[0]].clone();
            let channel = first.structure.matrix_channel().ok_or_else(|| {
                composition_error(composition::TensorCompositionError::NoMatrixChannel)
            })?;
            let mut term = first.chain_form(channel)?;
            for index in &order[1..] {
                let factor = factors[*index].clone();
                let channel = factor.structure.matrix_channel().ok_or_else(|| {
                    composition_error(composition::TensorCompositionError::NoMatrixChannel)
                })?;
                term = term.compose(
                    ConvertibleToSpensoNet(factor),
                    (0, 1),
                    (channel.input, channel.output),
                )?;
            }
            term = term.set_port_indices(&endpoints)?;
            let scalar = TensorExpression::from_structured(
                py,
                SymbolicTensor::new(coefficient, PartialStructure::from_logical_slots([])),
            )?;
            term = term.multiply_network(Self::from_arithmetic(
                ArithmeticStructure::Tensor(scalar),
                None,
            )?)?;
            result = Some(match result {
                Some(sum) => {
                    // Permuting matrix factors must not permute the public
                    // spectator axes. Their explicit identities determine the
                    // alignment before the component networks are added.
                    let slots = term.structure.structure().logical_slots();
                    let axes = sum
                        .structure
                        .structure()
                        .logical_slots()
                        .iter()
                        .map(|slot| {
                            slots.iter().position(|other| other == slot).ok_or_else(|| {
                                PyValueError::new_err(
                                    "projected factors have incompatible spectator ports",
                                )
                            })
                        })
                        .collect::<PyResult<Vec<_>>>()?;
                    sum.add_network(term.permute_axes(axes)?, false)?
                }
                None => term,
            });
        }
        let mut result = result.expect("a nonempty factor group has at least one permutation");
        // Computation is a sum of networks, while both semantic descriptions
        // retain a single projected channel that can be composed again.
        let replacements = result
            .materialized
            .structure()
            .logical_slots()
            .iter()
            .enumerate()
            .map(|(position, slot)| {
                let PartialIndex::Explicit(index) = slot.aind else {
                    unreachable!("network-facing ports are explicit")
                };
                (position, index)
            })
            .collect();
        result.materialized = structure
            .reindex_interface_ports(&replacements)
            .map_err(composition_error)?;
        result.structure = structure;
        result.descriptor = None;
        result.validate_graph_interface()?;
        Ok(result)
    }

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
        } else if let Ok(group) = ob.extract::<crate::projectors::PyFactorProjector>() {
            group.to_network(ob.py()).map(ConvertibleToSpensoNet)
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
        ArithmeticStructure::type_output()
            | SpensoNet::type_output()
            | Spensor::type_output()
            | crate::projectors::PyFactorProjector::type_output()
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
        intern: Option<crate::simplification::Intern>,
    ) -> PyResult<Self> {
        let replacements =
            TensorExpression::index_replacements(self.structure.structure(), indices, intern)?;
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
#[spenso_macros::track_usage(crate::record_usage)]
#[pymethods]
impl SpensoNet {
    /// Create a network from tensor components or symbolic algebra.
    ///
    /// Parameters
    /// ----------
    /// expr : Tensor, TensorNetwork, TensorExpression, or scalar expression
    ///     Computation to represent. A Tensor retains its data; a TensorNetwork
    ///     is copied with its execution progress.
    /// library : TensorLibrary, optional
    ///     Component definitions. Defaults to the built-in four-dimensional Dirac
    ///     and SU(3) library. Unregistered tensors receive symbolic components.
    ///
    /// Returns
    /// -------
    /// TensorNetwork
    ///     A new network with the corresponding external tensor interface.
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
    #[new]
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

    /// Construct the scalar 1 network.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> network = sp.TensorNetwork.one()
    /// >>> value = network.to_tensor().scalar()
    ///
    /// Returns
    /// -------
    /// TensorNetwork
    ///     A rank-zero multiplicative identity.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import TensorNetwork
    /// >>> TensorNetwork.one().result_scalar() == 1
    /// True
    #[staticmethod]
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

    /// Return the symbolic function used to group a tensor subexpression.
    ///
    /// Returns
    /// -------
    /// Expression
    ///     A callable Symbolica head. Its argument is kept as a grouped tensor
    ///     subexpression when parsed as a network.
    ///
    /// Notes
    /// -----
    /// This is useful when constructing raw Symbolica tensor syntax. Tensor-aware
    /// operations and simplifiers already understand this grouping.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica import S
    /// >>> from symbolica.community.tensor import TensorNetwork
    /// >>> grouped = TensorNetwork.bracket()(S("x") + 1)
    #[staticmethod]
    pub fn bracket() -> PythonExpression {
        PythonExpression {
            expr: Atom::var(SPENSO_TAG.bracket),
        }
    }

    /// Register a raw Symbolica function for elementwise tensor application.
    ///
    /// Parameters
    /// ----------
    /// str : str
    ///     Function name to register with the broadcast tag.
    ///
    /// Returns
    /// -------
    /// Expression
    ///     A callable Symbolica function head.
    ///
    /// Notes
    /// -----
    /// Prefer BroadcastFunction when applying the function to TensorExpression,
    /// Tensor, or TensorNetwork objects so the return type follows the operand.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica import S
    /// >>> from symbolica.community.tensor import TensorNetwork
    /// >>> applied = TensorNetwork.broadcast("raw_function")(S("x"))
    #[staticmethod]
    pub fn broadcast(str: &str) -> PythonExpression {
        PythonExpression {
            expr: Atom::var(symbol!(str, tag = SPENSO_TAG.broadcast)),
        }
    }
    /// Construct the scalar 0 network.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> network = sp.TensorNetwork.zero()
    /// >>> value = network.to_tensor().scalar()
    ///
    /// Returns
    /// -------
    /// TensorNetwork
    ///     A rank-zero additive zero.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import TensorNetwork
    /// >>> TensorNetwork.zero().result_scalar() == 0
    /// True
    #[staticmethod]
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

    /// Replace scalar expressions and symbolic component values in a copy.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> from symbolica import E, S
    /// >>> x = S("docs::x")
    /// >>> network = sp.TensorNetwork(sp.TensorExpression(x + 1))
    /// >>> replaced = network.replace(x, E("2"))
    ///
    /// Parameters
    /// ----------
    /// pattern : scalar expression
    ///     Symbolica pattern to find in scalar algebra and component expressions.
    /// rhs : scalar expression or HeldExpression
    ///     Replacement expression. Python replacement callbacks are not supported
    ///     by this network method; use TensorRule with TensorExpression for
    ///     whole-tensor replacement.
    /// cond : PatternRestriction or Condition, optional
    ///     Restriction applied to candidate matches.
    /// non_greedy_wildcards : sequence of Expression, optional
    ///     Wildcards that should prefer shorter matches.
    /// level_range : (int, int or None), optional
    ///     Minimum and maximum matching level; defaults to (0, None).
    /// level_is_tree_depth : bool, optional
    ///     Count tree depth rather than function nesting. Defaults to False.
    /// allow_new_wildcards_on_rhs : bool, optional
    ///     Permit unbound wildcard symbols in the replacement. Defaults to False.
    /// rhs_cache_size : int, optional
    ///     Maximum cached right-hand-side substitutions. Defaults to 100.
    /// repeat : bool, optional
    ///     Repeat replacements until unchanged. Defaults to False.
    ///
    /// Returns
    /// -------
    /// TensorNetwork
    ///     A copy with updated scalar/component expressions. Its external axes
    ///     and symbolic source descriptor are unchanged.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica import S
    /// >>> from symbolica.community.tensor import Tensor, TensorName, Representation
    /// >>> x = S("x")
    /// >>> values = Tensor.dense(TensorName.vector("eval_v")(Representation.euc(2)), [x, x**2])
    /// >>> evaluator = values.evaluator({}, {}, [x], iterations=1, n_cores=1)
    /// >>> from symbolica.community.tensor import TensorNetwork
    /// >>> network = TensorNetwork(values).replace(x, 2)
    /// >>> network.to_tensor()[0] == 2
    /// True
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

    /// Numerically evaluate symbolic component values in a network copy.
    ///
    /// Parameters
    /// ----------
    /// constants : mapping of Expression to float
    ///     Real values for symbolic parameters.
    /// functions : mapping of Expression to callable
    ///     Symbolica function heads and real-valued Python implementations.
    ///     A callback receives one list of evaluated argument values and
    ///     returns a float; for example, lambda args: args[0] ** 2.
    ///
    /// Returns
    /// -------
    /// TensorNetwork
    ///     Network with evaluated component values; contractions are not executed.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica import S
    /// >>> from symbolica.community.tensor import Tensor, TensorName, Representation
    /// >>> x = S("x")
    /// >>> values = Tensor.dense(TensorName.vector("eval_v")(Representation.euc(2)), [x, x**2])
    /// >>> evaluator = values.evaluator({}, {}, [x], iterations=1, n_cores=1)
    /// >>> from symbolica.community.tensor import TensorNetwork
    /// >>> network = TensorNetwork(values).evaluate({x: 2.0}, {})
    /// >>> network.to_tensor()[1]
    /// 4.0
    /// >>> f = S("f_network")
    /// >>> data = Tensor.dense(values.expression(), [f(x), x])
    /// >>> network = TensorNetwork(data).evaluate({x: 2.0}, {f: lambda args: args[0] + 1})
    /// >>> network.to_tensor()[0]
    /// 3.0
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

    /// Advance this network's computation in place.
    ///
    /// Parameters
    /// ----------
    /// library : TensorLibrary, optional
    ///     Component definitions. Defaults to the built-in four-dimensional Dirac
    ///     and SU(3) library. Unregistered tensors receive symbolic components.
    /// function_library : TensorFunctionLibrary, optional
    ///     Numerical implementations of broadcast functions. Defaults to the
    ///     built-in function library.
    /// n_steps : int, optional
    ///     Number of execution rounds. None runs the selected strategy to completion.
    /// mode : ExecutionMode, default ExecutionMode.All
    ///     All uses MinIntermediateCost, which estimates sparse overlap and
    ///     symbolic intermediate cost when choosing each contraction. Single
    ///     selects one minimum-degree contraction at a time. Scalar restricts
    ///     contraction to scalar work.
    ///     Preprocessing and ready operations may also run in an execution round.
    ///
    /// Returns
    /// -------
    /// None
    ///     This network's execution progress and stored intermediate values change.
    ///
    /// Notes
    /// -----
    /// Use ``to_tensor()`` to evaluate a copy, or ``step()`` for an intermediate
    /// copy. ``result_tensor()`` extracts the result after execution.
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
    /// >>> network.execute()
    /// >>> network.result_tensor()[1, 0]
    /// 3.0
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
                            .execute::<Steps<1>, MinIntermediateCost, _, _, _>(lib, fn_lib)
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
                        .execute::<Sequential, MinIntermediateCost, _, _, _>(lib, fn_lib)
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
    /// Read the component tensor from a completed network.
    ///
    /// Parameters
    /// ----------
    /// library : TensorLibrary, optional
    ///     Component definitions. Defaults to the built-in four-dimensional Dirac
    ///     and SU(3) library. Unregistered tensors receive symbolic components.
    ///
    /// Returns
    /// -------
    /// Tensor
    ///     Result components with the source computation's logical axis order.
    ///
    /// Raises
    /// ------
    /// RuntimeError
    ///     Work remains that prevents extracting one result tensor.
    ///
    /// Notes
    /// -----
    /// This does not run pending operations. Call execute() first, or use
    /// to_tensor() to execute a copy and extract the result in one operation.
    /// Existing component storage is preserved: exact Symbolica expressions
    /// are not converted to floating-point values.
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
    /// >>> network.execute()
    /// >>> network.result_tensor().shape
    /// (2, 2)
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
                ExecutionResult::One => Spensor::scalar_with_descriptor(
                    ConcreteOrParam::Param(Atom::one()),
                    descriptor,
                    name,
                    args,
                ),
                ExecutionResult::Zero => Spensor::scalar_with_descriptor(
                    ConcreteOrParam::Param(Atom::Zero),
                    descriptor,
                    name,
                    args,
                ),
                ExecutionResult::Val(value) => {
                    let value = value.into_owned();
                    // Network tensors are already stored in canonical Spenso axis
                    // order. The descriptor permutations only record how that
                    // storage maps back to the public logical interface.
                    Spensor::from_storage_with_descriptor(value, descriptor, name, args)
                }
            },
        )
    }

    /// Read the scalar value of a completed rank-zero network.
    ///
    /// Returns
    /// -------
    /// Expression
    ///     Scalar result as Symbolica algebra.
    ///
    /// Raises
    /// ------
    /// RuntimeError
    ///     The network has not reduced to a scalar result.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import TensorNetwork
    /// >>> network = TensorNetwork(3)
    /// >>> network.execute()
    /// >>> network.result_scalar() == 3
    /// True
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

    /// Return a readable object description for inspection.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> r = sp.Representation.euc(2)
    /// >>> A = sp.TensorName("docs::A")(r, r)
    /// >>> tensor = sp.Tensor.dense(A, [1.0, 2.0, 3.0, 4.0])
    /// >>> network = tensor("i", "j") * tensor("j", "k")
    /// >>> text = repr(network)
    fn __repr__(&self) -> String {
        format!(
            "TensorNetwork({})",
            display::format_structured(&self.structure, false)
        )
    }

    /// Return a readable text representation.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> r = sp.Representation.euc(2)
    /// >>> A = sp.TensorName("docs::A")(r, r)
    /// >>> tensor = sp.Tensor.dense(A, [1.0, 2.0, 3.0, 4.0])
    /// >>> network = tensor("i", "j") * tensor("j", "k")
    /// >>> text = str(network)
    fn __str__(&self) -> String {
        display::format_structured(&self.structure, false)
    }

    fn _repr_latex_(&self) -> String {
        display::structured_to_latex(&self.structure, &Default::default(), None)
    }

    fn _repr_pretty_(&self, pretty: &Bound<'_, PyAny>, cycle: bool) -> PyResult<()> {
        pretty.call_method1("text", (if cycle { "...".into() } else { self.__str__() },))?;
        Ok(())
    }

    /// Export the current network graph in Graphviz DOT syntax.
    ///
    /// Returns
    /// -------
    /// str
    ///     Graph description containing the current tensor nodes and connections.
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
    /// >>> source = network.to_dot()
    fn to_dot(&self) -> String {
        self.network.dot_pretty()
    }

    /// Generate Typst/Linnest source for the network graph.
    ///
    /// Parameters
    /// ----------
    /// config : dict or linnet.RenderConfig, optional
    ///     Graph layout and rendering options. Omit for the standard network view.
    ///
    /// Returns
    /// -------
    /// str
    ///     Typst graph source produced by Linnet.
    ///
    /// Notes
    /// -----
    /// This depicts the current graph, including execution progress. For the
    /// symbolic computation in tensor notation, use formatted() or to_svg().
    /// Graph rendering requires the Linnet/Typst rendering support.
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
    /// >>> output = network.to_linnest()
    #[pyo3(signature = (*, config=None))]
    fn to_linnest(
        &self,
        py: Python<'_>,
        #[gen_stub(override_type(type_repr="builtins.dict[builtins.str, typing.Any] | linnet.RenderConfig | None", imports=("builtins", "typing", "linnet")))]
        config: Option<&Bound<'_, PyAny>>,
    ) -> PyResult<String> {
        self.prepare_render(py, config)?
            .getattr("typst_source")?
            .extract()
    }

    /// Render the current network graph to SVG.
    ///
    /// Parameters
    /// ----------
    /// config : dict or linnet.RenderConfig, optional
    ///     Graph layout and rendering options. Omit for the standard network view.
    ///
    /// Returns
    /// -------
    /// str
    ///     SVG graph, with notebook-theme styling.
    ///
    /// Notes
    /// -----
    /// This depicts the current graph, including execution progress. For the
    /// symbolic computation in tensor notation, use formatted() or to_svg().
    /// Graph rendering requires the Linnet/Typst rendering support.
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
    /// >>> output = network.render()
    #[pyo3(signature = (*, config=None))]
    fn render(
        &self,
        py: Python<'_>,
        #[gen_stub(override_type(type_repr="builtins.dict[builtins.str, typing.Any] | linnet.RenderConfig | None", imports=("builtins", "typing", "linnet")))]
        config: Option<&Bound<'_, PyAny>>,
    ) -> PyResult<String> {
        let svg: String = self
            .prepare_render(py, config)?
            .call_method0("to_svg")?
            .extract()?;
        Ok(display::network::svg_theme(&svg))
    }

    /// Produce compact plain-text tensor notation.
    ///
    /// Parameters
    /// ----------
    /// show_dimensions : bool, optional
    ///     Override settings.show_dimensions for this call. None retains the
    ///     selected settings, whose default omits dimensions.
    /// settings : DisplaySettings, optional
    ///     Tensor notation and display choices. Defaults to DisplaySettings().
    ///
    /// Returns
    /// -------
    /// str
    ///     Readable tensor text, suitable for logs or terminal output.
    ///
    /// Notes
    /// -----
    /// Only ports layout and default spacing are supported in this source
    /// format. Use to_html() or to_svg() for other layouts and custom gaps.
    ///
    /// This formats the symbolic source expression, independent of execution
    /// progress. Use render() or to_html() for the network graph.
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
    /// >>> output = network.format_tensor()
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

    /// Produce static Typst math source.
    ///
    /// Parameters
    /// ----------
    /// show_dimensions : bool, optional
    ///     Override settings.show_dimensions for this call. None retains the
    ///     selected settings, whose default omits dimensions.
    /// settings : DisplaySettings, optional
    ///     Tensor notation and display choices. Defaults to DisplaySettings().
    ///
    /// Returns
    /// -------
    /// str
    ///     Typst source without an enclosing document; no compilation is performed.
    ///
    /// Notes
    /// -----
    /// Only ports layout and default spacing are supported in this source
    /// format. Use to_html() or to_svg() for other layouts and custom gaps.
    ///
    /// This formats the symbolic source expression, independent of execution
    /// progress. Use render() or to_html() for the network graph.
    ///
    /// Static Typst output remains available independently of the HTML explorer.
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
    /// >>> output = network.to_typst()
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

    /// Create a lazy rich display value for a notebook.
    ///
    /// Parameters
    /// ----------
    /// show_dimensions : bool, optional
    ///     Override settings.show_dimensions for this call. None retains the
    ///     selected settings, whose default omits dimensions.
    /// settings : DisplaySettings, optional
    ///     Tensor notation and display choices. Defaults to DisplaySettings().
    /// notation_source : str, optional
    ///     Custom Typst notation source used by the rich renderer. This is
    ///     Typst code, not a filename; omit it for the supplied tensor notation.
    ///
    /// Returns
    /// -------
    /// FormattedOutput
    ///     Each requested backend renders the complete source expression once
    ///     and caches it.
    ///
    /// Notes
    /// -----
    /// This formats the symbolic source expression, independent of execution
    /// progress, retaining a snapshot of that source. Use render() or to_html()
    /// for the network graph.
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
    /// >>> output = network.formatted()
    #[pyo3(signature = (show_dimensions = None, *, settings = None, notation_source = None))]
    fn formatted(
        &self,
        show_dimensions: Option<bool>,
        settings: Option<PyRef<'_, display::DisplaySettings>>,
        notation_source: Option<String>,
    ) -> PythonFormattedOutput {
        let settings = display::resolved_settings(show_dimensions, settings.as_deref());
        display::format_structured_output_rich(
            Arc::new(self.structure.clone()),
            settings,
            notation_source,
        )
    }

    /// Display the current network graph and execution status.
    ///
    /// Parameters
    /// ----------
    /// config : dict or linnet.RenderConfig, optional
    ///     Graph layout and rendering options. Omit for the standard network view.
    ///
    /// Returns
    /// -------
    /// str
    ///     Interactive HTML graph with its progress summary.
    ///
    /// Notes
    /// -----
    /// This depicts the current graph, including execution progress. For the
    /// symbolic computation in tensor notation, use formatted() or to_svg().
    /// Graph rendering requires the Linnet/Typst rendering support.
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
    /// >>> output = network.to_html()
    #[pyo3(signature = (*, config=None))]
    fn to_html(
        &self,
        py: Python<'_>,
        #[gen_stub(override_type(type_repr="builtins.dict[builtins.str, typing.Any] | linnet.RenderConfig | None", imports=("builtins", "typing", "linnet")))]
        config: Option<&Bound<'_, PyAny>>,
    ) -> PyResult<String> {
        Ok(display::network::html(
            &self.render(py, config)?,
            &self.status(),
        ))
    }

    /// Render static mathematical tensor notation to SVG.
    ///
    /// Parameters
    /// ----------
    /// show_dimensions : bool, optional
    ///     Override settings.show_dimensions for this call. None retains the
    ///     selected settings, whose default omits dimensions.
    /// settings : DisplaySettings, optional
    ///     Tensor notation and display choices. Defaults to DisplaySettings().
    /// notation_source : str, optional
    ///     Custom Typst notation source used by the rich renderer. This is
    ///     Typst code, not a filename; omit it for the supplied tensor notation.
    ///
    /// Returns
    /// -------
    /// str
    ///     SVG markup suitable for embedding or saving to a file.
    ///
    /// Notes
    /// -----
    /// Mathematical rendering uses the optional Typst runtime. The returned
    /// string is not automatically displayed; pass it to the notebook's HTML
    /// or SVG display facility.
    ///
    /// This formats the symbolic source expression, independent of execution
    /// progress. Use render() or to_html() for the network graph.
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
    /// >>> output = network.to_svg()
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

    /// Return the symbolic source computation of this network.
    ///
    /// Returns
    /// -------
    /// TensorExpression
    ///     A symbolic tensor with the same external axes.
    ///
    /// Notes
    /// -----
    /// Execution progress and replacement of component values do not rewrite
    /// this source expression. Use result_tensor after execution to inspect values.
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
    /// >>> network.expression().rank
    /// 2
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

    #[doc = python_doc!("Tensor.axes")]
    #[getter]
    #[gen_stub(override_return_type(type_repr = "tuple[Slot | Representation, ...]"))]
    fn axes(&self, py: Python<'_>) -> PyResult<Py<PyTuple>> {
        crate::metadata::axis_objects(py, self.structure.structure().logical_slots())
    }

    /// Canonical external signature of the whole tensor.
    ///
    /// Returns
    /// -------
    /// TensorStructure
    ///     Free axes in canonical order, independent of names and component layout.
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
    /// >>> network.structure.rank
    /// 2
    #[getter]
    fn structure(&self) -> crate::metadata::SpensoTensorStructure {
        crate::metadata::SpensoTensorStructure::from_interface(self.structure.structure())
    }

    /// Assign labels to the unresolved external axes.
    ///
    /// Parameters
    /// ----------
    /// *indices : int, str, Expression, Slot, or AUTO
    ///     One label per unresolved axis, in logical order. AUTO leaves the
    ///     corresponding axis unchanged. A Slot must have the matching
    ///     representation. Repeated compatible labels cause contraction.
    /// intern : {"indices", "flattened"} or None, optional
    ///     Encode compound index labels before checking the tensor structure.
    ///     "indices" retains reversible payloads; "flattened" creates readable names.
    ///     Both preserve index tags and leave scalar arguments and tensor heads intact.
    ///     None leaves labels unchanged and requires ordinary atomic index labels.
    ///
    /// Returns
    /// -------
    /// TensorNetwork
    ///     An indexed network retaining the component data; the original is unchanged.
    ///
    /// Notes
    /// -----
    /// Existing explicit labels are left in place. Use reindex to replace them.
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
    /// >>> indexed = network.index("i", "j")
    /// >>> indexed.rank
    /// 2
    #[pyo3(signature = (*indices, intern = None))]
    fn index(
        &self,
        indices: &Bound<'_, PyTuple>,
        intern: Option<crate::simplification::Intern>,
    ) -> PyResult<Self> {
        self.index_network(indices, intern)
    }

    /// Assign labels to all external axes.
    ///
    /// Parameters
    /// ----------
    /// *indices : int, str, Expression, Slot, or AUTO
    ///     One label per external axis, in logical order. AUTO leaves the
    ///     corresponding axis unchanged. A Slot must have the matching
    ///     representation. Repeated compatible labels cause contraction.
    /// intern : {"indices", "flattened"} or None, optional
    ///     Encode compound index labels before checking the tensor structure.
    ///     "indices" retains reversible payloads; "flattened" creates readable names.
    ///     Both preserve index tags and leave scalar arguments and tensor heads intact.
    ///     None leaves labels unchanged and requires ordinary atomic index labels.
    ///
    /// Returns
    /// -------
    /// TensorNetwork
    ///     An indexed network retaining the component data; the original is unchanged.
    ///
    /// Notes
    /// -----
    /// Existing labels are replaced as well as unresolved axes.
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
    /// >>> indexed = network.reindex("i", "j")
    /// >>> indexed.rank
    /// 2
    #[pyo3(signature = (*indices, intern=None))]
    pub(crate) fn reindex(
        &self,
        indices: &Bound<'_, PyTuple>,
        intern: Option<crate::simplification::Intern>,
    ) -> PyResult<Self> {
        let positions = (0..self.structure.rank()).collect::<Vec<_>>();
        let replacements = TensorExpression::port_replacements(
            self.structure.structure(),
            &positions,
            indices,
            intern,
        )?;
        self.set_port_indices(&replacements)
    }

    /// Rename selected external labels without changing the tensor rank.
    ///
    /// Parameters
    /// ----------
    /// mapping : dict
    ///     Old label to new label. A Slot key also restricts the representation.
    ///     Keys must identify existing external axes. Each axis may be renamed
    ///     once; renaming must preserve the external interface.
    /// intern : {"indices", "flattened"} or None, optional
    ///     Encode compound index labels before checking the tensor structure.
    ///     "indices" retains reversible payloads; "flattened" creates readable names.
    ///     Both preserve index tags and leave scalar arguments and tensor heads intact.
    ///     None leaves labels unchanged and requires ordinary atomic index labels.
    ///
    /// Returns
    /// -------
    /// TensorNetwork
    ///     A new value with renamed axes; the original is unchanged.
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
    /// >>> indexed = network.index("i", "j")
    /// >>> renamed = indexed.rename_indices({"i": "k"})
    /// >>> renamed.rank
    /// 2
    #[pyo3(signature = (mapping, *, intern=None))]
    pub(crate) fn rename_indices(
        &self,
        mapping: &Bound<'_, PyDict>,
        intern: Option<crate::simplification::Intern>,
    ) -> PyResult<Self> {
        let replacements =
            TensorExpression::named_replacements(self.structure.structure(), mapping, intern)?;
        let result = self.set_port_indices(&replacements)?;
        if result.structure.rank() != self.structure.rank() {
            return Err(PyValueError::new_err(
                "renaming would contract external ports; use reindex() to request a contraction",
            ));
        }
        Ok(result)
    }

    /// Reorder the external axes in logical component order.
    ///
    /// Parameters
    /// ----------
    /// axes : sequence of int
    ///     Each current axis position exactly once, in the desired new order.
    ///     For a matrix, [1, 0] exchanges the two axes.
    ///
    /// Returns
    /// -------
    /// TensorNetwork
    ///     A new tensor view with the reordered interface and correspondingly
    ///     reordered component access.
    ///
    /// Notes
    /// -----
    /// This permutes axes; it does not conjugate component values.
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
    /// >>> transposed = network.permute_axes([1, 0])
    /// >>> transposed.shape
    /// (2, 2)
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

    /// Implement ``copy.copy`` using an independent TensorNetwork copy.
    ///
    /// Returns
    /// -------
    /// TensorNetwork
    ///     Copy of the component data and execution progress.
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
    /// >>> import copy
    /// >>> copied = copy.copy(network)
    fn __copy__(&self) -> Self {
        self.clone()
    }

    /// Assign labels to the unresolved external axes.
    ///
    /// Parameters
    /// ----------
    /// *indices : int, str, Expression, Slot, or AUTO
    ///     One label per unresolved axis, in logical order. AUTO leaves the
    ///     corresponding axis unchanged. A Slot must have the matching
    ///     representation. Repeated compatible labels cause contraction.
    /// intern : {"indices", "flattened"} or None, optional
    ///     Encode compound index labels before checking the tensor structure.
    ///     "indices" retains reversible payloads; "flattened" creates readable names.
    ///     Both preserve index tags and leave scalar arguments and tensor heads intact.
    ///     None leaves labels unchanged and requires ordinary atomic index labels.
    ///
    /// Returns
    /// -------
    /// TensorNetwork
    ///     An indexed network retaining the component data; the original is unchanged.
    ///
    /// Notes
    /// -----
    /// Equivalent to ``index(*indices)``.
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
    /// >>> indexed = network("i", "j")
    /// >>> indexed.rank
    /// 2
    #[pyo3(signature = (*indices, intern = None))]
    fn __call__(
        &self,
        indices: &Bound<'_, PyTuple>,
        intern: Option<crate::simplification::Intern>,
    ) -> PyResult<Self> {
        self.index_network(indices, intern)
    }

    /// Negate every tensor component.
    ///
    /// Returns
    /// -------
    /// TensorNetwork
    ///     Negated symbolic algebra retaining component data lazily.
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
    /// >>> negated = -network
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

    /// Add tensors with compatible external interfaces.
    ///
    /// Parameters
    /// ----------
    /// rhs : scalar expression, TensorExpression, Tensor, or TensorNetwork
    ///     Other operand. Component data are retained in a lazy network.
    ///
    /// Returns
    /// -------
    /// TensorNetwork
    ///     New lazy calculation; operands are unchanged.
    ///
    /// Notes
    /// -----
    /// Scalar zero acts as the additive identity. Other scalars can be added
    /// only to rank-zero tensors.
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
    /// >>> result = network + network
    pub fn __add__(&self, rhs: ConvertibleToSpensoNet) -> PyResult<SpensoNet> {
        self.clone().add_network(rhs.to_net(), false)
    }

    /// Implement reflected addition.
    ///
    /// Parameters
    /// ----------
    /// rhs : scalar expression, TensorExpression, Tensor, or TensorNetwork
    ///     Other operand. Component data are retained in a lazy network.
    ///
    /// Returns
    /// -------
    /// TensorNetwork
    ///     New lazy calculation; operands are unchanged.
    ///
    /// Notes
    /// -----
    /// Scalar zero acts as the additive identity. Other scalars can be added
    /// only to rank-zero tensors.
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
    /// >>> result = network + network
    pub fn __radd__(&self, rhs: ConvertibleToSpensoNet) -> PyResult<SpensoNet> {
        self.__add__(rhs)
    }

    /// Subtract tensors with compatible external interfaces.
    ///
    /// Parameters
    /// ----------
    /// rhs : scalar expression, TensorExpression, Tensor, or TensorNetwork
    ///     Other operand. Component data are retained in a lazy network.
    ///
    /// Returns
    /// -------
    /// TensorNetwork
    ///     New lazy calculation; operands are unchanged.
    ///
    /// Notes
    /// -----
    /// Subtraction requires matching tensor axes; a nonzero scalar cannot
    /// be subtracted from a tensor with external axes.
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
    /// >>> result = network - network
    pub fn __sub__(&self, rhs: ConvertibleToSpensoNet) -> PyResult<SpensoNet> {
        self.clone().add_network(rhs.to_net(), true)
    }

    /// Implement reflected subtraction.
    ///
    /// Parameters
    /// ----------
    /// rhs : scalar expression, TensorExpression, Tensor, or TensorNetwork
    ///     Other operand. Component data are retained in a lazy network.
    ///
    /// Returns
    /// -------
    /// TensorNetwork
    ///     New lazy calculation; operands are unchanged.
    ///
    /// Notes
    /// -----
    /// Subtraction requires matching tensor axes; a nonzero scalar cannot
    /// be subtracted from a tensor with external axes.
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
    /// >>> result = network - network
    pub fn __rsub__(&self, rhs: ConvertibleToSpensoNet) -> PyResult<SpensoNet> {
        rhs.to_net().add_network(self.clone(), true)
    }

    /// Multiply tensors, contracting unambiguous compatible axes.
    ///
    /// Parameters
    /// ----------
    /// rhs : scalar expression, TensorExpression, Tensor, or TensorNetwork
    ///     Other operand. Component data are retained in a lazy network.
    ///
    /// Returns
    /// -------
    /// TensorNetwork
    ///     New lazy calculation; operands are unchanged.
    ///
    /// Notes
    /// -----
    /// Matching explicit labels contract. Compatible unresolved axes are
    /// paired only when the choice is unambiguous. Use outer(), contract_ports(),
    /// or compose() to make the intended pairing explicit.
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
    /// >>> result = network * 2
    pub fn __mul__(&self, rhs: ConvertibleToSpensoNet) -> PyResult<SpensoNet> {
        self.clone().multiply_network(rhs.to_net())
    }

    /// Implement reflected multiplication.
    ///
    /// Parameters
    /// ----------
    /// rhs : scalar expression, TensorExpression, Tensor, or TensorNetwork
    ///     Other operand. Component data are retained in a lazy network.
    ///
    /// Returns
    /// -------
    /// TensorNetwork
    ///     New lazy calculation; operands are unchanged.
    ///
    /// Notes
    /// -----
    /// Matching explicit labels contract. Compatible unresolved axes are
    /// paired only when the choice is unambiguous. Use outer(), contract_ports(),
    /// or compose() to make the intended pairing explicit.
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
    /// >>> result = 2 * network
    pub fn __rmul__(&self, rhs: ConvertibleToSpensoNet) -> PyResult<SpensoNet> {
        rhs.to_net().multiply_network(self.clone())
    }

    /// Divide tensor components by a scalar expression.
    ///
    /// Parameters
    /// ----------
    /// rhs : scalar expression, TensorExpression, Tensor, or TensorNetwork
    ///     Other operand. Component data are retained in a lazy network.
    ///
    /// Returns
    /// -------
    /// TensorNetwork
    ///     New lazy calculation; operands are unchanged.
    ///
    /// Notes
    /// -----
    /// The denominator must be scalar. For scalar divided by tensor, the
    /// tensor must also have rank zero.
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
    /// >>> result = network / 2
    pub fn __truediv__(&self, rhs: ConvertibleToSpensoNet) -> PyResult<SpensoNet> {
        self.clone().multiply_network(rhs.to_net().reciprocal()?)
    }

    /// Implement reflected division.
    ///
    /// Parameters
    /// ----------
    /// lhs : scalar expression, TensorExpression, Tensor, or TensorNetwork
    ///     Other operand. Component data are retained in a lazy network.
    ///
    /// Returns
    /// -------
    /// TensorNetwork
    ///     New lazy calculation; operands are unchanged.
    ///
    /// Notes
    /// -----
    /// The denominator must be scalar. For scalar divided by tensor, the
    /// tensor must also have rank zero.
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
    /// >>> result = 2 / TensorNetwork(3)
    pub fn __rtruediv__(&self, lhs: ConvertibleToSpensoNet) -> PyResult<SpensoNet> {
        if !lhs.0.structure.is_scalar() {
            return Err(PyTypeError::new_err(
                "tensor division requires a scalar numerator",
            ));
        }
        lhs.to_net().multiply_network(self.clone().reciprocal()?)
    }

    /// Form an outer product without implicit contractions between the operands.
    ///
    /// Parameters
    /// ----------
    /// rhs : TensorExpression, Tensor, or TensorNetwork
    ///     Tensor whose axes follow this tensor's axes.
    ///
    /// Returns
    /// -------
    /// TensorNetwork
    ///     A lazy outer product retaining the operands' component data.
    ///
    /// Notes
    /// -----
    /// Use this when compatible unresolved axes should remain independent.
    /// Choose explicit distinct labels if you need to refer to them separately.
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
    /// >>> product = network.outer(network)
    /// >>> product.rank
    /// 4
    pub fn outer(&self, rhs: ConvertibleToSpensoNet) -> PyResult<Self> {
        self.clone().apply_product(rhs.to_net(), ProductPlan::Outer)
    }

    /// Contract one chosen pair of axes between two tensors.
    ///
    /// Parameters
    /// ----------
    /// rhs : TensorExpression, Tensor, or TensorNetwork
    ///     Tensor to contract with this tensor.
    /// left, right : int
    ///     Zero-based axis positions in this tensor and rhs, respectively.
    ///     The representations must be compatible under contraction.
    ///
    /// Returns
    /// -------
    /// TensorNetwork
    ///     A lazy contraction retaining component data.
    ///
    /// Notes
    /// -----
    /// Other axes remain external. Axis numbers refer to ``axes``.
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
    /// >>> contracted = network.contract_ports(network, left=1, right=0)
    /// >>> contracted.rank
    /// 2
    #[pyo3(signature = (rhs, *, left, right))]
    pub fn contract_ports(
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

    /// Multiply two tensors along explicitly selected matrix channels.
    ///
    /// Parameters
    /// ----------
    /// rhs : TensorExpression, Tensor, or TensorNetwork
    ///     Next factor in the ordered matrix product.
    /// left, right : tuple of int and int
    ///     (input_axis, output_axis) in this tensor and rhs, respectively.
    ///     The left output contracts with the right input. Other axes are
    ///     retained as spectator axes.
    ///
    /// Returns
    /// -------
    /// TensorNetwork
    ///     A lazy ordered matrix product retaining component data.
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
    /// >>> product = network.compose(network, left=(0, 1), right=(0, 1))
    /// >>> product.rank
    /// 2
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

    /// Contract two rank-one tensors using their representation pairing.
    ///
    /// Parameters
    /// ----------
    /// rhs : TensorExpression, Tensor, or TensorNetwork
    ///     Rank-one tensor in a compatible space.
    ///
    /// Returns
    /// -------
    /// TensorNetwork
    ///     A lazy scalar contraction, including the metric signs of the space.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import TensorName, Tensor, TensorNetwork, Representation
    /// >>> vector = Tensor.dense(TensorName.vector("dot_v")(Representation.euc(2)), [1.0, 2.0])
    /// >>> TensorNetwork(vector).dot(vector).to_tensor().scalar() == 5
    /// True
    pub fn dot(&self, rhs: ConvertibleToSpensoNet) -> PyResult<Self> {
        if self.structure.rank() != 1 || rhs.0.structure.rank() != 1 {
            return Err(PyValueError::new_err(format!(
                "dot() requires rank-one operands, got ranks {} and {}",
                self.structure.rank(),
                rhs.0.structure.rank()
            )));
        }
        self.contract_ports(rhs, 0, 0)
    }

    /// Close a pair of matrix axes on this tensor.
    ///
    /// Parameters
    /// ----------
    /// channel : tuple of int and int, optional
    ///     (input_axis, output_axis) to contract. If omitted, the matrix channel
    ///     must be uniquely determined by the representations. Specify it when
    ///     more than one pairing is possible.
    ///
    /// Returns
    /// -------
    /// TensorNetwork
    ///     The traced tensor; spectator axes remain external.
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
    /// >>> network.trace(channel=(0, 1)).is_scalar
    /// True
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
            assert!(matches!(&result.tensor, ParamOrConcrete::Param(_)));
            let result = Py::new(py, result)?;
            for (indices, expected) in [((0, 0), -1), ((2, 2), 1), ((0, 1), 0)] {
                let component = result
                    .bind(py)
                    .get_item(indices)?
                    .extract::<PythonExpression>()?;
                assert_eq!(component.expr, Atom::num(expected));
            }
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
            let mut product = matrix.contract_ports(ConvertibleToSpensoNet(vector), 1, 0)?;

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
