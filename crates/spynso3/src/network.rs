use std::{collections::HashMap, ops::Deref};

use pyo3::{
    Borrowed,
    exceptions::{self, PyRuntimeError, PyTypeError, PyValueError},
    prelude::*,
    types::PyTuple,
};

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

use crate::{
    composition::{self, StructuredAtom},
    display,
    expression::TensorExpression,
    library::SpensorFunctionLibrary,
};

use super::{Spensor, fresh_open_owner, library::SpensorLibrary, structure::ArithmeticStructure};

use super::ModuleInit;

#[cfg(feature = "python_stubgen")]
use pyo3_stub_gen::{PyStubType, derive::*};

/// A tensor network representing computational graphs of tensor operations.
///
/// A tensor network is a graph-based representation of tensor computations where
/// nodes represent tensors and operations, edges represent tensor contractions and
/// data flow, and the network can be optimized and executed to compute results.
///
/// Tensor networks are particularly useful for symbolic manipulation of complex
/// tensor expressions, optimization of tensor contraction orders, efficient
/// evaluation of large tensor computations, and physics calculations involving
/// many-body systems.
///
/// A network retains the semantic source expression and its public tensor interface
/// separately from the executable graph and its stored values. Value specialization
/// and graph execution therefore do not rewrite the source expression returned by
/// `structure()` or used by the semantic display methods.
///
/// Examples
/// --------
/// >>> from symbolica.community.spenso import TensorNetwork, Tensor, TensorName, Representation
/// >>> rep = Representation.euc(2)
/// >>> T = TensorName("T")
/// >>> tensor = Tensor.dense(T(rep("mu"), rep("nu")), [1.0, 0.0, 0.0, 1.0])
/// >>> network = TensorNetwork(tensor)
/// >>> network.execute()
/// >>> result = network.result_tensor()
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
    pub(crate) structure: StructuredAtom,
    /// The source expression with all graph-facing ports made explicit.
    pub(crate) materialized: StructuredAtom,
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
/// Single : Execute one contraction at a time, useful for debugging
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

/// A Symbolica pattern restriction accepted by tensor-network replacement.
pub struct ReplacementCondition(Condition<PatternRestriction>);

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

fn non_scalar_zero_network(value: &StructuredAtom) -> Option<ParsingNet> {
    if !value.atom.as_view().is_zero() || value.is_scalar() {
        return None;
    }

    let replacements = value
        .interface
        .logical_slots()
        .into_iter()
        .filter_map(|slot| match slot.aind {
            PartialIndex::Explicit(_) => None,
            PartialIndex::Open(id) => Some((
                id,
                composition::fresh_dummy_index([&value.atom], [&value.interface]),
            )),
        })
        .collect();
    let structure: ShadowedStructure<AbstractIndex> = value
        .interface
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
            let structure = TensorExpression::structured(&expression);
            let materialized = TensorExpression::materialized(&expression)?;
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
                    materialized.atom.as_view(),
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
            let value = composition::normalize_closed_root_chain(StructuredAtom::new(
                atom.clone(),
                PartialStructure::from_logical_slots(
                    network
                        .graph
                        .dangling_indices()
                        .into_iter()
                        .map(|slot| slot.rep().slot(PartialIndex::Explicit(slot.aind()))),
                ),
            ))?;
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
            .interface
            .layout()
            .canonical_to_logical(&canonical);
        let replacements = logical
            .into_iter()
            .enumerate()
            .map(|(position, slot)| (position, slot.aind()))
            .collect::<HashMap<_, _>>();
        let materialized = composition::reindex_interface_ports(&value.descriptor, &replacements)
            .map_err(composition_error)?;
        let mut network = Self {
            network: Network::from_tensor(value.tensor),
            materialized,
            structure: value.descriptor,
            descriptor: value
                .descriptor_name
                .map(|name| (name, value.descriptor_args)),
        };
        for (position, slot) in network.semantic_slots().into_iter().enumerate() {
            if let PartialIndex::Explicit(index) = slot.aind {
                network.relabel_port(position, index)?;
            }
        }
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

#[derive(Clone, Copy)]
enum ProductPlan {
    Scalar,
    Outer,
    Contract(composition::PortPair),
    Compose(composition::MatrixChannel, composition::MatrixChannel),
}

fn composition_error(error: composition::TensorCompositionError) -> PyErr {
    PyValueError::new_err(error.to_string())
}

fn additive_interfaces_match(left: &PartialStructure, right: &PartialStructure) -> bool {
    // Unresolved ports remain positional; only explicit indices identify ports across terms.
    left.logical_slots() == right.logical_slots()
        || (left.open_positions().is_empty()
            && right.open_positions().is_empty()
            && left.canonical() == right.canonical())
}

fn product_plan(
    left: &StructuredAtom,
    right: &StructuredAtom,
) -> Result<(StructuredAtom, ProductPlan), composition::TensorCompositionError> {
    let result = composition::multiply(left, right)?;
    if left.is_scalar() || right.is_scalar() {
        return Ok((result, ProductPlan::Scalar));
    }

    let compatible = composition::compatible_pairs(&left.interface, &right.interface);
    if let (Some(left_channel), Some(right_channel)) = (
        composition::matrix_channel(left),
        composition::matrix_channel(right),
    ) {
        let channel_pair = composition::PortPair {
            left: left_channel.output,
            right: right_channel.input,
        };
        if compatible.contains(&channel_pair) {
            return Ok((result, ProductPlan::Compose(left_channel, right_channel)));
        }
    }

    Ok((
        result,
        match compatible.as_slice() {
            [] => ProductPlan::Outer,
            [pair] => ProductPlan::Contract(*pair),
            _ => unreachable!("composition::multiply rejected ambiguous compatible pairs"),
        },
    ))
}

impl SpensoNet {
    fn semantic_slots(&self) -> Vec<spenso::structure::partial::PartialSlot> {
        self.structure.interface.logical_slots()
    }

    fn relabel_port(&mut self, position: usize, index: AbstractIndex) -> PyResult<()> {
        let slots = self.materialized.interface.logical_slots();
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
        if current == index {
            return Ok(());
        }

        let materialized = composition::reindex_interface_ports(
            &self.materialized,
            &HashMap::from([(position, index)]),
        )
        .map_err(composition_error)?;
        let from = slot.rep().slot(current).to_lib();
        let to = slot.rep().slot(index).to_lib();
        let mut reindexed = Vec::new();
        for (tensor_index, tensor) in self.network.store.tensors.iter().enumerate() {
            let slots = tensor.external_structure();
            if !slots.contains(&from) {
                continue;
            }
            let indices = slots
                .into_iter()
                .map(|slot| if slot == from { index } else { slot.aind() })
                .collect::<Vec<_>>();
            let tensor = tensor
                .clone()
                .reindex_storage(&indices)
                .map_err(|error| PyRuntimeError::new_err(error.to_string()))?
                .apply();
            reindexed.push((tensor_index, tensor));
        }
        if !self.network.graph.relabel_dangling_slot(from, to) {
            return Err(PyRuntimeError::new_err(format!(
                "could not uniquely relabel network port {position}"
            )));
        }
        for (tensor_index, tensor) in reindexed {
            self.network.store.tensors[tensor_index] = tensor;
        }
        self.materialized = materialized;
        self.network.state = self.network.graph.state();
        Ok(())
    }

    fn freshen_open_ports(&mut self) -> PyResult<()> {
        let owner = fresh_open_owner();
        for (axis, position) in self
            .structure
            .interface
            .open_positions()
            .into_iter()
            .enumerate()
        {
            self.relabel_port(position, AbstractIndex::Open { owner, axis })?;
        }
        Ok(())
    }

    fn align_pair(left: &mut Self, right: &mut Self, pair: composition::PortPair) -> PyResult<()> {
        let left_slot = left
            .semantic_slots()
            .get(pair.left)
            .copied()
            .ok_or_else(|| {
                PyValueError::new_err(format!("invalid left port position {}", pair.left))
            })?;
        let right_slot = right
            .semantic_slots()
            .get(pair.right)
            .copied()
            .ok_or_else(|| {
                PyValueError::new_err(format!("invalid right port position {}", pair.right))
            })?;
        let target = match (left_slot.aind, right_slot.aind) {
            (PartialIndex::Explicit(left), PartialIndex::Explicit(right)) if left == right => left,
            (PartialIndex::Explicit(index), PartialIndex::Open(_))
            | (PartialIndex::Open(_), PartialIndex::Explicit(index)) => index,
            (PartialIndex::Open(_), PartialIndex::Open(_)) => composition::fresh_dummy_index(
                [&left.structure.atom, &right.structure.atom],
                [&left.structure.interface, &right.structure.interface],
            ),
            (PartialIndex::Explicit(_), PartialIndex::Explicit(_)) => {
                return Err(PyValueError::new_err(format!(
                    "explicit ports {} and {} do not carry the same abstract index",
                    pair.left, pair.right
                )));
            }
        };
        left.relabel_port(pair.left, target)?;
        right.relabel_port(pair.right, target)
    }

    fn align_self_pair(&mut self, pair: composition::PortPair) -> PyResult<()> {
        let slots = self.semantic_slots();
        let left = slots
            .get(pair.left)
            .copied()
            .ok_or_else(|| PyValueError::new_err(format!("invalid port position {}", pair.left)))?;
        let right = slots.get(pair.right).copied().ok_or_else(|| {
            PyValueError::new_err(format!("invalid port position {}", pair.right))
        })?;
        let target = match (left.aind, right.aind) {
            (PartialIndex::Explicit(left), PartialIndex::Explicit(right)) if left == right => left,
            (PartialIndex::Explicit(index), PartialIndex::Open(_))
            | (PartialIndex::Open(_), PartialIndex::Explicit(index)) => index,
            (PartialIndex::Open(_), PartialIndex::Open(_)) => {
                composition::fresh_dummy_index([&self.structure.atom], [&self.structure.interface])
            }
            (PartialIndex::Explicit(_), PartialIndex::Explicit(_)) => {
                return Err(PyValueError::new_err(format!(
                    "explicit ports {} and {} do not carry the same abstract index",
                    pair.left, pair.right
                )));
            }
        };
        self.relabel_port(pair.left, target)?;
        self.relabel_port(pair.right, target)
    }

    fn validate_graph_interface(&self) -> PyResult<()> {
        let mut expected = self
            .materialized
            .interface
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
        structure: StructuredAtom,
        materialized: StructuredAtom,
    ) -> PyResult<Self> {
        let structure =
            composition::normalize_closed_root_chain(structure).map_err(composition_error)?;
        let materialized =
            composition::normalize_closed_root_chain(materialized).map_err(composition_error)?;
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
        if !additive_interfaces_match(&self.structure.interface, &right.structure.interface) {
            return Err(PyValueError::new_err(if subtract {
                "subtraction requires compatible tensor interfaces"
            } else {
                "addition requires compatible tensor interfaces"
            }));
        }
        self.freshen_open_ports()?;
        right.freshen_open_ports()?;
        for position in self.structure.interface.open_positions() {
            Self::align_pair(
                &mut self,
                &mut right,
                composition::PortPair {
                    left: position,
                    right: position,
                },
            )?;
        }
        let structure = StructuredAtom::new(
            if subtract {
                self.structure.atom.as_ref() - right.structure.atom.as_ref()
            } else {
                self.structure.atom.as_ref() + right.structure.atom.as_ref()
            },
            self.structure.interface.clone(),
        );
        let materialized = StructuredAtom::new(
            if subtract {
                self.materialized.atom.as_ref() - right.materialized.atom.as_ref()
            } else {
                self.materialized.atom.as_ref() + right.materialized.atom.as_ref()
            },
            self.materialized.interface.clone(),
        );
        let network = if subtract {
            self.network - right.network
        } else {
            self.network + right.network
        };
        Self::finish(network, structure, materialized)
    }

    pub(crate) fn multiply_network(mut self, mut right: Self) -> PyResult<Self> {
        let (structure, plan) =
            product_plan(&self.structure, &right.structure).map_err(composition_error)?;
        self.freshen_open_ports()?;
        right.freshen_open_ports()?;
        let materialized = match plan {
            ProductPlan::Scalar => composition::multiply(&self.materialized, &right.materialized),
            ProductPlan::Outer => Ok(composition::outer(&self.materialized, &right.materialized)),
            ProductPlan::Contract(pair) => {
                Self::align_pair(&mut self, &mut right, pair)?;
                composition::contract(&self.materialized, &right.materialized, pair)
            }
            ProductPlan::Compose(left_channel, right_channel) => {
                Self::align_pair(
                    &mut self,
                    &mut right,
                    composition::PortPair {
                        left: left_channel.output,
                        right: right_channel.input,
                    },
                )?;
                composition::compose(
                    &self.materialized,
                    &right.materialized,
                    left_channel,
                    right_channel,
                )
            }
        }
        .map_err(composition_error)?;
        Self::finish(self.network * right.network, structure, materialized)
    }

    fn reciprocal(self) -> PyResult<Self> {
        if !self.structure.is_scalar() {
            return Err(PyValueError::new_err(
                "a non-scalar tensor cannot be used as a denominator",
            ));
        }
        let structure = StructuredAtom::new(
            Atom::num(1) / self.structure.atom.as_ref(),
            self.structure.interface.clone(),
        );
        let materialized = StructuredAtom::new(
            Atom::num(1) / self.materialized.atom.as_ref(),
            self.materialized.interface.clone(),
        );
        Self::finish(self.network.pow(-1), structure, materialized)
    }

    pub(crate) fn broadcast_function(self, function: Symbol) -> PyResult<Self> {
        let structure = StructuredAtom::new(
            FunctionBuilder::new(function)
                .add_arg(self.structure.atom.clone())
                .finish(),
            self.structure.interface.clone(),
        );
        let materialized = StructuredAtom::new(
            FunctionBuilder::new(function)
                .add_arg(self.materialized.atom.clone())
                .finish(),
            self.materialized.interface.clone(),
        );
        Self::finish(self.network.fun(function), structure, materialized)
    }

    pub(crate) fn index_network(
        &self,
        indices: &Bound<'_, PyTuple>,
        cook_indices: bool,
    ) -> PyResult<Self> {
        let replacements =
            TensorExpression::index_replacements(&self.structure.interface, indices, cook_indices)?;
        self.set_port_indices(&replacements)
    }

    pub(crate) fn set_port_indices(
        &self,
        replacements: &HashMap<usize, AbstractIndex>,
    ) -> PyResult<Self> {
        let structure = composition::reindex_interface_ports(&self.structure, replacements)
            .map_err(composition_error)?;
        let mut network = self.clone();
        network.freshen_open_ports()?;
        let owner = fresh_open_owner();
        for (axis, &position) in replacements.keys().enumerate() {
            network.relabel_port(position, AbstractIndex::Open { owner, axis })?;
        }
        for (&position, &index) in replacements {
            network.relabel_port(position, index)?;
        }
        network.network.graph.sew_dangling_slots();
        network.network.state = network.network.graph.state();

        let structure_slots = structure.interface.logical_slots();
        let materialized_slots = network.materialized.interface.logical_slots();
        if structure_slots.len() != materialized_slots.len() {
            return Err(PyRuntimeError::new_err(
                "semantic and materialized network interfaces have different ranks",
            ));
        }

        // The graph decides which assigned ports became contractions. Match its
        // dangling slots back to original logical positions to retain their order.
        let mut dangling = network.network.graph.dangling_indices();
        let mut retained = Vec::new();
        for (position, slot) in materialized_slots.iter().enumerate() {
            let PartialIndex::Explicit(index) = slot.aind else {
                return Err(PyRuntimeError::new_err(
                    "materialized network interface contains an unresolved port",
                ));
            };
            let slot = slot.rep().slot(index).to_lib();
            if let Some(found) = dangling.iter().position(|candidate| *candidate == slot) {
                dangling.swap_remove(found);
                retained.push(position);
            }
        }
        if !dangling.is_empty() {
            return Err(PyRuntimeError::new_err(format!(
                "tensor-network topology has dangling slots absent from its materialized interface: {dangling:?}"
            )));
        }

        let contracted = retained.len() < structure_slots.len();
        network.structure = composition::normalize_closed_root_chain(StructuredAtom::new(
            structure.atom,
            PartialStructure::from_logical_slots(
                retained.iter().map(|&position| structure_slots[position]),
            ),
        ))
        .map_err(composition_error)?;
        network.materialized = composition::normalize_closed_root_chain(StructuredAtom::new(
            network.materialized.atom.clone(),
            PartialStructure::from_logical_slots(
                retained
                    .iter()
                    .map(|&position| materialized_slots[position]),
            ),
        ))
        .map_err(composition_error)?;
        if contracted {
            network.descriptor = None;
        }
        network.validate_graph_interface()?;
        Ok(network)
    }

    pub(crate) fn chain_form(self, channel: composition::MatrixChannel) -> PyResult<Self> {
        fn convert(
            value: &StructuredAtom,
            channel: composition::MatrixChannel,
        ) -> Result<StructuredAtom, composition::TensorCompositionError> {
            let slots = value.interface.logical_slots();
            let input = *slots.get(channel.input).ok_or(
                composition::TensorCompositionError::InvalidPort {
                    position: channel.input,
                    rank: slots.len(),
                },
            )?;
            let output = *slots.get(channel.output).ok_or(
                composition::TensorCompositionError::InvalidPort {
                    position: channel.output,
                    rank: slots.len(),
                },
            )?;
            let atom = if value.atom.as_view().is_zero() {
                Atom::Zero
            } else {
                SPENSO_TAG.chain(
                    composition::port_atom(input),
                    composition::port_atom(output),
                    composition::chain_factors(value, channel)?,
                )
            };
            Ok(StructuredAtom::new(
                atom,
                PartialStructure::from_logical_slots(
                    [input, output].into_iter().chain(
                        slots
                            .into_iter()
                            .enumerate()
                            .filter(|(position, _)| {
                                *position != channel.input && *position != channel.output
                            })
                            .map(|(_, slot)| slot),
                    ),
                ),
            ))
        }

        let structure = convert(&self.structure, channel).map_err(composition_error)?;
        let materialized = convert(&self.materialized, channel).map_err(composition_error)?;
        Self::finish(self.network, structure, materialized)
    }
}

// #[gen_stub_pymethods]

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
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
        let value = StructuredAtom::new(
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
        let value = StructuredAtom::new(
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
                    .level_range(setting_level_range)
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
    /// n_steps : int, optional
    ///     Maximum number of execution steps (None for complete execution)
    /// mode : ExecutionMode, optional
    ///     Execution strategy (ExecutionMode.All, ExecutionMode.Scalar, or ExecutionMode.Single)
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
    fn execute(
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
    fn result_tensor(&self, library: Option<&SpensorLibrary>) -> PyResult<Spensor> {
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

    /// Return a DOT representation of the executable network graph.
    ///
    /// Generates a DOT format representation of the computational graph that can be
    /// visualized using graphviz or similar tools.
    fn __str__(&self) -> PyResult<String> {
        Ok(self.network.dot_pretty())
    }

    /// Format the semantic source structure using compact Spenso notation.
    #[pyo3(signature = (show_dimensions = false))]
    fn format_tensor(&self, show_dimensions: bool) -> String {
        display::format_structured(&self.structure, show_dimensions)
    }

    /// Format the semantic source structure as Typst math source.
    #[pyo3(signature = (show_dimensions = false))]
    fn to_typst(&self, show_dimensions: bool) -> String {
        display::structured_to_typst(&self.structure, show_dimensions)
    }

    /// Build Symbolica's rich display wrapper for the semantic source structure.
    #[pyo3(signature = (show_dimensions = false))]
    fn formatted(&self, show_dimensions: bool) -> PythonFormattedOutput {
        display::format_structured_output(&self.structure, show_dimensions)
    }

    /// Return the semantic source expression and its public tensor interface.
    ///
    /// This expression records the tensor-aware structure used for composition and
    /// provenance. It is not reconstructed from the current execution store, so
    /// `replace()`, `evaluate()`, and `execute()` leave it unchanged. Use
    /// `result_scalar()` or `result_tensor()` to inspect the current computed value.
    fn structure(&self, py: Python<'_>) -> PyResult<Py<TensorExpression>> {
        TensorExpression::from_atom_interface_descriptor(
            py,
            self.structure.atom.clone(),
            self.structure.interface.clone(),
            self.descriptor.as_ref().map(|(name, _)| *name),
            self.descriptor
                .as_ref()
                .map(|(_, args)| args.clone())
                .unwrap_or_default(),
        )
    }

    /// Fill the unresolved external ports with `indices` in interface order.
    ///
    /// Pass `AUTO` to leave a port unresolved. Set `cook_indices=True` to flatten nested
    /// symbolic index payloads before insertion.
    #[pyo3(signature = (*indices, cook_indices = false))]
    fn index(&self, indices: &Bound<'_, PyTuple>, cook_indices: bool) -> PyResult<Self> {
        self.index_network(indices, cook_indices)
    }

    /// Fill the unresolved external ports with `indices` in interface order.
    #[pyo3(signature = (*indices, cook_indices = false))]
    fn __call__(&self, indices: &Bound<'_, PyTuple>, cook_indices: bool) -> PyResult<Self> {
        self.index_network(indices, cook_indices)
    }

    pub fn __neg__(&self) -> PyResult<Self> {
        let structure = StructuredAtom::new(
            Atom::num(-1) * self.structure.atom.as_ref(),
            self.structure.interface.clone(),
        );
        let materialized = StructuredAtom::new(
            Atom::num(-1) * self.materialized.atom.as_ref(),
            self.materialized.interface.clone(),
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
        let mut left = self.clone();
        let mut right = rhs.to_net();
        let left_slots = left.semantic_slots();
        let right_slots = right.semantic_slots();
        let collisions =
            composition::compatible_pairs(&left.structure.interface, &right.structure.interface)
                .into_iter()
                .filter(|pair| {
                    matches!(left_slots[pair.left].aind, PartialIndex::Explicit(_))
                        && matches!(right_slots[pair.right].aind, PartialIndex::Explicit(_))
                })
                .collect::<Vec<_>>();
        if !collisions.is_empty() {
            return Err(PyValueError::new_err(format!(
                "outer product cannot preserve equal explicit indices on compatible port pairs {collisions:?}; use distinct indices"
            )));
        }
        let structure = composition::outer(&left.structure, &right.structure);
        left.freshen_open_ports()?;
        right.freshen_open_ports()?;
        let materialized = composition::outer(&left.materialized, &right.materialized);
        Self::finish(left.network * right.network, structure, materialized)
    }

    /// Contract one selected pair of public interface positions.
    #[pyo3(signature = (rhs, *, left, right))]
    pub fn contract(
        &self,
        rhs: ConvertibleToSpensoNet,
        left: usize,
        right: usize,
    ) -> PyResult<Self> {
        let mut lhs = self.clone();
        let mut rhs = rhs.to_net();
        let pair = composition::PortPair { left, right };
        let structure = composition::contract(&lhs.structure, &rhs.structure, pair)
            .map_err(composition_error)?;
        lhs.freshen_open_ports()?;
        rhs.freshen_open_ports()?;
        Self::align_pair(&mut lhs, &mut rhs, pair)?;
        let materialized = composition::contract(&lhs.materialized, &rhs.materialized, pair)
            .map_err(composition_error)?;
        Self::finish(lhs.network * rhs.network, structure, materialized)
    }

    /// Compose two explicitly selected matrix channels.
    #[pyo3(signature = (rhs, *, left, right))]
    pub fn compose(
        &self,
        rhs: ConvertibleToSpensoNet,
        left: (usize, usize),
        right: (usize, usize),
    ) -> PyResult<Self> {
        let mut lhs = self.clone();
        let mut rhs = rhs.to_net();
        let left_channel = composition::MatrixChannel {
            input: left.0,
            output: left.1,
        };
        let right_channel = composition::MatrixChannel {
            input: right.0,
            output: right.1,
        };
        let structure =
            composition::compose(&lhs.structure, &rhs.structure, left_channel, right_channel)
                .map_err(composition_error)?;
        lhs.freshen_open_ports()?;
        rhs.freshen_open_ports()?;
        Self::align_pair(
            &mut lhs,
            &mut rhs,
            composition::PortPair {
                left: left_channel.output,
                right: right_channel.input,
            },
        )?;
        let materialized = composition::compose(
            &lhs.materialized,
            &rhs.materialized,
            left_channel,
            right_channel,
        )
        .map_err(composition_error)?;
        Self::finish(lhs.network * rhs.network, structure, materialized)
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
            None => composition::matrix_channel(&self.structure)
                .ok_or_else(|| PyValueError::new_err("tensor has no unique matrix channel"))?,
        };
        let structure = composition::trace(&self.structure, selected).map_err(composition_error)?;
        let mut network = self.clone();
        network.freshen_open_ports()?;
        network.align_self_pair(composition::PortPair {
            left: selected.input,
            right: selected.output,
        })?;
        let materialized =
            composition::trace(&network.materialized, selected).map_err(composition_error)?;
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

    fn data_tensor(descriptor: StructuredAtom, name: Symbol) -> Spensor {
        let owner = fresh_open_owner();
        let storage = OrderedStructure::new(
            descriptor
                .interface
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
    ) -> (StructuredAtom, Symbol) {
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
        (StructuredAtom::new(atom, interface), name)
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
                network.structure.interface.logical_slots()[0].aind,
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
            let expression = TensorExpression::from_atom_interface_descriptor(
                py,
                descriptor.atom,
                descriptor.interface,
                Some(name),
                Vec::new(),
            )?;

            let indexed_expression = expression.bind(py).call1(("i", "i"))?;
            let indexed_expression = indexed_expression.extract::<PyRef<'_, TensorExpression>>()?;
            assert!(indexed_expression.interface.canonical().is_scalar());
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
                StructuredAtom::new(
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
                matches!(closed.structure.atom.as_view(), AtomView::Fun(function) if function.get_symbol() == SPENSO_TAG.trace)
            );
            assert!(
                matches!(closed.materialized.atom.as_view(), AtomView::Fun(function) if function.get_symbol() == SPENSO_TAG.trace)
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
            let expression = TensorExpression::from_atom_interface_descriptor(
                py,
                descriptor.atom,
                descriptor.interface,
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
            let matrix = TensorExpression::from_atom_interface_descriptor(
                py,
                matrix_descriptor.atom,
                matrix_descriptor.interface,
                Some(matrix_name),
                Vec::new(),
            )?;
            let vector = TensorExpression::from_atom_interface_descriptor(
                py,
                vector_descriptor.atom,
                vector_descriptor.interface,
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
            let expression = TensorExpression::from_atom_interface_descriptor(
                py,
                descriptor.atom,
                descriptor.interface,
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
            let expression = TensorExpression::from_atom_interface_descriptor(
                py,
                descriptor.atom,
                descriptor.interface,
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
            let expected_interface = left.bind(py).borrow().interface.logical_slots();

            let sum = left
                .bind(py)
                .call_method1("__add__", (right.clone_ref(py),))?;
            let difference = left
                .bind(py)
                .call_method1("__sub__", (right.clone_ref(py),))?;
            for value in [&sum, &difference] {
                let expression = value.extract::<PyRef<'_, TensorExpression>>()?;
                assert_eq!(expression.interface.logical_slots(), expected_interface);
            }

            let assert_result = |mut network: SpensoNet, expected: &[f64]| -> PyResult<()> {
                assert_eq!(
                    network.structure.interface.logical_slots(),
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
        let value = StructuredAtom::new(
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
        let value = StructuredAtom::new(
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
            product.structure.interface.logical_slots()[0].aind,
            PartialIndex::Explicit(AbstractIndex::Normal(41))
        );
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
                .interface
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
                .interface
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
