use std::{collections::HashMap, ops::Deref};

use pyo3::{
    exceptions::{self, PyRuntimeError, PyTypeError, PyValueError},
    prelude::*,
    types::{PyDict, PyTuple},
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
        HasStructure, PermutedStructure, ScalarTensor, TensorStructure,
        abstract_index::AbstractIndex,
        partial::{PartialIndex, PartialStructure, PartialStructureExt},
        permuted::Perm,
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
        ConvertibleToPatternRestriction, ConvertibleToReplaceWith, PythonExpression,
        PythonFormattedOutput,
    },
    atom::FunctionBuilder,
    domains::float::Complex as SymComplex,
    prelude::*,
};

use symbolica::api::python::ConvertibleToExpression;

use crate::{
    composition::{self, StructuredAtom},
    display,
    expression::TensorExpression,
    library::SpensorFunctionLibrary,
};

use super::{Spensor, library::SpensorLibrary, structure::ArithmeticStructure};

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
    /// Exact public tensor semantics, including unresolved ports.
    pub(crate) structure: StructuredAtom,
    /// The fully explicit counterpart represented by the execution graph.
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
        .structure
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
            let value = StructuredAtom::new(
                atom.clone(),
                PartialStructure::from_logical_slots(std::iter::empty()),
            );
            Ok(SpensoNet {
                network: ParsingNet::try_from_view(
                    atom.as_view(),
                    library,
                    &ParseSettings::default(),
                )?,
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
        let canonical = value.tensor.structure.external_structure();
        let representation_sorted = value.tensor.index_permutation.apply_slice_inv(&canonical);
        let logical = value
            .tensor
            .rep_permutation
            .apply_slice_inv(&representation_sorted);
        let replacements = logical
            .into_iter()
            .enumerate()
            .map(|(position, slot)| (position, slot.aind()))
            .collect::<HashMap<_, _>>();
        let materialized = composition::reindex_interface_ports(&value.descriptor, &replacements)
            .map_err(composition_error)?;
        let mut network = Self {
            network: Network::from_tensor(value.tensor.structure),
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
            let library = SpensorLibrary::construct();
            network_from_arithmetic(a, &library.library)
                .map(ConvertibleToSpensoNet)
                .map_err(|a| PyRuntimeError::new_err(a.to_string()))
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

fn exact_interface(left: &PartialStructure, right: &PartialStructure) -> bool {
    left.logical_slots() == right.logical_slots()
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
                .reindex(&indices)
                .map_err(|error| PyRuntimeError::new_err(error.to_string()))?
                .permute_reps_wrapped()
                .permute_inds();
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
        for position in self.structure.interface.open_positions() {
            let index = composition::fresh_dummy_index(
                [&self.structure.atom, &self.materialized.atom],
                [&self.structure.interface, &self.materialized.interface],
            );
            self.relabel_port(position, index)?;
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
        network: ParsingNet,
        structure: StructuredAtom,
        materialized: StructuredAtom,
    ) -> PyResult<Self> {
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
        if !exact_interface(&self.structure.interface, &right.structure.interface) {
            return Err(PyValueError::new_err(if subtract {
                "subtraction requires identical ordered tensor interfaces"
            } else {
                "addition requires identical ordered tensor interfaces"
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
        py: Python<'_>,
        indices: &Bound<'_, PyTuple>,
        cook_indices: bool,
    ) -> PyResult<Self> {
        let expression = TensorExpression::from_atom_interface_descriptor(
            py,
            self.structure.atom.clone(),
            self.structure.interface.clone(),
            self.descriptor.as_ref().map(|(name, _)| *name),
            self.descriptor
                .as_ref()
                .map(|(_, args)| args.clone())
                .unwrap_or_default(),
        )?;
        let kwargs = PyDict::new(py);
        kwargs.set_item("cook_indices", cook_indices)?;
        let indexed = expression
            .bind(py)
            .call_method("index", indices, Some(&kwargs))?;
        let indexed = indexed.cast::<TensorExpression>()?.borrow();
        let structure = TensorExpression::structured(&indexed);

        let mut network = self.clone();
        network.freshen_open_ports()?;
        let old_open = self.structure.interface.open_positions();
        let new_slots = structure.interface.logical_slots();
        for position in old_open {
            if let PartialIndex::Explicit(index) = new_slots[position].aind {
                network.relabel_port(position, index)?;
            }
        }
        network.structure = structure;
        network.validate_graph_interface()?;
        Ok(network)
    }

    pub(crate) fn set_port_indices(
        &self,
        replacements: &HashMap<usize, AbstractIndex>,
    ) -> PyResult<Self> {
        let structure = composition::reindex_interface_ports(&self.structure, replacements)
            .map_err(composition_error)?;
        let mut network = self.clone();
        network.freshen_open_ports()?;
        for (&position, &index) in replacements {
            network.relabel_port(position, index)?;
        }
        network.structure = structure;
        network.descriptor = None;
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
    ///     Optional tensor library for resolving named tensor references
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
    pub fn bracket() -> PythonExpression {
        PythonExpression {
            expr: Atom::var(SPENSO_TAG.bracket),
        }
    }

    #[staticmethod]
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

    /// Replace patterns in the tensor network using symbolic pattern matching.
    ///
    /// Applies pattern-based transformations to the network structure, allowing for
    /// symbolic simplifications, substitutions, and algebraic manipulations.
    ///
    /// Parameters
    /// ----------
    /// pattern : Expression
    ///     The symbolic pattern to match against
    /// rhs : Expression
    ///     The replacement expression or pattern
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
    ///     A new TensorNetwork with the replacements applied
    #[pyo3(signature = (pattern, rhs, _cond = None, non_greedy_wildcards = None, level_range = None, level_is_tree_depth = None, allow_new_wildcards_on_rhs = None, rhs_cache_size = None, repeat = None))]
    #[allow(clippy::too_many_arguments)]
    pub fn replace(
        &self,
        pattern: ConvertibleToExpression,
        rhs: ConvertibleToReplaceWith,
        _cond: Option<ConvertibleToPatternRestriction>,
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

        let cond = None;

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

    /// Evaluate symbolic expressions in the network with numerical values.
    ///
    /// Substitutes symbolic constants and functions with numerical values,
    /// converting symbolic parts of the network to concrete numerical tensors.
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
    ///     A new TensorNetwork with symbolic expressions evaluated
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
                    let tensor = PermutedStructure {
                        structure: value,
                        rep_permutation: descriptor.interface.rep_permutation.clone(),
                        index_permutation: descriptor.interface.index_permutation.clone(),
                    };
                    Spensor::from_storage_with_descriptor(tensor, descriptor, name, args)
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

    /// Return a string representation of the network structure.
    ///
    /// Generates a DOT format representation of the computational graph that can be
    /// visualized using graphviz or similar tools.
    fn __str__(&self) -> PyResult<String> {
        Ok(self.network.dot_pretty())
    }

    /// Format the exact semantic structure using compact Spenso notation.
    #[pyo3(signature = (show_dimensions = false))]
    fn format_tensor(&self, show_dimensions: bool) -> String {
        display::format_structured(&self.structure, show_dimensions)
    }

    /// Format the exact semantic structure as Typst math source.
    #[pyo3(signature = (show_dimensions = false))]
    fn to_typst(&self, show_dimensions: bool) -> String {
        display::structured_to_typst(&self.structure, show_dimensions)
    }

    /// Build Symbolica's rich display wrapper for the semantic structure.
    #[pyo3(signature = (show_dimensions = false))]
    fn formatted(&self, show_dimensions: bool) -> PythonFormattedOutput {
        display::format_structured_output(&self.structure, show_dimensions)
    }

    /// Return the exact structured expression represented by this network.
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

    #[pyo3(signature = (*indices, cook_indices = false))]
    fn index(
        &self,
        py: Python<'_>,
        indices: &Bound<'_, PyTuple>,
        cook_indices: bool,
    ) -> PyResult<Self> {
        self.index_network(py, indices, cook_indices)
    }

    #[pyo3(signature = (*indices, cook_indices = false))]
    fn __call__(
        &self,
        py: Python<'_>,
        indices: &Bound<'_, PyTuple>,
        cook_indices: bool,
    ) -> PyResult<Self> {
        self.index_network(py, indices, cook_indices)
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
    use idenso::representations::initialize;
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
        let storage = OrderedStructure::new(
            descriptor
                .interface
                .logical_slots()
                .into_iter()
                .enumerate()
                .map(|(position, slot)| slot.rep().slot(AbstractIndex::Normal(position)))
                .collect(),
        )
        .map_structure(|structure| ShadowedStructure {
            structure,
            global_name: Some(name),
            additional_args: None,
        })
        .map_structure(|structure| SparseTensor::<f64, _>::empty(structure, 0.0).into());
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
