use std::{collections::HashMap, sync::Arc};

use idenso::{
    IndexTooling,
    color::CS,
    dirac::AGS,
    tensor::{
        SymbolicTensor, composition,
        inference::{InterfaceInference, TensorInferenceError},
    },
};
use pyo3::{
    exceptions::{
        PyIndexError, PyKeyboardInterrupt, PyOverflowError, PyRuntimeError, PyTypeError,
        PyValueError,
    },
    prelude::*,
    types::{PyAny, PyDict, PySlice, PyTuple},
};
#[cfg(feature = "python_stubgen")]
use pyo3_stub_gen::{
    PyStubType, TypeInfo,
    generate::MethodType,
    inventory::submit,
    type_info::{MethodInfo, ParameterDefault, ParameterInfo, ParameterKind, PyMethodsInfo},
};
use spenso::{
    network::{
        library::symbolic::{ETS, ExplicitKey},
        tags::SPENSO_TAG,
    },
    structure::{
        Canonicalized, HasName, TensorStructure,
        abstract_index::AbstractIndex,
        partial::{PartialIndex, PartialStructure, PartialStructureExt},
        representation::{LibraryRep, Representation},
        slot::{IsAbstractSlot, Slot},
    },
};
use symbolica::{
    api::python::{
        ConvertibleToExpression, PythonExpression, PythonExpressionEvaluator, PythonFormattedOutput,
    },
    atom::{Atom, AtomView, Symbol},
    domains::float::Complex,
    evaluate::{EvaluatorBuilder, OptimizationSettings},
};

use crate::{
    ModuleInit, SliceOrIntOrExpanded, Spensor,
    aliases::AliasedTensorExpression,
    display,
    library::SpensorLibrary,
    network::{ConvertibleToSpensoNet, ExecutionMode, SpensoNet},
    simplification::{
        CanonicalizationError, CookingError, DiracAdjointError, NetworkToolingError, PyCookSettings,
    },
    structure::{
        ArithmeticStructure, ConvertibleToAbstractIndex, ConvertibleToDimension,
        ConvertibleToSpensoName, SpensoIndices, SpensoName, SpensoRepresentation, SpensoSlot,
    },
};

/// The local placeholder used to leave a tensor port unresolved.
#[cfg_attr(feature = "python_stubgen", pyo3_stub_gen::derive::gen_stub_pyclass)]
#[pyclass(frozen, name = "_AutoIndex", module = "symbolica.community.tensor")]
pub struct AutoIndex;

#[cfg(feature = "python_stubgen")]
pyo3_stub_gen::module_variable!("symbolica.community.tensor", "AUTO", AutoIndex);
#[cfg(feature = "python_stubgen")]
pyo3_stub_gen::module_variable!("symbolica.community.tensor", "_", AutoIndex);

#[spenso_macros::track_usage(crate::record_usage)]
#[pymethods]
impl AutoIndex {
    fn __repr__(&self) -> &'static str {
        "_"
    }

    fn __str__(&self) -> &'static str {
        "_"
    }
}

/// An immutable symbolic tensor with an ordered external interface.
///
/// `to_expression()` is the explicit boundary for generic Symbolica algebra.
///
/// Predefined tensors are constructed by typed factories and indexed afterward
/// in logical interface order. Dimensions are representation metadata, not
/// scalar tensor arguments.
///
/// Construct a tensor expression, preserving existing tensor metadata or inferring it.
///
/// Pass `cook_indices=sp.CookSettings.reversible()` to reversibly encode nested index payloads before
/// inferring the tensor interface for arbitrary nested index expressions.
/// Custom CookSettings control the index encoding, source filters, and output tags.
/// An explicit `structure` is checked against the resulting expression and
/// retains declared logical order and the interface of a symbolic zero.
///
/// Examples
/// --------
/// >>> from symbolica.community import tensor as sp
/// >>> from symbolica import E, S
/// >>> r = sp.Representation.euc(2)
/// >>> A = sp.TensorName("docs::A")(r, r)
/// >>> payload = S("docs::edge")(7, 2)
/// >>> label = S("docs::index")
/// >>> raw = A(label, "j").to_expression().replace(label, payload)
/// >>> value = sp.TensorExpression(raw, cook_indices=sp.CookSettings.reversible())
/// >>> mink = sp.Representation.mink(4)
/// >>> metric = sp.TensorExpression.g(mink)("mu", "nu")
/// >>> gamma = sp.TensorExpression.gamma(4)("i", "j", "mu")
/// >>> sigma = sp.TensorExpression.sigma(4)("mu", "nu", "i", "j")
///
/// Parameters
/// ----------
/// expression : TensorExpression | Expression or scalar
///     TensorExpression to copy, or scalar/Symbolica expression from which to infer the external tensor ports.
/// structure : TensorStructure, optional
///     Optional TensorStructure declaring the logical port order. It must agree with the expression and can retain the interface of a symbolic zero.
/// cook_indices : CookSettings, optional
///     None leaves index payloads unchanged. CookSettings selects encoding, filters, and output tags; CookSettings.indices() flattens nested index payloads.
#[cfg_attr(feature = "python_stubgen", pyo3_stub_gen::derive::gen_stub_pyclass)]
#[pyclass(
    frozen,
    from_py_object,
    name = "TensorExpression",
    module = "symbolica.community.tensor"
)]
#[derive(Clone)]
pub struct TensorExpression {
    pub(crate) value: Arc<SymbolicTensor<PartialStructure>>,
    /// An optional data identity, independent of the immutable Symbolica atom.
    ///
    /// In particular, naming a composite expression must not rewrite it into a
    /// single symbolic tensor leaf.
    pub(crate) name: Option<Symbol>,
    /// Scalar key arguments belonging to an inferred atomic data identity.
    /// Explicitly naming a composite starts a fresh identity with no key args.
    pub(crate) name_args: Vec<Atom>,
}

impl TensorExpression {
    pub(crate) fn inference_error(error: TensorInferenceError) -> PyErr {
        let message = match error {
            TensorInferenceError::Interrupted => return PyKeyboardInterrupt::new_err(()),
            TensorInferenceError::Network(error) => {
                return NetworkToolingError::new_err(error.to_string());
            }
            TensorInferenceError::InvalidBuiltinSignature { factory, signature } => format!(
                "predefined tensor `{factory}` requires {signature}; use TensorExpression.{factory}(...) to construct it"
            ),
            error => error.to_string(),
        };
        PyValueError::new_err(message)
    }

    /// Construct from untrusted symbolic syntax. Infer missing structure and
    /// validate supplied interfaces here, before tensor operations may reuse them.
    pub fn from_atom_interface(
        py: Python<'_>,
        atom: Atom,
        interface: impl Into<Option<PartialStructure>>,
    ) -> PyResult<Py<Self>> {
        if let Some(interface) = interface.into() {
            let original_rank = interface.canonical().order();
            let interface = InterfaceInference::merge_explicit_interface_sequence(&[interface])
                .map_err(Self::inference_error)?
                .canonicalize_open_ports();
            SymbolicTensor::validate_interface(&atom, &interface).map_err(Self::inference_error)?;
            let (name, arguments) = if interface.canonical().order() == original_rank {
                inferred_descriptor(atom.as_view())
            } else {
                (None, Vec::new())
            };
            return Self::from_known_parts(py, atom, interface, name, arguments);
        }
        let (name, name_args) = inferred_descriptor(atom.as_view());
        let value = if InterfaceInference::has_structured_syntax(atom.as_view()) {
            SymbolicTensor::infer(atom).map_err(Self::inference_error)?
        } else {
            SymbolicTensor::new(atom, PartialStructure::from_logical_slots([]))
        };
        let original_rank = value.rank();
        let value = value.finish_parts().map_err(Self::inference_error)?;
        let (name, name_args) = if value.rank() == original_rank {
            (name, name_args)
        } else {
            (None, Vec::new())
        };
        Self::from_shared(py, Arc::new(value), name, name_args)
    }

    /// Construct from tensor-aware parts whose interface was established by the
    /// caller (composition, indexing, storage, or validated syntax). Normalize
    /// contractions, but do not infer the same interface again.
    pub(crate) fn from_known_parts(
        py: Python<'_>,
        atom: Atom,
        interface: PartialStructure,
        name: Option<Symbol>,
        name_args: Vec<Atom>,
    ) -> PyResult<Py<Self>> {
        let original_rank = interface.canonical().order();
        let value =
            SymbolicTensor::checked_parts(atom, interface).map_err(Self::inference_error)?;
        // A contraction is a derived expression, not a declaration of the original stored data.
        let (name, name_args) = if value.rank() == original_rank {
            (name, name_args)
        } else {
            (None, Vec::new())
        };
        Self::from_shared(py, Arc::new(value), name, name_args)
    }

    /// Construct without inferring, validating, or normalizing the supplied parts.
    /// The caller is responsible for matching the atom, ordered slots, and identity;
    /// subsequent tensor operations trust this invariant.
    pub fn from_parts_unchecked(
        py: Python<'_>,
        atom: Atom,
        interface: PartialStructure,
        name: Option<Symbol>,
        name_args: Vec<Atom>,
    ) -> PyResult<Py<Self>> {
        Self::from_shared(
            py,
            Arc::new(SymbolicTensor::from_normalized_parts(atom, interface)),
            name,
            name_args,
        )
    }

    pub(crate) fn from_shared(
        py: Python<'_>,
        value: Arc<SymbolicTensor<PartialStructure>>,
        name: Option<Symbol>,
        name_args: Vec<Atom>,
    ) -> PyResult<Py<Self>> {
        Py::new(
            py,
            Self {
                value,
                name,
                name_args,
            },
        )
    }

    /// Both symbolic payload forms use Symbolica's evaluator state and options.
    pub(crate) fn build_evaluator(
        py: Python<'_>,
        builder: EvaluatorBuilder<'_>,
        iterations: usize,
        n_cores: usize,
    ) -> PyResult<PythonExpressionEvaluator> {
        let settings = OptimizationSettings::new()
            .horner_iterations(iterations)
            .cores(n_cores);
        let evaluator = py
            .detach(|| builder.optimization_settings(settings).build())
            .map_err(|error| PyValueError::new_err(error.to_string()))?;
        Ok(PythonExpressionEvaluator {
            rational_constants: evaluator.get_constants().to_vec(),
            eval_complex: evaluator.map_coeff(&|c| Complex::new(c.re.to_f64(), c.im.to_f64())),
            eval_real: None,
            #[cfg(feature = "native")]
            jit_real: None,
            #[cfg(feature = "native")]
            jit_complex: None,
            eval_double_float: None,
            eval_double_float_complex: None,
            eval_arb_prec: None,
            eval_arb_prec_complex: None,
            jit_compile: false,
            #[cfg(feature = "native")]
            jit_settings: symbolica::evaluate::JITCompilationSettings::new(),
        })
    }

    pub(crate) fn atom(&self) -> &Atom {
        self.value.expression()
    }

    pub(crate) fn interface(&self) -> &PartialStructure {
        self.value.structure()
    }

    pub(crate) fn from_transformed_atom(
        self_: &PyRef<'_, TensorExpression>,
        py: Python<'_>,
        atom: Atom,
    ) -> PyResult<Py<TensorExpression>> {
        let (value, contracted) = Self::structured(self_)
            .with_transformed_expression(atom)
            .map_err(Self::inference_error)?;
        let (name, args) = if contracted {
            (None, Vec::new())
        } else {
            Self::transformed_descriptor(self_, value.expression())
        };
        Self::from_shared(py, Arc::new(value), name, args)
    }

    fn with_cooked_indices(
        self_: &PyRef<'_, Self>,
        py: Python<'_>,
        settings: &idenso::CookSettings,
    ) -> PyResult<Py<Self>> {
        let value = Self::structured(self_)
            .with_cooked_indices(settings)
            .map_err(|error| CookingError::new_err(error.to_string()))?;
        let (name, args) = if value.rank() == self_.interface().canonical().order() {
            Self::transformed_descriptor(self_, value.expression())
        } else {
            (None, Vec::new())
        };
        Self::from_shared(py, Arc::new(value), name, args)
    }

    fn transformed_descriptor(self_: &PyRef<'_, Self>, atom: &Atom) -> (Option<Symbol>, Vec<Atom>) {
        let original = inferred_descriptor(self_.atom().as_view());
        if original == (self_.name, self_.name_args.clone()) {
            inferred_descriptor(atom.as_view())
        } else {
            // Explicit names assigned to composite expressions are metadata.
            (self_.name, self_.name_args.clone())
        }
    }

    /// Wrap a transformed atom after validating that it retains this tensor's
    /// interface. Preserve logical port order and the interface of tensor zeros.
    pub fn preserving_interface(
        self_: &PyRef<'_, Self>,
        py: Python<'_>,
        atom: Atom,
    ) -> PyResult<Py<Self>> {
        let (name, args) = Self::transformed_descriptor(self_, &atom);
        let value = Self::structured(self_)
            .with_checked_expression(atom)
            .map_err(Self::inference_error)?;
        let (name, args) = if value.rank() == self_.interface().canonical().order() {
            (name, args)
        } else {
            (None, Vec::new())
        };
        Self::from_shared(py, Arc::new(value), name, args)
    }

    pub(crate) fn from_structured(
        py: Python<'_>,
        value: SymbolicTensor<PartialStructure>,
    ) -> PyResult<Py<Self>> {
        Self::from_shared(py, Arc::new(value), None, Vec::new())
    }

    fn promoted_network(self_: &PyRef<'_, Self>, py: Python<'_>) -> PyResult<SpensoNet> {
        SpensoNet::from_arithmetic(
            ArithmeticStructure::Tensor(Py::new(py, Self::clone(self_))?),
            None,
        )
    }

    pub fn structured(&self) -> &SymbolicTensor<PartialStructure> {
        &self.value
    }

    pub(crate) fn index_replacements(
        interface: &PartialStructure,
        indices: &Bound<'_, PyTuple>,
        cook_indices: Option<&PyCookSettings>,
    ) -> PyResult<HashMap<usize, AbstractIndex>> {
        let open_positions = interface.open_positions();
        Self::port_replacements(interface, &open_positions, indices, cook_indices)
    }

    pub(crate) fn port_replacements(
        interface: &PartialStructure,
        positions: &[usize],
        indices: &Bound<'_, PyTuple>,
        cook_indices: Option<&PyCookSettings>,
    ) -> PyResult<HashMap<usize, AbstractIndex>> {
        if indices.len() != positions.len() {
            return Err(PyValueError::new_err(format!(
                "expected {} indices for selected ports, got {}",
                positions.len(),
                indices.len()
            )));
        }

        let logical = interface.logical_slots();
        let mut replacements = HashMap::new();
        for (argument, &position) in indices.iter().zip(positions) {
            if argument.is_instance_of::<AutoIndex>() {
                continue;
            }
            let index = if let Ok(slot) = argument.extract::<SpensoSlot>() {
                if logical[position].rep() != slot.slot.rep() {
                    return Err(PyValueError::new_err(format!(
                        "slot at unresolved port {position} has an incompatible representation"
                    )));
                }
                slot.slot.aind()
            } else {
                index_value(argument.extract()?, cook_indices)?
            };
            replacements.insert(position, index);
        }
        Ok(replacements)
    }

    pub(crate) fn named_replacements(
        interface: &PartialStructure,
        mapping: &Bound<'_, PyDict>,
        cook_indices: Option<&PyCookSettings>,
    ) -> PyResult<HashMap<usize, AbstractIndex>> {
        let slots = interface.logical_slots();
        let mut result = HashMap::new();
        for (key, value) in mapping.iter() {
            let typed = key.extract::<SpensoSlot>().ok();
            let index = if let Some(slot) = &typed {
                slot.slot.aind()
            } else {
                index_value(key.extract()?, cook_indices)?
            };
            let positions = slots
                .iter()
                .enumerate()
                .filter_map(|(position, slot)| {
                    (slot.aind == PartialIndex::Explicit(index)
                        && typed
                            .as_ref()
                            .is_none_or(|typed| typed.slot.rep() == slot.rep()))
                    .then_some(position)
                })
                .collect::<Vec<_>>();
            if positions.is_empty() {
                return Err(pyo3::exceptions::PyKeyError::new_err(format!(
                    "no external index matches {key}"
                )));
            }
            let arguments = PyTuple::new(mapping.py(), positions.iter().map(|_| value.clone()))?;
            for (position, index) in
                Self::port_replacements(interface, &positions, &arguments, cook_indices)?
            {
                if result.insert(position, index).is_some() {
                    return Err(PyValueError::new_err(format!(
                        "external port {position} is renamed more than once"
                    )));
                }
            }
        }
        Ok(result)
    }

    pub(crate) fn descriptor_name(self_: &PyRef<'_, Self>) -> Option<Symbol> {
        self_.name
    }

    pub(crate) fn descriptor_args(self_: &PyRef<'_, Self>) -> Vec<Atom> {
        self_.name_args.clone()
    }

    pub(crate) fn materialized_atom(self_: &PyRef<'_, Self>) -> PyResult<Atom> {
        Ok(Self::structured(self_)
            .materialized()
            .map_err(|error| PyValueError::new_err(error.to_string()))?
            .into_expression())
    }

    pub(crate) fn from_indices(py: Python<'_>, indices: &SpensoIndices) -> PyResult<Py<Self>> {
        let canonical = indices.structure.canonical().external_structure();
        let logical = indices.structure.layout().canonical_to_logical(&canonical);
        let interface = PartialStructure::from_logical_slots(
            logical
                .into_iter()
                .map(|slot| slot.rep().slot(PartialIndex::Explicit(slot.aind()))),
        );
        let atom = SymbolicTensor::from_canonicalized(&indices.structure)
            .ok_or_else(|| PyValueError::new_err("indexed tensor structure has no name"))?
            .into_expression();
        Self::from_known_parts(
            py,
            atom,
            interface,
            indices.structure.canonical().name(),
            indices.structure.canonical().args().unwrap_or_default(),
        )
    }

    pub(crate) fn from_structure(
        py: Python<'_>,
        structure: &Canonicalized<ExplicitKey<AbstractIndex>>,
    ) -> PyResult<Py<Self>> {
        let value = SymbolicTensor::from_signature(structure)
            .map_err(|error| PyRuntimeError::new_err(error.to_string()))?;
        Self::from_shared(
            py,
            Arc::new(value),
            structure.canonical().name(),
            structure.canonical().args().unwrap_or_default(),
        )
    }
}

impl ModuleInit for TensorExpression {}

// Operands are transient stack values; boxing this dispatch variant would add
// an allocation to every structured Python arithmetic operation.
#[allow(clippy::large_enum_variant)]
pub(crate) enum TensorOperand {
    Structured(SymbolicTensor<PartialStructure>),
    Scalar(Atom),
}

#[derive(IntoPyObject)]
enum TensorDispatch {
    Expression(Py<TensorExpression>),
    Network(Py<SpensoNet>),
}

#[cfg(feature = "python_stubgen")]
impl PyStubType for TensorDispatch {
    fn type_output() -> TypeInfo {
        TensorExpression::type_output() | SpensoNet::type_output()
    }
}

fn is_concrete_network(value: &Bound<'_, PyAny>) -> bool {
    value.is_instance_of::<SpensoNet>() || value.is_instance_of::<Spensor>()
}

fn concrete_network(value: &Bound<'_, PyAny>) -> PyResult<Option<SpensoNet>> {
    if let Ok(network) = value.extract::<SpensoNet>() {
        Ok(Some(network))
    } else if let Ok(tensor) = value.extract::<Spensor>() {
        SpensoNet::from_tensor(tensor).map(Some)
    } else {
        Ok(None)
    }
}

impl TensorOperand {
    fn into_structured(self) -> PyResult<SymbolicTensor<PartialStructure>> {
        match self {
            Self::Structured(value) => Ok(value),
            Self::Scalar(atom) => {
                SymbolicTensor::checked_parts(atom, PartialStructure::from_logical_slots([]))
                    .map_err(TensorExpression::inference_error)
            }
        }
    }

    pub(crate) fn extract(value: &Bound<'_, PyAny>) -> PyResult<Self> {
        if let Ok(value) = value.extract::<PyRef<'_, TensorExpression>>() {
            return Ok(Self::Structured(
                TensorExpression::structured(&value).clone(),
            ));
        }
        if let Ok(value) = value.extract::<ConvertibleToExpression>() {
            let atom = value.to_expression().expr;
            return if InterfaceInference::has_structured_syntax(atom.as_view()) {
                SymbolicTensor::infer(atom)
                    .map(Self::Structured)
                    .map_err(TensorExpression::inference_error)
            } else {
                Ok(Self::Scalar(atom))
            };
        }
        Err(PyTypeError::new_err(
            "expected a TensorExpression or Symbolica expression",
        ))
    }
}

pub(crate) fn index_value(
    value: ConvertibleToAbstractIndex,
    cook_indices: Option<&PyCookSettings>,
) -> PyResult<AbstractIndex> {
    match value {
        ConvertibleToAbstractIndex::Aind(index) => Ok(index),
        ConvertibleToAbstractIndex::Atom(expression) => {
            let atom = match (cook_indices, expression.expr.as_view()) {
                (Some(settings), AtomView::Fun(fun)) => {
                    settings.rust().cook_function(fun).map_err(|error| {
                        CookingError::new_err(format!("cannot cook index: {error:?}"))
                    })?
                }
                _ => expression.expr,
            };
            atom.as_view().try_into().map_err(|error| {
                let hint = if cook_indices.is_some() {
                    ""
                } else {
                    " Try passing cook_indices=CookSettings.indices()."
                };
                PyValueError::new_err(format!(
                    "cannot convert `{atom}` to an abstract index: {error}.{hint}"
                ))
            })
        }
        ConvertibleToAbstractIndex::Separator => Err(PyValueError::new_err(
            "the ';' separator is not valid in TensorExpression.index()",
        )),
    }
}

fn logical_dimensions(interface: &PartialStructure) -> PyResult<Vec<usize>> {
    interface
        .logical_slots()
        .into_iter()
        .enumerate()
        .map(|(position, slot)| {
            usize::try_from(slot.dim()).map_err(|_| {
                PyValueError::new_err(format!(
                    "tensor layout requires a concrete dimension at interface position {position}"
                ))
            })
        })
        .collect()
}

fn logical_size(dimensions: &[usize]) -> PyResult<usize> {
    dimensions.iter().try_fold(1usize, |size, dimension| {
        size.checked_mul(*dimension)
            .ok_or_else(|| PyOverflowError::new_err("tensor interface size overflows usize"))
    })
}

fn expanded_index(dimensions: &[usize], flat: usize) -> PyResult<Vec<usize>> {
    let size = logical_size(dimensions)?;
    if flat >= size {
        return Err(PyIndexError::new_err("flat index out of bounds"));
    }
    let mut remainder = flat;
    let mut expanded = vec![0; dimensions.len()];
    for (index, dimension) in dimensions.iter().enumerate().rev() {
        expanded[index] = remainder % dimension;
        remainder /= dimension;
    }
    Ok(expanded)
}

fn flat_index(dimensions: &[usize], expanded: &[usize]) -> PyResult<usize> {
    if expanded.len() != dimensions.len() {
        return Err(PyIndexError::new_err(format!(
            "expected {} coordinates, got {}",
            dimensions.len(),
            expanded.len()
        )));
    }
    expanded
        .iter()
        .zip(dimensions)
        .enumerate()
        .try_fold(0usize, |flat, (position, (index, dimension))| {
            if index >= dimension {
                return Err(PyIndexError::new_err(format!(
                    "coordinate {index} is out of bounds for dimension {dimension} at position {position}"
                )));
            }
            Ok(flat * dimension + index)
        })
}

fn inferred_descriptor(atom: AtomView<'_>) -> (Option<Symbol>, Vec<Atom>) {
    match atom {
        AtomView::Fun(function) if function.get_symbol().has_tag(&SPENSO_TAG.tensor) => {
            let mut inference = InterfaceInference::default();
            let name = SpensoName {
                name: function.get_symbol(),
            };
            // Generic bound leaves are derived contractions. Predefined compact
            // signatures, such as g(P(rep), Q(rep)), keep their data descriptors.
            if name.builtin_factory().is_none()
                && function
                    .iter()
                    .any(|argument| inference.compact_vector_port(argument))
            {
                return (None, Vec::new());
            }
            let args = function
                .iter()
                .filter(|argument| {
                    Slot::<LibraryRep, AbstractIndex>::try_from(*argument).is_err()
                        && Representation::<LibraryRep>::try_from(*argument).is_err()
                })
                .map(|argument| argument.to_owned())
                .collect();
            (Some(function.get_symbol()), args)
        }
        _ => (None, Vec::new()),
    }
}

#[cfg_attr(feature = "python_stubgen", pyo3_stub_gen::derive::gen_stub_pymethods)]
#[cfg_attr(not(feature = "python_stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[spenso_macros::track_usage(crate::record_usage)]
#[pymethods]
impl TensorExpression {
    /// Construct a tensor expression, preserving existing tensor metadata or inferring it.
    ///
    /// Pass `cook_indices=sp.CookSettings.reversible()` to reversibly encode nested index payloads before
    /// inferring the tensor interface for arbitrary nested index expressions.
    /// Custom CookSettings control the index encoding, source filters, and output tags.
    /// An explicit `structure` is checked against the resulting expression and
    /// retains declared logical order and the interface of a symbolic zero.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> from symbolica import E, S
    /// >>> r = sp.Representation.euc(2)
    /// >>> A = sp.TensorName("docs::A")(r, r)
    /// >>> payload = S("docs::edge")(7, 2)
    /// >>> label = S("docs::index")
    /// >>> raw = A(label, "j").to_expression().replace(label, payload)
    /// >>> value = sp.TensorExpression(raw, cook_indices=sp.CookSettings.reversible())
    ///
    /// Parameters
    /// ----------
    /// expression : TensorExpression | Expression or scalar
    ///     TensorExpression to copy, or scalar/Symbolica expression from which to infer the external tensor ports.
    /// structure : TensorStructure, optional
    ///     Optional TensorStructure declaring the logical port order. It must agree with the expression and can retain the interface of a symbolic zero.
    /// cook_indices : CookSettings, optional
    ///     None leaves index payloads unchanged. CookSettings selects encoding, filters, and output tags; CookSettings.indices() flattens nested index payloads.
    #[new]
    #[pyo3(signature = (expression, *, structure = None, cook_indices = None))]
    #[gen_stub(override_return_type(type_repr = "TensorExpression"))]
    fn new(
        py: Python<'_>,
        expression: &Bound<'_, PyAny>,
        structure: Option<&crate::metadata::SpensoTensorStructure>,
        cook_indices: Option<&PyCookSettings>,
    ) -> PyResult<Self> {
        if let Some(structure) = structure {
            let atom = if let Ok(value) = expression.extract::<PyRef<'_, Self>>() {
                value.atom().clone()
            } else {
                expression
                    .extract::<ConvertibleToExpression>()?
                    .to_expression()
                    .expr
            };
            let atom = if let Some(settings) = cook_indices {
                settings
                    .rust()
                    .try_cook_indices(atom.as_view())
                    .map_err(|error| {
                        CookingError::new_err(format!("cannot cook indices: {error:?}"))
                    })?
            } else {
                atom
            };
            let checked = Self::from_atom_interface(py, atom, Some(structure.interface.clone()))?;
            return Ok(checked.borrow(py).clone());
        }
        let tensor = if let Some(settings) = cook_indices {
            if let Ok(expression) = expression.extract::<PyRef<'_, Self>>() {
                Self::with_cooked_indices(&expression, py, &settings.rust())?
            } else {
                let atom = expression
                    .extract::<ConvertibleToExpression>()?
                    .to_expression()
                    .expr;
                let atom = settings
                    .rust()
                    .try_cook_indices(atom.as_view())
                    .map_err(|error| {
                        CookingError::new_err(format!("cannot cook indices: {error:?}"))
                    })?;
                Self::from_atom_interface(py, atom, None)?
            }
        } else {
            as_tensor(py, expression)?
        };
        let tensor = tensor.borrow(py);
        Ok(tensor.clone())
    }

    /// Create an unresolved metric with ports in `rep` and `other`.
    ///
    /// `other` defaults to `rep`; it may also be the dual of the same space.
    /// For example, `g(fund, fund.dual())("i", "j")` creates a fundamental identity.
    /// Call the result with two indices to fill its ports in logical order.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> r = sp.Representation.euc(2)
    /// >>> tensor = sp.TensorExpression.g(r)
    /// >>> ports = tensor.structure.slots
    ///
    /// Parameters
    /// ----------
    /// rep : Representation
    ///     Representation, including its dimension, for the tensor’s ports.
    /// other : Representation, optional
    ///     Representation paired with rep. None uses rep itself; a dualizable space can instead be paired with rep.dual().
    #[staticmethod]
    #[pyo3(signature = (rep, other = None))]
    fn g(
        py: Python<'_>,
        rep: &SpensoRepresentation,
        other: Option<&SpensoRepresentation>,
    ) -> PyResult<Py<TensorExpression>> {
        let left = rep.representation;
        let right = other.unwrap_or(rep).representation;
        if left != right && left.dual() != right {
            return Err(PyValueError::new_err(
                "metric ports must belong to the same representation space (possibly dual)",
            ));
        }
        let structure = ExplicitKey::<AbstractIndex>::from_iter([left, right], ETS.metric, None);
        Self::from_structure(py, &structure)
    }

    /// Create an unresolved musical-isomorphism tensor for `rep`.
    ///
    /// Call the result with two indices to fill its ports in logical order.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> r = sp.Representation.euc(2)
    /// >>> tensor = sp.TensorExpression.flat(r)
    /// >>> ports = tensor.structure.slots
    ///
    /// Parameters
    /// ----------
    /// rep : Representation
    ///     Representation, including its dimension, for the tensor’s ports.
    #[staticmethod]
    fn flat(py: Python<'_>, rep: &SpensoRepresentation) -> PyResult<Py<TensorExpression>> {
        let structure = ExplicitKey::<AbstractIndex>::from_iter(
            [rep.representation, rep.representation],
            ETS.flat,
            None,
        );
        Self::from_structure(py, &structure)
    }

    /// Construct a gamma matrix with two four-dimensional bispinor ports followed by
    /// one Lorentz port. Call the result as `gamma(i, j, mu)` to assign its indices.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> tensor = sp.TensorExpression.gamma(4)
    /// >>> ports = tensor.structure.slots
    ///
    /// Parameters
    /// ----------
    /// minkowski_dimension : int | Expression | str
    ///     Dimension of the Lorentz representation. The bispinor ports remain four-dimensional; a symbolic Lorentz dimension is allowed.
    #[staticmethod]
    fn gamma(
        py: Python<'_>,
        minkowski_dimension: ConvertibleToDimension,
    ) -> PyResult<Py<TensorExpression>> {
        Self::from_structure(py, &AGS.gamma_strct::<AbstractIndex>(minkowski_dimension.0))
    }

    /// Construct gamma-five with two bispinor ports. Call the result as `gamma5(i, j)`.
    /// Four-dimensional gamma-five identities require four-dimensional spinor and Lorentz spaces.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> tensor = sp.TensorExpression.gamma5(4)
    /// >>> ports = tensor.structure.slots
    ///
    /// Parameters
    /// ----------
    /// spinor_dimension : int | Expression | str
    ///     Dimension of the bispinor representation. Use 4 for the ordinary four-component Dirac matrices.
    #[staticmethod]
    fn gamma5(
        py: Python<'_>,
        spinor_dimension: ConvertibleToDimension,
    ) -> PyResult<Py<TensorExpression>> {
        Self::from_structure(py, &AGS.gamma5_strct::<AbstractIndex>(spinor_dimension.0))
    }

    /// Construct the left-chiral projector with two bispinor ports.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> tensor = sp.TensorExpression.projm(4)
    /// >>> ports = tensor.structure.slots
    ///
    /// Parameters
    /// ----------
    /// spinor_dimension : int | Expression | str
    ///     Dimension of the bispinor representation. Use 4 for the ordinary four-component Dirac matrices.
    #[staticmethod]
    fn projm(
        py: Python<'_>,
        spinor_dimension: ConvertibleToDimension,
    ) -> PyResult<Py<TensorExpression>> {
        Self::from_structure(py, &AGS.projm_strct::<AbstractIndex>(spinor_dimension.0))
    }

    /// Construct the right-chiral projector with two bispinor ports.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> tensor = sp.TensorExpression.projp(4)
    /// >>> ports = tensor.structure.slots
    ///
    /// Parameters
    /// ----------
    /// spinor_dimension : int | Expression | str
    ///     Dimension of the bispinor representation. Use 4 for the ordinary four-component Dirac matrices.
    #[staticmethod]
    fn projp(
        py: Python<'_>,
        spinor_dimension: ConvertibleToDimension,
    ) -> PyResult<Py<TensorExpression>> {
        Self::from_structure(py, &AGS.projp_strct::<AbstractIndex>(spinor_dimension.0))
    }

    /// Construct the sigma tensor with two Lorentz and two four-dimensional bispinor ports.
    /// Inspect `.structure.slots` before indexing when composing it with other tensor conventions.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> tensor = sp.TensorExpression.sigma(4)
    /// >>> ports = tensor.structure.slots
    ///
    /// Parameters
    /// ----------
    /// minkowski_dimension : int | Expression | str
    ///     Dimension of the Lorentz representation. The bispinor ports remain four-dimensional; a symbolic Lorentz dimension is allowed.
    #[staticmethod]
    fn sigma(
        py: Python<'_>,
        minkowski_dimension: ConvertibleToDimension,
    ) -> PyResult<Py<TensorExpression>> {
        Self::from_structure(py, &AGS.sigma_strct::<AbstractIndex>(minkowski_dimension.0))
    }

    /// Create an unresolved color structure constant.
    ///
    /// All three logical ports are adjoint representations of
    /// `adjoint_dimension`; call the result with `(a, b, c)`.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> tensor = sp.TensorExpression.f(8)
    /// >>> ports = tensor.structure.slots
    ///
    /// Parameters
    /// ----------
    /// adjoint_dimension : int | Expression | str
    ///     Dimension of the adjoint color representation, such as 8 for SU(3). Symbolic dimensions are also accepted.
    #[staticmethod]
    fn f(
        py: Python<'_>,
        adjoint_dimension: ConvertibleToDimension,
    ) -> PyResult<Py<TensorExpression>> {
        Self::from_structure(py, &CS.f_strct::<AbstractIndex>(adjoint_dimension.0))
    }

    /// Create an unresolved color generator.
    ///
    /// The logical ports are adjoint, fundamental, and antifundamental. Dimensions
    /// belong to those representations and are not tensor key arguments.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> tensor = sp.TensorExpression.t(8, 3)
    /// >>> ports = tensor.structure.slots
    ///
    /// Parameters
    /// ----------
    /// adjoint_dimension : int | Expression | str
    ///     Dimension of the adjoint color representation, such as 8 for SU(3). Symbolic dimensions are also accepted.
    /// fundamental_dimension : int | Expression | str
    ///     Dimension of the fundamental color representation, such as 3 for SU(3).
    #[staticmethod]
    fn t(
        py: Python<'_>,
        adjoint_dimension: ConvertibleToDimension,
        fundamental_dimension: ConvertibleToDimension,
    ) -> PyResult<Py<TensorExpression>> {
        Self::from_structure(
            py,
            &CS.t_strct::<AbstractIndex>(fundamental_dimension.0, adjoint_dimension.0),
        )
    }

    /// Number of external tensor ports in the ordered interface.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> r = sp.Representation.euc(2)
    /// >>> A = sp.TensorName("docs::A")(r, r)
    /// >>> result = A.rank
    #[getter]
    fn rank(&self) -> usize {
        self.interface().canonical().order()
    }

    /// Whether the expression has no external tensor ports.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> r = sp.Representation.euc(2)
    /// >>> A = sp.TensorName("docs::A")(r, r)
    /// >>> result = A.is_scalar
    #[getter]
    fn is_scalar(&self) -> bool {
        self.interface().canonical().is_scalar()
    }

    /// Tensor identity, scalar arguments, and external ports in logical order.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> r = sp.Representation.euc(2)
    /// >>> A = sp.TensorName("docs::A")(r, r)
    /// >>> result = A.shape
    #[getter]
    #[gen_stub(override_return_type(type_repr = "tuple[int | Expression, ...]"))]
    fn shape(&self, py: Python<'_>) -> PyResult<Py<PyTuple>> {
        self.structure().shape(py)
    }

    /// Tensor identity, scalar arguments, and external ports in logical order.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> r = sp.Representation.euc(2)
    /// >>> A = sp.TensorName("docs::A")(r, r)
    /// >>> result = A.structure
    #[getter]
    fn structure(&self) -> crate::metadata::SpensoTensorStructure {
        crate::metadata::SpensoTensorStructure {
            interface: self.interface().clone(),
            name: self.name,
            arguments: self.name_args.clone(),
        }
    }

    /// The optional identity used when this expression describes stored data.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> r = sp.Representation.euc(2)
    /// >>> A = sp.TensorName("docs::A")(r, r)
    /// >>> result = A.name
    #[getter]
    fn name(&self) -> Option<SpensoName> {
        self.name.map(|name| SpensoName { name })
    }

    /// Return a named descriptor without changing the underlying expression.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> r = sp.Representation.euc(2)
    /// >>> A = sp.TensorName("docs::A")(r, r)
    /// >>> named = (2 * A).with_name("docs::twice_A")
    ///
    /// Parameters
    /// ----------
    /// name : TensorName | str | Expression | TensorExpression
    ///     New data identity, as a TensorName, name string, symbol, or named TensorExpression. The tensor’s mathematical expression and logical interface are retained.
    fn with_name(
        self_: PyRef<'_, Self>,
        py: Python<'_>,
        name: ConvertibleToSpensoName,
    ) -> PyResult<Py<Self>> {
        let ConvertibleToSpensoName(name, args) = name;
        Self::from_shared(py, Arc::clone(&self_.value), Some(name.name), args)
    }

    /// Return the number of component positions in the tensor’s logical shape.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> r = sp.Representation.euc(2)
    /// >>> A = sp.TensorName("docs::A")(r, r)
    /// >>> size = len(A)
    fn __len__(&self) -> PyResult<usize> {
        logical_size(&logical_dimensions(self.interface())?)
    }

    /// Convert between flat positions and tensor coordinates in logical row-major order.
    /// An integer returns coordinates, a coordinate list returns a flat position, and a slice
    /// returns coordinate lists. This indexes the layout; use `.components()` for values.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> r = sp.Representation.euc(2)
    /// >>> A = sp.TensorName("docs::A")(r, r)
    /// >>> coordinates = A[3]
    /// >>> flat_index = A[[1, 1]]
    /// >>> all_coordinates = A[:]
    ///
    /// Parameters
    /// ----------
    /// item : int
    ///     Flat integer position, list of coordinates, or slice of flat positions. Coordinates follow .structure.slots and use row-major order.
    #[gen_stub(skip)]
    fn __getitem__(&self, py: Python<'_>, item: SliceOrIntOrExpanded) -> PyResult<Py<PyAny>> {
        let dimensions = logical_dimensions(self.interface())?;
        match item {
            SliceOrIntOrExpanded::Selection(_) => Err(PyTypeError::new_err(
                "coordinate conversion accepts nonnegative integers, integer coordinates, or a flat slice",
            )),
            SliceOrIntOrExpanded::Int(index) => expanded_index(&dimensions, index)?
                .into_pyobject(py)
                .map(|value| value.unbind().into_any()),
            SliceOrIntOrExpanded::Expanded(indices) => Ok(flat_index(&dimensions, &indices)?
                .into_pyobject(py)
                .map(|value| value.unbind().into_any())?),
            SliceOrIntOrExpanded::Slice(slice) => {
                let size = logical_size(&dimensions)?;
                let indices = slice.indices(size as isize)?;
                let expanded = (0..indices.slicelength)
                    .map(|offset| {
                        let index = indices.start + offset as isize * indices.step;
                        expanded_index(&dimensions, index as usize)
                    })
                    .collect::<PyResult<Vec<_>>>()?;
                expanded
                    .into_pyobject(py)
                    .map(|value| value.unbind().into_any())
            }
        }
    }

    /// Drop the structured interface and return an ordinary Symbolica `Expression`.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> r = sp.Representation.euc(2)
    /// >>> A = sp.TensorName("docs::A")(r, r)
    /// >>> scalar_algebra = A.to_expression()
    /// >>> restored = sp.TensorExpression(scalar_algebra)
    fn to_expression(self_: PyRef<'_, Self>) -> PythonExpression {
        PythonExpression {
            expr: self_.atom().clone(),
        }
    }

    /// Apply one tensor replacement rule or an ordered list of rules.
    /// The matched tensor must have explicit ports and a supported tensor structure.
    /// The replacement retains its external interface; index the tensor before creating
    /// a rule that targets unresolved ports.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> r = sp.Representation.euc(2)
    /// >>> A = sp.TensorName("docs::A")(r, r)
    /// >>> A = A("i", "j")
    /// >>> rule = sp.TensorRule(A, 2 * A)
    /// >>> replaced = A.replace(rule)
    ///
    /// Parameters
    /// ----------
    /// rules : TensorRule | list[TensorRule]
    ///     One TensorRule or an ordered list of rules. The first matching rule wins; the replacement is checked against the matched tensor interface.
    #[pyo3(signature = (rules))]
    fn replace(
        self_: PyRef<'_, Self>,
        py: Python<'_>,
        rules: &Bound<'_, PyAny>,
    ) -> PyResult<Py<Self>> {
        let rules = crate::tensor_rule::PyTensorRule::extract_rules(rules)?;
        let value = self_
            .value
            .replace_rules(&rules)
            .map_err(Self::inference_error)?;
        let (name, args) = Self::transformed_descriptor(&self_, value.expression());
        Self::from_shared(py, Arc::new(value), name, args)
    }

    /// Expand scalar algebra while preserving and validating the tensor interface.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> from symbolica import E, S
    /// >>> r = sp.Representation.euc(2)
    /// >>> A = sp.TensorName("docs::A")(r, r)
    /// >>> A = A("i", "j")
    /// >>> x = S("docs::x")
    /// >>> expanded = ((x + 1)**2 * A).expand()
    ///
    /// Parameters
    /// ----------
    /// var : Expression or scalar, optional
    ///     Optional scalar variable with respect to which to expand. None expands all supported scalar algebra.
    /// via_poly : bool, optional
    ///     Whether to use Symbolica’s polynomial expansion route. None uses direct expansion, as with False.
    #[pyo3(signature = (var = None, via_poly = None))]
    fn expand(
        self_: PyRef<'_, Self>,
        py: Python<'_>,
        var: Option<ConvertibleToExpression>,
        via_poly: Option<bool>,
    ) -> PyResult<Py<Self>> {
        let var = var.map(|value| value.to_expression());
        let value = Self::structured(&self_)
            .expanded(
                var.as_ref().map(|value| value.expr.as_view()),
                via_poly.unwrap_or(false),
            )
            .map_err(Self::inference_error)?;
        if value.expression() == self_.atom() {
            return Py::new(py, self_.clone());
        }
        let (name, args) = Self::transformed_descriptor(&self_, value.expression());
        Self::from_shared(py, Arc::new(value), name, args)
    }

    /// Build a numerical evaluator for repeated evaluations of a scalar tensor expression.
    /// Supply parameter values in the order given by `params`; each input row is one evaluation.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> from symbolica import E, S
    /// >>> x = S("docs::x")
    /// >>> evaluator = sp.TensorExpression((x + 1)**2).evaluator([x])
    /// >>> values = evaluator.evaluate([[2.0], [3.0]])
    ///
    /// Parameters
    /// ----------
    /// params : Sequence[Expression]
    ///     Symbolica variables in the exact order expected in each numerical input row. Include every free numerical parameter that remains after substitutions.
    /// iterations : int
    ///     Number of optimization iterations when constructing the evaluator. More iterations can improve repeated evaluations at a higher setup cost.
    /// n_cores : int
    ///     Worker count used while optimizing the evaluator. Use 1 for serial setup.
    #[pyo3(signature = (params, *, iterations=1, n_cores=1))]
    fn evaluator(
        &self,
        py: Python<'_>,
        params: Vec<PythonExpression>,
        iterations: usize,
        n_cores: usize,
    ) -> PyResult<PythonExpressionEvaluator> {
        let params = params
            .into_iter()
            .map(|value| value.expr)
            .collect::<Vec<_>>();
        Self::build_evaluator(py, self.value.evaluator(&params), iterations, n_cores)
    }

    /// Copy the value and its tensor metadata without mutating the original.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> r = sp.Representation.euc(2)
    /// >>> A = sp.TensorName("docs::A")(r, r)
    /// >>> import copy
    /// >>> independent = copy.copy(A)
    fn __copy__(self_: PyRef<'_, Self>, py: Python<'_>) -> PyResult<Py<Self>> {
        Py::new(py, self_.clone())
    }

    /// Replace four-dimensional Lorentz slots and compact representations by dimension `D`.
    ///
    /// Spinor and color dimensions, scalar coefficients, and other Lorentz dimensions
    /// are unchanged. Apply this before contracting four-dimensional Lorentz indices.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> from symbolica import E, S
    /// >>> D = S("docs::D")
    /// >>> gamma = sp.TensorExpression.gamma(4).with_lorentz_dimension(D)
    ///
    /// Parameters
    /// ----------
    /// dimension : int | Expression | str
    ///     New Lorentz dimension, as an integer, string, or atomic Symbolica expression. This changes Lorentz slots; spinor and color dimensions are retained.
    fn with_lorentz_dimension(
        self_: PyRef<'_, Self>,
        py: Python<'_>,
        dimension: ConvertibleToDimension,
    ) -> PyResult<Py<Self>> {
        let value = Self::structured(&self_)
            .with_lorentz_dimension(dimension.0)
            .map_err(Self::inference_error)?;
        let (name, args) = if value.rank() == self_.interface().canonical().order() {
            (self_.name, self_.name_args.clone())
        } else {
            (None, Vec::new())
        };
        Self::from_shared(py, Arc::new(value), name, args)
    }

    /// Apply the selected tensor identities while keeping intermediate expressions factored.
    /// Choose `AlgebraSettings.hep()` for gamma, color, and epsilon identities.
    /// The result retains shared subexpressions; call `.expand()` when an explicit result is needed.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> gamma = sp.TensorExpression.gamma(4)
    /// >>> trace = gamma("a", "b", "mu") * gamma("b", "a", "nu")
    /// >>> reduced = trace.reduce_algebra(sp.AlgebraSettings.hep())
    /// >>> expanded = reduced.expand()
    ///
    /// Parameters
    /// ----------
    /// settings : AlgebraSettings, optional
    ///     AlgebraSettings selecting identity families. None uses the default settings with no optional gamma, color, or epsilon family enabled.
    #[pyo3(signature=(settings=None))]
    fn reduce_algebra(
        &self,
        settings: Option<&crate::simplification::PyAlgebraSettings>,
    ) -> PyResult<AliasedTensorExpression> {
        let settings = settings.copied().unwrap_or_default();
        let value = self
            .value
            .reduce_algebra(settings.rust())
            .map_err(Self::inference_error)?;
        Ok(AliasedTensorExpression::from_parts(
            value,
            self.name,
            self.name_args.clone(),
        ))
    }

    /// Place explicit indices in a named scope so independently constructed expressions
    /// can be combined without accidental dummy-index collisions. The tensor’s logical
    /// port order is retained. Use `dummies_only=True` to keep free labels shared.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> from symbolica import E, S
    /// >>> r = sp.Representation.euc(2)
    /// >>> A = sp.TensorName("docs::A")(r, r)
    /// >>> scoped = A("i", "j").wrap_indices(S("docs::diagram1"))
    ///
    /// Parameters
    /// ----------
    /// header : Expression
    ///     Symbolica symbol naming the index scope. Use different symbols for independently built expressions whose dummy labels must not collide.
    /// dummies_only : bool
    ///     If True, scope only contracted indices, leaving the external labels unchanged. False scopes all explicit indices.
    #[pyo3(signature = (header, *, dummies_only=false))]
    fn wrap_indices(
        self_: PyRef<'_, Self>,
        py: Python<'_>,
        header: Symbol,
        dummies_only: bool,
    ) -> PyResult<Py<Self>> {
        let value = self_
            .value
            .wrap_indices(header, dummies_only)
            .map_err(Self::inference_error)?;
        Self::from_shared(py, Arc::new(value), self_.name, self_.name_args.clone())
    }

    /// Return the ordered external, uncontracted indices as Symbolica expressions.
    ///
    /// Raises `ValueError` when the expression cannot be parsed as a tensor network.
    ///
    /// Uses the stored logical slot order when all external ports are indexed.
    ///
    /// Identifies and returns all indices that are not summed over (i.e., not dummy
    /// indices). These are the "free" indices that appear in the final result and
    /// determine the tensor rank of the expression. For dualizable representations,
    /// downstairs indices are represented wrapped in `dind(...)`.
    ///
    /// This is essential for:
    /// - Verifying index conservation in tensor equations
    /// - Determining the rank and structure of tensor expressions
    /// - Debugging index contractions
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> r = sp.Representation.euc(2)
    /// >>> A = sp.TensorName("docs::A")(r, r)
    /// >>> free_indices = A("i", "j").list_dangling()
    ///
    /// Returns
    /// -------
    /// A list of expressions, each representing a free (dangling) index.
    ///
    /// Raises
    /// ------
    /// `ValueError` when the expression cannot be parsed as a tensor network.
    fn list_dangling(self_: PyRef<'_, Self>) -> PyResult<Vec<PythonExpression>> {
        if self_.interface().open_positions().is_empty() {
            return Ok(self_
                .interface()
                .logical_slots()
                .into_iter()
                .map(|slot| composition::port_atom(slot).into())
                .collect());
        }
        self_
            .atom()
            .list_dangling::<AbstractIndex>()
            .map(|indices| indices.into_iter().map(Into::into).collect())
            .map_err(|error| PyValueError::new_err(error.to_string()))
    }

    /// Canonically order tensor factors and deterministically rename contracted indices.
    ///
    /// Raises `CanonicalizationError` when the expression cannot be parsed or canonicalized.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> r = sp.Representation.euc(2)
    /// >>> A = sp.TensorName("docs::A")(r, r)
    /// >>> canonical = (A("i", "j") * A("j", "k")).canonize()
    fn canonize(self_: PyRef<'_, Self>, py: Python<'_>) -> PyResult<Py<Self>> {
        let value = Self::structured(&self_)
            .canonize()
            .map_err(|error| match error {
                idenso::CanonicalizationError::Tensor(error) => Self::inference_error(error),
                error => CanonicalizationError::new_err(error.to_string()),
            })?;
        let (name, args) = Self::transformed_descriptor(&self_, value.expression());
        Self::from_shared(py, Arc::new(value), name, args)
    }

    /// Take the complex-conjugate, reversed Dirac chain, including the gamma-zero
    /// factors needed at open bispinor endpoints. Use `preserve_indices=True` when
    /// external labels must remain attached to the original amplitude legs.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> gamma = sp.TensorExpression.gamma(4)("a", "b", "mu")
    /// >>> adjoint = gamma.dirac_adjoint(preserve_indices=True)
    ///
    /// Parameters
    /// ----------
    /// preserve_indices : bool
    ///     Keep external labels attached to their original amplitude legs. False uses the usual endpoint exchange when taking the Dirac adjoint.
    ///
    /// Returns
    /// -------
    /// The representation-aware Dirac adjoint.
    #[pyo3(signature = (*, preserve_indices=false))]
    fn dirac_adjoint(
        self_: PyRef<'_, Self>,
        py: Python<'_>,
        preserve_indices: bool,
    ) -> PyResult<Py<Self>> {
        let result = self_
            .atom()
            .dirac_adjoint::<AbstractIndex>(preserve_indices)
            .map_err(|error| DiracAdjointError::new_err(error.to_string()))?;
        Self::from_transformed_atom(&self_, py, result)
    }

    /// Normalize index-free scalar-product notation without contracting indexed factors.
    /// Use `contract()` first when explicit vector indices should be contracted.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> r = sp.Representation.euc(2)
    /// >>> p = sp.TensorName.vector("docs::p")(r)
    /// >>> q = sp.TensorName.vector("docs::q")(r)
    /// >>> value = sp.dot(p, q)
    /// >>> rewritten = value.to_dots()
    fn to_dots(self_: PyRef<'_, Self>, py: Python<'_>) -> PyResult<Py<Self>> {
        let value = self_.value.to_dots().map_err(Self::inference_error)?;
        let (name, args) = Self::transformed_descriptor(&self_, value.expression());
        Self::from_shared(py, Arc::new(value), name, args)
    }

    /// Open dots into symbolic indexed contractions, retaining other shorthands.
    /// This does not expand the numerator or evaluate finite tensor components.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> r = sp.Representation.euc(2)
    /// >>> p = sp.TensorName.vector("docs::p")(r)
    /// >>> q = sp.TensorName.vector("docs::q")(r)
    /// >>> value = sp.dot(p, q)
    /// >>> rewritten = value.undo_dots()
    fn undo_dots(self_: PyRef<'_, Self>, py: Python<'_>) -> PyResult<Py<Self>> {
        let value = self_.value.undo_dots().map_err(Self::inference_error)?;
        let (name, args) = Self::transformed_descriptor(&self_, value.expression());
        Self::from_shared(py, Arc::new(value), name, args)
    }

    /// Expose collected chain factors with fresh compatible dummy indices.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> r = sp.Representation.euc(2)
    /// >>> A = sp.TensorName("docs::A")(r, r)
    /// >>> value = sp.chain(r("i"), r("j"), A, A)
    /// >>> opened = value.undo_chain()
    fn undo_chain(self_: PyRef<'_, Self>, py: Python<'_>) -> PyResult<Py<Self>> {
        let value = self_.value.undo_chain().map_err(Self::inference_error)?;
        let (name, args) = Self::transformed_descriptor(&self_, value.expression());
        Self::from_shared(py, Arc::new(value), name, args)
    }

    /// Open traces as closed chains without evaluating them.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> r = sp.Representation.euc(2)
    /// >>> A = sp.TensorName("docs::A")(r, r)
    /// >>> value = sp.trace(r, A, A)
    /// >>> opened = value.undo_trace()
    fn undo_trace(self_: PyRef<'_, Self>, py: Python<'_>) -> PyResult<Py<Self>> {
        let value = self_.value.undo_trace().map_err(Self::inference_error)?;
        let (name, args) = Self::transformed_descriptor(&self_, value.expression());
        Self::from_shared(py, Arc::new(value), name, args)
    }

    /// Materialize unresolved ports and parse this expression as a tensor network.
    ///
    /// When `library` is omitted, use the built-in four-dimensional HEP and SU(3) library.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> r = sp.Representation.euc(2)
    /// >>> A = sp.TensorName("docs::A")(r, r)
    /// >>> tensor = sp.Tensor.dense(A, [1.0, 2.0, 3.0, 4.0])
    /// >>> library = sp.TensorLibrary()
    /// >>> library.register(tensor)
    /// >>> network = A.to_network(library)
    ///
    /// Parameters
    /// ----------
    /// library : TensorLibrary, optional
    ///     Component definitions for named tensors. None uses the built-in four-dimensional HEP and SU(3) library; supply a TensorLibrary for your own tensors.
    #[pyo3(signature = (library = None))]
    fn to_network(
        self_: PyRef<'_, Self>,
        py: Python<'_>,
        library: Option<&SpensorLibrary>,
    ) -> PyResult<SpensoNet> {
        let expression = Self::from_known_parts(
            py,
            self_.atom().clone(),
            self_.interface().clone(),
            self_.name,
            self_.name_args.clone(),
        )?;
        SpensoNet::from_arithmetic(ArithmeticStructure::Tensor(expression), library)
            .map_err(|error| PyRuntimeError::new_err(error.to_string()))
    }

    /// Execute the expression and return its component Tensor in logical axis order.
    ///
    /// Unregistered names acquire symbolic components. This uses the same execution
    /// path as to_network().execute() and does not mutate the expression or library.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> r = sp.Representation.euc(2)
    /// >>> A = sp.TensorName("docs::A")(r, r)
    /// >>> tensor = sp.Tensor.dense(A, [1.0, 2.0, 3.0, 4.0])
    /// >>> library = sp.TensorLibrary()
    /// >>> library.register(tensor)
    /// >>> result = A.to_tensor(library)
    ///
    /// Parameters
    /// ----------
    /// library : TensorLibrary, optional
    ///     Component definitions for named tensors. None uses the built-in four-dimensional HEP and SU(3) library; supply a TensorLibrary for your own tensors.
    /// function_library : TensorFunctionLibrary, optional
    ///     Callbacks for custom elementwise functions during component evaluation. None uses only the available built-in functions.
    #[pyo3(signature = (library=None, *, function_library=None))]
    fn to_tensor(
        self_: PyRef<'_, Self>,
        py: Python<'_>,
        library: Option<&SpensorLibrary>,
        function_library: Option<&crate::library::SpensorFunctionLibrary>,
    ) -> PyResult<Spensor> {
        let mut network = Self::to_network(self_, py, library)?;
        network.execute(library, function_library, None, ExecutionMode::All)?;
        network.result_tensor(library)
    }

    /// Return component values in logical row-major interface order.
    ///
    /// Parse and execute through the same path as `to_network()`. Tensors absent
    /// from the library acquire symbolic components, with the same identities
    /// used when they occur in a larger network. Registered tensors use their
    /// library values; compound expressions are contracted before extraction.
    /// Abstract index labels do not change a named tensor's component identities.
    ///
    /// The result is a flat list (one element for a scalar). Dimensions must be
    /// concrete. This does not mutate the expression or register generated data.
    /// When `library` is omitted, use the default HEP library; pass
    /// `TensorLibrary.hep_lib_atom()` to select its atom-valued variant.
    /// Values retain the Python types returned by `TensorNetwork.result_tensor()`.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> r = sp.Representation.euc(2)
    /// >>> A = sp.TensorName("docs::A")(r, r)
    /// >>> tensor = sp.Tensor.dense(A, [1.0, 2.0, 3.0, 4.0])
    /// >>> library = sp.TensorLibrary()
    /// >>> library.register(tensor)
    /// >>> entries = A.components(library)
    ///
    /// Parameters
    /// ----------
    /// library : TensorLibrary, optional
    ///     Component definitions for named tensors. None uses the built-in four-dimensional HEP and SU(3) library; supply a TensorLibrary for your own tensors.
    #[pyo3(signature = (library = None))]
    #[gen_stub(override_return_type(
        type_repr = "builtins.list[Expression | builtins.complex | builtins.float]"
    ))]
    fn components(
        self_: PyRef<'_, Self>,
        py: Python<'_>,
        library: Option<&SpensorLibrary>,
    ) -> PyResult<Py<PyAny>> {
        Self::to_tensor(self_, py, library, None)?
            .__getitem__(SliceOrIntOrExpanded::Slice(PySlice::full(py)))
    }

    /// Fill the unresolved external ports with `indices` in interface order.
    ///
    /// Pass `AUTO` to leave a port unresolved. Pass `cook_indices=sp.CookSettings.reversible()` to encode nested
    /// symbolic index payloads before insertion. Repeated compatible indices contract their
    /// ports, in which case the result no longer carries the original stored-data identity.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> r = sp.Representation.euc(2)
    /// >>> A = sp.TensorName("docs::A")(r, r)
    /// >>> indexed = A.index("i", "j")
    ///
    /// Parameters
    /// ----------
    /// indices : index label, Slot, or AUTO
    ///     One index per selected port, in logical axis order. Accepts integers, strings, expressions, typed Slots, and AUTO; AUTO leaves that port unchanged.
    /// cook_indices : CookSettings, optional
    ///     None leaves index payloads unchanged. CookSettings selects encoding, filters, and output tags; CookSettings.indices() flattens nested index payloads.
    #[pyo3(signature = (*indices, cook_indices = None))]
    fn index(
        self_: PyRef<'_, Self>,
        py: Python<'_>,
        indices: &Bound<'_, PyTuple>,
        cook_indices: Option<&PyCookSettings>,
    ) -> PyResult<Py<Self>> {
        let replacements = Self::index_replacements(self_.interface(), indices, cook_indices)?;
        let value = Self::structured(&self_)
            .indexed(&replacements)
            .map_err(|error| PyValueError::new_err(error.to_string()))?;
        let (name, args) = if value.rank() == self_.interface().canonical().order() {
            (self_.name, self_.name_args.clone())
        } else {
            (None, Vec::new())
        };
        Self::from_shared(py, Arc::new(value), name, args)
    }

    /// Assign all external ports in logical order, including already indexed ports.
    /// AUTO keeps the current port. Repeated compatible indices contract.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> r = sp.Representation.euc(2)
    /// >>> A = sp.TensorName("docs::A")(r, r)
    /// >>> indexed = A("i", "j")
    /// >>> renamed = indexed.reindex("mu", sp.AUTO)
    ///
    /// Parameters
    /// ----------
    /// indices : index label, Slot, or AUTO
    ///     One index per selected port, in logical axis order. Accepts integers, strings, expressions, typed Slots, and AUTO; AUTO leaves that port unchanged.
    /// cook_indices : CookSettings, optional
    ///     None leaves index payloads unchanged. CookSettings selects encoding, filters, and output tags; CookSettings.indices() flattens nested index payloads.
    #[pyo3(signature = (*indices, cook_indices=None))]
    fn reindex(
        self_: PyRef<'_, Self>,
        py: Python<'_>,
        indices: &Bound<'_, PyTuple>,
        cook_indices: Option<&PyCookSettings>,
    ) -> PyResult<Py<Self>> {
        let positions = (0..self_.interface().canonical().order()).collect::<Vec<_>>();
        let replacements =
            Self::port_replacements(self_.interface(), &positions, indices, cook_indices)?;
        let value = Self::structured(&self_)
            .reindex_interface_ports(&replacements)
            .map_err(|e| PyValueError::new_err(e.to_string()))?;
        let (atom, interface) = value.into_parts();
        Self::from_known_parts(py, atom, interface, self_.name, self_.name_args.clone())
    }

    /// Rename external indices simultaneously without changing rank or capturing dummy indices.
    /// Keys are index labels or typed Slots; typed keys disambiguate representations.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> r = sp.Representation.euc(2)
    /// >>> A = sp.TensorName("docs::A")(r, r)
    /// >>> indexed = A("i", "j")
    /// >>> swapped = indexed.rename_indices({"i": "j", "j": "i"})
    ///
    /// Parameters
    /// ----------
    /// mapping : dict[int | str | Expression | Slot, index label, Slot, or AUTO]
    ///     Old-to-new external index labels, applied simultaneously. Typed Slot keys disambiguate representations. Renaming cannot contract ports or capture dummy indices; use reindex() for contractions.
    /// cook_indices : CookSettings, optional
    ///     None leaves index payloads unchanged. CookSettings selects encoding, filters, and output tags; CookSettings.indices() flattens nested index payloads.
    #[pyo3(signature = (mapping, *, cook_indices=None))]
    fn rename_indices(
        self_: PyRef<'_, Self>,
        py: Python<'_>,
        mapping: &Bound<'_, PyDict>,
        cook_indices: Option<&PyCookSettings>,
    ) -> PyResult<Py<Self>> {
        let replacements = Self::named_replacements(self_.interface(), mapping, cook_indices)?;
        let value = Self::structured(&self_)
            .reindex_interface_ports(&replacements)
            .map_err(|e| PyValueError::new_err(e.to_string()))?;
        let (atom, interface) = value.into_parts();
        let result =
            Self::from_known_parts(py, atom, interface, self_.name, self_.name_args.clone())?;
        if result.borrow(py).interface().canonical().order()
            != self_.interface().canonical().order()
        {
            return Err(PyValueError::new_err(
                "renaming would contract external ports; use reindex() to request a contraction",
            ));
        }
        Ok(result)
    }

    /// Return a tensor with external axes in the specified logical order.
    /// Unresolved ports acquire fresh indices so identical representations remain distinct.
    /// Use reindex() to assign preferred labels after permutation.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> r = sp.Representation.euc(2)
    /// >>> A = sp.TensorName("docs::A")(r, r)
    /// >>> transposed = A.permute_axes([1, 0])
    ///
    /// Parameters
    /// ----------
    /// axes : Sequence[int]
    ///     Permutation of all axis positions, starting at zero. For a matrix, [1, 0] transposes the logical axes.
    fn permute_axes(
        self_: PyRef<'_, Self>,
        py: Python<'_>,
        axes: Vec<usize>,
    ) -> PyResult<Py<Self>> {
        let value = Self::structured(&self_)
            .permuted(&axes)
            .map_err(|error| PyValueError::new_err(error.to_string()))?;
        let (atom, interface) = value.into_parts();
        Self::from_known_parts(py, atom, interface, self_.name, self_.name_args.clone())
    }

    /// Fill the unresolved external ports with `indices` in interface order.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> r = sp.Representation.euc(2)
    /// >>> A = sp.TensorName("docs::A")(r, r)
    /// >>> indexed = A("i", "j")
    ///
    /// Parameters
    /// ----------
    /// indices : index label, Slot, or AUTO
    ///     One index per selected port, in logical axis order. Accepts integers, strings, expressions, typed Slots, and AUTO; AUTO leaves that port unchanged.
    /// cook_indices : CookSettings, optional
    ///     None leaves index payloads unchanged. CookSettings selects encoding, filters, and output tags; CookSettings.indices() flattens nested index payloads.
    #[pyo3(signature = (*indices, cook_indices = None))]
    fn __call__(
        self_: PyRef<'_, Self>,
        py: Python<'_>,
        indices: &Bound<'_, PyTuple>,
        cook_indices: Option<&PyCookSettings>,
    ) -> PyResult<Py<Self>> {
        Self::index(self_, py, indices, cook_indices)
    }

    /// Negate every tensor coefficient, preserving the external axes.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> r = sp.Representation.euc(2)
    /// >>> A = sp.TensorName("docs::A")(r, r)
    /// >>> negative = -A
    fn __neg__(self_: PyRef<'_, Self>, py: Python<'_>) -> PyResult<Py<Self>> {
        Self::from_known_parts(
            py,
            Atom::num(-1) * self_.atom(),
            self_.interface().clone(),
            None,
            Vec::new(),
        )
    }

    /// Add compatible tensors, preserving their logical external axes.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> r = sp.Representation.euc(2)
    /// >>> A = sp.TensorName("docs::A")(r, r)
    /// >>> result = A + A
    ///
    /// Parameters
    /// ----------
    /// rhs : TensorExpression | Expression or scalar
    ///     Right operand, with a compatible tensor interface for addition or subtraction. A denominator must be scalar.
    #[gen_stub(skip)]
    fn __add__(
        self_: PyRef<'_, Self>,
        py: Python<'_>,
        rhs: &Bound<'_, PyAny>,
    ) -> PyResult<TensorDispatch> {
        if let Some(right) = concrete_network(rhs)? {
            return Self::promoted_network(&self_, py)?
                .add_network(right, false)
                .and_then(|network| Py::new(py, network))
                .map(TensorDispatch::Network);
        }
        let left = Self::structured(&self_);
        let right = match TensorOperand::extract(rhs)? {
            TensorOperand::Structured(right) => {
                return left
                    .try_add(&right)
                    .map_err(|error| PyValueError::new_err(error.to_string()))
                    .and_then(|value| Self::from_structured(py, value))
                    .map(TensorDispatch::Expression);
            }
            TensorOperand::Scalar(right) if right.as_view().is_zero() => {
                return Self::from_shared(
                    py,
                    Arc::clone(&self_.value),
                    self_.name,
                    self_.name_args.clone(),
                )
                .map(TensorDispatch::Expression);
            }
            TensorOperand::Scalar(right) if left.is_scalar() => right,
            TensorOperand::Scalar(_) => {
                return Err(PyValueError::new_err(
                    "cannot add a scalar expression to a non-scalar tensor",
                ));
            }
        };
        let (left, interface) = (left.expression(), left.structure().clone());
        Self::from_known_parts(py, left + right.as_ref(), interface, None, Vec::new())
            .map(TensorDispatch::Expression)
    }

    /// Add compatible tensors, preserving their logical external axes.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> r = sp.Representation.euc(2)
    /// >>> A = sp.TensorName("docs::A")(r, r)
    /// >>> result = A + A
    ///
    /// Parameters
    /// ----------
    /// lhs : TensorExpression | Expression or scalar
    ///     Left operand. Addition and subtraction require compatible tensor interfaces.
    #[gen_stub(skip)]
    fn __radd__(
        self_: PyRef<'_, Self>,
        py: Python<'_>,
        lhs: &Bound<'_, PyAny>,
    ) -> PyResult<TensorDispatch> {
        if let Some(left) = concrete_network(lhs)? {
            return left
                .add_network(Self::promoted_network(&self_, py)?, false)
                .and_then(|network| Py::new(py, network))
                .map(TensorDispatch::Network);
        }
        Self::__add__(self_, py, lhs)
    }

    /// Subtract compatible tensors, preserving their logical external axes.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> r = sp.Representation.euc(2)
    /// >>> A = sp.TensorName("docs::A")(r, r)
    /// >>> result = A - A
    ///
    /// Parameters
    /// ----------
    /// rhs : TensorExpression | Expression or scalar
    ///     Right operand, with a compatible tensor interface for addition or subtraction. A denominator must be scalar.
    #[gen_stub(skip)]
    fn __sub__(
        self_: PyRef<'_, Self>,
        py: Python<'_>,
        rhs: &Bound<'_, PyAny>,
    ) -> PyResult<TensorDispatch> {
        if let Some(right) = concrete_network(rhs)? {
            return Self::promoted_network(&self_, py)?
                .add_network(right, true)
                .and_then(|network| Py::new(py, network))
                .map(TensorDispatch::Network);
        }
        let left = Self::structured(&self_);
        let right = match TensorOperand::extract(rhs)? {
            TensorOperand::Structured(right) => {
                return left
                    .try_sub(&right)
                    .map_err(|error| PyValueError::new_err(error.to_string()))
                    .and_then(|value| Self::from_structured(py, value))
                    .map(TensorDispatch::Expression);
            }
            TensorOperand::Scalar(right) if right.as_view().is_zero() => {
                return Self::from_shared(
                    py,
                    Arc::clone(&self_.value),
                    self_.name,
                    self_.name_args.clone(),
                )
                .map(TensorDispatch::Expression);
            }
            TensorOperand::Scalar(right) if left.is_scalar() => right,
            TensorOperand::Scalar(_) => {
                return Err(PyValueError::new_err(
                    "cannot subtract a scalar expression from a non-scalar tensor",
                ));
            }
        };
        let (left, interface) = (left.expression(), left.structure().clone());
        Self::from_known_parts(py, left - right.as_ref(), interface, None, Vec::new())
            .map(TensorDispatch::Expression)
    }

    /// Subtract compatible tensors, preserving their logical external axes.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> r = sp.Representation.euc(2)
    /// >>> A = sp.TensorName("docs::A")(r, r)
    /// >>> result = A - A
    ///
    /// Parameters
    /// ----------
    /// lhs : TensorExpression | Expression or scalar
    ///     Left operand. Addition and subtraction require compatible tensor interfaces.
    #[gen_stub(skip)]
    fn __rsub__(
        self_: PyRef<'_, Self>,
        py: Python<'_>,
        lhs: &Bound<'_, PyAny>,
    ) -> PyResult<TensorDispatch> {
        if let Some(left) = concrete_network(lhs)? {
            return left
                .add_network(Self::promoted_network(&self_, py)?, true)
                .and_then(|network| Py::new(py, network))
                .map(TensorDispatch::Network);
        }
        let right = Self::structured(&self_);
        let left = match TensorOperand::extract(lhs)? {
            TensorOperand::Structured(left) => {
                return left
                    .try_sub(right)
                    .map_err(|error| PyValueError::new_err(error.to_string()))
                    .and_then(|value| Self::from_structured(py, value))
                    .map(TensorDispatch::Expression);
            }
            TensorOperand::Scalar(left) if left.as_view().is_zero() => {
                return Self::__neg__(self_, py).map(TensorDispatch::Expression);
            }
            TensorOperand::Scalar(left) if right.is_scalar() => left,
            TensorOperand::Scalar(_) => {
                return Err(PyValueError::new_err(
                    "cannot subtract a non-scalar tensor from a scalar expression",
                ));
            }
        };
        let (right, interface) = (right.expression(), right.structure().clone());
        Self::from_known_parts(py, left.as_ref() - right, interface, None, Vec::new())
            .map(TensorDispatch::Expression)
    }

    /// Multiply tensor operands, contracting repeated compatible explicit indices.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> r = sp.Representation.euc(2)
    /// >>> A = sp.TensorName("docs::A")(r, r)
    /// >>> result = A * 2
    ///
    /// Parameters
    /// ----------
    /// rhs : TensorExpression | Expression or scalar
    ///     Right operand, with a compatible tensor interface for addition or subtraction. A denominator must be scalar.
    #[gen_stub(skip)]
    fn __mul__(
        self_: PyRef<'_, Self>,
        py: Python<'_>,
        rhs: &Bound<'_, PyAny>,
    ) -> PyResult<TensorDispatch> {
        if let Some(right) = concrete_network(rhs)? {
            return Self::promoted_network(&self_, py)?
                .multiply_network(right)
                .and_then(|network| Py::new(py, network))
                .map(TensorDispatch::Network);
        }
        let left = Self::structured(&self_);
        match TensorOperand::extract(rhs)? {
            TensorOperand::Structured(right) => left
                .multiply(&right)
                .map_err(|error| PyValueError::new_err(error.to_string()))
                .and_then(|value| Self::from_structured(py, value))
                .map(TensorDispatch::Expression),
            TensorOperand::Scalar(right) => Self::from_known_parts(
                py,
                left.expression() * right.as_ref(),
                left.structure().clone(),
                None,
                Vec::new(),
            )
            .map(TensorDispatch::Expression),
        }
    }

    /// Multiply tensor operands, contracting repeated compatible explicit indices.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> r = sp.Representation.euc(2)
    /// >>> A = sp.TensorName("docs::A")(r, r)
    /// >>> result = 2 * A
    ///
    /// Parameters
    /// ----------
    /// lhs : TensorExpression | Expression or scalar
    ///     Left operand. Addition and subtraction require compatible tensor interfaces.
    #[gen_stub(skip)]
    fn __rmul__(
        self_: PyRef<'_, Self>,
        py: Python<'_>,
        lhs: &Bound<'_, PyAny>,
    ) -> PyResult<TensorDispatch> {
        if let Some(left) = concrete_network(lhs)? {
            return left
                .multiply_network(Self::promoted_network(&self_, py)?)
                .and_then(|network| Py::new(py, network))
                .map(TensorDispatch::Network);
        }
        let right = Self::structured(&self_);
        match TensorOperand::extract(lhs)? {
            TensorOperand::Structured(left) => left
                .multiply(right)
                .map_err(|error| PyValueError::new_err(error.to_string()))
                .and_then(|value| Self::from_structured(py, value))
                .map(TensorDispatch::Expression),
            TensorOperand::Scalar(left) => Self::from_known_parts(
                py,
                left.as_ref() * right.expression(),
                right.structure().clone(),
                None,
                Vec::new(),
            )
            .map(TensorDispatch::Expression),
        }
    }

    /// Divide tensor coefficients by a scalar denominator.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> r = sp.Representation.euc(2)
    /// >>> A = sp.TensorName("docs::A")(r, r)
    /// >>> result = A / 2
    ///
    /// Parameters
    /// ----------
    /// rhs : TensorExpression | Expression or scalar
    ///     Right operand, with a compatible tensor interface for addition or subtraction. A denominator must be scalar.
    #[gen_stub(skip)]
    fn __truediv__(
        self_: PyRef<'_, Self>,
        py: Python<'_>,
        rhs: &Bound<'_, PyAny>,
    ) -> PyResult<TensorDispatch> {
        if let Some(right) = concrete_network(rhs)? {
            return Self::promoted_network(&self_, py)?
                .__truediv__(ConvertibleToSpensoNet(right))
                .and_then(|network| Py::new(py, network))
                .map(TensorDispatch::Network);
        }
        let denominator = TensorOperand::extract(rhs)?.into_structured()?;
        let value = self_
            .value
            .divide(&denominator)
            .map_err(Self::inference_error)?;
        Self::from_structured(py, value).map(TensorDispatch::Expression)
    }

    /// Divide a scalar by a scalar tensor expression.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> scalar = sp.TensorExpression(2)
    /// >>> result = 1 / scalar
    ///
    /// Parameters
    /// ----------
    /// lhs : TensorExpression | Expression or scalar
    ///     Left operand. Addition and subtraction require compatible tensor interfaces.
    #[gen_stub(skip)]
    fn __rtruediv__(
        self_: PyRef<'_, Self>,
        py: Python<'_>,
        lhs: &Bound<'_, PyAny>,
    ) -> PyResult<TensorDispatch> {
        if let Some(left) = concrete_network(lhs)? {
            return Self::promoted_network(&self_, py)?
                .__rtruediv__(ConvertibleToSpensoNet(left))
                .and_then(|network| Py::new(py, network))
                .map(TensorDispatch::Network);
        }
        let numerator = TensorOperand::extract(lhs)?.into_structured()?;
        let value = numerator
            .divide(&self_.value)
            .map_err(Self::inference_error)?;
        Self::from_structured(py, value).map(TensorDispatch::Expression)
    }

    /// Raise a scalar tensor expression to a power. Non-scalar tensor powers are rejected.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> from symbolica import E, S
    /// >>> x = S("docs::x")
    /// >>> scalar = sp.TensorExpression(x)
    /// >>> result = scalar**2
    ///
    /// Parameters
    /// ----------
    /// exponent : TensorExpression | Expression or scalar
    ///     Scalar exponent. Tensor-valued exponents are not supported.
    /// modulo : Any, optional
    ///     Must be None; modular exponentiation is not supported.
    fn __pow__(
        self_: PyRef<'_, Self>,
        py: Python<'_>,
        exponent: &Bound<'_, PyAny>,
        modulo: Option<&Bound<'_, PyAny>>,
    ) -> PyResult<Py<Self>> {
        if modulo.is_some() {
            return Err(PyTypeError::new_err(
                "modular tensor powers are not supported",
            ));
        }
        let exponent = TensorOperand::extract(exponent)?.into_structured()?;
        let value = self_.value.pow(&exponent).map_err(Self::inference_error)?;
        Self::from_structured(py, value)
    }

    /// Raise a scalar tensor expression to a power. Non-scalar tensor powers are rejected.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> from symbolica import E, S
    /// >>> x = S("docs::x")
    /// >>> scalar = sp.TensorExpression(x)
    /// >>> result = 2**scalar
    ///
    /// Parameters
    /// ----------
    /// base : TensorExpression | Expression or scalar
    ///     Scalar base of the power.
    /// modulo : Any, optional
    ///     Must be None; modular exponentiation is not supported.
    fn __rpow__(
        self_: PyRef<'_, Self>,
        py: Python<'_>,
        base: &Bound<'_, PyAny>,
        modulo: Option<&Bound<'_, PyAny>>,
    ) -> PyResult<Py<Self>> {
        if modulo.is_some() {
            return Err(PyTypeError::new_err(
                "modular tensor powers are not supported",
            ));
        }
        let base = TensorOperand::extract(base)?.into_structured()?;
        let value = base.pow(&self_.value).map_err(Self::inference_error)?;
        Self::from_structured(py, value)
    }

    /// Compare normalized expressions and their ordered tensor interfaces for equality.
    /// For a comparison to a plain Symbolica expression, first use `.to_expression()`.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> r = sp.Representation.euc(2)
    /// >>> A = sp.TensorName("docs::A")(r, r)
    /// >>> assert A == sp.TensorExpression(A)
    ///
    /// Parameters
    /// ----------
    /// other : object
    ///     TensorExpression to compare. Other types can provide their own reverse comparison.
    fn __eq__(&self, py: Python<'_>, other: &Bound<'_, PyAny>) -> PyResult<Py<PyAny>> {
        let Ok(other) = other.extract::<PyRef<'_, Self>>() else {
            return Ok(py.NotImplemented());
        };
        Ok(self
            .value
            .same_value(&other.value)
            .into_pyobject(py)?
            .to_owned()
            .into_any()
            .unbind())
    }

    /// Test whether normalized expressions or their ordered tensor interfaces differ.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> r = sp.Representation.euc(2)
    /// >>> A = sp.TensorName("docs::A")(r, r)
    /// >>> assert A != 2 * A
    ///
    /// Parameters
    /// ----------
    /// other : object
    ///     TensorExpression to compare. Other types can provide their own reverse comparison.
    fn __ne__(&self, py: Python<'_>, other: &Bound<'_, PyAny>) -> PyResult<Py<PyAny>> {
        let Ok(other) = other.extract::<PyRef<'_, Self>>() else {
            return Ok(py.NotImplemented());
        };
        Ok((!self.value.same_value(&other.value))
            .into_pyobject(py)?
            .to_owned()
            .into_any()
            .unbind())
    }

    /// Hash the immutable tensor value for use in dictionaries and sets.
    /// Equal tensor expressions have equal hashes, including their ordered interfaces.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> r = sp.Representation.euc(2)
    /// >>> A = sp.TensorName("docs::A")(r, r)
    /// >>> labels = {A: "matrix A"}
    /// >>> assert labels[sp.TensorExpression(A)] == "matrix A"
    fn __hash__(&self) -> isize {
        use std::hash::{DefaultHasher, Hasher};
        let mut state = DefaultHasher::new();
        self.value.hash_value(&mut state);
        let hash = state.finish() as isize;
        if hash == -1 { -2 } else { hash }
    }

    /// Test whether the symbolic tensor expression is nonzero.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> r = sp.Representation.euc(2)
    /// >>> A = sp.TensorName("docs::A")(r, r)
    /// >>> nonzero = bool(A)
    fn __bool__(&self) -> bool {
        !self.value.is_zero()
    }

    /// Form an outer product, keeping both operands’ ports free even when their labels coincide.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> r = sp.Representation.euc(2)
    /// >>> A = sp.TensorName("docs::A")(r, r)
    /// >>> result = A.outer(A)
    ///
    /// Parameters
    /// ----------
    /// rhs : TensorExpression | Expression
    ///     Other tensor, network, or scalar operand. Concrete component data or a network makes the result a TensorNetwork.
    #[gen_stub(skip)]
    fn outer(
        self_: PyRef<'_, Self>,
        py: Python<'_>,
        rhs: &Bound<'_, PyAny>,
    ) -> PyResult<TensorDispatch> {
        if let Some(right) = concrete_network(rhs)? {
            return Self::promoted_network(&self_, py)?
                .outer(ConvertibleToSpensoNet(right))
                .and_then(|network| Py::new(py, network))
                .map(TensorDispatch::Network);
        }
        let left = Self::structured(&self_);
        let TensorOperand::Structured(right) = TensorOperand::extract(rhs)? else {
            return Err(PyTypeError::new_err("outer() requires a tensor operand"));
        };
        left.outer(&right)
            .map_err(|error| PyValueError::new_err(error.to_string()))
            .and_then(|value| Self::from_structured(py, value))
            .map(TensorDispatch::Expression)
    }

    /// Contract compatible indices, retaining shared subexpressions in the result.
    /// Use `ContractSettings` to select representations or control chain and trace collection.
    /// Call `.expand()` to materialize the factored result.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> r = sp.Representation.euc(2)
    /// >>> metric = sp.TensorExpression.g(r)("i", "j")
    /// >>> vector = sp.TensorName.vector("docs::p")(r("j"))
    /// >>> contracted = (metric * vector).contract()
    /// >>> explicit = contracted.expand()
    ///
    /// Parameters
    /// ----------
    /// settings : ContractSettings, optional
    ///     ContractSettings selecting contractions and collection. None enables the default contraction policy.
    #[pyo3(signature=(settings=None))]
    fn contract(
        &self,
        settings: Option<&crate::simplification::PyContractSettings>,
    ) -> PyResult<AliasedTensorExpression> {
        let settings = settings.cloned().unwrap_or_default();
        let value = self
            .value
            .contract(settings.rust())
            .map_err(Self::inference_error)?;
        Ok(AliasedTensorExpression::from_parts(
            Arc::new(value),
            self.name,
            self.name_args.clone(),
        ))
    }

    /// Contract one chosen pair of axes and retain the other external ports in logical order.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> r = sp.Representation.euc(2)
    /// >>> A = sp.TensorName("docs::A")(r, r)
    /// >>> result = A.contract_ports(A, left=1, right=0)
    ///
    /// Parameters
    /// ----------
    /// rhs : TensorExpression | Expression
    ///     Other tensor, network, or scalar operand. Concrete component data or a network makes the result a TensorNetwork.
    /// left : int
    ///     Zero-based axis position in this operand to contract.
    /// right : int
    ///     Zero-based axis position in rhs to contract with the selected left axis.
    #[pyo3(signature = (rhs, *, left, right))]
    #[gen_stub(skip)]
    fn contract_ports(
        self_: PyRef<'_, Self>,
        py: Python<'_>,
        rhs: &Bound<'_, PyAny>,
        left: usize,
        right: usize,
    ) -> PyResult<TensorDispatch> {
        if let Some(rhs) = concrete_network(rhs)? {
            return Self::promoted_network(&self_, py)?
                .contract_ports(ConvertibleToSpensoNet(rhs), left, right)
                .and_then(|network| Py::new(py, network))
                .map(TensorDispatch::Network);
        }
        let value = Self::structured(&self_);
        let TensorOperand::Structured(rhs) = TensorOperand::extract(rhs)? else {
            return Err(PyTypeError::new_err(
                "contract_ports() requires a tensor operand",
            ));
        };
        value
            .contract_ports(&rhs, &[composition::PortPair { left, right }])
            .map_err(|error| PyValueError::new_err(error.to_string()))
            .and_then(|value| Self::from_structured(py, value))
            .map(TensorDispatch::Expression)
    }

    /// Compose two explicitly selected matrix channels, retaining the uncontracted endpoints.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> r = sp.Representation.euc(2)
    /// >>> A = sp.TensorName("docs::A")(r, r)
    /// >>> result = A.compose(A, left=(0, 1), right=(0, 1))
    ///
    /// Parameters
    /// ----------
    /// rhs : TensorExpression | Expression
    ///     Other tensor, network, or scalar operand. Concrete component data or a network makes the result a TensorNetwork.
    /// left : tuple[int, int]
    ///     (input_axis, output_axis) channel in this operand.
    /// right : tuple[int, int]
    ///     (input_axis, output_axis) channel in rhs. The selected connecting ports must have compatible representations.
    #[pyo3(signature = (rhs, *, left, right))]
    #[gen_stub(skip)]
    fn compose(
        self_: PyRef<'_, Self>,
        py: Python<'_>,
        rhs: &Bound<'_, PyAny>,
        left: (usize, usize),
        right: (usize, usize),
    ) -> PyResult<TensorDispatch> {
        if let Some(rhs) = concrete_network(rhs)? {
            return Self::promoted_network(&self_, py)?
                .compose(ConvertibleToSpensoNet(rhs), left, right)
                .and_then(|network| Py::new(py, network))
                .map(TensorDispatch::Network);
        }
        let value = Self::structured(&self_);
        let TensorOperand::Structured(rhs) = TensorOperand::extract(rhs)? else {
            return Err(PyTypeError::new_err("compose() requires a tensor operand"));
        };
        SymbolicTensor::compose(
            value,
            &rhs,
            composition::MatrixChannel {
                input: left.0,
                output: left.1,
            },
            composition::MatrixChannel {
                input: right.0,
                output: right.1,
            },
        )
        .map_err(|error| PyValueError::new_err(error.to_string()))
        .and_then(|value| Self::from_structured(py, value))
        .map(TensorDispatch::Expression)
    }

    /// Close `channel`, or the unique matrix channel when it is omitted.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> r = sp.Representation.euc(2)
    /// >>> A = sp.TensorName("docs::A")(r, r)
    /// >>> trace = A.trace(channel=(0, 1))
    ///
    /// Parameters
    /// ----------
    /// channel : tuple[int, int], optional
    ///     Pair (input_axis, output_axis) identifying the matrix channel to close. None selects the unique compatible channel; ambiguous tensors require an explicit pair.
    #[pyo3(signature = (*, channel = None))]
    fn trace(
        self_: PyRef<'_, Self>,
        py: Python<'_>,
        channel: Option<(usize, usize)>,
    ) -> PyResult<Py<Self>> {
        let value = Self::structured(&self_);
        let traced = match channel {
            Some((input, output)) => {
                value.trace_ports(composition::MatrixChannel { input, output })
            }
            None => value.trace_unique(),
        };
        traced
            .map_err(|error| PyValueError::new_err(error.to_string()))
            .and_then(|value| Self::from_structured(py, value))
    }

    /// Format this structured expression using compact Spenso notation.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> r = sp.Representation.euc(2)
    /// >>> A = sp.TensorName("docs::A")(r, r)
    /// >>> text = A.format_tensor(settings=sp.DisplaySettings(index_style="raw", show_dimensions=True))
    ///
    /// Parameters
    /// ----------
    /// show_dimensions : bool, optional
    ///     Override whether representation dimensions are displayed. None uses settings; the default display settings hide dimensions.
    /// settings : DisplaySettings, optional
    ///     Optional DisplaySettings controlling notation, indices, and component layout. None uses the default presentation.
    #[pyo3(signature = (show_dimensions = None, *, settings = None))]
    fn format_tensor(
        self_: PyRef<'_, Self>,
        show_dimensions: Option<bool>,
        settings: Option<PyRef<'_, display::DisplaySettings>>,
    ) -> PyResult<String> {
        let settings = display::resolved_settings(show_dimensions, settings.as_deref());
        display::validate_plain_source_settings(&settings)?;
        Ok(display::format_structured_settings(
            Self::structured(&self_),
            &settings,
        ))
    }

    /// Export LaTeX with tensor notation, index alphabets, and model parameter labels.
    /// ``max_line_length`` wraps top-level sums using Symbolica's usual rules.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> r = sp.Representation.euc(2)
    /// >>> A = sp.TensorName("docs::A")(r, r)
    /// >>> latex_source = A("i", "j").to_latex(max_line_length=80)
    ///
    /// Parameters
    /// ----------
    /// max_line_length : int, optional
    ///     Optional target line length for the exported LaTeX source. None leaves line wrapping to the formatter.
    #[pyo3(signature = (max_line_length = None))]
    fn to_latex(self_: PyRef<'_, Self>, max_line_length: Option<usize>) -> String {
        display::structured_to_latex(Self::structured(&self_), false, max_line_length)
    }

    /// Format this structured expression as Typst math source.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> r = sp.Representation.euc(2)
    /// >>> A = sp.TensorName("docs::A")(r, r)
    /// >>> typst_source = A.to_typst(settings=sp.DisplaySettings(index_style="raw", show_dimensions=True))
    ///
    /// Parameters
    /// ----------
    /// show_dimensions : bool, optional
    ///     Override whether representation dimensions are displayed. None uses settings; the default display settings hide dimensions.
    /// settings : DisplaySettings, optional
    ///     Optional DisplaySettings controlling notation, indices, and component layout. None uses the default presentation.
    #[pyo3(signature = (show_dimensions = None, *, settings = None))]
    fn to_typst(
        self_: PyRef<'_, Self>,
        show_dimensions: Option<bool>,
        settings: Option<PyRef<'_, display::DisplaySettings>>,
    ) -> PyResult<String> {
        let settings = display::resolved_settings(show_dimensions, settings.as_deref());
        display::validate_typst_source_settings(&settings)?;
        Ok(display::structured_to_typst_with_settings(
            Self::structured(&self_),
            &settings,
        ))
    }

    /// Build Symbolica's rich display value, including semantic HTML when the
    /// optional ``gammaloop[typst-display]`` renderer is installed.
    ///
    /// ``notation_source`` is a trusted complete replacement for the bundled
    /// ``notation.typ`` module, not a style fragment. Typst executes it. When
    /// HTML cannot be rendered, the result retains its LaTeX and text forms.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> r = sp.Representation.euc(2)
    /// >>> A = sp.TensorName("docs::A")(r, r)
    /// >>> display_value = A.formatted(settings=sp.DisplaySettings(index_style="raw", show_dimensions=True))
    ///
    /// Parameters
    /// ----------
    /// show_dimensions : bool, optional
    ///     Override whether representation dimensions are displayed. None uses settings; the default display settings hide dimensions.
    /// settings : DisplaySettings, optional
    ///     Optional DisplaySettings controlling notation, indices, and component layout. None uses the default presentation.
    /// notation_source : str, optional
    ///     Optional custom Typst notation definitions. None uses the bundled tensor notation. Pass trusted Typst source, not untrusted user text.
    #[pyo3(signature = (show_dimensions = None, *, settings = None, notation_source = None))]
    fn formatted(
        self_: PyRef<'_, Self>,
        py: Python<'_>,
        show_dimensions: Option<bool>,
        settings: Option<PyRef<'_, display::DisplaySettings>>,
        notation_source: Option<String>,
    ) -> PythonFormattedOutput {
        let settings = display::resolved_settings(show_dimensions, settings.as_deref());
        display::format_structured_output_rich(
            py,
            Self::structured(&self_),
            &settings,
            notation_source.as_deref(),
        )
    }

    /// Compile this expression to semantic HTML with the optional Typst
    /// renderer.
    ///
    /// Install ``gammaloop[typst-display]`` to enable this method.
    /// ``notation_source``, when supplied, is trusted Typst code replacing the
    /// complete bundled ``notation.typ`` module.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> r = sp.Representation.euc(2)
    /// >>> A = sp.TensorName("docs::A")(r, r)
    /// >>> html = A.to_html(settings=sp.DisplaySettings(index_style="raw", show_dimensions=True))
    ///
    /// Parameters
    /// ----------
    /// show_dimensions : bool, optional
    ///     Override whether representation dimensions are displayed. None uses settings; the default display settings hide dimensions.
    /// settings : DisplaySettings, optional
    ///     Optional DisplaySettings controlling notation, indices, and component layout. None uses the default presentation.
    /// notation_source : str, optional
    ///     Optional custom Typst notation definitions. None uses the bundled tensor notation. Pass trusted Typst source, not untrusted user text.
    #[pyo3(signature = (show_dimensions = None, *, settings = None, notation_source = None))]
    fn to_html(
        self_: PyRef<'_, Self>,
        py: Python<'_>,
        show_dimensions: Option<bool>,
        settings: Option<PyRef<'_, display::DisplaySettings>>,
        notation_source: Option<String>,
    ) -> PyResult<String> {
        let settings = display::resolved_settings(show_dimensions, settings.as_deref());
        display::structured_to_html(
            py,
            Self::structured(&self_),
            &settings,
            notation_source.as_deref(),
        )
    }

    /// Compile this expression to SVG with the optional Typst renderer.
    ///
    /// Install ``gammaloop[typst-display]`` to enable this method.
    /// ``notation_source``, when supplied, is trusted Typst code replacing the
    /// complete bundled ``notation.typ`` module.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> r = sp.Representation.euc(2)
    /// >>> A = sp.TensorName("docs::A")(r, r)
    /// >>> svg = A.to_svg(settings=sp.DisplaySettings(index_style="raw", show_dimensions=True))
    ///
    /// Parameters
    /// ----------
    /// show_dimensions : bool, optional
    ///     Override whether representation dimensions are displayed. None uses settings; the default display settings hide dimensions.
    /// settings : DisplaySettings, optional
    ///     Optional DisplaySettings controlling notation, indices, and component layout. None uses the default presentation.
    /// notation_source : str, optional
    ///     Optional custom Typst notation definitions. None uses the bundled tensor notation. Pass trusted Typst source, not untrusted user text.
    #[pyo3(signature = (show_dimensions = None, *, settings = None, notation_source = None))]
    fn to_svg(
        self_: PyRef<'_, Self>,
        py: Python<'_>,
        show_dimensions: Option<bool>,
        settings: Option<PyRef<'_, display::DisplaySettings>>,
        notation_source: Option<String>,
    ) -> PyResult<String> {
        let settings = display::resolved_settings(show_dimensions, settings.as_deref());
        display::structured_to_svg(
            py,
            Self::structured(&self_),
            &settings,
            notation_source.as_deref(),
        )
    }

    /// Return a readable object description for inspection.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> r = sp.Representation.euc(2)
    /// >>> A = sp.TensorName("docs::A")(r, r)
    /// >>> text = repr(A)
    fn __repr__(self_: PyRef<'_, Self>) -> String {
        display::format_structured(Self::structured(&self_), false)
    }

    /// Return a readable text representation.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> r = sp.Representation.euc(2)
    /// >>> A = sp.TensorName("docs::A")(r, r)
    /// >>> text = str(A)
    fn __str__(self_: PyRef<'_, Self>) -> String {
        display::format_structured(Self::structured(&self_), false)
    }

    fn _repr_pretty_(
        self_: PyRef<'_, Self>,
        pretty: &Bound<'_, PyAny>,
        cycle: bool,
    ) -> PyResult<()> {
        let text = if cycle {
            "...".to_string()
        } else {
            display::format_structured(Self::structured(&self_), false)
        };
        pretty.call_method1("text", (text,))?;
        Ok(())
    }

    fn _repr_html_(self_: PyRef<'_, Self>, py: Python<'_>) -> Option<String> {
        let settings = display::DisplaySettings::default();
        display::structured_to_html(py, Self::structured(&self_), &settings, None).ok()
    }

    fn _repr_latex_(self_: PyRef<'_, Self>) -> String {
        Self::to_latex(self_, None)
    }
}

/// Restore tensor-aware dispatch after a base Symbolica transformation.
///
/// Examples
/// --------
/// >>> from symbolica.community import tensor as sp
/// >>> from symbolica import E, S
/// >>> r = sp.Representation.euc(2)
/// >>> A = sp.TensorName("docs::A")(r, r)
/// >>> x = S("docs::x")
/// >>> raw = (x * A).to_expression().derivative(x)
/// >>> tensor = sp.as_tensor(raw)
///
/// Parameters
/// ----------
/// expression : Any
///     TensorExpression to copy, or a Symbolica/scalar expression from which to infer the tensor interface. Useful after scalar algebra operations return a base Expression.
#[cfg_attr(
    feature = "python_stubgen",
    pyo3_stub_gen::derive::gen_stub_pyfunction(module = "symbolica.community.tensor")
)]
#[spenso_macros::track_usage(crate::record_usage)]
#[pyfunction]
pub fn as_tensor(py: Python<'_>, expression: &Bound<'_, PyAny>) -> PyResult<Py<TensorExpression>> {
    if let Ok(expression) = expression.extract::<PyRef<'_, TensorExpression>>() {
        return Py::new(py, expression.clone());
    }
    let expression = expression
        .extract::<ConvertibleToExpression>()
        .map_err(|_| PyTypeError::new_err("expected an Expression or TensorExpression"))?;
    TensorExpression::from_atom_interface(py, expression.to_expression().expr, None)
}

fn structured_operand(
    value: &Bound<'_, PyAny>,
    operation: &str,
) -> PyResult<SymbolicTensor<PartialStructure>> {
    match TensorOperand::extract(value)? {
        TensorOperand::Structured(value) => Ok(value),
        TensorOperand::Scalar(_) => Err(PyTypeError::new_err(format!(
            "{operation} requires structured tensor operands"
        ))),
    }
}

// Preserve actual Rust type names/imports for the symbolic function overloads.
#[cfg(feature = "python_stubgen")]
macro_rules! symbolic_function_overload {
    ($name:literal, $doc:literal, [$($argument:literal: $input:ty),*], [$($value:literal: $kind:ident),+]) => {
        submit! {
            pyo3_stub_gen::type_info::PyFunctionInfo {
                name: $name,
                parameters: &[
                    $(ParameterInfo {
                        name: $argument,
                        kind: ParameterKind::PositionalOrKeyword,
                        default: ParameterDefault::None,
                        type_info: <$input>::type_input,
                    },)*
                    $(ParameterInfo {
                        name: $value,
                        kind: ParameterKind::$kind,
                        default: ParameterDefault::None,
                        type_info: || TensorExpression::type_input() | PythonExpression::type_input(),
                    },)+
                ],
                r#return: TensorExpression::type_output,
                doc: $doc, module: Some("symbolica.community.tensor"),
                is_async: false, deprecated: None, type_ignored: None, is_overload: true,
                file: file!(), line: line!(), column: column!(), index: 0,
            }
        }
    };
}

#[cfg(feature = "python_stubgen")]
symbolic_function_overload!("dot", r##"Contract two rank-one tensors into the canonical dot form.

Examples
--------
>>> from symbolica.community import tensor as sp
>>> r = sp.Representation.euc(2)
>>> p = sp.TensorName.vector("docs::p")(r)
>>> q = sp.TensorName.vector("docs::q")(r)
>>> scalar_product = sp.dot(p, q)

Parameters
----------
left : TensorExpression | Expression
    Left rank-one tensor operand.
right : TensorExpression | Expression
    Right rank-one tensor operand with a compatible representation."##, [], ["left": PositionalOrKeyword, "right": PositionalOrKeyword]);

/// Contract two rank-one tensors into the canonical dot form.
///
/// Examples
/// --------
/// >>> from symbolica.community import tensor as sp
/// >>> r = sp.Representation.euc(2)
/// >>> p = sp.TensorName.vector("docs::p")(r)
/// >>> q = sp.TensorName.vector("docs::q")(r)
/// >>> scalar_product = sp.dot(p, q)
///
/// Parameters
/// ----------
/// left : TensorExpression | Expression
///     Left rank-one tensor operand.
/// right : TensorExpression | Expression
///     Right rank-one tensor operand with a compatible representation.
#[cfg_attr(
    feature = "python_stubgen",
    pyo3_stub_gen::derive::gen_stub_pyfunction(
        module = "symbolica.community.tensor",
        no_default_overload = true,
        python_overload = r#"
        import typing

        @overload
        def dot(left: typing.Union[Tensor, TensorNetwork], right: pyo3_stub_gen.RustType["ConvertibleToSpensoNet"]) -> pyo3_stub_gen.RustType["SpensoNet"]:
            """Contract two rank-one tensors into the canonical dot form.

            Examples
            --------
            >>> from symbolica.community import tensor as sp
            >>> r = sp.Representation.euc(2)
            >>> p = sp.TensorName.vector("docs::p")(r)
            >>> q = sp.TensorName.vector("docs::q")(r)
            >>> scalar_product = sp.dot(p, q)

            Parameters
            ----------
            left : Union[Tensor, TensorNetwork]
                Left rank-one tensor operand.
            right : pyo3_stub_gen.RustType['ConvertibleToSpensoNet']
                Right rank-one tensor operand with a compatible representation."""
        @overload
        def dot(left: pyo3_stub_gen.RustType["ConvertibleToSpensoNet"], right: typing.Union[Tensor, TensorNetwork]) -> pyo3_stub_gen.RustType["SpensoNet"]:
            """Contract two rank-one tensors into the canonical dot form.

            Examples
            --------
            >>> from symbolica.community import tensor as sp
            >>> r = sp.Representation.euc(2)
            >>> p = sp.TensorName.vector("docs::p")(r)
            >>> q = sp.TensorName.vector("docs::q")(r)
            >>> scalar_product = sp.dot(p, q)

            Parameters
            ----------
            left : pyo3_stub_gen.RustType['ConvertibleToSpensoNet']
                Left rank-one tensor operand.
            right : Union[Tensor, TensorNetwork]
                Right rank-one tensor operand with a compatible representation."""
        "#,
    )
)]
#[spenso_macros::track_usage(crate::record_usage)]
#[pyfunction]
fn dot(
    py: Python<'_>,
    left: &Bound<'_, PyAny>,
    right: &Bound<'_, PyAny>,
) -> PyResult<TensorDispatch> {
    if is_concrete_network(left) || is_concrete_network(right) {
        let left = left.extract::<ConvertibleToSpensoNet>()?.to_net();
        let right = right.extract::<ConvertibleToSpensoNet>()?;
        return left
            .dot(right)
            .and_then(|network| Py::new(py, network))
            .map(TensorDispatch::Network);
    }
    let left = structured_operand(left, "dot()")?;
    let right = structured_operand(right, "dot()")?;
    left.dot(&right)
        .map_err(|error| PyValueError::new_err(error.to_string()))
        .and_then(|value| TensorExpression::from_structured(py, value))
        .map(TensorDispatch::Expression)
}

// Python typing permits only one unbounded tuple unpack. A concrete first or last
// factor guarantees a network; an arbitrary mixed sequence needs the union fallback
// because it may be empty or contain only symbolic factors at runtime.
#[cfg(feature = "python_stubgen")]
symbolic_function_overload!("chain", r##"Build an explicitly-ended ordered tensor chain.

Examples
--------
>>> from symbolica.community import tensor as sp
>>> r = sp.Representation.euc(2)
>>> A = sp.TensorName("docs::A")(r, r)
>>> chain = sp.chain(r("i"), r("j"), A, A)

Parameters
----------
start_slot : Slot
    Typed Slot labelling the first free endpoint of the chain.
end_slot : Slot
    Typed Slot labelling the last free endpoint. Matching compatible endpoints close the chain.
factors : TensorExpression | Expression
    Ordered tensor factors. Symbolic operands produce a TensorExpression; any component tensor or network promotes the result to a TensorNetwork."##, ["start_slot": SpensoSlot, "end_slot": SpensoSlot], ["factors": VarPositional]);

/// Build an explicitly-ended ordered tensor chain.
///
/// Examples
/// --------
/// >>> from symbolica.community import tensor as sp
/// >>> r = sp.Representation.euc(2)
/// >>> A = sp.TensorName("docs::A")(r, r)
/// >>> chain = sp.chain(r("i"), r("j"), A, A)
///
/// Parameters
/// ----------
/// start_slot : Slot
///     Typed Slot labelling the first free endpoint of the chain.
/// end_slot : Slot
///     Typed Slot labelling the last free endpoint. Matching compatible endpoints close the chain.
/// factors : TensorExpression | Expression
///     Ordered tensor factors. Symbolic operands produce a TensorExpression; any component tensor or network promotes the result to a TensorNetwork.
#[cfg_attr(
    feature = "python_stubgen",
    pyo3_stub_gen::derive::gen_stub_pyfunction(
        module = "symbolica.community.tensor",
        no_default_overload = true,
        python_overload = r#"
        import typing

        @overload
        def chain(start_slot: pyo3_stub_gen.RustType["SpensoSlot"], end_slot: pyo3_stub_gen.RustType["SpensoSlot"], factor: typing.Union[Tensor, TensorNetwork], /, *factors: pyo3_stub_gen.RustType["ConvertibleToSpensoNet"]) -> pyo3_stub_gen.RustType["SpensoNet"]:
            """Build an explicitly-ended ordered tensor chain.

            Examples
            --------
            >>> from symbolica.community import tensor as sp
            >>> r = sp.Representation.euc(2)
            >>> A = sp.TensorName("docs::A")(r, r)
            >>> chain = sp.chain(r("i"), r("j"), A, A)

            Parameters
            ----------
            start_slot : pyo3_stub_gen.RustType['SpensoSlot']
                Typed Slot labelling the first free endpoint of the chain.
            end_slot : pyo3_stub_gen.RustType['SpensoSlot']
                Typed Slot labelling the last free endpoint. Matching compatible endpoints close the chain.
            factor : Union[Tensor, TensorNetwork]
                First tensor factor in the ordered product.
            factors : pyo3_stub_gen.RustType['ConvertibleToSpensoNet']
                Ordered tensor factors. Symbolic operands produce a TensorExpression; any component tensor or network promotes the result to a TensorNetwork."""
        @overload
        def chain(start_slot: pyo3_stub_gen.RustType["SpensoSlot"], end_slot: pyo3_stub_gen.RustType["SpensoSlot"], first: pyo3_stub_gen.RustType["ConvertibleToSpensoNet"], second: typing.Union[Tensor, TensorNetwork], /, *factors: pyo3_stub_gen.RustType["ConvertibleToSpensoNet"]) -> pyo3_stub_gen.RustType["SpensoNet"]:
            """Build an explicitly-ended ordered tensor chain.

            Examples
            --------
            >>> from symbolica.community import tensor as sp
            >>> r = sp.Representation.euc(2)
            >>> A = sp.TensorName("docs::A")(r, r)
            >>> chain = sp.chain(r("i"), r("j"), A, A)

            Parameters
            ----------
            start_slot : pyo3_stub_gen.RustType['SpensoSlot']
                Typed Slot labelling the first free endpoint of the chain.
            end_slot : pyo3_stub_gen.RustType['SpensoSlot']
                Typed Slot labelling the last free endpoint. Matching compatible endpoints close the chain.
            first : pyo3_stub_gen.RustType['ConvertibleToSpensoNet']
                First tensor factor in the ordered product.
            second : Union[Tensor, TensorNetwork]
                Second tensor factor in the ordered product.
            factors : pyo3_stub_gen.RustType['ConvertibleToSpensoNet']
                Ordered tensor factors. Symbolic operands produce a TensorExpression; any component tensor or network promotes the result to a TensorNetwork."""
        @overload
        def chain(start_slot: pyo3_stub_gen.RustType["SpensoSlot"], end_slot: pyo3_stub_gen.RustType["SpensoSlot"], *factors: pyo3_stub_gen.RustType["ConvertibleToSpensoNet"]) -> pyo3_stub_gen.RustType["TensorDispatch"]:
            """Build an explicitly-ended ordered tensor chain.

            Examples
            --------
            >>> from symbolica.community import tensor as sp
            >>> r = sp.Representation.euc(2)
            >>> A = sp.TensorName("docs::A")(r, r)
            >>> chain = sp.chain(r("i"), r("j"), A, A)

            Parameters
            ----------
            start_slot : pyo3_stub_gen.RustType['SpensoSlot']
                Typed Slot labelling the first free endpoint of the chain.
            end_slot : pyo3_stub_gen.RustType['SpensoSlot']
                Typed Slot labelling the last free endpoint. Matching compatible endpoints close the chain.
            factors : pyo3_stub_gen.RustType['ConvertibleToSpensoNet']
                Ordered tensor factors. Symbolic operands produce a TensorExpression; any component tensor or network promotes the result to a TensorNetwork."""
        "#,
    )
)]
#[spenso_macros::track_usage(crate::record_usage)]
#[pyfunction]
#[pyo3(signature = (start_slot, end_slot, *factors))]
fn chain(
    py: Python<'_>,
    start_slot: SpensoSlot,
    end_slot: SpensoSlot,
    factors: &Bound<'_, PyTuple>,
) -> PyResult<TensorDispatch> {
    if factors.iter().any(|factor| is_concrete_network(&factor)) {
        let mut factors = factors.iter();
        let first = factors
            .next()
            .expect("the empty factor sequence was handled above")
            .extract::<ConvertibleToSpensoNet>()?
            .to_net();
        let first_channel = first
            .structure
            .matrix_channel()
            .ok_or_else(|| PyValueError::new_err("chain factors have no unique matrix channel"))?;
        let mut value = first.chain_form(first_channel)?;
        for factor in factors {
            let factor = factor.extract::<ConvertibleToSpensoNet>()?.to_net();
            let factor_channel = factor.structure.matrix_channel().ok_or_else(|| {
                PyValueError::new_err("chain factors have no unique matrix channel")
            })?;
            value = value.compose(
                ConvertibleToSpensoNet(factor),
                (0, 1),
                (factor_channel.input, factor_channel.output),
            )?;
        }
        let slots = value.structure.structure().logical_slots();
        if start_slot.slot.rep() != slots[0].rep() || end_slot.slot.rep() != slots[1].rep() {
            return Err(PyValueError::new_err(
                "chain endpoints are incompatible with the factor channel",
            ));
        }
        let value = if start_slot.slot.aind() == end_slot.slot.aind() {
            value.trace(Some((0, 1)))
        } else {
            value.set_port_indices(&HashMap::from([
                (0, start_slot.slot.aind()),
                (1, end_slot.slot.aind()),
            ]))
        }?;
        return Py::new(py, value).map(TensorDispatch::Network);
    }

    let factors = factors
        .iter()
        .map(|factor| structured_operand(&factor, "chain()"))
        .collect::<PyResult<Vec<_>>>()?;
    SymbolicTensor::chain(start_slot.slot, end_slot.slot, &factors)
        .map_err(|error| PyValueError::new_err(error.to_string()))
        .and_then(|value| TensorExpression::from_structured(py, value))
        .map(TensorDispatch::Expression)
}

#[cfg(feature = "python_stubgen")]
symbolic_function_overload!("trace", r##"Close an ordered factor sequence into a canonical cyclic trace.

Examples
--------
>>> from symbolica.community import tensor as sp
>>> r = sp.Representation.euc(2)
>>> A = sp.TensorName("docs::A")(r, r)
>>> trace = sp.trace(r, A, A)

Parameters
----------
representation : Representation
    Representation of the matrix channel to trace over.
factors : TensorExpression | Expression
    Ordered tensor factors. Symbolic operands produce a TensorExpression; any component tensor or network promotes the result to a TensorNetwork."##, ["representation": SpensoRepresentation], ["factors": VarPositional]);

/// Close an ordered factor sequence into a canonical cyclic trace.
///
/// Examples
/// --------
/// >>> from symbolica.community import tensor as sp
/// >>> r = sp.Representation.euc(2)
/// >>> A = sp.TensorName("docs::A")(r, r)
/// >>> trace = sp.trace(r, A, A)
///
/// Parameters
/// ----------
/// representation : Representation
///     Representation of the matrix channel to trace over.
/// factors : TensorExpression | Expression
///     Ordered tensor factors. Symbolic operands produce a TensorExpression; any component tensor or network promotes the result to a TensorNetwork.
#[cfg_attr(
    feature = "python_stubgen",
    pyo3_stub_gen::derive::gen_stub_pyfunction(
        module = "symbolica.community.tensor",
        no_default_overload = true,
        python_overload = r#"
        import typing

        @overload
        def trace(representation: pyo3_stub_gen.RustType["SpensoRepresentation"], factor: typing.Union[Tensor, TensorNetwork], /, *factors: pyo3_stub_gen.RustType["ConvertibleToSpensoNet"]) -> pyo3_stub_gen.RustType["SpensoNet"]:
            """Close an ordered factor sequence into a canonical cyclic trace.

            Examples
            --------
            >>> from symbolica.community import tensor as sp
            >>> r = sp.Representation.euc(2)
            >>> A = sp.TensorName("docs::A")(r, r)
            >>> trace = sp.trace(r, A, A)

            Parameters
            ----------
            representation : pyo3_stub_gen.RustType['SpensoRepresentation']
                Representation of the matrix channel to trace over.
            factor : Union[Tensor, TensorNetwork]
                First tensor factor in the ordered product.
            factors : pyo3_stub_gen.RustType['ConvertibleToSpensoNet']
                Ordered tensor factors. Symbolic operands produce a TensorExpression; any component tensor or network promotes the result to a TensorNetwork."""
        @overload
        def trace(representation: pyo3_stub_gen.RustType["SpensoRepresentation"], first: pyo3_stub_gen.RustType["ConvertibleToSpensoNet"], second: typing.Union[Tensor, TensorNetwork], /, *factors: pyo3_stub_gen.RustType["ConvertibleToSpensoNet"]) -> pyo3_stub_gen.RustType["SpensoNet"]:
            """Close an ordered factor sequence into a canonical cyclic trace.

            Examples
            --------
            >>> from symbolica.community import tensor as sp
            >>> r = sp.Representation.euc(2)
            >>> A = sp.TensorName("docs::A")(r, r)
            >>> trace = sp.trace(r, A, A)

            Parameters
            ----------
            representation : pyo3_stub_gen.RustType['SpensoRepresentation']
                Representation of the matrix channel to trace over.
            first : pyo3_stub_gen.RustType['ConvertibleToSpensoNet']
                First tensor factor in the ordered product.
            second : Union[Tensor, TensorNetwork]
                Second tensor factor in the ordered product.
            factors : pyo3_stub_gen.RustType['ConvertibleToSpensoNet']
                Ordered tensor factors. Symbolic operands produce a TensorExpression; any component tensor or network promotes the result to a TensorNetwork."""
        @overload
        def trace(representation: pyo3_stub_gen.RustType["SpensoRepresentation"], *factors: pyo3_stub_gen.RustType["ConvertibleToSpensoNet"]) -> pyo3_stub_gen.RustType["TensorDispatch"]:
            """Close an ordered factor sequence into a canonical cyclic trace.

            Examples
            --------
            >>> from symbolica.community import tensor as sp
            >>> r = sp.Representation.euc(2)
            >>> A = sp.TensorName("docs::A")(r, r)
            >>> trace = sp.trace(r, A, A)

            Parameters
            ----------
            representation : pyo3_stub_gen.RustType['SpensoRepresentation']
                Representation of the matrix channel to trace over.
            factors : pyo3_stub_gen.RustType['ConvertibleToSpensoNet']
                Ordered tensor factors. Symbolic operands produce a TensorExpression; any component tensor or network promotes the result to a TensorNetwork."""
        "#,
    )
)]
#[spenso_macros::track_usage(crate::record_usage)]
#[pyfunction]
#[pyo3(signature = (representation, *factors))]
fn trace(
    py: Python<'_>,
    representation: SpensoRepresentation,
    factors: &Bound<'_, PyTuple>,
) -> PyResult<TensorDispatch> {
    if factors.iter().any(|factor| is_concrete_network(&factor)) {
        let mut factors = factors.iter();
        let mut value = factors
            .next()
            .expect("the empty factor sequence was handled above")
            .extract::<ConvertibleToSpensoNet>()?
            .to_net();
        for factor in factors {
            value = value.multiply_network(factor.extract::<ConvertibleToSpensoNet>()?.to_net())?;
        }
        let channel = value
            .structure
            .matrix_channel()
            .ok_or_else(|| PyValueError::new_err("trace factors have no unique matrix channel"))?;
        let input = value.structure.structure().logical_slots()[channel.input];
        if representation.representation != input.rep() {
            return Err(PyValueError::new_err(
                "trace representation does not match the factor channel",
            ));
        }
        return value
            .trace(Some((channel.input, channel.output)))
            .and_then(|network| Py::new(py, network))
            .map(TensorDispatch::Network);
    }

    let factors = factors
        .iter()
        .map(|factor| structured_operand(&factor, "trace()"))
        .collect::<PyResult<Vec<_>>>()?;
    SymbolicTensor::trace(representation.representation, &factors)
        .map_err(|error| PyValueError::new_err(error.to_string()))
        .and_then(|value| TensorExpression::from_structured(py, value))
        .map(TensorDispatch::Expression)
}

pub(crate) fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    TensorExpression::init(m)?;
    m.add_class::<AutoIndex>()?;
    m.add_function(wrap_pyfunction!(as_tensor, m)?)?;
    m.add_function(wrap_pyfunction!(dot, m)?)?;
    m.add_function(wrap_pyfunction!(chain, m)?)?;
    m.add_function(wrap_pyfunction!(trace, m)?)?;
    let auto = Py::new(m.py(), AutoIndex)?;
    m.add("AUTO", auto.clone_ref(m.py()))?;
    m.add("_", auto)?;
    Ok(())
}

// Nested RustType markers in Python unions lose their trait mapping. Register
// symbolic overloads with the same typed inventory used by __getitem__ below.
#[cfg(feature = "python_stubgen")]
macro_rules! symbolic_method_overload {
    ($name:literal, $doc:literal, $argument:literal, $input:ty $(, $keyword:literal: $kind:ty)*) => {
        submit! {
            PyMethodsInfo {
                struct_id: std::any::TypeId::of::<TensorExpression>,
                attrs: &[], getters: &[], setters: &[],
                file: file!(), line: line!(), column: column!(),
                methods: &[MethodInfo {
                    name: $name,
                    parameters: &[
                        ParameterInfo {
                            name: $argument,
                            kind: ParameterKind::PositionalOrKeyword,
                            default: ParameterDefault::None,
                            type_info: || TensorExpression::type_input() | <$input>::type_input(),
                        },
                        $(ParameterInfo {
                            name: $keyword,
                            kind: ParameterKind::KeywordOnly,
                            default: ParameterDefault::None,
                            type_info: <$kind>::type_input,
                        },)*
                    ],
                    r#type: MethodType::Instance,
                    r#return: TensorExpression::type_output,
                    doc: $doc, is_async: false, deprecated: None, type_ignored: None,
                    is_overload: true,
                }],
            }
        }
    };
}

#[cfg(feature = "python_stubgen")]
symbolic_method_overload!(
    "__add__",
    r##"Add compatible tensors, preserving their logical external axes.

Examples
--------
>>> from symbolica.community import tensor as sp
>>> r = sp.Representation.euc(2)
>>> A = sp.TensorName("docs::A")(r, r)
>>> result = A + A

Parameters
----------
rhs : TensorExpression | Expression or scalar
    Right operand, with a compatible tensor interface for addition or subtraction. A denominator must be scalar."##,
    "rhs",
    ConvertibleToExpression
);
#[cfg(feature = "python_stubgen")]
symbolic_method_overload!(
    "__radd__",
    r##"Add compatible tensors, preserving their logical external axes.

Examples
--------
>>> from symbolica.community import tensor as sp
>>> r = sp.Representation.euc(2)
>>> A = sp.TensorName("docs::A")(r, r)
>>> result = A + A

Parameters
----------
lhs : TensorExpression | Expression or scalar
    Left operand. Addition and subtraction require compatible tensor interfaces."##,
    "lhs",
    ConvertibleToExpression
);
#[cfg(feature = "python_stubgen")]
symbolic_method_overload!(
    "__sub__",
    r##"Subtract compatible tensors, preserving their logical external axes.

Examples
--------
>>> from symbolica.community import tensor as sp
>>> r = sp.Representation.euc(2)
>>> A = sp.TensorName("docs::A")(r, r)
>>> result = A - A

Parameters
----------
rhs : TensorExpression | Expression or scalar
    Right operand, with a compatible tensor interface for addition or subtraction. A denominator must be scalar."##,
    "rhs",
    ConvertibleToExpression
);
#[cfg(feature = "python_stubgen")]
symbolic_method_overload!(
    "__rsub__",
    r##"Subtract compatible tensors, preserving their logical external axes.

Examples
--------
>>> from symbolica.community import tensor as sp
>>> r = sp.Representation.euc(2)
>>> A = sp.TensorName("docs::A")(r, r)
>>> result = A - A

Parameters
----------
lhs : TensorExpression | Expression or scalar
    Left operand. Addition and subtraction require compatible tensor interfaces."##,
    "lhs",
    ConvertibleToExpression
);
#[cfg(feature = "python_stubgen")]
symbolic_method_overload!(
    "__mul__",
    r##"Multiply tensor operands, contracting repeated compatible explicit indices.

Examples
--------
>>> from symbolica.community import tensor as sp
>>> r = sp.Representation.euc(2)
>>> A = sp.TensorName("docs::A")(r, r)
>>> result = A * 2

Parameters
----------
rhs : TensorExpression | Expression or scalar
    Right operand, with a compatible tensor interface for addition or subtraction. A denominator must be scalar."##,
    "rhs",
    ConvertibleToExpression
);
#[cfg(feature = "python_stubgen")]
symbolic_method_overload!(
    "__rmul__",
    r##"Multiply tensor operands, contracting repeated compatible explicit indices.

Examples
--------
>>> from symbolica.community import tensor as sp
>>> r = sp.Representation.euc(2)
>>> A = sp.TensorName("docs::A")(r, r)
>>> result = 2 * A

Parameters
----------
lhs : TensorExpression | Expression or scalar
    Left operand. Addition and subtraction require compatible tensor interfaces."##,
    "lhs",
    ConvertibleToExpression
);
#[cfg(feature = "python_stubgen")]
symbolic_method_overload!(
    "__truediv__",
    r##"Divide tensor coefficients by a scalar denominator.

Examples
--------
>>> from symbolica.community import tensor as sp
>>> r = sp.Representation.euc(2)
>>> A = sp.TensorName("docs::A")(r, r)
>>> result = A / 2

Parameters
----------
rhs : TensorExpression | Expression or scalar
    Right operand, with a compatible tensor interface for addition or subtraction. A denominator must be scalar."##,
    "rhs",
    ConvertibleToExpression
);
#[cfg(feature = "python_stubgen")]
symbolic_method_overload!(
    "__rtruediv__",
    r##"Divide a scalar by a scalar tensor expression.

Examples
--------
>>> from symbolica.community import tensor as sp
>>> scalar = sp.TensorExpression(2)
>>> result = 1 / scalar

Parameters
----------
lhs : TensorExpression | Expression or scalar
    Left operand. Addition and subtraction require compatible tensor interfaces."##,
    "lhs",
    ConvertibleToExpression
);
#[cfg(feature = "python_stubgen")]
symbolic_method_overload!(
    "outer",
    r##"Form an outer product, keeping both operands’ ports free even when their labels coincide.

Examples
--------
>>> from symbolica.community import tensor as sp
>>> r = sp.Representation.euc(2)
>>> A = sp.TensorName("docs::A")(r, r)
>>> result = A.outer(A)

Parameters
----------
rhs : TensorExpression | Expression
    Other tensor, network, or scalar operand. Concrete component data or a network makes the result a TensorNetwork."##,
    "rhs",
    PythonExpression
);
#[cfg(feature = "python_stubgen")]
symbolic_method_overload!("contract_ports", r##"Contract one chosen pair of axes and retain the other external ports in logical order.

Examples
--------
>>> from symbolica.community import tensor as sp
>>> r = sp.Representation.euc(2)
>>> A = sp.TensorName("docs::A")(r, r)
>>> result = A.contract_ports(A, left=1, right=0)

Parameters
----------
rhs : TensorExpression | Expression
    Other tensor, network, or scalar operand. Concrete component data or a network makes the result a TensorNetwork.
left : int
    Zero-based axis position in this operand to contract.
right : int
    Zero-based axis position in rhs to contract with the selected left axis."##, "rhs", PythonExpression, "left": usize, "right": usize);
#[cfg(feature = "python_stubgen")]
symbolic_method_overload!("compose", r##"Compose two explicitly selected matrix channels, retaining the uncontracted endpoints.

Examples
--------
>>> from symbolica.community import tensor as sp
>>> r = sp.Representation.euc(2)
>>> A = sp.TensorName("docs::A")(r, r)
>>> result = A.compose(A, left=(0, 1), right=(0, 1))

Parameters
----------
rhs : TensorExpression | Expression
    Other tensor, network, or scalar operand. Concrete component data or a network makes the result a TensorNetwork.
left : tuple[int, int]
    (input_axis, output_axis) channel in this operand.
right : tuple[int, int]
    (input_axis, output_axis) channel in rhs. The selected connecting ports must have compatible representations."##, "rhs", PythonExpression, "left": (usize, usize), "right": (usize, usize));

// The stub generator's Python parser requires typing.Union instead of `|`.
#[cfg(feature = "python_stubgen")]
submit! {
    pyo3_stub_gen::derive::gen_methods_from_python! {
        r#"
        import typing

        class TensorExpression:
            @overload
            def __add__(self, rhs: typing.Union[Tensor, TensorNetwork]) -> pyo3_stub_gen.RustType["SpensoNet"]:
                """Add compatible tensors, preserving their logical external axes.

                Examples
                --------
                >>> from symbolica.community import tensor as sp
                >>> r = sp.Representation.euc(2)
                >>> A = sp.TensorName("docs::A")(r, r)
                >>> result = A + A

                Parameters
                ----------
                rhs : Union[Tensor, TensorNetwork]
                    Right operand, with a compatible tensor interface for addition or subtraction. A denominator must be scalar."""
            @overload
            def __radd__(self, lhs: typing.Union[Tensor, TensorNetwork]) -> pyo3_stub_gen.RustType["SpensoNet"]:
                """Add compatible tensors, preserving their logical external axes.

                Examples
                --------
                >>> from symbolica.community import tensor as sp
                >>> r = sp.Representation.euc(2)
                >>> A = sp.TensorName("docs::A")(r, r)
                >>> result = A + A

                Parameters
                ----------
                lhs : Union[Tensor, TensorNetwork]
                    Left operand. Addition and subtraction require compatible tensor interfaces."""
            @overload
            def __sub__(self, rhs: typing.Union[Tensor, TensorNetwork]) -> pyo3_stub_gen.RustType["SpensoNet"]:
                """Subtract compatible tensors, preserving their logical external axes.

                Examples
                --------
                >>> from symbolica.community import tensor as sp
                >>> r = sp.Representation.euc(2)
                >>> A = sp.TensorName("docs::A")(r, r)
                >>> result = A - A

                Parameters
                ----------
                rhs : Union[Tensor, TensorNetwork]
                    Right operand, with a compatible tensor interface for addition or subtraction. A denominator must be scalar."""
            @overload
            def __rsub__(self, lhs: typing.Union[Tensor, TensorNetwork]) -> pyo3_stub_gen.RustType["SpensoNet"]:
                """Subtract compatible tensors, preserving their logical external axes.

                Examples
                --------
                >>> from symbolica.community import tensor as sp
                >>> r = sp.Representation.euc(2)
                >>> A = sp.TensorName("docs::A")(r, r)
                >>> result = A - A

                Parameters
                ----------
                lhs : Union[Tensor, TensorNetwork]
                    Left operand. Addition and subtraction require compatible tensor interfaces."""
            @overload
            def __mul__(self, rhs: typing.Union[Tensor, TensorNetwork]) -> pyo3_stub_gen.RustType["SpensoNet"]:
                """Multiply tensor operands, contracting repeated compatible explicit indices.

                Examples
                --------
                >>> from symbolica.community import tensor as sp
                >>> r = sp.Representation.euc(2)
                >>> A = sp.TensorName("docs::A")(r, r)
                >>> result = A * 2

                Parameters
                ----------
                rhs : Union[Tensor, TensorNetwork]
                    Right operand, with a compatible tensor interface for addition or subtraction. A denominator must be scalar."""
            @overload
            def __rmul__(self, lhs: typing.Union[Tensor, TensorNetwork]) -> pyo3_stub_gen.RustType["SpensoNet"]:
                """Multiply tensor operands, contracting repeated compatible explicit indices.

                Examples
                --------
                >>> from symbolica.community import tensor as sp
                >>> r = sp.Representation.euc(2)
                >>> A = sp.TensorName("docs::A")(r, r)
                >>> result = 2 * A

                Parameters
                ----------
                lhs : Union[Tensor, TensorNetwork]
                    Left operand. Addition and subtraction require compatible tensor interfaces."""
            @overload
            def __truediv__(self, rhs: typing.Union[Tensor, TensorNetwork]) -> pyo3_stub_gen.RustType["SpensoNet"]:
                """Divide tensor coefficients by a scalar denominator.

                Examples
                --------
                >>> from symbolica.community import tensor as sp
                >>> r = sp.Representation.euc(2)
                >>> A = sp.TensorName("docs::A")(r, r)
                >>> result = A / 2

                Parameters
                ----------
                rhs : Union[Tensor, TensorNetwork]
                    Right operand, with a compatible tensor interface for addition or subtraction. A denominator must be scalar."""
            @overload
            def __rtruediv__(self, lhs: typing.Union[Tensor, TensorNetwork]) -> pyo3_stub_gen.RustType["SpensoNet"]:
                """Divide a scalar by a scalar tensor expression.

                Examples
                --------
                >>> from symbolica.community import tensor as sp
                >>> scalar = sp.TensorExpression(2)
                >>> result = 1 / scalar

                Parameters
                ----------
                lhs : Union[Tensor, TensorNetwork]
                    Left operand. Addition and subtraction require compatible tensor interfaces."""
            @overload
            def outer(self, rhs: typing.Union[Tensor, TensorNetwork]) -> pyo3_stub_gen.RustType["SpensoNet"]:
                """Form an outer product, keeping both operands’ ports free even when their labels coincide.

                Examples
                --------
                >>> from symbolica.community import tensor as sp
                >>> r = sp.Representation.euc(2)
                >>> A = sp.TensorName("docs::A")(r, r)
                >>> result = A.outer(A)

                Parameters
                ----------
                rhs : Union[Tensor, TensorNetwork]
                    Other tensor, network, or scalar operand. Concrete component data or a network makes the result a TensorNetwork."""
            @overload
            def contract_ports(self, rhs: typing.Union[Tensor, TensorNetwork], *, left: int, right: int) -> pyo3_stub_gen.RustType["SpensoNet"]:
                """Contract one chosen pair of axes and retain the other external ports in logical order.

                Examples
                --------
                >>> from symbolica.community import tensor as sp
                >>> r = sp.Representation.euc(2)
                >>> A = sp.TensorName("docs::A")(r, r)
                >>> result = A.contract_ports(A, left=1, right=0)

                Parameters
                ----------
                rhs : Union[Tensor, TensorNetwork]
                    Other tensor, network, or scalar operand. Concrete component data or a network makes the result a TensorNetwork.
                left : int
                    Zero-based axis position in this operand to contract.
                right : int
                    Zero-based axis position in rhs to contract with the selected left axis."""
            @overload
            def compose(self, rhs: typing.Union[Tensor, TensorNetwork], *, left: tuple[int, int], right: tuple[int, int]) -> pyo3_stub_gen.RustType["SpensoNet"]:
                """Compose two explicitly selected matrix channels, retaining the uncontracted endpoints.

                Examples
                --------
                >>> from symbolica.community import tensor as sp
                >>> r = sp.Representation.euc(2)
                >>> A = sp.TensorName("docs::A")(r, r)
                >>> result = A.compose(A, left=(0, 1), right=(0, 1))

                Parameters
                ----------
                rhs : Union[Tensor, TensorNetwork]
                    Other tensor, network, or scalar operand. Concrete component data or a network makes the result a TensorNetwork.
                left : tuple[int, int]
                    (input_axis, output_axis) channel in this operand.
                right : tuple[int, int]
                    (input_axis, output_axis) channel in rhs. The selected connecting ports must have compatible representations."""
        "#
    }
}

#[cfg(feature = "python_stubgen")]
submit! {
    PyMethodsInfo {
        struct_id: std::any::TypeId::of::<TensorExpression>,
        attrs: &[],
        getters: &[],
        setters: &[],
        file: file!(),
        line: line!(),
        column: column!(),
        methods: &[
            MethodInfo {
                name: "__getitem__",
                parameters: &[ParameterInfo {
                    name: "item",
                    kind: ParameterKind::PositionalOrKeyword,
                    default: ParameterDefault::None,
                    type_info: usize::type_input,
                }],
                r#type: MethodType::Instance,
                r#return: Vec::<usize>::type_output,
                doc: r##"Convert between flat positions and tensor coordinates in logical row-major order.
An integer returns coordinates, a coordinate list returns a flat position, and a slice
returns coordinate lists. This indexes the layout; use `.components()` for values.

Examples
--------
>>> from symbolica.community import tensor as sp
>>> r = sp.Representation.euc(2)
>>> A = sp.TensorName("docs::A")(r, r)
>>> coordinates = A[3]
>>> flat_index = A[[1, 1]]
>>> all_coordinates = A[:]

Parameters
----------
item : int
    Flat integer position, list of coordinates, or slice of flat positions. Coordinates follow .structure.slots and use row-major order."##,
                is_async: false,
                deprecated: None,
                type_ignored: Some(pyo3_stub_gen::type_info::IgnoreTarget::Specified(&["override"])),
                is_overload: true,
            },
            MethodInfo {
                name: "__getitem__",
                parameters: &[ParameterInfo {
                    name: "item",
                    kind: ParameterKind::PositionalOrKeyword,
                    default: ParameterDefault::None,
                    type_info: Vec::<usize>::type_input,
                }],
                r#type: MethodType::Instance,
                r#return: usize::type_output,
                doc: r##"Convert between flat positions and tensor coordinates in logical row-major order.
An integer returns coordinates, a coordinate list returns a flat position, and a slice
returns coordinate lists. This indexes the layout; use `.components()` for values.

Examples
--------
>>> from symbolica.community import tensor as sp
>>> r = sp.Representation.euc(2)
>>> A = sp.TensorName("docs::A")(r, r)
>>> coordinates = A[3]
>>> flat_index = A[[1, 1]]
>>> all_coordinates = A[:]

Parameters
----------
item : int
    Flat integer position, list of coordinates, or slice of flat positions. Coordinates follow .structure.slots and use row-major order."##,
                is_async: false,
                deprecated: None,
                type_ignored: Some(pyo3_stub_gen::type_info::IgnoreTarget::Specified(&["override"])),
                is_overload: true,
            },
            MethodInfo {
                name: "__getitem__",
                parameters: &[ParameterInfo {
                    name: "item",
                    kind: ParameterKind::PositionalOrKeyword,
                    default: ParameterDefault::None,
                    type_info: || TypeInfo::builtin("slice"),
                }],
                r#type: MethodType::Instance,
                r#return: Vec::<Vec<usize>>::type_output,
                doc: r##"Convert between flat positions and tensor coordinates in logical row-major order.
An integer returns coordinates, a coordinate list returns a flat position, and a slice
returns coordinate lists. This indexes the layout; use `.components()` for values.

Examples
--------
>>> from symbolica.community import tensor as sp
>>> r = sp.Representation.euc(2)
>>> A = sp.TensorName("docs::A")(r, r)
>>> coordinates = A[3]
>>> flat_index = A[[1, 1]]
>>> all_coordinates = A[:]

Parameters
----------
item : int
    Flat integer position, list of coordinates, or slice of flat positions. Coordinates follow .structure.slots and use row-major order."##,
                is_async: false,
                deprecated: None,
                type_ignored: Some(pyo3_stub_gen::type_info::IgnoreTarget::Specified(&["override", "invalid-method-override"])),
                is_overload: true,
            },
        ],
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use idenso::{
        IndexTooling, IndexToolingError,
        representations::{Bispinor, ColorAdjoint, ColorFundamental},
    };
    use spenso::structure::{
        dimension::Dimension,
        representation::{ExtendibleReps, Minkowski, RepName},
    };
    use symbolica::{
        atom::{AtomCore, FunctionBuilder},
        symbol,
    };

    fn infer_interface(atom: &Atom) -> PyResult<PartialStructure> {
        InterfaceInference::default()
            .infer_validated(atom.as_view())
            .map_err(TensorExpression::inference_error)
    }

    fn explicit_interface(
        representation: Representation<LibraryRep>,
        index: AbstractIndex,
    ) -> PartialStructure {
        PartialStructure::from_logical_slots([representation.slot(PartialIndex::Explicit(index))])
    }

    #[cfg(feature = "python_stubgen")]
    #[test]
    fn dispatch_stubs_preserve_symbolic_results_and_concrete_promotion() {
        let methods = pyo3_stub_gen::inventory::iter::<PyMethodsInfo>
            .into_iter()
            .filter(|info| (info.struct_id)() == std::any::TypeId::of::<TensorExpression>())
            .flat_map(|info| info.methods);
        for name in [
            "__add__",
            "__radd__",
            "__sub__",
            "__rsub__",
            "__mul__",
            "__rmul__",
            "__truediv__",
            "__rtruediv__",
            "outer",
            "contract_ports",
            "compose",
        ] {
            let overloads = methods
                .clone()
                .filter(|method| method.name == name)
                .collect::<Vec<_>>();
            assert_eq!(
                overloads.len(),
                2,
                "{name} must not retain the Any overload"
            );
            assert!(overloads.iter().all(|method| method.is_overload));
            let symbolic = if name.starts_with("__") {
                TensorExpression::type_input() | ConvertibleToExpression::type_input()
            } else {
                TensorExpression::type_input() | PythonExpression::type_input()
            };
            assert!(overloads.iter().any(|method| {
                (method.parameters[0].type_info)() == symbolic
                    && (method.r#return)() == TensorExpression::type_output()
            }));
            assert!(overloads.iter().any(|method| {
                (method.parameters[0].type_info)().to_string()
                    == "typing.Union[Tensor, TensorNetwork]"
                    && (method.r#return)() == SpensoNet::type_output()
            }));
            for method in overloads {
                assert!(
                    method.parameters[1..]
                        .iter()
                        .all(|parameter| { matches!(parameter.kind, ParameterKind::KeywordOnly) })
                );
            }
        }

        let functions = pyo3_stub_gen::inventory::iter::<pyo3_stub_gen::type_info::PyFunctionInfo>
            .into_iter()
            .filter(|info| info.module == Some("symbolica.community.tensor"));
        for name in ["chain", "trace"] {
            let overloads = functions
                .clone()
                .filter(|function| function.name == name)
                .collect::<Vec<_>>();
            assert_eq!(overloads.len(), 4);
            // Concrete prefixes let type checkers resolve common mixed calls without
            // relying on overlapping variadic-tuple inference.
            for concrete_parameter in ["factor", "second"] {
                assert!(overloads.iter().any(|function| {
                    function.parameters.iter().any(|parameter| {
                        parameter.name == concrete_parameter
                            && matches!(parameter.kind, ParameterKind::PositionalOnly)
                            && (parameter.type_info)().to_string()
                                == "typing.Union[Tensor, TensorNetwork]"
                    }) && (function.r#return)() == SpensoNet::type_output()
                }));
            }
            assert!(overloads.iter().all(|function| function.is_overload));
            assert!(overloads.iter().any(|function| {
                let factors = function.parameters.last().unwrap();
                matches!(factors.kind, ParameterKind::VarPositional)
                    && (factors.type_info)()
                        == (TensorExpression::type_input() | PythonExpression::type_input())
                    && (function.r#return)() == TensorExpression::type_output()
            }));
            // A possibly empty, mixed sequence cannot promise concrete promotion.
            assert!(overloads.iter().any(|function| {
                (function.r#return)() == TensorDispatch::type_output()
                    && !(function.parameters.last().unwrap().type_info)()
                        .to_string()
                        .contains("Any")
            }));
        }
    }

    #[test]
    fn default_color_simplification_preserves_external_tensor_ports() {
        idenso::representations::initialize();
        Python::initialize();
        let mink = Minkowski {}.new_rep(Dimension::from(symbol!("color_ports_D")));
        let adjoint = ColorAdjoint {}.new_rep(8);
        let fundamental = ColorFundamental {}.new_rep(3);
        let [mu, nu] = [1, 2].map(|index| {
            mink.slot::<AbstractIndex, _>(AbstractIndex::Normal(index))
                .to_atom()
        });
        let [a, b, c, d] = [3, 4, 5, 6].map(|index| {
            adjoint
                .slot::<AbstractIndex, _>(AbstractIndex::Normal(index))
                .to_atom()
        });
        let color_factors = [
            idenso::color_f!(&a, &c, &d) * idenso::color_f!(&b, &c, &d),
            SPENSO_TAG.trace(
                fundamental.to_symbolic([]),
                [idenso::color_t!(&a), idenso::color_t!(&b)],
            ),
        ];

        for color in color_factors {
            let numerator = ETS.metric(&mu, &nu) * color;
            let expected = infer_interface(&numerator).unwrap();
            assert_eq!(expected.canonical().order(), 4);

            // Casimir/index representation arguments are scalar metadata, not ports.
            let simplified = SymbolicTensor::infer(numerator)
                .unwrap()
                .reduce_algebra(&idenso::tensor::simplification::AlgebraSettings {
                    color: Some(idenso::color::ColorSimplifySettings::default()),
                    ..Default::default()
                })
                .unwrap()
                .resolved()
                .unwrap()
                .into_expression();
            let actual = infer_interface(&simplified).unwrap();
            assert_eq!(actual.canonical(), expected.canonical(), "{simplified}");
        }
    }

    #[test]
    fn compact_epsilon_preserves_only_uncontracted_ports() {
        idenso::representations::initialize();
        Python::initialize();
        let mink = Minkowski {}.new_rep(4);
        let [mu, nu] = [1, 2].map(|index| {
            mink.slot::<AbstractIndex, _>(AbstractIndex::Normal(index))
                .to_atom()
        });
        let [p, q] = [
            spenso::vector_symbol!("epsilon_validation_p"),
            spenso::vector_symbol!("epsilon_validation_q"),
        ]
        .map(|symbol| {
            FunctionBuilder::new(symbol)
                .add_arg(mink.to_symbolic([]))
                .finish()
        });
        let epsilon = idenso::epsilon::epsilon4(&mu, &nu, &p, &q);
        let interface = infer_interface(&epsilon).unwrap();
        assert_eq!(interface.canonical().order(), 2);
        let scalar = epsilon.clone() * epsilon;
        assert!(infer_interface(&scalar).unwrap().canonical().is_scalar());
        assert!(
            infer_interface(
                SymbolicTensor::infer(scalar)
                    .unwrap()
                    .reduce_algebra(&idenso::tensor::AlgebraSettings {
                        epsilon: true,
                        ..Default::default()
                    })
                    .unwrap()
                    .resolved()
                    .unwrap()
                    .expression()
            )
            .unwrap()
            .canonical()
            .is_scalar()
        );
        assert!(infer_interface(&idenso::epsilon::epsilon4(&mu, &nu, &p, Atom::num(1))).is_err());
    }

    #[test]
    fn explicit_dyad_projection_is_independent_of_metric_and_operand_order() {
        idenso::representations::initialize();
        Python::initialize();
        let p = spenso::vector_symbol!("dyad_projection_p");
        let q = spenso::vector_symbol!("dyad_projection_q");
        for dimension in [
            Dimension::Concrete(4),
            Dimension::from(symbol!("dyad_projection_D")),
        ] {
            let mink = Minkowski {}.new_rep(dimension);
            let [mu, nu] = [201, 202].map(|index| {
                mink.slot::<AbstractIndex, _>(AbstractIndex::Normal(index))
                    .to_atom()
            });
            let vector = |name, slot: &Atom| FunctionBuilder::new(name).add_arg(slot).finish();
            let dyad = vector(p, &mu) * vector(q, &nu);
            let expected_dot = FunctionBuilder::new(SPENSO_TAG.dot)
                .add_arg(vector(p, &mink.to_symbolic([])))
                .add_arg(vector(q, &mink.to_symbolic([])))
                .finish();
            for (input, output) in [(&mu, &nu), (&nu, &mu)] {
                let metric = FunctionBuilder::new(ETS.metric)
                    .add_arg(input)
                    .add_arg(output)
                    .finish();
                for (numerator, expected) in [
                    (dyad.clone(), expected_dot.clone()),
                    (
                        dyad.clone() + &metric,
                        expected_dot.clone() + dimension.to_symbolic(),
                    ),
                ] {
                    let numerator = SymbolicTensor::new(
                        numerator.clone(),
                        infer_interface(&numerator).unwrap(),
                    );
                    let metric =
                        SymbolicTensor::new(metric.clone(), infer_interface(&metric).unwrap());
                    for (left, right) in [(&numerator, &metric), (&metric, &numerator)] {
                        let product = left.multiply(right).unwrap();
                        assert!(product.is_scalar());
                        assert!(
                            infer_interface(product.expression())
                                .unwrap()
                                .canonical()
                                .is_scalar()
                        );
                        assert_eq!(
                            product
                                .contract(Default::default())
                                .unwrap()
                                .expanded()
                                .unwrap()
                                .to_dots()
                                .unwrap()
                                .into_expression(),
                            expected
                        );
                    }
                }
            }
        }
    }

    #[test]
    fn explicit_dyad_contraction_preserves_uncontracted_slots() {
        idenso::representations::initialize();
        Python::initialize();
        let mink = Minkowski {}.new_rep(4);
        let [mu, nu, rho] = [211, 212, 213].map(|index| {
            mink.slot::<AbstractIndex, _>(AbstractIndex::Normal(index))
                .to_atom()
        });
        let p = spenso::vector_symbol!("dyad_partial_p");
        let q = spenso::vector_symbol!("dyad_partial_q");
        let vector = |name, slot: &Atom| FunctionBuilder::new(name).add_arg(slot).finish();
        let numerator = vector(p, &mu) * vector(q, &nu);
        let numerator =
            SymbolicTensor::new(numerator.clone(), infer_interface(&numerator).unwrap());
        let expected = vector(p, &mu) * vector(q, &rho);
        for (input, output) in [(&nu, &rho), (&rho, &nu)] {
            let metric = FunctionBuilder::new(ETS.metric)
                .add_arg(input)
                .add_arg(output)
                .finish();
            let metric = SymbolicTensor::new(metric.clone(), infer_interface(&metric).unwrap());
            for (left, right) in [(&numerator, &metric), (&metric, &numerator)] {
                let product = left.multiply(right).unwrap();
                assert_eq!(product.rank(), 2);
                assert_eq!(
                    infer_interface(product.expression()).unwrap().canonical(),
                    infer_interface(&expected).unwrap().canonical()
                );
                assert_eq!(
                    product
                        .contract(Default::default())
                        .unwrap()
                        .resolved()
                        .unwrap()
                        .into_expression(),
                    expected
                );
            }
        }
    }

    #[test]
    fn raw_gamma_trace_accepts_longitudinal_projection() {
        idenso::representations::initialize();
        Python::initialize();
        let mink = Minkowski {}.new_rep(Dimension::from(symbol!("gamma_projection_D")));
        let bis = Bispinor {}.new_rep(4);
        let [mu, nu] = [1, 2].map(|index| {
            mink.slot::<AbstractIndex, _>(AbstractIndex::Normal(index))
                .to_atom()
        });
        let momentum = spenso::vector_symbol!("gamma_projection_p");
        let numerator = SPENSO_TAG.trace(
            bis.to_symbolic([]),
            [idenso::gamma!(&mu), idenso::gamma!(&nu)],
        );
        let projector = FunctionBuilder::new(momentum).add_arg(&mu).finish()
            * FunctionBuilder::new(momentum).add_arg(&nu).finish();
        let numerator =
            SymbolicTensor::new(numerator.clone(), infer_interface(&numerator).unwrap());
        let projector =
            SymbolicTensor::new(projector.clone(), infer_interface(&projector).unwrap());
        let contracted = numerator
            .contract_ports(&projector, &[composition::PortPair { left: 0, right: 0 }])
            .unwrap();
        assert!(contracted.is_scalar());

        let projected = contracted
            .contract(Default::default())
            .unwrap()
            .resolved()
            .unwrap()
            .into_expression();
        // Network contraction produces valid slashed gamma factors inside the trace.
        let interface = infer_interface(&projected)
            .expect("projecting a valid gamma trace must retain a valid tensor expression");
        assert!(interface.canonical().is_scalar());
    }

    #[test]
    fn compact_gamma_collection_preserves_tensor_ports_and_traces() {
        idenso::representations::initialize();
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let bis = Bispinor {}.new_rep(4);
            let [a, b, c] = [101, 102, 103].map(|index| {
                bis.slot::<AbstractIndex, _>(AbstractIndex::Normal(index)).to_atom()
            });
            let momentum = spenso::vector_symbol!("compact_gamma_momentum");
            for dimension in [Dimension::Concrete(4), Dimension::from(symbol!("compact_gamma_D"))] {
                let mink = Minkowski {}.new_rep(dimension);
                let mu = mink.slot::<AbstractIndex, _>(AbstractIndex::Normal(104)).to_atom();
                let nu = mink.slot::<AbstractIndex, _>(AbstractIndex::Normal(105)).to_atom();
                let vector = FunctionBuilder::new(momentum).add_arg(&mu).finish();
                for endpoint in [&a, &c] {
                    let original = idenso::gamma!(&a, &b, &mu)
                        * vector.clone()
                        * idenso::gamma!(&b, endpoint.as_view(), &nu);
                    let tensor = TensorExpression::from_atom_interface(py, original.clone(), None)?;
                    let compact = TensorExpression::from_structured(py, tensor.borrow(py).value.contract(Default::default()).map_err(TensorExpression::inference_error)?.resolved().map_err(TensorExpression::inference_error)?)?;
                    let settings = crate::simplification::PyGammaSimplifySettings::from_rust(idenso::dirac::GammaSimplifySettings {
                        output: idenso::dirac::GammaOutput::Chains,
                        ..Default::default()
                    });
                    let collected = compact.borrow(py).reduce_algebra(Some(&crate::simplification::PyAlgebraSettings::new(Some(settings), None, false, 16)))?;
                    let collected = TensorExpression::from_structured(py, collected.value.resolved().map_err(TensorExpression::inference_error)?)?;
                    let collected = collected.borrow(py);
                    let expected = infer_interface(&original)?;
                    assert_eq!(collected.interface().canonical(), expected.canonical());
                    let atom = &collected.atom();
                    let head = if endpoint == &a { SPENSO_TAG.trace } else { SPENSO_TAG.chain };
                    assert!(matches!(atom.as_view(), AtomView::Fun(function) if function.get_symbol() == head));
                    let slash = idenso::gamma!(
                        FunctionBuilder::new(momentum).add_arg(mink.to_symbolic([])).finish()
                    );
                    let mut contains_slash = false;
                    atom.visitor(&mut |value| {
                        contains_slash |= value == slash.as_view();
                        true
                    });
                    assert!(contains_slash);
                }
            }
            Ok(())
        }).unwrap();
    }

    #[test]
    fn compact_gamma_arguments_retain_signature_validation() {
        idenso::representations::initialize();
        Python::initialize();
        let bis = Bispinor {}.new_rep(4);
        let wrong_bis = Bispinor {}.new_rep(2);
        let mink = Minkowski {}.new_rep(4);
        let wrong_rep = ExtendibleReps::EUCLIDEAN.new_rep(Dimension::Concrete(4));
        let momentum = spenso::vector_symbol!("compact_gamma_validation_p");
        let compact = FunctionBuilder::new(momentum)
            .add_arg(mink.to_symbolic([]))
            .finish();
        let invalid = [
            FunctionBuilder::new(AGS.gamma)
                .add_arg(wrong_bis.to_symbolic([]))
                .add_arg(bis.to_symbolic([]))
                .add_arg(&compact)
                .finish(),
            FunctionBuilder::new(AGS.gamma)
                .add_arg(bis.to_symbolic([]))
                .add_arg(bis.to_symbolic([]))
                .add_arg(
                    FunctionBuilder::new(momentum)
                        .add_arg(wrong_rep.to_symbolic([]))
                        .finish(),
                )
                .finish(),
            FunctionBuilder::new(AGS.gamma)
                .add_arg(bis.to_symbolic([]))
                .add_arg(bis.to_symbolic([]))
                .add_arg(
                    FunctionBuilder::new(momentum)
                        .add_arg(
                            mink.slot::<AbstractIndex, _>(AbstractIndex::Normal(106))
                                .to_atom(),
                        )
                        .finish(),
                )
                .finish(),
            FunctionBuilder::new(AGS.gamma)
                .add_arg(bis.to_symbolic([]))
                .add_arg(bis.to_symbolic([]))
                .add_arg(ETS.metric(mink.to_symbolic([]), mink.to_symbolic([])))
                .finish(),
            FunctionBuilder::new(AGS.gamma)
                .add_arg(bis.to_symbolic([]))
                .add_arg(bis.to_symbolic([]))
                .add_arg(Atom::num(1))
                .finish(),
        ];
        for atom in invalid {
            assert!(
                infer_interface(&atom).is_err(),
                "invalid compact gamma accepted: {atom}"
            );
        }
        // A trace supplies its formal spin channel; this is distinct from
        // the fixed standalone gamma factory signature checked above.
        let formal_trace = SPENSO_TAG.trace(wrong_bis.to_symbolic([]), [idenso::gamma!(&compact)]);
        assert!(
            infer_interface(&formal_trace)
                .unwrap()
                .canonical()
                .is_scalar()
        );
    }

    #[test]
    fn tensor_expression_public_surface_has_runtime_docstrings() {
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let expression_type = py.get_type::<TensorExpression>();
            for name in [
                "g",
                "flat",
                "gamma",
                "gamma5",
                "projm",
                "projp",
                "sigma",
                "f",
                "t",
                "rank",
                "is_scalar",
                "structure",
                "name",
                "with_name",
                "to_expression",
                "expand",
                "reduce_algebra",
                "wrap_indices",
                "list_dangling",
                "canonize",
                "dirac_adjoint",
                "to_dots",
                "undo_dots",
                "to_network",
                "index",
                "outer",
                "contract",
                "compose",
                "trace",
                "format_tensor",
                "to_typst",
                "formatted",
            ] {
                let documentation = expression_type
                    .getattr(name)?
                    .getattr("__doc__")?
                    .extract::<Option<String>>()?;
                assert!(
                    documentation.is_some_and(|documentation| !documentation.trim().is_empty()),
                    "TensorExpression.{name} is missing its Python docstring"
                );
            }
            Ok(())
        })
        .unwrap();
    }

    #[test]
    fn builtin_factories_define_logical_interfaces_without_key_arguments() {
        use idenso::representations::ColorAntiFundamental;

        idenso::representations::initialize();
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let expression_type = py.get_type::<TensorExpression>();
            let euc = ExtendibleReps::EUCLIDEAN.new_rep(Dimension::Concrete(4));
            let euc_object = Py::new(
                py,
                SpensoRepresentation {
                    representation: euc,
                },
            )?;
            let adjoint: Representation<LibraryRep> =
                ColorAdjoint {}.new_rep(Dimension::Concrete(8)).cast();
            let fundamental: Representation<LibraryRep> =
                ColorFundamental {}.new_rep(Dimension::Concrete(3)).cast();
            let antifundamental: Representation<LibraryRep> = ColorAntiFundamental {}
                .new_rep(Dimension::Concrete(3))
                .cast();
            let minkowski: Representation<LibraryRep> =
                Minkowski {}.new_rep(Dimension::Concrete(4)).cast();
            let bispinor: Representation<LibraryRep> =
                Bispinor {}.new_rep(Dimension::Concrete(4)).cast();

            let cases = [
                (
                    expression_type.call_method1("g", (euc_object.clone_ref(py),))?,
                    vec![euc, euc],
                ),
                (
                    expression_type.call_method1("flat", (euc_object,))?,
                    vec![euc, euc],
                ),
                (
                    expression_type.call_method1("gamma", (4,))?,
                    vec![bispinor, bispinor, minkowski],
                ),
                (
                    expression_type.call_method1("gamma5", (4,))?,
                    vec![bispinor, bispinor],
                ),
                (
                    expression_type.call_method1("projm", (4,))?,
                    vec![bispinor, bispinor],
                ),
                (
                    expression_type.call_method1("projp", (4,))?,
                    vec![bispinor, bispinor],
                ),
                (
                    expression_type.call_method1("sigma", (4,))?,
                    vec![minkowski, minkowski, bispinor, bispinor],
                ),
                (
                    expression_type.call_method1("f", (8,))?,
                    vec![adjoint, adjoint, adjoint],
                ),
                (
                    expression_type.call_method1("t", (8, 3))?,
                    vec![adjoint, fundamental, antifundamental],
                ),
            ];
            let structures = [
                ExplicitKey::<AbstractIndex>::from_iter([euc, euc], ETS.metric, None),
                ExplicitKey::<AbstractIndex>::from_iter([euc, euc], ETS.flat, None),
                AGS.gamma_strct::<AbstractIndex>(4),
                AGS.gamma5_strct::<AbstractIndex>(4),
                AGS.projm_strct::<AbstractIndex>(4),
                AGS.projp_strct::<AbstractIndex>(4),
                AGS.sigma_strct::<AbstractIndex>(4),
                CS.f_strct::<AbstractIndex>(8),
                CS.t_strct::<AbstractIndex>(3, 8),
            ];

            for ((value, expected), structure) in cases.iter().zip(&structures) {
                let tensor = value.extract::<PyRef<'_, TensorExpression>>()?;
                assert_eq!(
                    tensor
                        .interface()
                        .logical_slots()
                        .into_iter()
                        .map(|slot| slot.rep())
                        .collect::<Vec<_>>(),
                    *expected
                );
                assert!(tensor.name_args.is_empty());
                assert_eq!(tensor.interface().open_positions().len(), expected.len());
                let AtomView::Fun(function) = tensor.atom().as_view() else {
                    panic!("a predefined tensor factory must remain an atomic leaf")
                };
                let canonical_slots = function
                    .iter()
                    .map(|argument| {
                        Slot::<LibraryRep, AbstractIndex>::try_from(argument)
                            .expect("factory ports must carry distinct open indices")
                    })
                    .collect::<Vec<_>>();
                assert_eq!(
                    canonical_slots
                        .iter()
                        .map(IsAbstractSlot::rep)
                        .collect::<Vec<_>>(),
                    structure
                        .canonical()
                        .external_reps_iter()
                        .collect::<Vec<_>>()
                );
                assert_eq!(function.get_nargs(), expected.len());
                drop(tensor);

                let logical_slots = expected
                    .iter()
                    .enumerate()
                    .map(|(position, representation)| {
                        representation.slot(AbstractIndex::Normal(100 + position))
                    })
                    .collect::<Vec<_>>();
                let indices = logical_slots
                    .iter()
                    .map(|slot| Py::new(py, SpensoSlot { slot: *slot }))
                    .collect::<PyResult<Vec<_>>>()?;
                let full = PyTuple::new(py, indices.iter())?;
                let indexed = value
                    .call(&full, None)?
                    .extract::<PyRef<'_, TensorExpression>>()?;
                assert_eq!(
                    indexed.interface().logical_slots(),
                    logical_slots
                        .iter()
                        .map(|slot| slot.rep().slot(PartialIndex::Explicit(slot.aind())))
                        .collect::<Vec<_>>()
                );
                let AtomView::Fun(function) = indexed.atom().as_view() else {
                    panic!("an indexed predefined tensor must remain an atomic leaf")
                };
                assert_eq!(
                    function
                        .iter()
                        .map(
                            |argument| Slot::<LibraryRep, AbstractIndex>::try_from(argument)
                                .unwrap()
                        )
                        .collect::<Vec<_>>(),
                    structure.layout().logical_to_canonical(&logical_slots)
                );
                drop(indexed);

                let short = PyTuple::new(py, indices.iter().take(indices.len() - 1))?;
                let error = value
                    .call(&short, None)
                    .expect_err("every predefined factory must enforce its indexing arity");
                assert!(error.is_instance_of::<PyValueError>(py));
            }
            Ok(())
        })
        .unwrap();
    }

    #[test]
    fn builtin_factories_keep_symbolic_dimensions_in_ports_only() {
        idenso::representations::initialize();
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let dimension = Dimension::from(symbol!("factory_symbolic_dimension"));
            let adjoint_dimension = Dimension::from(symbol!("factory_symbolic_adjoint"));
            let fundamental_dimension = Dimension::from(symbol!("factory_symbolic_fundamental"));
            let concrete_four = Dimension::Concrete(4);
            let euc_object = SpensoRepresentation {
                representation: ExtendibleReps::EUCLIDEAN.new_rep(dimension),
            };
            let cases = vec![
                (
                    TensorExpression::g(py, &euc_object, None)?,
                    vec![dimension; 2],
                ),
                (TensorExpression::flat(py, &euc_object)?, vec![dimension; 2]),
                (
                    TensorExpression::gamma(py, ConvertibleToDimension(dimension))?,
                    vec![concrete_four, concrete_four, dimension],
                ),
                (
                    TensorExpression::gamma5(py, ConvertibleToDimension(dimension))?,
                    vec![dimension; 2],
                ),
                (
                    TensorExpression::projm(py, ConvertibleToDimension(dimension))?,
                    vec![dimension; 2],
                ),
                (
                    TensorExpression::projp(py, ConvertibleToDimension(dimension))?,
                    vec![dimension; 2],
                ),
                (
                    TensorExpression::sigma(py, ConvertibleToDimension(dimension))?,
                    vec![dimension, dimension, concrete_four, concrete_four],
                ),
                (
                    TensorExpression::f(py, ConvertibleToDimension(adjoint_dimension))?,
                    vec![adjoint_dimension; 3],
                ),
                (
                    TensorExpression::t(
                        py,
                        ConvertibleToDimension(adjoint_dimension),
                        ConvertibleToDimension(fundamental_dimension),
                    )?,
                    vec![
                        adjoint_dimension,
                        fundamental_dimension,
                        fundamental_dimension,
                    ],
                ),
            ];

            for (expression, expected_dimensions) in cases {
                let expression = expression.bind(py).borrow();
                assert_eq!(
                    expression
                        .interface()
                        .logical_slots()
                        .iter()
                        .map(IsAbstractSlot::dim)
                        .collect::<Vec<_>>(),
                    expected_dimensions
                );
                assert!(expression.name_args.is_empty());
                let AtomView::Fun(function) = expression.atom().as_view() else {
                    panic!("a symbolic-dimension factory must produce one tensor leaf")
                };
                assert_eq!(function.get_nargs(), expected_dimensions.len());
                assert!(
                    function
                        .iter()
                        .all(
                            |argument| Slot::<LibraryRep, AbstractIndex>::try_from(argument)
                                .is_ok()
                        )
                );
            }
            Ok(())
        })
        .unwrap();
    }

    #[test]
    fn raw_predefined_tensor_leaves_use_their_shared_structure() {
        Python::initialize();
        idenso::representations::initialize();
        let minkowski = Minkowski {}.new_rep(Dimension::Concrete(4));
        let malformed = FunctionBuilder::new(CS.t)
            .add_args([
                minkowski.to_symbolic([]),
                minkowski.to_symbolic([]),
                minkowski.to_symbolic([]),
            ])
            .finish();
        let error = infer_interface(&malformed).expect_err("invalid t ports must be rejected");
        assert!(error.to_string().contains("TensorExpression.t"));

        let generator_structure = CS.t_strct::<AbstractIndex>(3, 8);
        let expected = generator_structure.layout().canonical_to_logical(
            &generator_structure
                .canonical()
                .external_reps_iter()
                .collect::<Vec<_>>(),
        );
        let generator = SymbolicTensor::from_signature(&generator_structure).unwrap();
        let inferred = infer_interface(generator.expression()).unwrap();
        assert_eq!(
            inferred
                .logical_slots()
                .into_iter()
                .map(|slot| slot.rep())
                .collect::<Vec<_>>(),
            expected
        );

        let bispinor: Representation<LibraryRep> =
            Bispinor {}.new_rep(Dimension::Concrete(4)).cast();
        let gamma = FunctionBuilder::new(AGS.gamma)
            .add_args([
                minkowski
                    .slot::<AbstractIndex, _>(AbstractIndex::Normal(1))
                    .to_atom(),
                bispinor
                    .slot::<AbstractIndex, _>(AbstractIndex::Normal(2))
                    .to_atom(),
                bispinor
                    .slot::<AbstractIndex, _>(AbstractIndex::Normal(3))
                    .to_atom(),
            ])
            .finish();
        let error = infer_interface(&gamma)
            .expect_err("raw predefined tensors must use canonical atom ordering");
        assert!(error.to_string().contains("TensorExpression.gamma"));
    }

    #[test]
    fn tensor_operand_reinfers_tensor_bearing_plain_expressions() {
        idenso::representations::initialize();
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let representation = ExtendibleReps::MINKOWSKI.new_rep(Dimension::Concrete(4));
            let index = AbstractIndex::Normal(61);
            let slot = representation.slot::<AbstractIndex, _>(index);
            let atom = FunctionBuilder::new(SPENSO_TAG.tensor_symbol("plain_tensor_operand"))
                .add_arg(slot.to_atom())
                .finish();
            let expression = PythonExpression { expr: atom }.into_pyobject(py)?;

            let TensorOperand::Structured(value) = TensorOperand::extract(expression.as_any())?
            else {
                panic!("tagged tensor expression was classified as scalar")
            };
            assert_eq!(
                value.structure().logical_slots(),
                explicit_interface(representation, index).logical_slots()
            );
            Ok(())
        })
        .unwrap();
    }

    #[test]
    fn indexing_preserves_data_identity_only_without_contraction() {
        idenso::representations::initialize();
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let representation = ExtendibleReps::EUCLIDEAN.new_rep(Dimension::Concrete(3));
            let interface = PartialStructure::from_logical_slots([
                representation.slot(PartialIndex::open(0)),
                representation.slot(PartialIndex::open(1)),
            ]);
            let name = SPENSO_TAG.tensor_symbol("indexed_data_identity");
            let argument = Atom::num(7);
            let atom = FunctionBuilder::new(name)
                .add_arg(argument.clone())
                .add_args(
                    interface
                        .logical_slots()
                        .into_iter()
                        .map(composition::port_atom),
                )
                .finish();
            let expression = TensorExpression::from_known_parts(
                py,
                atom,
                interface,
                Some(name),
                vec![argument.clone()],
            )?;

            let relabeled = expression.bind(py).call1(("i", "j"))?;
            let relabeled = relabeled.extract::<PyRef<'_, TensorExpression>>()?;
            assert_eq!(relabeled.interface().canonical().order(), 2);
            assert_eq!(relabeled.name, Some(name));
            assert_eq!(relabeled.name_args, vec![argument]);
            drop(relabeled);

            let contracted = expression.bind(py).call1(("i", "i"))?;
            let contracted = contracted.extract::<PyRef<'_, TensorExpression>>()?;
            assert!(contracted.interface().canonical().is_scalar());
            assert_eq!(contracted.name, None);
            assert!(contracted.name_args.is_empty());
            Ok(())
        })
        .unwrap();
    }

    #[test]
    fn tensor_constructor_cooks_wrapped_indices_and_restores_encoded_tensors() {
        idenso::representations::initialize();
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let representation = SpensoRepresentation {
                representation: ExtendibleReps::EUCLIDEAN
                    .new_rep(Dimension::Concrete(3))
                    .cast(),
            };
            let metric = TensorExpression::g(py, &representation, None)?;
            let tensor_type = py.get_type::<TensorExpression>();
            let kwargs = pyo3::types::PyDict::new(py);
            kwargs.set_item("cook_indices", PyCookSettings::indices())?;
            for (right, rank) in [("mu", 0), ("nu", 2)] {
                let indexed = metric.bind(py).call1(("mu", right))?;
                let wrapped = indexed.call_method1(
                    "wrap_indices",
                    (PythonExpression {
                        expr: Atom::var(symbol!("unified_index_wrapper")),
                    },),
                )?;
                assert!(wrapped.is_instance_of::<TensorExpression>());
                let restored = tensor_type.call1((&wrapped,))?;
                let restored = restored.extract::<PyRef<'_, TensorExpression>>()?;
                let scoped = wrapped.extract::<PyRef<'_, TensorExpression>>()?;
                assert_eq!(restored.atom(), scoped.atom());
                assert_eq!(restored.interface(), scoped.interface());
                assert_eq!(scoped.interface().canonical().order(), rank);
                let cooked = tensor_type.call((&wrapped,), Some(&kwargs))?;
                let cooked = cooked.extract::<Py<TensorExpression>>()?;
                assert_eq!(cooked.borrow(py).interface().canonical().order(), rank);
                if rank == 0 {
                    let simplified = cooked
                        .bind(py)
                        .call_method0("contract")?
                        .call_method0("to_expression")?;
                    let simplified = simplified.extract::<PyRef<'_, TensorExpression>>()?;
                    assert_eq!(*simplified.atom(), Atom::num(3));
                }
                let indexed = indexed.extract::<PyRef<'_, TensorExpression>>()?;
                let settings = idenso::CookSettings::reversible();
                let encoded = settings.try_cook(indexed.atom().as_view()).unwrap();
                let restored = TensorExpression::from_atom_interface(
                    py,
                    settings.uncook(encoded.as_view()),
                    None,
                )?;
                let restored = restored.borrow(py);
                assert_eq!(restored.atom(), indexed.atom());
                // The symmetric metric can normalize its arguments into a different
                // order. Bare cooking cannot retain the separately declared layout.
                assert_eq!(
                    restored.interface().canonical(),
                    indexed.interface().canonical()
                );
                let declared_layout = pyo3::types::PyDict::new(py);
                declared_layout.set_item("structure", indexed.structure())?;
                let restored = tensor_type.call(
                    (PythonExpression {
                        expr: restored.atom().clone(),
                    },),
                    Some(&declared_layout),
                )?;
                let restored = restored.extract::<PyRef<'_, TensorExpression>>()?;
                assert_eq!(restored.atom(), indexed.atom());
                assert_eq!(restored.interface(), indexed.interface());
            }
            Ok(())
        })
        .unwrap();
    }

    #[test]
    fn nonsymmetric_tensor_cooking_preserves_argument_and_port_order() {
        idenso::representations::initialize();
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let representation = ExtendibleReps::EUCLIDEAN
                .new_rep(Dimension::Concrete(3))
                .cast::<LibraryRep>();
            let structure = ExplicitKey::<AbstractIndex>::from_iter(
                [representation, representation],
                SPENSO_TAG.tensor_symbol("nonsymmetric_cooking_tensor"),
                None,
            );
            let tensor = TensorExpression::from_structure(py, &structure)?;
            let mut atoms = Vec::new();
            for indices in [("mu", "nu"), ("nu", "mu")] {
                let indexed = tensor.bind(py).call1(indices)?;
                let indexed = indexed.extract::<PyRef<'_, TensorExpression>>()?;
                let settings = idenso::CookSettings::reversible();
                let encoded = settings.try_cook(indexed.atom().as_view()).unwrap();
                let restored = TensorExpression::from_atom_interface(
                    py,
                    settings.uncook(encoded.as_view()),
                    None,
                )?;
                let restored = restored.borrow(py);
                assert_eq!(restored.atom(), indexed.atom());
                assert_eq!(restored.interface(), indexed.interface());
                atoms.push(indexed.atom().clone());
            }
            assert_ne!(atoms[0], atoms[1]);
            Ok(())
        })
        .unwrap();
    }

    #[test]
    fn typed_dummy_scopes_preserve_free_indices() {
        idenso::representations::initialize();
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let representation = Minkowski {}.new_rep(4);
            let mu = representation
                .slot::<AbstractIndex, _>(AbstractIndex::Normal(17))
                .to_atom();
            let nu = representation
                .slot::<AbstractIndex, _>(AbstractIndex::Normal(19))
                .to_atom();
            let momentum = spenso::vector_symbol!("wrapped_dummy_momentum");
            let atom = ETS.metric(&mu, &nu) * FunctionBuilder::new(momentum).add_arg(&mu).finish();
            let tensor = TensorExpression::from_atom_interface(py, atom, None)?;
            let header = symbol!("wrapped_dummy_payload");
            let wrapped = TensorExpression::wrap_indices(tensor.borrow(py), py, header, true)?;
            let wrapped_mu = representation
                .slot::<AbstractIndex, _>(AbstractIndex::Normal(17).scoped(header))
                .to_atom();
            assert_eq!(
                *wrapped.borrow(py).atom(),
                ETS.metric(&wrapped_mu, &nu)
                    * FunctionBuilder::new(momentum).add_arg(&wrapped_mu).finish()
            );
            assert_eq!(
                wrapped.borrow(py).interface(),
                tensor.borrow(py).interface()
            );
            let simplified = wrapped
                .bind(py)
                .call_method0("contract")?
                .call_method0("to_expression")?;
            let simplified = simplified.extract::<PyRef<'_, TensorExpression>>()?;
            assert_eq!(
                *simplified.atom(),
                FunctionBuilder::new(momentum).add_arg(nu).finish()
            );
            Ok(())
        })
        .unwrap();
    }

    #[test]
    fn wrap_indices_preserves_tensor_expression() {
        idenso::representations::initialize();
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let representation = ExtendibleReps::EUCLIDEAN.new_rep(Dimension::Concrete(3));
            let index = AbstractIndex::Normal(67);
            let slot = representation.slot::<AbstractIndex, _>(index);
            let atom = FunctionBuilder::new(SPENSO_TAG.tensor_symbol("wrapped_index_expression"))
                .add_arg(slot.to_atom())
                .add_arg(slot.to_atom())
                .finish();
            let interface = PartialStructure::from_logical_slots([
                representation.slot(PartialIndex::Explicit(index)),
                representation.slot(PartialIndex::Explicit(index)),
            ]);
            let expression =
                TensorExpression::from_known_parts(py, atom.clone(), interface, None, Vec::new())?;
            assert!(
                expression
                    .bind(py)
                    .borrow()
                    .interface()
                    .canonical()
                    .is_scalar()
            );

            let header = symbolica::symbol!("wrapped_index_header");
            let wrapped =
                TensorExpression::wrap_indices(expression.bind(py).borrow(), py, header, false)?;
            assert_eq!(*wrapped.borrow(py).atom(), atom.wrap_indices(header));
            assert!(wrapped.bind(py).is_instance_of::<TensorExpression>());
            Ok(())
        })
        .unwrap();
    }

    #[test]
    fn explicit_algebra_preserves_open_indexed_zero_and_scalar_tensors() {
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let open = TensorExpression::gamma(py, ConvertibleToDimension(4.into()))?;
            let indexed = open
                .bind(py)
                .call1(("j", "i", "mu"))?
                .extract::<Py<TensorExpression>>()?;
            let zero = TensorExpression::from_atom_interface(
                py,
                Atom::Zero,
                open.borrow(py).interface().clone(),
            )?;
            let scalar = TensorExpression::from_atom_interface(py, Atom::one(), None)?;
            let x = Atom::var(symbol!("algebra_x"));
            let y = Atom::var(symbol!("algebra_y"));
            let denominator = &x + Atom::one();
            let coefficient =
                (Atom::num(2) * &x + Atom::num(4) * x.clone().pow(2) + Atom::num(6) * &y)
                    / &denominator
                    + (x.clone().pow(2) - Atom::one()) / denominator;
            let x_py = PythonExpression { expr: x }.into_pyobject(py)?;
            let y_py = PythonExpression { expr: y }.into_pyobject(py)?;
            let empty = PyTuple::empty(py);
            let one = PyTuple::new(py, [x_py.as_any()])?;
            let two = PyTuple::new(py, [x_py.as_any(), y_py.as_any()])?;
            let variables = PyTuple::new(py, [&two])?;
            let name = symbol!("algebra_tensor_data");
            for tensor in [open, indexed, zero, scalar] {
                let tensor = tensor.borrow(py);
                let expression = TensorExpression::from_known_parts(
                    py,
                    &coefficient * tensor.atom(),
                    tensor.interface().clone(),
                    Some(name),
                    vec![Atom::num(7)],
                )?;
                let original = expression.borrow(py);
                for (method, args) in [
                    ("expand_num", &empty),
                    ("collect_num", &empty),
                    ("collect_factors", &empty),
                    ("collect_by_coefficient", &empty),
                    ("collect_symbol", &one),
                    ("collect_horner", &empty),
                    ("collect_horner", &variables),
                    ("together", &empty),
                    ("cancel", &empty),
                    ("apart", &empty),
                    ("apart", &one),
                    ("apart", &two),
                ] {
                    let raw = TensorExpression::to_expression(expression.borrow(py))
                        .into_pyobject(py)?
                        .call_method1(method, args)?
                        .extract::<PythonExpression>()?;
                    let result = TensorExpression::preserving_interface(&original, py, raw.expr)?;
                    let result = result.borrow(py);
                    assert_eq!(result.interface(), original.interface(), "{method}");
                    assert_eq!(result.name, original.name, "{method}");
                    assert_eq!(result.name_args, original.name_args, "{method}");
                    assert!(
                        (result.atom().clone() - original.atom().clone())
                            .together()
                            .cancel()
                            .expand()
                            .is_zero(),
                        "{method}"
                    );
                }
                let copied = py
                    .import("copy")?
                    .call_method1("copy", (expression.bind(py),))?;
                let copied = copied.extract::<PyRef<'_, TensorExpression>>()?;
                assert_eq!(copied.interface(), original.interface());
                assert_eq!(copied.name, original.name);
                assert_eq!(copied.name_args, original.name_args);
                assert_eq!(copied.atom(), original.atom());
            }
            Ok(())
        })
        .unwrap();
    }

    #[test]
    fn collect_symbol_rejects_callbacks_that_strip_the_tensor_interface() {
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let tensor = TensorExpression::gamma(py, ConvertibleToDimension(4.into()))?;
            let x = PythonExpression {
                expr: Atom::var(symbol!("collect_symbol_x")),
            };
            let expression = TensorExpression::from_atom_interface(
                py,
                &x.expr * tensor.borrow(py).atom(),
                tensor.borrow(py).interface().clone(),
            )?;
            let kwargs = pyo3::types::PyDict::new(py);
            let callback = pyo3::types::PyCFunction::new_closure(py, None, None, |_, _| {
                Ok::<_, PyErr>(PythonExpression { expr: Atom::one() })
            })?;
            kwargs.set_item("coeff_map", callback)?;
            let raw = TensorExpression::to_expression(expression.borrow(py))
                .into_pyobject(py)?
                .call_method("collect_symbol", (x,), Some(&kwargs))?
                .extract::<PythonExpression>()?;
            let error =
                TensorExpression::preserving_interface(&expression.borrow(py), py, raw.expr)
                    .unwrap_err();
            assert!(error.is_instance_of::<PyValueError>(py));
            Ok(())
        })
        .unwrap();
    }

    #[test]
    fn tensor_constructor_factor_collect_and_dimension_preserve_metadata() {
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let tensor = TensorExpression::gamma(py, ConvertibleToDimension(4.into()))?;
            let tensor = TensorExpression::with_name(
                tensor.borrow(py),
                py,
                ConvertibleToSpensoName(
                    SpensoName {
                        name: symbol!("tensor_convenience_name"),
                    },
                    Vec::new(),
                ),
            )?;
            let constructed = py
                .get_type::<TensorExpression>()
                .call1((tensor.bind(py),))?;
            let constructed = constructed.extract::<PyRef<'_, TensorExpression>>()?;
            assert_eq!(constructed.interface(), tensor.borrow(py).interface());
            assert_eq!(constructed.name, tensor.borrow(py).name);
            let ordinary = TensorExpression::to_expression(tensor.borrow(py)).into_pyobject(py)?;
            let inferred = py.get_type::<TensorExpression>().call1((&ordinary,))?;
            assert_eq!(
                inferred
                    .extract::<PyRef<'_, TensorExpression>>()?
                    .interface(),
                tensor.borrow(py).interface()
            );

            let x = PythonExpression {
                expr: Atom::var(symbol!("tensor_convenience_x")),
            };
            let atom = (x.expr.clone() + Atom::one()).pow(2) * tensor.borrow(py).atom().clone();
            let expression = TensorExpression::from_known_parts(
                py,
                atom.expand(),
                tensor.borrow(py).interface().clone(),
                tensor.borrow(py).name,
                vec![Atom::num(7)],
            )?;
            let x = x.into_pyobject(py)?;
            let factored = TensorExpression::to_expression(expression.borrow(py))
                .into_pyobject(py)?
                .call_method0("factor")?;
            let collected = factored
                .call_method1("collect", (&x,))?
                .extract::<PythonExpression>()?;
            let collected =
                TensorExpression::preserving_interface(&expression.borrow(py), py, collected.expr)?;
            let collected = collected.borrow(py);
            assert_eq!(collected.interface(), expression.borrow(py).interface());
            assert_eq!(collected.name, expression.borrow(py).name);
            assert_eq!(collected.name_args, vec![Atom::num(7)]);
            assert!(
                (collected.atom().clone() - expression.borrow(py).atom().clone())
                    .expand()
                    .is_zero()
            );

            let kwargs = pyo3::types::PyDict::new(py);
            let strip_tensor = pyo3::types::PyCFunction::new_closure(py, None, None, |_, _| {
                Ok::<_, PyErr>(PythonExpression { expr: Atom::one() })
            })?;
            kwargs.set_item("coeff_map", strip_tensor)?;
            let raw = factored
                .call_method("collect", (&x,), Some(&kwargs))?
                .extract::<PythonExpression>()?;
            let error =
                TensorExpression::preserving_interface(&expression.borrow(py), py, raw.expr)
                    .unwrap_err();
            assert!(error.is_instance_of::<PyValueError>(py));
            let zero = TensorExpression::from_atom_interface(
                py,
                Atom::Zero,
                tensor.borrow(py).interface().clone(),
            )?;
            let raw = TensorExpression::to_expression(zero.borrow(py))
                .into_pyobject(py)?
                .call_method0("factor")?
                .extract::<PythonExpression>()?;
            let zero = TensorExpression::preserving_interface(&zero.borrow(py), py, raw.expr)?;
            assert_eq!(zero.borrow(py).interface(), tensor.borrow(py).interface());

            let d = Dimension::from(symbol!("tensor_convenience_D"));
            let promoted = TensorExpression::with_lorentz_dimension(
                tensor.borrow(py),
                py,
                ConvertibleToDimension(d),
            )?;
            let expected = TensorExpression::gamma(py, ConvertibleToDimension(d))?;
            let promoted_indexed = promoted.bind(py).call1(("i", "j", "mu"))?;
            let expected_indexed = expected.bind(py).call1(("i", "j", "mu"))?;
            assert_eq!(
                promoted_indexed
                    .extract::<PyRef<'_, TensorExpression>>()?
                    .atom(),
                expected_indexed
                    .extract::<PyRef<'_, TensorExpression>>()?
                    .atom(),
            );
            let again = TensorExpression::with_lorentz_dimension(
                promoted.borrow(py),
                py,
                ConvertibleToDimension(d),
            )?;
            assert_eq!(again.borrow(py).atom(), promoted.borrow(py).atom());
            assert_eq!(promoted.borrow(py).name, tensor.borrow(py).name);
            Ok(())
        })
        .unwrap();
    }

    #[test]
    fn expand_and_zero_network_preserve_structured_metadata() {
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let representation = ExtendibleReps::EUCLIDEAN.new_rep(Dimension::Concrete(4));
            let interface =
                PartialStructure::from_logical_slots([representation.slot(PartialIndex::open(0))]);
            let tensor = FunctionBuilder::new(SPENSO_TAG.tensor_symbol("expand_tensor"))
                .add_arg(representation.to_symbolic([]))
                .finish();
            let variable = Atom::var(symbolica::symbol!("expand_variable"));
            let name = SPENSO_TAG.tensor_symbol("expanded_data_name");
            let expression = TensorExpression::from_known_parts(
                py,
                (variable + Atom::num(1)) * tensor,
                interface.clone(),
                Some(name),
                vec![Atom::num(7)],
            )?;

            let expanded = TensorExpression::expand(expression.bind(py).borrow(), py, None, None)?;
            let expanded_ref = expanded.bind(py).borrow();
            assert!(matches!(expanded_ref.atom(), Atom::Add(_)));
            assert_eq!(
                expanded_ref.interface().logical_slots(),
                interface.logical_slots()
            );
            assert_eq!(expanded_ref.name, Some(name));
            assert_eq!(expanded_ref.name_args, vec![Atom::num(7)]);
            drop(expanded_ref);
            let unchanged = TensorExpression::expand(expanded.borrow(py), py, None, None)?;
            assert_eq!(unchanged.borrow(py).atom(), expanded.borrow(py).atom());
            assert_eq!(
                unchanged.borrow(py).interface(),
                expanded.borrow(py).interface()
            );
            assert_eq!(unchanged.borrow(py).name, Some(name));
            assert_eq!(unchanged.borrow(py).name_args, vec![Atom::num(7)]);

            let ordinary = TensorExpression::to_expression(expression.bind(py).borrow());
            let ordinary = ordinary.expand(None, None)?;
            let ordinary = ordinary.into_pyobject(py)?;
            assert!(!ordinary.is_instance_of::<TensorExpression>());

            let zero = TensorExpression::from_atom_interface(py, Atom::Zero, interface.clone())?;
            let expanded_zero = TensorExpression::expand(zero.borrow(py), py, None, None)?;
            assert!(expanded_zero.borrow(py).atom().as_view().is_zero());
            assert_eq!(*expanded_zero.borrow(py).interface(), interface);
            let network = TensorExpression::to_network(zero.bind(py).borrow(), py, None)?;
            assert_eq!(
                network.structure.structure().logical_slots(),
                interface.logical_slots()
            );
            assert!(network.materialized.structure().open_positions().is_empty());
            assert!(network.network.state.is_tensor());
            Ok(())
        })
        .unwrap();
    }

    #[test]
    fn metric_identities_retain_axis_order_and_tensor_zeros() {
        idenso::representations::initialize();
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let rep = ExtendibleReps::EUCLIDEAN.new_rep(Dimension::Concrete(3));
            let other = ExtendibleReps::MINKOWSKI.new_rep(Dimension::Concrete(4));
            let [i, j] = [31, 37].map(AbstractIndex::Normal);
            let k = AbstractIndex::Normal(41);
            let tensor = |index| {
                FunctionBuilder::new(spenso::tensor_symbol!("retained_metric_tensor"))
                    .add_arg(rep.slot::<AbstractIndex, _>(index).to_atom())
                    .add_arg(other.slot::<AbstractIndex, _>(k).to_atom())
                    .finish()
            };
            let metric = |left, right| {
                FunctionBuilder::new(ETS.metric)
                    .add_arg(rep.slot::<AbstractIndex, _>(left).to_atom())
                    .add_arg(rep.slot::<AbstractIndex, _>(right).to_atom())
                    .finish()
            };
            let interface = PartialStructure::from_logical_slots([
                other.slot(PartialIndex::Explicit(k)),
                rep.slot(PartialIndex::Explicit(i)),
            ]);
            let name = symbol!("retained_metric_data");
            for (atom, expected) in [
                (metric(i, j) * tensor(j), tensor(i)),
                ((metric(j, j) - Atom::num(3)) * tensor(i), Atom::Zero),
            ] {
                let expression = TensorExpression::from_known_parts(
                    py,
                    atom,
                    interface.clone(),
                    Some(name),
                    vec![Atom::num(7)],
                )?;
                let result = expression
                    .bind(py)
                    .call_method0("contract")?
                    .call_method0("to_expression")?;
                let result = result.extract::<PyRef<'_, TensorExpression>>()?;
                assert_eq!(*result.atom(), expected, "contract");
                assert_eq!(
                    result.interface().logical_slots(),
                    interface.logical_slots(),
                    "contract"
                );
                assert_eq!(result.name, Some(name));
                assert_eq!(result.name_args, vec![Atom::num(7)]);
                let dangling = TensorExpression::list_dangling(result)?;
                assert_eq!(
                    dangling
                        .into_iter()
                        .map(|value| value.expr)
                        .collect::<Vec<_>>(),
                    interface
                        .logical_slots()
                        .into_iter()
                        .map(composition::port_atom)
                        .collect::<Vec<_>>()
                );
            }
            Ok(())
        })
        .unwrap();
    }

    #[test]
    fn gamma_and_dot_rewrites_retain_permuted_interfaces() {
        idenso::representations::initialize();
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let mink = ExtendibleReps::MINKOWSKI.new_rep(Dimension::Concrete(4));
            let bis = Bispinor {}.new_rep(4).cast::<LibraryRep>();
            let [a, b, mu, nu] = [43, 47, 53, 59].map(AbstractIndex::Normal);
            let gamma = |left, right, lorentz| {
                FunctionBuilder::new(AGS.gamma)
                    .add_arg(bis.slot::<AbstractIndex, _>(left).to_atom())
                    .add_arg(bis.slot::<AbstractIndex, _>(right).to_atom())
                    .add_arg(mink.slot::<AbstractIndex, _>(lorentz).to_atom())
                    .finish()
            };
            let gamma_product = gamma(a, b, mu) * gamma(b, a, nu);
            let ports = [nu, mu].map(|index| mink.slot(PartialIndex::Explicit(index)));
            let gamma_expression = TensorExpression::from_atom_interface(
                py,
                gamma_product.clone(),
                PartialStructure::from_logical_slots(ports),
            )?;
            for output in [
                idenso::dirac::GammaOutput::Reduced,
                idenso::dirac::GammaOutput::Chains,
            ] {
                let settings = crate::simplification::PyGammaSimplifySettings::from_rust(
                    idenso::dirac::GammaSimplifySettings {
                        output,
                        ..Default::default()
                    },
                );
                let result = gamma_expression.bind(py).call_method1(
                    "reduce_algebra",
                    (Py::new(
                        py,
                        crate::simplification::PyAlgebraSettings::new(
                            Some(settings),
                            None,
                            output == idenso::dirac::GammaOutput::Reduced,
                            16,
                        ),
                    )?,),
                )?;
                let result = result.call_method0("to_expression")?;
                let result = result.extract::<PyRef<'_, TensorExpression>>()?;
                assert_ne!(
                    *result.atom(),
                    gamma_product,
                    "fixture must exercise {output:?}"
                );
                assert_eq!(result.interface().logical_slots(), ports, "{output:?}");
            }
            let vector = |name| {
                FunctionBuilder::new(name)
                    .add_arg(mink.to_symbolic([]))
                    .finish()
            };
            let dot = FunctionBuilder::new(SPENSO_TAG.dot)
                .add_arg(vector(spenso::vector_symbol!("retained_dot_p")))
                .add_arg(vector(spenso::vector_symbol!("retained_dot_q")))
                .finish();
            let tensor = FunctionBuilder::new(spenso::tensor_symbol!("retained_dot_tensor"))
                .add_args([mu, nu].map(|index| mink.slot::<AbstractIndex, _>(index).to_atom()))
                .finish();
            let atom = dot * tensor;
            let expression = TensorExpression::from_atom_interface(
                py,
                atom.clone(),
                PartialStructure::from_logical_slots(ports),
            )?;
            let result = expression.bind(py).call_method0("undo_dots")?;
            let result = result.extract::<PyRef<'_, TensorExpression>>()?;
            assert_ne!(*result.atom(), atom, "fixture must exercise undo_dots");
            assert_eq!(result.interface().logical_slots(), ports);
            Ok(())
        })
        .unwrap();
    }

    #[test]
    fn metric_only_contract_keeps_explicit_vectors_until_full_contraction() {
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let mu = spenso::mink!(4, 94431);
            let nu = spenso::mink!(4, 94432);
            let source = spenso::g!(&mu, &nu) * spenso::p!(&nu) * spenso::q!(&mu);
            let tensor = TensorExpression::from_atom_interface(py, source, None)?;
            let kwargs = PyDict::new(py);
            kwargs.set_item("rank_one", false)?;
            kwargs.set_item("collect_chains", false)?;
            kwargs.set_item("collect_traces", false)?;
            let policy = py
                .get_type::<crate::simplification::PyContractSettings>()
                .call((), Some(&kwargs))?;
            let metric_only = tensor.bind(py).call_method1("contract", (&policy,))?;
            let value = metric_only.extract::<PyRef<'_, AliasedTensorExpression>>()?;
            assert_eq!(
                value
                    .value
                    .resolved()
                    .map_err(TensorExpression::inference_error)?
                    .expression(),
                &(spenso::p!(&nu) * spenso::q!(&nu))
            );
            assert!(!value.value.contraction_complete());
            let again = metric_only.call_method1("contract", (&policy,))?;
            let again = again.extract::<PyRef<'_, AliasedTensorExpression>>()?;
            assert_eq!(again.value.root(), value.value.root());
            let full = metric_only.call_method0("contract")?;
            let full = full.extract::<PyRef<'_, AliasedTensorExpression>>()?;
            assert_eq!(
                full.value
                    .resolved()
                    .map_err(TensorExpression::inference_error)?
                    .expression(),
                &spenso::g!(spenso::p!(spenso::mink!(4)), spenso::q!(spenso::mink!(4)))
            );
            assert!(full.value.contraction_complete());
            Ok(())
        })
        .unwrap();
    }

    #[test]
    fn cooking_maps_stored_slots_and_contracts_collisions() {
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let rep = ExtendibleReps::EUCLIDEAN.new_rep(Dimension::Concrete(3));
            let named = symbol!("retained_cooking_index", tags = [&SPENSO_TAG.index]);
            let [i, j] = [71, 73].map(|owner| AbstractIndex::Named(named.into(), owner, 2));
            let slots = [i, j].map(|index| rep.slot(PartialIndex::Explicit(index)));
            let atom = FunctionBuilder::new(spenso::tensor_symbol!("retained_cooking_tensor"))
                .add_args(slots.map(composition::port_atom))
                .finish();
            let original = TensorExpression::from_atom_interface(
                py,
                atom.clone(),
                PartialStructure::from_logical_slots(slots.into_iter().rev()),
            )?;
            let settings = PyCookSettings::indices();
            let expected = slots
                .into_iter()
                .rev()
                .map(|slot| {
                    settings
                        .rust()
                        .try_cook_indices(composition::port_atom(slot).as_view())
                        .unwrap()
                })
                .collect::<Vec<_>>();
            let kwargs = PyDict::new(py);
            kwargs.set_item("cook_indices", settings.clone())?;
            let result = py
                .get_type::<TensorExpression>()
                .call((&original,), Some(&kwargs))?;
            let result = result.extract::<PyRef<'_, TensorExpression>>()?;
            assert_ne!(*result.atom(), atom);
            assert_eq!(
                result
                    .interface()
                    .logical_slots()
                    .into_iter()
                    .map(composition::port_atom)
                    .collect::<Vec<_>>(),
                expected
            );
            // A typed zero has no index syntax to rewrite, but its interface
            // still needs the same cooking operation, including mixed open ports.
            let mut zero_slots = slots.into_iter().rev().collect::<Vec<_>>();
            zero_slots.push(rep.slot(PartialIndex::open(2)));
            let zero = TensorExpression::from_atom_interface(
                py,
                Atom::Zero,
                PartialStructure::from_logical_slots(zero_slots),
            )?;
            let cooked_zero = Py::new(
                py,
                TensorExpression::new(py, zero.bind(py).as_any(), None, Some(&settings))?,
            )?;
            let cooked_zero = cooked_zero.borrow(py);
            assert!(cooked_zero.atom().is_zero());
            assert_eq!(cooked_zero.interface().open_positions(), vec![2]);
            assert_eq!(
                cooked_zero.interface().logical_slots()[..2]
                    .iter()
                    .copied()
                    .map(composition::port_atom)
                    .collect::<Vec<_>>(),
                expected
            );
            // Flattening may identify a previously named index with an existing
            // scalar label. Derive that contraction from the mapped slots.
            let original_slot = rep.slot::<AbstractIndex, _>(i).to_atom();
            let cooked_slot = settings
                .rust()
                .try_cook_indices(original_slot.as_view())
                .unwrap();
            let metric = FunctionBuilder::new(ETS.metric)
                .add_arg(original_slot)
                .add_arg(cooked_slot)
                .finish();
            let metric = TensorExpression::from_atom_interface(py, metric, None)?;
            assert_eq!(metric.borrow(py).interface().canonical().order(), 2);
            let cooked = Py::new(
                py,
                TensorExpression::new(py, metric.bind(py).as_any(), None, Some(&settings))?,
            )?;
            assert!(cooked.borrow(py).interface().canonical().is_scalar());
            let contraction = PyDict::new(py);
            contraction.set_item("collect_chains", false)?;
            contraction.set_item("collect_traces", false)?;
            let contraction = py
                .get_type::<crate::simplification::PyContractSettings>()
                .call((), Some(&contraction))?;
            let contraction =
                contraction.extract::<PyRef<'_, crate::simplification::PyContractSettings>>()?;
            let simplified = TensorExpression::contract(&cooked.borrow(py), Some(&contraction))?;
            assert_eq!(
                simplified
                    .value
                    .resolved()
                    .map_err(TensorExpression::inference_error)?
                    .expression(),
                &Atom::num(3)
            );
            Ok(())
        })
        .unwrap();
    }

    #[test]
    fn simplify_expansion_keeps_the_same_order_as_expansion() {
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let rep = ExtendibleReps::EUCLIDEAN.new_rep(Dimension::Concrete(2));
            let other = ExtendibleReps::MINKOWSKI.new_rep(Dimension::Concrete(3));
            let slots = [
                rep.slot(PartialIndex::Explicit(AbstractIndex::Normal(61))),
                other.slot(PartialIndex::Explicit(AbstractIndex::Normal(67))),
            ];
            let tensor = FunctionBuilder::new(spenso::tensor_symbol!("retained_pipeline_tensor"))
                .add_args(slots.map(composition::port_atom))
                .finish();
            let x = Atom::var(symbol!("retained_pipeline_x"));
            let atom = (x + Atom::one()) * tensor;
            let interface = PartialStructure::from_logical_slots(slots.into_iter().rev());
            let original =
                TensorExpression::from_atom_interface(py, atom.clone(), interface.clone())?;
            let expanded = TensorExpression::expand(original.borrow(py), py, None, None)?;
            let settings = crate::simplification::PyAlgebraSettings::new(None, None, false, 16);
            let simplified = original.borrow(py).reduce_algebra(Some(&settings))?;
            let simplified = simplified
                .value
                .expanded()
                .map_err(TensorExpression::inference_error)?;
            assert_ne!(simplified.expression(), &atom);
            assert_eq!(simplified.expression(), expanded.borrow(py).atom());
            assert_eq!(
                simplified.structure().logical_slots(),
                interface.logical_slots()
            );
            Ok(())
        })
        .unwrap();
    }

    #[test]
    fn constructor_rejects_inconsistent_supplied_parts() {
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let rep = ExtendibleReps::EUCLIDEAN.new_rep(Dimension::Concrete(3));
            let slot = rep.slot(PartialIndex::Explicit(AbstractIndex::Normal(101)));
            // Duplicate indices contract; a declared extra port cannot be invented.
            let structure = crate::metadata::SpensoTensorStructure {
                interface: PartialStructure::from_logical_slots([
                    slot,
                    slot,
                    rep.slot(PartialIndex::Explicit(AbstractIndex::Normal(103))),
                ]),
                name: Some(symbol!("unsafe_supplied_data")),
                arguments: vec![Atom::num(7)],
            };
            let kwargs = PyDict::new(py);
            kwargs.set_item("structure", structure.clone())?;
            let malformed = FunctionBuilder::new(SPENSO_TAG.dot)
                .add_arg(Atom::one())
                .finish();
            for atom in [Atom::one(), malformed] {
                assert!(
                    TensorExpression::from_atom_interface(
                        py,
                        atom.clone(),
                        structure.interface.clone()
                    )
                    .is_err()
                );
                let expression = PythonExpression { expr: atom.clone() }.into_pyobject(py)?;
                assert!(
                    py.get_type::<TensorExpression>()
                        .call((expression,), Some(&kwargs))
                        .is_err()
                );
            }
            Ok(())
        })
        .unwrap();
    }

    #[test]
    fn algebra_preserves_validated_interfaces() {
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let rep = ExtendibleReps::EUCLIDEAN.new_rep(Dimension::Concrete(3));
            let slots = [79, 83]
                .map(|index| rep.slot(PartialIndex::Explicit(AbstractIndex::Normal(index))));
            let tensor = FunctionBuilder::new(spenso::tensor_symbol!("constructed_algebra_tensor"))
                .add_args(slots.map(composition::port_atom))
                .finish();
            let interface = PartialStructure::from_logical_slots(slots.into_iter().rev());
            let x = Atom::var(symbol!("constructed_algebra_x"));
            let y = Atom::var(symbol!("constructed_algebra_y"));
            let coefficient =
                (Atom::num(2) * &x + Atom::num(4) * x.clone().pow(2) + Atom::num(6) * &y)
                    / (&x + Atom::one())
                    + (x.clone().pow(2) - Atom::one()) / (&x + Atom::one());
            let x_py = PythonExpression { expr: x.clone() }.into_pyobject(py)?;
            let empty = PyTuple::empty(py);
            let one = PyTuple::new(py, [x_py.as_any()])?;
            for atom in [&coefficient * &tensor, tensor, Atom::Zero] {
                let original =
                    TensorExpression::from_atom_interface(py, atom.clone(), interface.clone())?;
                for (method, args) in [
                    ("factor", &empty),
                    ("expand", &empty),
                    ("expand_num", &empty),
                    ("collect", &one),
                    ("collect_num", &empty),
                    ("collect_factors", &empty),
                    ("collect_by_coefficient", &empty),
                    ("collect_symbol", &one),
                    ("collect_horner", &empty),
                    ("together", &empty),
                    ("cancel", &empty),
                    ("apart", &one),
                ] {
                    let raw = TensorExpression::to_expression(original.borrow(py))
                        .into_pyobject(py)?
                        .call_method1(method, args)?
                        .extract::<PythonExpression>()?;
                    let result =
                        TensorExpression::preserving_interface(&original.borrow(py), py, raw.expr)?;
                    let result = result.borrow(py);
                    assert_eq!(
                        result.interface().logical_slots(),
                        interface.logical_slots(),
                        "{method}"
                    );
                    assert!(
                        (result.atom() - &atom)
                            .together()
                            .cancel()
                            .expand()
                            .is_zero(),
                        "{method}"
                    );
                }
                for via_poly in [false, true] {
                    let expanded =
                        TensorExpression::expand(original.borrow(py), py, None, Some(via_poly))?;
                    assert_eq!(
                        expanded.borrow(py).interface().logical_slots(),
                        interface.logical_slots()
                    );
                }
                let settings = crate::simplification::PyAlgebraSettings::new(None, None, false, 16);
                let simplified = original.borrow(py).reduce_algebra(Some(&settings))?;
                let simplified = simplified
                    .value
                    .expanded()
                    .map_err(TensorExpression::inference_error)?;
                assert_eq!(
                    simplified.structure().logical_slots(),
                    interface.logical_slots()
                );
            }
            Ok(())
        })
        .unwrap();
    }

    #[test]
    fn tensor_name_validates_structure_changing_normalizers_at_construction() {
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let name = spenso::tensor_symbol!(
                "constructor_strips_ports",
                norm = |_, output| {
                    **output = Atom::one();
                }
            );
            let tensor_name = Py::new(py, SpensoName { name })?;
            let representation = SpensoRepresentation {
                representation: ExtendibleReps::EUCLIDEAN.new_rep(Dimension::Concrete(3)),
            };
            let error = tensor_name
                .bind(py)
                .call1((representation,))
                .expect_err("a normalizer cannot silently strip a tensor's ports");
            assert!(error.is_instance_of::<PyValueError>(py));
            assert!(error.to_string().contains("interface"));
            Ok(())
        })
        .unwrap();
    }

    #[test]
    fn untrusted_interfaces_are_checked_and_known_composition_preserves_structure() {
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let rep = ExtendibleReps::EUCLIDEAN.new_rep(Dimension::Concrete(3));
            let slot = rep.slot(PartialIndex::Explicit(AbstractIndex::Normal(89)));
            let atom = FunctionBuilder::new(spenso::tensor_symbol!("constructor_checked_tensor"))
                .add_arg(composition::port_atom(slot))
                .finish();
            let interface = PartialStructure::from_logical_slots([slot]);
            let tensor =
                TensorExpression::from_atom_interface(py, atom.clone(), interface.clone())?;
            assert!(
                TensorExpression::from_atom_interface(
                    py,
                    atom,
                    PartialStructure::from_logical_slots([])
                )
                .is_err()
            );
            let mismatched = rep.slot(PartialIndex::Explicit(AbstractIndex::Normal(97)));
            assert!(
                TensorExpression::from_atom_interface(
                    py,
                    tensor.borrow(py).atom().clone(),
                    PartialStructure::from_logical_slots([mismatched])
                )
                .is_err()
            );
            let number = PythonExpression { expr: Atom::num(2) }.into_pyobject(py)?;
            let doubled = tensor.bind(py).call_method1("__mul__", (number,))?;
            let doubled = doubled.extract::<PyRef<'_, TensorExpression>>()?;
            assert_eq!(
                doubled.interface().logical_slots(),
                interface.logical_slots()
            );
            let original = tensor.borrow(py);
            let invalid_rewrite =
                TensorExpression::preserving_interface(&original, py, Atom::one());
            assert!(invalid_rewrite.is_err());
            let unchanged =
                TensorExpression::preserving_interface(&original, py, original.atom().clone())?;
            assert_eq!(
                unchanged.borrow(py).interface().logical_slots(),
                interface.logical_slots()
            );
            Ok(())
        })
        .unwrap();
    }

    #[test]
    fn expansion_preserves_checked_logical_ports_names_and_scalar_dots() {
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let rep = ExtendibleReps::MINKOWSKI.new_rep(Dimension::Concrete(4));
            let slots = [31, 37]
                .map(|index| rep.slot(PartialIndex::Explicit(AbstractIndex::Normal(index))));
            let tensor = FunctionBuilder::new(spenso::tensor_symbol!("expansion_checked_tensor"))
                .add_args(slots.map(composition::port_atom))
                .finish();
            let p = FunctionBuilder::new(spenso::vector_symbol!("expansion_checked_p"))
                .add_arg(rep.to_symbolic([]))
                .finish();
            let q = FunctionBuilder::new(spenso::vector_symbol!("expansion_checked_q"))
                .add_arg(rep.to_symbolic([]))
                .finish();
            let x = Atom::var(symbol!("expansion_checked_x"));
            let coefficient = (&x + Atom::one()) * ETS.metric(p, q).pow(2);
            let name = symbol!("expansion_checked_data");
            let metric = FunctionBuilder::new(ETS.metric)
                .add_args(slots.map(composition::port_atom))
                .finish();
            for leaf in [tensor, metric] {
                let atom = &coefficient * leaf;
                assert!(
                    InterfaceInference::default().algebra_preserves_leaf_interfaces(atom.as_view())
                );
                // An explicit interface may carry a different logical ordering
                // from the encoded leaf. Expansion must retain the supplied order.
                let interface = PartialStructure::from_logical_slots(slots.into_iter().rev());
                let original = TensorExpression::from_known_parts(
                    py,
                    atom.clone(),
                    interface.clone(),
                    Some(name),
                    vec![Atom::num(7)],
                )?;
                let var = PythonExpression { expr: x.clone() }
                    .into_pyobject(py)?
                    .extract::<ConvertibleToExpression>()?;
                for (var, via_poly) in [(None, None), (None, Some(true)), (Some(var), None)] {
                    let expanded =
                        TensorExpression::expand(original.borrow(py), py, var, via_poly)?;
                    let expanded = expanded.borrow(py);
                    assert_eq!(*expanded.atom(), atom.expand());
                    assert_eq!(*expanded.interface(), interface);
                    assert_eq!(expanded.name, Some(name));
                    assert_eq!(expanded.name_args, vec![Atom::num(7)]);
                }
            }
            Ok(())
        })
        .unwrap();
    }

    #[test]
    fn expansion_retains_checked_errors_for_powers_open_ports_and_cancellation() {
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let rep = ExtendibleReps::MINKOWSKI.new_rep(Dimension::Concrete(4));
            let explicit = rep
                .slot::<AbstractIndex, _>(AbstractIndex::Normal(41))
                .to_atom();
            let directed = ColorFundamental {}
                .new_rep(3)
                .slot::<AbstractIndex, _>(AbstractIndex::Normal(43))
                .to_atom();
            let vector = |name: &str, port: &Atom| {
                FunctionBuilder::new(SPENSO_TAG.tensor_symbol(name))
                    .add_arg(port)
                    .finish()
            };
            let [p, q, r, s] = [
                "expansion_error_p",
                "expansion_error_q",
                "expansion_error_r",
                "expansion_error_s",
            ]
            .map(|name| vector(name, &explicit));
            let open = rep.to_symbolic([]);
            let [op, oq, or] = [
                "expansion_error_p",
                "expansion_error_q",
                "expansion_error_r",
            ]
            .map(|name| vector(name, &open));
            let [dp, dq, dr] = [
                "expansion_error_p",
                "expansion_error_q",
                "expansion_error_r",
            ]
            .map(|name| vector(name, &directed));
            let x = Atom::var(symbol!("expansion_error_x"));
            // Reject invalid supplied metadata at construction, even if scalar
            // expansion could later cancel incompatible terms to zero.
            let cancellation = (&p + &x) * (&p - &x) - p.clone().pow(2) + x.clone().pow(2);
            assert!(
                TensorExpression::from_atom_interface(
                    py,
                    cancellation.clone(),
                    PartialStructure::from_logical_slots([])
                )
                .is_err()
            );
            let canceled = TensorExpression::from_atom_interface(
                py,
                cancellation.expand(),
                PartialStructure::from_logical_slots([]),
            )?;
            assert!(canceled.borrow(py).atom().is_zero());
            assert!(
                TensorExpression::from_atom_interface(
                    py,
                    x.clone(),
                    PartialStructure::from_logical_slots([rep.slot(PartialIndex::open(0))])
                )
                .is_err()
            );
            let cases = [
                ((&p * &q + &r * &s).pow(2), None, false),
                ((&op + &oq) * (&op + &or), None, false),
                ((&dp + &dq) * (&dp + &dr), None, false),
                ((p.clone().pow(2) + &x) * &p, None, false),
                (
                    Atom::Zero,
                    Some(PartialStructure::from_logical_slots([
                        rep.slot(PartialIndex::open(0))
                    ])),
                    true,
                ),
            ];
            for (atom, interface, succeeds) in cases {
                let original = TensorExpression::from_atom_interface(py, atom.clone(), interface)?;
                let expected =
                    TensorExpression::preserving_interface(&original.borrow(py), py, atom.expand());
                let actual = TensorExpression::expand(original.borrow(py), py, None, None);
                assert_eq!(actual.is_ok(), succeeds, "{atom}: {actual:?}");
                match (actual, expected) {
                    (Ok(actual), Ok(expected)) => {
                        let actual = actual.borrow(py);
                        let expected = expected.borrow(py);
                        assert_eq!(actual.atom(), expected.atom());
                        assert_eq!(actual.interface(), expected.interface());
                        assert_eq!(actual.name, expected.name);
                        assert_eq!(actual.name_args, expected.name_args);
                    }
                    (Err(actual), Err(expected)) => {
                        assert_eq!(actual.to_string(), expected.to_string())
                    }
                    _ => panic!("expansion changed checked behavior for {atom}"),
                }
            }
            Ok(())
        })
        .unwrap();
    }

    #[test]
    fn expansion_does_not_replay_materialization_callbacks() {
        use std::sync::{Arc, Mutex};
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let calls = Arc::new(Mutex::new(Vec::new()));
            let observed = Arc::clone(&calls);
            let callback = spenso::tensor_symbol!(
                "expansion_materialization_callback",
                norm = move |node, _| {
                    observed.lock().unwrap().push(node.to_owned());
                }
            );
            let rep = ExtendibleReps::MINKOWSKI.new_rep(Dimension::Concrete(4));
            let tensor = FunctionBuilder::new(callback)
                .add_arg(rep.to_symbolic([]))
                .finish();
            // Establish that the observer records actual normalization calls.
            assert_eq!(*calls.lock().unwrap(), vec![tensor.clone()]);
            let x = Atom::var(symbol!("expansion_callback_x"));
            let atom = (x + Atom::one()) * tensor;
            let original = TensorExpression::from_atom_interface(
                py,
                atom.clone(),
                PartialStructure::from_logical_slots([rep.slot(PartialIndex::open(0))]),
            )?;
            calls.lock().unwrap().clear();
            let expanded = atom.expand();
            let expansion_calls = std::mem::take(&mut *calls.lock().unwrap());
            let expected =
                TensorExpression::preserving_interface(&original.borrow(py), py, expanded)?;
            assert!(calls.lock().unwrap().is_empty());
            let actual = TensorExpression::expand(original.borrow(py), py, None, None)?;
            assert_eq!(*calls.lock().unwrap(), expansion_calls);
            let actual = actual.borrow(py);
            let expected = expected.borrow(py);
            assert_eq!(actual.atom(), expected.atom());
            assert_eq!(actual.interface(), expected.interface());
            Ok(())
        })
        .unwrap();
    }

    #[test]
    fn addition_accepts_permuted_logical_interfaces_with_equal_canonical_slots() {
        idenso::representations::initialize();
        let representation = ExtendibleReps::EUCLIDEAN.new_rep(Dimension::Concrete(2));
        let i = AbstractIndex::Normal(31);
        let j = AbstractIndex::Normal(37);
        let tensor = |name: &str, indices: [AbstractIndex; 2]| {
            let interface = PartialStructure::from_logical_slots(
                indices.map(|index| representation.slot(PartialIndex::Explicit(index))),
            );
            let atom = FunctionBuilder::new(SPENSO_TAG.tensor_symbol(name))
                .add_args(
                    interface
                        .logical_slots()
                        .into_iter()
                        .map(composition::port_atom),
                )
                .finish();
            SymbolicTensor::new(atom, interface)
        };
        let left = tensor("canonical_addition_left", [i, j]);
        let right = tensor("canonical_addition_right", [j, i]);
        let expected_sum = left.expression() + right.expression();
        let expected_difference = left.expression() - right.expression();

        assert_ne!(
            left.structure().logical_slots(),
            right.structure().logical_slots()
        );
        assert!(InterfaceInference::additive_interfaces_match(
            left.structure(),
            right.structure()
        ));
        let mixed_left = PartialStructure::from_logical_slots([
            representation.slot(PartialIndex::Explicit(i)),
            representation.slot(PartialIndex::open(0)),
        ]);
        let mixed_right = PartialStructure::from_logical_slots([
            representation.slot(PartialIndex::open(0)),
            representation.slot(PartialIndex::Explicit(i)),
        ]);
        assert_eq!(mixed_left.canonical(), mixed_right.canonical());
        assert!(!InterfaceInference::additive_interfaces_match(
            &mixed_left,
            &mixed_right
        ));
        let open_interface = PartialStructure::from_logical_slots([
            representation.slot(PartialIndex::open(0)),
            representation.slot(PartialIndex::open(1)),
        ]);
        let open_tensor = |name| {
            FunctionBuilder::new(SPENSO_TAG.tensor_symbol(name))
                .add_arg(representation.to_symbolic([]))
                .add_arg(representation.to_symbolic([]))
                .finish()
        };
        let open_sum = open_tensor("canonical_open_addition_left")
            + open_tensor("canonical_open_addition_right");
        assert_eq!(
            InterfaceInference::default()
                .infer_validated(open_sum.as_view())
                .unwrap()
                .logical_slots(),
            open_interface.logical_slots()
        );
        let inferred = InterfaceInference::default()
            .infer_validated(expected_sum.as_view())
            .unwrap();
        assert!(InterfaceInference::additive_interfaces_match(
            left.structure(),
            &inferred
        ));

        let standalone = SymbolicTensor::new(left.expression().clone(), right.structure().clone());
        assert_eq!(
            standalone.presentation_atom(),
            tensor("canonical_addition_left", [j, i]).into_expression()
        );

        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let expected = left.structure().logical_slots();
            let left = TensorExpression::from_structured(py, left)?;
            let right = TensorExpression::from_structured(py, right)?;
            let TensorDispatch::Expression(sum) =
                TensorExpression::__add__(left.bind(py).borrow(), py, right.bind(py).as_any())?
            else {
                panic!("adding symbolic tensor expressions produced a network")
            };
            let sum = sum.bind(py).borrow();
            assert_eq!(sum.interface().logical_slots(), expected);
            assert_eq!(*sum.atom(), expected_sum);
            assert_eq!(
                TensorExpression::structured(&sum).presentation_atom(),
                expected_sum
            );
            drop(sum);

            let TensorDispatch::Expression(difference) =
                TensorExpression::__sub__(left.bind(py).borrow(), py, right.bind(py).as_any())?
            else {
                panic!("subtracting symbolic tensor expressions produced a network")
            };
            let difference = difference.bind(py).borrow();
            assert_eq!(difference.interface().logical_slots(), expected);
            assert_eq!(*difference.atom(), expected_difference);
            assert_eq!(
                TensorExpression::structured(&difference).presentation_atom(),
                expected_difference
            );
            drop(difference);

            let TensorDispatch::Expression(reverse_difference) =
                TensorExpression::__rsub__(right.bind(py).borrow(), py, left.bind(py).as_any())?
            else {
                panic!("reverse-subtracting symbolic tensor expressions produced a network")
            };
            let reverse_difference = reverse_difference.bind(py).borrow();
            assert_eq!(reverse_difference.interface().logical_slots(), expected);
            assert_eq!(*reverse_difference.atom(), expected_difference);
            Ok(())
        })
        .unwrap();
    }

    #[test]
    fn reflected_division_retains_tensor_numerator_interface_and_typed_zero() {
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let ports = [
                ExtendibleReps::MINKOWSKI
                    .new_rep(Dimension::Concrete(4))
                    .slot(PartialIndex::Explicit(AbstractIndex::Normal(73001))),
                ExtendibleReps::EUCLIDEAN
                    .new_rep(Dimension::Concrete(3))
                    .slot(PartialIndex::Explicit(AbstractIndex::Normal(73003))),
            ];
            let interface = PartialStructure::from_logical_slots(ports);
            let expression =
                FunctionBuilder::new(spenso::tensor_symbol!("reflected_division_numerator"))
                    .add_args(ports.into_iter().rev().map(composition::port_atom))
                    .finish();
            let divisor = Atom::num(2);
            let denominator = TensorExpression::from_atom_interface(
                py,
                divisor.clone(),
                PartialStructure::from_logical_slots([]),
            )?;
            for expression in [expression, Atom::Zero] {
                let numerator = TensorExpression::from_atom_interface(
                    py,
                    expression.clone(),
                    interface.clone(),
                )?;
                let TensorDispatch::Expression(quotient) = TensorExpression::__rtruediv__(
                    denominator.bind(py).borrow(),
                    py,
                    numerator.bind(py).as_any(),
                )?
                else {
                    panic!("symbolic reflected division produced a network")
                };
                let quotient = quotient.borrow(py);
                assert_eq!(*quotient.atom(), &expression / &divisor);
                assert_eq!(quotient.interface().logical_slots(), ports);
                assert_eq!(numerator.borrow(py).interface().logical_slots(), ports);
                let error = TensorExpression::__rtruediv__(
                    numerator.bind(py).borrow(),
                    py,
                    denominator.bind(py).as_any(),
                )
                .err()
                .expect("a tensor denominator must be rejected, including typed zero");
                assert!(error.is_instance_of::<PyValueError>(py));
                assert!(error.to_string().contains("non-scalar tensor"));
            }
            Ok(())
        })
        .unwrap();
    }

    #[test]
    fn scalar_zero_is_a_polymorphic_additive_identity() {
        idenso::representations::initialize();
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let representation = ExtendibleReps::EUCLIDEAN.new_rep(Dimension::Concrete(2));
            let interface =
                PartialStructure::from_logical_slots([representation.slot(PartialIndex::open(0))]);
            let name = SPENSO_TAG.tensor_symbol("zero_additive_identity");
            let argument = Atom::num(7);
            let atom = FunctionBuilder::new(name)
                .add_arg(argument.clone())
                .add_arg(representation.to_symbolic([]))
                .finish();
            let expression = TensorExpression::from_known_parts(
                py,
                atom.clone(),
                interface.clone(),
                Some(name),
                vec![argument.clone()],
            )?;
            let assert_identity = |value: &Bound<'_, PyAny>| -> PyResult<()> {
                let value = value.extract::<PyRef<'_, TensorExpression>>()?;
                assert_eq!(*value.atom(), atom);
                assert_eq!(value.interface().logical_slots(), interface.logical_slots());
                assert_eq!(value.name, Some(name));
                assert_eq!(value.name_args, vec![argument.clone()]);
                Ok(())
            };

            for method in ["__add__", "__radd__", "__sub__"] {
                let identity = expression.bind(py).call_method1(method, (0,))?;
                assert_identity(&identity)?;
            }

            let singleton = PyTuple::new(py, [expression.clone_ref(py)])?;
            let sum = py.import("builtins")?.call_method1("sum", (singleton,))?;
            assert_identity(&sum)?;

            let pair = PyTuple::new(py, [expression.clone_ref(py), expression.clone_ref(py)])?;
            let sum = py.import("builtins")?.call_method1("sum", (pair,))?;
            let sum = sum.extract::<PyRef<'_, TensorExpression>>()?;
            assert_eq!(*sum.atom(), atom.as_ref() + atom.as_ref());
            assert_eq!(sum.interface().logical_slots(), interface.logical_slots());
            assert_eq!(sum.name, None);
            assert!(sum.name_args.is_empty());
            drop(sum);

            let negated = expression.bind(py).call_method1("__rsub__", (0,))?;
            let negated = negated.extract::<PyRef<'_, TensorExpression>>()?;
            assert_eq!(*negated.atom(), Atom::num(-1) * atom.as_ref());
            assert_eq!(
                negated.interface().logical_slots(),
                interface.logical_slots()
            );
            assert_eq!(negated.name, None);
            assert!(negated.name_args.is_empty());
            drop(negated);

            for method in ["__add__", "__radd__", "__sub__", "__rsub__"] {
                let error = expression
                    .bind(py)
                    .call_method1(method, (1,))
                    .expect_err("nonzero scalar/tensor addition must remain invalid");
                assert!(error.is_instance_of::<PyValueError>(py));
            }

            let scalar = TensorExpression::from_atom_interface(
                py,
                FunctionBuilder::new(SPENSO_TAG.tensor_symbol("rank_zero_addition")).finish(),
                PartialStructure::from_logical_slots(std::iter::empty()),
            )?;
            assert!(scalar.bind(py).call_method1("__add__", (1,)).is_ok());
            Ok(())
        })
        .unwrap();
    }

    #[test]
    fn rank_zero_tensor_calls_reinfer_as_structured_scalars() {
        let name = SPENSO_TAG.tensor_symbol("rank_zero_tensor_call");
        for atom in [
            FunctionBuilder::new(name).finish(),
            FunctionBuilder::new(name).add_arg(Atom::num(7)).finish(),
        ] {
            assert!(infer_interface(&atom).unwrap().canonical().is_scalar());
        }
    }

    #[test]
    fn tensor_scalar_arguments_must_precede_mixed_ports() {
        let representation = ExtendibleReps::EUCLIDEAN.new_rep(Dimension::Concrete(4));
        let slot = representation
            .slot::<AbstractIndex, _>(AbstractIndex::Normal(29))
            .to_atom();
        let name = SPENSO_TAG.tensor_symbol("mixed_tensor_ports");
        let valid = FunctionBuilder::new(name)
            .add_arg(Atom::num(3))
            .add_arg(&slot)
            .add_arg(representation.to_symbolic([]))
            .finish();
        let invalid = FunctionBuilder::new(name)
            .add_arg(slot)
            .add_arg(Atom::num(3))
            .finish();

        let interface = infer_interface(&valid).unwrap();
        assert_eq!(interface.canonical().order(), 2);
        assert!(matches!(
            interface.logical_slots()[0].aind,
            PartialIndex::Explicit(AbstractIndex::Normal(29))
        ));
        assert!(matches!(
            interface.logical_slots()[1].aind,
            PartialIndex::Open(_)
        ));
        assert!(infer_interface(&invalid).is_err());
    }

    #[test]
    fn explicit_scalar_functions_keep_tensor_metadata_opaque() {
        let scalar = symbol!("opaque_tensor_metadata"; Scalar);
        let vector = spenso::vector_symbol!("opaque_metadata_vector");
        let representation = ExtendibleReps::EUCLIDEAN.new_rep(Dimension::Concrete(4));
        let metadata = FunctionBuilder::new(vector).add_arg(Atom::num(3)).finish();
        let coefficient = FunctionBuilder::new(scalar)
            .add_arg(&metadata)
            .add_arg(representation.to_symbolic([]))
            .finish();
        let tensor = FunctionBuilder::new(vector)
            .add_arg(Atom::num(3))
            .add_arg(representation.to_symbolic([]))
            .finish();

        assert!(!InterfaceInference::has_structured_syntax(
            coefficient.as_view()
        ));
        assert!(infer_interface(&metadata).is_err());
        let product = tensor * coefficient.pow(Atom::num(-2));
        let inferred = SymbolicTensor::infer(product.clone()).unwrap();
        assert_eq!(inferred.structure().canonical().order(), 1);
        assert_eq!(inferred.expression(), &product);
    }

    #[test]
    fn rank_one_tensor_calls_require_one_final_structural_port() {
        let name = spenso::vector_symbol!("rank_one_tensor_call");
        let representation = ExtendibleReps::EUCLIDEAN.new_rep(Dimension::Concrete(4));
        let valid = FunctionBuilder::new(name)
            .add_arg(Atom::num(5))
            .add_arg(representation.to_symbolic([]))
            .finish();
        let missing = FunctionBuilder::new(name).add_arg(Atom::num(5)).finish();
        let extra = FunctionBuilder::new(name)
            .add_arg(representation.to_symbolic([]))
            .add_arg(representation.to_symbolic([]))
            .finish();

        assert_eq!(infer_interface(&valid).unwrap().canonical().order(), 1);
        assert!(infer_interface(&missing).is_err());
        assert!(infer_interface(&extra).is_err());
    }

    #[test]
    fn explicit_color_tensor_squares_contract_all_index_pairs() {
        idenso::representations::initialize();
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let adjoint = ColorAdjoint {}.new_rep(8);
            let [a, b, c] = [1, 2, 3].map(|index| {
                adjoint
                    .slot::<AbstractIndex, _>(AbstractIndex::Normal(index))
                    .to_atom()
            });
            let atom = idenso::color_f!(&a, &b, &c) * idenso::color_f!(&c, &b, &a);
            let expression = PythonExpression { expr: atom.clone() }.into_pyobject(py)?;
            let tensor = py.get_type::<TensorExpression>().call1((expression,))?;
            let tensor = tensor.extract::<Py<TensorExpression>>()?;
            assert!(tensor.borrow(py).interface().canonical().is_scalar());
            assert_eq!(*tensor.borrow(py).atom(), atom);
            let simplified = tensor
                .bind(py)
                .call_method1(
                    "reduce_algebra",
                    (Py::new(
                        py,
                        crate::simplification::PyAlgebraSettings::new(
                            None,
                            Some(crate::simplification::PyColorSimplifySettings::from_rust(
                                Default::default(),
                            )),
                            false,
                            16,
                        ),
                    )?,),
                )?
                .call_method0("to_expression")?;
            let simplified = simplified.extract::<PyRef<'_, TensorExpression>>()?;
            assert!(simplified.interface().canonical().is_scalar());
            assert_eq!(
                *simplified.atom(),
                Atom::num(-8) * CS.cas(Atom::num(2), adjoint.to_symbolic([]))
            );

            let open = TensorExpression::f(py, ConvertibleToDimension(8.into()))?;
            let square = open.borrow(py).atom().as_ref().pow(Atom::num(2));
            let error = TensorExpression::from_atom_interface(py, square, None)
                .expect_err("a bare rank-three square still has ambiguous port pairings");
            assert!(error.to_string().contains("ambiguous"));
            Ok(())
        })
        .unwrap();
    }

    #[test]
    fn chain_placeholders_require_a_factor_scope() {
        Python::initialize();
        let placeholder = Atom::var(SPENSO_TAG.chain_in);
        assert!(
            InterfaceInference::validate_placeholder_scope(placeholder.as_view(), false).is_err()
        );

        let representation = ExtendibleReps::EUCLIDEAN.new_rep(Dimension::Concrete(4));
        let factor = FunctionBuilder::new(ETS.metric)
            .add_arg(Atom::var(SPENSO_TAG.chain_in))
            .add_arg(Atom::var(SPENSO_TAG.chain_out))
            .finish();
        let chain = FunctionBuilder::new(SPENSO_TAG.chain)
            .add_arg(representation.to_symbolic([]))
            .add_arg(representation.to_symbolic([]))
            .add_arg(factor)
            .finish();

        assert!(InterfaceInference::validate_placeholder_scope(chain.as_view(), false).is_ok());

        let bispinor: Representation<LibraryRep> =
            Bispinor {}.new_rep(Dimension::Concrete(4)).cast();
        let minkowski: Representation<LibraryRep> =
            Minkowski {}.new_rep(Dimension::Concrete(4)).cast();
        let gamma_factor = |input: Symbol, output: Option<Symbol>| {
            let mut factor = FunctionBuilder::new(AGS.gamma).add_arg(Atom::var(input));
            if let Some(output) = output {
                factor = factor.add_arg(Atom::var(output));
            }
            factor
                .add_arg(
                    minkowski
                        .slot::<AbstractIndex, _>(AbstractIndex::Normal(23))
                        .to_atom(),
                )
                .finish()
        };
        let gamma_chain =
            |factor| SPENSO_TAG.chain(bispinor.to_symbolic([]), bispinor.to_symbolic([]), [factor]);
        assert!(
            infer_interface(&gamma_chain(gamma_factor(
                SPENSO_TAG.chain_in,
                Some(SPENSO_TAG.chain_out),
            )))
            .is_ok()
        );
        assert!(infer_interface(&gamma_chain(gamma_factor(SPENSO_TAG.chain_in, None))).is_err());
        assert!(
            infer_interface(&gamma_chain(gamma_factor(
                SPENSO_TAG.chain_in,
                Some(SPENSO_TAG.chain_in),
            )))
            .is_err()
        );
    }

    #[test]
    fn placeholder_only_tensor_leaves_remain_structured_chain_factors() {
        Python::initialize();
        idenso::representations::initialize();
        let representation = ExtendibleReps::EUCLIDEAN.new_rep(Dimension::Concrete(4));
        let placeholder_leaf = |name| {
            FunctionBuilder::new(name)
                .add_arg(Atom::var(SPENSO_TAG.chain_in))
                .add_arg(Atom::var(SPENSO_TAG.chain_out))
                .finish()
        };
        let generic_name = SPENSO_TAG.tensor_symbol("placeholder_only_factor");
        let generic_with_key = FunctionBuilder::new(generic_name)
            .add_arg(Atom::num(7))
            .add_arg(Atom::var(SPENSO_TAG.chain_in))
            .add_arg(Atom::var(SPENSO_TAG.chain_out))
            .finish();

        for factor in [generic_with_key, placeholder_leaf(ETS.metric)] {
            let chain = SPENSO_TAG.chain(
                representation.to_symbolic([]),
                representation.to_symbolic([]),
                [factor.clone()],
            );
            assert_eq!(infer_interface(&chain).unwrap().canonical().order(), 2);

            let trace = SPENSO_TAG.trace(representation.to_symbolic([]), [factor]);
            assert!(infer_interface(&trace).unwrap().canonical().is_scalar());
        }

        let invalid_factors = [
            FunctionBuilder::new(generic_name)
                .add_arg(Atom::var(SPENSO_TAG.chain_in))
                .finish(),
            FunctionBuilder::new(generic_name)
                .add_arg(Atom::var(SPENSO_TAG.chain_in))
                .add_arg(Atom::var(SPENSO_TAG.chain_in))
                .finish(),
            FunctionBuilder::new(generic_name)
                .add_arg(Atom::var(SPENSO_TAG.chain_in))
                .add_arg(Atom::var(SPENSO_TAG.chain_out))
                .add_arg(Atom::num(7))
                .finish(),
        ];
        for factor in invalid_factors {
            let chain = SPENSO_TAG.chain(
                representation.to_symbolic([]),
                representation.to_symbolic([]),
                [factor],
            );
            assert!(infer_interface(&chain).is_err());
        }

        Python::attach(|py| -> PyResult<()> {
            let metric = TensorExpression::g(py, &SpensoRepresentation { representation }, None)?;
            let right = metric.clone_ref(py);
            let TensorDispatch::Expression(composed) = TensorExpression::compose(
                metric.bind(py).borrow(),
                py,
                right.bind(py).as_any(),
                (0, 1),
                (0, 1),
            )?
            else {
                panic!("composing symbolic metrics produced a network")
            };
            let reinferred = TensorExpression::from_atom_interface(
                py,
                composed.borrow(py).atom().clone(),
                None,
            )?;
            assert_eq!(
                reinferred.bind(py).borrow().interface().canonical().order(),
                2
            );
            Ok(())
        })
        .unwrap();
    }

    #[test]
    fn tensor_expression_boundary_rejects_nested_chain_scopes() {
        Python::initialize();
        Python::attach(|py| {
            let representation = ExtendibleReps::EUCLIDEAN.new_rep(Dimension::Concrete(4));
            let factor = FunctionBuilder::new(SPENSO_TAG.tensor_symbol("nested_python_factor"))
                .add_arg(Atom::var(SPENSO_TAG.chain_in))
                .add_arg(Atom::var(SPENSO_TAG.chain_out))
                .finish();
            let inner = SPENSO_TAG.trace(representation.to_symbolic([]), [factor]);
            let outer = SPENSO_TAG.chain(
                representation
                    .slot::<AbstractIndex, _>(AbstractIndex::Normal(41))
                    .to_atom(),
                representation
                    .slot::<AbstractIndex, _>(AbstractIndex::Normal(43))
                    .to_atom(),
                [inner],
            );
            let error = TensorExpression::from_atom_interface(
                py,
                outer,
                PartialStructure::from_logical_slots([]),
            )
            .expect_err("nested chain/trace syntax must be rejected");

            assert!(error.to_string().contains("cannot be nested"));
        });
    }

    #[test]
    fn indexing_root_chain_endpoints_equally_produces_a_root_trace() {
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let representation = ExtendibleReps::EUCLIDEAN.new_rep(Dimension::Concrete(4));
            let factor = FunctionBuilder::new(SPENSO_TAG.tensor_symbol("indexed_chain_factor"))
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
            let indexed = expression.bind(py).call1(("i", "i"))?;
            let indexed = indexed.extract::<PyRef<'_, TensorExpression>>()?;

            assert!(indexed.interface().canonical().is_scalar());
            assert!(
                matches!(indexed.atom().as_view(), AtomView::Fun(function) if function.get_symbol() == SPENSO_TAG.trace)
            );
            Ok(())
        })
        .unwrap();
    }

    #[test]
    fn empty_chain_with_equal_endpoints_is_a_trace() {
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let representation = ExtendibleReps::EUCLIDEAN.new_rep(Dimension::Concrete(4));
            let index = AbstractIndex::Normal(47);
            let factors = PyTuple::empty(py);
            let TensorDispatch::Expression(expression) = chain(
                py,
                SpensoSlot {
                    slot: representation.slot(index),
                },
                SpensoSlot {
                    slot: representation.slot(index),
                },
                &factors,
            )?
            else {
                panic!("an empty symbolic chain unexpectedly produced a network")
            };
            let expression = expression.bind(py).borrow();

            assert!(expression.interface().canonical().is_scalar());
            assert!(
                matches!(expression.atom().as_view(), AtomView::Fun(function) if function.get_symbol() == SPENSO_TAG.trace)
            );
            Ok(())
        })
        .unwrap();
    }

    #[test]
    fn index_tooling_reports_unresolved_composed_chain_errors() {
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let representation = ExtendibleReps::EUCLIDEAN.new_rep(Dimension::Concrete(4));
            let matrix = |name| {
                let ports = [
                    representation.slot(PartialIndex::open(0)),
                    representation.slot(PartialIndex::open(1)),
                ];
                let atom = FunctionBuilder::new(SPENSO_TAG.tensor_symbol(name))
                    .add_arg(representation.to_symbolic([]))
                    .add_arg(representation.to_symbolic([]))
                    .finish();
                SymbolicTensor::new(atom, PartialStructure::from_logical_slots(ports))
            };
            let composed = matrix("index_tooling_left")
                .multiply(&matrix("index_tooling_right"))
                .expect("compatible unresolved matrices should compose");
            assert!(matches!(
                composed.expression().as_view(),
                AtomView::Fun(function) if function.get_symbol() == SPENSO_TAG.chain
            ));
            let rust_list_error = composed
                .expression()
                .list_dangling::<AbstractIndex>()
                .expect_err("unresolved chain endpoints should return a parse error");
            assert!(matches!(
                rust_list_error,
                IndexToolingError::ListDangling { reason } if reason.contains("invalid chain start")
            ));
            let expression = TensorExpression::from_structured(py, composed)?;
            let list_error = TensorExpression::list_dangling(expression.bind(py).borrow())
                .err()
                .expect("unresolved chain endpoints should return a parse error");
            assert!(list_error.is_instance_of::<PyValueError>(py));
            assert!(
                list_error
                    .to_string()
                    .contains("cannot list dangling indices")
            );
            for dummies_only in [false, true] {
                let wrapped = TensorExpression::wrap_indices(
                    expression.borrow(py),
                    py,
                    symbolica::symbol!("index_tooling_wrap"),
                    dummies_only,
                )?;
                assert_eq!(wrapped.borrow(py).atom(), expression.borrow(py).atom());
                assert_eq!(
                    wrapped.borrow(py).interface(),
                    expression.borrow(py).interface()
                );
            }
            Ok(())
        })
        .unwrap();
    }

    #[test]
    fn network_tooling_failures_are_python_value_errors() {
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let malformed = FunctionBuilder::new(SPENSO_TAG.dot)
                .add_arg(Atom::var(symbolica::symbol!("malformed_dot_operand")))
                .finish();
            // Deliberately inject malformed internals to test tooling error translation.
            let expression = TensorExpression::from_parts_unchecked(
                py,
                malformed,
                PartialStructure::from_logical_slots(std::iter::empty()),
                None,
                Vec::new(),
            )?;

            for name in ["undo_dots", "dirac_adjoint"] {
                let error = expression
                    .bind(py)
                    .call_method0(name)
                    .expect_err("malformed dot notation should return an error");
                assert!(error.is_instance_of::<PyValueError>(py));
                assert!(error.to_string().contains("cannot parse tensor network"));
                assert!(error.to_string().contains("Invalid dot function"));
            }
            Ok(())
        })
        .unwrap();
    }

    #[test]
    fn dimension_rewrite_clears_descriptor_only_when_ports_contract() {
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let index = PartialIndex::Explicit(AbstractIndex::Normal(73331));
            let ports = [4, 6].map(|dimension| {
                ExtendibleReps::MINKOWSKI
                    .new_rep(Dimension::Concrete(dimension))
                    .slot(index)
            });
            let expression = FunctionBuilder::new(spenso::tensor_symbol!("dimension_named_tensor"))
                .add_args(ports.map(idenso::tensor::composition::port_atom))
                .finish();
            let name = symbol!("dimension_explicit_data_name");
            let original = TensorExpression::from_known_parts(
                py,
                expression,
                PartialStructure::from_logical_slots(ports),
                Some(name),
                vec![Atom::num(7)],
            )?;
            assert_eq!(original.borrow(py).interface().canonical().order(), 2);
            let result = TensorExpression::with_lorentz_dimension(
                original.borrow(py),
                py,
                ConvertibleToDimension(Dimension::Concrete(6)),
            )?;
            assert!(result.borrow(py).interface().canonical().is_scalar());
            assert!(result.borrow(py).name.is_none());
            assert!(result.borrow(py).name_args.is_empty());

            for dimension in [4, 5] {
                let result = TensorExpression::with_lorentz_dimension(
                    original.borrow(py),
                    py,
                    ConvertibleToDimension(Dimension::Concrete(dimension)),
                )?;
                assert_eq!(result.borrow(py).interface().canonical().order(), 2);
                assert_eq!(result.borrow(py).name, Some(name));
                assert_eq!(result.borrow(py).name_args, vec![Atom::num(7)]);
            }
            Ok(())
        })
        .unwrap();
    }

    #[test]
    fn bound_vector_operands_are_composite_descriptors() {
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let representation = ExtendibleReps::MINKOWSKI.new_rep(Dimension::Concrete(4));
            let slot = representation
                .slot::<AbstractIndex, _>(AbstractIndex::Normal(73211))
                .to_atom();
            let name = spenso::tensor_symbol!("bound_descriptor_tensor");
            let metadata = Atom::num(7);
            for vector in [
                spenso::vector_symbol!("bound_descriptor_p"),
                spenso::vector_symbol!("bound_descriptor_q"),
            ] {
                let bound = FunctionBuilder::new(vector)
                    .add_arg(representation.to_symbolic([]))
                    .finish();
                for arguments in [[&bound, &slot], [&slot, &bound]] {
                    let atom = FunctionBuilder::new(name)
                        .add_arg(&metadata)
                        .add_args(arguments)
                        .finish();
                    let tensor = TensorExpression::from_atom_interface(py, atom.clone(), None)?;
                    let tensor = tensor.borrow(py);
                    assert_eq!(*tensor.atom(), atom);
                    assert_eq!(tensor.interface().canonical().order(), 1);
                    assert!(tensor.name.is_none());
                    assert!(tensor.name_args.is_empty());
                }
            }
            let atom = FunctionBuilder::new(name)
                .add_arg(&metadata)
                .add_arg(&slot)
                .finish();
            let tensor = TensorExpression::from_atom_interface(py, atom, None)?;
            let tensor = tensor.borrow(py);
            assert_eq!(tensor.name, Some(name));
            assert_eq!(tensor.name_args, [metadata]);
            Ok(())
        })
        .unwrap();
    }

    #[test]
    fn compact_metric_retains_builtin_descriptor() {
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let representation = ExtendibleReps::MINKOWSKI.new_rep(Dimension::Concrete(4));
            let vectors = [
                spenso::vector_symbol!("builtin_descriptor_p"),
                spenso::vector_symbol!("builtin_descriptor_q"),
            ]
            .map(|name| {
                FunctionBuilder::new(name)
                    .add_arg(representation.to_symbolic([]))
                    .finish()
            });
            let atom = FunctionBuilder::new(ETS.metric).add_args(&vectors).finish();
            let tensor = TensorExpression::from_atom_interface(py, atom.clone(), None)?;
            let tensor = tensor.borrow(py);
            assert_eq!(*tensor.atom(), atom);
            assert!(tensor.interface().canonical().is_scalar());
            assert_eq!(tensor.name, Some(ETS.metric));
            assert_eq!(tensor.name_args, vectors);
            Ok(())
        })
        .unwrap();
    }

    #[test]
    fn transformed_diagonal_leaf_clears_descriptor_only_when_raw_ports_contract() {
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let x = Atom::var(symbol!("diagonal_descriptor_x"));
            let name = symbol!("diagonal_descriptor_data");
            let original = TensorExpression::from_known_parts(
                py,
                &x + Atom::one(),
                PartialStructure::from_logical_slots([]),
                Some(name),
                vec![Atom::num(7)],
            )?;
            let representation = ExtendibleReps::EUCLIDEAN.new_rep(Dimension::Concrete(3));
            let slot = representation
                .slot::<AbstractIndex, _>(AbstractIndex::Normal(73201))
                .to_atom();
            let diagonal =
                FunctionBuilder::new(spenso::tensor_symbol!("diagonal_descriptor_tensor"))
                    .add_arg(&slot)
                    .add_arg(&slot)
                    .finish();
            let raw = InterfaceInference::default()
                .infer_validated(diagonal.as_view())
                .unwrap();
            assert_eq!(
                raw.canonical().order(),
                2,
                "the fixture must expose an intermediate port pair"
            );
            let result =
                TensorExpression::from_transformed_atom(&original.borrow(py), py, diagonal)?;
            let result = result.borrow(py);
            assert_eq!(result.interface().canonical().order(), 0);
            assert!(result.name.is_none());
            assert!(result.name_args.is_empty());

            // A scalar rewrite with no port contraction retains explicitly assigned identity.
            let scalar = TensorExpression::from_transformed_atom(
                &original.borrow(py),
                py,
                x + Atom::num(2),
            )?;
            let scalar = scalar.borrow(py);
            assert_eq!(scalar.name, Some(name));
            assert_eq!(scalar.name_args, vec![Atom::num(7)]);
            Ok(())
        })
        .unwrap();
    }
}
