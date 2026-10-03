use std::{collections::HashMap, sync::Arc};

use idenso::{
    IndexTooling,
    color::CS,
    dirac::{AGS, spinor_matrix_structure},
    epsilon::EPSILON_SYMBOL,
    tensor::{
        SymbolicTensor, composition,
        inference::{InterfaceInference, TensorInferenceError},
    },
};
use pyo3::{
    IntoPyObjectExt,
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
        ConvertibleToExpression, ConvertibleToPatternRestriction, ConvertibleToReplaceWith,
        FormattedOutputBackend, PythonExpression, PythonFormattedOutput, PythonReplacement,
        PythonTransformer,
    },
    atom::{Atom, AtomCore, AtomView, Symbol},
};

use crate::{
    ModuleInit, SliceOrIntOrExpanded, Spensor, display,
    library::SpensorLibrary,
    network::{ConvertibleToSpensoNet, ExecutionMode, SpensoNet},
    simplification::{
        CanonicalizationError, CookingError, DiracAdjointError, Intern, NetworkToolingError,
    },
    structure::{
        ArithmeticStructure, ConvertibleToAbstractIndex, ConvertibleToDimension,
        ConvertibleToSpensoName, SpensoIndices, SpensoName, SpensoRepresentation, SpensoSlot,
        SpensoSlotOrArgOrRep,
    },
};

/// Sentinel type for an unresolved tensor axis.
///
/// Use the exported singleton ``AUTO`` (also available as ``_``), rather than
/// constructing this type. During index() or reindex(), AUTO leaves that axis
/// unchanged. It is not a wildcard and not a component coordinate.
/// Unresolved axes display as hollow squares. When several axes of a tensor
/// are unresolved, their squares carry zero-based axis numbers, excluding
/// scalar parameters. These local numbers do not imply contractions. Hover
/// over a square in HTML output to see its axis and representation.
///
/// Examples
/// --------
/// >>> from symbolica.community.tensor import Representation, TensorName, TensorExpression
/// >>> space = Representation.euc(2)
/// >>> A = TensorName("M")(space, space)
/// >>> from symbolica.community.tensor import AUTO
/// >>> partly_indexed = A(AUTO, "j")
/// >>> partly_indexed.rank
/// 2
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

/// A symbolic tensor expression with an ordered set of external axes.
///
/// Construct tensors from TensorName and Representation objects, or wrap a
/// Symbolica expression whose tensor structure can be inferred. Scalar tensors
/// have rank zero. Components are obtained separately through ``to_tensor``.
///
/// Products contract matching explicit indices and unambiguous compatible
/// unresolved axes. Use ``outer`` to keep axes independent, ``contract_ports``
/// to choose a contraction, or ``compose`` to choose matrix channels.
/// Arithmetic with Tensor or TensorNetwork operands returns a TensorNetwork.
///
/// Reductions preserve this tensor type and keep unrelated sums factorized.
/// Ordinary algebra, including ``expand_num()``, ``collect_factors()``, and
/// ``factor()``, also returns a tensor expression with its ordered interface.
/// This is a Symbolica ``Expression`` subclass: inspection, matching, evaluators,
/// and other unoverridden methods are inherited directly. Inherited transformations
/// return ordinary expressions; the tensor-aware overrides retain the interface.
/// Use ``to_expression()`` to explicitly drop the tensor interface.
///
/// Pass ``intern="indices"`` to reversibly encode compound index labels, or
/// ``intern="flattened"`` to create readable names. Both modes operate only on
/// index payloads. An explicit ``structure`` retains the declared axis order
/// and can describe the interface of a symbolic zero.
///
/// Examples
/// --------
/// >>> from symbolica.community.tensor import Representation, TensorName, TensorExpression
/// >>> space = Representation.euc(3)
/// >>> A = TensorName("A")(space, space)
/// >>> indexed = A("i", "j")
/// >>> indexed.rank
/// 2
#[cfg_attr(feature = "python_stubgen", pyo3_stub_gen::derive::gen_stub_pyclass)]
#[pyclass(
    frozen,
    from_py_object,
    extends = PythonExpression,
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
            TensorInferenceError::InvalidBuiltinSignature { factory, signature } => {
                let constructor = match factory {
                    "gamma" => "dirac_gamma",
                    "f" => "color_f",
                    "t" => "color_t",
                    name => name,
                };
                format!(
                    "predefined tensor `{factory}` requires {signature}; use TensorExpression.{constructor}(...) to construct it"
                )
            }
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
            }
            .into_initializer(),
        )
    }

    /// Initialize the inherited immutable expression from the shared tensor owner.
    /// All typed results use this boundary, keeping inherited Symbolica methods
    /// on exactly the same atom without re-inferring the tensor interface.
    pub(crate) fn into_initializer(self) -> PyClassInitializer<Self> {
        let expression = PythonExpression {
            expr: self.atom().clone(),
        };
        (self, expression).into()
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

    /// Reuse the shared proof for algebra that treats tensor leaves as opaque.
    fn from_algebra_atom(
        self_: &PyRef<'_, Self>,
        py: Python<'_>,
        atom: Atom,
    ) -> PyResult<Py<Self>> {
        if atom == self_.atom() {
            return Py::new(py, Self::clone(self_).into_initializer());
        }
        let value = self_
            .value
            .with_algebra_result(atom)
            .map_err(Self::inference_error)?;
        let (name, args) = Self::transformed_descriptor(self_, value.expression());
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
            ArithmeticStructure::Tensor(Py::new(py, Self::clone(self_).into_initializer())?),
            None,
        )
    }

    pub fn structured(&self) -> &SymbolicTensor<PartialStructure> {
        &self.value
    }

    pub(crate) fn index_replacements(
        interface: &PartialStructure,
        indices: &Bound<'_, PyTuple>,
        intern: Option<Intern>,
    ) -> PyResult<HashMap<usize, AbstractIndex>> {
        let open_positions = interface.open_positions();
        Self::port_replacements(interface, &open_positions, indices, intern)
    }

    pub(crate) fn port_replacements(
        interface: &PartialStructure,
        positions: &[usize],
        indices: &Bound<'_, PyTuple>,
        intern: Option<Intern>,
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
                index_value(argument.extract()?, intern)?
            };
            replacements.insert(position, index);
        }
        Ok(replacements)
    }

    pub(crate) fn named_replacements(
        interface: &PartialStructure,
        mapping: &Bound<'_, PyDict>,
        intern: Option<Intern>,
    ) -> PyResult<HashMap<usize, AbstractIndex>> {
        let slots = interface.logical_slots();
        let mut result = HashMap::new();
        for (key, value) in mapping.iter() {
            let typed = key.extract::<SpensoSlot>().ok();
            let index = if let Some(slot) = &typed {
                slot.slot.aind()
            } else {
                index_value(key.extract()?, intern)?
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
                Self::port_replacements(interface, &positions, &arguments, intern)?
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
    value.is_instance_of::<SpensoNet>()
        || value.is_instance_of::<Spensor>()
        || value
            .extract::<crate::projectors::PyFactorProjector>()
            .is_ok_and(|group| group.network.is_some())
}

fn concrete_network(value: &Bound<'_, PyAny>) -> PyResult<Option<SpensoNet>> {
    if let Ok(network) = value.extract::<SpensoNet>() {
        Ok(Some(network))
    } else if let Ok(tensor) = value.extract::<Spensor>() {
        SpensoNet::from_tensor(tensor).map(Some)
    } else if let Ok(group) = value.extract::<crate::projectors::PyFactorProjector>() {
        Ok(group.network)
    } else {
        Ok(None)
    }
}

impl TensorOperand {
    pub(crate) fn into_structured(self) -> PyResult<SymbolicTensor<PartialStructure>> {
        match self {
            Self::Structured(value) => Ok(value),
            Self::Scalar(atom) => {
                SymbolicTensor::checked_parts(atom, PartialStructure::from_logical_slots([]))
                    .map_err(TensorExpression::inference_error)
            }
        }
    }

    pub(crate) fn extract(value: &Bound<'_, PyAny>) -> PyResult<Self> {
        if let Ok(group) = value.extract::<crate::projectors::PyFactorProjector>() {
            if group.network.is_none() {
                return Ok(Self::Structured(group.structure));
            }
            return Err(PyTypeError::new_err(
                "component factor projectors require tensor-network dispatch",
            ));
        }
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
    intern: Option<Intern>,
) -> PyResult<AbstractIndex> {
    match value {
        ConvertibleToAbstractIndex::Aind(index) => Ok(index),
        ConvertibleToAbstractIndex::Atom(expression) => {
            let atom = match (intern, expression.expr.as_view()) {
                (Some(settings), AtomView::Fun(fun)) => {
                    settings.rust().cook_function(fun).map_err(|error| {
                        CookingError::new_err(format!("cannot cook index: {error:?}"))
                    })?
                }
                _ => expression.expr,
            };
            atom.as_view().try_into().map_err(|error| {
                let hint = if intern.is_some() {
                    ""
                } else {
                    " Try passing intern=\"indices\"."
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
    /// Wrap symbolic algebra as a tensor and check its external axes.
    ///
    /// Parameters
    /// ----------
    /// expression : TensorExpression or scalar expression
    ///     Existing tensor, Symbolica Expression, or scalar convertible to one.
    ///     Tensor syntax is inferred; an ordinary scalar has no external axes.
    /// structure : TensorStructure, optional
    ///     Canonical free-axis signature to check against the expression. It does
    ///     not reorder arguments or component axes. For zero, it supplies the rank
    ///     and axes that cannot be inferred from the scalar expression alone.
    /// intern : {"indices", "flattened"} or None, optional
    ///     Encode compound index labels before checking the tensor structure.
    ///     "indices" retains reversible payloads; "flattened" creates readable names.
    ///     Both preserve index tags and leave scalar arguments and tensor heads intact.
    ///     None leaves labels unchanged and requires ordinary atomic index labels.
    ///
    /// Returns
    /// -------
    /// TensorExpression
    ///     A checked tensor expression, preserving existing metadata where possible.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import TensorExpression, TensorStructure, Representation
    /// >>> zero = TensorExpression(0, structure=TensorStructure([Representation.euc(3)]))
    /// >>> zero.rank
    /// 1
    /// >>> from symbolica import S
    /// >>> from symbolica.community.tensor import TensorName
    /// >>> vector = TensorName.vector("v")(Representation.euc(3))
    /// >>> label, payload = S("i"), S("edge")(7, 2)
    /// >>> raw = vector(label).to_expression().replace(label, payload)
    /// >>> TensorExpression(raw, intern="indices").rank
    /// 1
    #[new]
    #[pyo3(signature = (expression, *, structure = None, intern = None))]
    #[gen_stub(override_return_type(type_repr = "TensorExpression"))]
    fn new(
        py: Python<'_>,
        expression: &Bound<'_, PyAny>,
        structure: Option<&crate::metadata::SpensoTensorStructure>,
        intern: Option<Intern>,
    ) -> PyResult<PyClassInitializer<Self>> {
        if let Some(structure) = structure {
            let atom = if let Ok(value) = expression.extract::<PyRef<'_, Self>>() {
                value.atom().clone()
            } else {
                expression
                    .extract::<ConvertibleToExpression>()?
                    .to_expression()
                    .expr
            };
            let atom = if let Some(settings) = intern {
                settings
                    .rust()
                    .try_cook_indices(atom.as_view())
                    .map_err(|error| {
                        CookingError::new_err(format!("cannot cook indices: {error:?}"))
                    })?
            } else {
                atom
            };
            let zero_interface = atom
                .as_view()
                .is_zero()
                .then(|| Canonicalized::identity(structure.interface.clone()));
            let checked = Self::from_atom_interface(py, atom, zero_interface)?;
            if checked.borrow(py).structure() != *structure {
                return Err(PyValueError::new_err(
                    "supplied structure does not match the expression's external interface",
                ));
            }
            return Ok(checked.borrow(py).clone().into_initializer());
        }
        let tensor = if let Some(settings) = intern {
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
        Ok(tensor.clone().into_initializer())
    }

    /// Construct the metric pairing, optionally with explicit indices.
    ///
    /// Parameters
    /// ----------
    /// rep : Representation or Slot
    ///     First axis. A representation leaves it unresolved; a slot supplies
    ///     both its representation and its index.
    /// other : Representation or Slot, optional
    ///     Second axis. Defaults to an unresolved axis in the first axis's
    ///     representation, even when the first argument is a slot.
    ///
    /// Returns
    /// -------
    /// TensorExpression
    ///     The metric with supplied indices in argument order. Representations
    ///     remain unresolved for later indexing. Matching compatible indices
    ///     contract, just as when indexing an existing metric.
    ///
    /// Raises
    /// ------
    /// ValueError
    ///     If the two spaces are neither equal nor dual partners, including
    ///     when their dimensions differ.
    /// TypeError
    ///     If an argument is not a Representation or Slot.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import TensorExpression, Representation
    /// >>> r = Representation.mink(4)
    /// >>> metric = TensorExpression.g(r("mu"), r("nu"))
    /// >>> metric == TensorExpression.g(r)("mu", "nu")
    /// True
    /// >>> partly_indexed = TensorExpression.g(r("mu"), r)
    /// >>> partly_indexed("nu") == metric
    /// True
    /// >>> f = Representation.cof(3)
    /// >>> identity = TensorExpression.g(f("i"), f.dual()("j"))
    /// >>> identity.rank
    /// 2
    #[staticmethod]
    #[pyo3(signature = (rep, other = None))]
    fn g(
        py: Python<'_>,
        #[gen_stub(override_type(type_repr = "Representation | Slot"))] rep: SpensoSlotOrArgOrRep,
        #[gen_stub(override_type(type_repr = "Representation | Slot | None"))] other: Option<
            SpensoSlotOrArgOrRep,
        >,
    ) -> PyResult<Py<TensorExpression>> {
        let port = |argument, position| match argument {
            SpensoSlotOrArgOrRep::Slot(slot) => Ok(slot
                .slot
                .rep()
                .slot(PartialIndex::Explicit(slot.slot.aind()))),
            SpensoSlotOrArgOrRep::Rep(rep) => {
                Ok(rep.representation.slot(PartialIndex::open(position)))
            }
            SpensoSlotOrArgOrRep::Arg(_) => Err(PyTypeError::new_err(
                "metric ports must be Representation or Slot objects",
            )),
        };
        let left = port(rep, 0)?;
        let right = other
            .map(|other| port(other, 1))
            .transpose()?
            .unwrap_or_else(|| left.rep().slot(PartialIndex::open(1)));
        if left.rep() != right.rep() && left.rep().dual() != right.rep() {
            return Err(PyValueError::new_err(format!(
                "metric ports must belong to the same representation space (possibly dual); got {} and {}",
                left.rep().to_symbolic([]),
                right.rep().to_symbolic([]),
            )));
        }
        let atom = ETS.metric(composition::port_atom(left), composition::port_atom(right));
        Self::from_known_parts(
            py,
            atom,
            PartialStructure::from_logical_slots([left, right]),
            Some(ETS.metric),
            Vec::new(),
        )
    }

    /// Construct the metric map for raising or lowering an index.
    ///
    /// Parameters
    /// ----------
    /// rep : Representation
    ///     Space of both axes.
    ///
    /// Returns
    /// -------
    /// TensorExpression
    ///     A rank-2 symbolic tensor with distinct unresolved axes. Call the
    ///     result with index labels to fill those axes in the order described below.
    ///
    /// Notes
    /// -----
    /// This map accounts for the signs in the representation metric.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import TensorExpression, Representation
    /// >>> tensor = TensorExpression.flat(Representation.mink(4))
    /// >>> tensor.rank
    /// 2
    #[staticmethod]
    fn flat(py: Python<'_>, rep: &SpensoRepresentation) -> PyResult<Py<TensorExpression>> {
        let structure = ExplicitKey::<AbstractIndex>::from_iter(
            [rep.representation, rep.representation],
            ETS.flat,
            None,
        );
        Self::from_structure(py, &structure)
    }

    /// Construct the Dirac gamma matrices for Clifford algebra.
    ///
    /// Parameters
    /// ----------
    /// minkowski_dimension : int, Expression, or str
    ///     Dimension of the Minkowski vector index, for example 4 or a symbol D.
    ///
    /// Returns
    /// -------
    /// TensorExpression
    ///     A rank-3 symbolic tensor with distinct unresolved axes. Call the
    ///     result with index labels to fill those axes in the order described below.
    ///
    /// Notes
    /// -----
    /// The axes are (spinor row, spinor column, Lorentz index), with four
    /// components on each spinor axis. The anticommutator obeys
    /// {gamma^mu, gamma^nu} = 2 g^{mu,nu} I. This helper is unrelated to
    /// the scalar gamma special function. Explicit library matrices are four-dimensional.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import TensorExpression, Representation
    /// >>> tensor = TensorExpression.dirac_gamma(4)
    /// >>> tensor.rank
    /// 3
    #[staticmethod]
    fn dirac_gamma(
        py: Python<'_>,
        minkowski_dimension: ConvertibleToDimension,
    ) -> PyResult<Py<TensorExpression>> {
        Self::from_structure(py, &AGS.gamma_strct::<AbstractIndex>(minkowski_dimension.0))
    }

    /// Construct the time-component Dirac matrix gamma^0.
    ///
    /// Parameters
    /// ----------
    /// spinor_dimension : int, Expression, or str
    ///     Number of spinor components.
    ///
    /// Returns
    /// -------
    /// TensorExpression
    ///     A rank-2 symbolic tensor with distinct unresolved axes. Call the
    ///     result with index labels to fill those axes in the order described below.
    ///
    /// Notes
    /// -----
    /// The two axes are spinor row and column. The HEP libraries provide
    /// four-component matrices in the Weyl basis.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import TensorExpression, Representation
    /// >>> tensor = TensorExpression.gamma0(4)
    /// >>> tensor.rank
    /// 2
    #[staticmethod]
    fn gamma0(
        py: Python<'_>,
        spinor_dimension: ConvertibleToDimension,
    ) -> PyResult<Py<TensorExpression>> {
        Self::from_structure(py, &AGS.gamma0_strct::<AbstractIndex>(spinor_dimension.0))
    }

    /// Construct the Dirac charge-conjugation matrix.
    ///
    /// Parameters
    /// ----------
    /// spinor_dimension : int, Expression, or str
    ///     Number of spinor components.
    ///
    /// Returns
    /// -------
    /// TensorExpression
    ///     A rank-2 symbolic tensor with distinct unresolved axes. Call the
    ///     result with index labels to fill those axes in the order described below.
    ///
    /// Notes
    /// -----
    /// The axes are spinor row and column. The head is antisymmetric.
    /// The HEP libraries provide the four-component Weyl-basis convention.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import TensorExpression, Representation
    /// >>> tensor = TensorExpression.charge_conjugation(4)
    /// >>> tensor.rank
    /// 2
    #[staticmethod]
    fn charge_conjugation(
        py: Python<'_>,
        spinor_dimension: ConvertibleToDimension,
    ) -> PyResult<Py<TensorExpression>> {
        Self::from_structure(
            py,
            &spinor_matrix_structure::<AbstractIndex>(AGS.charge_conjugation, spinor_dimension.0),
        )
    }

    /// Construct the totally antisymmetric Levi-Civita tensor.
    ///
    /// Parameters
    /// ----------
    /// rep : Representation
    ///     Space of every axis.
    /// rank : int, default 4
    ///     Number of axes, independent of the representation dimension.
    ///
    /// Returns
    /// -------
    /// TensorExpression
    ///     A symbolic tensor with ``rank`` distinct unresolved axes.
    ///
    /// Notes
    /// -----
    /// This constructs symbolic epsilon syntax. It does not populate a component
    /// library; rank and space dimension need not be equal for symbolic work.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import TensorExpression, Representation
    /// >>> tensor = TensorExpression.levi_civita(Representation.euc(3), rank=3)
    /// >>> tensor.rank
    /// 3
    #[staticmethod]
    #[pyo3(signature = (rep, rank=4))]
    fn levi_civita(
        py: Python<'_>,
        rep: &SpensoRepresentation,
        rank: usize,
    ) -> PyResult<Py<TensorExpression>> {
        let structure = ExplicitKey::<AbstractIndex>::from_iter(
            std::iter::repeat_n(rep.representation, rank),
            *EPSILON_SYMBOL,
            None,
        );
        Self::from_structure(py, &structure)
    }

    /// Expand normalized factor groups into their permutation sums.
    ///
    /// Returns
    /// -------
    /// TensorExpression
    ///     Tensor algebra with symmetric, antisymmetric, and cyclic groups
    ///     expanded. Other sums and products remain factored.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import Representation, TensorName, TensorExpression
    /// >>> space = Representation.euc(2)
    /// >>> A = TensorName("M")(space, space)
    /// >>> from symbolica.community.tensor import FactorProjector, trace
    /// >>> expression = trace(space, FactorProjector.symmetric(A, A))
    /// >>> expression.expand_projectors().is_scalar
    /// True
    fn expand_projectors(self_: PyRef<'_, Self>, py: Python<'_>) -> PyResult<Py<Self>> {
        let value = Self::structured(&self_)
            .expand_projectors()
            .map_err(|error| PyValueError::new_err(error.to_string()))?;
        Self::from_structured(py, value)
    }

    /// Construct the Dirac chirality matrix gamma^5.
    ///
    /// Parameters
    /// ----------
    /// spinor_dimension : int, Expression, or str
    ///     Number of spinor components.
    ///
    /// Returns
    /// -------
    /// TensorExpression
    ///     A rank-2 symbolic tensor with distinct unresolved axes. Call the
    ///     result with index labels to fill those axes in the order described below.
    ///
    /// Notes
    /// -----
    /// The axes are spinor row and column. The HEP libraries provide the
    /// four-component Weyl-basis chirality matrix.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import TensorExpression, Representation
    /// >>> tensor = TensorExpression.gamma5(4)
    /// >>> tensor.rank
    /// 2
    #[staticmethod]
    fn gamma5(
        py: Python<'_>,
        spinor_dimension: ConvertibleToDimension,
    ) -> PyResult<Py<TensorExpression>> {
        Self::from_structure(py, &AGS.gamma5_strct::<AbstractIndex>(spinor_dimension.0))
    }

    /// Construct the left-chiral Dirac projector (I - gamma^5)/2.
    ///
    /// Parameters
    /// ----------
    /// spinor_dimension : int, Expression, or str
    ///     Number of spinor components.
    ///
    /// Returns
    /// -------
    /// TensorExpression
    ///     A rank-2 symbolic tensor with distinct unresolved axes. Call the
    ///     result with index labels to fill those axes in the order described below.
    ///
    /// Notes
    /// -----
    /// The two axes are spinor row and column. This is a projector on spinors,
    /// not a FactorProjector for permuting tensor factors.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import TensorExpression, Representation
    /// >>> tensor = TensorExpression.projm(4)
    /// >>> tensor.rank
    /// 2
    #[staticmethod]
    fn projm(
        py: Python<'_>,
        spinor_dimension: ConvertibleToDimension,
    ) -> PyResult<Py<TensorExpression>> {
        Self::from_structure(py, &AGS.projm_strct::<AbstractIndex>(spinor_dimension.0))
    }

    /// Construct the right-chiral Dirac projector (I + gamma^5)/2.
    ///
    /// Parameters
    /// ----------
    /// spinor_dimension : int, Expression, or str
    ///     Number of spinor components.
    ///
    /// Returns
    /// -------
    /// TensorExpression
    ///     A rank-2 symbolic tensor with distinct unresolved axes. Call the
    ///     result with index labels to fill those axes in the order described below.
    ///
    /// Notes
    /// -----
    /// The two axes are spinor row and column. This is a projector on spinors,
    /// not a FactorProjector for permuting tensor factors.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import TensorExpression, Representation
    /// >>> tensor = TensorExpression.projp(4)
    /// >>> tensor.rank
    /// 2
    #[staticmethod]
    fn projp(
        py: Python<'_>,
        spinor_dimension: ConvertibleToDimension,
    ) -> PyResult<Py<TensorExpression>> {
        Self::from_structure(py, &AGS.projp_strct::<AbstractIndex>(spinor_dimension.0))
    }

    /// Construct the antisymmetric Dirac sigma tensor.
    ///
    /// Parameters
    /// ----------
    /// minkowski_dimension : int, Expression, or str
    ///     Dimension of each Minkowski vector index.
    ///
    /// Returns
    /// -------
    /// TensorExpression
    ///     A rank-4 symbolic tensor with distinct unresolved axes. Call the
    ///     result with index labels to fill those axes in the order described below.
    ///
    /// Notes
    /// -----
    /// The logical axes are (mu, nu, spinor row, spinor column). Both spinor
    /// axes have dimension four.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import TensorExpression, Representation
    /// >>> tensor = TensorExpression.sigma(4)
    /// >>> tensor.rank
    /// 4
    #[staticmethod]
    fn sigma(
        py: Python<'_>,
        minkowski_dimension: ConvertibleToDimension,
    ) -> PyResult<Py<TensorExpression>> {
        Self::from_structure(py, &AGS.sigma_strct::<AbstractIndex>(minkowski_dimension.0))
    }

    /// Construct the antisymmetric color structure constants f^{abc}.
    ///
    /// Parameters
    /// ----------
    /// adjoint_dimension : int, Expression, or str
    ///     Dimension of the adjoint representation, for example 8 for SU(3).
    ///
    /// Returns
    /// -------
    /// TensorExpression
    ///     A rank-3 symbolic tensor with distinct unresolved axes. Call the
    ///     result with index labels to fill those axes in the order described below.
    ///
    /// Notes
    /// -----
    /// The axes are three adjoint indices (a, b, c). This dimension is the
    /// number of generators, not the dimension of the fundamental representation.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import TensorExpression, Representation
    /// >>> tensor = TensorExpression.color_f(8)
    /// >>> tensor.rank
    /// 3
    #[staticmethod]
    fn color_f(
        py: Python<'_>,
        adjoint_dimension: ConvertibleToDimension,
    ) -> PyResult<Py<TensorExpression>> {
        Self::from_structure(py, &CS.f_strct::<AbstractIndex>(adjoint_dimension.0))
    }

    /// Construct the fundamental color generators T^a.
    ///
    /// Parameters
    /// ----------
    /// adjoint_dimension : int, Expression, or str
    ///     Adjoint dimension, for example 8 for SU(3).
    /// fundamental_dimension : int, Expression, or str
    ///     Fundamental dimension, for example 3 for SU(3).
    ///
    /// Returns
    /// -------
    /// TensorExpression
    ///     A rank-3 symbolic tensor with distinct unresolved axes. Call the
    ///     result with index labels to fill those axes in the order described below.
    ///
    /// Notes
    /// -----
    /// The axes are (adjoint, fundamental, antifundamental). The HEP libraries
    /// provide SU(3) generators; symbolic color identities retain representation
    /// invariants unless their substitution is requested.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import TensorExpression, Representation
    /// >>> tensor = TensorExpression.color_t(8, 3)
    /// >>> tensor.rank
    /// 3
    #[staticmethod]
    fn color_t(
        py: Python<'_>,
        adjoint_dimension: ConvertibleToDimension,
        fundamental_dimension: ConvertibleToDimension,
    ) -> PyResult<Py<TensorExpression>> {
        Self::from_structure(
            py,
            &CS.t_strct::<AbstractIndex>(fundamental_dimension.0, adjoint_dimension.0),
        )
    }

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
    /// >>> space = Representation.euc(3)
    /// >>> A = TensorName("A")(space, space)
    /// >>> A.rank
    /// 2
    #[getter]
    fn rank(&self) -> usize {
        self.interface().canonical().order()
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
    /// >>> space = Representation.euc(3)
    /// >>> A = TensorName("A")(space, space)
    /// >>> A.is_scalar
    /// False
    #[getter]
    fn is_scalar(&self) -> bool {
        self.interface().canonical().is_scalar()
    }

    /// Dimensions in logical axis order.
    ///
    /// Returns
    /// -------
    /// tuple of int or Expression
    ///     Concrete dimensions are integers; symbolic dimensions remain expressions.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import Representation, TensorName, TensorExpression
    /// >>> space = Representation.euc(3)
    /// >>> A = TensorName("A")(space, space)
    /// >>> A.shape
    /// (3, 3)
    #[getter]
    #[gen_stub(override_return_type(type_repr = "tuple[int | Expression, ...]"))]
    fn shape(&self, py: Python<'_>) -> PyResult<Py<PyTuple>> {
        crate::metadata::axis_shape(py, self.interface().logical_slots())
    }

    #[doc = python_doc!("Tensor.axes")]
    #[getter]
    #[gen_stub(override_return_type(type_repr = "tuple[Slot | Representation, ...]"))]
    fn axes(&self, py: Python<'_>) -> PyResult<Py<PyTuple>> {
        crate::metadata::axis_objects(py, self.interface().logical_slots())
    }

    /// Canonical external signature of the whole expression.
    ///
    /// Returns
    /// -------
    /// TensorStructure
    ///     Immutable signature in canonical order. Inspect this to determine
    ///     rank, representations, or existing index labels.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import Representation, TensorName, TensorExpression
    /// >>> space = Representation.euc(3)
    /// >>> A = TensorName("A")(space, space)
    /// >>> A.structure.shape
    /// (3, 3)
    #[getter]
    fn structure(&self) -> crate::metadata::SpensoTensorStructure {
        crate::metadata::SpensoTensorStructure::from_interface(self.interface())
    }

    /// Optional name identifying stored tensor data.
    ///
    /// Returns
    /// -------
    /// TensorName or None
    ///     Composite expressions may have no name.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import Representation, TensorName, TensorExpression
    /// >>> space = Representation.euc(3)
    /// >>> A = TensorName("A")(space, space)
    /// >>> name = A.name
    #[getter(name)]
    fn data_name(&self) -> Option<SpensoName> {
        self.name.map(|name| SpensoName { name })
    }

    /// Scalar arguments identifying the named tensor.
    ///
    /// Returns
    /// -------
    /// tuple of Expression
    ///     Scalar key arguments in their original order; empty for an unnamed
    ///     expression. These belong to the expression, not its external signature.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica import S
    /// >>> from symbolica.community.tensor import Representation, TensorName
    /// >>> x = S("x")
    /// >>> A = TensorName("A")(x, 7, Representation.euc(3))
    /// >>> A.arguments == (x, 7)
    /// True
    #[getter]
    #[gen_stub(override_return_type(type_repr = "tuple[Expression, ...]"))]
    fn arguments(&self, py: Python<'_>) -> PyResult<Py<PyTuple>> {
        Ok(PyTuple::new(
            py,
            self.name_args.iter().cloned().map(PythonExpression::from),
        )?
        .unbind())
    }

    /// Assign a data identity without rewriting the symbolic computation.
    ///
    /// Parameters
    /// ----------
    /// name : TensorName, str, Expression, or TensorExpression
    ///     New tensor name. A rank-zero atomic TensorExpression may also supply
    ///     scalar key arguments, for example TensorName("B")(7).
    ///
    /// Returns
    /// -------
    /// TensorExpression
    ///     A new value with the requested identity and the same axes.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import Representation, TensorName, TensorExpression
    /// >>> space = Representation.euc(2)
    /// >>> A = TensorName("M")(space, space)
    /// >>> named = A.with_name("named_matrix")
    /// >>> named.rank
    /// 2
    fn with_name(
        self_: PyRef<'_, Self>,
        py: Python<'_>,
        name: ConvertibleToSpensoName,
    ) -> PyResult<Py<Self>> {
        let ConvertibleToSpensoName(name, args) = name;
        Self::from_shared(py, Arc::clone(&self_.value), Some(name.name), args)
    }

    /// Count all components in the tensor shape.
    ///
    /// Returns
    /// -------
    /// int
    ///     Product of the axis dimensions, including implicit zeros; dimensions
    ///     must be concrete.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import Representation, TensorName, TensorExpression
    /// >>> space = Representation.euc(2)
    /// >>> A = TensorName("M")(space, space)
    /// >>> len(A)
    /// 4
    fn __len__(&self) -> PyResult<usize> {
        logical_size(&logical_dimensions(self.interface())?)
    }

    #[doc = python_doc!("TensorExpression.__getitem__")]
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

    /// Return ordinary Symbolica algebra without the tensor interface.
    ///
    /// Returns
    /// -------
    /// Expression
    ///     The underlying expression. Axis metadata is no longer carried by
    ///     the result; operations on it are ordinary Symbolica operations.
    ///
    /// Notes
    /// -----
    /// TensorExpression already inherits Expression. Convert when deliberately
    /// requesting ordinary Symbolica behavior, including replacements that change
    /// the tensor interface, rather than its tensor-aware overrides.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import Representation, TensorName, TensorExpression
    /// >>> space = Representation.euc(2)
    /// >>> A = TensorName("M")(space, space)
    /// >>> raw = A.to_expression()
    fn to_expression(&self) -> PythonExpression {
        PythonExpression {
            expr: self.atom().clone(),
        }
    }

    /// Test whether the expression is exactly zero, including zeros with nonzero tensor rank.
    /// This does not perform algebraic simplification.
    ///
    /// Returns
    /// -------
    /// bool
    ///     True for an exact symbolic zero, irrespective of rank.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import TensorName, Representation
    /// >>> A = TensorName("A")(Representation.euc(3))("i")
    /// >>> assert (0 * A).is_zero()
    /// >>> assert (0 * A).rank == 1
    fn is_zero(&self) -> bool {
        self.atom().is_zero()
    }

    /// Test whether the expression is exactly the scalar one.
    ///
    /// Returns
    /// -------
    /// bool
    ///     True for scalar one; an identity tensor is not the scalar one.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import TensorExpression
    /// >>> assert TensorExpression(1).is_one()
    fn is_one(&self) -> bool {
        self.atom().is_one()
    }

    /// Test whether products and powers are expanded, optionally with respect to var.
    /// The expression and its tensor interface are unchanged.
    ///
    /// Parameters
    /// ----------
    /// var : Expression, optional
    ///     Restrict the expansion check to this indeterminate.
    ///
    /// Returns
    /// -------
    /// bool
    ///     True if no further arithmetic distribution is needed for the requested scope.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica import S
    /// >>> from symbolica.community.tensor import TensorExpression
    /// >>> x = S("x")
    /// >>> value = TensorExpression((x + 1)**2)
    /// >>> assert not value.is_expanded()
    /// >>> assert value.expand().is_expanded(x)
    #[pyo3(signature = (var = None))]
    fn is_expanded(&self, var: Option<ConvertibleToExpression>) -> bool {
        self.atom()
            .is_expanded(var.map(|value| value.to_expression().expr))
    }

    /// Count top-level additive terms without expanding the expression.
    /// Unlike len(tensor), this counts symbolic terms rather than tensor components.
    ///
    /// Returns
    /// -------
    /// int
    ///     Number of terms in the outermost sum, or one for a non-sum.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica import S
    /// >>> from symbolica.community.tensor import TensorExpression
    /// >>> x = S("x")
    /// >>> value = TensorExpression((x + 1)**2)
    /// >>> assert value.nterms() == 1
    /// >>> assert value.expand().nterms() == 3
    fn nterms(&self) -> usize {
        self.atom().nterms()
    }

    /// Replace symbolic patterns or apply whole-tensor replacement rules.
    ///
    /// Parameters
    /// ----------
    /// pattern : Expression, TensorRule, or list of TensorRule
    ///     Symbolica pattern when rhs is supplied. Without rhs, apply prepared
    ///     TensorRules in priority order, leaving scalar key arguments opaque.
    /// rhs : Expression, HeldExpression, or callable, optional
    ///     Symbolica replacement, including callbacks. The result must preserve
    ///     the tensor's ordered interface. Use reindex() to change external labels.
    /// cond : PatternRestriction or Condition, optional
    ///     Restriction on a symbolic match.
    /// non_greedy_wildcards : list of Expression, optional
    ///     Wildcards whose matches should be selected non-greedily.
    /// min_level, max_level : int, optional
    ///     Inclusive traversal bounds, starting at zero.
    /// level_is_tree_depth : bool, default False
    ///     Count expression-tree levels rather than function nesting.
    /// partial : bool, default True
    ///     Allow matches inside larger expressions.
    /// allow_new_wildcards_on_rhs : bool, default False
    ///     Allow wildcard symbols that do not occur in the pattern.
    /// rhs_cache_size : int, optional
    ///     Cache repeated symbolic replacements; zero disables caching.
    /// repeat, once, bottom_up, nested : bool, default False
    ///     Symbolica's replacement traversal controls. Prepared TensorRules own
    ///     their conditions and cache settings and do not accept these options.
    ///
    /// Returns
    /// -------
    /// TensorExpression
    ///     Rewritten value retaining its checked tensor interface and metadata.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica import S
    /// >>> from symbolica.community.tensor import Representation, TensorName, TensorRule
    /// >>> A = TensorName("M")(Representation.euc(2))("i")
    /// >>> x, y = S("x", "y")
    /// >>> assert (x * A).replace(x, y) == y * A
    /// >>> assert A.replace(TensorRule(A, 2 * A)) == 2 * A
    #[allow(clippy::too_many_arguments)]
    #[pyo3(signature = (pattern, rhs=None, cond=None, non_greedy_wildcards=None, min_level=0, max_level=None, level_is_tree_depth=false, partial=true, allow_new_wildcards_on_rhs=false, rhs_cache_size=None, repeat=false, once=false, bottom_up=false, nested=false))]
    fn replace(
        self_: PyRef<'_, Self>,
        py: Python<'_>,
        #[gen_stub(override_type(
            type_repr = "typing.Union[symbolica.core.Expression, int, float, complex, TensorRule, list[TensorRule]]"
        ))]
        pattern: &Bound<'_, PyAny>,
        rhs: Option<ConvertibleToReplaceWith>,
        cond: Option<ConvertibleToPatternRestriction>,
        non_greedy_wildcards: Option<Vec<PythonExpression>>,
        min_level: usize,
        max_level: Option<usize>,
        level_is_tree_depth: bool,
        partial: bool,
        allow_new_wildcards_on_rhs: bool,
        rhs_cache_size: Option<usize>,
        repeat: bool,
        once: bool,
        bottom_up: bool,
        nested: bool,
    ) -> PyResult<Py<Self>> {
        if let Some(rhs) = rhs {
            let result = self_.as_super().replace(
                pattern.extract()?,
                rhs,
                cond,
                non_greedy_wildcards,
                min_level,
                max_level,
                level_is_tree_depth,
                partial,
                allow_new_wildcards_on_rhs,
                rhs_cache_size,
                repeat,
                once,
                bottom_up,
                nested,
            )?;
            return Self::preserving_interface(&self_, py, result.expr);
        }
        if cond.is_some()
            || non_greedy_wildcards.is_some()
            || min_level != 0
            || max_level.is_some()
            || level_is_tree_depth
            || !partial
            || allow_new_wildcards_on_rhs
            || rhs_cache_size.is_some()
            || repeat
            || once
            || bottom_up
            || nested
        {
            return Err(PyTypeError::new_err(
                "TensorRules own their replacement settings; traversal options require a symbolic pattern and rhs",
            ));
        }
        let rules = crate::tensor_rule::PyTensorRule::extract_rules(pattern)?;
        let value = self_
            .value
            .replace_rules(&rules)
            .map_err(Self::inference_error)?;
        let (name, args) = Self::transformed_descriptor(&self_, value.expression());
        Self::from_shared(py, Arc::new(value), name, args)
    }

    /// Apply simultaneous Symbolica replacements and validate the resulting tensor.
    ///
    /// Parameters
    /// ----------
    /// replacements : list of Replacement
    ///     Symbolica replacements, including their conditions and callbacks.
    /// repeat, once, bottom_up, nested : bool, default False
    ///     Traversal controls with the same meaning as Expression.replace_multiple.
    ///
    /// Returns
    /// -------
    /// TensorExpression
    ///     Rewritten value preserving ordered ports, metadata, and typed zeros.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica import S, Replacement
    /// >>> from symbolica.community.tensor import Representation, TensorName
    /// >>> x, y = S("x", "y")
    /// >>> A = TensorName("A")(Representation.euc(3))("i")
    /// >>> assert (x*A).replace_multiple([Replacement(x, y)]) == y*A
    #[pyo3(signature = (replacements, repeat=false, once=false, bottom_up=false, nested=false))]
    fn replace_multiple(
        self_: PyRef<'_, Self>,
        py: Python<'_>,
        replacements: Vec<PythonReplacement>,
        repeat: bool,
        once: bool,
        bottom_up: bool,
        nested: bool,
    ) -> PyResult<Py<Self>> {
        let result =
            self_
                .as_super()
                .replace_multiple(replacements, repeat, once, bottom_up, nested)?;
        Self::preserving_interface(&self_, py, result.expr)
    }

    /// Apply a Symbolica Transformer and retain the validated tensor interface.
    ///
    /// Parameters
    /// ----------
    /// op : Transformer
    ///     Transformations to apply, including user callbacks.
    /// n_cores : int, optional
    ///     Number of cores for Symbolica's term mapping.
    /// stats_to_file : str, optional
    ///     JSON output file for the transformer's statistics.
    ///
    /// Returns
    /// -------
    /// TensorExpression
    ///     Transformed expression with its original ordered tensor interface.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica import S, T
    /// >>> from symbolica.community.tensor import Representation, TensorName
    /// >>> x = S("x")
    /// >>> A = TensorName("A")(Representation.euc(3))("i")
    /// >>> result = ((x + 1)**2 * A).map(T().expand())
    /// >>> assert result == ((x + 1)**2 * A).expand()
    #[pyo3(signature = (op, n_cores=None, stats_to_file=None))]
    fn map(
        self_: PyRef<'_, Self>,
        py: Python<'_>,
        op: PythonTransformer,
        n_cores: Option<usize>,
        stats_to_file: Option<String>,
    ) -> PyResult<Py<Self>> {
        let result = self_.as_super().map(op, py, n_cores, stats_to_file)?;
        Self::preserving_interface(&self_, py, result.expr)
    }

    /// Differentiate scalar coefficients, preserving tensor zeros and ports.
    /// Formal derivatives of unknown tensor functions are not supported by the
    /// network executor. Differentiate their component expressions with
    /// Tensor.map_components() instead, or supply a TensorName derivative callback.
    ///
    /// Parameters
    /// ----------
    /// x : Expression
    ///     Scalar differentiation variable.
    ///
    /// Returns
    /// -------
    /// TensorExpression
    ///     Derivative with the same ordered external ports, even when zero.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica import S
    /// >>> from symbolica.community.tensor import Representation, TensorName
    /// >>> x = S("x")
    /// >>> A = TensorName("A")(Representation.euc(3))("i")
    /// >>> assert (x**2 * A).derivative(x) == 2*x*A
    /// >>> assert A.derivative(x).rank == A.rank
    fn derivative(
        self_: PyRef<'_, Self>,
        py: Python<'_>,
        x: ConvertibleToExpression,
    ) -> PyResult<Py<Self>> {
        let result = self_.as_super().derivative(x)?;
        let mut formal_tensor_derivative = false;
        result.expr.visitor(&mut |value| {
            if let AtomView::Fun(function) = value
                && function.get_symbol() == Symbol::DERIVATIVE
                && function.get_nargs() % 2 == 1
                && matches!(function.iter().nth(function.get_nargs() / 2), Some(AtomView::Var(head)) if head.get_symbol().has_tag(&SPENSO_TAG.tensor))
            {
                formal_tensor_derivative = true;
            }
            !formal_tensor_derivative
        });
        if formal_tensor_derivative {
            return Err(PyValueError::new_err(
                "formal tensor derivatives are not supported; differentiate component expressions with Tensor.map_components()",
            ));
        }
        Self::preserving_interface(&self_, py, result.expr)
    }

    /// Distribute scalar products and powers over sums.
    ///
    /// Parameters
    /// ----------
    /// var : scalar expression, optional
    ///     Restrict expansion to terms involving this expression. None expands
    ///     throughout the expression.
    /// via_poly : bool, optional
    ///     Use Symbolica's polynomial-based expansion when True. Defaults to False.
    ///
    /// Returns
    /// -------
    /// TensorExpression
    ///     Expanded algebra with a checked tensor interface.
    ///
    /// Notes
    /// -----
    /// Label unresolved axes before expansion when their identities cannot be
    /// tracked through a rewritten sum. Expansion can enlarge expressions. Tensor simplifiers
    /// do not generally require expanding the whole expression first.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import Representation, TensorName, TensorExpression
    /// >>> space = Representation.euc(2)
    /// >>> A = TensorName("M")(space, space)
    /// >>> from symbolica import S
    /// >>> x = S("x")
    /// >>> expanded = ((x + 1)**2 * A("i", "j")).expand()
    /// >>> expanded.rank
    /// 2
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
            return Py::new(py, self_.clone().into_initializer());
        }
        let (name, args) = Self::transformed_descriptor(&self_, value.expression());
        Self::from_shared(py, Arc::new(value), name, args)
    }

    /// Factor scalar algebra while retaining ordered tensor ports and data metadata.
    /// Tensor leaves are treated as symbolic indeterminates; no tensor identity is applied.
    ///
    /// Parameters
    /// ----------
    /// complex : bool, default False
    ///     Factor over the complex rationals rather than the rationals.
    /// extension : list of Expression, optional
    ///     Algebraic numbers defining the coefficient field, as in Expression.factor.
    ///
    /// Returns
    /// -------
    /// TensorExpression
    ///     Factored expression with the original tensor interface, including typed zeros.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica import S
    /// >>> from symbolica.community.tensor import TensorName, Representation
    /// >>> x = S("x")
    /// >>> A = TensorName("A")(Representation.euc(3))("i")
    /// >>> factored = ((x**2 + 2*x + 1) * A).factor()
    /// >>> assert factored == (x + 1)**2 * A
    #[pyo3(signature = (complex = false, extension = None))]
    fn factor(
        self_: PyRef<'_, Self>,
        py: Python<'_>,
        complex: bool,
        extension: Option<Vec<ConvertibleToExpression>>,
    ) -> PyResult<Py<Self>> {
        let result = self_.as_super().factor(complex, extension)?;
        Self::from_algebra_atom(&self_, py, result.expr)
    }

    /// Distribute numerical coefficients over sums without expanding symbolic products.
    /// The result remains a TensorExpression with the same ordered ports and metadata.
    ///
    /// Returns
    /// -------
    /// TensorExpression
    ///     Rewritten scalar algebra with the original tensor interface and metadata.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica import S
    /// >>> from symbolica.community.tensor import TensorName, Representation
    /// >>> x, y = S("x", "y")
    /// >>> A = TensorName("A")(Representation.euc(3))("i")
    /// >>> result = (2 * (x + y) * A).expand_num()
    /// >>> assert result == (2*x + 2*y) * A
    fn expand_num(self_: PyRef<'_, Self>, py: Python<'_>) -> PyResult<Py<Self>> {
        let result = self_.as_super().expand_num();
        Self::from_algebra_atom(&self_, py, result.expr)
    }

    /// Collect terms in variables or functions, retaining the tensor interface.
    ///
    /// Parameters
    /// ----------
    /// *x : Expression
    ///     Variables or functions to collect in, with Expression.collect semantics.
    /// key_map, coeff_map : callable, optional
    ///     Map each collected key or coefficient. Callbacks receive and return ordinary
    ///     Expressions. Their results must preserve this tensor's ordered external ports;
    ///     otherwise the operation raises ValueError. Zero retains the original interface.
    ///
    /// Returns
    /// -------
    /// TensorExpression
    ///     Collected expression with the original interface and metadata.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica import S
    /// >>> from symbolica.community.tensor import TensorName, Representation
    /// >>> x, y = S("x", "y")
    /// >>> A = TensorName("A")(Representation.euc(3))("i")
    /// >>> result = (x*A + x*y*A).collect(x)
    /// >>> assert result.expand() == (x*(1+y)*A).expand()
    #[pyo3(signature = (*x, key_map = None, coeff_map = None))]
    fn collect(
        self_: PyRef<'_, Self>,
        py: Python<'_>,
        #[gen_stub(override_type(type_repr = "symbolica.core.Expression"))] x: &Bound<'_, PyTuple>,
        #[gen_stub(override_type(
            type_repr = "typing.Optional[typing.Callable[[symbolica.core.Expression], symbolica.core.Expression]]"
        ))]
        key_map: Option<Py<PyAny>>,
        #[gen_stub(override_type(
            type_repr = "typing.Optional[typing.Callable[[symbolica.core.Expression], symbolica.core.Expression]]"
        ))]
        coeff_map: Option<Py<PyAny>>,
    ) -> PyResult<Py<Self>> {
        let has_callbacks = key_map.is_some() || coeff_map.is_some();
        let result = self_.as_super().collect(x, key_map, coeff_map)?;
        if has_callbacks {
            Self::preserving_interface(&self_, py, result.expr)
        } else {
            Self::from_algebra_atom(&self_, py, result.expr)
        }
    }

    /// Extract common numerical coefficients from sums, retaining ordered ports and metadata.
    /// This reverses numerical distribution performed by expand_num().
    ///
    /// Returns
    /// -------
    /// TensorExpression
    ///     Expression with common numerical coefficients extracted.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica import S
    /// >>> from symbolica.community.tensor import TensorName, Representation
    /// >>> x, y = S("x", "y")
    /// >>> A = TensorName("A")(Representation.euc(3))("i")
    /// >>> assert ((2*x + 4*y)*A).collect_num() == 2*(x + 2*y)*A
    fn collect_num(self_: PyRef<'_, Self>, py: Python<'_>) -> PyResult<Py<Self>> {
        let result = self_.as_super().collect_num();
        Self::from_algebra_atom(&self_, py, result.expr)
    }

    /// Collect common factors from nested sums without applying tensor identities.
    /// Ordered ports, unresolved axes, tensor zeros, and data metadata are retained.
    ///
    /// Returns
    /// -------
    /// TensorExpression
    ///     Rewritten scalar algebra with the original tensor interface and metadata.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica import S
    /// >>> from symbolica.community.tensor import TensorName, Representation
    /// >>> x, y = S("x", "y")
    /// >>> A = TensorName("A")(Representation.euc(3))("i")
    /// >>> result = (x*A + y*A).collect_factors()
    /// >>> assert result == (x+y)*A
    fn collect_factors(self_: PyRef<'_, Self>, py: Python<'_>) -> PyResult<Py<Self>> {
        let result = self_.as_super().collect_factors();
        Self::from_algebra_atom(&self_, py, result.expr)
    }

    /// Group terms with equal numerical coefficients, retaining ordered ports and metadata.
    ///
    /// Returns
    /// -------
    /// TensorExpression
    ///     Expression grouped by numerical coefficient with the original tensor interface.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica import S
    /// >>> from symbolica.community.tensor import TensorExpression
    /// >>> x, y = S("x", "y")
    /// >>> result = TensorExpression(2*x + 2*y).collect_by_coefficient()
    /// >>> assert result.to_expression() == 2*(x + y)
    fn collect_by_coefficient(self_: PyRef<'_, Self>, py: Python<'_>) -> PyResult<Py<Self>> {
        let result = self_.as_super().collect_by_coefficient();
        Self::from_algebra_atom(&self_, py, result.expr)
    }

    /// Collect terms by powers of a variable or calls with the same function head.
    ///
    /// Parameters
    /// ----------
    /// x : Expression
    ///     Variable or function symbol to collect in.
    /// key_map, coeff_map : callable, optional
    ///     Map ordinary Expressions. As with collect(), the result is checked for
    ///     preservation of the tensor's ordered external ports.
    ///
    /// Returns
    /// -------
    /// TensorExpression
    ///     Collected expression with the original interface and metadata.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica import S
    /// >>> from symbolica.community.tensor import TensorName, Representation
    /// >>> x, f = S("x", "f")
    /// >>> A = TensorName("A")(Representation.euc(3))("i")
    /// >>> result = (f(x)*A + 2*f(x)*A).collect_symbol(f)
    /// >>> assert result.expand() == (3*f(x)*A).expand()
    #[pyo3(signature = (x, key_map = None, coeff_map = None))]
    fn collect_symbol(
        self_: PyRef<'_, Self>,
        py: Python<'_>,
        x: PythonExpression,
        #[gen_stub(override_type(
            type_repr = "typing.Optional[typing.Callable[[symbolica.core.Expression], symbolica.core.Expression]]"
        ))]
        key_map: Option<Py<PyAny>>,
        #[gen_stub(override_type(
            type_repr = "typing.Optional[typing.Callable[[symbolica.core.Expression], symbolica.core.Expression]]"
        ))]
        coeff_map: Option<Py<PyAny>>,
    ) -> PyResult<Py<Self>> {
        let has_callbacks = key_map.is_some() || coeff_map.is_some();
        let result = self_.as_super().collect_symbol(x, key_map, coeff_map)?;
        if has_callbacks {
            Self::preserving_interface(&self_, py, result.expr)
        } else {
            Self::from_algebra_atom(&self_, py, result.expr)
        }
    }

    /// Rewrite polynomial algebra in Horner form, retaining ordered ports and metadata.
    ///
    /// Parameters
    /// ----------
    /// vars : list of Expression, optional
    ///     Variable order. Omit it to use Symbolica's heuristic ordering.
    ///
    /// Returns
    /// -------
    /// TensorExpression
    ///     Horner-form expression with the original tensor interface.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica import S
    /// >>> from symbolica.community.tensor import TensorExpression
    /// >>> x, y = S("x", "y")
    /// >>> result = TensorExpression(x**2 + x*y).collect_horner([x, y])
    /// >>> assert result.to_expression() == x*(x + y)
    #[pyo3(signature = (vars = None))]
    fn collect_horner(
        self_: PyRef<'_, Self>,
        py: Python<'_>,
        vars: Option<Vec<PythonExpression>>,
    ) -> PyResult<Py<Self>> {
        let result = self_.as_super().collect_horner(vars)?;
        Self::from_algebra_atom(&self_, py, result.expr)
    }

    /// Combine scalar denominators into one fraction, retaining ordered ports and metadata.
    /// This is ordinary rational algebra and performs no tensor contractions.
    ///
    /// Returns
    /// -------
    /// TensorExpression
    ///     Single rational expression with the original tensor interface.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica import S
    /// >>> from symbolica.community.tensor import TensorExpression
    /// >>> x, y = S("x", "y")
    /// >>> assert TensorExpression(1/x + 1/y).together().to_expression() == (x + y)/(x*y)
    fn together(self_: PyRef<'_, Self>, py: Python<'_>) -> PyResult<Py<Self>> {
        let result = self_.as_super().together()?;
        Self::from_algebra_atom(&self_, py, result.expr)
    }

    /// Cancel common numerator and denominator factors, retaining ordered ports and metadata.
    /// This is ordinary rational algebra and performs no tensor contractions.
    ///
    /// Returns
    /// -------
    /// TensorExpression
    ///     Rational expression with common factors cancelled and the original tensor interface.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica import S
    /// >>> from symbolica.community.tensor import TensorName, Representation
    /// >>> x = S("x")
    /// >>> A = TensorName("A")(Representation.euc(3))("i")
    /// >>> assert ((x**2 - 1)/(x - 1)*A).cancel() == (x + 1)*A
    fn cancel(self_: PyRef<'_, Self>, py: Python<'_>) -> PyResult<Py<Self>> {
        let result = self_.as_super().cancel()?;
        Self::from_algebra_atom(&self_, py, result.expr)
    }

    /// Decompose scalar denominators into partial fractions, retaining the tensor interface.
    ///
    /// Parameters
    /// ----------
    /// *variables : Expression
    ///     Variables for partial fractioning. Omit them to use all indeterminates,
    ///     as in Expression.apart().
    ///
    /// Returns
    /// -------
    /// TensorExpression
    ///     Partial fractions with the original ordered ports and metadata.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica import S
    /// >>> from symbolica.community.tensor import TensorName, Representation
    /// >>> x = S("x")
    /// >>> A = TensorName("A")(Representation.euc(3))("i")
    /// >>> assert (A/(x*(x + 1))).apart(x) == A/x - A/(x + 1)
    #[pyo3(signature = (*variables))]
    fn apart(
        self_: PyRef<'_, Self>,
        py: Python<'_>,
        #[gen_stub(override_type(type_repr = "symbolica.core.Expression"))] variables: &Bound<
            '_,
            PyTuple,
        >,
    ) -> PyResult<Py<Self>> {
        // Zero has no partial fractions; Symbolica's multivariate path cannot
        // lift an empty numerator polynomial into its coefficient ring.
        if self_.atom().is_zero() {
            return Py::new(py, self_.clone().into_initializer());
        }
        let result = self_.as_super().apart(variables)?;
        Self::from_algebra_atom(&self_, py, result.expr)
    }

    /// Copy the expression and its tensor metadata.
    ///
    /// Returns
    /// -------
    /// TensorExpression
    ///     A copy with the same ordered axes and optional data name.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import Representation, TensorName, TensorExpression
    /// >>> space = Representation.euc(2)
    /// >>> A = TensorName("M")(space, space)
    /// >>> import copy
    /// >>> copy.copy(A).shape
    /// (2, 2)
    fn __copy__(self_: PyRef<'_, Self>, py: Python<'_>) -> PyResult<Py<Self>> {
        Py::new(py, self_.clone().into_initializer())
    }

    /// Copy the tensor's immutable expression, ordered axes, and data metadata.
    ///
    /// Parameters
    /// ----------
    /// memo : dict
    ///     Copy bookkeeping maintained by Python's copy module.
    ///
    /// Returns
    /// -------
    /// TensorExpression
    ///     Independent Python wrapper sharing the immutable tensor value.
    ///
    /// Examples
    /// --------
    /// >>> import copy
    /// >>> from symbolica.community.tensor import Representation, TensorName
    /// >>> A = TensorName("A")(Representation.euc(3))
    /// >>> zero = 0 * A
    /// >>> assert copy.deepcopy(zero).structure == zero.structure
    #[pyo3(signature = (memo))]
    fn __deepcopy__(
        self_: PyRef<'_, Self>,
        py: Python<'_>,
        memo: &Bound<'_, PyDict>,
    ) -> PyResult<Py<Self>> {
        let pointer = self_.as_ptr() as usize;
        let copied = Self::__copy__(self_, py)?;
        memo.set_item(pointer, copied.bind(py))?;
        Ok(copied)
    }

    /// Reject pickling that would silently discard the ordered tensor interface.
    /// Use to_expression() when explicitly serializing only the symbolic algebra.
    ///
    /// Raises
    /// ------
    /// TypeError
    ///     Tensor metadata serialization is not supported.
    #[gen_stub(override_return_type(type_repr = "typing.NoReturn"))]
    fn __reduce__(&self) -> PyResult<()> {
        Err(PyTypeError::new_err(
            "TensorExpression pickling is not supported; use to_expression() to serialize only its symbolic algebra",
        ))
    }

    /// Replace four-dimensional Minkowski index spaces by a supplied dimension.
    ///
    /// Parameters
    /// ----------
    /// dimension : int, Expression, or str
    ///     New Lorentz dimension, such as a symbol D. An Expression must be
    ///     a single symbol.
    ///
    /// Returns
    /// -------
    /// TensorExpression
    ///     Updated explicit slots and compact representation arguments.
    ///
    /// Notes
    /// -----
    /// Spinor and color dimensions, scalar coefficients, and Minkowski spaces
    /// with dimensions other than four are unchanged. Apply this before performing
    /// contractions that have already replaced four-dimensional traces by numbers.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica import S
    /// >>> from symbolica.community.tensor import TensorExpression
    /// >>> gamma = TensorExpression.dirac_gamma(4).with_lorentz_dimension(S("D"))
    /// >>> gamma.shape[:2]
    /// (4, 4)
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

    /// Reduce enabled tensor identities and their required index contractions.
    ///
    /// Gamma and color identities are enabled by default; epsilon identities are
    /// opt-in. The default also completes supported structural contractions, so
    /// a second contract() call is normally unnecessary. The input is unchanged.
    ///
    /// Parameters
    /// ----------
    /// gamma, color : bool, default True
    ///     Enable Dirac/Clifford and color identities, respectively. Each family
    ///     authorizes the structural contractions required to apply its identities.
    ///     Set the other families to False when requesting only one family.
    /// epsilon : bool, default False
    ///     Enable supported Levi-Civita pair contractions and their prerequisites.
    ///     These can generate determinant sums. Epsilon tensors emitted by gamma
    ///     reduction do not implicitly enable this option. Intrinsic antisymmetry
    ///     normalization during construction is independent of this setting.
    /// contract : {"fully", "dots", "minimal", "selected", "none"}, default "fully"
    ///     Structural work in addition to enabled identities and their prerequisites:
    ///
    ///     - "fully": contract compatible metrics/identities and vectors, collect
    ///       matrix chains, and close compatible chains into trace notation.
    ///     - "dots": the same work, plus canonical dot notation for scalar products.
    ///     - "minimal": perform structural contractions that do not multiply out
    ///       independent sum alternatives. Boundary substitutions through existing
    ///       sums and contractions within their branches are allowed. Products of
    ///       sums stay factored when contracting them would require distribution.
    ///     - "selected": the same structural operations, restricted to connections
    ///       in representations. This does not restrict identity prerequisites.
    ///     - "none": no additional structural pass. Enabled identities still get
    ///       their prerequisite contractions; this does not disable all contraction.
    ///
    ///     These modes do not enable identity families. In particular, collecting
    ///     a trace does not evaluate it, and contracting a color delta does not
    ///     apply a generator or Fierz identity. "minimal" restricts additional
    ///     structural work; enabled algebra identities retain their prerequisites
    ///     and can still introduce sums. It is not a global ban on expansion.
    /// representations : sequence of Representation or RepresentationName, optional
    ///     Required with contract="selected" and rejected with the other modes.
    ///     Select connecting representation families, including their duals, rather
    ///     than whole tensor factors. For example, Representation.mink(4) selects
    ///     the Minkowski family, including symbolic-D instances; it does not change
    ///     dimensions. Each connection must still have compatible dimensions and
    ///     orientation. [] permits no additional index contractions or trace closure.
    ///     Identity prerequisites and existing constructor normalization are unaffected.
    /// collect_coefficients : bool, default True
    ///     Collect generated exact coefficients in contracted Dirac traces and
    ///     their bound vector factors so equal terms combine and cancel. Includes
    ///     Symbolica's expand_num() distribution of numerical factors over sums
    ///     inside those generated coefficients. False retains their nested,
    ///     factored output. Other family kernels retain their
    ///     existing coefficient normalization. This is not global expansion:
    ///     unrelated input sums and scalar prefactors remain factorized. It neither
    ///     expands free Dirac traces nor reverses simplification already performed.
    /// gamma_output : {"reduced", "chains"}, optional
    ///     Defaults to "reduced": apply supported Clifford and trace identities.
    ///     "chains" collects Dirac chains/traces without Clifford reduction or
    ///     trace evaluation, even if gamma_evaluate_traces=True. Explicit gamma0
    ///     and gamma_conjugate preparation still applies when requested.
    /// gamma_ordering : {"repeated_pairs", "canonical"}, optional
    ///     Ordering for reduced open Dirac chains. The default "repeated_pairs"
    ///     brings repeated indices together without fully sorting unrelated
    ///     matrices. "canonical" orders the chain and can introduce metric-term
    ///     sums through anticommutation. It is inactive for gamma_output="chains".
    /// gamma_evaluate_traces : bool, optional
    ///     Evaluate supported closed Dirac traces with reduced output; defaults to
    ///     True. False retains trace notation while other enabled identities remain
    ///     available. Trace normalization follows the spinor representation dimension.
    /// gamma0 : bool, optional
    ///     Apply gamma-zero factoring before chain collection; defaults to False.
    ///     This controls the preparation pass, not every gamma-zero identity already
    ///     included in reduced Dirac algebra.
    /// gamma_conjugate : bool, optional
    ///     Rewrite explicitly complex-conjugated Dirac matrices before collection;
    ///     defaults to False. This does not conjugate the input tensor.
    /// gamma_expand_three_gamma_epsilon : bool, optional
    ///     With reduced output, expand three compatible four-dimensional gammas
    ///     in a gamma5/epsilon basis; defaults to False. This can introduce sums
    ///     and epsilon tensors, but does not enable epsilon identities.
    /// color_evaluate_traces : bool, optional
    ///     Evaluate supported closed color traces; defaults to True. False retains
    ///     unevaluated traces, but enabled Fierz identities may still join or change
    ///     them. It does not disable the rest of color algebra.
    /// color_expand_fierz : bool, optional
    ///     Apply fundamental Fierz identities between different open chains or
    ///     traces; defaults to True. This may introduce sums. False disables these
    ///     cross-chain expansions, not every color identity or source of sums.
    /// color_substitute_cof_dimension_invariants : bool, optional
    ///     Write supported color invariants in terms of the fundamental SU(N)
    ///     dimension, with T_F=1/2; defaults to False, retaining symbolic invariants.
    ///     In the supported fundamental/adjoint spaces this gives
    ///     C_F=(N**2-1)/(2*N) and C_A=N. The dimension comes from the representation:
    ///     this does not set every color space to SU(3) or reduce arbitrary groups.
    /// max_steps_per_domain : int or None, default None
    ///     Optional nonnegative limit on successful planner transformations in each
    ///     active expression domain. None imposes no automatic work limit. This is
    ///     not a bound on time, memory, terms, or work inside one kernel. Zero leaves
    ///     eligible work untouched and reports Capped; an already completed request
    ///     can still report Complete. Unfinished results remain exact.
    ///
    /// Returns
    /// -------
    /// TensorExpression
    ///     A typed result carrying its checked interface and reduction_status.
    ///     No alias wrapper or separate materialization step is required.
    ///
    /// Raises
    /// ------
    /// ValueError
    ///     An option value is invalid, representations is missing or used with the
    ///     wrong mode, a family-specific option is supplied for a disabled family,
    ///     or a rewrite cannot produce a valid tensor interface.
    ///
    /// Notes
    /// -----
    /// Select the families together; the shared planner schedules identities and
    /// contractions until no further permitted work is found. Under "fully" and
    /// "dots", a closed color delta produces its dimension even with color=False.
    /// Under "none", an unrelated Lorentz numerator can remain indexed and factored
    /// after its color sector is reduced. Use contract() for structural work alone.
    ///
    /// This is not a global arithmetic expansion. Planning starts with a shallow
    /// view of factor interfaces and opens selected regions when needed. Necessary
    /// local distribution and identity-generated sums are allowed; unrelated sums
    /// and scalar coefficients stay factored. Request expand() explicitly when the
    /// final output must be a fully expanded polynomial.
    ///
    /// Dimensions belong to the tensor representations. Set them before reduction,
    /// for example with with_lorentz_dimension(D). Ordinary Clifford contractions
    /// support compatible symbolic Lorentz dimensions. Four-dimensional gamma5,
    /// axial-trace and epsilon-basis rules are not a general D-dimensional gamma5
    /// prescription; changing D after evaluation cannot recover lost dimension factors.
    ///
    /// Family-specific None values use the defaults documented above. Even an
    /// explicitly supplied False option requires its family to be enabled. Reusable
    /// configurations are ordinary dictionaries passed with **options.
    ///
    /// Inspect reduction_status on the returned reduction result. Complete means
    /// no eligible work remains under these settings, not under disabled families
    /// or other output requirements. Deferred retains exact work the current
    /// implementation could not finish. Capped means an explicit limit stopped it.
    /// Rerunning a Complete result with the same settings is stable. Capped or
    /// Deferred results can be passed back for further work, without a guarantee
    /// that an unsupported Deferred case will progress. contraction_complete is a
    /// separate structural certificate and can be False for a Complete "none" result.
    ///
    /// See Also
    /// --------
    /// contract : Structural contractions without enabling identities.
    /// to_dots : Change surviving scalar-product notation without contracting indices.
    /// expand : Explicit arithmetic distribution.
    /// undo_trace : Unfold an existing trace without evaluating it.
    ///
    /// Examples
    /// --------
    /// A closed Dirac word evaluates to 4*D even with no additional structural
    /// pass. The spinor dimension remains four while the Lorentz dimension is D.
    ///
    /// >>> from symbolica import E, S
    /// >>> from symbolica.community import tensor as sp
    /// >>> D = S("algebra_docs::D")
    /// >>> gamma = sp.TensorExpression.dirac_gamma(4).with_lorentz_dimension(D)
    /// >>> word = gamma("a", "b", "mu") * gamma("b", "a", "mu")
    /// >>> result = word.simplify_algebra(color=False, contract="none")
    /// >>> assert result.to_expression() == 4 * D
    /// >>> assert result.reduction_status == sp.ReductionStatus.Complete
    ///
    /// Retain nested identity output when that form is preferable. Both requests
    /// represent the same tensor and leave unrelated scalar prefactors factored.
    ///
    /// >>> nested = word.simplify_algebra(color=False, collect_coefficients=False)
    /// >>> assert (nested - result).expand().to_expression() == 0
    ///
    /// Color-only reduction with explicit SU(3) invariants gives T^a T^a = 4/3 I.
    /// Omit the invariant-substitution option to retain symbolic Casimirs.
    ///
    /// >>> fundamental = sp.Representation.cof(3)
    /// >>> identity = sp.TensorExpression.g(fundamental, fundamental.dual())
    /// >>> T = sp.TensorExpression.color_t(8, 3)
    /// >>> pair = T("a", "i", "j") * T("a", "j", "k")
    /// >>> color_only = dict(gamma=False, color=True, contract="none",
    /// ...                   color_substitute_cof_dimension_invariants=True)
    /// >>> assert pair.simplify_algebra(**color_only) == E("4/3") * identity("i", "k")
    ///
    /// Structural color-delta contraction is independent of color identities.
    /// With both identity families disabled, "none" preserves this product;
    /// "fully" and a matching "selected" filter produce its dimension.
    ///
    /// >>> cycle = identity("i", "j") * identity("j", "i")
    /// >>> structural = dict(gamma=False, color=False)
    /// >>> assert cycle.simplify_algebra(**structural, contract="none") == cycle
    /// >>> assert cycle.simplify_algebra(**structural).to_expression() == 3
    /// >>> selected = cycle.simplify_algebra(**structural, contract="selected",
    /// ...                                    representations=[fundamental.name])
    /// >>> assert selected.to_expression() == 3
    /// >>> assert cycle.simplify_algebra(**structural, contract="selected",
    /// ...                               representations=[]) == cycle
    ///
    /// "dots" includes scalar-product notation conversion. Epsilon reduction is
    /// separately enabled; here the Euclidean identity epsilon_ijk epsilon_ijk = 6
    /// is evaluated without gamma or color identities.
    ///
    /// >>> lorentz = sp.Representation.mink(4)
    /// >>> p = sp.TensorName.vector("algebra_docs::p")(lorentz)
    /// >>> q = sp.TensorName.vector("algebra_docs::q")(lorentz)
    /// >>> dotted = (p("mu") * q("mu")).simplify_algebra(**structural, contract="dots")
    /// >>> assert dotted == sp.dot(p, q)
    /// >>> epsilon = sp.TensorExpression.levi_civita(sp.Representation.euc(3), rank=3)
    /// >>> product = epsilon("i", "j", "k") * epsilon("i", "j", "k")
    /// >>> assert product.simplify_algebra(**structural, epsilon=True).to_expression() == 6
    ///
    /// Minimal contraction leaves a product of tensor sums exact and factorized.
    /// It can report Complete even though full contraction would do more work.
    /// The equivalent structural-only request is contract(expand=False).
    ///
    /// >>> metric = sp.TensorExpression.g(lorentz)
    /// >>> A = sp.TensorName("algebra_docs::A")(lorentz, lorentz)
    /// >>> factored = (metric("mu", "nu") + A("mu", "nu")) * (p("mu") + q("mu"))
    /// >>> minimal = factored.simplify_algebra(**structural, contract="minimal")
    /// >>> assert minimal == factored
    /// >>> assert minimal.reduction_status == sp.ReductionStatus.Complete
    /// >>> assert minimal == factored.contract(expand=False)
    /// >>> assert minimal != factored.contract()
    ///
    /// An explicit limit preserves the exact remaining expression. Continue by
    /// calling the same operation without a limit; no manual expansion is needed.
    ///
    /// >>> pending = cycle.simplify_algebra(**structural, max_steps_per_domain=0)
    /// >>> assert pending.reduction_status == sp.ReductionStatus.Capped
    /// >>> assert pending == cycle
    /// >>> assert pending.simplify_algebra(**structural).to_expression() == 3
    #[pyo3(signature=(*, gamma=true, color=true, epsilon=false, contract="fully", representations=None, collect_coefficients=true, gamma_output=None, gamma_ordering=None, gamma_evaluate_traces=None, gamma0=None, gamma_conjugate=None, gamma_expand_three_gamma_epsilon=None, color_evaluate_traces=None, color_expand_fierz=None, color_substitute_cof_dimension_invariants=None, max_steps_per_domain=None))]
    #[allow(clippy::too_many_arguments)]
    fn simplify_algebra(
        &self,
        py: Python<'_>,
        gamma: bool,
        color: bool,
        epsilon: bool,
        #[gen_stub(override_type(type_repr="typing.Literal['fully', 'minimal', 'dots', 'selected', 'none']", imports=("typing")))]
        contract: &str,
        #[gen_stub(override_type(type_repr="typing.Optional[typing.Sequence[Representation | RepresentationName]]", imports=("typing")))]
        representations: Option<Vec<Bound<'_, PyAny>>>,
        collect_coefficients: bool,
        #[gen_stub(override_type(type_repr="typing.Optional[typing.Literal['reduced', 'chains']]", imports=("typing")))]
        gamma_output: Option<&str>,
        #[gen_stub(override_type(type_repr="typing.Optional[typing.Literal['repeated_pairs', 'canonical']]", imports=("typing")))]
        gamma_ordering: Option<&str>,
        gamma_evaluate_traces: Option<bool>,
        gamma0: Option<bool>,
        gamma_conjugate: Option<bool>,
        gamma_expand_three_gamma_epsilon: Option<bool>,
        color_evaluate_traces: Option<bool>,
        color_expand_fierz: Option<bool>,
        color_substitute_cof_dimension_invariants: Option<bool>,
        max_steps_per_domain: Option<usize>,
    ) -> PyResult<Py<Self>> {
        let settings = crate::simplification::algebra_settings(
            gamma,
            color,
            epsilon,
            contract,
            representations,
            collect_coefficients,
            gamma_output,
            gamma_ordering,
            gamma_evaluate_traces,
            gamma0,
            gamma_conjugate,
            gamma_expand_three_gamma_epsilon,
            color_evaluate_traces,
            color_expand_fierz,
            color_substitute_cof_dimension_invariants,
            max_steps_per_domain,
        )?;
        let value = self
            .value
            .simplify_algebra(&settings)
            .map_err(Self::inference_error)?;
        Self::from_shared(py, Arc::new(value), self.name, self.name_args.clone())
    }

    /// Give explicit indices a scope that separates them from other copies.
    ///
    /// Parameters
    /// ----------
    /// header : Expression
    ///     A single Symbolica symbol naming the scope.
    /// dummies_only : bool, default False
    ///     Scope only contracted indices, leaving external labels unchanged.
    ///
    /// Returns
    /// -------
    /// TensorExpression
    ///     Scoped indices with the same internal contraction relationships.
    ///
    /// Notes
    /// -----
    /// Applying the same outer scope twice is idempotent; different scopes nest.
    /// Alphabet display uses primed labels for scoped copies. This is useful when
    /// multiplying independently summed expressions or an amplitude and its adjoint.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import Representation, TensorName, TensorExpression
    /// >>> space = Representation.euc(2)
    /// >>> A = TensorName("M")(space, space)
    /// >>> from symbolica import S
    /// >>> scoped = A("i", "j").wrap_indices(S("bra"))
    /// >>> scoped.rank
    /// 2
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

    /// List the external indices that are not summed over.
    ///
    /// Returns
    /// -------
    /// list of Expression
    ///     Representation-aware index expressions, in logical order when all
    ///     axes are labeled. Dual indices retain their dual wrapper.
    ///
    /// Notes
    /// -----
    /// For unresolved axes, prefer ``structure.axes`` to inspect the interface.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import Representation, TensorName, TensorExpression
    /// >>> space = Representation.euc(2)
    /// >>> A = TensorName("M")(space, space)
    /// >>> len(A("i", "j").list_dangling())
    /// 2
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

    /// Canonicalize tensor factors and dummy-index labels.
    ///
    /// Returns
    /// -------
    /// TensorExpression
    ///     An equivalent expression with deterministically named contracted indices.
    ///
    /// Raises
    /// ------
    /// CanonicalizationError
    ///     The tensor expression cannot be canonicalized.
    ///
    /// Notes
    /// -----
    /// This identifies different spellings of the same dummy-index contractions;
    /// it is not a general polynomial simplifier.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import Representation, TensorName, TensorExpression
    /// >>> space = Representation.euc(2)
    /// >>> A = TensorName("M")(space, space)
    /// >>> canonical = A("i", "j").canonize()
    /// >>> canonical.rank
    /// 2
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

    /// Construct the Dirac adjoint of a spinor tensor expression.
    ///
    /// Parameters
    /// ----------
    /// preserve_indices : bool, default False
    ///     Keep external labels on their original physical legs. False uses
    ///     matrix-adjoint ordering and exchanges open-chain endpoints.
    ///
    /// Returns
    /// -------
    /// TensorExpression
    ///     The complex conjugate with reversed spinor chains and the gamma^0
    ///     boundary factors required by a Dirac adjoint.
    ///
    /// Raises
    /// ------
    /// DiracAdjointError
    ///     The index structure does not admit a consistent Dirac adjoint.
    ///
    /// Notes
    /// -----
    /// The expression must use the registered Dirac representation and tensor
    /// forms. This operation is specialized to Dirac spinors; it is not an ordinary
    /// transpose of arbitrary tensor components.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import TensorName, Representation
    /// >>> spinor = TensorName("psi")(Representation.bis(4)("a"))
    /// >>> spinor.dirac_adjoint().rank
    /// 1
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

    /// Write compact metric products using dot notation.
    ///
    /// Returns
    /// -------
    /// TensorExpression
    ///     An equivalent value in the requested notation.
    ///
    /// Notes
    /// -----
    /// This changes notation; it does not itself contract explicit vector indices.
    /// Use contract() first when those contractions are needed.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import Representation, TensorName, TensorExpression
    /// >>> space = Representation.euc(2)
    /// >>> A = TensorName("M")(space, space)
    /// >>> converted = A.to_dots()
    fn to_dots(self_: PyRef<'_, Self>, py: Python<'_>) -> PyResult<Py<Self>> {
        let value = self_.value.to_dots().map_err(Self::inference_error)?;
        let (name, args) = Self::transformed_descriptor(&self_, value.expression());
        Self::from_shared(py, Arc::new(value), name, args)
    }

    /// Write dot products as explicit indexed contractions.
    ///
    /// Returns
    /// -------
    /// TensorExpression
    ///     An equivalent value in the requested notation.
    ///
    /// Notes
    /// -----
    /// Other compact forms and factored sums are retained. This does not
    /// evaluate component data or expand a polynomial.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import Representation, TensorName, TensorExpression
    /// >>> space = Representation.euc(2)
    /// >>> A = TensorName("M")(space, space)
    /// >>> converted = A.undo_dots()
    fn undo_dots(self_: PyRef<'_, Self>, py: Python<'_>) -> PyResult<Py<Self>> {
        let value = self_.value.undo_dots().map_err(Self::inference_error)?;
        let (name, args) = Self::transformed_descriptor(&self_, value.expression());
        Self::from_shared(py, Arc::new(value), name, args)
    }

    /// Expose collected chain factors with fresh compatible dummy indices.
    ///
    /// Returns
    /// -------
    /// TensorExpression
    ///     Equivalent algebra in the requested notation.
    ///
    /// Notes
    /// -----
    /// This changes notation without distributing polynomial sums or requesting
    /// Dirac/color identities. Use undo_trace().undo_chain() to expose indexed trace factors.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import Representation, TensorName, chain, trace
    /// >>> space = Representation.euc(2)
    /// >>> A = TensorName("M")(space, space)
    /// >>> expression = chain(space("i"), space("j"), A, A)
    /// >>> opened = expression.undo_chain()
    fn undo_chain(self_: PyRef<'_, Self>, py: Python<'_>) -> PyResult<Py<Self>> {
        let value = self_.value.undo_chain().map_err(Self::inference_error)?;
        let (name, args) = Self::transformed_descriptor(&self_, value.expression());
        Self::from_shared(py, Arc::new(value), name, args)
    }

    /// Open compact traces as closed chains without evaluating them.
    ///
    /// Returns
    /// -------
    /// TensorExpression
    ///     Equivalent algebra in the requested notation.
    ///
    /// Notes
    /// -----
    /// This changes notation without distributing polynomial sums or requesting
    /// Dirac/color identities. Use undo_trace().undo_chain() to expose indexed trace factors.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import Representation, TensorName, chain, trace
    /// >>> space = Representation.euc(2)
    /// >>> A = TensorName("M")(space, space)
    /// >>> expression = trace(space, A, A)
    /// >>> opened = expression.undo_trace()
    fn undo_trace(self_: PyRef<'_, Self>, py: Python<'_>) -> PyResult<Py<Self>> {
        let value = self_.value.undo_trace().map_err(Self::inference_error)?;
        let (name, args) = Self::transformed_descriptor(&self_, value.expression());
        Self::from_shared(py, Arc::new(value), name, args)
    }

    /// Build an executable network for this symbolic tensor.
    ///
    /// Parameters
    /// ----------
    /// library : TensorLibrary, optional
    ///     Component definitions. Defaults to the built-in four-dimensional Dirac
    ///     and SU(3) library. Unregistered tensors receive symbolic components.
    ///
    /// Returns
    /// -------
    /// TensorNetwork
    ///     A new network with the same external axes. Construction does not
    ///     execute its contractions.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import Representation, TensorName, TensorExpression
    /// >>> space = Representation.euc(2)
    /// >>> A = TensorName("M")(space, space)
    /// >>> network = A.to_network()
    /// >>> network.rank
    /// 2
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

    /// Evaluate this expression and return its components.
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
    /// The original expression and the supplied libraries are unchanged.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import Representation, TensorName, TensorExpression
    /// >>> space = Representation.euc(2)
    /// >>> A = TensorName("M")(space, space)
    /// >>> A.to_tensor().shape
    /// (2, 2)
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

    /// Evaluate this tensor and return its flat component list.
    ///
    /// Parameters
    /// ----------
    /// library : TensorLibrary, optional
    ///     Component definitions. Defaults to the built-in four-dimensional Dirac
    ///     and SU(3) library. Unregistered tensors receive symbolic components.
    ///
    /// Returns
    /// -------
    /// list of Expression, float, or complex
    ///     Components in logical row-major order; one element for a scalar.
    ///     Dimensions must be concrete.
    ///
    /// Notes
    /// -----
    /// Unregistered tensors receive symbolic components. Relabeling abstract
    /// indices does not change those component identities.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import Representation, TensorName, TensorExpression
    /// >>> space = Representation.euc(2)
    /// >>> A = TensorName("M")(space, space)
    /// >>> len(A.components())
    /// 4
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
    /// TensorExpression
    ///     An indexed expression.
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
    /// >>> indexed = A.index("i", "j")
    /// >>> indexed.rank
    /// 2
    #[pyo3(signature = (*indices, intern = None))]
    fn index(
        self_: PyRef<'_, Self>,
        py: Python<'_>,
        indices: &Bound<'_, PyTuple>,
        intern: Option<Intern>,
    ) -> PyResult<Py<Self>> {
        let replacements = Self::index_replacements(self_.interface(), indices, intern)?;
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
    /// TensorExpression
    ///     An indexed expression.
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
    /// >>> indexed = A.reindex("i", "j")
    /// >>> indexed.rank
    /// 2
    #[pyo3(signature = (*indices, intern=None))]
    fn reindex(
        self_: PyRef<'_, Self>,
        py: Python<'_>,
        indices: &Bound<'_, PyTuple>,
        intern: Option<Intern>,
    ) -> PyResult<Py<Self>> {
        let positions = (0..self_.interface().canonical().order()).collect::<Vec<_>>();
        let replacements = Self::port_replacements(self_.interface(), &positions, indices, intern)?;
        let value = Self::structured(&self_)
            .reindex_interface_ports(&replacements)
            .map_err(|e| PyValueError::new_err(e.to_string()))?;
        let (atom, interface) = value.into_parts();
        Self::from_known_parts(py, atom, interface, self_.name, self_.name_args.clone())
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
    /// TensorExpression
    ///     A new value with renamed axes; the original is unchanged.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import Representation, TensorName, TensorExpression
    /// >>> space = Representation.euc(2)
    /// >>> A = TensorName("M")(space, space)
    /// >>> indexed = A.index("i", "j")
    /// >>> renamed = indexed.rename_indices({"i": "k"})
    /// >>> renamed.rank
    /// 2
    #[pyo3(signature = (mapping, *, intern=None))]
    fn rename_indices(
        self_: PyRef<'_, Self>,
        py: Python<'_>,
        mapping: &Bound<'_, PyDict>,
        intern: Option<Intern>,
    ) -> PyResult<Py<Self>> {
        let replacements = Self::named_replacements(self_.interface(), mapping, intern)?;
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
    /// TensorExpression
    ///     A new tensor view with the reordered interface.
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
    /// >>> transposed = A.permute_axes([1, 0])
    /// >>> transposed.shape
    /// (2, 2)
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
    /// TensorExpression
    ///     An indexed expression.
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
    /// >>> indexed = A("i", "j")
    /// >>> indexed.rank
    /// 2
    #[pyo3(signature = (*indices, intern = None))]
    fn __call__(
        self_: PyRef<'_, Self>,
        py: Python<'_>,
        indices: &Bound<'_, PyTuple>,
        intern: Option<Intern>,
    ) -> PyResult<Py<Self>> {
        Self::index(self_, py, indices, intern)
    }

    /// Negate every tensor component.
    ///
    /// Returns
    /// -------
    /// TensorExpression
    ///     Negated symbolic algebra with the same external axes.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import Representation, TensorName, TensorExpression
    /// >>> space = Representation.euc(2)
    /// >>> A = TensorName("M")(space, space)
    /// >>> negated = -A
    fn __neg__(self_: PyRef<'_, Self>, py: Python<'_>) -> PyResult<Py<Self>> {
        Self::from_known_parts(
            py,
            Atom::num(-1) * self_.atom(),
            self_.interface().clone(),
            None,
            Vec::new(),
        )
    }

    #[doc = python_doc!("TensorExpression.__add__")]
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

    #[doc = python_doc!("TensorExpression.__radd__")]
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

    #[doc = python_doc!("TensorExpression.__sub__")]
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

    #[doc = python_doc!("TensorExpression.__rsub__")]
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

    #[doc = python_doc!("TensorExpression.__mul__")]
    #[gen_stub(skip)]
    fn __mul__(
        self_: PyRef<'_, Self>,
        py: Python<'_>,
        rhs: &Bound<'_, PyAny>,
    ) -> PyResult<Py<PyAny>> {
        match rhs.getattr("__symbolica_rmul__") {
            Ok(product) => return product.call1((&self_,)).map(Bound::unbind),
            Err(error) if error.is_instance_of::<pyo3::exceptions::PyAttributeError>(py) => {}
            Err(error) => return Err(error),
        }
        if let Some(right) = concrete_network(rhs)? {
            return Self::promoted_network(&self_, py)?
                .multiply_network(right)
                .and_then(|network| Py::new(py, network))
                .map(TensorDispatch::Network)
                .and_then(|result| result.into_py_any(py));
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
        .and_then(|result| result.into_py_any(py))
    }

    /// Preserve the tensor type when a Symbolica expression multiplies this tensor.
    ///
    /// Parameters
    /// ----------
    /// lhs : Expression
    ///     Left operand; a tensor expression retains its tensor structure.
    ///
    /// Returns
    /// -------
    /// TensorExpression
    ///     The same ordered product as ``lhs * self``.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica import S
    /// >>> from symbolica.community.tensor import Representation, TensorName
    /// >>> vector = TensorName.vector("v")(Representation.euc(2))
    /// >>> (S("x") * vector).rank
    /// 1
    #[gen_stub(override_return_type(type_repr = "TensorExpression"))]
    fn __symbolica_rmul__(
        self_: PyRef<'_, Self>,
        py: Python<'_>,
        lhs: &Bound<'_, PythonExpression>,
    ) -> PyResult<TensorDispatch> {
        Self::__rmul__(self_, py, lhs.as_any())
    }

    #[doc = python_doc!("TensorExpression.__rmul__")]
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

    #[doc = python_doc!("TensorExpression.__truediv__")]
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

    #[doc = python_doc!("TensorExpression.__rtruediv__")]
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

    /// Raise tensor algebra to a scalar power.
    ///
    /// Parameters
    /// ----------
    /// exponent : TensorExpression or scalar expression
    ///     Scalar exponent.
    /// modulo : object, optional
    ///     Must be None; modular tensor powers are not supported.
    ///
    /// Returns
    /// -------
    /// TensorExpression
    ///     Powered tensor expression. A non-scalar base requires a nonnegative
    ///     integer exponent and uses repeated tensor multiplication, including
    ///     its contraction and ambiguity rules. This is not an elementwise power.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import TensorExpression
    /// >>> squared = TensorExpression(3)**2
    /// >>> squared.to_expression() == 9
    /// True
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

    /// Use a scalar tensor expression as an exponent.
    ///
    /// Parameters
    /// ----------
    /// base : TensorExpression or scalar expression
    ///     Scalar base.
    /// modulo : None, optional
    ///     Must be None; modular powers are not supported.
    ///
    /// Returns
    /// -------
    /// TensorExpression
    ///     Scalar power with a checked tensor interface.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import TensorExpression
    /// >>> result = 2**TensorExpression(3)
    /// >>> result.to_expression() == 8
    /// True
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
    ///     TensorExpression to compare. Values without an ordered tensor interface are unequal.
    fn __eq__(&self, other: &Bound<'_, PyAny>) -> bool {
        // Returning NotImplemented lets Expression's reflected comparison ignore
        // ordered tensor ports, violating equality transitivity and our hash.
        other
            .extract::<PyRef<'_, Self>>()
            .is_ok_and(|other| self.value.same_value(&other.value))
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
    ///     TensorExpression to compare. Values without an ordered tensor interface are unequal.
    fn __ne__(&self, other: &Bound<'_, PyAny>) -> bool {
        !self.__eq__(other)
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

    #[doc = python_doc!("TensorExpression.outer")]
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

    /// Report completion under the last reduction settings and budget.
    ///
    /// Returns
    /// -------
    /// ReductionStatus
    ///     Complete, Deferred, or Capped. Unfinished results remain exact and can be
    ///     passed to contract() or simplify_algebra() again.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import TensorExpression, ReductionStatus
    /// >>> TensorExpression(3).contract().reduction_status == ReductionStatus.Complete
    /// True
    #[getter]
    fn reduction_status(&self) -> crate::simplification::PyReductionStatus {
        self.value.reduction_status().into()
    }

    /// Whether symbolic metric/vector contraction is certified complete.
    ///
    /// Returns
    /// -------
    /// bool
    ///     False also covers a value not yet contracted or stopped before completion.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import TensorExpression
    /// >>> TensorExpression(3).contract().contraction_complete
    /// True
    #[getter]
    fn contraction_complete(&self) -> bool {
        self.value.contraction_complete()
    }

    /// Contract compatible indices without applying tensor-algebra identities.
    ///
    /// Eliminate metrics and identity tensors, including closed dimension factors;
    /// substitute vectors into compatible ports; form scalar products; and collect
    /// ordered matrix chains and unevaluated traces. Matrix order and orientation
    /// are preserved. The input is unchanged.
    ///
    /// Parameters
    /// ----------
    /// representations : sequence of Representation or RepresentationName, optional
    ///     Filter connections by representation family: None permits all, [] none.
    ///     Families include dual representations. For example, Representation.mink(4)
    ///     permits the Minkowski family, including symbolic-D slots; the 4 is not
    ///     a dimension filter. Each connection must still satisfy dimension and
    ///     orientation compatibility. Matching representation names alone is insufficient.
    ///     The filter acts on connecting ports, not whole tensor factors: a Lorentz
    ///     connection can enter a gamma matrix while its spinor ports survive.
    ///     Chain assembly and trace closure also require a permitted connection.
    /// metrics : bool, default True
    ///     Eliminate compatible metric/identity tensors by index substitution and
    ///     evaluate closed metric loops to their representation dimension.
    /// rank_one : bool, default True
    ///     Substitute vectors into compatible tensor ports (Schoonschip notation)
    ///     and contract vector pairs into scalar products. This can enter a gamma
    ///     matrix without applying a gamma identity. False retains those connections.
    /// collect_chains : bool, default True
    ///     Represent connected matrix products as ordered chains without evaluating
    ///     their matrices, including words with supplied rank-one endpoints.
    ///     Supplied endpoints do not add external axes. False retains the indexed
    ///     factors for open products.
    /// collect_traces : bool, default True
    ///     Represent compatible closed matrix products as unevaluated traces.
    ///     False retains indexed closure, including when collect_chains=True.
    ///     It does not unfold trace notation already present in the input.
    /// expand : bool, default True
    ///     Permit local arithmetic distribution required by a selected contraction.
    ///     False performs the minimal structural work: substitute at boundaries
    ///     and contract within existing sum branches, but do not multiply out
    ///     independent sum alternatives. A product of tensor sums can therefore
    ///     remain indexed and still be Complete under expand=False. Neither value
    ///     requests global polynomial expansion; expand() is a separate operation.
    /// order : sequence of int, optional
    ///     Optional zero-based permutation of every normalized top-level factor,
    ///     used for the initial contraction round. Factors follow Symbolica's
    ///     normalized order, not necessarily the order in the Python source. None
    ///     lets the planner choose. This affects scheduling, never matrix order;
    ///     omit it unless you need to prescribe a factor order explicitly.
    /// max_steps_per_domain : int or None, default None
    ///     Optional nonnegative limit on successful planner transformations in each
    ///     active expression domain. None imposes no automatic work limit. This is
    ///     not a limit on time, memory, generated terms, or work inside one kernel.
    ///     Zero preserves eligible work and reports Capped; a request with no eligible
    ///     work can still report Complete. All unfinished results remain exact.
    ///
    /// Returns
    /// -------
    /// TensorExpression
    ///     A typed result with its checked external interface and reduction_status.
    ///     Unrelated scalar coefficients and sums remain factored. No alias wrapper
    ///     or separate materialization step is required.
    ///
    /// Raises
    /// ------
    /// ValueError
    ///     order is not a permutation of all normalized top-level factors, or a
    ///     substitution cannot produce a valid tensor interface.
    ///
    /// Notes
    /// -----
    /// Contracting a color delta does not enable generator or Fierz identities.
    /// Collecting a Dirac or color trace does not evaluate it. Use simplify_algebra()
    /// for identities, contract_ports(rhs, left=..., right=...) for an explicit
    /// binary port contraction, and to_tensor() for component-library evaluation.
    /// With rank_one=False, collect_chains=False and collect_traces=False, only
    /// metric/identity contraction is requested.
    ///
    /// Planning starts with a depth-one view of factor interfaces. It opens a
    /// selected region only when the contraction needs its interior, reusing
    /// unchanged leaves and known interfaces. With expand=True, necessary local
    /// distribution is allowed; an unrelated scalar sum is not routinely expanded.
    /// Internal dummy contractions can be processed even inside a scalar-interface
    /// leaf. This is distinct from expanding the whole expression with expand().
    ///
    /// Settings do not undo intrinsic normalization performed during construction.
    /// A metric with identical compatible ports can already be its dimension,
    /// typed vector multiplication can already form a scalar product, and symmetric
    /// dot arguments are already ordered. representations=[] and disabled flags
    /// preserve these existing normalizations. Assign explicit compatible indices
    /// before requesting a particular Einstein contraction; unresolved axes are
    /// distinct port identities, not dummy indices to pair by guesswork.
    ///
    /// Scalar products can remain in compact metric notation g(p(R), q(R)).
    /// to_dots() converts that surviving notation to dot(p(R), q(R)) without
    /// searching for repeated indices. undo_dots(), undo_chain() and undo_trace()
    /// unfold existing notation with fresh compatible dummy indices, preserving
    /// matrix order. They do not recover expressions already changed by identities.
    ///
    /// Complete means no eligible work remains under the requested settings;
    /// disabled contractions and unevaluated traces may remain. Deferred retains
    /// exact work the current implementation could not finish. Capped means an
    /// explicit work limit stopped the operation. Completed reruns with the same
    /// settings are stable. Unfinished results can be passed back to contract(),
    /// possibly with different settings; unsupported cases need not make progress.
    /// Inspect this operation's reduction_status before other transformations.
    /// contraction_complete is a separate structural certificate, not a substitute
    /// for completion under a restricted representation filter or expand=False.
    ///
    /// See Also
    /// --------
    /// simplify_algebra : Gamma, color and epsilon identities plus their contractions.
    /// contract_ports : Explicit binary contraction of selected logical ports.
    /// to_dots : Scalar-product notation conversion without index contraction.
    /// undo_trace : Turn trace notation into a closed indexed chain.
    /// expand : Explicit arithmetic distribution.
    ///
    /// Examples
    /// --------
    /// A metric relabels a vector. A Lorentz-only filter also permits a metric to
    /// enter a gamma matrix without touching its two spinor ports.
    ///
    /// >>> from symbolica import S
    /// >>> from symbolica.community import tensor as sp
    /// >>> lorentz = sp.Representation.mink(4)
    /// >>> metric = sp.TensorExpression.g(lorentz)
    /// >>> p = sp.TensorName.vector("contract_docs::p")(lorentz)
    /// >>> source = metric("mu", "nu") * p("nu")
    /// >>> assert source.contract() == p("mu")
    /// >>> gamma = sp.TensorExpression.dirac_gamma(4)
    /// >>> mixed = gamma("a", "b", "nu") * metric("mu", "nu")
    /// >>> assert mixed.contract(representations=[lorentz]) == gamma("a", "b", "mu")
    ///
    /// Trace collection and trace evaluation are separate operations. The collected
    /// trace represents the same tensor as 4*g(mu,nu), but contract() leaves it
    /// unevaluated. Unfolding and contracting it recovers the collected notation.
    ///
    /// >>> word = gamma("a", "b", "mu") * gamma("b", "a", "nu")
    /// >>> collected = word.contract()
    /// >>> assert collected.simplify_algebra(color=False) == 4 * metric("mu", "nu")
    /// >>> assert collected.undo_trace().undo_chain().contract() == collected
    /// >>> uncollected = word.contract(collect_traces=False)
    /// >>> assert uncollected.contract() == collected
    ///
    /// A supplied spinor remains attached when collecting an open gamma word.
    /// Only the unresolved spinor slot remains an external axis.
    ///
    /// >>> spinor = sp.Representation.bis(4)
    /// >>> psi = sp.TensorName.vector("contract_docs::psi")(spinor)
    /// >>> word = psi("a") * gamma("a", "b", "mu") * p("mu") * gamma("b", sp.AUTO, "nu") * p("nu")
    /// >>> collected = word.contract()
    /// >>> assert collected.rank == 1
    /// >>> assert collected.contract() == collected
    ///
    /// A closed color delta gives N_c without any color-algebra identity. An empty
    /// filter retains this product, and an unrelated scalar power remains factored.
    ///
    /// >>> Nc = S("contract_docs::Nc")
    /// >>> fundamental = sp.Representation.cof(Nc)
    /// >>> delta = sp.TensorExpression.g(fundamental, fundamental.dual())
    /// >>> cycle = delta("i", "j") * delta("j", "i")
    /// >>> assert cycle.contract().to_expression() == Nc
    /// >>> assert cycle.contract(representations=[]) == cycle
    /// >>> x, y = S("contract_docs::x", "contract_docs::y")
    /// >>> assert ((x + y)**8 * cycle).contract().to_expression() == (x + y)**8 * Nc
    ///
    /// Minimal contraction still relabels ports through a sum. It leaves products
    /// of independent tensor sums factorized instead of enumerating their pairings.
    ///
    /// >>> q = sp.TensorName.vector("contract_docs::q")(lorentz)
    /// >>> A = sp.TensorName("contract_docs::A")(lorentz, lorentz)
    /// >>> relabel = metric("mu", "nu") * (p("nu") + q("nu"))
    /// >>> assert relabel != p("mu") + q("mu")
    /// >>> assert relabel.contract(expand=False) == p("mu") + q("mu")
    /// >>> product = (metric("mu", "nu") + A("mu", "nu")) * (p("mu") + q("mu"))
    /// >>> minimal = product.contract(expand=False)
    /// >>> assert minimal == product
    /// >>> assert minimal.reduction_status == sp.ReductionStatus.Complete
    /// >>> full = product.contract()
    /// >>> assert minimal != full
    /// >>> assert minimal.contract() == full
    ///
    /// An explicit zero limit retains eligible work exactly. Continuing without
    /// that limit completes it, while the original input remains available.
    ///
    /// >>> pending = cycle.contract(max_steps_per_domain=0)
    /// >>> assert pending.reduction_status == sp.ReductionStatus.Capped
    /// >>> assert pending == cycle
    /// >>> assert pending.contract().to_expression() == Nc
    #[pyo3(signature=(*, representations=None, metrics=true, rank_one=true, collect_chains=true, collect_traces=true, expand=true, order=None, max_steps_per_domain=None))]
    #[allow(clippy::too_many_arguments)]
    fn contract(
        &self,
        py: Python<'_>,
        #[gen_stub(override_type(type_repr="typing.Optional[typing.Sequence[Representation | RepresentationName]]", imports=("typing")))]
        representations: Option<Vec<Bound<'_, PyAny>>>,
        metrics: bool,
        rank_one: bool,
        collect_chains: bool,
        collect_traces: bool,
        expand: bool,
        order: Option<Vec<usize>>,
        max_steps_per_domain: Option<usize>,
    ) -> PyResult<Py<Self>> {
        let representations = crate::simplification::contraction_representations(representations)?;
        let settings = idenso::tensor::ContractSettings {
            representations: representations.as_deref(),
            metrics,
            rank_one,
            collect_chains,
            collect_traces,
            expand,
            order: order.as_deref(),
            max_passes: max_steps_per_domain,
        };
        let value = self
            .value
            .contract(settings)
            .map_err(Self::inference_error)?;
        Self::from_shared(py, Arc::new(value), self.name, self.name_args.clone())
    }

    #[doc = python_doc!("TensorExpression.contract_ports")]
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

    #[doc = python_doc!("TensorExpression.compose")]
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
    /// TensorExpression
    ///     The traced tensor; spectator axes remain external.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import Representation, TensorName, TensorExpression
    /// >>> space = Representation.euc(2)
    /// >>> A = TensorName("M")(space, space)
    /// >>> A.trace(channel=(0, 1)).is_scalar
    /// True
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
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import Representation, TensorName, TensorExpression
    /// >>> space = Representation.euc(2)
    /// >>> A = TensorName("M")(space, space)
    /// >>> output = A.format_tensor()
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

    /// Produce LaTeX for the symbolic tensor expression.
    ///
    /// Parameters
    /// ----------
    /// max_line_length : int, optional
    ///     Requested line-length bound for the printer. None leaves it unbounded.
    /// settings : DisplaySettings, optional
    ///     Index labels, dimensions and invariant notation. Defaults to
    ///     DisplaySettings(). Only ports layout and default spacing are supported.
    ///
    /// Returns
    /// -------
    /// str
    ///     LaTeX tensor notation, including the enclosing math delimiters.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import Representation, TensorName, TensorExpression
    /// >>> space = Representation.euc(2)
    /// >>> A = TensorName("M")(space, space)
    /// >>> source = A.to_latex()
    #[pyo3(signature = (max_line_length = None, *, settings = None))]
    fn to_latex(
        self_: PyRef<'_, Self>,
        max_line_length: Option<usize>,
        settings: Option<PyRef<'_, display::DisplaySettings>>,
    ) -> PyResult<String> {
        let settings = display::resolved_settings(None, settings.as_deref());
        display::reject_renderer_only_settings(&settings, "to_latex")?;
        Ok(display::structured_to_latex(
            Self::structured(&self_),
            &settings,
            max_line_length,
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
    /// Static Typst output remains available independently of the HTML explorer.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import Representation, TensorName, TensorExpression
    /// >>> space = Representation.euc(2)
    /// >>> A = TensorName("M")(space, space)
    /// >>> output = A.to_typst()
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
    ///     Small expressions render once per backend; large notebook displays use bounded pages.
    ///
    /// Notes
    /// -----
    /// Automatic representations of large expressions are bounded previews.
    /// Paging retains factorization and ordered ports; explicit exports remain complete.
    /// Settings are captured here; custom printers run on the first request
    /// for each backend. Create a new formatted value after changing printers.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import Representation, TensorName, TensorExpression
    /// >>> space = Representation.euc(2)
    /// >>> A = TensorName("M")(space, space)
    /// >>> output = A.formatted()
    #[pyo3(signature = (show_dimensions = None, *, settings = None, notation_source = None))]
    fn formatted(
        self_: PyRef<'_, Self>,
        show_dimensions: Option<bool>,
        settings: Option<PyRef<'_, display::DisplaySettings>>,
        notation_source: Option<String>,
    ) -> PythonFormattedOutput {
        let settings = display::resolved_settings(show_dimensions, settings.as_deref());
        display::format_structured_output_rich(Arc::clone(&self_.value), settings, notation_source)
    }

    /// Create a bounded, interactive MathML viewer in marimo or Jupyter.
    ///
    /// Only the current page is rendered. Nested expressions retain their
    /// factorization, and omitted subexpressions can be opened independently.
    /// Requires the optional notebook-display dependencies for live navigation.
    ///
    /// Parameters
    /// ----------
    /// page_size : int, default 25
    ///     Maximum entries shown on a page: 25, 100, 250, or 500. Rendering may
    ///     show fewer entries when an individual expression is large.
    /// settings : DisplaySettings, optional
    ///     Tensor notation, index labels, and other presentation preferences.
    /// notation_source : str, optional
    ///     Custom Typst notation source, with the same meaning as in to_html().
    ///
    /// Returns
    /// -------
    /// object
    ///     Notebook viewer retaining the expression and rendering pages on demand.
    ///     Display it in marimo or Jupyter to navigate its widget.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import TensorExpression
    /// >>> viewer = TensorExpression(1).paged(page_size=25)
    #[pyo3(signature = (page_size=25, *, settings=None, notation_source=None))]
    fn paged(
        self_: PyRef<'_, Self>,
        py: Python<'_>,
        page_size: usize,
        settings: Option<PyRef<'_, display::DisplaySettings>>,
        notation_source: Option<String>,
    ) -> PyResult<Py<PyAny>> {
        display::paging::viewer(
            py,
            Arc::clone(&self_.value),
            settings.as_deref().cloned().unwrap_or_default(),
            notation_source,
            page_size,
        )
    }

    /// Return marimo's notebook presentation with bounded pages for large expressions.
    ///
    /// Returns
    /// -------
    /// object
    ///     Rich notebook output or an interactive viewer, using default settings.
    fn _display_(self_: PyRef<'_, Self>, py: Python<'_>) -> PyResult<Py<PyAny>> {
        display::format_structured_output_rich(
            Arc::clone(&self_.value),
            display::DisplaySettings::default(),
            None,
        )
        ._display_(py)
    }

    /// Return the rich representations requested by a Jupyter frontend.
    ///
    /// Parameters
    /// ----------
    /// include, exclude : list of str, optional
    ///     MIME types to include or exclude from the returned representations.
    ///
    /// Returns
    /// -------
    /// object
    ///     MIME bundle for the complete small expression or paged large expression.
    #[pyo3(signature = (include=None, exclude=None))]
    fn _repr_mimebundle_(
        self_: PyRef<'_, Self>,
        py: Python<'_>,
        include: Option<Vec<String>>,
        exclude: Option<Vec<String>>,
    ) -> PyResult<Py<PyAny>> {
        display::format_structured_output_rich(
            Arc::clone(&self_.value),
            display::DisplaySettings::default(),
            None,
        )
        ._repr_mimebundle_(py, include, exclude)
    }

    /// Render the tensor as an HTML fragment.
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
    ///     Self-contained output fragment for an HTML-capable notebook or page.
    ///
    /// Notes
    /// -----
    /// Mathematical rendering uses the embedded Typst compiler. The returned
    /// string is not automatically displayed; pass it to the notebook's HTML
    /// or SVG display facility.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import Representation, TensorName, TensorExpression
    /// >>> space = Representation.euc(2)
    /// >>> A = TensorName("M")(space, space)
    /// >>> output = A.to_html()
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
    /// Mathematical rendering uses the embedded Typst compiler. The returned
    /// string is not automatically displayed; pass it to the notebook's HTML
    /// or SVG display facility.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import Representation, TensorName, TensorExpression
    /// >>> space = Representation.euc(2)
    /// >>> A = TensorName("M")(space, space)
    /// >>> output = A.to_svg()
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
        if display::paging::is_large(&self_.value) {
            "Tensor expression (large; use the paged notebook viewer or explicit export)".into()
        } else {
            display::format_structured(Self::structured(&self_), false)
        }
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
        if display::paging::is_large(&self_.value) {
            "Tensor expression (large; use the paged notebook viewer or explicit export)".into()
        } else {
            display::format_structured(Self::structured(&self_), false)
        }
    }

    fn _repr_pretty_(
        self_: PyRef<'_, Self>,
        pretty: &Bound<'_, PyAny>,
        cycle: bool,
    ) -> PyResult<()> {
        let text = if cycle {
            "...".to_string()
        } else if display::paging::is_large(&self_.value) {
            "Tensor expression (large; use the paged notebook viewer or explicit export)".into()
        } else {
            display::format_structured_backend(
                pretty.py(),
                Self::structured(&self_),
                &display::DisplaySettings::default(),
                None,
                FormattedOutputBackend::Text,
            )?
        };
        pretty.call_method1("text", (text,))?;
        Ok(())
    }

    fn _repr_html_(self_: PyRef<'_, Self>, py: Python<'_>) -> Option<String> {
        if display::paging::is_large(&self_.value) {
            return display::paging::static_html(
                py,
                Arc::clone(&self_.value),
                display::DisplaySettings::default(),
                None,
            )
            .ok();
        }
        let settings = display::DisplaySettings::default();
        display::format_structured_backend(
            py,
            Self::structured(&self_),
            &settings,
            None,
            FormattedOutputBackend::Html,
        )
        .ok()
    }

    fn _repr_latex_(self_: PyRef<'_, Self>) -> PyResult<String> {
        if display::paging::is_large(&self_.value) {
            Ok(r"$$\text{Large tensor expression: use the paged viewer}$$".into())
        } else {
            Self::to_latex(self_, None, None)
        }
    }
}

/// Convert a symbolic expression to a TensorExpression.
///
/// Parameters
/// ----------
/// expression : TensorExpression or scalar expression
///     Tensor syntax to infer, or an ordinary scalar. An existing
///     TensorExpression retains its tensor metadata.
///
/// Returns
/// -------
/// TensorExpression
///     Tensor-aware algebra, including rank-zero scalars.
///
/// Notes
/// -----
/// Equivalent to ``TensorExpression(expression)`` without its optional
/// structure and index-cooking settings.
///
/// Examples
/// --------
/// >>> from symbolica.community.tensor import as_tensor
/// >>> as_tensor(2).is_scalar
/// True
#[cfg_attr(
    feature = "python_stubgen",
    pyo3_stub_gen::derive::gen_stub_pyfunction(module = "symbolica.community.tensor")
)]
#[spenso_macros::track_usage(crate::record_usage)]
#[pyfunction]
pub fn as_tensor(py: Python<'_>, expression: &Bound<'_, PyAny>) -> PyResult<Py<TensorExpression>> {
    if let Ok(expression) = expression.extract::<PyRef<'_, TensorExpression>>() {
        return Py::new(py, expression.clone().into_initializer());
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
    ($name:literal, $doc:expr, [$($argument:literal: $input:ty),*], [$($value:literal: $kind:ident),+]) => {
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
symbolic_function_overload!("dot", python_doc!("dot"), [], ["left": PositionalOrKeyword, "right": PositionalOrKeyword]);

#[doc = python_doc!("dot")]
#[cfg_attr(
    feature = "python_stubgen",
    pyo3_stub_gen::derive::gen_stub_pyfunction(
        module = "symbolica.community.tensor",
        no_default_overload = true,
        python_overload = r#"
        import typing

        @overload
        def dot(left: typing.Union[Tensor, TensorNetwork], right: pyo3_stub_gen.RustType["ConvertibleToSpensoNet"]) -> pyo3_stub_gen.RustType["SpensoNet"]:
            ...
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
symbolic_function_overload!("chain", python_doc!("chain"), ["start_slot": SpensoSlot, "end_slot": SpensoSlot], ["factors": VarPositional]);

#[doc = python_doc!("chain")]
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
symbolic_function_overload!("trace", python_doc!("trace"), ["representation": SpensoRepresentation], ["factors": VarPositional]);

#[doc = python_doc!("trace")]
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

/// Return the canonical real color-count symbol, whose default numerical value is 3.
///
/// Construct it on demand so importing Symbolica leaves time to set a license key.
/// Use a separate dimension symbol for formal SU(N) calculations when the default
/// numerical value is not appropriate.
///
/// Examples
/// --------
/// >>> from symbolica.community.tensor import Nc
/// >>> Nc().evaluate({})
/// 3
#[cfg_attr(
    feature = "python_stubgen",
    pyo3_stub_gen::derive::gen_stub_pyfunction(module = "symbolica.community.tensor")
)]
#[spenso_macros::track_usage(crate::record_usage)]
#[pyfunction(name = "Nc")]
fn nc() -> PythonExpression {
    PythonExpression::from(Atom::var(CS.nc))
}

pub(crate) fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    TensorExpression::init(m)?;
    m.add_function(wrap_pyfunction!(nc, m)?)?;
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
                doc: python_doc!("TensorExpression.__getitem__"),
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
                doc: python_doc!("TensorExpression.__getitem__"),
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
                doc: python_doc!("TensorExpression.__getitem__"),
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
                .simplify_algebra(&idenso::tensor::simplification::AlgebraSettings {
                    color: Some(idenso::color::ColorSimplifySettings::default()),
                    ..Default::default()
                })
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
                    .simplify_algebra(&idenso::tensor::AlgebraSettings {
                        epsilon: true,
                        ..Default::default()
                    })
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
                                .expanded(None, false)
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
                    let compact = TensorExpression::from_structured(py, tensor.borrow(py).value.contract(Default::default()).map_err(TensorExpression::inference_error)?)?;
                    let kwargs = pyo3::types::PyDict::new(py);
                    kwargs.set_item("gamma", true)?;
                    kwargs.set_item("gamma_output", "chains")?;
                    let collected = compact.bind(py).call_method("simplify_algebra", (), Some(&kwargs))?;
                    let collected = collected.extract::<PyRef<'_, TensorExpression>>()?;
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
                "dirac_gamma",
                "gamma5",
                "projm",
                "projp",
                "sigma",
                "color_f",
                "color_t",
                "rank",
                "is_scalar",
                "structure",
                "name",
                "with_name",
                "to_expression",
                "expand",
                "simplify_algebra",
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
                    expression_type.call_method1("dirac_gamma", (4,))?,
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
                    expression_type.call_method1("color_f", (8,))?,
                    vec![adjoint, adjoint, adjoint],
                ),
                (
                    expression_type.call_method1("color_t", (8, 3))?,
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
                let canonical_representations = function
                    .iter()
                    .map(|argument| {
                        if function.get_symbol() == ETS.metric {
                            // The slot-aware metric factory keeps unresolved
                            // arguments compact; its interface owns their
                            // distinct occurrence-local identities.
                            assert!(Slot::<LibraryRep, AbstractIndex>::try_from(argument).is_err());
                            Representation::<LibraryRep>::try_from(argument)
                                .expect("an unresolved metric argument must be a representation")
                        } else {
                            Slot::<LibraryRep, AbstractIndex>::try_from(argument)
                                .expect("factory ports must carry distinct open indices")
                                .rep()
                        }
                    })
                    .collect::<Vec<_>>();
                assert_eq!(
                    canonical_representations,
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
                    TensorExpression::g(py, SpensoSlotOrArgOrRep::Rep(euc_object.clone()), None)?,
                    vec![dimension; 2],
                ),
                (TensorExpression::flat(py, &euc_object)?, vec![dimension; 2]),
                (
                    TensorExpression::dirac_gamma(py, ConvertibleToDimension(dimension))?,
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
                    TensorExpression::color_f(py, ConvertibleToDimension(adjoint_dimension))?,
                    vec![adjoint_dimension; 3],
                ),
                (
                    TensorExpression::color_t(
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
                        .all(|argument| if function.get_symbol() == ETS.metric {
                            Slot::<LibraryRep, AbstractIndex>::try_from(argument).is_err()
                                && Representation::<LibraryRep>::try_from(argument).is_ok()
                        } else {
                            Slot::<LibraryRep, AbstractIndex>::try_from(argument).is_ok()
                        })
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
        assert!(error.to_string().contains("TensorExpression.color_t"));

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
        assert!(error.to_string().contains("TensorExpression.dirac_gamma"));
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
            let metric =
                TensorExpression::g(py, SpensoSlotOrArgOrRep::Rep(representation.clone()), None)?;
            let tensor_type = py.get_type::<TensorExpression>();
            let kwargs = pyo3::types::PyDict::new(py);
            kwargs.set_item("intern", "indices")?;
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
                    let simplified = cooked.bind(py).call_method0("contract")?;
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
    fn tensor_constructor_interns_only_indices_and_preserves_tags() {
        use symbolica::atom::UserData;

        idenso::representations::initialize();
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let payload = FunctionBuilder::new(symbol!(
                "constructor_index_source",
                tags = ["idenso::constructor_input_index"]
            ))
            .add_arg(Atom::num(7))
            .finish();
            let argument = FunctionBuilder::new(symbol!("constructor_scalar_argument"))
                .add_arg(Atom::num(9))
                .finish();
            let representation = Minkowski {}.new_rep(4).cast::<LibraryRep>();
            let raw = FunctionBuilder::new(SPENSO_TAG.tensor_symbol("constructor_interned_tensor"))
                .add_arg(argument.clone())
                .add_arg(representation.to_symbolic([payload.clone()]))
                .finish();
            let expression = PythonExpression { expr: raw.clone() }.into_pyobject(py)?;
            let tensor_type = py.get_type::<TensorExpression>();
            for (mode, policy) in [
                ("indices", Intern::Indices),
                ("flattened", Intern::Flattened),
            ] {
                let kwargs = PyDict::new(py);
                kwargs.set_item("intern", mode)?;
                let result = tensor_type.call((&expression,), Some(&kwargs))?;
                let result = result.extract::<PyRef<'_, TensorExpression>>()?;
                let slots = result.interface().logical_slots();
                assert_eq!(slots.len(), 1);
                assert_eq!(slots[0].rep(), representation);
                assert_eq!(result.name_args, vec![argument.clone()]);
                let PartialIndex::Explicit(AbstractIndex::Symbol(symbol)) = slots[0].aind else {
                    panic!("interning must retain an explicit external index");
                };
                let symbol: Symbol = symbol.into();
                assert!(symbol.has_tag("idenso::constructor_input_index"));
                if mode == "indices" {
                    assert_eq!(symbol.get_data(), &UserData::Atom(payload.clone()));
                    assert_eq!(policy.rust().uncook(result.atom().as_view()), raw);
                } else {
                    assert_ne!(symbol.get_data(), &UserData::Atom(payload.clone()));
                    assert_eq!(
                        *result.atom(),
                        policy.rust().try_cook_indices(raw.as_view()).unwrap()
                    );
                }
            }
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
            let simplified = wrapped.bind(py).call_method0("contract")?;
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
            let open = TensorExpression::dirac_gamma(py, ConvertibleToDimension(4.into()))?;
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
                    let result = expression
                        .bind(py)
                        .call_method1(method, args)?
                        .extract::<Py<TensorExpression>>()?;
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
            let tensor = TensorExpression::dirac_gamma(py, ConvertibleToDimension(4.into()))?;
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
            let error = expression
                .bind(py)
                .call_method("collect_symbol", (x,), Some(&kwargs))
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
            let tensor = TensorExpression::dirac_gamma(py, ConvertibleToDimension(4.into()))?;
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
            assert_eq!(&constructed.as_super().expr, constructed.atom());
            assert!(tensor.bind(py).is_instance_of::<PythonExpression>());
            let inherited = tensor.bind(py).call_method0("get_all_symbols")?;
            let inherited = inherited.extract::<Vec<PythonExpression>>()?;
            assert!(inherited == tensor.borrow(py).to_expression().get_all_symbols(true));
            let ordinary = tensor.borrow(py).to_expression().into_pyobject(py)?;
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
            let factored = expression.bind(py).call_method0("factor")?;
            let collected = factored
                .call_method1("collect", (&x,))?
                .extract::<Py<TensorExpression>>()?;
            let collected = collected.borrow(py);
            assert_eq!(&collected.as_super().expr, collected.atom());
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
            let error = factored
                .call_method("collect", (&x,), Some(&kwargs))
                .unwrap_err();
            assert!(error.is_instance_of::<PyValueError>(py));
            let zero = TensorExpression::from_atom_interface(
                py,
                Atom::Zero,
                tensor.borrow(py).interface().clone(),
            )?;
            let zero = zero
                .bind(py)
                .call_method0("factor")?
                .extract::<Py<TensorExpression>>()?;
            assert_eq!(zero.borrow(py).interface(), tensor.borrow(py).interface());

            let d = Dimension::from(symbol!("tensor_convenience_D"));
            let promoted = TensorExpression::with_lorentz_dimension(
                tensor.borrow(py),
                py,
                ConvertibleToDimension(d),
            )?;
            let expected = TensorExpression::dirac_gamma(py, ConvertibleToDimension(d))?;
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

            let ordinary = expression.bind(py).borrow().to_expression();
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
                let result = expression.bind(py).call_method0("contract")?;
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
                let kwargs = pyo3::types::PyDict::new(py);
                kwargs.set_item("gamma", true)?;
                kwargs.set_item(
                    "gamma_output",
                    if output == idenso::dirac::GammaOutput::Reduced {
                        "reduced"
                    } else {
                        "chains"
                    },
                )?;
                kwargs.set_item("epsilon", output == idenso::dirac::GammaOutput::Reduced)?;
                let result =
                    gamma_expression
                        .bind(py)
                        .call_method("simplify_algebra", (), Some(&kwargs))?;
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
            let metric_only = tensor.bind(py).call_method("contract", (), Some(&kwargs))?;
            let value = metric_only.extract::<PyRef<'_, TensorExpression>>()?;
            assert_eq!(
                value.value.expression(),
                &(spenso::p!(&nu) * spenso::q!(&nu))
            );
            assert!(!value.value.contraction_complete());
            let again = metric_only.call_method("contract", (), Some(&kwargs))?;
            let again = again.extract::<PyRef<'_, TensorExpression>>()?;
            assert_eq!(again.value, value.value);
            let full = metric_only.call_method0("contract")?;
            let full = full.extract::<PyRef<'_, TensorExpression>>()?;
            assert_eq!(
                full.value.expression(),
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
            let settings = Intern::Flattened;
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
            kwargs.set_item("intern", "flattened")?;
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
                TensorExpression::new(py, zero.bind(py).as_any(), None, Some(settings))?,
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
                TensorExpression::new(py, metric.bind(py).as_any(), None, Some(settings))?,
            )?;
            assert!(cooked.borrow(py).interface().canonical().is_scalar());
            let contraction = PyDict::new(py);
            contraction.set_item("collect_chains", false)?;
            contraction.set_item("collect_traces", false)?;
            let simplified = cooked
                .bind(py)
                .call_method("contract", (), Some(&contraction))?;
            let simplified = simplified.extract::<PyRef<'_, TensorExpression>>()?;
            assert_eq!(simplified.value.expression(), &Atom::num(3));
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
            let simplified = original.bind(py).call_method0("simplify_algebra")?;
            let simplified = simplified.extract::<PyRef<'_, TensorExpression>>()?;
            let simplified = simplified
                .value
                .expanded(None, false)
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
                ])
                .into_canonical(),
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
                        Canonicalized::identity(structure.interface.clone())
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
                    let result = original
                        .bind(py)
                        .call_method1(method, args)?
                        .extract::<Py<TensorExpression>>()?;
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
                let simplified = original.bind(py).call_method0("simplify_algebra")?;
                let simplified = simplified.extract::<PyRef<'_, TensorExpression>>()?;
                let simplified = simplified
                    .value
                    .expanded(None, false)
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
            let kwargs = pyo3::types::PyDict::new(py);
            kwargs.set_item("color", true)?;
            let simplified = tensor
                .bind(py)
                .call_method("simplify_algebra", (), Some(&kwargs))?;
            let simplified = simplified.extract::<PyRef<'_, TensorExpression>>()?;
            assert!(simplified.interface().canonical().is_scalar());
            assert_eq!(
                *simplified.atom(),
                Atom::num(-8) * CS.cas(Atom::num(2), adjoint.to_symbolic([]))
            );

            let open = TensorExpression::color_f(py, ConvertibleToDimension(8.into()))?;
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
            let metric = TensorExpression::g(
                py,
                SpensoSlotOrArgOrRep::Rep(SpensoRepresentation { representation }),
                None,
            )?;
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
