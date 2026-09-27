use std::collections::HashMap;

use idenso::{
    Cookable, IndexTooling,
    color::{CS, ColorCasimirSettings, ColorSimplifier},
    dirac::{AGS, GammaSimplifier},
    epsilon::EpsilonSimplifier,
    selective_expand::SelectiveExpand,
    shorthands::{
        UndoShorthands,
        chain::Chain,
        metric::MetricSimplifier,
        schoonschip::{Schoonschip, SchoonschipSettings},
    },
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
        ConvertibleToExpression, ConvertibleToPatternRestriction, ConvertibleToReplaceWith,
        PythonExpression, PythonFormattedOutput, PythonReplacement, PythonTransformer,
    },
    atom::{Atom, AtomCore, AtomView, Symbol},
};

use crate::{
    ModuleInit, SliceOrIntOrExpanded, Spensor, display,
    library::SpensorLibrary,
    network::{ConvertibleToSpensoNet, ExecutionMode, ReplacementCondition, SpensoNet},
    simplification::{
        CanonicalizationError, CookingError, DiracAdjointError, DotExpansionError,
        GammaConjugationError, NetworkToolingError, PyColorCasimirSettings,
        PyColorSimplifySettings, PyCookSettings, PyGammaSimplifySettings, PySchoonschipSettings,
        expansion::{PythonTerm, python_terms},
    },
    structure::{
        ArithmeticStructure, ConvertibleToAbstractIndex, ConvertibleToDimension,
        ConvertibleToSpensoName, SpensoIndices, SpensoName, SpensoRepresentation, SpensoSlot,
    },
};

/// The local placeholder used to leave a tensor port unresolved.
#[cfg_attr(feature = "python_stubgen", pyo3_stub_gen::derive::gen_stub_pyclass)]
#[pyclass(frozen, name = "_AutoIndex", module = "symbolica.community.spenso")]
pub struct AutoIndex;

#[cfg(feature = "python_stubgen")]
pyo3_stub_gen::module_variable!("symbolica.community.spenso", "AUTO", AutoIndex);
#[cfg(feature = "python_stubgen")]
pyo3_stub_gen::module_variable!("symbolica.community.spenso", "_", AutoIndex);

#[pymethods]
impl AutoIndex {
    fn __repr__(&self) -> &'static str {
        "_"
    }

    fn __str__(&self) -> &'static str {
        "_"
    }
}

/// A Symbolica expression with an ordered external tensor interface.
///
/// Predefined tensors are constructed by typed factories and indexed afterward
/// in logical interface order. Dimensions are representation metadata, not
/// scalar tensor arguments.
///
/// Examples
/// --------
/// >>> from symbolica.community.spenso import Representation, TensorExpression
/// >>> mink = Representation.mink(4)
/// >>> metric = TensorExpression.g(mink)("mu", "nu")
/// >>> flat = TensorExpression.flat(mink)("mu", "nu")
/// >>> gamma = TensorExpression.gamma(4)("i", "j", "mu")
/// >>> gamma5 = TensorExpression.gamma5(4)("i", "j")
/// >>> projm = TensorExpression.projm(4)("i", "j")
/// >>> projp = TensorExpression.projp(4)("i", "j")
/// >>> sigma = TensorExpression.sigma(4)("mu", "nu", "i", "j")
/// >>> structure_constant = TensorExpression.f(8)("a", "b", "c")
/// >>> generator = TensorExpression.t(8, 3)("a", "i", "j")
#[cfg_attr(feature = "python_stubgen", pyo3_stub_gen::derive::gen_stub_pyclass)]
#[pyclass(
    frozen,
    extends = PythonExpression,
    from_py_object,
    name = "TensorExpression",
    module = "symbolica.community.spenso"
)]
#[derive(Clone)]
pub struct TensorExpression {
    pub(crate) interface: PartialStructure,
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
            let interface = InterfaceInference::merge_explicit_interface_sequence(&[interface])
                .map_err(Self::inference_error)?
                .canonicalize_open_ports();
            SymbolicTensor::validate_interface(&atom, &interface).map_err(Self::inference_error)?;
            return Self::from_known_parts(py, atom, interface, None, Vec::new());
        }
        let (name, name_args) = inferred_descriptor(atom.as_view());
        let value = if InterfaceInference::has_structured_syntax(atom.as_view()) {
            SymbolicTensor::infer(atom).map_err(Self::inference_error)?
        } else {
            SymbolicTensor::new(atom, PartialStructure::from_logical_slots([]))
        };
        Self::from_known_parts(py, value.expression, value.structure, name, name_args)
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
        let atom = value.expression;
        let interface = value.structure;
        // A contraction is a derived expression, not a declaration of the original stored data.
        let (name, name_args) = if interface.canonical().order() == original_rank {
            (name, name_args)
        } else {
            (None, Vec::new())
        };
        Self::from_parts_unchecked(
            py,
            atom,
            interface.canonicalize_open_ports(),
            name,
            name_args,
        )
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
        Py::new(
            py,
            (
                Self {
                    interface,
                    name,
                    name_args,
                },
                PythonExpression { expr: atom },
            ),
        )
    }

    pub(crate) fn from_transformed_atom(
        self_: &PyRef<'_, TensorExpression>,
        py: Python<'_>,
        atom: Atom,
    ) -> PyResult<Py<TensorExpression>> {
        if atom == self_.as_super().expr {
            return Py::new(py, (Self::clone(self_), PythonExpression { expr: atom }));
        }
        let interface = if atom.as_view().is_zero() {
            self_.interface.clone()
        } else if InterfaceInference::has_structured_syntax(atom.as_view()) {
            SymbolicTensor::infer(atom.clone())
                .map_err(Self::inference_error)?
                .structure
        } else if self_.interface.canonical().is_scalar() {
            PartialStructure::from_logical_slots([])
        } else {
            return Err(PyValueError::new_err(
                "Tensor transformation removed the tensor syntax of a non-scalar expression",
            ));
        };
        // Descriptor identity depends on contractions of the inferred raw ports,
        // even when the input and final result happen to have the same rank.
        let (name, args) = Self::transformed_descriptor(self_, &atom);
        Self::from_known_parts(py, atom, interface, name, args)
    }

    /// Wrap a tensor identity using the shared callback-aware interface proof.
    fn from_preserved_atom(
        self_: &PyRef<'_, Self>,
        py: Python<'_>,
        atom: Atom,
    ) -> PyResult<Py<Self>> {
        if atom == self_.as_super().expr {
            return Py::new(py, (Self::clone(self_), PythonExpression { expr: atom }));
        }
        let value = Self::structured(self_)
            .with_rewritten_expression(atom)
            .map_err(Self::inference_error)?;
        let atom = value.expression;
        let (name, args) = Self::transformed_descriptor(self_, &atom);
        Self::from_parts_unchecked(py, atom, value.structure, name, args)
    }

    /// Wrap scalar algebra using the shared tensor-interface proof.
    fn from_algebra_atom(
        self_: &PyRef<'_, Self>,
        py: Python<'_>,
        atom: Atom,
    ) -> PyResult<Py<Self>> {
        if atom == self_.as_super().expr {
            return Py::new(py, (Self::clone(self_), PythonExpression { expr: atom }));
        }
        let value = Self::structured(self_)
            .with_algebra_result(atom)
            .map_err(Self::inference_error)?;
        let (name, args) = Self::transformed_descriptor(self_, &value.expression);
        Self::from_parts_unchecked(py, value.expression, value.structure, name, args)
    }

    fn with_cooked_indices(
        self_: &PyRef<'_, Self>,
        py: Python<'_>,
        settings: &idenso::CookSettings,
    ) -> PyResult<Py<Self>> {
        let value = Self::structured(self_)
            .with_cooked_indices(settings)
            .map_err(|error| CookingError::new_err(error.to_string()))?;
        let (name, args) = if value.rank() == self_.interface.canonical().order() {
            Self::transformed_descriptor(self_, &value.expression)
        } else {
            (None, Vec::new())
        };
        Self::from_parts_unchecked(py, value.expression, value.structure, name, args)
    }

    fn transformed_descriptor(self_: &PyRef<'_, Self>, atom: &Atom) -> (Option<Symbol>, Vec<Atom>) {
        let original = inferred_descriptor(self_.as_super().expr.as_view());
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
        if atom == self_.as_super().expr {
            return Self::from_preserved_atom(self_, py, atom);
        }
        let (name, args) = Self::transformed_descriptor(self_, &atom);
        let value = Self::structured(self_)
            .with_checked_expression(atom)
            .map_err(Self::inference_error)?;
        let (name, args) = if value.rank() == self_.interface.canonical().order() {
            (name, args)
        } else {
            (None, Vec::new())
        };
        Self::from_parts_unchecked(py, value.expression, value.structure, name, args)
    }

    pub(crate) fn from_structured(
        py: Python<'_>,
        value: SymbolicTensor<PartialStructure>,
    ) -> PyResult<Py<Self>> {
        Self::from_known_parts(py, value.expression, value.structure, None, Vec::new())
    }

    fn promoted_network(self_: &PyRef<'_, Self>, py: Python<'_>) -> PyResult<SpensoNet> {
        let expression = Self::from_known_parts(
            py,
            self_.as_super().expr.clone(),
            self_.interface.clone(),
            self_.name,
            self_.name_args.clone(),
        )?;
        SpensoNet::from_arithmetic(ArithmeticStructure::Tensor(expression), None)
    }

    pub(crate) fn structured(self_: &PyRef<'_, Self>) -> SymbolicTensor<PartialStructure> {
        SymbolicTensor::from_normalized_parts(
            self_.as_super().expr.clone(),
            self_.interface.clone(),
        )
    }

    pub(crate) fn index_replacements(
        interface: &PartialStructure,
        indices: &Bound<'_, PyTuple>,
        cook_indices: Option<&crate::simplification::PyCookSettings>,
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
            .expression)
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
            .expression;
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
        Self::from_known_parts(
            py,
            value.expression,
            value.structure,
            structure.canonical().name(),
            structure.canonical().args().unwrap_or_default(),
        )
    }
}

impl ModuleInit for TensorExpression {}

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
    pub(crate) fn extract(value: &Bound<'_, PyAny>) -> PyResult<Self> {
        if let Ok(value) = value.extract::<PyRef<'_, TensorExpression>>() {
            return Ok(Self::Structured(TensorExpression::structured(&value)));
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
    cook: Option<&PyCookSettings>,
) -> PyResult<AbstractIndex> {
    match value {
        ConvertibleToAbstractIndex::Aind(index) => Ok(index),
        ConvertibleToAbstractIndex::Atom(expression) => {
            let atom = match (cook, expression.expr.as_view()) {
                (Some(settings), AtomView::Fun(fun)) => {
                    settings.rust().cook_function(fun).map_err(|error| {
                        CookingError::new_err(format!("cannot cook index: {error:?}"))
                    })?
                }
                _ => expression.expr,
            };
            atom.as_view().try_into().map_err(|error| {
                let hint = if cook.is_some() {
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
#[pymethods]
impl TensorExpression {
    /// Construct a tensor expression, preserving existing tensor metadata or inferring it.
    ///
    /// Pass `cook_indices=CookSettings.indices()` to flatten nested index payloads before
    /// inferring the tensor interface for arbitrary nested index expressions.
    /// Custom settings control the index encoding, source filters, and output tags.
    #[new]
    #[pyo3(signature = (expression, *, cook_indices = None))]
    #[gen_stub(override_return_type(type_repr = "TensorExpression"))]
    fn new(
        py: Python<'_>,
        expression: &Bound<'_, PyAny>,
        cook_indices: Option<&PyCookSettings>,
    ) -> PyResult<(TensorExpression, PythonExpression)> {
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
        Ok((
            tensor.clone(),
            PythonExpression {
                expr: tensor.as_super().expr.clone(),
            },
        ))
    }

    /// Construct without tensor inference, validation, or contraction normalization.
    ///
    /// Copy the supplied TensorStructure exactly, including logical slot order,
    /// identity, and arguments. The caller must ensure it matches the expression;
    /// subsequent tensor operations trust it. Prefer TensorExpression(expression)
    /// when that correspondence has not already been established.
    #[staticmethod]
    #[pyo3(signature = (expression, *, structure))]
    pub fn unsafe_from_expression(
        py: Python<'_>,
        expression: ConvertibleToExpression,
        structure: &crate::metadata::SpensoTensorStructure,
    ) -> PyResult<Py<Self>> {
        Self::from_parts_unchecked(
            py,
            expression.to_expression().expr,
            structure.interface.clone(),
            structure.name,
            structure.arguments.clone(),
        )
    }

    /// Create an unresolved metric with ports in `rep` and `other`.
    ///
    /// `other` defaults to `rep`; it may also be the dual of the same space.
    /// For example, `g(fund, fund.dual())("i", "j")` creates a fundamental identity.
    /// Call the result with two indices to fill its ports in logical order.
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
    #[staticmethod]
    fn flat(py: Python<'_>, rep: &SpensoRepresentation) -> PyResult<Py<TensorExpression>> {
        let structure = ExplicitKey::<AbstractIndex>::from_iter(
            [rep.representation, rep.representation],
            ETS.flat,
            None,
        );
        Self::from_structure(py, &structure)
    }

    /// Create an unresolved gamma matrix with Minkowski dimension
    /// `minkowski_dimension`.
    ///
    /// Its public ports use storage order: bispinor, bispinor, then Minkowski.
    /// The bispinor dimension is four. Call the result with `(i, j, mu)` to
    /// index it.
    #[staticmethod]
    // Tensor semantics intentionally specialize the inherited Expression API.
    #[gen_stub(type_ignore = ["override", "invalid-method-override"])]
    fn gamma(
        py: Python<'_>,
        minkowski_dimension: ConvertibleToDimension,
    ) -> PyResult<Py<TensorExpression>> {
        Self::from_structure(py, &AGS.gamma_strct::<AbstractIndex>(minkowski_dimension.0))
    }

    /// Create an unresolved gamma-five matrix with two bispinor ports of
    /// `spinor_dimension`; call the result with `(i, j)`.
    #[staticmethod]
    fn gamma5(
        py: Python<'_>,
        spinor_dimension: ConvertibleToDimension,
    ) -> PyResult<Py<TensorExpression>> {
        Self::from_structure(py, &AGS.gamma5_strct::<AbstractIndex>(spinor_dimension.0))
    }

    /// Create an unresolved left-chiral projector with two bispinor ports of
    /// `spinor_dimension`; call the result with `(i, j)`.
    #[staticmethod]
    fn projm(
        py: Python<'_>,
        spinor_dimension: ConvertibleToDimension,
    ) -> PyResult<Py<TensorExpression>> {
        Self::from_structure(py, &AGS.projm_strct::<AbstractIndex>(spinor_dimension.0))
    }

    /// Create an unresolved right-chiral projector with two bispinor ports of
    /// `spinor_dimension`; call the result with `(i, j)`.
    #[staticmethod]
    fn projp(
        py: Python<'_>,
        spinor_dimension: ConvertibleToDimension,
    ) -> PyResult<Py<TensorExpression>> {
        Self::from_structure(py, &AGS.projp_strct::<AbstractIndex>(spinor_dimension.0))
    }

    /// Create an unresolved sigma matrix with Minkowski dimension
    /// `minkowski_dimension`.
    ///
    /// Its logical ports are two Minkowski ports followed by two four-dimensional
    /// bispinor ports. Call the result with `(mu, nu, i, j)` to index it.
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
    /// >>> from symbolica.community.spenso import TensorExpression
    /// >>> generator = TensorExpression.t(8, 3)
    /// >>> indexed = generator("a", "i", "j")
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
    #[getter]
    fn rank(&self) -> usize {
        self.interface.canonical().order()
    }

    /// Whether the expression has no external tensor ports.
    #[getter]
    fn is_scalar(&self) -> bool {
        self.interface.canonical().is_scalar()
    }

    /// Tensor identity, scalar arguments, and external ports in logical order.
    #[getter]
    #[gen_stub(override_return_type(type_repr = "tuple[int | Expression, ...]"))]
    fn shape(&self, py: Python<'_>) -> PyResult<Py<PyTuple>> {
        self.structure().shape(py)
    }

    /// Tensor identity, scalar arguments, and external ports in logical order.
    #[getter]
    fn structure(&self) -> crate::metadata::SpensoTensorStructure {
        crate::metadata::SpensoTensorStructure {
            interface: self.interface.clone(),
            name: self.name,
            arguments: self.name_args.clone(),
        }
    }

    /// The optional identity used when this expression describes stored data.
    #[getter]
    fn name(&self) -> Option<SpensoName> {
        self.name.map(|name| SpensoName { name })
    }

    /// Return a named descriptor without changing the underlying expression.
    fn with_name(
        self_: PyRef<'_, Self>,
        py: Python<'_>,
        name: ConvertibleToSpensoName,
    ) -> PyResult<Py<Self>> {
        let ConvertibleToSpensoName(name, args) = name;
        Self::from_known_parts(
            py,
            self_.as_super().expr.clone(),
            self_.interface.clone(),
            Some(name.name),
            args,
        )
    }

    fn __len__(&self) -> PyResult<usize> {
        logical_size(&logical_dimensions(&self.interface)?)
    }

    #[gen_stub(skip)]
    fn __getitem__(&self, py: Python<'_>, item: SliceOrIntOrExpanded) -> PyResult<Py<PyAny>> {
        let dimensions = logical_dimensions(&self.interface)?;
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
    fn to_expression(self_: PyRef<'_, Self>) -> PythonExpression {
        PythonExpression {
            expr: self_.as_super().expr.clone(),
        }
    }

    /// Re-parse the underlying symbolic expression and rebuild its ordered tensor interface.
    fn reinfer(self_: PyRef<'_, Self>, py: Python<'_>) -> PyResult<Py<Self>> {
        Self::from_atom_interface(py, self_.as_super().expr.clone(), None)
    }

    /// Replace scalar coefficients or tensor factors and validate the external interface.
    ///
    /// Supports Symbolica's patterns, callbacks, conditions and traversal options.
    /// Index-changing rewrites must use reindex/rename_indices, or explicitly drop
    /// the tensor interface with to_expression() before rebuilding it.
    #[allow(clippy::too_many_arguments)]
    #[pyo3(signature = (pattern, rhs, cond=None, non_greedy_wildcards=None, min_level=0, max_level=None, level_is_tree_depth=false, partial=true, allow_new_wildcards_on_rhs=false, rhs_cache_size=None, repeat=false, once=false, bottom_up=false, nested=false))]
    fn replace(
        self_: PyRef<'_, Self>,
        py: Python<'_>,
        pattern: ConvertibleToExpression,
        rhs: ConvertibleToReplaceWith,
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
        let result = self_.as_super().replace(
            pattern,
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
        Self::preserving_interface(&self_, py, result.expr)
    }

    /// Apply a TensorRule, or construct one from a pattern and right-hand side.
    ///
    /// Traverses sums and products, leaving scalar metadata and tensor arguments
    /// opaque. Each replacement must preserve the factor's explicit interface
    /// without introducing additional index occurrences. Unresolved ports and
    /// unsupported tensor powers are rejected. Conditions are evaluated at each
    /// match; repeated right-hand sides and their checks share a bounded cache.
    /// Remaining user normalization or evaluation hooks are rejected; callbacks
    /// may produce an admitted normalized RHS. Use `rhs_cache_size=0` when RHS
    /// callbacks have side effects.
    #[pyo3(signature = (pattern, rhs=None, *, cond=None, rhs_cache_size=None))]
    fn replace_tensor(
        self_: PyRef<'_, Self>,
        py: Python<'_>,
        pattern: Bound<'_, PyAny>,
        rhs: Option<ConvertibleToReplaceWith>,
        cond: Option<ReplacementCondition>,
        rhs_cache_size: Option<usize>,
    ) -> PyResult<Py<Self>> {
        let borrowed = pattern.extract::<PyRef<'_, crate::tensor_rule::PyTensorRule>>();
        let prepared;
        let rule = match &borrowed {
            Ok(rule) => {
                if rhs.is_some() || cond.is_some() || rhs_cache_size.is_some() {
                    return Err(PyTypeError::new_err(
                        "a TensorRule already owns its RHS, conditions and cache settings",
                    ));
                }
                &rule.rule
            }
            Err(_) => {
                prepared = crate::tensor_rule::PyTensorRule::new(
                    pattern.extract()?,
                    rhs.ok_or_else(|| {
                        PyTypeError::new_err("a pattern requires a replacement RHS")
                    })?,
                    cond,
                    rhs_cache_size.unwrap_or(100),
                )?;
                &prepared.rule
            }
        };
        let value = Self::structured(&self_)
            .replace(rule)
            .map_err(Self::inference_error)?;
        let (name, args) = Self::transformed_descriptor(&self_, &value.expression);
        Self::from_parts_unchecked(py, value.expression, value.structure, name, args)
    }

    /// Apply simultaneous Symbolica replacements and validate the resulting tensor.
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

    /// Expand scalar algebra while preserving and validating the tensor interface.
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
        if value.expression == self_.as_super().expr {
            return Py::new(
                py,
                (
                    Self::clone(&self_),
                    PythonExpression {
                        expr: value.expression,
                    },
                ),
            );
        }
        let (name, args) = Self::transformed_descriptor(&self_, &value.expression);
        Self::from_parts_unchecked(py, value.expression, value.structure, name, args)
    }

    /// Factor scalar algebra while preserving and validating the tensor interface.
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

    /// Collect terms while preserving the tensor interface. Callback results are validated.
    #[pyo3(signature = (*x, key_map = None, coeff_map = None))]
    fn collect(
        self_: PyRef<'_, Self>,
        py: Python<'_>,
        x: &Bound<'_, PyTuple>,
        key_map: Option<Py<PyAny>>,
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

    /// Distribute numerical factors while preserving the tensor interface.
    fn expand_num(self_: PyRef<'_, Self>, py: Python<'_>) -> PyResult<Py<Self>> {
        let result = self_.as_super().expand_num();
        Self::from_preserved_atom(&self_, py, result.expr)
    }

    /// Extract common numerical factors while preserving the tensor interface.
    fn collect_num(self_: PyRef<'_, Self>, py: Python<'_>) -> PyResult<Py<Self>> {
        let result = self_.as_super().collect_num();
        Self::from_preserved_atom(&self_, py, result.expr)
    }

    /// Collect common factors from nested sums while preserving the tensor interface.
    fn collect_factors(self_: PyRef<'_, Self>, py: Python<'_>) -> PyResult<Py<Self>> {
        let result = self_.as_super().collect_factors();
        Self::from_algebra_atom(&self_, py, result.expr)
    }

    /// Group terms with equal numerical coefficients while preserving the tensor interface.
    fn collect_by_coefficient(self_: PyRef<'_, Self>, py: Python<'_>) -> PyResult<Py<Self>> {
        let result = self_.as_super().collect_by_coefficient();
        Self::from_preserved_atom(&self_, py, result.expr)
    }

    /// Combine scalar denominators while preserving the tensor interface.
    fn together(self_: PyRef<'_, Self>, py: Python<'_>) -> PyResult<Py<Self>> {
        let result = self_.as_super().together()?;
        Self::from_algebra_atom(&self_, py, result.expr)
    }

    /// Cancel common numerator and denominator factors while preserving the tensor interface.
    fn cancel(self_: PyRef<'_, Self>, py: Python<'_>) -> PyResult<Py<Self>> {
        let result = self_.as_super().cancel()?;
        Self::from_algebra_atom(&self_, py, result.expr)
    }

    /// Collect powers of a symbol while preserving the tensor interface.
    /// Callback results are validated, as for `collect`.
    #[pyo3(signature = (x, key_map = None, coeff_map = None))]
    fn collect_symbol(
        self_: PyRef<'_, Self>,
        py: Python<'_>,
        x: PythonExpression,
        key_map: Option<Py<PyAny>>,
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

    /// Rewrite in Horner form while preserving the tensor interface.
    /// Omit `vars` to choose the variable order heuristically.
    #[pyo3(signature = (vars = None))]
    fn collect_horner(
        self_: PyRef<'_, Self>,
        py: Python<'_>,
        vars: Option<Vec<PythonExpression>>,
    ) -> PyResult<Py<Self>> {
        let result = self_.as_super().collect_horner(vars)?;
        Self::from_algebra_atom(&self_, py, result.expr)
    }

    /// Decompose scalar denominators into partial fractions while preserving the tensor interface.
    /// Omit `variables` to decompose in all variables.
    #[pyo3(signature = (*variables))]
    fn apart(
        self_: PyRef<'_, Self>,
        py: Python<'_>,
        variables: &Bound<'_, PyTuple>,
    ) -> PyResult<Py<Self>> {
        // Zero has no partial fractions; Symbolica's multivariate path cannot
        // lift an empty numerator polynomial into its coefficient ring.
        if self_.as_super().expr.is_zero() {
            return Self::__copy__(self_, py);
        }
        let result = self_.as_super().apart(variables)?;
        Self::from_algebra_atom(&self_, py, result.expr)
    }

    /// Copy the expression together with its ordered interface and data identity.
    fn __copy__(self_: PyRef<'_, Self>, py: Python<'_>) -> PyResult<Py<Self>> {
        Self::from_known_parts(
            py,
            self_.as_super().expr.clone(),
            self_.interface.clone(),
            self_.name,
            self_.name_args.clone(),
        )
    }

    /// Replace four-dimensional Lorentz slots and compact representations by dimension `D`.
    ///
    /// Spinor and color dimensions, scalar coefficients, and other Lorentz dimensions
    /// are unchanged. Apply this before contracting four-dimensional Lorentz indices.
    fn with_lorentz_dimension(
        self_: PyRef<'_, Self>,
        py: Python<'_>,
        dimension: ConvertibleToDimension,
    ) -> PyResult<Py<Self>> {
        let value = Self::structured(&self_)
            .with_lorentz_dimension(dimension.0)
            .map_err(Self::inference_error)?;
        let (name, args) = if value.rank() == self_.interface.canonical().order() {
            (self_.name, self_.name_args.clone())
        } else {
            (None, Vec::new())
        };
        Self::from_parts_unchecked(py, value.expression, value.structure, name, args)
    }

    /// Apply the selected algebra identities while retaining the ordered external interface.
    ///
    /// Apply selected algebra passes to a fixed point and return a tensor expression.
    ///
    /// The default contracts metrics only. Use SimplifySettings.hep() for gamma,
    /// color, and epsilon algebra too, or construct explicit settings. Dimensions
    /// are taken from slots; dimension changes remain an explicit operation.
    /// No full polynomial expansion occurs unless settings.expand is true.
    /// Raises ValueError if the selected passes have not stabilized by max_passes.
    #[pyo3(signature=(settings=None))]
    fn simplify(
        self_: PyRef<'_, Self>,
        py: Python<'_>,
        settings: Option<&crate::simplification::PySimplifySettings>,
    ) -> PyResult<Py<Self>> {
        let settings = settings.copied().unwrap_or_default();
        let mut atom = self_.as_super().expr.clone();
        for _ in 0..settings.max_passes {
            let before = atom.clone();
            if settings.metrics {
                atom = atom.simplify_metrics();
            }
            if let Some(gamma) = settings.gamma {
                atom = atom.simplify_gamma_with(gamma.rust());
            }
            if let Some(color) = settings.color {
                atom = atom.simplify_color_with(color.rust());
            }
            if settings.epsilon {
                atom = atom.simplify_epsilon();
            }
            if settings.expand {
                atom = atom.expand();
            }
            if atom == before {
                return if settings.expand {
                    Self::from_algebra_atom(&self_, py, atom)
                } else {
                    Self::from_preserved_atom(&self_, py, atom)
                };
            }
        }
        Err(PyValueError::new_err(format!(
            "simplification did not stabilize within {} passes",
            settings.max_passes
        )))
    }

    /// Simplify registered Spenso gamma chains and traces with Idenso's default rules.
    ///
    /// The dimension-generic part applies compatible Clifford anticommutation, adjacent
    /// contractions, and ordinary odd/even trace recursion. Chisholm identities, gamma-five
    /// anticommutation and traces, gamma-zero conjugation, and chiral-projector rules are applied
    /// only when the expression carries explicit four-dimensional Minkowski and bispinor
    /// representations.
    ///
    /// This method does not select or implement a dimensional-regularization gamma-five scheme.
    /// Its gamma-five rules are strictly four-dimensional, and the default method does
    /// not enable the optional three-gamma epsilon expansion available through `GammaSimplifySettings`.
    /// Gamma factors must use the Spenso representation-aware forms registered on import;
    /// unrecognized plain Symbolica functions are left unchanged.
    ///
    /// # Examples
    /// ```python
    /// >>> from symbolica import E
    /// >>> from symbolica.community.spenso import TensorExpression
    /// >>> # Built-in representations are registered automatically on import.
    /// >>> trace = E('''
    /// ...     gamma(bis(4,a),bis(4,b),mink(4,mu))
    /// ...     * gamma(bis(4,b),bis(4,a),mink(4,nu))
    /// ... ''', default_namespace="spenso")
    /// >>> simplified = TensorExpression(trace).simplify_gamma()
    /// >>> "gamma(" not in str(simplified) and "g(" in str(simplified)
    /// True
    /// ```
    ///
    /// The native gamma argument order is `bis(dim,alpha), bis(dim,beta), mink(dim,mu)`:
    /// `alpha` and `beta` are spinor indices, followed by the Lorentz index `mu`.
    /// These forms can also be constructed through the HEP tensor library.
    ///
    /// # Arguments
    /// - `self`: expression containing gamma matrix products and traces
    ///
    /// # Returns
    /// The simplified expression with gamma algebra applied.
    ///
    /// # Examples:
    /// ```python
    /// from symbolica.community.spenso import TensorLibrary, TensorName
    /// from symbolica.community.spenso import TensorExpression
    /// from symbolica import S, Expression
    /// # Get HEP library with standard tensors
    /// hep_lib = TensorLibrary.hep_lib()
    /// # Access standard tensors like gamma matrices
    /// gamma_structure = hep_lib[S("spenso::gamma")].expression()
    /// print(gamma_structure)
    /// print(TensorExpression(gamma_structure(3, 4, 7) * gamma_structure(7, 4, 3)).simplify_gamma())
    /// ```
    #[pyo3(signature = (settings = None))]
    fn simplify_gamma(
        self_: PyRef<'_, Self>,
        py: Python<'_>,
        settings: Option<&PyGammaSimplifySettings>,
    ) -> PyResult<Py<Self>> {
        let result = match settings {
            Some(settings) => self_.as_super().expr.simplify_gamma_with(settings.rust()),
            None => self_.as_super().expr.simplify_gamma(),
        };
        Self::from_preserved_atom(&self_, py, result)
    }

    /// Convert bispinor tensors into chain/trace shorthands and join adjacent gamma chains.
    fn collect_gamma_chains(self_: PyRef<'_, Self>, py: Python<'_>) -> PyResult<Py<Self>> {
        let result = self_.as_super().expr.collect_gamma_chains();
        Self::from_preserved_atom(&self_, py, result)
    }

    /// Simplify products and linear combinations involving the time-like gamma matrix `gamma0`.
    fn simplify_gamma0(self_: PyRef<'_, Self>, py: Python<'_>) -> PyResult<Py<Self>> {
        let result = self_.as_super().expr.simplify_gamma0();
        Self::from_preserved_atom(&self_, py, result)
    }

    /// Rewrite conjugated Dirac matrices using gamma0 sandwiches and fresh spinor indices.
    ///
    /// Includes ordinary gamma matrices and Hermitian four-dimensional gamma0,
    /// gamma5, and chiral projectors. The special-matrix rules leave other spinor
    /// dimensions unchanged; this does not choose a gamma-five regularization scheme.
    ///
    /// Raises `GammaConjugationError` when the expression cannot be rewritten consistently.
    fn simplify_gamma_conjugate(self_: PyRef<'_, Self>, py: Python<'_>) -> PyResult<Py<Self>> {
        let result = self_
            .as_super()
            .expr
            .simplify_gamma_conj::<AbstractIndex>()
            .map_err(|error| GammaConjugationError::new_err(error.to_string()))?;
        Self::from_preserved_atom(&self_, py, result)
    }

    /// Simplify Levi-Civita/metric contractions and pairs of Levi-Civita tensors.
    ///
    /// Simplify Levi-Civita/metric contractions and pairs of Levi-Civita tensors to a fixed point.
    fn simplify_epsilon(self_: PyRef<'_, Self>, py: Python<'_>) -> PyResult<Py<Self>> {
        let result = self_.as_super().expr.simplify_epsilon();
        Self::from_preserved_atom(&self_, py, result)
    }

    /// Contract metric and identity tensors while retaining the ordered external interface.
    ///
    /// Simplifies contractions involving metric tensors and identity tensors.
    ///
    /// Applies fundamental tensor algebra rules for metric and identity tensors:
    ///
    /// **Metric tensor rules:**
    /// - `gᵘᵛ pᵥ → pᵘ` (index raising/lowering)
    /// - `gᵘᵛ gᵥρ → gᵘρ` or `δᵘρ` (metric composition)
    /// - `gᵘᵤ → D` (dimension of spacetime)
    /// - `ηᵘᵛ pᵥ → pᵘ` (flat metric contractions)
    ///
    /// **Identity tensor rules:**
    /// - `δᵘᵛ pᵥ → pᵘ` (Kronecker delta contraction)
    /// - `δᵘᵤ → D` (trace of identity)
    ///
    /// The function recognizes metrics as `spenso::g(...)`
    ///
    /// # Arguments
    /// - `self`: expression containing metric/identity tensor contractions
    ///
    /// # Returns
    /// The simplified expression with metric rules applied.
    ///
    /// # Examples:
    /// ```python
    /// from symbolica.community.spenso import TensorExpression
    /// from symbolica.community.spenso import Representation, TensorExpression, TensorName
    /// q = TensorName("q")
    /// rep = Representation.euc(3)
    /// g = TensorExpression.g(rep)
    /// # With slots (creates TensorExpression)
    /// mu = rep("mu")
    /// nu = rep("nu")
    /// print(TensorExpression(g('mu', 'nu') * q(mu)).simplify_metrics())
    /// ```
    fn simplify_metrics(self_: PyRef<'_, Self>, py: Python<'_>) -> PyResult<Py<Self>> {
        let result = self_.as_super().expr.simplify_metrics();
        Self::from_preserved_atom(&self_, py, result)
    }

    /// Apply Idenso's SU(N) color-algebra simplifier while retaining the external interface.
    ///
    /// Simplify registered Spenso color chains, traces, generators, and structure constants.
    ///
    /// With the default Python settings, the simplifier evaluates supported closed traces and
    /// expands contractions between generators on separate fundamental chains. Its normalization
    /// conventions are
    ///
    /// - `Tr(T^a T^b) = TR δ^{ab}`;
    /// - `Σ_a (T^a)_i^j (T^a)_k^l = TR (δ_i^l δ_k^j - δ_i^j δ_k^l/Nc)`;
    /// - `Σ_a (T^a)_i^j (T^a)_j^k = CF δ_i^k`;
    /// - `Σ_{c,d} f^{acd} f^{bcd} = CA δ^{ab}`.
    /// Antisymmetry and Jacobi identities apply to the registered structure constants.
    ///
    /// `CA = Nc`, `CF = (Nc² - 1)/(2Nc)`, and `TR = 1/2` are the conventional fundamental
    /// SU(Nc) specialization, not identities imposed on every input. The default simplifier keeps
    /// representation invariants symbolic where possible; explicit dimension substitution is a
    /// separate `ColorSimplifySettings` option.
    ///
    /// # Examples
    /// ```python
    /// >>> from symbolica import E
    /// >>> from symbolica.community.spenso import TensorExpression
    /// >>> # Built-in representations are registered automatically on import.
    /// >>> generators = E('''
    /// ...     t(coad(Nc^2-1,a),cof(Nc,i),dind(cof(Nc,j)))
    /// ...     * t(coad(Nc^2-1,a),cof(Nc,k),dind(cof(Nc,l)))
    /// ... ''', default_namespace="spenso")
    /// >>> simplified = TensorExpression(generators).simplify_color()
    /// >>> "t(" not in str(simplified) and "g(" in str(simplified)
    /// True
    /// ```
    ///
    /// **Representation invariants:**
    /// Use `Representation.dimension`, `.casimir()`, `.dynkin_index()`, and `.gram(...)`
    /// to construct the scalar invariants associated with explicitly typed color structures.
    ///
    /// # Arguments
    /// - `self`: expression containing SU(N) color structures
    ///
    /// # Returns
    /// The simplified expression, reduced to representation-owned scalar invariants when possible.
    /// Unsupported or open indexed structures may remain explicitly in the result; their presence
    /// is not an error.
    ///
    /// # Notes
    /// Only representation-aware Spenso color forms are recognized. Plain Symbolica functions with
    /// similar names are left unchanged.
    #[pyo3(signature = (settings = None))]
    fn simplify_color(
        self_: PyRef<'_, Self>,
        py: Python<'_>,
        settings: Option<&PyColorSimplifySettings>,
    ) -> PyResult<Py<Self>> {
        let result = match settings {
            Some(settings) => self_.as_super().expr.simplify_color_with(settings.rust()),
            None => self_.as_super().expr.simplify_color(),
        };
        Self::from_preserved_atom(&self_, py, result)
    }

    /// Factor around tensors carrying fundamental, antifundamental, or adjoint color.
    ///
    /// Factor an expression around tensors carrying fundamental, antifundamental, or adjoint color.
    fn collect_color(self_: PyRef<'_, Self>, py: Python<'_>) -> PyResult<Py<Self>> {
        let result = self_.as_super().expr.collect_color();
        Self::from_preserved_atom(&self_, py, result)
    }

    /// Factor around recognized scalar color invariants such as Casimirs and indices.
    ///
    /// Factor an expression around recognized scalar color invariants such as Casimirs and indices.
    fn collect_color_constants(self_: PyRef<'_, Self>, py: Python<'_>) -> PyResult<Py<Self>> {
        let result = self_.as_super().expr.collect_color_constants();
        Self::from_preserved_atom(&self_, py, result)
    }

    /// Rewrite supplied color dimensions and invariants into a representation-aware Casimir basis.
    ///
    /// Rewrite the supplied color-representation dimensions and invariants into a Casimir basis.
    ///
    /// `fundamental` and `adjoint` are registered Spenso representations supplied by the
    /// tensor-construction API. Only scalar coefficient positions are rewritten.
    #[pyo3(signature = (*, fundamental, adjoint, settings = None))]
    fn to_color_casimir(
        self_: PyRef<'_, Self>,
        py: Python<'_>,
        fundamental: &SpensoRepresentation,
        adjoint: &SpensoRepresentation,
        settings: Option<&PyColorCasimirSettings>,
    ) -> PyResult<Py<Self>> {
        let fundamental = fundamental.representation.to_symbolic([]);
        let adjoint = adjoint.representation.to_symbolic([]);
        let result = self_.as_super().expr.to_color_casimir_with(
            fundamental.as_view(),
            adjoint.as_view(),
            settings.map_or_else(ColorCasimirSettings::default, PyColorCasimirSettings::rust),
        );
        Self::from_preserved_atom(&self_, py, result)
    }

    /// Replace supported `cof(N)` invariants by explicit dimension formulas.
    ///
    /// Replace supported `cof(N)` Casimir, Dynkin-index, and Gram invariants by dimension formulas.
    fn to_cof_dimension_invariants(self_: PyRef<'_, Self>, py: Python<'_>) -> PyResult<Py<Self>> {
        let result = self_.as_super().expr.to_cof_dimension_invariants();
        Self::from_preserved_atom(&self_, py, result)
    }

    /// Expand around color structures and wrap each scalar coefficient with `symbol`.
    ///
    /// Expand around color structures and wrap each resulting scalar coefficient with `symbol`.
    fn wrap_color(self_: PyRef<'_, Self>, py: Python<'_>, symbol: Symbol) -> PyResult<Py<Self>> {
        let result = self_.as_super().expr.wrap_color(symbol);
        Self::from_transformed_atom(&self_, py, result)
    }

    /// Selectively expand around Minkowski structures into `(structure, coefficient)` pairs.
    ///
    /// Expand products around factors carrying registered Minkowski indices.
    ///
    /// This is a selective symbolic expansion: Minkowski-bearing factors become polynomial
    /// variables while unrelated sectors remain coefficients. It does not substitute explicit
    /// four-vector components or choose a metric signature.
    ///
    /// # Arguments
    /// - `self`: a factorized Spenso-compatible expression.
    ///
    /// # Returns
    /// `(structure, coefficient)` pairs distributed around Minkowski-bearing factors.
    ///
    /// # Examples
    /// ```python
    /// >>> from symbolica.community.spenso import TensorExpression
    /// >>> from symbolica.community.spenso import Representation, TensorName
    /// >>> # Built-in representations are registered automatically on import.
    /// >>> minkowski = Representation.mink(4)
    /// >>> mu, nu = minkowski("mu"), minkowski("nu")
    /// >>> p, q, r = TensorName("p"), TensorName("q"), TensorName("r")
    /// >>> p_mu = p(mu).to_expression()
    /// >>> q_nu, r_nu = q(nu).to_expression(), r(nu).to_expression()
    /// >>> factorized = p_mu * (q_nu + r_nu)
    /// >>> terms = TensorExpression(factorized).expand_mink()
    /// >>> sum(structure * coefficient for structure, coefficient in terms) == p_mu * q_nu + p_mu * r_nu
    /// True
    /// ```
    fn expand_mink(self_: PyRef<'_, Self>) -> Vec<PythonTerm> {
        python_terms(self_.as_super().expr.expand_mink())
    }

    /// Selectively expand around bispinor structures into `(structure, coefficient)` pairs.
    ///
    /// Expand products around factors carrying registered bispinor indices.
    ///
    /// # Arguments
    /// - `self`: a factorized Spenso-compatible expression.
    ///
    /// # Returns
    /// `(structure, coefficient)` pairs distributed around bispinor-bearing factors. No explicit
    /// spinor components are substituted.
    ///
    /// # Examples
    /// ```python
    /// >>> from symbolica.community.spenso import TensorExpression
    /// >>> from symbolica.community.spenso import Representation, TensorName
    /// >>> # Built-in representations are registered automatically on import.
    /// >>> bispinor = Representation.bis(4)
    /// >>> alpha, beta = bispinor("alpha"), bispinor("beta")
    /// >>> u, v, w = TensorName("u"), TensorName("v"), TensorName("w")
    /// >>> u_alpha = u(alpha).to_expression()
    /// >>> v_beta, w_beta = v(beta).to_expression(), w(beta).to_expression()
    /// >>> factorized = u_alpha * (v_beta + w_beta)
    /// >>> terms = TensorExpression(factorized).expand_bis()
    /// >>> sum(structure * coefficient for structure, coefficient in terms) == u_alpha * v_beta + u_alpha * w_beta
    /// True
    /// ```
    fn expand_bis(self_: PyRef<'_, Self>) -> Vec<PythonTerm> {
        python_terms(self_.as_super().expr.expand_bis())
    }

    /// Selectively expand around Minkowski and bispinor structures into factorized pairs.
    ///
    /// Expand products around factors carrying Minkowski or bispinor indices.
    ///
    /// This combines the selection patterns of `expand_mink()` and `expand_bis()` in one
    /// coefficient pass. Other representation families remain in the coefficient sector.
    ///
    /// # Arguments
    /// - `self`: a factorized Spenso-compatible expression.
    ///
    /// # Returns
    /// `(structure, coefficient)` pairs distributed around both selected representation families.
    ///
    /// # Examples
    /// ```python
    /// >>> from symbolica.community.spenso import TensorExpression
    /// >>> from symbolica.community.spenso import Representation, TensorName
    /// >>> # Built-in representations are registered automatically on import.
    /// >>> minkowski, bispinor = Representation.mink(4), Representation.bis(4)
    /// >>> p_mu = TensorName("p")(minkowski("mu")).to_expression()
    /// >>> q_mu = TensorName("q")(minkowski("mu")).to_expression()
    /// >>> u_a = TensorName("u")(bispinor("a")).to_expression()
    /// >>> v_a = TensorName("v")(bispinor("a")).to_expression()
    /// >>> factorized = (p_mu + q_mu) * (u_a + v_a)
    /// >>> expected = p_mu * u_a + p_mu * v_a + q_mu * u_a + q_mu * v_a
    /// >>> terms = TensorExpression(factorized).expand_mink_bis()
    /// >>> sum(structure * coefficient for structure, coefficient in terms) == expected
    /// True
    /// ```
    fn expand_mink_bis(self_: PyRef<'_, Self>) -> Vec<PythonTerm> {
        python_terms(self_.as_super().expr.expand_mink_bis())
    }

    /// Selectively expand around metric tensors into `(structure, coefficient)` pairs.
    ///
    /// Expand products around registered metric tensors.
    ///
    /// This is a structural expansion only. It neither contracts the metrics nor substitutes a
    /// dimension or signature; call `simplify_metrics()` separately for supported contractions.
    ///
    /// # Arguments
    /// - `self`: a factorized Spenso-compatible expression.
    ///
    /// # Returns
    /// `(structure, coefficient)` pairs distributed around metric factors.
    ///
    /// # Examples
    /// ```python
    /// >>> from symbolica.community.spenso import TensorExpression
    /// >>> from symbolica.community.spenso import Representation, TensorExpression
    /// >>> # Built-in representations are registered automatically on import.
    /// >>> minkowski = Representation.mink(4)
    /// >>> metric = TensorExpression.g(minkowski)
    /// >>> g_mn = metric(minkowski("mu"), minkowski("nu")).to_expression()
    /// >>> g_rs = metric(minkowski("rho"), minkowski("sigma")).to_expression()
    /// >>> g_ab = metric(minkowski("alpha"), minkowski("beta")).to_expression()
    /// >>> factorized = g_mn * (g_rs + g_ab)
    /// >>> terms = TensorExpression(factorized).expand_metrics()
    /// >>> sum(structure * coefficient for structure, coefficient in terms) == g_mn * g_rs + g_mn * g_ab
    /// True
    /// ```
    fn expand_metrics(self_: PyRef<'_, Self>) -> Vec<PythonTerm> {
        python_terms(self_.as_super().expr.expand_metrics())
    }

    /// Selectively expand around color structures into `(structure, coefficient)` pairs.
    ///
    /// Expand products around registered color factors.
    ///
    /// Fundamental, antifundamental, adjoint, color-chain, color-trace, and supported invariant
    /// factors form the selected sector. This only distributes the symbolic expression; use
    /// `simplify_color()` separately to apply SU(N) identities.
    ///
    /// # Arguments
    /// - `self`: a factorized Spenso-compatible expression.
    ///
    /// # Returns
    /// `(structure, coefficient)` pairs distributed around color-bearing factors.
    ///
    /// # Examples
    /// ```python
    /// >>> from symbolica.community.spenso import TensorExpression
    /// >>> from symbolica.community.spenso import Representation, TensorExpression
    /// >>> # Built-in representations are registered automatically on import.
    /// >>> adjoint, fundamental = Representation.coad(8), Representation.cof(3)
    /// >>> antifundamental = fundamental.dual()
    /// >>> generator = TensorExpression.t(8, 3)
    /// >>> t_a = generator(
    /// ...     adjoint("a"), fundamental("i"), antifundamental("j")
    /// ... ).to_expression()
    /// >>> t_b = generator(
    /// ...     adjoint("b"), fundamental("k"), antifundamental("l")
    /// ... ).to_expression()
    /// >>> t_c = generator(
    /// ...     adjoint("c"), fundamental("m"), antifundamental("n")
    /// ... ).to_expression()
    /// >>> factorized = t_a * (t_b + t_c)
    /// >>> terms = TensorExpression(factorized).expand_color()
    /// >>> sum(structure * coefficient for structure, coefficient in terms) == t_a * t_b + t_a * t_c
    /// True
    /// ```
    fn expand_color(self_: PyRef<'_, Self>) -> Vec<PythonTerm> {
        python_terms(ColorSimplifier::expand_color(&self_.as_super().expr))
    }

    /// Selectively expand around the supplied Symbolica patterns into factorized pairs.
    ///
    /// Selectively expand around the supplied expression patterns.
    ///
    /// Results retain the Rust API's `(structure, coefficient)` factorization.
    fn expand_in_patterns(
        self_: PyRef<'_, Self>,
        patterns: Vec<PythonExpression>,
    ) -> Vec<PythonTerm> {
        let patterns = patterns
            .iter()
            .map(|pattern| pattern.expr.to_pattern())
            .collect::<Vec<_>>();
        python_terms(self_.as_super().expr.expand_in_patterns(&patterns))
    }

    /// Put all explicit indices in a named scope, preserving the tensor interface.
    ///
    /// Scoped copies contract internally as before, but their indices are distinct
    /// from the original. Alphabet display uses primed labels for scoped indices.
    /// Applying the same outer scope twice is idempotent; different scopes nest.
    ///
    /// >>> conjugate = tensor.dirac_adjoint().wrap_indices(S("bra"))
    fn wrap_indices(self_: PyRef<'_, Self>, py: Python<'_>, header: Symbol) -> PyResult<Py<Self>> {
        let slots = self_.interface.logical_slots().into_iter().map(|mut slot| {
            if let PartialIndex::Explicit(index) = slot.aind {
                slot.aind = PartialIndex::Explicit(index.scoped(header));
            }
            slot
        });
        Self::from_known_parts(
            py,
            self_.as_super().expr.wrap_indices(header),
            PartialStructure::from_logical_slots(slots),
            self_.name,
            self_.name_args.clone(),
        )
    }

    /// Flatten nested representation-index payloads using index cooking by default.
    ///
    /// Transform both the expression and its stored explicit slots, retaining logical order.
    ///
    /// Transforms hierarchical index expressions within tensor function arguments
    /// into simplified, flat symbolic representations. This "cooking" process is
    /// essential for pattern matching, simplification, and computational efficiency
    /// when dealing with complex tensor expressions.
    ///
    /// **Index Cooking Transformation:**
    /// - Nested structure: `mink(4, f(g(h(μ))))` → `mink(4, f_g_h_mu)`
    /// - Function chains: `lorentz(up(mu))` → `lorentz(up_mu)`
    /// - Complex arguments: `tensor(rep(dim,type(idx)))` → `tensor(rep(dim,type_idx))`
    ///
    /// **Scope:**
    /// - Only affects indices appearing as function arguments
    /// - Preserves top-level function structure
    /// # Arguments
    /// - `self`: expression containing complex nested index structures
    ///
    /// # Returns
    /// Expression with flattened, simplified index names.
    ///
    /// # Examples:
    /// ```python
    /// from symbolica import S
    /// from symbolica.community.spenso import CookSettings, Representation, TensorName
    ///
    /// rep = Representation.euc(3)
    /// template = TensorName("T")(rep, rep)
    /// nested = S("outer")(S("mu"))
    /// # Cook an arbitrary payload when filling the tensor's open ports.
    /// tensor = template(nested, "nu", cook_indices=CookSettings.indices())
    /// ```
    #[pyo3(signature = (settings = None))]
    fn cook_indices(
        self_: PyRef<'_, Self>,
        py: Python<'_>,
        settings: Option<&PyCookSettings>,
    ) -> PyResult<Py<Self>> {
        Self::with_cooked_indices(&self_, py, &PyCookSettings::indices_or(settings))
    }

    /// Encode one function call as an ordinary Symbolica expression.
    ///
    /// Convert a single function call into a flattened variable symbol.
    ///
    /// Transforms a function expression with arguments into a single symbolic variable
    /// whose name encodes both the function name and its arguments. This is the
    /// atomic version of `cook_indices()`, operating on individual function calls
    /// rather than complete expressions.
    ///
    /// **Function Cooking Transform:**
    /// - Simple function: `f(a, b)` → `f_a_b`
    /// - Nested arguments: `tensor(rep(mu))` → `tensor_rep_mu`
    /// - Multiple arguments: `gamma(alpha, beta, mu)` → `gamma_alpha_beta_mu`
    /// - Complex names: `my_function(x, y)` → `my_function_x_y`
    ///
    ///
    /// **Constraints:**
    /// - Input must be a single function call (not sum, product, etc.)
    /// - Arguments must be cookable (symbols, numbers, simple functions)
    /// - Cannot cook expressions containing polynomials or complex structures
    ///
    /// # Arguments
    /// - `self`: expression representing a single function call to cook
    ///
    /// # Returns
    /// Expression containing the flattened variable symbol.
    ///
    /// # Raises
    /// `TypeError` if input is not a cookable function or contains invalid argument types.
    ///
    /// # Examples:
    /// ```python
    /// import symbolica as sp
    /// from symbolica.community.spenso import TensorExpression
    ///
    /// # Simple function cooking
    /// f = sp.S('f')
    /// a, b = sp.S('a','b')
    ///
    /// cooked = TensorExpression(f(a, b)).cook_function()
    /// print(cooked)  # Outputs: f_a_b
    /// ```
    #[pyo3(signature = (settings = None))]
    fn cook_function(
        self_: PyRef<'_, Self>,
        settings: Option<&PyCookSettings>,
    ) -> PyResult<PythonExpression> {
        let settings = PyCookSettings::flattened_or(settings);
        self_
            .as_super()
            .expr
            .cook_function_with_settings(&settings)
            .map(Into::into)
            .map_err(|error| CookingError::new_err(format!("cannot cook: {error:?}")))
    }

    /// Wrap only contracted-index payloads with `header` and return an ordinary expression.
    ///
    /// Wrapped payloads remain available for symbolic matching and custom cooking. Restore
    /// tensor inference with `TensorExpression(result, cook_indices=CookSettings.indices())`.
    /// Raises `ValueError` when the expression cannot be parsed as a tensor network.
    ///
    /// Wraps only the dummy (contracted) indices within the expression using a header symbol.
    ///
    /// Similar to `wrap_indices`, but selectively identifies and wraps only contracted
    /// indices (those appearing once upstairs and once downstairs, or twice in a
    /// self-dual representation), leaving external (dangling) indices untouched.
    /// This is crucial for proper index management in tensor calculations.
    ///
    /// Contracted indices are those that:
    /// - Appear in both upper and lower positions (for dualizable reps)
    /// - Appear twice in the same position (for self-dual reps)
    /// - Are summed over (Einstein summation convention)
    ///
    /// # Arguments
    /// - `self`: input expression containing both dummy and free indices
    /// - `header`: symbol to use as wrapper function name for dummy indices only
    ///
    /// # Returns
    /// A new expression with only contracted indices wrapped.
    ///
    /// # Raises
    /// `ValueError` when the expression cannot be parsed as a tensor network.
    ///
    /// # Examples:
    /// ```python
    /// from symbolica.community.spenso import TensorName, Slot, Representation
    /// import symbolica as sp
    /// from symbolica.community.spenso import TensorExpression
    ///
    /// T = TensorName("T")
    /// rep = Representation.euc(3)
    /// # With slots (creates TensorExpression)
    /// mu = rep("mu")
    /// nu = rep("nu")
    /// x = sp.S("x")
    /// tensor_with_args = T(x, mu, nu, nu)  # T(x; mu, nu, nu)
    /// # print(tensor_with_args)
    /// print(tensor_with_args.wrap_dummies(sp.S('wrap')))
    ///
    /// ```
    fn wrap_dummies(self_: PyRef<'_, Self>, header: Symbol) -> PyResult<PythonExpression> {
        self_
            .as_super()
            .expr
            .wrap_dummies::<AbstractIndex>(header)
            .map(Into::into)
            .map_err(|error| PyValueError::new_err(error.to_string()))
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
    /// # Arguments
    /// - `self`: tensor expression to analyze
    ///
    /// # Returns
    /// A list of expressions, each representing a free (dangling) index.
    ///
    /// # Raises
    /// `ValueError` when the expression cannot be parsed as a tensor network.
    ///
    /// # Examples:
    /// ```python
    /// from symbolica.community.spenso import TensorName, Slot, Representation
    /// import symbolica as sp
    /// from symbolica.community.spenso import TensorExpression
    ///
    /// T = TensorName("T")
    /// rep = Representation.euc(3)
    /// # With slots (creates TensorExpression)
    /// mu = rep("mu")
    /// nu = rep("nu")
    /// x = sp.S("x")
    /// tensor_with_args = T(x, mu, nu, nu)  # T(x; mu, nu, nu)
    /// # print(tensor_with_args)
    /// print(tensor_with_args.list_dangling())
    /// ```
    fn list_dangling(self_: PyRef<'_, Self>) -> PyResult<Vec<PythonExpression>> {
        if self_.interface.open_positions().is_empty() {
            return Ok(self_
                .interface
                .logical_slots()
                .into_iter()
                .map(|slot| composition::port_atom(slot).into())
                .collect());
        }
        self_
            .as_super()
            .expr
            .list_dangling::<AbstractIndex>()
            .map(|indices| indices.into_iter().map(Into::into).collect())
            .map_err(|error| PyValueError::new_err(error.to_string()))
    }

    /// Canonically order tensor factors and deterministically rename contracted indices.
    ///
    /// Raises `CanonicalizationError` when the expression cannot be parsed or canonicalized.
    fn canonize(self_: PyRef<'_, Self>, py: Python<'_>) -> PyResult<Py<Self>> {
        let result = self_
            .as_super()
            .expr
            .canonize::<AbstractIndex>(AbstractIndex::Dummy)
            .map_err(|error| CanonicalizationError::new_err(error.to_string()))?;
        Self::from_preserved_atom(&self_, py, result)
    }

    /// Replace nested tensor subexpressions by aliases and return the root plus alias mappings.
    ///
    /// Replace nested tensor subexpressions by generated aliases.
    ///
    /// Returns the rewritten root followed by sorted `(alias, original)` pairs.
    fn alias_subtensors(
        self_: PyRef<'_, Self>,
        tensor_name: &str,
    ) -> (PythonExpression, Vec<(PythonExpression, PythonExpression)>) {
        let (root, aliases) = self_
            .as_super()
            .expr
            .alias_subtensors(tensor_name)
            .into_inner_with_aliases();
        let mut aliases = aliases.into_iter().collect::<Vec<_>>();
        aliases.sort_by_cached_key(|(alias, _)| alias.to_string());
        (
            root.into(),
            aliases
                .into_iter()
                .map(|(alias, original)| (alias.into(), original.into()))
                .collect(),
        )
    }

    /// Complex-conjugate this tensor expression and re-infer its interface.
    ///
    /// Complex-conjugate an expression while keeping unevaluated conjugations explicit.
    fn spenso_conjugate(self_: PyRef<'_, Self>, py: Python<'_>) -> PyResult<Py<Self>> {
        let result = self_.as_super().expr.spenso_conj();
        Self::from_transformed_atom(&self_, py, result)
    }

    /// Complex-conjugate this expression and transpose slots in `representation`.
    ///
    /// Complex-conjugate an expression and transpose tensor slots in `representation`.
    fn conjugate_transpose(
        self_: PyRef<'_, Self>,
        py: Python<'_>,
        representation: &SpensoRepresentation,
    ) -> PyResult<Py<Self>> {
        let result = self_
            .as_super()
            .expr
            .conjugate_transpose(representation.representation.rep);
        Self::from_transformed_atom(&self_, py, result)
    }

    /// Construct the physics-aware Dirac adjoint and re-infer the tensor interface.
    ///
    /// Construct the physics-aware Dirac adjoint of a tensor expression.
    ///
    /// Idenso takes the symbolic complex conjugate, reverses compatible open bispinor chains, and
    /// inserts the registered `gamma0` factors required at dangling bispinor slots.
    /// With `preserve_indices=True`, external labels remain attached to the same physical
    /// amplitude legs instead of exchanging the two endpoints of each open chain. Use this
    /// when summing diagrams with different fermion pairings before squaring an amplitude.
    /// The gamma-zero boundary factors still apply. The default matrix-adjoint convention
    /// exchanges endpoints and is unchanged.
    /// The input must use the representation-aware Spenso forms registered on import.
    /// Raises `DiracAdjointError` when the tensor network does not define a consistent adjoint.
    ///
    /// # Examples
    /// ```python
    /// >>> from symbolica.community.spenso import TensorExpression
    /// >>> from symbolica.community.spenso import Representation, TensorName
    /// >>> # Built-in representations are registered automatically on import.
    /// >>> # Re-importing the module does not require explicit re-registration.
    /// >>> bispinor = Representation.bis(4)
    /// >>> spinor = TensorName("u")(bispinor("alpha")).to_expression()
    /// >>> adjoint = TensorExpression(spinor).dirac_adjoint()
    /// >>> len(TensorExpression(adjoint).list_dangling()) == 1
    /// True
    /// >>> "gamma0" in str(adjoint)
    /// True
    /// ```
    ///
    /// # Arguments
    /// - `self`: a Spenso-compatible tensor expression.
    ///
    /// # Returns
    /// The representation-aware Dirac adjoint.
    #[pyo3(signature = (*, preserve_indices=false))]
    fn dirac_adjoint(
        self_: PyRef<'_, Self>,
        py: Python<'_>,
        preserve_indices: bool,
    ) -> PyResult<Py<Self>> {
        let result = self_
            .as_super()
            .expr
            .dirac_adjoint::<AbstractIndex>(preserve_indices)
            .map_err(|error| DiracAdjointError::new_err(error.to_string()))?;
        Self::from_transformed_atom(&self_, py, result)
    }

    /// Encode this expression using Idenso's reversible cooking format.
    ///
    /// Encode selected functions as symbols using reversible cooking by default.
    ///
    /// Raises `CookingError` when a selected function payload cannot be encoded.
    #[pyo3(signature = (settings = None))]
    fn cook(
        self_: PyRef<'_, Self>,
        settings: Option<&PyCookSettings>,
    ) -> PyResult<PythonExpression> {
        PyCookSettings::reversible_or(settings)
            .try_cook(self_.as_super().expr.as_view())
            .map(Into::into)
            .map_err(|error| CookingError::new_err(format!("cannot cook: {error:?}")))
    }

    /// Restore a reversibly cooked expression and re-infer its tensor interface.
    ///
    /// Restore symbols produced by matching reversible cooking settings.
    #[pyo3(signature = (settings = None))]
    fn uncook(
        self_: PyRef<'_, Self>,
        py: Python<'_>,
        settings: Option<&PyCookSettings>,
    ) -> PyResult<Py<Self>> {
        let result =
            PyCookSettings::reversible_or(settings).uncook(self_.as_super().expr.as_view());
        Self::from_transformed_atom(&self_, py, result)
    }

    /// Simplify tensor shorthands using the configured Schoonschip traversal.
    #[pyo3(signature = (settings = None))]
    fn schoonschip(
        self_: PyRef<'_, Self>,
        py: Python<'_>,
        settings: Option<&PySchoonschipSettings>,
    ) -> PyResult<Py<Self>> {
        let result = match settings {
            Some(settings) => self_
                .as_super()
                .expr
                .schoonschip_with_settings(&settings.rust()),
            None => self_.as_super().expr.schoonschip(),
        };
        Self::from_preserved_atom(&self_, py, result)
    }

    /// Parse and contract this expression as a symbolic tensor network.
    ///
    /// Parse and contract a symbolic tensor network using Schoonschip rules.
    ///
    /// Set `expand_contracted_sums=True` to distribute sums before contracted products are executed.
    ///
    /// Raises `NetworkToolingError` when the expression is not a valid tensor network or contraction
    /// fails.
    #[pyo3(signature = (settings = None, *, expand_contracted_sums = false))]
    fn schoonschip_net(
        self_: PyRef<'_, Self>,
        py: Python<'_>,
        settings: Option<&PySchoonschipSettings>,
        expand_contracted_sums: bool,
    ) -> PyResult<Py<Self>> {
        let mut settings = settings
            .map(PySchoonschipSettings::rust)
            .unwrap_or_else(SchoonschipSettings::default_network);
        settings.expand_contracted_sums |= expand_contracted_sums;
        let result = self_
            .as_super()
            .expr
            .schoonschip_with_net::<false, AbstractIndex>(&settings)
            .map_err(|error| NetworkToolingError::new_err(error.to_string()))?;
        Self::from_preserved_atom(&self_, py, result)
    }

    /// Convert contracted rank-one tensors into compact dot-product notation.
    ///
    /// Converts contracted Lorentz/Minkowski indices into dot product notation.
    ///
    /// Automatically identifies and converts patterns like `p(mink(D, mu)) * q(mink(D, mu))`
    /// into the compact, representation-carrying `dot(p(mink(D)), q(mink(D)))` notation. This
    /// simplification is essential for physics calculations involving four-vectors.
    ///
    /// The function recognizes:
    /// - Contracted vector indices: `pᵘqᵤ → p·q`
    /// - Multiple contractions: `pᵘqᵤrᵛsᵥ → (p·q)(r·s)`
    /// - Self-contractions: `pᵘpᵤ → p²`
    ///
    /// # Arguments
    /// - `self`: expression containing contracted Minkowski vector indices
    ///
    /// # Returns
    /// The expression with vector contractions converted to dot products.
    ///
    /// # Examples:
    /// ```python
    /// from symbolica.community.spenso import TensorExpression
    /// from symbolica.community.spenso import Representation, TensorName
    /// p = TensorName("p")
    /// q = TensorName("q")
    /// rep = Representation.euc(3)
    /// # With slots (creates TensorExpression)
    /// mu = rep("mu")
    /// nu = rep("nu")
    ///
    /// print(TensorExpression(p(mu) * q(mu)).to_dots())
    /// ```
    fn to_dots(self_: PyRef<'_, Self>, py: Python<'_>) -> PyResult<Py<Self>> {
        let result = self_.as_super().expr.to_dots();
        Self::from_preserved_atom(&self_, py, result)
    }

    /// Canonicalize compact dot-product shorthands without expanding them.
    ///
    /// Canonicalize compact dot-product shorthands without expanding their tensor structure.
    fn normalize_dots(self_: PyRef<'_, Self>, py: Python<'_>) -> PyResult<Py<Self>> {
        let result = self_.as_super().expr.normalize_dots();
        Self::from_preserved_atom(&self_, py, result)
    }

    /// Expand dot products into explicit metric and indexed-vector contractions.
    ///
    /// Expand compact dot products into explicit metric and indexed-vector contractions.
    ///
    /// Raises `DotExpansionError` when a dot product does not define a valid tensor contraction.
    fn expand_dots(self_: PyRef<'_, Self>, py: Python<'_>) -> PyResult<Py<Self>> {
        let result = self_
            .as_super()
            .expr
            .expand_dots()
            .map_err(|error| DotExpansionError::new_err(error.to_string()))?;
        Self::from_preserved_atom(&self_, py, result)
    }

    /// Replace metric shorthand such as `g(p(rep), q(rep))` by a compact dot product.
    ///
    /// Replace metric shorthand such as `g(p(rep), q(rep.dual()))` by
    /// `dot(p(rep), q(rep.dual()))`.
    fn metric_shorthand_to_dot(self_: PyRef<'_, Self>, py: Python<'_>) -> PyResult<Py<Self>> {
        let result = self_.as_super().expr.metric_shorthand_to_dot();
        Self::from_preserved_atom(&self_, py, result)
    }

    /// Expand every Idenso tensor shorthand into explicit tensor syntax.
    ///
    /// Raises `NetworkToolingError` when the expression is not a valid tensor network or evaluation fails.
    fn undo_all(self_: PyRef<'_, Self>, py: Python<'_>) -> PyResult<Py<Self>> {
        let result = self_
            .as_super()
            .expr
            .undo_all::<AbstractIndex>()
            .map_err(|error| NetworkToolingError::new_err(error.to_string()))?;
        Self::from_preserved_atom(&self_, py, result)
    }

    /// Expand Schoonschip shorthands while leaving dots, chains, and traces compact.
    ///
    /// Raises `NetworkToolingError` when the expression is not a valid tensor network or evaluation fails.
    fn undo_schoonschip(self_: PyRef<'_, Self>, py: Python<'_>) -> PyResult<Py<Self>> {
        let result = self_
            .as_super()
            .expr
            .undo_schoonschip::<AbstractIndex>()
            .map_err(|error| NetworkToolingError::new_err(error.to_string()))?;
        Self::from_preserved_atom(&self_, py, result)
    }

    /// Expand dot-product shorthands while leaving other shorthands compact.
    ///
    /// Raises `NetworkToolingError` when the expression is not a valid tensor network or evaluation fails.
    fn undo_dots(self_: PyRef<'_, Self>, py: Python<'_>) -> PyResult<Py<Self>> {
        let result = self_
            .as_super()
            .expr
            .undo_dots::<AbstractIndex>()
            .map_err(|error| NetworkToolingError::new_err(error.to_string()))?;
        Self::from_preserved_atom(&self_, py, result)
    }

    /// Expand open-chain shorthands while leaving other shorthands compact.
    ///
    /// Raises `NetworkToolingError` when the expression is not a valid tensor network or evaluation fails.
    fn undo_chain(self_: PyRef<'_, Self>, py: Python<'_>) -> PyResult<Py<Self>> {
        let result = self_
            .as_super()
            .expr
            .undo_chain::<AbstractIndex>()
            .map_err(|error| NetworkToolingError::new_err(error.to_string()))?;
        Self::from_preserved_atom(&self_, py, result)
    }

    /// Expand trace shorthands while leaving other shorthands compact.
    ///
    /// Raises `NetworkToolingError` when the expression is not a valid tensor network or evaluation fails.
    fn undo_trace(self_: PyRef<'_, Self>, py: Python<'_>) -> PyResult<Py<Self>> {
        let result = self_
            .as_super()
            .expr
            .undo_trace::<AbstractIndex>()
            .map_err(|error| NetworkToolingError::new_err(error.to_string()))?;
        Self::from_preserved_atom(&self_, py, result)
    }

    /// Join adjacent open chains for `representation`, retaining the external interface.
    ///
    /// Join adjacent open chains for the supplied representation.
    fn collect_chains(
        self_: PyRef<'_, Self>,
        py: Python<'_>,
        representation: &SpensoRepresentation,
    ) -> PyResult<Py<Self>> {
        let result = self_
            .as_super()
            .expr
            .collect_chains(representation.representation.rep);
        Self::from_preserved_atom(&self_, py, result)
    }

    /// Rewrite tensors with two `representation` slots as open-chain factors.
    ///
    /// Rewrite tensors with two slots in `representation` as explicit open-chain factors.
    fn chainify(
        self_: PyRef<'_, Self>,
        py: Python<'_>,
        representation: &SpensoRepresentation,
    ) -> PyResult<Py<Self>> {
        let result = self_
            .as_super()
            .expr
            .chainify(representation.representation.rep);
        Self::from_preserved_atom(&self_, py, result)
    }

    /// Convert chains whose endpoints coincide into trace shorthands.
    fn normalize_chains(self_: PyRef<'_, Self>, py: Python<'_>) -> PyResult<Py<Self>> {
        let result = self_.as_super().expr.normalize_chains();
        Self::from_preserved_atom(&self_, py, result)
    }

    /// Replace one-factor chain shorthands by their underlying tensor factor.
    fn undo_single_length(self_: PyRef<'_, Self>, py: Python<'_>) -> PyResult<Py<Self>> {
        let result = self_.as_super().expr.undo_single_length();
        Self::from_preserved_atom(&self_, py, result)
    }

    /// Materialize unresolved ports and parse this expression as a tensor network.
    ///
    /// When `library` is omitted, use the built-in four-dimensional HEP and SU(3) library.
    #[pyo3(signature = (library = None))]
    fn to_network(
        self_: PyRef<'_, Self>,
        py: Python<'_>,
        library: Option<&SpensorLibrary>,
    ) -> PyResult<SpensoNet> {
        let expression = Self::from_known_parts(
            py,
            self_.as_super().expr.clone(),
            self_.interface.clone(),
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
    /// Pass `AUTO` to leave a port unresolved. Pass `cook_indices=CookSettings.indices()` to flatten nested
    /// symbolic index payloads before insertion. Repeated compatible indices contract their
    /// ports, in which case the result no longer carries the original stored-data identity.
    #[pyo3(signature = (*indices, cook_indices = None))]
    fn index(
        self_: PyRef<'_, Self>,
        py: Python<'_>,
        indices: &Bound<'_, PyTuple>,
        cook_indices: Option<&crate::simplification::PyCookSettings>,
    ) -> PyResult<Py<Self>> {
        let replacements = Self::index_replacements(&self_.interface, indices, cook_indices)?;
        let value = Self::structured(&self_)
            .indexed(&replacements)
            .map_err(|error| PyValueError::new_err(error.to_string()))?;
        let (name, args) = if value.rank() == self_.interface.canonical().order() {
            (self_.name, self_.name_args.clone())
        } else {
            (None, Vec::new())
        };
        Self::from_parts_unchecked(py, value.expression, value.structure, name, args)
    }

    /// Assign all external ports in logical order, including already indexed ports.
    /// AUTO keeps the current port. Repeated compatible indices contract.
    #[pyo3(signature = (*indices, cook_indices=None))]
    fn reindex(
        self_: PyRef<'_, Self>,
        py: Python<'_>,
        indices: &Bound<'_, PyTuple>,
        cook_indices: Option<&PyCookSettings>,
    ) -> PyResult<Py<Self>> {
        let positions = (0..self_.interface.canonical().order()).collect::<Vec<_>>();
        let replacements =
            Self::port_replacements(&self_.interface, &positions, indices, cook_indices)?;
        let value = Self::structured(&self_)
            .reindex_interface_ports(&replacements)
            .map_err(|e| PyValueError::new_err(e.to_string()))?;
        Self::from_known_parts(
            py,
            value.expression,
            value.structure,
            self_.name,
            self_.name_args.clone(),
        )
    }

    /// Rename external indices simultaneously without changing rank or capturing dummy indices.
    /// Keys are index labels or typed Slots; typed keys disambiguate representations.
    #[pyo3(signature = (mapping, *, cook_indices=None))]
    fn rename_indices(
        self_: PyRef<'_, Self>,
        py: Python<'_>,
        mapping: &Bound<'_, PyDict>,
        cook_indices: Option<&PyCookSettings>,
    ) -> PyResult<Py<Self>> {
        let replacements = Self::named_replacements(&self_.interface, mapping, cook_indices)?;
        let value = Self::structured(&self_)
            .reindex_interface_ports(&replacements)
            .map_err(|e| PyValueError::new_err(e.to_string()))?;
        let result = Self::from_known_parts(
            py,
            value.expression,
            value.structure,
            self_.name,
            self_.name_args.clone(),
        )?;
        if result.borrow(py).interface.canonical().order() != self_.interface.canonical().order() {
            return Err(PyValueError::new_err(
                "renaming would contract external ports; use reindex() to request a contraction",
            ));
        }
        Ok(result)
    }

    /// Return a tensor with external axes in the specified logical order.
    /// Unresolved ports acquire fresh indices so identical representations remain distinct.
    /// Use reindex() to assign preferred labels after permutation.
    fn permute_axes(
        self_: PyRef<'_, Self>,
        py: Python<'_>,
        axes: Vec<usize>,
    ) -> PyResult<Py<Self>> {
        let value = Self::structured(&self_)
            .permuted(&axes)
            .map_err(|error| PyValueError::new_err(error.to_string()))?;
        Self::from_known_parts(
            py,
            value.expression,
            value.structure,
            self_.name,
            self_.name_args.clone(),
        )
    }

    /// Fill the unresolved external ports with `indices` in interface order.
    #[pyo3(signature = (*indices, cook_indices = None))]
    // Tensor semantics intentionally specialize the inherited Expression API.
    #[gen_stub(type_ignore = ["override", "invalid-method-override"])]
    fn __call__(
        self_: PyRef<'_, Self>,
        py: Python<'_>,
        indices: &Bound<'_, PyTuple>,
        cook_indices: Option<&crate::simplification::PyCookSettings>,
    ) -> PyResult<Py<Self>> {
        Self::index(self_, py, indices, cook_indices)
    }

    fn __neg__(self_: PyRef<'_, Self>, py: Python<'_>) -> PyResult<Py<Self>> {
        Self::from_known_parts(
            py,
            Atom::num(-1) * self_.as_super().expr.as_ref(),
            self_.interface.clone(),
            None,
            Vec::new(),
        )
    }

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
        let (right, interface) = match TensorOperand::extract(rhs)? {
            TensorOperand::Structured(right) => {
                return left
                    .try_add(&right)
                    .map_err(|error| PyValueError::new_err(error.to_string()))
                    .and_then(|value| Self::from_structured(py, value))
                    .map(TensorDispatch::Expression);
            }
            TensorOperand::Scalar(right) if right.as_view().is_zero() => {
                return Self::from_known_parts(
                    py,
                    left.expression,
                    left.structure,
                    self_.name,
                    self_.name_args.clone(),
                )
                .map(TensorDispatch::Expression);
            }
            TensorOperand::Scalar(right) if left.is_scalar() => (right, left.structure),
            TensorOperand::Scalar(_) => {
                return Err(PyValueError::new_err(
                    "cannot add a scalar expression to a non-scalar tensor",
                ));
            }
        };
        Self::from_known_parts(
            py,
            left.expression.as_ref() + right.as_ref(),
            interface,
            None,
            Vec::new(),
        )
        .map(TensorDispatch::Expression)
    }

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
        let (right, interface) = match TensorOperand::extract(rhs)? {
            TensorOperand::Structured(right) => {
                return left
                    .try_sub(&right)
                    .map_err(|error| PyValueError::new_err(error.to_string()))
                    .and_then(|value| Self::from_structured(py, value))
                    .map(TensorDispatch::Expression);
            }
            TensorOperand::Scalar(right) if right.as_view().is_zero() => {
                return Self::from_known_parts(
                    py,
                    left.expression,
                    left.structure,
                    self_.name,
                    self_.name_args.clone(),
                )
                .map(TensorDispatch::Expression);
            }
            TensorOperand::Scalar(right) if left.is_scalar() => (right, left.structure),
            TensorOperand::Scalar(_) => {
                return Err(PyValueError::new_err(
                    "cannot subtract a scalar expression from a non-scalar tensor",
                ));
            }
        };
        Self::from_known_parts(
            py,
            left.expression.as_ref() - right.as_ref(),
            interface,
            None,
            Vec::new(),
        )
        .map(TensorDispatch::Expression)
    }

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
        let (left, interface) = match TensorOperand::extract(lhs)? {
            TensorOperand::Structured(left) => {
                return left
                    .try_sub(&right)
                    .map_err(|error| PyValueError::new_err(error.to_string()))
                    .and_then(|value| Self::from_structured(py, value))
                    .map(TensorDispatch::Expression);
            }
            TensorOperand::Scalar(left) if left.as_view().is_zero() => {
                return Self::__neg__(self_, py).map(TensorDispatch::Expression);
            }
            TensorOperand::Scalar(left) if right.is_scalar() => (left, right.structure),
            TensorOperand::Scalar(_) => {
                return Err(PyValueError::new_err(
                    "cannot subtract a non-scalar tensor from a scalar expression",
                ));
            }
        };
        Self::from_known_parts(
            py,
            left.as_ref() - right.expression.as_ref(),
            interface,
            None,
            Vec::new(),
        )
        .map(TensorDispatch::Expression)
    }

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
                left.expression.as_ref() * right.as_ref(),
                left.structure,
                None,
                Vec::new(),
            )
            .map(TensorDispatch::Expression),
        }
    }

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
                .multiply(&right)
                .map_err(|error| PyValueError::new_err(error.to_string()))
                .and_then(|value| Self::from_structured(py, value))
                .map(TensorDispatch::Expression),
            TensorOperand::Scalar(left) => Self::from_known_parts(
                py,
                left.as_ref() * right.expression.as_ref(),
                right.structure,
                None,
                Vec::new(),
            )
            .map(TensorDispatch::Expression),
        }
    }

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
        let value = Self::structured(&self_);
        let rhs = match TensorOperand::extract(rhs)? {
            TensorOperand::Scalar(rhs) => rhs,
            TensorOperand::Structured(rhs) if rhs.is_scalar() => rhs.expression,
            TensorOperand::Structured(_) => {
                return Err(PyTypeError::new_err(
                    "tensor division requires a scalar denominator",
                ));
            }
        };
        Self::from_known_parts(
            py,
            value.expression.as_ref() / rhs.as_ref(),
            value.structure,
            None,
            Vec::new(),
        )
        .map(TensorDispatch::Expression)
    }

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
        let value = Self::structured(&self_);
        if !value.is_scalar() {
            return Err(PyValueError::new_err(
                "a non-scalar tensor cannot be used as a denominator",
            ));
        }
        let (lhs, structure) = match TensorOperand::extract(lhs)? {
            TensorOperand::Scalar(lhs) => (lhs, value.structure),
            TensorOperand::Structured(lhs) => (lhs.expression, lhs.structure),
        };
        Self::from_known_parts(
            py,
            lhs.as_ref() / value.expression.as_ref(),
            structure,
            None,
            Vec::new(),
        )
        .map(TensorDispatch::Expression)
    }

    /// Form an outer product without contracting compatible ports.
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

    /// Contract one selected pair of ordered interface positions.
    #[pyo3(signature = (rhs, *, left, right))]
    #[gen_stub(skip)]
    fn contract(
        self_: PyRef<'_, Self>,
        py: Python<'_>,
        rhs: &Bound<'_, PyAny>,
        left: usize,
        right: usize,
    ) -> PyResult<TensorDispatch> {
        if let Some(rhs) = concrete_network(rhs)? {
            return Self::promoted_network(&self_, py)?
                .contract(ConvertibleToSpensoNet(rhs), left, right)
                .and_then(|network| Py::new(py, network))
                .map(TensorDispatch::Network);
        }
        let value = Self::structured(&self_);
        let TensorOperand::Structured(rhs) = TensorOperand::extract(rhs)? else {
            return Err(PyTypeError::new_err("contract() requires a tensor operand"));
        };
        value
            .contract_ports(&rhs, &[composition::PortPair { left, right }])
            .map_err(|error| PyValueError::new_err(error.to_string()))
            .and_then(|value| Self::from_structured(py, value))
            .map(TensorDispatch::Expression)
    }

    /// Compose two selected `(input, output)` matrix channels.
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
            &value,
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
    #[pyo3(signature = (show_dimensions = None, *, settings = None))]
    fn format_tensor(
        self_: PyRef<'_, Self>,
        show_dimensions: Option<bool>,
        settings: Option<PyRef<'_, display::DisplaySettings>>,
    ) -> PyResult<String> {
        let settings = display::resolved_settings(show_dimensions, settings.as_deref());
        display::validate_plain_source_settings(&settings)?;
        Ok(display::format_structured_settings(
            &Self::structured(&self_),
            &settings,
        ))
    }

    /// Export LaTeX with tensor notation, index alphabets, and model parameter labels.
    /// ``max_line_length`` wraps top-level sums using Symbolica's usual rules.
    #[pyo3(signature = (max_line_length = None))]
    fn to_latex(self_: PyRef<'_, Self>, max_line_length: Option<usize>) -> String {
        display::structured_to_latex(&Self::structured(&self_), false, max_line_length)
    }

    /// Format this structured expression as Typst math source.
    #[pyo3(signature = (show_dimensions = None, *, settings = None))]
    // Tensor semantics intentionally specialize the inherited Expression API.
    #[gen_stub(type_ignore = ["override", "invalid-method-override"])]
    fn to_typst(
        self_: PyRef<'_, Self>,
        show_dimensions: Option<bool>,
        settings: Option<PyRef<'_, display::DisplaySettings>>,
    ) -> PyResult<String> {
        let settings = display::resolved_settings(show_dimensions, settings.as_deref());
        display::validate_typst_source_settings(&settings)?;
        Ok(display::structured_to_typst_with_settings(
            &Self::structured(&self_),
            &settings,
        ))
    }

    #[pyo3(signature = (show_dimensions = None, *, settings = None, notation_source = None))]
    /// Build Symbolica's rich display value, including semantic HTML when the
    /// optional ``gammaloop[typst-display]`` renderer is installed.
    ///
    /// ``notation_source`` is a trusted complete replacement for the bundled
    /// ``notation.typ`` module, not a style fragment. Typst executes it. When
    /// HTML cannot be rendered, the result retains its LaTeX and text forms.
    // Tensor semantics intentionally specialize the inherited Expression API.
    #[gen_stub(type_ignore = ["override", "invalid-method-override"])]
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
            &Self::structured(&self_),
            &settings,
            notation_source.as_deref(),
        )
    }

    #[pyo3(signature = (show_dimensions = None, *, settings = None, notation_source = None))]
    /// Compile this expression to semantic HTML with the optional Typst
    /// renderer.
    ///
    /// Install ``gammaloop[typst-display]`` to enable this method.
    /// ``notation_source``, when supplied, is trusted Typst code replacing the
    /// complete bundled ``notation.typ`` module.
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
            &Self::structured(&self_),
            &settings,
            notation_source.as_deref(),
        )
    }

    #[pyo3(signature = (show_dimensions = None, *, settings = None, notation_source = None))]
    /// Compile this expression to SVG with the optional Typst renderer.
    ///
    /// Install ``gammaloop[typst-display]`` to enable this method.
    /// ``notation_source``, when supplied, is trusted Typst code replacing the
    /// complete bundled ``notation.typ`` module.
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
            &Self::structured(&self_),
            &settings,
            notation_source.as_deref(),
        )
    }

    fn __repr__(self_: PyRef<'_, Self>) -> String {
        display::format_structured(&Self::structured(&self_), false)
    }

    fn __str__(self_: PyRef<'_, Self>) -> String {
        display::format_structured(&Self::structured(&self_), false)
    }

    fn _repr_pretty_(
        self_: PyRef<'_, Self>,
        pretty: &Bound<'_, PyAny>,
        cycle: bool,
    ) -> PyResult<()> {
        let text = if cycle {
            "...".to_string()
        } else {
            display::format_structured(&Self::structured(&self_), false)
        };
        pretty.call_method1("text", (text,))?;
        Ok(())
    }

    // Tensor semantics intentionally specialize the inherited Expression API.
    #[gen_stub(type_ignore = ["override", "invalid-method-override"])]
    fn _repr_html_(self_: PyRef<'_, Self>, py: Python<'_>) -> Option<String> {
        let settings = display::DisplaySettings::default();
        display::structured_to_html(py, &Self::structured(&self_), &settings, None).ok()
    }

    fn _repr_latex_(self_: PyRef<'_, Self>) -> String {
        Self::to_latex(self_, None)
    }
}

/// Restore tensor-aware dispatch after a base Symbolica transformation.
#[cfg_attr(
    feature = "python_stubgen",
    pyo3_stub_gen::derive::gen_stub_pyfunction(module = "symbolica.community.spenso")
)]
#[pyfunction]
pub fn as_tensor(py: Python<'_>, expression: &Bound<'_, PyAny>) -> PyResult<Py<TensorExpression>> {
    if let Ok(expression) = expression.extract::<PyRef<'_, TensorExpression>>() {
        return TensorExpression::from_known_parts(
            py,
            expression.as_super().expr.clone(),
            expression.interface.clone(),
            expression.name,
            expression.name_args.clone(),
        );
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

/// Contract two rank-one tensors into the canonical dot form.
#[cfg_attr(
    feature = "python_stubgen",
    pyo3_stub_gen::derive::gen_stub_pyfunction(
        module = "symbolica.community.spenso",
        no_default_overload = true,
        python_overload = r#"
        import typing

        @overload
        def dot(left: typing.Union[Tensor, TensorNetwork], right: pyo3_stub_gen.RustType["ConvertibleToSpensoNet"]) -> pyo3_stub_gen.RustType["SpensoNet"]:
            """Contract two rank-one tensors into the canonical dot form."""
        @overload
        def dot(left: pyo3_stub_gen.RustType["ConvertibleToSpensoNet"], right: typing.Union[Tensor, TensorNetwork]) -> pyo3_stub_gen.RustType["SpensoNet"]: ...
        @overload
        def dot(left: pyo3_stub_gen.RustType["PythonExpression"], right: pyo3_stub_gen.RustType["PythonExpression"]) -> pyo3_stub_gen.RustType["TensorExpression"]: ...
        "#,
    )
)]
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
/// Build an explicitly-ended ordered tensor chain.
#[cfg_attr(
    feature = "python_stubgen",
    pyo3_stub_gen::derive::gen_stub_pyfunction(
        module = "symbolica.community.spenso",
        no_default_overload = true,
        python_overload = r#"
        import typing

        @overload
        def chain(start_slot: pyo3_stub_gen.RustType["SpensoSlot"], end_slot: pyo3_stub_gen.RustType["SpensoSlot"], factor: typing.Union[Tensor, TensorNetwork], /, *factors: pyo3_stub_gen.RustType["ConvertibleToSpensoNet"]) -> pyo3_stub_gen.RustType["SpensoNet"]: ...
        @overload
        def chain(start_slot: pyo3_stub_gen.RustType["SpensoSlot"], end_slot: pyo3_stub_gen.RustType["SpensoSlot"], first: pyo3_stub_gen.RustType["ConvertibleToSpensoNet"], second: typing.Union[Tensor, TensorNetwork], /, *factors: pyo3_stub_gen.RustType["ConvertibleToSpensoNet"]) -> pyo3_stub_gen.RustType["SpensoNet"]: ...
        @overload
        def chain(start_slot: pyo3_stub_gen.RustType["SpensoSlot"], end_slot: pyo3_stub_gen.RustType["SpensoSlot"], *factors: pyo3_stub_gen.RustType["PythonExpression"]) -> pyo3_stub_gen.RustType["TensorExpression"]:
            """Build an explicitly-ended ordered tensor chain."""
        @overload
        def chain(start_slot: pyo3_stub_gen.RustType["SpensoSlot"], end_slot: pyo3_stub_gen.RustType["SpensoSlot"], *factors: pyo3_stub_gen.RustType["ConvertibleToSpensoNet"]) -> pyo3_stub_gen.RustType["TensorDispatch"]: ...
        "#,
    )
)]
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
        let slots = value.structure.structure.logical_slots();
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

/// Close an ordered factor sequence into a canonical cyclic trace.
#[cfg_attr(
    feature = "python_stubgen",
    pyo3_stub_gen::derive::gen_stub_pyfunction(
        module = "symbolica.community.spenso",
        no_default_overload = true,
        python_overload = r#"
        import typing

        @overload
        def trace(representation: pyo3_stub_gen.RustType["SpensoRepresentation"], factor: typing.Union[Tensor, TensorNetwork], /, *factors: pyo3_stub_gen.RustType["ConvertibleToSpensoNet"]) -> pyo3_stub_gen.RustType["SpensoNet"]: ...
        @overload
        def trace(representation: pyo3_stub_gen.RustType["SpensoRepresentation"], first: pyo3_stub_gen.RustType["ConvertibleToSpensoNet"], second: typing.Union[Tensor, TensorNetwork], /, *factors: pyo3_stub_gen.RustType["ConvertibleToSpensoNet"]) -> pyo3_stub_gen.RustType["SpensoNet"]: ...
        @overload
        def trace(representation: pyo3_stub_gen.RustType["SpensoRepresentation"], *factors: pyo3_stub_gen.RustType["PythonExpression"]) -> pyo3_stub_gen.RustType["TensorExpression"]:
            """Close an ordered factor sequence into a canonical cyclic trace."""
        @overload
        def trace(representation: pyo3_stub_gen.RustType["SpensoRepresentation"], *factors: pyo3_stub_gen.RustType["ConvertibleToSpensoNet"]) -> pyo3_stub_gen.RustType["TensorDispatch"]: ...
        "#,
    )
)]
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
        let input = value.structure.structure.logical_slots()[channel.input];
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

// The stub generator's Python parser requires typing.Union instead of `|`.
#[cfg(feature = "python_stubgen")]
submit! {
    pyo3_stub_gen::derive::gen_methods_from_python! {
        r#"
        import typing

        class TensorExpression:
            @overload
            def __add__(self, rhs: typing.Union[Tensor, TensorNetwork]) -> pyo3_stub_gen.RustType["SpensoNet"]: ...
            @overload
            def __add__(self, rhs: pyo3_stub_gen.RustType["ConvertibleToExpression"]) -> pyo3_stub_gen.RustType["TensorExpression"]: ...
            @overload
            def __radd__(self, lhs: typing.Union[Tensor, TensorNetwork]) -> pyo3_stub_gen.RustType["SpensoNet"]: ...
            @overload
            def __radd__(self, lhs: pyo3_stub_gen.RustType["ConvertibleToExpression"]) -> pyo3_stub_gen.RustType["TensorExpression"]: ...
            @overload
            def __sub__(self, rhs: typing.Union[Tensor, TensorNetwork]) -> pyo3_stub_gen.RustType["SpensoNet"]: ...
            @overload
            def __sub__(self, rhs: pyo3_stub_gen.RustType["ConvertibleToExpression"]) -> pyo3_stub_gen.RustType["TensorExpression"]: ...
            @overload
            def __rsub__(self, lhs: typing.Union[Tensor, TensorNetwork]) -> pyo3_stub_gen.RustType["SpensoNet"]: ...
            @overload
            def __rsub__(self, lhs: pyo3_stub_gen.RustType["ConvertibleToExpression"]) -> pyo3_stub_gen.RustType["TensorExpression"]: ...
            @overload
            def __mul__(self, rhs: typing.Union[Tensor, TensorNetwork]) -> pyo3_stub_gen.RustType["SpensoNet"]: ...
            @overload
            def __mul__(self, rhs: pyo3_stub_gen.RustType["ConvertibleToExpression"]) -> pyo3_stub_gen.RustType["TensorExpression"]: ...
            @overload
            def __rmul__(self, lhs: typing.Union[Tensor, TensorNetwork]) -> pyo3_stub_gen.RustType["SpensoNet"]: ...
            @overload
            def __rmul__(self, lhs: pyo3_stub_gen.RustType["ConvertibleToExpression"]) -> pyo3_stub_gen.RustType["TensorExpression"]: ...
            @overload
            def __truediv__(self, rhs: typing.Union[Tensor, TensorNetwork]) -> pyo3_stub_gen.RustType["SpensoNet"]: ...
            @overload
            def __truediv__(self, rhs: pyo3_stub_gen.RustType["ConvertibleToExpression"]) -> pyo3_stub_gen.RustType["TensorExpression"]: ...
            @overload
            def __rtruediv__(self, lhs: typing.Union[Tensor, TensorNetwork]) -> pyo3_stub_gen.RustType["SpensoNet"]: ...
            @overload
            def __rtruediv__(self, lhs: pyo3_stub_gen.RustType["ConvertibleToExpression"]) -> pyo3_stub_gen.RustType["TensorExpression"]: ...
            @overload
            def outer(self, rhs: typing.Union[Tensor, TensorNetwork]) -> pyo3_stub_gen.RustType["SpensoNet"]:
                """Form an outer product without contracting compatible ports."""
            @overload
            def outer(self, rhs: pyo3_stub_gen.RustType["PythonExpression"]) -> pyo3_stub_gen.RustType["TensorExpression"]:
                """Form an outer product without contracting compatible ports."""
            @overload
            def contract(self, rhs: typing.Union[Tensor, TensorNetwork], *, left: int, right: int) -> pyo3_stub_gen.RustType["SpensoNet"]:
                """Contract one selected pair of ordered interface positions."""
            @overload
            def contract(self, rhs: pyo3_stub_gen.RustType["PythonExpression"], *, left: int, right: int) -> pyo3_stub_gen.RustType["TensorExpression"]:
                """Contract one selected pair of ordered interface positions."""
            @overload
            def compose(self, rhs: typing.Union[Tensor, TensorNetwork], *, left: tuple[int, int], right: tuple[int, int]) -> pyo3_stub_gen.RustType["SpensoNet"]:
                """Compose two selected `(input, output)` matrix channels."""
            @overload
            def compose(self, rhs: pyo3_stub_gen.RustType["PythonExpression"], *, left: tuple[int, int], right: tuple[int, int]) -> pyo3_stub_gen.RustType["TensorExpression"]:
                """Compose two selected `(input, output)` matrix channels."""
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
                doc: "Convert a logical row-major flat index to tensor coordinates.",
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
                doc: "Convert tensor coordinates to a logical row-major flat index.",
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
                doc: "Expand a slice of logical row-major flat indices to tensor coordinates.",
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
    use symbolica::{atom::FunctionBuilder, symbol};

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
            "contract",
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
                ConvertibleToExpression::type_input()
            } else {
                PythonExpression::type_input()
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
            .filter(|info| info.module == Some("symbolica.community.spenso"));
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
                    && (factors.type_info)() == PythonExpression::type_input()
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
        use idenso::color::ColorSimplifier;

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
            let simplified = numerator.simplify_color();
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
            infer_interface(&scalar.simplify_epsilon())
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
        let settings = SchoonschipSettings::new(None).with_expanded_contracted_sums();
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
                            infer_interface(&product.expression)
                                .unwrap()
                                .canonical()
                                .is_scalar()
                        );
                        assert_eq!(
                            product
                                .expression
                                .schoonschip_with_net::<false, AbstractIndex>(&settings)
                                .unwrap()
                                .to_dots()
                                .expand_in(Atom::var(ETS.metric).as_view())
                                .schoonschip_with_net::<false, AbstractIndex>(&settings)
                                .unwrap()
                                .to_dots(),
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
                    infer_interface(&product.expression).unwrap().canonical(),
                    infer_interface(&expected).unwrap().canonical()
                );
                assert_eq!(
                    product
                        .expression
                        .schoonschip_net::<AbstractIndex>()
                        .unwrap(),
                    expected
                );
            }
        }
    }

    #[test]
    fn raw_gamma_trace_accepts_longitudinal_projection() {
        use idenso::shorthands::schoonschip::{Schoonschip, SchoonschipSettings};

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
            .expression
            .schoonschip_with_net::<false, AbstractIndex>(
                &SchoonschipSettings::default_network().with_expanded_contracted_sums(),
            )
            .unwrap();
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
                    let compact = TensorExpression::schoonschip(tensor.borrow(py), py, None)?;
                    let collected = TensorExpression::collect_gamma_chains(compact.borrow(py), py)?;
                    let collected = collected.borrow(py);
                    let expected = infer_interface(&original)?;
                    assert_eq!(collected.interface.canonical(), expected.canonical());
                    let atom = &collected.as_super().expr;
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
            SPENSO_TAG.trace(wrong_bis.to_symbolic([]), [idenso::gamma!(&compact)]),
        ];
        for atom in invalid {
            assert!(
                infer_interface(&atom).is_err(),
                "invalid compact gamma accepted: {atom}"
            );
        }
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
                "reinfer",
                "expand",
                "simplify_gamma",
                "collect_gamma_chains",
                "simplify_gamma0",
                "simplify_gamma_conjugate",
                "simplify_epsilon",
                "simplify_metrics",
                "simplify_color",
                "collect_color",
                "collect_color_constants",
                "to_color_casimir",
                "to_cof_dimension_invariants",
                "wrap_color",
                "expand_mink",
                "expand_bis",
                "expand_mink_bis",
                "expand_metrics",
                "expand_color",
                "expand_in_patterns",
                "wrap_indices",
                "cook_indices",
                "cook_function",
                "wrap_dummies",
                "list_dangling",
                "canonize",
                "alias_subtensors",
                "spenso_conjugate",
                "conjugate_transpose",
                "dirac_adjoint",
                "cook",
                "uncook",
                "schoonschip",
                "schoonschip_net",
                "to_dots",
                "normalize_dots",
                "expand_dots",
                "metric_shorthand_to_dot",
                "undo_all",
                "undo_schoonschip",
                "undo_dots",
                "undo_chain",
                "undo_trace",
                "collect_chains",
                "chainify",
                "normalize_chains",
                "undo_single_length",
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
                        .interface
                        .logical_slots()
                        .into_iter()
                        .map(|slot| slot.rep())
                        .collect::<Vec<_>>(),
                    *expected
                );
                assert!(tensor.name_args.is_empty());
                assert_eq!(tensor.interface.open_positions().len(), expected.len());
                let AtomView::Fun(function) = tensor.as_super().expr.as_view() else {
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
                    indexed.interface.logical_slots(),
                    logical_slots
                        .iter()
                        .map(|slot| slot.rep().slot(PartialIndex::Explicit(slot.aind())))
                        .collect::<Vec<_>>()
                );
                let AtomView::Fun(function) = indexed.as_super().expr.as_view() else {
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
                        .interface
                        .logical_slots()
                        .iter()
                        .map(IsAbstractSlot::dim)
                        .collect::<Vec<_>>(),
                    expected_dimensions
                );
                assert!(expression.name_args.is_empty());
                let AtomView::Fun(function) = expression.as_super().expr.as_view() else {
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
        let inferred = infer_interface(&generator.expression).unwrap();
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
                value.structure.logical_slots(),
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
            assert_eq!(relabeled.interface.canonical().order(), 2);
            assert_eq!(relabeled.name, Some(name));
            assert_eq!(relabeled.name_args, vec![argument]);
            drop(relabeled);

            let contracted = expression.bind(py).call1(("i", "i"))?;
            let contracted = contracted.extract::<PyRef<'_, TensorExpression>>()?;
            assert!(contracted.interface.canonical().is_scalar());
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
                assert_eq!(restored.as_super().expr, scoped.as_super().expr);
                assert_eq!(restored.interface, scoped.interface);
                assert_eq!(scoped.interface.canonical().order(), rank);
                let cooked = tensor_type.call((&wrapped,), Some(&kwargs))?;
                let cooked = cooked.extract::<Py<TensorExpression>>()?;
                assert_eq!(cooked.borrow(py).interface.canonical().order(), rank);
                if rank == 0 {
                    let simplified = cooked.bind(py).call_method0("simplify_metrics")?;
                    let simplified = simplified.extract::<PyRef<'_, TensorExpression>>()?;
                    assert_eq!(simplified.as_super().expr, Atom::num(3));
                }
                let encoded = indexed.call_method0("cook")?;
                assert!(!encoded.is_instance_of::<TensorExpression>());
                let restored = tensor_type.call1((encoded,))?.call_method0("uncook")?;
                let restored = restored.extract::<PyRef<'_, TensorExpression>>()?;
                let indexed = indexed.extract::<PyRef<'_, TensorExpression>>()?;
                assert_eq!(restored.as_super().expr, indexed.as_super().expr);
                assert_eq!(restored.interface, indexed.interface);
            }
            Ok(())
        })
        .unwrap();
    }

    #[test]
    fn tensor_constructor_preserves_custom_index_encoding_and_tags() {
        use symbolica::atom::UserData;

        idenso::representations::initialize();
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let module = pyo3::types::PyModule::new(py, "spenso")?;
            crate::simplification::register(&module)?;
            let mode = module
                .getattr("CookMode")?
                .getattr("ReversibleEncoding")?
                .extract()?;
            let payload = FunctionBuilder::new(symbol!(
                "constructor_index_source",
                tags = ["idenso::constructor_input_index"]
            ))
            .add_arg(Atom::num(7))
            .finish();
            let representation = Minkowski {}.new_rep(4).cast::<LibraryRep>();
            let raw = FunctionBuilder::new(SPENSO_TAG.tensor_symbol("constructor_cooked_tensor"))
                .add_arg(representation.to_symbolic([payload.clone()]))
                .finish();
            let expression = PythonExpression { expr: raw.clone() }.into_pyobject(py)?;
            let tensor_type = py.get_type::<TensorExpression>();
            for preserve_tags in [false, true] {
                let settings = PyCookSettings::new(
                    mode,
                    None,
                    Some(vec!["idenso::constructor_output_index".into()]),
                    preserve_tags,
                );
                let kwargs = pyo3::types::PyDict::new(py);
                kwargs.set_item("cook_indices", settings.clone())?;
                let cooked = tensor_type.call((&expression,), Some(&kwargs))?;
                let cooked = cooked.extract::<PyRef<'_, TensorExpression>>()?;
                let slots = cooked.interface.logical_slots();
                assert_eq!(slots.len(), 1);
                assert_eq!(slots[0].rep(), representation);
                let PartialIndex::Explicit(AbstractIndex::Symbol(symbol)) = slots[0].aind else {
                    panic!("cooking must retain an explicit external index");
                };
                let symbol: Symbol = symbol.into();
                let expected_tag = if preserve_tags {
                    "idenso::constructor_input_index"
                } else {
                    "idenso::constructor_output_index"
                };
                assert!(symbol.has_tag(expected_tag));
                assert_eq!(symbol.get_data(), &UserData::Atom(payload.clone()));
                assert_eq!(
                    settings.rust().uncook(cooked.as_super().expr.as_view()),
                    raw
                );
            }
            Ok(())
        })
        .unwrap();
    }

    #[test]
    fn wrap_dummies_preserves_raw_payloads_and_free_indices() {
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
            let wrapped = TensorExpression::wrap_dummies(tensor.borrow(py), header)?;
            let payload = FunctionBuilder::new(header).add_arg(Atom::num(17)).finish();
            let wrapped_mu = representation.to_symbolic([payload]);
            assert_eq!(
                wrapped.expr,
                ETS.metric(&wrapped_mu, &nu)
                    * FunctionBuilder::new(momentum).add_arg(&wrapped_mu).finish()
            );
            let wrapped = wrapped.into_pyobject(py)?;
            assert!(!wrapped.is_instance_of::<TensorExpression>());
            let kwargs = pyo3::types::PyDict::new(py);
            kwargs.set_item("cook_indices", PyCookSettings::indices())?;
            let cooked = py
                .get_type::<TensorExpression>()
                .call((wrapped,), Some(&kwargs))?;
            let cooked = cooked.extract::<Py<TensorExpression>>()?;
            assert_eq!(cooked.borrow(py).interface, tensor.borrow(py).interface);
            let simplified = cooked.bind(py).call_method0("simplify_metrics")?;
            let simplified = simplified.extract::<PyRef<'_, TensorExpression>>()?;
            assert_eq!(
                simplified.as_super().expr,
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
            let expression = TensorExpression::from_structured(
                py,
                SymbolicTensor::new(atom.clone(), interface),
            )?;
            assert!(
                expression
                    .bind(py)
                    .borrow()
                    .interface
                    .canonical()
                    .is_scalar()
            );

            let header = symbolica::symbol!("wrapped_index_header");
            let wrapped = TensorExpression::wrap_indices(expression.bind(py).borrow(), py, header)?;
            assert_eq!(
                wrapped.borrow(py).as_super().expr,
                atom.wrap_indices(header)
            );
            assert!(wrapped.bind(py).is_instance_of::<TensorExpression>());
            Ok(())
        })
        .unwrap();
    }

    #[test]
    fn inherited_algebra_preserves_open_indexed_zero_and_scalar_tensors() {
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
                open.borrow(py).interface.clone(),
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
                    &coefficient * &tensor.as_super().expr,
                    tensor.interface.clone(),
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
                    let result = expression.bind(py).call_method1(method, args)?;
                    let result = result.extract::<PyRef<'_, TensorExpression>>()?;
                    assert_eq!(result.interface, original.interface, "{method}");
                    assert_eq!(result.name, original.name, "{method}");
                    assert_eq!(result.name_args, original.name_args, "{method}");
                    assert!(
                        (result.as_super().expr.clone() - original.as_super().expr.clone())
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
                assert_eq!(copied.interface, original.interface);
                assert_eq!(copied.name, original.name);
                assert_eq!(copied.name_args, original.name_args);
                assert_eq!(copied.as_super().expr, original.as_super().expr);
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
                &x.expr * &tensor.borrow(py).as_super().expr,
                tensor.borrow(py).interface.clone(),
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
            assert_eq!(constructed.interface, tensor.borrow(py).interface);
            assert_eq!(constructed.name, tensor.borrow(py).name);
            let ordinary = TensorExpression::to_expression(tensor.borrow(py)).into_pyobject(py)?;
            let inferred = py.get_type::<TensorExpression>().call1((&ordinary,))?;
            assert_eq!(
                inferred.extract::<PyRef<'_, TensorExpression>>()?.interface,
                tensor.borrow(py).interface
            );

            let x = PythonExpression {
                expr: Atom::var(symbol!("tensor_convenience_x")),
            };
            let atom =
                (x.expr.clone() + Atom::one()).pow(2) * tensor.borrow(py).as_super().expr.clone();
            let expression = TensorExpression::from_known_parts(
                py,
                atom.expand(),
                tensor.borrow(py).interface.clone(),
                tensor.borrow(py).name,
                vec![Atom::num(7)],
            )?;
            let x = x.into_pyobject(py)?;
            let factored = expression.bind(py).call_method0("factor")?;
            let collected = factored.call_method1("collect", (&x,))?;
            let collected = collected.extract::<PyRef<'_, TensorExpression>>()?;
            assert_eq!(collected.interface, expression.borrow(py).interface);
            assert_eq!(collected.name, expression.borrow(py).name);
            assert_eq!(collected.name_args, vec![Atom::num(7)]);
            assert!(
                (collected.as_super().expr.clone() - expression.borrow(py).as_super().expr.clone())
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
                tensor.borrow(py).interface.clone(),
            )?;
            let zero = zero.bind(py).call_method0("factor")?;
            assert_eq!(
                zero.extract::<PyRef<'_, TensorExpression>>()?.interface,
                tensor.borrow(py).interface
            );

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
                    .as_super()
                    .expr,
                expected_indexed
                    .extract::<PyRef<'_, TensorExpression>>()?
                    .as_super()
                    .expr,
            );
            let again = TensorExpression::with_lorentz_dimension(
                promoted.borrow(py),
                py,
                ConvertibleToDimension(d),
            )?;
            assert_eq!(
                again.borrow(py).as_super().expr,
                promoted.borrow(py).as_super().expr
            );
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
            assert!(matches!(expanded_ref.as_super().expr, Atom::Add(_)));
            assert_eq!(
                expanded_ref.interface.logical_slots(),
                interface.logical_slots()
            );
            assert_eq!(expanded_ref.name, Some(name));
            assert_eq!(expanded_ref.name_args, vec![Atom::num(7)]);
            drop(expanded_ref);
            let unchanged = TensorExpression::expand(expanded.borrow(py), py, None, None)?;
            assert_eq!(
                unchanged.borrow(py).as_super().expr,
                expanded.borrow(py).as_super().expr
            );
            assert_eq!(
                unchanged.borrow(py).interface,
                expanded.borrow(py).interface
            );
            assert_eq!(unchanged.borrow(py).name, Some(name));
            assert_eq!(unchanged.borrow(py).name_args, vec![Atom::num(7)]);

            let ordinary = TensorExpression::to_expression(expression.bind(py).borrow());
            let ordinary = ordinary.expand(None, None)?;
            let ordinary = ordinary.into_pyobject(py)?;
            assert!(!ordinary.is_instance_of::<TensorExpression>());

            let zero = TensorExpression::from_atom_interface(py, Atom::Zero, interface.clone())?;
            let expanded_zero = TensorExpression::expand(zero.borrow(py), py, None, None)?;
            assert!(expanded_zero.borrow(py).as_super().expr.as_view().is_zero());
            assert_eq!(expanded_zero.borrow(py).interface, interface);
            let network = TensorExpression::to_network(zero.bind(py).borrow(), py, None)?;
            assert_eq!(
                network.structure.structure.logical_slots(),
                interface.logical_slots()
            );
            assert!(network.materialized.structure.open_positions().is_empty());
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
                for method in ["simplify", "simplify_metrics"] {
                    let result = expression.bind(py).call_method0(method)?;
                    let result = result.extract::<PyRef<'_, TensorExpression>>()?;
                    assert_eq!(result.as_super().expr, expected, "{method}");
                    assert_eq!(
                        result.interface.logical_slots(),
                        interface.logical_slots(),
                        "{method}"
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
            for method in ["simplify_gamma", "collect_gamma_chains"] {
                let result = gamma_expression.bind(py).call_method0(method)?;
                let result = result.extract::<PyRef<'_, TensorExpression>>()?;
                assert_ne!(
                    result.as_super().expr,
                    gamma_product,
                    "fixture must exercise {method}"
                );
                assert_eq!(result.interface.logical_slots(), ports, "{method}");
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
            for method in ["undo_dots", "undo_all"] {
                let result = expression.bind(py).call_method0(method)?;
                let result = result.extract::<PyRef<'_, TensorExpression>>()?;
                assert_ne!(
                    result.as_super().expr,
                    atom,
                    "fixture must exercise {method}"
                );
                assert_eq!(result.interface.logical_slots(), ports, "{method}");
            }
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
            for constructor in [false, true] {
                let result = if constructor {
                    py.get_type::<TensorExpression>()
                        .call((&original,), Some(&kwargs))?
                } else {
                    original.bind(py).call_method0("cook_indices")?
                };
                let result = result.extract::<PyRef<'_, TensorExpression>>()?;
                assert_ne!(result.as_super().expr, atom);
                assert_eq!(
                    result
                        .interface
                        .logical_slots()
                        .into_iter()
                        .map(composition::port_atom)
                        .collect::<Vec<_>>(),
                    expected
                );
            }
            // A typed zero has no index syntax to rewrite, but its interface
            // still needs the same cooking operation, including mixed open ports.
            let mut zero_slots = slots.into_iter().rev().collect::<Vec<_>>();
            zero_slots.push(rep.slot(PartialIndex::open(2)));
            let zero = TensorExpression::from_atom_interface(
                py,
                Atom::Zero,
                PartialStructure::from_logical_slots(zero_slots),
            )?;
            let cooked_zero = TensorExpression::cook_indices(zero.borrow(py), py, Some(&settings))?;
            let cooked_zero = cooked_zero.borrow(py);
            assert!(cooked_zero.as_super().expr.is_zero());
            assert_eq!(cooked_zero.interface.open_positions(), vec![2]);
            assert_eq!(
                cooked_zero.interface.logical_slots()[..2]
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
            assert_eq!(metric.borrow(py).interface.canonical().order(), 2);
            let cooked = TensorExpression::cook_indices(metric.borrow(py), py, Some(&settings))?;
            assert!(cooked.borrow(py).interface.canonical().is_scalar());
            let simplified = TensorExpression::simplify_metrics(cooked.borrow(py), py)?;
            assert_eq!(simplified.borrow(py).as_super().expr, Atom::num(3));
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
            let settings = crate::simplification::PySimplifySettings {
                metrics: false,
                expand: true,
                ..Default::default()
            };
            let simplified = TensorExpression::simplify(original.borrow(py), py, Some(&settings))?;
            assert_ne!(simplified.borrow(py).as_super().expr, atom);
            assert_eq!(
                simplified.borrow(py).as_super().expr,
                expanded.borrow(py).as_super().expr
            );
            assert_eq!(
                simplified.borrow(py).interface.logical_slots(),
                interface.logical_slots()
            );
            Ok(())
        })
        .unwrap();
    }

    #[test]
    fn unsafe_constructor_copies_even_inconsistent_parts_without_validation() {
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let rep = ExtendibleReps::EUCLIDEAN.new_rep(Dimension::Concrete(3));
            let slot = rep.slot(PartialIndex::Explicit(AbstractIndex::Normal(101)));
            // Duplicate indices would normally contract. The unchecked API must
            // copy this metadata rather than normalize or "repair" it.
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
                let result = py.get_type::<TensorExpression>().call_method(
                    "unsafe_from_expression",
                    (expression,),
                    Some(&kwargs),
                )?;
                let result = result.extract::<PyRef<'_, TensorExpression>>()?;
                assert_eq!(result.as_super().expr, atom);
                assert_eq!(
                    result.interface.logical_slots(),
                    structure.interface.logical_slots()
                );
                assert_eq!(result.name, structure.name);
                assert_eq!(result.name_args, structure.arguments);
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
                    let result = original.bind(py).call_method1(method, args)?;
                    let result = result.extract::<PyRef<'_, TensorExpression>>()?;
                    assert_eq!(
                        result.interface.logical_slots(),
                        interface.logical_slots(),
                        "{method}"
                    );
                    assert!(
                        (&result.as_super().expr - &atom)
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
                        expanded.borrow(py).interface.logical_slots(),
                        interface.logical_slots()
                    );
                }
                let settings = crate::simplification::PySimplifySettings {
                    metrics: false,
                    expand: true,
                    ..Default::default()
                };
                let simplified =
                    TensorExpression::simplify(original.borrow(py), py, Some(&settings))?;
                assert_eq!(
                    simplified.borrow(py).interface.logical_slots(),
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
                    tensor.borrow(py).as_super().expr.clone(),
                    PartialStructure::from_logical_slots([mismatched])
                )
                .is_err()
            );
            let number = PythonExpression { expr: Atom::num(2) }.into_pyobject(py)?;
            let doubled = tensor.bind(py).call_method1("__mul__", (number,))?;
            let doubled = doubled.extract::<PyRef<'_, TensorExpression>>()?;
            assert_eq!(doubled.interface.logical_slots(), interface.logical_slots());
            let original = tensor.borrow(py);
            let invalid_rewrite =
                TensorExpression::preserving_interface(&original, py, Atom::one());
            assert!(invalid_rewrite.is_err());
            let unchanged = TensorExpression::preserving_interface(
                &original,
                py,
                original.as_super().expr.clone(),
            )?;
            assert_eq!(
                unchanged.borrow(py).interface.logical_slots(),
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
                    assert_eq!(expanded.as_super().expr, atom.expand());
                    assert_eq!(expanded.interface, interface);
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
            assert!(canceled.borrow(py).as_super().expr.is_zero());
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
                        assert_eq!(actual.as_super().expr, expected.as_super().expr);
                        assert_eq!(actual.interface, expected.interface);
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
            assert_eq!(actual.as_super().expr, expected.as_super().expr);
            assert_eq!(actual.interface, expected.interface);
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
        let expected_sum = left.expression.as_ref() + right.expression.as_ref();
        let expected_difference = left.expression.as_ref() - right.expression.as_ref();

        assert_ne!(
            left.structure.logical_slots(),
            right.structure.logical_slots()
        );
        assert!(InterfaceInference::additive_interfaces_match(
            &left.structure,
            &right.structure
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
            &left.structure,
            &inferred
        ));

        let standalone = SymbolicTensor::new(left.expression.clone(), right.structure.clone());
        assert_eq!(
            standalone.presentation_atom(),
            tensor("canonical_addition_left", [j, i]).expression
        );

        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let expected = left.structure.logical_slots();
            let left = TensorExpression::from_structured(py, left)?;
            let right = TensorExpression::from_structured(py, right)?;
            let TensorDispatch::Expression(sum) =
                TensorExpression::__add__(left.bind(py).borrow(), py, right.bind(py).as_any())?
            else {
                panic!("adding symbolic tensor expressions produced a network")
            };
            let sum = sum.bind(py).borrow();
            assert_eq!(sum.interface.logical_slots(), expected);
            assert_eq!(sum.as_super().expr, expected_sum);
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
            assert_eq!(difference.interface.logical_slots(), expected);
            assert_eq!(difference.as_super().expr, expected_difference);
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
            assert_eq!(reverse_difference.interface.logical_slots(), expected);
            assert_eq!(reverse_difference.as_super().expr, expected_difference);
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
                assert_eq!(quotient.as_super().expr, &expression / &divisor);
                assert_eq!(quotient.interface.logical_slots(), ports);
                assert_eq!(numerator.borrow(py).interface.logical_slots(), ports);
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
                assert_eq!(value.as_super().expr, atom);
                assert_eq!(value.interface.logical_slots(), interface.logical_slots());
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
            assert_eq!(sum.as_super().expr, atom.as_ref() + atom.as_ref());
            assert_eq!(sum.interface.logical_slots(), interface.logical_slots());
            assert_eq!(sum.name, None);
            assert!(sum.name_args.is_empty());
            drop(sum);

            let negated = expression.bind(py).call_method1("__rsub__", (0,))?;
            let negated = negated.extract::<PyRef<'_, TensorExpression>>()?;
            assert_eq!(negated.as_super().expr, Atom::num(-1) * atom.as_ref());
            assert_eq!(negated.interface.logical_slots(), interface.logical_slots());
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
        assert_eq!(inferred.structure.canonical().order(), 1);
        assert_eq!(inferred.expression, product);
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
            assert!(tensor.borrow(py).interface.canonical().is_scalar());
            assert_eq!(tensor.borrow(py).as_super().expr, atom);
            let simplified = tensor.bind(py).call_method0("simplify_color")?;
            let simplified = simplified.extract::<PyRef<'_, TensorExpression>>()?;
            assert!(simplified.interface.canonical().is_scalar());
            assert_eq!(
                simplified.as_super().expr,
                Atom::num(-8) * CS.cas(Atom::num(2), adjoint.to_symbolic([]))
            );

            let open = TensorExpression::f(py, ConvertibleToDimension(8.into()))?;
            let square = open.borrow(py).as_super().expr.as_ref().pow(Atom::num(2));
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
            let reinferred = TensorExpression::reinfer(composed.bind(py).borrow(), py)?;
            assert_eq!(
                reinferred.bind(py).borrow().interface.canonical().order(),
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

            assert!(indexed.interface.canonical().is_scalar());
            assert!(
                matches!(indexed.as_super().expr.as_view(), AtomView::Fun(function) if function.get_symbol() == SPENSO_TAG.trace)
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

            assert!(expression.interface.canonical().is_scalar());
            assert!(
                matches!(expression.as_super().expr.as_view(), AtomView::Fun(function) if function.get_symbol() == SPENSO_TAG.trace)
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
                composed.expression.as_view(),
                AtomView::Fun(function) if function.get_symbol() == SPENSO_TAG.chain
            ));
            let rust_list_error = composed
                .expression
                .list_dangling::<AbstractIndex>()
                .expect_err("unresolved chain endpoints should return a parse error");
            let rust_wrap_error = composed
                .expression
                .wrap_dummies::<AbstractIndex>(symbolica::symbol!("index_tooling_wrap"))
                .expect_err("unresolved chain endpoints should return a parse error");
            assert!(matches!(
                rust_list_error,
                IndexToolingError::ListDangling { reason }
                    if reason.contains("invalid chain start")
            ));
            assert!(matches!(
                rust_wrap_error,
                IndexToolingError::WrapDummies { reason }
                    if reason.contains("invalid chain start")
            ));

            let expression = TensorExpression::from_structured(py, composed)?;
            let list_error = TensorExpression::list_dangling(expression.bind(py).borrow())
                .err()
                .expect("unresolved chain endpoints should return a parse error");
            let wrap_error = TensorExpression::wrap_dummies(
                expression.bind(py).borrow(),
                symbolica::symbol!("index_tooling_wrap"),
            )
            .err()
            .expect("unresolved chain endpoints should return a parse error");

            for (error, context) in [
                (list_error, "cannot list dangling indices"),
                (wrap_error, "cannot wrap dummy indices"),
            ] {
                assert!(error.is_instance_of::<PyValueError>(py));
                assert!(error.to_string().contains(context));
                assert!(error.to_string().contains("invalid chain start"));
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

            for name in ["undo_dots", "schoonschip_net", "dirac_adjoint"] {
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
            assert_eq!(original.borrow(py).interface.canonical().order(), 2);
            let result = TensorExpression::with_lorentz_dimension(
                original.borrow(py),
                py,
                ConvertibleToDimension(Dimension::Concrete(6)),
            )?;
            assert!(result.borrow(py).interface.canonical().is_scalar());
            assert!(result.borrow(py).name.is_none());
            assert!(result.borrow(py).name_args.is_empty());

            for dimension in [4, 5] {
                let result = TensorExpression::with_lorentz_dimension(
                    original.borrow(py),
                    py,
                    ConvertibleToDimension(Dimension::Concrete(dimension)),
                )?;
                assert_eq!(result.borrow(py).interface.canonical().order(), 2);
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
                    assert_eq!(tensor.as_super().expr, atom);
                    assert_eq!(tensor.interface.canonical().order(), 1);
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
            assert_eq!(tensor.as_super().expr, atom);
            assert!(tensor.interface.canonical().is_scalar());
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
            assert_eq!(result.interface.canonical().order(), 0);
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
