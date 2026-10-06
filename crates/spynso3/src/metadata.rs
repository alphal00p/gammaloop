//! Python views of tensor metadata, independent of expressions and component data.
use pyo3::{
    exceptions::{PyTypeError, PyValueError},
    prelude::*,
    types::PyTuple,
};
#[cfg(feature = "python_stubgen")]
use pyo3_stub_gen::derive::*;
#[cfg(not(feature = "python_stubgen"))]
use pyo3_stub_gen_derive::remove_gen_stub;
use spenso::structure::{
    OrderedStructure, TensorStructure,
    abstract_index::AbstractIndex,
    dimension::Dimension,
    partial::{PartialIndex, PartialSlot, PartialStructure, PartialStructureExt},
    representation::{LibraryRep, RepName},
    slot::IsAbstractSlot,
};
use symbolica::api::python::PythonExpression;

use crate::{
    ModuleInit, display,
    structure::{SpensoRepresentation, SpensoSlot},
};

/// A representation identity independent of its dimension.
///
/// Obtain this read-only metadata through ``Representation.name``. It distinguishes
/// self-dual spaces from the two partners of a dualizable space.
///
/// Examples
/// --------
/// >>> from symbolica.community.tensor import Representation
/// >>> identity = Representation.mink(4).name
/// >>> identity.duality
/// 'self-dual'
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    frozen,
    from_py_object,
    eq,
    name = "RepresentationName",
    module = "symbolica.community.tensor"
)]
#[derive(Clone, PartialEq, Eq)]
pub struct SpensoRepresentationName {
    pub(crate) rep: LibraryRep,
}

impl ModuleInit for SpensoRepresentationName {}

impl SpensoRepresentationName {
    pub(crate) fn metric_rule(&self) -> &'static str {
        if !self.rep.is_self_dual() {
            "dual pairing"
        } else if self.rep == spenso::structure::representation::ExtendibleReps::MINKOWSKI {
            "(+, −, …)"
        } else if matches!(self.rep, LibraryRep::InlineMetric(_)) {
            "index-dependent"
        } else {
            "(+, +, …)"
        }
    }
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[spenso_macros::track_usage(crate::record_usage, on_success)]
#[pymethods]
impl SpensoRepresentationName {
    /// The exact registered name of the representation.
    ///
    /// Returns
    /// -------
    /// str
    ///     Name without dimension or index arguments.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import Representation
    /// >>> identity = Representation.mink(4).name
    /// >>> name = identity.name
    #[getter]
    pub fn name(&self) -> String {
        self.rep.symbol().get_name().to_owned()
    }

    /// How this representation pairs with another space.
    ///
    /// Returns
    /// -------
    /// str
    ///     One of "self-dual", "base", or "dual".
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import Representation
    /// >>> identity = Representation.mink(4).name
    /// >>> identity.duality
    /// 'self-dual'
    #[getter]
    pub fn duality(&self) -> &'static str {
        if self.rep.is_self_dual() {
            "self-dual"
        } else if self.rep.is_dual() {
            "dual"
        } else {
            "base"
        }
    }

    /// Return the dimension-independent identity of the paired space.
    ///
    /// Returns
    /// -------
    /// RepresentationName
    ///     The same identity for a self-dual space; the dual partner otherwise.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import Representation
    /// >>> identity = Representation.cof(3).name
    /// >>> identity.dual().dual() == identity
    /// True
    fn dual(&self) -> Self {
        Self {
            rep: self.rep.dual(),
        }
    }

    /// Return the sign of the canonical pairing at a component coordinate.
    ///
    /// Parameters
    /// ----------
    /// index : int
    ///     Nonnegative component coordinate. This identity has no dimension, so
    ///     the coordinate is not checked against a particular tensor shape.
    ///
    /// Returns
    /// -------
    /// int
    ///     +1 or -1. In Minkowski space, zero gives +1 and spatial coordinates -1.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import Representation
    /// >>> metric = Representation.mink(4).name
    /// >>> [metric.metric_sign(i) for i in range(4)]
    /// [1, -1, -1, -1]
    fn metric_sign(&self, index: usize) -> i8 {
        if self.rep.is_neg(index) { -1 } else { 1 }
    }

    /// Return a readable text representation.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> name = sp.Representation.mink(4).name
    /// >>> text = str(name)
    fn __str__(&self) -> String {
        self.rep.to_string()
    }

    /// Return a readable object description for inspection.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> name = sp.Representation.mink(4).name
    /// >>> text = repr(name)
    fn __repr__(&self) -> String {
        format!(
            "RepresentationName({}, duality={:?})",
            self.name(),
            self.duality()
        )
    }

    /// Render a compact HTML view of this tensor metadata.
    ///
    /// Returns
    /// -------
    /// str
    ///     HTML fragment for display in a notebook or page.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import Representation
    /// >>> value = Representation.mink(4).name
    /// >>> html = value.to_html()
    fn to_html(&self) -> String {
        display::metadata::representation_name(self)
    }

    fn _repr_html_(&self) -> String {
        self.to_html()
    }
}

/// The canonical external signature of an opaque tensor.
///
/// A Slot denotes a labeled axis; a Representation denotes an unresolved axis.
/// This immutable signature contains only free axes: no tensor name, scalar
/// arguments, expression, component data, or layout permutation. Axes follow
/// Spenso's canonical representation/index ordering. A(i,j), A(j,i), and their
/// sum have the same structure. The expression retains argument order; tensor
/// and network component views retain their own axis order.
///
/// Examples
/// --------
/// >>> from symbolica.community.tensor import Representation, TensorName, TensorExpression
/// >>> space = Representation.euc(3)
/// >>> A = TensorName("A")(space, space)
/// >>> metadata = A.structure
/// >>> metadata.rank
/// 2
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    frozen,
    from_py_object,
    eq,
    name = "TensorStructure",
    module = "symbolica.community.tensor"
)]
#[derive(Clone, PartialEq, Eq)]
pub struct SpensoTensorStructure {
    pub(crate) interface: OrderedStructure<LibraryRep, PartialIndex<AbstractIndex>>,
}

impl ModuleInit for SpensoTensorStructure {}

impl SpensoTensorStructure {
    pub(crate) fn from_interface(interface: &PartialStructure) -> Self {
        // Open identifiers are occurrence-local. Renumber them in canonical order
        // so two equal signatures do not retain the source argument positions.
        let mut interface = interface.canonical().clone();
        let mut next_open = 0;
        for slot in &mut interface {
            if matches!(slot.aind, PartialIndex::Open(_)) {
                slot.aind = PartialIndex::open(next_open);
                next_open += 1;
            }
        }
        Self { interface }
    }
}

pub(crate) fn axis_objects(py: Python<'_>, slots: Vec<PartialSlot>) -> PyResult<Py<PyTuple>> {
    let slots = slots
        .into_iter()
        .map(|slot| match slot.aind {
            PartialIndex::Explicit(index) => Py::new(
                py,
                SpensoSlot {
                    slot: slot.rep().slot(index),
                },
            )
            .map(Py::into_any),
            PartialIndex::Open(_) => Py::new(
                py,
                SpensoRepresentation {
                    representation: slot.rep(),
                },
            )
            .map(Py::into_any),
        })
        .collect::<PyResult<Vec<_>>>()?;
    Ok(PyTuple::new(py, slots)?.unbind())
}

pub(crate) fn axis_shape(py: Python<'_>, slots: Vec<PartialSlot>) -> PyResult<Py<PyTuple>> {
    let shape = slots
        .into_iter()
        .map(|slot| match slot.rep().dim {
            Dimension::Concrete(n) => n
                .into_pyobject(py)
                .map(|v| v.into_any().unbind())
                .map_err(Into::into),
            dim => Py::new(py, PythonExpression::from(dim.to_symbolic())).map(Py::into_any),
        })
        .collect::<PyResult<Vec<_>>>()?;
    Ok(PyTuple::new(py, shape)?.unbind())
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "python_stubgen"), remove_gen_stub)]
#[spenso_macros::track_usage(crate::record_usage, on_success)]
#[pymethods]
impl SpensoTensorStructure {
    /// Describe an opaque tensor by its canonical free axes.
    ///
    /// Parameters
    /// ----------
    /// slots : sequence of Slot or Representation
    ///     Free axes, sorted into canonical order. Repeated unresolved spaces
    ///     remain distinct axes. This describes a signature without contracting it.
    ///
    /// Returns
    /// -------
    /// TensorStructure
    ///     Canonical signature, independent of the supplied sequence order.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import Representation, TensorStructure
    /// >>> metadata = TensorStructure([Representation.euc(2), Representation.euc(3)])
    /// >>> metadata.shape
    /// (2, 3)
    #[new]
    #[pyo3(signature = (slots))]
    fn new(
        #[gen_stub(override_type(type_repr = "typing.Sequence[Slot | Representation]", imports = ("typing")))]
        slots: Vec<Bound<'_, PyAny>>,
    ) -> PyResult<Self> {
        let slots = slots
            .into_iter()
            .enumerate()
            .map(|(position, slot)| {
                if let Ok(slot) = slot.extract::<SpensoSlot>() {
                    Ok(slot.slot.rep().slot(PartialIndex::Explicit(slot.slot.aind)))
                } else if let Ok(rep) = slot.extract::<SpensoRepresentation>() {
                    Ok(rep.representation.slot(PartialIndex::open(position)))
                } else {
                    Err(PyTypeError::new_err(
                        "structure slots must be Slot or Representation objects",
                    ))
                }
            })
            .collect::<PyResult<Vec<_>>>()?;
        Ok(Self::from_interface(&PartialStructure::from_logical_slots(
            slots,
        )))
    }

    /// Free axes in canonical representation/index order.
    ///
    /// The signature is independent of argument order, factor order, and summand
    /// order. It describes the whole expression as an opaque tensor. Component
    /// coordinates belong to the tensor's current view, not to this signature.
    ///
    /// Returns
    /// -------
    /// tuple of Slot or Representation
    ///     Labeled and unresolved axes, respectively.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import Representation, TensorName, TensorExpression
    /// >>> space = Representation.euc(3)
    /// >>> A = TensorName("A")(space, space)
    /// >>> A.structure.axes == (space, space)
    /// True
    #[getter]
    #[gen_stub(override_return_type(type_repr = "tuple[Slot | Representation, ...]"))]
    fn axes(&self, py: Python<'_>) -> PyResult<Py<PyTuple>> {
        axis_objects(py, self.interface.external_structure())
    }

    /// Return every external axis as an explicitly indexed slot.
    ///
    /// Returns
    /// -------
    /// list of Slot
    ///     Slots in canonical order; an empty list for a scalar. Representations
    ///     and dimensions may differ between axes.
    ///
    /// Raises
    /// ------
    /// ValueError
    ///     If any axis is unresolved. No indices are invented or axes omitted.
    ///     Use ``axes`` to inspect a mixed structure.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import Representation, TensorStructure
    /// >>> space = Representation.euc(3)
    /// >>> metadata = TensorStructure([space("j"), space("i")])
    /// >>> metadata.slots() == TensorStructure([space("i"), space("j")]).slots()
    /// True
    fn slots(&self) -> PyResult<Vec<SpensoSlot>> {
        self.interface
            .slots()
            .map(|slots| slots.into_iter().map(|slot| SpensoSlot { slot }).collect())
            .map_err(|error| PyValueError::new_err(error.to_string()))
    }

    /// Return every external axis as an unresolved representation.
    ///
    /// Returns
    /// -------
    /// list of Representation
    ///     Representations in canonical order; an empty list for a scalar.
    ///     Different spaces, dimensions, and dualities are allowed.
    ///
    /// Raises
    /// ------
    /// ValueError
    ///     If any axis has an explicit index. Labels are never discarded.
    ///     Use ``axes`` to inspect a mixed structure.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import Representation, TensorStructure
    /// >>> spaces = [Representation.mink(4), Representation.euc(3)]
    /// >>> TensorStructure(spaces) == TensorStructure(spaces[::-1])
    /// True
    fn representations(&self) -> PyResult<Vec<SpensoRepresentation>> {
        self.interface
            .representations()
            .map(|representations| {
                representations
                    .into_iter()
                    .map(|representation| SpensoRepresentation { representation })
                    .collect()
            })
            .map_err(|error| PyValueError::new_err(error.to_string()))
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
    /// >>> A = A.structure
    /// >>> A.rank
    /// 2
    #[getter]
    pub(crate) fn rank(&self) -> usize {
        self.interface.order()
    }

    /// Dimensions in canonical signature order.
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
    /// >>> A = A.structure
    /// >>> A.shape
    /// (3, 3)
    #[getter]
    #[gen_stub(override_return_type(type_repr = "tuple[int | Expression, ...]"))]
    pub(crate) fn shape(&self, py: Python<'_>) -> PyResult<Py<PyTuple>> {
        axis_shape(py, self.interface.external_structure())
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
    /// >>> A = A.structure
    /// >>> A.is_scalar
    /// False
    #[getter]
    fn is_scalar(&self) -> bool {
        self.rank() == 0
    }

    /// Return a readable object description for inspection.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> r = sp.Representation.euc(2)
    /// >>> structure = sp.TensorStructure([r("i"), r("j")])
    /// >>> text = repr(structure)
    fn __repr__(&self) -> String {
        let slots = self
            .interface
            .external_structure()
            .iter()
            .map(|slot| idenso::tensor::composition::port_atom(*slot).to_string())
            .collect::<Vec<_>>()
            .join(", ");
        format!("TensorStructure([{slots}])")
    }

    /// Render a compact HTML view of this tensor metadata.
    ///
    /// Parameters
    /// ----------
    /// settings : DisplaySettings, optional
    ///     Index and tensor presentation choices; defaults to DisplaySettings().
    ///
    /// Returns
    /// -------
    /// str
    ///     HTML fragment for display in a notebook or page.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import Representation, TensorName, TensorExpression
    /// >>> space = Representation.euc(2)
    /// >>> A = TensorName("M")(space, space)
    /// >>> value = A.structure
    /// >>> html = value.to_html()
    #[pyo3(signature = (*, settings=None))]
    fn to_html(
        &self,
        py: Python<'_>,
        settings: Option<display::DisplaySettings>,
    ) -> PyResult<String> {
        display::metadata::structure(py, self, &settings.unwrap_or_default())
    }

    fn _repr_html_(&self, py: Python<'_>) -> PyResult<String> {
        self.to_html(py, None)
    }
}
