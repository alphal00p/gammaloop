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
    TensorStructure,
    dimension::Dimension,
    partial::{PartialIndex, PartialStructure, PartialStructureExt},
    representation::{LibraryRep, RepName},
    slot::IsAbstractSlot,
};
use symbolica::{
    api::python::{ConvertibleToExpression, PythonExpression},
    atom::{Atom, Symbol},
};

use crate::{
    ModuleInit, display,
    structure::{SpensoName, SpensoRepresentation, SpensoSlot},
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
#[spenso_macros::track_usage(crate::record_usage)]
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

/// The identity, scalar arguments, and ordered axes of a tensor.
///
/// A Slot denotes a labeled axis; a Representation denotes an unresolved axis.
/// This immutable object stores metadata, without an expression or components.
/// Axis positions and component shapes follow the order of ``slots``.
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
    pub(crate) interface: PartialStructure,
    pub(crate) name: Option<Symbol>,
    pub(crate) arguments: Vec<Atom>,
}

impl ModuleInit for SpensoTensorStructure {}

pub(crate) fn validate_axis_permutation(rank: usize, axes: &[usize]) -> PyResult<()> {
    let mut sorted = axes.to_vec();
    sorted.sort_unstable();
    if sorted != (0..rank).collect::<Vec<_>>() {
        return Err(PyValueError::new_err(format!(
            "axes must be a permutation of 0..{rank}"
        )));
    }
    Ok(())
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "python_stubgen"), remove_gen_stub)]
#[spenso_macros::track_usage(crate::record_usage)]
#[pymethods]
impl SpensoTensorStructure {
    /// Describe a tensor interface in logical axis order.
    ///
    /// Parameters
    /// ----------
    /// slots : sequence of Slot or Representation
    ///     External axes, in the order used for indexing and component coordinates.
    /// name : TensorName, optional
    ///     Name identifying stored tensor data.
    /// arguments : sequence of scalar expressions, optional
    ///     Non-index arguments that distinguish named tensors, such as A(x, 7).
    ///
    /// Returns
    /// -------
    /// TensorStructure
    ///     Immutable metadata for the supplied axes and optional identity.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import Representation, TensorStructure
    /// >>> metadata = TensorStructure([Representation.euc(2), Representation.euc(3)])
    /// >>> metadata.shape
    /// (2, 3)
    #[new]
    #[pyo3(signature = (slots, *, name=None, arguments=None))]
    fn new(
        #[gen_stub(override_type(type_repr = "typing.Sequence[Slot | Representation]", imports = ("typing")))]
        slots: Vec<Bound<'_, PyAny>>,
        name: Option<SpensoName>,
        arguments: Option<Vec<ConvertibleToExpression>>,
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
        Ok(Self {
            interface: PartialStructure::from_logical_slots(slots),
            name: name.map(|value| value.name),
            arguments: arguments
                .unwrap_or_default()
                .into_iter()
                .map(|arg| arg.to_expression().expr)
                .collect(),
        })
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
    /// >>> A = A.structure
    /// >>> name = A.name
    #[getter]
    fn name(&self) -> Option<SpensoName> {
        self.name.map(|name| SpensoName { name })
    }

    /// Scalar arguments that distinguish the named tensor.
    ///
    /// Returns
    /// -------
    /// tuple of Expression
    ///     Non-index arguments, in their original order; empty if there are none.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import Representation, TensorName, TensorExpression
    /// >>> space = Representation.euc(3)
    /// >>> A = TensorName("A")(space, space)
    /// >>> A.structure.arguments
    /// ()
    #[getter]
    #[gen_stub(override_return_type(type_repr = "tuple[Expression, ...]"))]
    fn arguments(&self, py: Python<'_>) -> PyResult<Py<PyTuple>> {
        Ok(PyTuple::new(
            py,
            self.arguments.iter().cloned().map(PythonExpression::from),
        )?
        .unbind())
    }

    /// External axes, in the order used by tensor indexing.
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
    /// >>> A.structure.slots == (space, space)
    /// True
    #[getter]
    #[gen_stub(override_return_type(type_repr = "tuple[Slot | Representation, ...]"))]
    fn slots(&self, py: Python<'_>) -> PyResult<Py<PyTuple>> {
        let slots = self
            .interface
            .logical_slots()
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
        self.interface.canonical().order()
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
    /// >>> A = A.structure
    /// >>> A.shape
    /// (3, 3)
    #[getter]
    #[gen_stub(override_return_type(type_repr = "tuple[int | Expression, ...]"))]
    pub(crate) fn shape(&self, py: Python<'_>) -> PyResult<Py<PyTuple>> {
        let shape = self
            .interface
            .logical_slots()
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

    /// Reorder the metadata axes.
    ///
    /// Parameters
    /// ----------
    /// axes : sequence of int
    ///     Each old axis position exactly once, in the desired new order.
    ///
    /// Returns
    /// -------
    /// TensorStructure
    ///     New metadata with the same name and scalar arguments.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import Representation, TensorStructure
    /// >>> metadata = TensorStructure([Representation.euc(2), Representation.euc(3)])
    /// >>> metadata.permute_axes([1, 0]).shape
    /// (3, 2)
    fn permute_axes(&self, axes: Vec<usize>) -> PyResult<Self> {
        validate_axis_permutation(self.rank(), &axes)?;
        let slots = self.interface.logical_slots();
        Ok(Self {
            interface: PartialStructure::from_logical_slots(axes.iter().map(|&axis| slots[axis])),
            name: self.name,
            arguments: self.arguments.clone(),
        })
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
            .logical_slots()
            .iter()
            .map(|slot| idenso::tensor::composition::port_atom(*slot).to_string())
            .collect::<Vec<_>>()
            .join(", ");
        let name = self
            .name
            .map(|s| format!(", name={:?}", s.get_name()))
            .unwrap_or_default();
        let args = self
            .arguments
            .iter()
            .map(ToString::to_string)
            .collect::<Vec<_>>()
            .join(", ");
        format!("TensorStructure([{slots}]{name}, arguments=({args}))")
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
