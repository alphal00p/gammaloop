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

/// A dimension-independent representation identity, including its duality.
///
/// Examples
/// --------
/// >>> from symbolica.community import tensor as sp
/// >>> name = sp.Representation.mink(4).name
/// >>> text = repr(name)
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
    /// Exact registered symbol name, without dimension or index arguments.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> name = sp.Representation.mink(4).name
    /// >>> result = name.name
    #[getter]
    pub fn name(&self) -> String {
        self.rep.symbol().get_name().to_owned()
    }

    /// Whether this is a self-dual, base, or dual representation.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> name = sp.Representation.mink(4).name
    /// >>> result = name.duality
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

    /// Return the representation identity paired with this one.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> fundamental = sp.Representation.cof(3).name
    /// >>> antifundamental = fundamental.dual()
    fn dual(&self) -> Self {
        Self {
            rep: self.rep.dual(),
        }
    }

    /// Sign of the canonical contraction pairing at a nonnegative component index.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> name = sp.Representation.mink(4).name
    /// >>> time_sign = name.metric_sign(0)
    /// >>> space_sign = name.metric_sign(1)
    /// >>> assert (time_sign, space_sign) == (1, -1)
    ///
    /// Parameters
    /// ----------
    /// index : int
    ///     Nonnegative component position. Minkowski position 0 has positive metric sign; spatial positions have negative sign.
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

    /// Compact metadata display with the registered name, duality, and metric rule.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> name = sp.Representation.mink(4).name
    /// >>> html = name.to_html()
    fn to_html(&self) -> String {
        display::metadata::representation_name(self)
    }

    fn _repr_html_(&self) -> String {
        self.to_html()
    }
}

/// Describe a tensor’s identity and ordered external axes without storing component data.
/// Use Slots for indexed axes and Representations for ports that will be indexed later.
///
/// Examples
/// --------
/// >>> from symbolica.community import tensor as sp
/// >>> r = sp.Representation.euc(2)
/// >>> structure = sp.TensorStructure([r("i"), r], name=sp.TensorName("docs::A"))
/// >>> assert structure.rank == 2
///
/// Parameters
/// ----------
/// slots : Sequence[Slot | Representation]
///     Ordered external ports. A Slot fixes an index label; a Representation leaves the port unresolved.
/// name : TensorName, optional
///     Optional TensorName recording the data identity. Omit it for anonymous tensor metadata.
/// arguments : Sequence[Expression or scalar], optional
///     Optional scalar arguments belonging to the tensor identity. These are metadata, not extra tensor axes.
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
    /// Describe a tensor’s identity and ordered external axes without storing component data.
    /// Use Slots for indexed axes and Representations for ports that will be indexed later.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> r = sp.Representation.euc(2)
    /// >>> structure = sp.TensorStructure([r("i"), r], name=sp.TensorName("docs::A"))
    /// >>> assert structure.rank == 2
    ///
    /// Parameters
    /// ----------
    /// slots : Sequence[Slot | Representation]
    ///     Ordered external ports. A Slot fixes an index label; a Representation leaves the port unresolved.
    /// name : TensorName, optional
    ///     Optional TensorName recording the data identity. Omit it for anonymous tensor metadata.
    /// arguments : Sequence[Expression or scalar], optional
    ///     Optional scalar arguments belonging to the tensor identity. These are metadata, not extra tensor axes.
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

    /// Optional tensor data identity; unnamed expressions have no name.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> r = sp.Representation.euc(2)
    /// >>> structure = sp.TensorStructure([r("i"), r("j")])
    /// >>> result = structure.name
    #[getter]
    fn name(&self) -> Option<SpensoName> {
        self.name.map(|name| SpensoName { name })
    }

    /// Scalar arguments belonging to the tensor identity, in their original order.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> r = sp.Representation.euc(2)
    /// >>> structure = sp.TensorStructure([r("i"), r("j")])
    /// >>> result = structure.arguments
    #[getter]
    #[gen_stub(override_return_type(type_repr = "tuple[Expression, ...]"))]
    fn arguments(&self, py: Python<'_>) -> PyResult<Py<PyTuple>> {
        Ok(PyTuple::new(
            py,
            self.arguments.iter().cloned().map(PythonExpression::from),
        )?
        .unbind())
    }

    /// External slots in logical order, including unresolved representations.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> r = sp.Representation.euc(2)
    /// >>> structure = sp.TensorStructure([r("i"), r("j")])
    /// >>> result = structure.slots
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

    /// Number of external slots, including unresolved ports.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> r = sp.Representation.euc(2)
    /// >>> structure = sp.TensorStructure([r("i"), r("j")])
    /// >>> result = structure.rank
    #[getter]
    pub(crate) fn rank(&self) -> usize {
        self.interface.canonical().order()
    }

    /// Dimensions in logical order; symbolic dimensions remain symbolic.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> r = sp.Representation.euc(2)
    /// >>> structure = sp.TensorStructure([r("i"), r("j")])
    /// >>> result = structure.shape
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

    /// Whether there are no external axes.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> r = sp.Representation.euc(2)
    /// >>> structure = sp.TensorStructure([r("i"), r("j")])
    /// >>> result = structure.is_scalar
    #[getter]
    fn is_scalar(&self) -> bool {
        self.rank() == 0
    }

    /// Reorder metadata axes without changing their labels, dimensions, or identity.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> r = sp.Representation.euc(2)
    /// >>> structure = sp.TensorStructure([r("i"), r("j")])
    /// >>> transposed = structure.permute_axes([1, 0])
    /// >>> assert transposed.slots == structure.slots[::-1]
    ///
    /// Parameters
    /// ----------
    /// axes : Sequence[int]
    ///     Permutation of all axis positions, starting at zero. For a matrix, [1, 0] transposes the logical axes.
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

    /// Inspect the tensor identity and ordered slots without changing its algebra.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> r = sp.Representation.euc(2)
    /// >>> structure = sp.TensorStructure([r("i"), r("j")])
    /// >>> html = structure.to_html()
    ///
    /// Parameters
    /// ----------
    /// settings : DisplaySettings, optional
    ///     Optional DisplaySettings controlling notation, indices, and component layout. None uses the default presentation.
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
