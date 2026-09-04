use std::collections::{HashMap, HashSet};

use idenso::{
    Cookable,
    color::CS,
    dirac::AGS,
    python::{
        PyColorCasimirSettings, PyColorSimplifySettings, PyCookSettings, PyGammaSimplifySettings,
        PySchoonschipSettings, RegisteredRepresentation,
    },
    representations::{Bispinor, ColorAdjoint, ColorAntiFundamental, ColorFundamental},
    tensor::SymbolicTensor,
};
use pyo3::{
    exceptions::{PyIndexError, PyOverflowError, PyRuntimeError, PyTypeError, PyValueError},
    prelude::*,
    types::{PyAny, PyTuple},
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
        parsing::{AtomStructureExt, StrictTensorFilter, StructureInferenceMode},
        tags::SPENSO_TAG,
    },
    shadowing,
    structure::{
        Canonicalized, HasName, OrderedStructure, StructureError, TensorStructure,
        abstract_index::AbstractIndex,
        dimension::Dimension,
        partial::{PartialIndex, PartialStructure, PartialStructureExt},
        representation::{LibraryRep, Minkowski, RepName, Representation},
        slot::{IsAbstractSlot, Slot},
    },
};
use symbolica::{
    api::python::{ConvertibleToExpression, PythonExpression, PythonFormattedOutput},
    atom::{Atom, AtomCore, AtomView, FunctionBuilder, Symbol},
    domains::rational::Rational,
};

use crate::{
    ModuleInit, SliceOrIntOrExpanded, Spensor,
    composition::{self, StructuredAtom},
    display,
    library::SpensorLibrary,
    network::{ConvertibleToSpensoNet, SpensoNet},
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
/// >>> gamma = TensorExpression.gamma(4)("mu", "i", "j")
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
    pub(crate) fn from_atom_interface(
        py: Python<'_>,
        atom: Atom,
        interface: PartialStructure,
    ) -> PyResult<Py<Self>> {
        Self::from_atom_interface_named(py, atom, interface, None)
    }

    pub(crate) fn from_atom_interface_named(
        py: Python<'_>,
        atom: Atom,
        interface: PartialStructure,
        name: Option<Symbol>,
    ) -> PyResult<Py<Self>> {
        Self::from_atom_interface_descriptor(py, atom, interface, name, Vec::new())
    }

    pub(crate) fn from_atom_interface_descriptor(
        py: Python<'_>,
        atom: Atom,
        interface: PartialStructure,
        name: Option<Symbol>,
        name_args: Vec<Atom>,
    ) -> PyResult<Py<Self>> {
        let original_rank = interface.canonical().order();
        atom.validate_chain_like_nesting()
            .map_err(|error| PyValueError::new_err(error.to_string()))?;
        validate_placeholder_scope(atom.as_view(), false).map_err(PyValueError::new_err)?;
        composition::validate_explicit_index_occurrences(&atom)
            .map_err(|error| PyValueError::new_err(error.to_string()))?;
        let value = composition::normalize_closed_root_chain(StructuredAtom::new(atom, interface))
            .map_err(|error| PyValueError::new_err(error.to_string()))?;
        let atom = value.atom;
        let interface = value.interface;
        let interface = merge_explicit_interface_sequence(&[interface])?;
        // A contraction is a derived expression, not a declaration of the original stored data.
        let (name, name_args) = if interface.canonical().order() == original_rank {
            (name, name_args)
        } else {
            (None, Vec::new())
        };
        Py::new(
            py,
            (
                Self {
                    interface: interface.canonicalize_open_ports(),
                    name,
                    name_args,
                },
                PythonExpression { expr: atom },
            ),
        )
    }

    pub(crate) fn from_structured(py: Python<'_>, value: StructuredAtom) -> PyResult<Py<Self>> {
        Self::from_atom_interface(py, value.atom, value.interface)
    }

    fn promoted_network(self_: &PyRef<'_, Self>, py: Python<'_>) -> PyResult<SpensoNet> {
        let expression = Self::from_atom_interface_descriptor(
            py,
            self_.as_super().expr.clone(),
            self_.interface.clone(),
            self_.name,
            self_.name_args.clone(),
        )?;
        SpensoNet::from_arithmetic(ArithmeticStructure::Tensor(expression), None)
    }

    pub(crate) fn structured(self_: &PyRef<'_, Self>) -> StructuredAtom {
        StructuredAtom::new(self_.as_super().expr.clone(), self_.interface.clone())
    }

    pub(crate) fn index_replacements(
        interface: &PartialStructure,
        indices: &Bound<'_, PyTuple>,
        cook_indices: bool,
    ) -> PyResult<HashMap<usize, AbstractIndex>> {
        let open_positions = interface.open_positions();
        if indices.len() != open_positions.len() {
            return Err(PyValueError::new_err(format!(
                "expected {} indices for unresolved ports, got {}",
                open_positions.len(),
                indices.len()
            )));
        }

        let logical = interface.logical_slots();
        let mut replacements = HashMap::new();
        for (argument, position) in indices.iter().zip(open_positions) {
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

    pub(crate) fn descriptor_name(self_: &PyRef<'_, Self>) -> Option<Symbol> {
        self_.name
    }

    pub(crate) fn descriptor_args(self_: &PyRef<'_, Self>) -> Vec<Atom> {
        self_.name_args.clone()
    }

    pub(crate) fn materialized_atom(self_: &PyRef<'_, Self>) -> PyResult<Atom> {
        Ok(Self::materialized(self_)?.atom)
    }

    pub(crate) fn materialized(self_: &PyRef<'_, Self>) -> PyResult<StructuredAtom> {
        let value = Self::structured(self_);
        let replacements = self_
            .interface
            .open_positions()
            .into_iter()
            .map(|position| {
                (
                    position,
                    composition::fresh_dummy_index([&value.atom], [&value.interface]),
                )
            })
            .collect::<HashMap<_, _>>();
        let atom = composition::materialize_interface_ports(&value, &replacements)
            .map_err(|error| PyValueError::new_err(error.to_string()))?;
        let mut logical = self_.interface.logical_slots();
        for (position, index) in replacements {
            logical[position].set_aind(PartialIndex::Explicit(index));
        }
        Ok(StructuredAtom::new(
            atom,
            PartialStructure::from_logical_slots(logical),
        ))
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
        Self::from_atom_interface_descriptor(
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
        let value = value_to_structured_atom(structure)?;
        Self::from_atom_interface_descriptor(
            py,
            value.atom,
            value.interface,
            structure.canonical().name(),
            structure.canonical().args().unwrap_or_default(),
        )
    }
}

impl ModuleInit for TensorExpression {}

fn additive_interfaces_match(left: &PartialStructure, right: &PartialStructure) -> bool {
    // Unresolved ports remain positional; only explicit indices identify ports across terms.
    left.logical_slots() == right.logical_slots()
        || (left.open_positions().is_empty()
            && right.open_positions().is_empty()
            && left.canonical() == right.canonical())
}

fn builtin_tensor_structure(
    symbol: Symbol,
    arguments: &[AtomView<'_>],
) -> PyResult<Option<Canonicalized<ExplicitKey<AbstractIndex>>>> {
    let builtin = symbol == ETS.metric
        || symbol == ETS.flat
        || symbol == AGS.gamma
        || symbol == AGS.gamma5
        || symbol == AGS.projm
        || symbol == AGS.projp
        || symbol == AGS.sigma
        || symbol == CS.f
        || symbol == CS.t;
    if !builtin {
        return Ok(None);
    }

    let compact_metric = symbol == ETS.metric
        && arguments.len() == 2
        && arguments.iter().all(|argument| {
            Slot::<LibraryRep, AbstractIndex>::try_from(*argument).is_err()
                && Representation::<LibraryRep>::try_from(*argument).is_err()
                && !is_chain_placeholder(*argument)
                && has_structured_syntax(*argument)
        });
    if compact_metric {
        return Ok(None);
    }

    let ports = arguments
        .iter()
        .map(|argument| {
            if is_chain_placeholder(*argument) {
                Ok(None)
            } else {
                Slot::<LibraryRep, AbstractIndex>::try_from(*argument)
                    .map(|slot| Some(slot.rep()))
                    .or_else(|_| Representation::<LibraryRep>::try_from(*argument).map(Some))
                    .map_err(|_| ())
            }
        })
        .collect::<Result<Vec<_>, _>>();
    let dimension = |representation: LibraryRep| {
        ports
            .as_ref()
            .ok()
            .and_then(|ports| {
                ports
                    .iter()
                    .flatten()
                    .find(|port| port.rep == representation)
            })
            .map(|port| port.dim)
    };
    let invalid = |factory: &str, signature: &str| {
        PyValueError::new_err(format!(
            "predefined tensor `{factory}` requires {signature}; use \
             TensorExpression.{factory}(...) to construct it"
        ))
    };

    let (factory, signature, expected): (_, _, Option<Canonicalized<ExplicitKey<AbstractIndex>>>) =
        if symbol == ETS.metric || symbol == ETS.flat {
            let factory = if symbol == ETS.metric { "g" } else { "flat" };
            let ports = ports
                .as_ref()
                .map_err(|_| invalid(factory, "two compatible representation ports"))?;
            if arguments.len() != 2 {
                return Err(invalid(factory, "two compatible representation ports"));
            }
            let visible = ports.iter().flatten().copied().collect::<Vec<_>>();
            let compatible = visible.len() < 2
                || if symbol == ETS.metric {
                    visible[0] == visible[1] || visible[0].matches(&visible[1])
                } else {
                    visible[0] == visible[1]
                };
            if !compatible {
                return Err(invalid(factory, "two compatible representation ports"));
            }
            let expected = match visible.as_slice() {
                [left, right] => Some(ExplicitKey::from_iter([*left, *right], symbol, None)),
                _ => None,
            };
            (factory, "two compatible representation ports", expected)
        } else if symbol == AGS.gamma {
            (
                "gamma",
                "one Minkowski and two four-dimensional bispinor ports",
                Some(
                    AGS.gamma_strct(
                        dimension(Minkowski {}.into()).unwrap_or(Dimension::Concrete(1)),
                    ),
                ),
            )
        } else if symbol == AGS.gamma5 || symbol == AGS.projm || symbol == AGS.projp {
            let dimension = dimension(Bispinor {}.into()).unwrap_or(Dimension::Concrete(1));
            let (factory, expected) = if symbol == AGS.gamma5 {
                ("gamma5", AGS.gamma5_strct(dimension))
            } else if symbol == AGS.projm {
                ("projm", AGS.projm_strct(dimension))
            } else {
                ("projp", AGS.projp_strct(dimension))
            };
            (factory, "two equal bispinor ports", Some(expected))
        } else if symbol == AGS.sigma {
            (
                "sigma",
                "two equal Minkowski and two four-dimensional bispinor ports",
                Some(
                    AGS.sigma_strct(
                        dimension(Minkowski {}.into()).unwrap_or(Dimension::Concrete(1)),
                    ),
                ),
            )
        } else if symbol == CS.f {
            (
                "f",
                "three equal adjoint ports",
                Some(
                    CS.f_strct(dimension(ColorAdjoint {}.into()).unwrap_or(Dimension::Concrete(1))),
                ),
            )
        } else if symbol == CS.t {
            let adjoint_dimension =
                dimension(ColorAdjoint {}.into()).unwrap_or(Dimension::Concrete(1));
            let fundamental_dimension = dimension(ColorFundamental {}.into())
                .or_else(|| dimension(ColorAntiFundamental {}.into()))
                .unwrap_or(Dimension::Concrete(1));
            (
                "t",
                "adjoint, fundamental, and antifundamental ports with matching color dimensions",
                Some(CS.t_strct(fundamental_dimension, adjoint_dimension)),
            )
        } else {
            unreachable!("every predefined tensor symbol was handled")
        };

    let ports = ports.map_err(|_| invalid(factory, signature))?;
    let placeholders = arguments
        .iter()
        .enumerate()
        .filter_map(|(position, argument)| {
            let AtomView::Var(variable) = *argument else {
                return None;
            };
            let symbol = variable.get_symbol();
            is_chain_placeholder(*argument).then_some((position, symbol))
        })
        .collect::<Vec<_>>();
    if !placeholders.is_empty()
        && (!matches!(placeholders.as_slice(), [(_, input), (_, output)] if *input != *output)
            || !placeholders
                .iter()
                .any(|(_, symbol)| *symbol == SPENSO_TAG.chain_in)
            || !placeholders
                .iter()
                .any(|(_, symbol)| *symbol == SPENSO_TAG.chain_out))
    {
        return Err(invalid(factory, signature));
    }
    let Some(expected) = expected else {
        return Ok(None);
    };
    let expected_ports = expected
        .canonical()
        .external_reps_iter()
        .collect::<Vec<_>>();
    if ports.len() != expected_ports.len()
        || ports
            .iter()
            .zip(&expected_ports)
            .any(|(actual, expected)| actual.is_some_and(|actual| actual != *expected))
    {
        return Err(invalid(factory, signature));
    }
    if !placeholders.is_empty() {
        let input = placeholders
            .iter()
            .find(|(_, symbol)| *symbol == SPENSO_TAG.chain_in)
            .expect("placeholder names were validated")
            .0;
        let output = placeholders
            .iter()
            .find(|(_, symbol)| *symbol == SPENSO_TAG.chain_out)
            .expect("placeholder names were validated")
            .0;
        let input = expected_ports[input];
        let output = expected_ports[output];
        let oriented = input.rep.is_self_dual() || (input.rep.is_base() && output.rep.is_dual());
        if !input.matches(&output) || !oriented {
            return Err(invalid(factory, signature));
        }
    }
    Ok(Some(expected))
}

fn validate_builtin_placeholder_channels(
    value: AtomView<'_>,
    input: Representation<LibraryRep>,
    output: Representation<LibraryRep>,
) -> PyResult<()> {
    match value {
        AtomView::Add(add) => {
            for term in add.iter() {
                validate_builtin_placeholder_channels(term, input, output)?;
            }
        }
        AtomView::Mul(mul) => {
            for factor in mul.iter() {
                validate_builtin_placeholder_channels(factor, input, output)?;
            }
        }
        AtomView::Pow(power) => {
            let (base, exponent) = power.get_base_exp();
            validate_builtin_placeholder_channels(base, input, output)?;
            validate_builtin_placeholder_channels(exponent, input, output)?;
        }
        AtomView::Fun(function) => {
            let arguments = function.iter().collect::<Vec<_>>();
            if arguments
                .iter()
                .any(|argument| is_chain_placeholder(*argument))
            {
                let resolved = arguments
                    .iter()
                    .map(|argument| match *argument {
                        AtomView::Var(variable) if variable.get_symbol() == SPENSO_TAG.chain_in => {
                            input.to_symbolic([])
                        }
                        AtomView::Var(variable)
                            if variable.get_symbol() == SPENSO_TAG.chain_out =>
                        {
                            output.to_symbolic([])
                        }
                        argument => argument.to_owned(),
                    })
                    .collect::<Vec<_>>();
                let resolved = resolved.iter().map(Atom::as_view).collect::<Vec<_>>();
                builtin_tensor_structure(function.get_symbol(), &resolved)?;
            }
            for argument in arguments {
                if !is_chain_placeholder(argument) {
                    validate_builtin_placeholder_channels(argument, input, output)?;
                }
            }
        }
        _ => {}
    }
    Ok(())
}

fn from_idenso_atom(
    self_: &PyRef<'_, TensorExpression>,
    py: Python<'_>,
    atom: Atom,
) -> PyResult<Py<TensorExpression>> {
    let interface = if atom.as_view().is_zero() {
        self_.interface.clone()
    } else if has_structured_syntax(atom.as_view()) {
        reinfer_structured(atom.clone())?.interface
    } else if self_.interface.canonical().is_scalar() {
        PartialStructure::from_logical_slots(std::iter::empty())
    } else {
        return Err(PyValueError::new_err(
            "Idenso transformation removed the tensor syntax of a non-scalar expression",
        ));
    };
    TensorExpression::from_atom_interface_descriptor(
        py,
        atom,
        interface,
        self_.name,
        self_.name_args.clone(),
    )
}

fn idenso_input(self_: &PyRef<'_, TensorExpression>) -> PythonExpression {
    PythonExpression {
        expr: self_.as_super().expr.clone(),
    }
}

pub(crate) enum TensorOperand {
    Structured(StructuredAtom),
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
            return if has_structured_syntax(atom.as_view()) {
                reinfer_structured(atom).map(Self::Structured)
            } else {
                Ok(Self::Scalar(atom))
            };
        }
        Err(PyTypeError::new_err(
            "expected a TensorExpression or Symbolica expression",
        ))
    }
}

pub(crate) fn value_to_structured_atom(
    value: &Canonicalized<ExplicitKey<AbstractIndex>>,
) -> PyResult<StructuredAtom> {
    let structure = value.canonical();
    let name = structure
        .global_name
        .ok_or_else(|| PyRuntimeError::new_err("tensor structure has no name"))?;
    let args = structure.additional_args.as_deref().unwrap_or_default();
    let canonical = value.canonical().external_reps_iter().collect::<Vec<_>>();
    let logical = value.layout().canonical_to_logical(&canonical);
    let logical_axes = (0..logical.len()).collect::<Vec<_>>();
    let canonical_axes = value.layout().logical_to_canonical(&logical_axes);
    let owner = crate::fresh_open_owner();
    let atom = FunctionBuilder::new(name)
        .add_args(args)
        .add_args(
            canonical
                .iter()
                .zip(canonical_axes)
                .map(|(representation, axis)| {
                    representation
                        .slot::<AbstractIndex, _>(AbstractIndex::Open { owner, axis })
                        .to_atom()
                }),
        )
        .finish();
    let interface = PartialStructure::from_logical_slots(
        logical
            .into_iter()
            .enumerate()
            .map(|(position, representation)| representation.slot(PartialIndex::open(position))),
    );
    Ok(StructuredAtom::new(atom, interface))
}

fn index_value(value: ConvertibleToAbstractIndex, cook: bool) -> PyResult<AbstractIndex> {
    match value {
        ConvertibleToAbstractIndex::Aind(index) => Ok(index),
        ConvertibleToAbstractIndex::Atom(expression) => {
            let atom = if cook {
                expression.expr.cook_indices()
            } else {
                expression.expr
            };
            atom.as_view().try_into().map_err(|error| {
                let hint = if cook {
                    ""
                } else {
                    " Try setting cook_indices=True."
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

fn has_structured_syntax(value: AtomView<'_>) -> bool {
    value.is_tensorial(StrictTensorFilter::Tagged)
        || is_unmaterialized_tensor_leaf(value)
        || matches!(
            value,
            AtomView::Fun(projector)
                if projector.get_symbol() == *shadowing::SYM
                    || projector.get_symbol() == *shadowing::ANTISYM
                    || projector.get_symbol() == *shadowing::CYCLIC
        )
}

fn is_unmaterialized_tensor_leaf(value: AtomView<'_>) -> bool {
    let AtomView::Fun(function) = value else {
        return false;
    };
    function.get_symbol().has_tag(&SPENSO_TAG.tensor)
        && function.iter().any(is_chain_placeholder)
        && function.iter().all(|argument| {
            Slot::<LibraryRep, AbstractIndex>::try_from(argument).is_err()
                && Representation::<LibraryRep>::try_from(argument).is_err()
        })
}

fn is_chain_placeholder(value: AtomView<'_>) -> bool {
    matches!(
        value,
        AtomView::Var(variable)
            if variable.get_symbol() == SPENSO_TAG.chain_in
                || variable.get_symbol() == SPENSO_TAG.chain_out
    )
}

/// Namespaced `in`/`out` placeholders are wiring labels, not free abstract indices. They may
/// only occur as direct ports of tensor leaves inside a chain or trace factor.
fn validate_placeholder_scope(value: AtomView<'_>, inside_factor: bool) -> Result<(), String> {
    match value {
        AtomView::Var(_) if is_chain_placeholder(value) => Err(
            "Spenso chain placeholders are only valid inside chain or trace tensor factors".into(),
        ),
        AtomView::Add(sum) => {
            for term in sum.iter() {
                validate_placeholder_scope(term, inside_factor)?;
            }
            Ok(())
        }
        AtomView::Mul(product) => {
            for factor in product.iter() {
                validate_placeholder_scope(factor, inside_factor)?;
            }
            Ok(())
        }
        AtomView::Pow(power) => {
            let (base, exponent) = power.get_base_exp();
            validate_placeholder_scope(base, inside_factor)?;
            validate_placeholder_scope(exponent, inside_factor)
        }
        AtomView::Fun(function) => {
            let symbol = function.get_symbol();
            if symbol == SPENSO_TAG.chain {
                for (position, argument) in function.iter().enumerate() {
                    validate_placeholder_scope(argument, position >= 2)?;
                }
                return Ok(());
            }
            if symbol == SPENSO_TAG.trace {
                for (position, argument) in function.iter().enumerate() {
                    validate_placeholder_scope(argument, position >= 1)?;
                }
                return Ok(());
            }

            let tensor_leaf = symbol.has_tag(&SPENSO_TAG.tensor);
            let mut inputs = 0;
            let mut outputs = 0;
            for argument in function.iter() {
                if is_chain_placeholder(argument) {
                    if !(inside_factor && tensor_leaf) {
                        return Err(
                            "Spenso chain placeholders are only valid as tensor ports inside chain or trace factors"
                                .into(),
                        );
                    }
                    let AtomView::Var(variable) = argument else {
                        unreachable!("chain placeholders are variables")
                    };
                    if variable.get_symbol() == SPENSO_TAG.chain_in {
                        inputs += 1;
                    } else {
                        outputs += 1;
                    }
                } else {
                    validate_placeholder_scope(argument, inside_factor)?;
                }
            }
            if inputs + outputs > 0 && (inputs != 1 || outputs != 1) {
                return Err(format!(
                    "tensor factor `{symbol}` requires exactly one direct `in` and one direct `out` placeholder"
                ));
            }
            Ok(())
        }
        _ => Ok(()),
    }
}

fn infer_interface(atom: &Atom) -> PyResult<PartialStructure> {
    atom.validate_chain_like_nesting()
        .map_err(|error| PyValueError::new_err(error.to_string()))?;
    validate_placeholder_scope(atom.as_view(), false).map_err(PyValueError::new_err)?;

    let mut syntax_error = None;
    let _ = atom.replace_map(|value, _, _| {
        if syntax_error.is_some() {
            return;
        }

        match value {
            AtomView::Var(variable) if variable.get_symbol().has_tag(&SPENSO_TAG.tensor) => {
                syntax_error = Some(format!(
                    "tensor symbol `{}` must be called with at least one structural port",
                    variable.get_symbol()
                ));
            }
            AtomView::Fun(function) => {
                let symbol = function.get_symbol();
                let arguments = function.iter().collect::<Vec<_>>();
                if symbol == SPENSO_TAG.bracket && arguments.is_empty() {
                    syntax_error = Some("ordered tensor products require at least one operand".into());
                    return;
                }
                if symbol == SPENSO_TAG.dot && arguments.len() != 2 {
                    syntax_error = Some("dot requires exactly two operands".into());
                    return;
                }
                if symbol == SPENSO_TAG.chain && arguments.len() < 2 {
                    syntax_error = Some("chain requires explicit start and end ports".into());
                    return;
                }
                if symbol == SPENSO_TAG.trace {
                    let Some(representation) = arguments.first() else {
                        syntax_error = Some("trace requires a representation".into());
                        return;
                    };
                    if Representation::<LibraryRep>::try_from(*representation).is_err() {
                        syntax_error =
                            Some("trace metadata is not a Spenso representation".into());
                        return;
                    }
                }
                if symbol.has_tag(&SPENSO_TAG.broadcast) && arguments.len() != 1 {
                    syntax_error = Some(format!(
                        "broadcast function `{symbol}` requires exactly one argument"
                    ));
                    return;
                }
                if !symbol.has_tag(&SPENSO_TAG.tensor) {
                    return;
                }

                let compact_metric = symbol == ETS.metric
                    && arguments.len() == 2
                    && arguments.iter().all(|argument| {
                        Slot::<LibraryRep, AbstractIndex>::try_from(*argument).is_err()
                            && Representation::<LibraryRep>::try_from(*argument).is_err()
                            && has_structured_syntax(*argument)
                    });
                if let Err(error) = builtin_tensor_structure(symbol, &arguments) {
                    syntax_error = Some(error.to_string());
                    return;
                }

                let ports = arguments
                    .iter()
                    .enumerate()
                    .filter(|(_, argument)| {
                        Slot::<LibraryRep, AbstractIndex>::try_from(**argument).is_ok()
                            || Representation::<LibraryRep>::try_from(**argument).is_ok()
                            || matches!(
                                **argument,
                                AtomView::Var(variable)
                                    if variable.get_symbol() == SPENSO_TAG.chain_in
                                        || variable.get_symbol() == SPENSO_TAG.chain_out
                            )
                    })
                    .map(|(position, _)| position)
                    .collect::<Vec<_>>();
                let ports_are_final = ports
                    .iter()
                    .copied()
                    .eq(arguments.len().saturating_sub(ports.len())..arguments.len());
                if !ports_are_final && !compact_metric {
                    syntax_error = Some(format!(
                        "tensor function `{symbol}` requires scalar arguments before structural ports"
                    ));
                } else if symbol.has_tag(&SPENSO_TAG.rank1) && ports.len() != 1 {
                    syntax_error = Some(format!(
                        "rank-one tensor function `{symbol}` requires exactly one final structural port"
                    ));
                }
            }
            _ => {}
        }
    });
    if let Some(error) = syntax_error {
        return Err(PyValueError::new_err(error));
    }

    infer_validated_interface(atom)
}

fn lower_tensor_powers(value: AtomView<'_>) -> PyResult<Atom> {
    match value {
        AtomView::Add(sum) => {
            let mut result = Atom::Zero;
            for term in sum.iter() {
                result += lower_tensor_powers(term)?;
            }
            Ok(result)
        }
        AtomView::Mul(product) => {
            let mut result = Atom::num(1);
            for factor in product.iter() {
                result *= lower_tensor_powers(factor)?;
            }
            Ok(result)
        }
        AtomView::Pow(power) => {
            let (base, exponent) = power.get_base_exp();
            let base = lower_tensor_powers(base)?;
            let exponent = lower_tensor_powers(exponent)?;
            if !has_structured_syntax(base.as_view()) {
                return Ok(base.pow(exponent));
            }

            // Validate multiplicity and exponent constraints against the
            // unexpanded spelling before lowering it through tensor-aware
            // multiplication.
            infer_validated_interface(&base.as_ref().pow(exponent.as_ref()))?;
            let raw_interface = infer_interface(&base)?;
            let interface = merge_explicit_interface_sequence(&[raw_interface])?;
            if interface.canonical().is_scalar() {
                return Ok(base.pow(exponent));
            }

            let exponent = Rational::try_from(exponent.as_view()).map_err(|_| {
                PyValueError::new_err("a non-scalar tensor power requires an integer exponent")
            })?;
            if exponent.denominator() != 1 || exponent.numerator().is_negative() {
                return Err(PyValueError::new_err(
                    "a non-scalar tensor power requires a non-negative integer exponent",
                ));
            }
            let repetitions = usize::try_from(exponent.numerator().clone()).map_err(|_| {
                PyValueError::new_err("tensor power exponent is too large to materialize")
            })?;
            if repetitions == 0 {
                return Ok(Atom::num(1));
            }

            let factor = StructuredAtom::new(base, interface);
            let mut result = factor.clone();
            for _ in 1..repetitions {
                result = composition::multiply(&result, &factor)
                    .map_err(|error| PyValueError::new_err(error.to_string()))?;
            }
            Ok(result.atom)
        }
        AtomView::Fun(function) => {
            let mut result = FunctionBuilder::new(function.get_symbol());
            for argument in function.iter() {
                result = result.add_arg(lower_tensor_powers(argument)?);
            }
            Ok(result.finish())
        }
        _ => Ok(value.to_owned()),
    }
}

fn reinfer_structured(atom: Atom) -> PyResult<StructuredAtom> {
    let atom = lower_tensor_powers(atom.as_view())?;
    let interface = infer_interface(&atom)?;
    Ok(StructuredAtom::new(atom, interface))
}

fn infer_validated_interface(atom: &Atom) -> PyResult<PartialStructure> {
    if let AtomView::Add(sum) = atom.as_view() {
        let mut expected = None;
        let mut has_scalar_term = false;
        for term in sum.iter() {
            if !has_structured_syntax(term) {
                has_scalar_term = true;
                continue;
            }
            let actual =
                merge_explicit_interface_sequence(&[infer_validated_interface(&term.to_owned())?])?;
            let Some(current) = &expected else {
                expected = Some(actual);
                continue;
            };
            if !additive_interfaces_match(current, &actual) {
                return Err(PyValueError::new_err(
                    "tensor summands do not have compatible tensor interfaces",
                ));
            }
        }
        let Some(expected) = expected else {
            return Err(PyValueError::new_err(
                "expression does not contain valid tagged Spenso tensor syntax",
            ));
        };
        if has_scalar_term && !expected.canonical().is_scalar() {
            return Err(PyValueError::new_err(
                "cannot add a scalar expression to a non-scalar tensor",
            ));
        }
        return Ok(expected);
    }

    if let AtomView::Mul(product) = atom.as_view() {
        let mut interfaces = Vec::new();
        for factor in product.iter() {
            if !has_structured_syntax(factor) {
                continue;
            }
            interfaces.push(infer_validated_interface(&factor.to_owned())?);
        }
        if interfaces.is_empty() {
            return Err(PyValueError::new_err(
                "expression does not contain valid tagged Spenso tensor syntax",
            ));
        }
        return merge_explicit_interface_sequence(&interfaces);
    }

    if let AtomView::Pow(power) = atom.as_view() {
        let (base, exponent) = power.get_base_exp();
        if !has_structured_syntax(base) {
            return Err(PyValueError::new_err(
                "expression does not contain valid tagged Spenso tensor syntax",
            ));
        }
        let interface = infer_validated_interface(&base.to_owned())?;
        if has_structured_syntax(exponent)
            && !infer_validated_interface(&exponent.to_owned())?
                .canonical()
                .is_scalar()
        {
            return Err(PyValueError::new_err("a tensor exponent must be scalar"));
        }
        if interface.canonical().is_scalar() {
            return Ok(interface);
        }
        if !interface
            .logical_slots()
            .iter()
            .all(|slot| slot.rep().rep.is_self_dual())
        {
            return Err(PyValueError::new_err(format!(
                "invalid power of non-self-dual tensor `{atom}`"
            )));
        }
        let exponent = Rational::try_from(exponent)
            .map_err(|_| PyValueError::new_err(format!("invalid tensor power `{atom}`")))?;
        if exponent.denominator() != 1 {
            return Err(PyValueError::new_err(format!(
                "fractional tensor power `{atom}` has no well-defined interface"
            )));
        }
        let repetitions = exponent.numerator().abs();
        let slots = interface.logical_slots();
        for slot in &slots {
            let PartialIndex::Explicit(index) = slot.aind else {
                continue;
            };
            let compatible = slots
                .iter()
                .filter(|candidate| {
                    matches!(candidate.aind, PartialIndex::Explicit(candidate_index) if candidate_index == index)
                        && slot.rep().matches(&candidate.rep())
                })
                .count();
            let exceeds_einstein_multiplicity = match repetitions.to_i64() {
                Some(0) => false,
                Some(1) => compatible > 2,
                Some(2) => compatible > 1,
                _ => compatible > 0,
            };
            if exceeds_einstein_multiplicity {
                return Err(PyValueError::new_err(format!(
                    "explicit index `{index}` occurs on more than two compatible ports in `{atom}`"
                )));
            }
        }
        if exponent.numerator() % 2 == 0 {
            return Ok(PartialStructure::from_logical_slots([]));
        }
        return Ok(interface);
    }

    if let AtomView::Fun(function) = atom.as_view() {
        let symbol = function.get_symbol();
        let arguments = function.iter().collect::<Vec<_>>();
        let unmaterialized_tensor_leaf = is_unmaterialized_tensor_leaf(atom.as_view());
        if symbol == *shadowing::SYM
            || symbol == *shadowing::ANTISYM
            || symbol == *shadowing::CYCLIC
        {
            let mut interfaces = Vec::new();
            for argument in arguments {
                if !has_structured_syntax(argument) {
                    continue;
                }
                interfaces.push(infer_validated_interface(&argument.to_owned())?);
            }
            if interfaces.is_empty() {
                return Err(PyValueError::new_err(
                    "tensor projectors require at least one tensor operand",
                ));
            }
            return merge_explicit_interface_sequence(&interfaces);
        }

        if symbol == SPENSO_TAG.bracket {
            let mut interfaces = Vec::new();
            for argument in arguments {
                if !has_structured_syntax(argument) {
                    continue;
                }
                interfaces.push(infer_validated_interface(&argument.to_owned())?);
            }
            if interfaces.is_empty() {
                return Err(PyValueError::new_err(
                    "ordered tensor products require at least one tensor operand",
                ));
            }
            return merge_explicit_interface_sequence(&interfaces);
        }

        let compact_metric = symbol == ETS.metric
            && arguments.len() == 2
            && arguments.iter().all(|argument| {
                Slot::<LibraryRep, AbstractIndex>::try_from(*argument).is_err()
                    && Representation::<LibraryRep>::try_from(*argument).is_err()
                    && !is_chain_placeholder(*argument)
                    && has_structured_syntax(*argument)
            });
        if symbol == SPENSO_TAG.dot || compact_metric {
            let [left, right] = arguments.as_slice() else {
                return Err(PyValueError::new_err(
                    "inner products require exactly two operands",
                ));
            };
            if !has_structured_syntax(*left) || !has_structured_syntax(*right) {
                return Err(PyValueError::new_err(
                    "inner products require rank-one tensor operands",
                ));
            }
            let left = infer_validated_interface(&left.to_owned())?;
            let right = infer_validated_interface(&right.to_owned())?;
            if left.canonical().order() != 1 || right.canonical().order() != 1 {
                return Err(PyValueError::new_err(format!(
                    "inner products require rank-one operands, got ranks {} and {}",
                    left.canonical().order(),
                    right.canonical().order()
                )));
            }
            let left = left.logical_slots()[0];
            let right = right.logical_slots()[0];
            if !left.rep().matches(&right.rep()) {
                return Err(PyValueError::new_err(
                    "inner-product operands carry incompatible representations",
                ));
            }
            if matches!(
                (left.aind, right.aind),
                (PartialIndex::Explicit(left), PartialIndex::Explicit(right)) if left != right
            ) {
                return Err(PyValueError::new_err(
                    "inner-product operands carry unequal explicit indices",
                ));
            }
            return Ok(PartialStructure::from_logical_slots([]));
        }

        if symbol.has_tag(&SPENSO_TAG.broadcast) {
            let argument = arguments[0];
            if !has_structured_syntax(argument) {
                return Err(PyValueError::new_err(format!(
                    "broadcast function `{symbol}` does not contain a structured tensor argument"
                )));
            }
            return infer_validated_interface(&argument.to_owned());
        }

        if let Some(structure) = builtin_tensor_structure(symbol, &arguments)? {
            let canonical_ports = arguments
                .iter()
                .enumerate()
                .map(|(position, argument)| {
                    if let Ok(slot) = Slot::<LibraryRep, AbstractIndex>::try_from(*argument) {
                        let index = match slot.aind() {
                            AbstractIndex::Open { axis, .. } => PartialIndex::open(axis),
                            index => PartialIndex::Explicit(index),
                        };
                        Some(slot.rep().slot(index))
                    } else {
                        Representation::<LibraryRep>::try_from(*argument)
                            .ok()
                            .map(|representation| representation.slot(PartialIndex::open(position)))
                    }
                })
                .collect::<Option<Vec<_>>>();
            if let Some(canonical_ports) = canonical_ports {
                return Ok(PartialStructure::from_logical_slots(
                    structure.layout().canonical_to_logical(&canonical_ports),
                ));
            }
        }

        if symbol == SPENSO_TAG.chain {
            let endpoints = arguments[..2]
                .iter()
                .map(|argument| {
                    if let Ok(slot) = Slot::<LibraryRep, AbstractIndex>::try_from(*argument) {
                        Ok(slot.rep().slot(PartialIndex::Explicit(slot.aind())))
                    } else if let Ok(representation) =
                        Representation::<LibraryRep>::try_from(*argument)
                    {
                        Ok(representation.slot(PartialIndex::open(0)))
                    } else {
                        Err(PyValueError::new_err(
                            "chain endpoints must be Spenso slots or representations",
                        ))
                    }
                })
                .collect::<PyResult<Vec<_>>>()?;
            let input = endpoints[0].rep();
            let output = endpoints[1].rep();
            if !input.matches(&output)
                || !(input.rep.is_self_dual() || (input.rep.is_base() && output.rep.is_dual()))
            {
                return Err(PyValueError::new_err(
                    "chain endpoints do not form a compatible input-to-output channel",
                ));
            }
            let mut interfaces = Vec::new();
            for factor in &arguments[2..] {
                if !has_structured_syntax(*factor) {
                    return Err(PyValueError::new_err(
                        "chain factors must contain structured tensors",
                    ));
                }
                validate_builtin_placeholder_channels(*factor, input, output)?;
                interfaces.push(infer_validated_interface(&factor.to_owned())?);
            }
            let spectators = merge_explicit_interface_sequence(&interfaces)?;
            return Ok(PartialStructure::from_logical_slots(
                endpoints.into_iter().chain(spectators.logical_slots()),
            ));
        } else if symbol == SPENSO_TAG.trace {
            let representation =
                Representation::<LibraryRep>::try_from(arguments[0]).map_err(|_| {
                    PyValueError::new_err("trace metadata is not a Spenso representation")
                })?;
            let dual = representation.dual();
            let mut interfaces = Vec::new();
            for factor in shadowing::trace_factor_views(&arguments[1..]) {
                if !has_structured_syntax(factor) {
                    return Err(PyValueError::new_err(
                        "trace factors must contain structured tensors",
                    ));
                }
                validate_builtin_placeholder_channels(factor, representation, dual)?;
                interfaces.push(infer_validated_interface(&factor.to_owned())?);
            }
            return merge_explicit_interface_sequence(&interfaces);
        } else if symbol.has_tag(&SPENSO_TAG.tensor) {
            for argument in arguments {
                if Slot::<LibraryRep, AbstractIndex>::try_from(argument).is_ok()
                    || Representation::<LibraryRep>::try_from(argument).is_ok()
                    || matches!(
                        argument,
                        AtomView::Var(variable)
                            if variable.get_symbol() == SPENSO_TAG.chain_in
                                || variable.get_symbol() == SPENSO_TAG.chain_out
                    )
                    || !has_structured_syntax(argument)
                {
                    continue;
                }
                if !infer_validated_interface(&argument.to_owned())?
                    .canonical()
                    .is_scalar()
                {
                    return Err(PyValueError::new_err(format!(
                        "tensor metadata for `{symbol}` must be scalar"
                    )));
                }
            }
            if unmaterialized_tensor_leaf {
                return Ok(PartialStructure::from_logical_slots([]));
            }
        }
    }

    if !atom.is_tensorial(StrictTensorFilter::Tagged) {
        return Err(PyValueError::new_err(
            "expression does not contain valid tagged Spenso tensor syntax",
        ));
    }

    let mut open_markers = HashSet::new();
    let materialized = atom.replace_map(|value, _, output| {
        if Slot::<LibraryRep, AbstractIndex>::try_from(value).is_ok() {
            return;
        }
        if let Ok(representation) = Representation::<LibraryRep>::try_from(value) {
            let marker = loop {
                let marker = composition::fresh_dummy_index([atom], std::iter::empty());
                if open_markers.insert(marker) {
                    break marker;
                }
            };
            **output = representation.slot::<AbstractIndex, _>(marker).to_atom();
        }
    });
    let inferred = match materialized
        .infer_structure::<OrderedStructure<LibraryRep, AbstractIndex>>(
            StructureInferenceMode::Fast,
        ) {
        Ok(inferred) => inferred,
        Err(StructureError::EmptyStructure(_)) => {
            return Ok(PartialStructure::from_logical_slots([]));
        }
        Err(error) => {
            return Err(PyValueError::new_err(format!(
                "invalid Spenso expression: {error}"
            )));
        }
    };
    // `OrderedStructure` canonicalizes slots and its fast function inference
    // does not retain that permutation. Read the direct structural arguments
    // back from the atom so ordinary non-built-in leaves follow their encoded
    // syntax order.
    let logical = syntactic_leaf_slots(materialized.as_view());
    if logical.len() != inferred.canonical().order() {
        return Err(PyValueError::new_err(format!(
            "invalid Spenso expression: inferred {} ports but found {} direct structural arguments",
            inferred.canonical().order(),
            logical.len()
        )));
    }
    Ok(PartialStructure::from_logical_slots(
        logical.into_iter().map(|slot| {
            let index = match slot.aind() {
                AbstractIndex::Open { axis, .. } => PartialIndex::open(axis),
                index if open_markers.contains(&index) => PartialIndex::open(0),
                index => PartialIndex::Explicit(index),
            };
            slot.rep().slot(index)
        }),
    ))
}

fn syntactic_leaf_slots(value: AtomView<'_>) -> Vec<Slot<LibraryRep, AbstractIndex>> {
    if let Ok(slot) = Slot::<LibraryRep, AbstractIndex>::try_from(value) {
        return vec![slot];
    }
    let AtomView::Fun(function) = value else {
        return Vec::new();
    };
    function
        .iter()
        .filter_map(|argument| Slot::<LibraryRep, AbstractIndex>::try_from(argument).ok())
        .collect()
}

fn merge_explicit_interface_sequence(
    interfaces: &[PartialStructure],
) -> PyResult<PartialStructure> {
    let slots = interfaces
        .iter()
        .flat_map(PartialStructureExt::logical_slots)
        .collect::<Vec<_>>();
    let mut contracted = HashSet::new();
    for (position, slot) in slots.iter().enumerate() {
        let PartialIndex::Explicit(index) = slot.aind else {
            continue;
        };
        let compatible = slots
            .iter()
            .enumerate()
            .filter(|(candidate_position, candidate)| {
                if *candidate_position == position {
                    return true;
                }
                matches!(candidate.aind, PartialIndex::Explicit(candidate_index) if candidate_index == index)
                    && slot.rep().matches(&candidate.rep())
            })
            .map(|(candidate_position, _)| candidate_position)
            .collect::<Vec<_>>();
        if compatible.len() > 2 {
            return Err(PyValueError::new_err(format!(
                "explicit index `{index}` occurs on more than two compatible ports {compatible:?}"
            )));
        }
        if compatible.len() == 2 {
            contracted.extend(compatible);
        }
    }

    Ok(PartialStructure::from_logical_slots(
        slots
            .into_iter()
            .enumerate()
            .filter(|(position, _)| !contracted.contains(position))
            .map(|(_, slot)| slot),
    ))
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
    /// Create an unresolved metric tensor for `rep`.
    ///
    /// Call the result with two indices to fill its ports in logical order.
    #[staticmethod]
    fn g(py: Python<'_>, rep: &SpensoRepresentation) -> PyResult<Py<TensorExpression>> {
        let structure = ExplicitKey::<AbstractIndex>::from_iter(
            [rep.representation, rep.representation],
            ETS.metric,
            None,
        );
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
    /// Its logical ports are Minkowski, bispinor, and bispinor. The bispinor
    /// dimension is four. Call the result with `(mu, i, j)` to index it.
    #[staticmethod]
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

    /// Ordered external interface as concrete `Slot` objects or unresolved `Representation`s.
    #[getter]
    fn interface(&self, py: Python<'_>) -> PyResult<Py<PyTuple>> {
        let values = self
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
        Ok(PyTuple::new(py, values)?.unbind())
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
        Self::from_atom_interface_descriptor(
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
        let atom = self_.as_super().expr.clone();
        let (name, name_args) = inferred_descriptor(atom.as_view());
        let value = reinfer_structured(atom)?;
        Self::from_atom_interface_descriptor(py, value.atom, value.interface, name, name_args)
    }

    /// Expand scalar algebra while preserving and validating the tensor interface.
    #[pyo3(signature = (var = None, via_poly = None))]
    fn expand(
        self_: PyRef<'_, Self>,
        py: Python<'_>,
        var: Option<ConvertibleToExpression>,
        via_poly: Option<bool>,
    ) -> PyResult<Py<Self>> {
        let atom = if let Some(var) = var {
            let var = var.to_expression();
            if !matches!(var.expr, Atom::Var(_) | Atom::Fun(_)) {
                return Err(PyValueError::new_err(
                    "Expansion must be done wrt an indeterminate",
                ));
            }
            if via_poly.unwrap_or(false) {
                self_
                    .as_super()
                    .expr
                    .as_view()
                    .expand_via_poly::<i16, Atom>(Some(var.expr.clone()))
            } else {
                self_
                    .as_super()
                    .expr
                    .as_view()
                    .expand_in(var.expr.as_view())
            }
        } else if via_poly.unwrap_or(false) {
            self_
                .as_super()
                .expr
                .as_view()
                .expand_via_poly::<i16, Atom>(None)
        } else {
            self_.as_super().expr.expand()
        };
        let inferred = if atom.as_view().is_zero() {
            self_.interface.clone()
        } else {
            merge_explicit_interface_sequence(&[infer_interface(&atom)?])?.canonicalize_open_ports()
        };
        if !additive_interfaces_match(&self_.interface, &inferred) {
            return Err(PyValueError::new_err(
                "expanded expression does not preserve a compatible tensor interface",
            ));
        }
        Self::from_atom_interface_descriptor(
            py,
            atom,
            self_.interface.clone(),
            self_.name,
            self_.name_args.clone(),
        )
    }

    /// Apply Idenso's gamma-algebra simplifier and re-infer the tensor interface.
    #[pyo3(signature = (settings = None))]
    fn simplify_gamma(
        self_: PyRef<'_, Self>,
        py: Python<'_>,
        settings: Option<&PyGammaSimplifySettings>,
    ) -> PyResult<Py<Self>> {
        let result = idenso::python::simplify_gamma(&idenso_input(&self_), settings);
        from_idenso_atom(&self_, py, result.expr)
    }

    /// Convert bispinor tensors into chain/trace shorthands and join adjacent gamma chains.
    fn collect_gamma_chains(self_: PyRef<'_, Self>, py: Python<'_>) -> PyResult<Py<Self>> {
        let result = idenso::python::collect_gamma_chains(&idenso_input(&self_));
        from_idenso_atom(&self_, py, result.expr)
    }

    /// Simplify products and linear combinations involving the time-like gamma matrix `gamma0`.
    fn simplify_gamma0(self_: PyRef<'_, Self>, py: Python<'_>) -> PyResult<Py<Self>> {
        let result = idenso::python::simplify_gamma0(&idenso_input(&self_));
        from_idenso_atom(&self_, py, result.expr)
    }

    /// Rewrite conjugated gamma matrices as gamma0-sandwiched matrices.
    fn simplify_gamma_conjugate(self_: PyRef<'_, Self>, py: Python<'_>) -> PyResult<Py<Self>> {
        let result = idenso::python::simplify_gamma_conjugate(&idenso_input(&self_))?;
        from_idenso_atom(&self_, py, result.expr)
    }

    /// Simplify Levi-Civita/metric contractions and pairs of Levi-Civita tensors.
    fn simplify_epsilon(self_: PyRef<'_, Self>, py: Python<'_>) -> PyResult<Py<Self>> {
        let result = idenso::python::simplify_epsilon(&idenso_input(&self_));
        from_idenso_atom(&self_, py, result.expr)
    }

    /// Contract metric and identity tensors and re-infer the external interface.
    fn simplify_metrics(self_: PyRef<'_, Self>, py: Python<'_>) -> PyResult<Py<Self>> {
        let result = idenso::python::simplify_metrics(&idenso_input(&self_));
        from_idenso_atom(&self_, py, result.expr)
    }

    /// Apply Idenso's SU(N) color-algebra simplifier and re-infer the interface.
    #[pyo3(signature = (settings = None))]
    fn simplify_color(
        self_: PyRef<'_, Self>,
        py: Python<'_>,
        settings: Option<&PyColorSimplifySettings>,
    ) -> PyResult<Py<Self>> {
        let result = idenso::python::simplify_color(&idenso_input(&self_), settings);
        from_idenso_atom(&self_, py, result.expr)
    }

    /// Factor around tensors carrying fundamental, antifundamental, or adjoint color.
    fn collect_color(self_: PyRef<'_, Self>, py: Python<'_>) -> PyResult<Py<Self>> {
        let result = idenso::python::collect_color(&idenso_input(&self_));
        from_idenso_atom(&self_, py, result.expr)
    }

    /// Factor around recognized scalar color invariants such as Casimirs and indices.
    fn collect_color_constants(self_: PyRef<'_, Self>, py: Python<'_>) -> PyResult<Py<Self>> {
        let result = idenso::python::collect_color_constants(&idenso_input(&self_));
        from_idenso_atom(&self_, py, result.expr)
    }

    /// Rewrite supplied color dimensions and invariants into a representation-aware Casimir basis.
    #[pyo3(signature = (*, fundamental, adjoint, settings = None))]
    fn to_color_casimir(
        self_: PyRef<'_, Self>,
        py: Python<'_>,
        fundamental: &SpensoRepresentation,
        adjoint: &SpensoRepresentation,
        settings: Option<&PyColorCasimirSettings>,
    ) -> PyResult<Py<Self>> {
        let result = idenso::python::to_color_casimir(
            &idenso_input(&self_),
            RegisteredRepresentation(fundamental.representation),
            RegisteredRepresentation(adjoint.representation),
            settings,
        );
        from_idenso_atom(&self_, py, result.expr)
    }

    /// Replace supported `cof(N)` invariants by explicit dimension formulas.
    fn to_cof_dimension_invariants(self_: PyRef<'_, Self>, py: Python<'_>) -> PyResult<Py<Self>> {
        let result = idenso::python::to_cof_dimension_invariants(&idenso_input(&self_));
        from_idenso_atom(&self_, py, result.expr)
    }

    /// Expand around color structures and wrap each scalar coefficient with `symbol`.
    fn wrap_color(self_: PyRef<'_, Self>, py: Python<'_>, symbol: Symbol) -> PyResult<Py<Self>> {
        let result = idenso::python::wrap_color(&idenso_input(&self_), symbol);
        from_idenso_atom(&self_, py, result.expr)
    }

    /// Selectively expand around Minkowski structures into `(structure, coefficient)` pairs.
    fn expand_mink(self_: PyRef<'_, Self>) -> Vec<idenso::python::PythonTerm> {
        idenso::python::expand_mink(&idenso_input(&self_))
    }

    /// Selectively expand around bispinor structures into `(structure, coefficient)` pairs.
    fn expand_bis(self_: PyRef<'_, Self>) -> Vec<idenso::python::PythonTerm> {
        idenso::python::expand_bis(&idenso_input(&self_))
    }

    /// Selectively expand around Minkowski and bispinor structures into factorized pairs.
    fn expand_mink_bis(self_: PyRef<'_, Self>) -> Vec<idenso::python::PythonTerm> {
        idenso::python::expand_mink_bis(&idenso_input(&self_))
    }

    /// Selectively expand around metric tensors into `(structure, coefficient)` pairs.
    fn expand_metrics(self_: PyRef<'_, Self>) -> Vec<idenso::python::PythonTerm> {
        idenso::python::expand_metrics(&idenso_input(&self_))
    }

    /// Selectively expand around color structures into `(structure, coefficient)` pairs.
    fn expand_color(self_: PyRef<'_, Self>) -> Vec<idenso::python::PythonTerm> {
        idenso::python::expand_color(&idenso_input(&self_))
    }

    /// Selectively expand around the supplied Symbolica patterns into factorized pairs.
    fn expand_in_patterns(
        self_: PyRef<'_, Self>,
        patterns: Vec<PythonExpression>,
    ) -> Vec<idenso::python::PythonTerm> {
        idenso::python::expand_in_patterns(&idenso_input(&self_), patterns)
    }

    /// Wrap every abstract-index payload with `header` and return an ordinary expression.
    ///
    /// Wrapped payloads are not Spenso abstract indices until they are cooked, so the result
    /// intentionally has no `TensorExpression` interface.
    fn wrap_indices(self_: PyRef<'_, Self>, header: Symbol) -> PythonExpression {
        idenso::python::wrap_indices(&idenso_input(&self_), header)
    }

    /// Flatten nested representation-index payloads using index cooking by default.
    #[pyo3(signature = (settings = None))]
    fn cook_indices(
        self_: PyRef<'_, Self>,
        py: Python<'_>,
        settings: Option<&PyCookSettings>,
    ) -> PyResult<Py<Self>> {
        let result = idenso::python::cook_indices(&idenso_input(&self_), settings)?;
        from_idenso_atom(&self_, py, result.expr)
    }

    /// Encode one function call as an ordinary Symbolica expression.
    #[pyo3(signature = (settings = None))]
    fn cook_function(
        self_: PyRef<'_, Self>,
        settings: Option<&PyCookSettings>,
    ) -> PyResult<PythonExpression> {
        idenso::python::cook_function(&idenso_input(&self_), settings)
    }

    /// Wrap only contracted-index payloads with `header` and re-infer the interface.
    ///
    /// Raises `ValueError` when the expression cannot be parsed as a tensor network.
    fn wrap_dummies(self_: PyRef<'_, Self>, py: Python<'_>, header: Symbol) -> PyResult<Py<Self>> {
        let result = idenso::python::wrap_dummies(&idenso_input(&self_), header)?;
        from_idenso_atom(&self_, py, result.expr)
    }

    /// Return the ordered external, uncontracted indices as Symbolica expressions.
    ///
    /// Raises `ValueError` when the expression cannot be parsed as a tensor network.
    fn list_dangling(self_: PyRef<'_, Self>) -> PyResult<Vec<PythonExpression>> {
        idenso::python::list_dangling(&idenso_input(&self_))
    }

    /// Canonically order tensor factors and deterministically rename contracted indices.
    fn canonize(self_: PyRef<'_, Self>, py: Python<'_>) -> PyResult<Py<Self>> {
        let result = idenso::python::canonize(&idenso_input(&self_))?;
        from_idenso_atom(&self_, py, result.expr)
    }

    /// Replace nested tensor subexpressions by aliases and return the root plus alias mappings.
    fn alias_subtensors(
        self_: PyRef<'_, Self>,
        tensor_name: &str,
    ) -> (PythonExpression, Vec<(PythonExpression, PythonExpression)>) {
        idenso::python::alias_subtensors(&idenso_input(&self_), tensor_name)
    }

    /// Complex-conjugate this tensor expression and re-infer its interface.
    fn spenso_conjugate(self_: PyRef<'_, Self>, py: Python<'_>) -> PyResult<Py<Self>> {
        let result = idenso::python::spenso_conjugate(&idenso_input(&self_));
        from_idenso_atom(&self_, py, result.expr)
    }

    /// Complex-conjugate this expression and transpose slots in `representation`.
    fn conjugate_transpose(
        self_: PyRef<'_, Self>,
        py: Python<'_>,
        representation: &SpensoRepresentation,
    ) -> PyResult<Py<Self>> {
        let result = idenso::python::conjugate_transpose(
            &idenso_input(&self_),
            RegisteredRepresentation(representation.representation),
        );
        from_idenso_atom(&self_, py, result.expr)
    }

    /// Construct the physics-aware Dirac adjoint and re-infer the tensor interface.
    fn dirac_adjoint(self_: PyRef<'_, Self>, py: Python<'_>) -> PyResult<Py<Self>> {
        let result = idenso::python::dirac_adjoint(&idenso_input(&self_))?;
        from_idenso_atom(&self_, py, result.expr)
    }

    /// Encode this expression using Idenso's reversible cooking format.
    #[pyo3(signature = (settings = None))]
    fn cook(
        self_: PyRef<'_, Self>,
        settings: Option<&PyCookSettings>,
    ) -> PyResult<PythonExpression> {
        idenso::python::cook(&idenso_input(&self_), settings)
    }

    /// Restore a reversibly cooked expression and re-infer its tensor interface.
    #[pyo3(signature = (settings = None))]
    fn uncook(
        self_: PyRef<'_, Self>,
        py: Python<'_>,
        settings: Option<&PyCookSettings>,
    ) -> PyResult<Py<Self>> {
        let result = idenso::python::uncook(&idenso_input(&self_), settings);
        from_idenso_atom(&self_, py, result.expr)
    }

    /// Simplify tensor shorthands using the configured Schoonschip traversal.
    #[pyo3(signature = (settings = None))]
    fn schoonschip(
        self_: PyRef<'_, Self>,
        py: Python<'_>,
        settings: Option<&PySchoonschipSettings>,
    ) -> PyResult<Py<Self>> {
        let result = idenso::python::schoonschip(&idenso_input(&self_), settings);
        from_idenso_atom(&self_, py, result.expr)
    }

    /// Parse and contract this expression as a symbolic tensor network.
    #[pyo3(signature = (settings = None, *, expand_contracted_sums = false))]
    fn schoonschip_net(
        self_: PyRef<'_, Self>,
        py: Python<'_>,
        settings: Option<&PySchoonschipSettings>,
        expand_contracted_sums: bool,
    ) -> PyResult<Py<Self>> {
        let result = idenso::python::schoonschip_net(
            &idenso_input(&self_),
            settings,
            expand_contracted_sums,
        )?;
        from_idenso_atom(&self_, py, result.expr)
    }

    /// Convert contracted rank-one tensors into compact dot-product notation.
    fn to_dots(self_: PyRef<'_, Self>, py: Python<'_>) -> PyResult<Py<Self>> {
        let result = idenso::python::to_dots(&idenso_input(&self_));
        from_idenso_atom(&self_, py, result.expr)
    }

    /// Canonicalize compact dot-product shorthands without expanding them.
    fn normalize_dots(self_: PyRef<'_, Self>, py: Python<'_>) -> PyResult<Py<Self>> {
        let result = idenso::python::normalize_dots(&idenso_input(&self_));
        from_idenso_atom(&self_, py, result.expr)
    }

    /// Expand dot products into explicit metric and indexed-vector contractions.
    fn expand_dots(self_: PyRef<'_, Self>, py: Python<'_>) -> PyResult<Py<Self>> {
        let result = idenso::python::expand_dots(&idenso_input(&self_))?;
        from_idenso_atom(&self_, py, result.expr)
    }

    /// Replace metric shorthand such as `g(p(rep), q(rep))` by a compact dot product.
    fn metric_shorthand_to_dot(self_: PyRef<'_, Self>, py: Python<'_>) -> PyResult<Py<Self>> {
        let result = idenso::python::metric_shorthand_to_dot(&idenso_input(&self_));
        from_idenso_atom(&self_, py, result.expr)
    }

    /// Expand every Idenso tensor shorthand into explicit tensor syntax.
    fn undo_all(self_: PyRef<'_, Self>, py: Python<'_>) -> PyResult<Py<Self>> {
        let result = idenso::python::undo_all(&idenso_input(&self_))?;
        from_idenso_atom(&self_, py, result.expr)
    }

    /// Expand Schoonschip shorthands while leaving dots, chains, and traces compact.
    fn undo_schoonschip(self_: PyRef<'_, Self>, py: Python<'_>) -> PyResult<Py<Self>> {
        let result = idenso::python::undo_schoonschip(&idenso_input(&self_))?;
        from_idenso_atom(&self_, py, result.expr)
    }

    /// Expand dot-product shorthands while leaving other shorthands compact.
    fn undo_dots(self_: PyRef<'_, Self>, py: Python<'_>) -> PyResult<Py<Self>> {
        let result = idenso::python::undo_dots(&idenso_input(&self_))?;
        from_idenso_atom(&self_, py, result.expr)
    }

    /// Expand open-chain shorthands while leaving other shorthands compact.
    fn undo_chain(self_: PyRef<'_, Self>, py: Python<'_>) -> PyResult<Py<Self>> {
        let result = idenso::python::undo_chain(&idenso_input(&self_))?;
        from_idenso_atom(&self_, py, result.expr)
    }

    /// Expand trace shorthands while leaving other shorthands compact.
    fn undo_trace(self_: PyRef<'_, Self>, py: Python<'_>) -> PyResult<Py<Self>> {
        let result = idenso::python::undo_trace(&idenso_input(&self_))?;
        from_idenso_atom(&self_, py, result.expr)
    }

    /// Join adjacent open chains for `representation` and re-infer the interface.
    fn collect_chains(
        self_: PyRef<'_, Self>,
        py: Python<'_>,
        representation: &SpensoRepresentation,
    ) -> PyResult<Py<Self>> {
        let result = idenso::python::collect_chains(
            &idenso_input(&self_),
            RegisteredRepresentation(representation.representation),
        );
        from_idenso_atom(&self_, py, result.expr)
    }

    /// Rewrite tensors with two `representation` slots as open-chain factors.
    fn chainify(
        self_: PyRef<'_, Self>,
        py: Python<'_>,
        representation: &SpensoRepresentation,
    ) -> PyResult<Py<Self>> {
        let result = idenso::python::chainify(
            &idenso_input(&self_),
            RegisteredRepresentation(representation.representation),
        );
        from_idenso_atom(&self_, py, result.expr)
    }

    /// Convert chains whose endpoints coincide into trace shorthands.
    fn normalize_chains(self_: PyRef<'_, Self>, py: Python<'_>) -> PyResult<Py<Self>> {
        let result = idenso::python::normalize_chains(&idenso_input(&self_));
        from_idenso_atom(&self_, py, result.expr)
    }

    /// Replace one-factor chain shorthands by their underlying tensor factor.
    fn undo_single_length(self_: PyRef<'_, Self>, py: Python<'_>) -> PyResult<Py<Self>> {
        let result = idenso::python::undo_single_length(&idenso_input(&self_));
        from_idenso_atom(&self_, py, result.expr)
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
        let expression = Self::from_atom_interface_descriptor(
            py,
            self_.as_super().expr.clone(),
            self_.interface.clone(),
            self_.name,
            self_.name_args.clone(),
        )?;
        SpensoNet::from_arithmetic(ArithmeticStructure::Tensor(expression), library)
            .map_err(|error| PyRuntimeError::new_err(error.to_string()))
    }

    /// Fill the unresolved external ports with `indices` in interface order.
    ///
    /// Pass `AUTO` to leave a port unresolved. Set `cook_indices=True` to flatten nested
    /// symbolic index payloads before insertion. Repeated compatible indices contract their
    /// ports, in which case the result no longer carries the original stored-data identity.
    #[pyo3(signature = (*indices, cook_indices = false))]
    fn index(
        self_: PyRef<'_, Self>,
        py: Python<'_>,
        indices: &Bound<'_, PyTuple>,
        cook_indices: bool,
    ) -> PyResult<Py<Self>> {
        let replacements = Self::index_replacements(&self_.interface, indices, cook_indices)?;
        let mut logical = self_.interface.logical_slots();
        for (&position, &index) in &replacements {
            logical[position].set_aind(PartialIndex::Explicit(index));
        }
        let value = Self::structured(&self_);
        let atom = composition::materialize_interface_ports(&value, &replacements)
            .map_err(|error| PyValueError::new_err(error.to_string()))?;
        Self::from_atom_interface_descriptor(
            py,
            atom,
            PartialStructure::from_logical_slots(logical),
            self_.name,
            self_.name_args.clone(),
        )
    }

    /// Fill the unresolved external ports with `indices` in interface order.
    #[pyo3(signature = (*indices, cook_indices = false))]
    fn __call__(
        self_: PyRef<'_, Self>,
        py: Python<'_>,
        indices: &Bound<'_, PyTuple>,
        cook_indices: bool,
    ) -> PyResult<Py<Self>> {
        Self::index(self_, py, indices, cook_indices)
    }

    fn __neg__(self_: PyRef<'_, Self>, py: Python<'_>) -> PyResult<Py<Self>> {
        Self::from_atom_interface(
            py,
            Atom::num(-1) * self_.as_super().expr.as_ref(),
            self_.interface.clone(),
        )
    }

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
                if !additive_interfaces_match(&left.interface, &right.interface) {
                    return Err(PyValueError::new_err(
                        "addition requires compatible tensor interfaces",
                    ));
                }
                (right.atom, left.interface)
            }
            TensorOperand::Scalar(right) if right.as_view().is_zero() => {
                return Self::from_atom_interface_descriptor(
                    py,
                    left.atom,
                    left.interface,
                    self_.name,
                    self_.name_args.clone(),
                )
                .map(TensorDispatch::Expression);
            }
            TensorOperand::Scalar(right) if left.is_scalar() => (right, left.interface),
            TensorOperand::Scalar(_) => {
                return Err(PyValueError::new_err(
                    "cannot add a scalar expression to a non-scalar tensor",
                ));
            }
        };
        Self::from_atom_interface(py, left.atom.as_ref() + right.as_ref(), interface)
            .map(TensorDispatch::Expression)
    }

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
                if !additive_interfaces_match(&left.interface, &right.interface) {
                    return Err(PyValueError::new_err(
                        "subtraction requires compatible tensor interfaces",
                    ));
                }
                (right.atom, left.interface)
            }
            TensorOperand::Scalar(right) if right.as_view().is_zero() => {
                return Self::from_atom_interface_descriptor(
                    py,
                    left.atom,
                    left.interface,
                    self_.name,
                    self_.name_args.clone(),
                )
                .map(TensorDispatch::Expression);
            }
            TensorOperand::Scalar(right) if left.is_scalar() => (right, left.interface),
            TensorOperand::Scalar(_) => {
                return Err(PyValueError::new_err(
                    "cannot subtract a scalar expression from a non-scalar tensor",
                ));
            }
        };
        Self::from_atom_interface(py, left.atom.as_ref() - right.as_ref(), interface)
            .map(TensorDispatch::Expression)
    }

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
                if !additive_interfaces_match(&left.interface, &right.interface) {
                    return Err(PyValueError::new_err(
                        "subtraction requires compatible tensor interfaces",
                    ));
                }
                (left.atom, left.interface)
            }
            TensorOperand::Scalar(left) if left.as_view().is_zero() => {
                return Self::__neg__(self_, py).map(TensorDispatch::Expression);
            }
            TensorOperand::Scalar(left) if right.is_scalar() => (left, right.interface),
            TensorOperand::Scalar(_) => {
                return Err(PyValueError::new_err(
                    "cannot subtract a non-scalar tensor from a scalar expression",
                ));
            }
        };
        Self::from_atom_interface(py, left.as_ref() - right.atom.as_ref(), interface)
            .map(TensorDispatch::Expression)
    }

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
            TensorOperand::Structured(right) => composition::multiply(&left, &right)
                .map_err(|error| PyValueError::new_err(error.to_string()))
                .and_then(|value| Self::from_structured(py, value))
                .map(TensorDispatch::Expression),
            TensorOperand::Scalar(right) => {
                Self::from_atom_interface(py, left.atom.as_ref() * right.as_ref(), left.interface)
                    .map(TensorDispatch::Expression)
            }
        }
    }

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
            TensorOperand::Structured(left) => composition::multiply(&left, &right)
                .map_err(|error| PyValueError::new_err(error.to_string()))
                .and_then(|value| Self::from_structured(py, value))
                .map(TensorDispatch::Expression),
            TensorOperand::Scalar(left) => {
                Self::from_atom_interface(py, left.as_ref() * right.atom.as_ref(), right.interface)
                    .map(TensorDispatch::Expression)
            }
        }
    }

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
            TensorOperand::Structured(rhs) if rhs.is_scalar() => rhs.atom,
            TensorOperand::Structured(_) => {
                return Err(PyTypeError::new_err(
                    "tensor division requires a scalar denominator",
                ));
            }
        };
        Self::from_atom_interface(py, value.atom.as_ref() / rhs.as_ref(), value.interface)
            .map(TensorDispatch::Expression)
    }

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
        let lhs = match TensorOperand::extract(lhs)? {
            TensorOperand::Scalar(lhs) => lhs,
            TensorOperand::Structured(lhs) if lhs.is_scalar() => lhs.atom,
            TensorOperand::Structured(_) => {
                return Err(PyTypeError::new_err(
                    "tensor division requires a scalar numerator",
                ));
            }
        };
        Self::from_atom_interface(py, lhs.as_ref() / value.atom.as_ref(), value.interface)
            .map(TensorDispatch::Expression)
    }

    /// Form an outer product without contracting compatible ports.
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
        let left_slots = left.interface.logical_slots();
        let right_slots = right.interface.logical_slots();
        let collisions = composition::compatible_pairs(&left.interface, &right.interface)
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
        Self::from_structured(py, composition::outer(&left, &right)).map(TensorDispatch::Expression)
    }

    /// Contract one selected pair of ordered interface positions.
    #[pyo3(signature = (rhs, *, left, right))]
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
        composition::contract(&value, &rhs, composition::PortPair { left, right })
            .map_err(|error| PyValueError::new_err(error.to_string()))
            .and_then(|value| Self::from_structured(py, value))
            .map(TensorDispatch::Expression)
    }

    /// Compose two selected `(input, output)` matrix channels.
    #[pyo3(signature = (rhs, *, left, right))]
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
        composition::compose(
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
                composition::trace(&value, composition::MatrixChannel { input, output })
            }
            None => composition::trace_unique(&value),
        };
        traced
            .map_err(|error| PyValueError::new_err(error.to_string()))
            .and_then(|value| Self::from_structured(py, value))
    }

    /// Format this structured expression using compact Spenso notation.
    #[pyo3(signature = (show_dimensions = false))]
    fn format_tensor(self_: PyRef<'_, Self>, show_dimensions: bool) -> String {
        display::format_structured(&Self::structured(&self_), show_dimensions)
    }

    /// Format this structured expression as Typst math source.
    #[pyo3(signature = (show_dimensions = false))]
    fn to_typst(self_: PyRef<'_, Self>, show_dimensions: bool) -> String {
        display::structured_to_typst(&Self::structured(&self_), show_dimensions)
    }

    /// Build Symbolica's rich display wrapper for this structured expression.
    #[pyo3(signature = (show_dimensions = false))]
    fn formatted(self_: PyRef<'_, Self>, show_dimensions: bool) -> PythonFormattedOutput {
        display::format_structured_output(&Self::structured(&self_), show_dimensions)
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

    fn _repr_html_(self_: PyRef<'_, Self>) -> Option<String> {
        display::format_structured_output(&Self::structured(&self_), false).html
    }

    fn _repr_latex_(self_: PyRef<'_, Self>) -> String {
        display::structured_to_latex(&Self::structured(&self_), false)
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
        return TensorExpression::from_atom_interface_descriptor(
            py,
            expression.as_super().expr.clone(),
            expression.interface.clone(),
            expression.name,
            expression.name_args.clone(),
        );
    }
    let expression = expression
        .extract::<ConvertibleToExpression>()
        .map_err(|_| {
            PyTypeError::new_err("as_tensor() expects an Expression or TensorExpression")
        })?;
    let atom = expression.to_expression().expr;
    let (name, name_args) = inferred_descriptor(atom.as_view());
    let value = reinfer_structured(atom)?;
    TensorExpression::from_atom_interface_descriptor(
        py,
        value.atom,
        value.interface,
        name,
        name_args,
    )
}

fn structured_operand(value: &Bound<'_, PyAny>, operation: &str) -> PyResult<StructuredAtom> {
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
    pyo3_stub_gen::derive::gen_stub_pyfunction(module = "symbolica.community.spenso")
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
    if left.rank() != 1 || right.rank() != 1 {
        return Err(PyValueError::new_err(format!(
            "dot() requires rank-one operands, got ranks {} and {}",
            left.rank(),
            right.rank()
        )));
    }
    composition::contract(&left, &right, composition::PortPair { left: 0, right: 0 })
        .map_err(|error| PyValueError::new_err(error.to_string()))
        .and_then(|value| TensorExpression::from_structured(py, value))
        .map(TensorDispatch::Expression)
}

/// Build an explicitly-ended ordered tensor chain.
#[cfg_attr(
    feature = "python_stubgen",
    pyo3_stub_gen::derive::gen_stub_pyfunction(module = "symbolica.community.spenso")
)]
#[pyfunction]
#[pyo3(signature = (start_slot, end_slot, *factors))]
fn chain(
    py: Python<'_>,
    start_slot: SpensoSlot,
    end_slot: SpensoSlot,
    factors: &Bound<'_, PyTuple>,
) -> PyResult<TensorDispatch> {
    if factors.is_empty() {
        let start_representation = start_slot.slot.rep();
        let end_representation = end_slot.slot.rep();
        if !start_representation.matches(&end_representation) {
            return Err(PyValueError::new_err(
                "chain endpoints carry incompatible representations",
            ));
        }
        if !start_representation.rep.is_self_dual()
            && !(start_representation.rep.is_base() && end_representation.rep.is_dual())
        {
            return Err(PyValueError::new_err(
                "chain endpoints do not form an input-to-output propagation channel",
            ));
        }
        let closed = start_slot.slot.aind() == end_slot.slot.aind();
        let atom = if closed {
            SPENSO_TAG.trace(
                start_representation.to_symbolic([]),
                std::iter::empty::<Atom>(),
            )
        } else {
            SPENSO_TAG.chain(
                start_slot.slot.to_atom(),
                end_slot.slot.to_atom(),
                std::iter::empty::<Atom>(),
            )
        };
        let interface = if closed {
            PartialStructure::from_logical_slots([])
        } else {
            PartialStructure::from_logical_slots([
                start_slot
                    .slot
                    .rep()
                    .slot(PartialIndex::Explicit(start_slot.slot.aind())),
                end_slot
                    .slot
                    .rep()
                    .slot(PartialIndex::Explicit(end_slot.slot.aind())),
            ])
        };
        return TensorExpression::from_atom_interface(py, atom, interface)
            .map(TensorDispatch::Expression);
    }

    if factors.iter().any(|factor| is_concrete_network(&factor)) {
        let mut factors = factors.iter();
        let first = factors
            .next()
            .expect("the empty factor sequence was handled above")
            .extract::<ConvertibleToSpensoNet>()?
            .to_net();
        let first_channel = composition::matrix_channel(&first.structure)
            .ok_or_else(|| PyValueError::new_err("chain factors have no unique matrix channel"))?;
        let mut value = first.chain_form(first_channel)?;
        for factor in factors {
            let factor = factor.extract::<ConvertibleToSpensoNet>()?.to_net();
            let factor_channel =
                composition::matrix_channel(&factor.structure).ok_or_else(|| {
                    PyValueError::new_err("chain factors have no unique matrix channel")
                })?;
            value = value.compose(
                ConvertibleToSpensoNet(factor),
                (0, 1),
                (factor_channel.input, factor_channel.output),
            )?;
        }
        let slots = value.structure.interface.logical_slots();
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

    let mut factors = factors.iter();
    let first = factors
        .next()
        .expect("the empty factor sequence was handled above");
    let first = structured_operand(&first, "chain()")?;
    let first_channel = composition::matrix_channel(&first)
        .ok_or_else(|| PyValueError::new_err("chain factors have no unique matrix channel"))?;
    let first_slots = first.interface.logical_slots();
    let input = first_slots[first_channel.input];
    let output = first_slots[first_channel.output];
    let atom = if first.atom.as_view().is_zero() {
        Atom::Zero
    } else {
        SPENSO_TAG.chain(
            composition::port_atom(input),
            composition::port_atom(output),
            composition::chain_factors(&first, first_channel)
                .map_err(|error| PyValueError::new_err(error.to_string()))?,
        )
    };
    let interface = PartialStructure::from_logical_slots(
        [input, output].into_iter().chain(
            first_slots
                .into_iter()
                .enumerate()
                .filter(|(position, _)| {
                    *position != first_channel.input && *position != first_channel.output
                })
                .map(|(_, slot)| slot),
        ),
    );
    let mut value = StructuredAtom::new(atom, interface);
    for factor in factors {
        let factor = structured_operand(&factor, "chain()")?;
        let factor_channel = composition::matrix_channel(&factor)
            .ok_or_else(|| PyValueError::new_err("chain factors have no unique matrix channel"))?;
        value = composition::compose(
            &value,
            &factor,
            composition::MatrixChannel {
                input: 0,
                output: 1,
            },
            factor_channel,
        )
        .map_err(|error| PyValueError::new_err(error.to_string()))?;
    }
    let channel = composition::MatrixChannel {
        input: 0,
        output: 1,
    };
    let slots = value.interface.logical_slots();
    if start_slot.slot.rep() != slots[channel.input].rep()
        || end_slot.slot.rep() != slots[channel.output].rep()
    {
        return Err(PyValueError::new_err(
            "chain endpoints are incompatible with the factor channel",
        ));
    }
    if start_slot.slot.aind() == end_slot.slot.aind() {
        return composition::trace(&value, channel)
            .map_err(|error| PyValueError::new_err(error.to_string()))
            .and_then(|value| TensorExpression::from_structured(py, value))
            .map(TensorDispatch::Expression);
    }
    let replacements = HashMap::from([
        (channel.input, start_slot.slot.aind()),
        (channel.output, end_slot.slot.aind()),
    ]);
    let atom = composition::materialize_interface_ports(&value, &replacements)
        .map_err(|error| PyValueError::new_err(error.to_string()))?;
    let mut interface = slots;
    interface[channel.input] = start_slot
        .slot
        .rep()
        .slot(PartialIndex::Explicit(start_slot.slot.aind()));
    interface[channel.output] = end_slot
        .slot
        .rep()
        .slot(PartialIndex::Explicit(end_slot.slot.aind()));
    TensorExpression::from_atom_interface(py, atom, PartialStructure::from_logical_slots(interface))
        .map(TensorDispatch::Expression)
}

/// Close an ordered factor sequence into a canonical cyclic trace.
#[cfg_attr(
    feature = "python_stubgen",
    pyo3_stub_gen::derive::gen_stub_pyfunction(module = "symbolica.community.spenso")
)]
#[pyfunction]
#[pyo3(signature = (representation, *factors))]
fn trace(
    py: Python<'_>,
    representation: SpensoRepresentation,
    factors: &Bound<'_, PyTuple>,
) -> PyResult<TensorDispatch> {
    if factors.is_empty() {
        let atom = SPENSO_TAG.trace(
            representation.representation.to_symbolic([]),
            std::iter::empty::<Atom>(),
        );
        return TensorExpression::from_atom_interface(
            py,
            atom,
            PartialStructure::from_logical_slots(std::iter::empty()),
        )
        .map(TensorDispatch::Expression);
    }

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
        let channel = composition::matrix_channel(&value.structure)
            .ok_or_else(|| PyValueError::new_err("trace factors have no unique matrix channel"))?;
        let input = value.structure.interface.logical_slots()[channel.input];
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

    let mut factors = factors.iter();
    let first = factors
        .next()
        .expect("the empty factor sequence was handled above");
    let mut value = structured_operand(&first, "trace()")?;
    for factor in factors {
        let factor = structured_operand(&factor, "trace()")?;
        value = composition::multiply(&value, &factor)
            .map_err(|error| PyValueError::new_err(error.to_string()))?;
    }
    let channel = composition::matrix_channel(&value)
        .ok_or_else(|| PyValueError::new_err("trace factors have no unique matrix channel"))?;
    let input = value.interface.logical_slots()[channel.input];
    if representation.representation != input.rep() {
        return Err(PyValueError::new_err(
            "trace representation does not match the factor channel",
        ));
    }
    composition::trace(&value, channel)
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
                type_ignored: None,
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
                type_ignored: None,
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
                type_ignored: None,
                is_overload: true,
            },
        ],
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use idenso::{IndexTooling, IndexToolingError};
    use spenso::structure::{dimension::Dimension, representation::ExtendibleReps};
    use symbolica::{atom::FunctionBuilder, symbol};

    fn explicit_interface(
        representation: Representation<LibraryRep>,
        index: AbstractIndex,
    ) -> PartialStructure {
        PartialStructure::from_logical_slots([representation.slot(PartialIndex::Explicit(index))])
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
                "interface",
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
                    vec![minkowski, bispinor, bispinor],
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
                (TensorExpression::g(py, &euc_object)?, vec![dimension; 2]),
                (TensorExpression::flat(py, &euc_object)?, vec![dimension; 2]),
                (
                    TensorExpression::gamma(py, ConvertibleToDimension(dimension))?,
                    vec![dimension, concrete_four, concrete_four],
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
        let generator = value_to_structured_atom(&generator_structure).unwrap();
        let inferred = infer_interface(&generator.atom).unwrap();
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
                value.interface.logical_slots(),
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
            let expression = TensorExpression::from_atom_interface_descriptor(
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
    fn wrap_indices_returns_an_ordinary_expression() {
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
                StructuredAtom::new(atom.clone(), interface),
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
            let wrapped = TensorExpression::wrap_indices(expression.bind(py).borrow(), header);
            assert_eq!(wrapped.expr, atom.wrap_indices(header));
            let wrapped = wrapped.into_pyobject(py)?;
            assert!(!wrapped.is_instance_of::<TensorExpression>());
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
            let expression = TensorExpression::from_atom_interface_descriptor(
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

            let ordinary = TensorExpression::to_expression(expression.bind(py).borrow());
            let ordinary = ordinary.expand(None, None)?;
            let ordinary = ordinary.into_pyobject(py)?;
            assert!(!ordinary.is_instance_of::<TensorExpression>());

            let zero = TensorExpression::from_atom_interface(py, Atom::Zero, interface.clone())?;
            let network = TensorExpression::to_network(zero.bind(py).borrow(), py, None)?;
            assert_eq!(
                network.structure.interface.logical_slots(),
                interface.logical_slots()
            );
            assert!(network.materialized.interface.open_positions().is_empty());
            assert!(network.network.state.is_tensor());
            Ok(())
        })
        .unwrap();
    }

    #[test]
    fn recursive_merging_rejects_three_compatible_explicit_ports() {
        let representation = ExtendibleReps::EUCLIDEAN.new_rep(Dimension::Concrete(4));
        let index = AbstractIndex::Normal(11);
        let interfaces = vec![
            explicit_interface(representation, index),
            explicit_interface(representation, index),
            explicit_interface(representation, index),
        ];

        assert!(merge_explicit_interface_sequence(&interfaces).is_err());
    }

    #[test]
    fn recursive_merging_contracts_one_explicit_pair() {
        let representation = ExtendibleReps::EUCLIDEAN.new_rep(Dimension::Concrete(4));
        let index = AbstractIndex::Normal(13);
        let interfaces = vec![
            explicit_interface(representation, index),
            explicit_interface(representation, index),
        ];

        assert!(
            merge_explicit_interface_sequence(&interfaces)
                .unwrap()
                .canonical()
                .is_scalar()
        );
    }

    #[test]
    fn recursive_merging_normalizes_explicit_indices_within_one_interface() {
        let representation = ExtendibleReps::EUCLIDEAN.new_rep(Dimension::Concrete(4));
        let index = AbstractIndex::Normal(15);
        let pair = PartialStructure::from_logical_slots([
            representation.slot(PartialIndex::Explicit(index)),
            representation.slot(PartialIndex::Explicit(index)),
        ]);
        let triple = PartialStructure::from_logical_slots([
            representation.slot(PartialIndex::Explicit(index)),
            representation.slot(PartialIndex::Explicit(index)),
            representation.slot(PartialIndex::Explicit(index)),
        ]);

        assert!(
            merge_explicit_interface_sequence(&[pair])
                .unwrap()
                .canonical()
                .is_scalar()
        );
        assert!(merge_explicit_interface_sequence(&[triple]).is_err());
    }

    #[test]
    fn additive_reinference_compares_canonical_scalar_interfaces() {
        let representation = ExtendibleReps::EUCLIDEAN.new_rep(Dimension::Concrete(4));
        let diagonal = |name: &str, index| {
            let slot = representation.slot::<AbstractIndex, _>(index).to_atom();
            FunctionBuilder::new(SPENSO_TAG.tensor_symbol(name))
                .add_arg(&slot)
                .add_arg(&slot)
                .finish()
        };
        let left = diagonal("additive_scalar_left", AbstractIndex::Normal(21));
        let right = diagonal("additive_scalar_right", AbstractIndex::Normal(23));

        assert!(
            infer_validated_interface(&(left.as_ref() + Atom::num(1)))
                .unwrap()
                .canonical()
                .is_scalar()
        );
        assert!(
            infer_validated_interface(&(left + right))
                .unwrap()
                .canonical()
                .is_scalar()
        );
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
            StructuredAtom::new(atom, interface)
        };
        let left = tensor("canonical_addition_left", [i, j]);
        let right = tensor("canonical_addition_right", [j, i]);
        let expected_sum = left.atom.as_ref() + right.atom.as_ref();
        let expected_difference = left.atom.as_ref() - right.atom.as_ref();

        assert_ne!(
            left.interface.logical_slots(),
            right.interface.logical_slots()
        );
        assert!(additive_interfaces_match(&left.interface, &right.interface));
        let mixed_left = PartialStructure::from_logical_slots([
            representation.slot(PartialIndex::Explicit(i)),
            representation.slot(PartialIndex::open(0)),
        ]);
        let mixed_right = PartialStructure::from_logical_slots([
            representation.slot(PartialIndex::open(0)),
            representation.slot(PartialIndex::Explicit(i)),
        ]);
        assert_eq!(mixed_left.canonical(), mixed_right.canonical());
        assert!(!additive_interfaces_match(&mixed_left, &mixed_right));
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
            infer_validated_interface(&open_sum)
                .unwrap()
                .logical_slots(),
            open_interface.logical_slots()
        );
        let inferred = infer_validated_interface(&expected_sum).unwrap();
        assert!(additive_interfaces_match(&left.interface, &inferred));

        let standalone = StructuredAtom::new(left.atom.clone(), right.interface.clone());
        assert_eq!(
            standalone.presentation_atom(),
            tensor("canonical_addition_left", [j, i]).atom
        );

        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let expected = left.interface.logical_slots();
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
            let expression = TensorExpression::from_atom_interface_descriptor(
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
    fn powers_reject_more_than_two_explicit_index_occurrences() {
        let representation = ExtendibleReps::EUCLIDEAN.new_rep(Dimension::Concrete(4));
        let index = AbstractIndex::Normal(17);
        let tensor = FunctionBuilder::new(SPENSO_TAG.tensor_symbol("powered_explicit_index"))
            .add_arg(representation.slot::<AbstractIndex, _>(index).to_atom())
            .finish();

        assert!(
            infer_validated_interface(&tensor.as_ref().pow(Atom::num(2)))
                .unwrap()
                .canonical()
                .is_scalar()
        );
        assert!(infer_validated_interface(&tensor.as_ref().pow(Atom::num(3))).is_err());
    }

    #[test]
    fn powers_count_repeated_explicit_indices_in_the_base() {
        let representation = ExtendibleReps::EUCLIDEAN.new_rep(Dimension::Concrete(4));
        let index = AbstractIndex::Normal(19);
        let slot = representation.slot::<AbstractIndex, _>(index).to_atom();
        let tensor = FunctionBuilder::new(SPENSO_TAG.tensor_symbol("powered_repeated_index"))
            .add_arg(&slot)
            .add_arg(&slot)
            .finish();

        assert!(infer_validated_interface(&tensor.as_ref().pow(Atom::num(2))).is_err());
    }

    #[test]
    fn open_ports_retain_power_parity_semantics() {
        let representation = ExtendibleReps::EUCLIDEAN.new_rep(Dimension::Concrete(4));
        let tensor = FunctionBuilder::new(SPENSO_TAG.tensor_symbol("powered_open_index"))
            .add_arg(representation.to_symbolic([]))
            .finish();

        assert!(
            infer_validated_interface(&tensor.as_ref().pow(Atom::num(2)))
                .unwrap()
                .canonical()
                .is_scalar()
        );
        assert_eq!(
            infer_validated_interface(&tensor.as_ref().pow(Atom::num(3)))
                .unwrap()
                .canonical()
                .order(),
            1
        );
    }

    #[test]
    fn leaf_slots_follow_syntactic_argument_order() {
        let euclidean = ExtendibleReps::EUCLIDEAN.new_rep(Dimension::Concrete(4));
        let minkowski = ExtendibleReps::MINKOWSKI.new_rep(Dimension::Concrete(4));
        let atom = FunctionBuilder::new(SPENSO_TAG.bracket)
            .add_arg(
                euclidean
                    .slot::<AbstractIndex, _>(AbstractIndex::Normal(17))
                    .to_atom(),
            )
            .add_arg(
                minkowski
                    .slot::<AbstractIndex, _>(AbstractIndex::Normal(19))
                    .to_atom(),
            )
            .finish();

        assert_eq!(
            syntactic_leaf_slots(atom.as_view())
                .into_iter()
                .map(|slot| slot.rep())
                .collect::<Vec<_>>(),
            vec![euclidean, minkowski]
        );
    }

    #[test]
    fn chain_placeholders_require_a_factor_scope() {
        Python::initialize();
        let placeholder = Atom::var(SPENSO_TAG.chain_in);
        assert!(validate_placeholder_scope(placeholder.as_view(), false).is_err());

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

        assert!(validate_placeholder_scope(chain.as_view(), false).is_ok());

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
            let metric = TensorExpression::g(py, &SpensoRepresentation { representation })?;
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
                StructuredAtom::new(atom, PartialStructure::from_logical_slots(ports))
            };
            let composed = composition::multiply(
                &matrix("index_tooling_left"),
                &matrix("index_tooling_right"),
            )
            .expect("compatible unresolved matrices should compose");
            assert!(matches!(
                composed.atom.as_view(),
                AtomView::Fun(function) if function.get_symbol() == SPENSO_TAG.chain
            ));
            let rust_list_error = composed
                .atom
                .list_dangling::<AbstractIndex>()
                .expect_err("unresolved chain endpoints should return a parse error");
            let rust_wrap_error = composed
                .atom
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
                py,
                symbolica::symbol!("index_tooling_wrap"),
            )
            .expect_err("unresolved chain endpoints should return a parse error");

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
            let expression = TensorExpression::from_atom_interface(
                py,
                malformed,
                PartialStructure::from_logical_slots(std::iter::empty()),
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
}
