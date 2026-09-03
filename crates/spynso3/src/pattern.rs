//! Python builders for Spenso's tagged Symbolica pattern syntax.

#[cfg(not(feature = "python_stubgen"))]
use pyo3_stub_gen_derive::remove_gen_stub;

use idenso::{
    color::CS,
    dirac::AGS,
    representations::{Bispinor, ColorAdjoint, ColorAntiFundamental, ColorFundamental},
};
use pyo3::{
    exceptions::{PyTypeError, PyValueError},
    prelude::*,
};
#[cfg(feature = "python_stubgen")]
use pyo3_stub_gen::{PyStubType, TypeInfo, derive::*};
use spenso::{
    network::{
        library::symbolic::{ETS, ExplicitKey},
        tags::{SPENSO_TAG, SymbolAtomExt},
    },
    structure::{
        Canonicalized,
        abstract_index::{AIND_SYMBOLS, AbstractIndex},
        representation::{Minkowski, RepName},
        slot::IsAbstractSlot,
    },
};
use symbolica::{
    api::python::{ConvertibleToExpression, PythonExpression},
    atom::{Atom, AtomView, FunctionBuilder, Symbol},
};

use crate::{
    ModuleInit,
    structure::{SpensoName, SpensoRepresentation, SpensoSlot},
};

/// An expression accepted in a tensor pattern's scalar-argument section.
struct PatternArgument(Atom);

impl<'a, 'py> FromPyObject<'a, 'py> for PatternArgument {
    type Error = PyErr;

    fn extract(value: pyo3::Borrowed<'a, 'py, PyAny>) -> PyResult<Self> {
        value
            .extract::<ConvertibleToExpression>()
            .map(|value| Self(value.to_expression().expr))
    }
}

impl<'py> IntoPyObject<'py> for PatternArgument {
    type Target = PythonExpression;
    type Output = Bound<'py, PythonExpression>;
    type Error = PyErr;

    fn into_pyobject(self, py: Python<'py>) -> Result<Self::Output, Self::Error> {
        PythonExpression::from(self.0).into_pyobject(py)
    }
}

#[cfg(feature = "python_stubgen")]
impl PyStubType for PatternArgument {
    fn type_input() -> TypeInfo {
        ConvertibleToExpression::type_input()
    }

    fn type_output() -> TypeInfo {
        PythonExpression::type_output()
    }
}

/// An expression accepted in a tensor pattern's structural-port section.
struct PatternPort(Atom);

impl<'a, 'py> FromPyObject<'a, 'py> for PatternPort {
    type Error = PyErr;

    fn extract(value: pyo3::Borrowed<'a, 'py, PyAny>) -> PyResult<Self> {
        if let Ok(slot) = value.extract::<SpensoSlot>() {
            Ok(Self(slot.slot.to_atom()))
        } else if let Ok(representation) = value.extract::<SpensoRepresentation>() {
            Ok(Self(representation.representation.to_symbolic([])))
        } else if let Ok(expression) = value.extract::<ConvertibleToExpression>() {
            Ok(Self(expression.to_expression().expr))
        } else {
            Err(PyTypeError::new_err(
                "pattern ports must be Slots, Representations, or Symbolica expressions",
            ))
        }
    }
}

impl<'py> IntoPyObject<'py> for PatternPort {
    type Target = PythonExpression;
    type Output = Bound<'py, PythonExpression>;
    type Error = PyErr;

    fn into_pyobject(self, py: Python<'py>) -> Result<Self::Output, Self::Error> {
        PythonExpression::from(self.0).into_pyobject(py)
    }
}

#[cfg(feature = "python_stubgen")]
impl PyStubType for PatternPort {
    fn type_input() -> TypeInfo {
        SpensoSlot::type_input()
            | SpensoRepresentation::type_input()
            | ConvertibleToExpression::type_input()
    }

    fn type_output() -> TypeInfo {
        Self::type_input()
    }
}

/// A concrete or symbolic stripped representation used by metric patterns.
struct RepresentationPattern(Atom);

impl<'a, 'py> FromPyObject<'a, 'py> for RepresentationPattern {
    type Error = PyErr;

    fn extract(value: pyo3::Borrowed<'a, 'py, PyAny>) -> PyResult<Self> {
        if let Ok(representation) = value.extract::<SpensoRepresentation>() {
            Ok(Self(representation.representation.to_symbolic([])))
        } else if let Ok(expression) = value.extract::<PythonExpression>() {
            Ok(Self(expression.expr))
        } else {
            Err(PyTypeError::new_err(
                "representation pattern must be a Representation or a PortPattern",
            ))
        }
    }
}

#[cfg(feature = "python_stubgen")]
impl PyStubType for RepresentationPattern {
    fn type_input() -> TypeInfo {
        SpensoRepresentation::type_input() | PythonExpression::type_input()
    }

    fn type_output() -> TypeInfo {
        Self::type_input()
    }
}

fn expression_atoms(values: Vec<PatternArgument>) -> Vec<Atom> {
    values.into_iter().map(|value| value.0).collect()
}

fn port_atoms(values: Vec<PatternPort>) -> Vec<Atom> {
    values.into_iter().map(|value| value.0).collect()
}

fn fixed_tensor(head: Symbol, args: Vec<Atom>, ports: Vec<Atom>) -> Atom {
    FunctionBuilder::new(head)
        .add_args(&args)
        .add_args(&ports)
        .finish()
}

fn wildcard_head<E: std::fmt::Display>(
    name: &str,
    make_symbol: impl FnOnce(&str) -> Result<Symbol, E>,
) -> PyResult<Symbol> {
    let Some(base) = name.strip_suffix('_') else {
        return Err(PyValueError::new_err(format!(
            "wildcard tensor and representation heads must end in exactly one underscore, got {name:?}",
        )));
    };
    if base.ends_with('_') {
        return Err(PyValueError::new_err(format!(
            "wildcard tensor and representation heads must end in exactly one underscore, got {name:?}",
        )));
    }
    let symbol = make_symbol(name).map_err(|error| {
        PyValueError::new_err(format!("invalid wildcard head {name:?}: {error}"))
    })?;
    let wildcard = Atom::var(symbol);
    let AtomView::Var(variable) = wildcard.as_view() else {
        unreachable!("a symbol atom is always a variable")
    };
    if variable.get_wildcard_level() != 1 {
        return Err(PyValueError::new_err(format!(
            "wildcard tensor and representation heads must end in exactly one underscore, got {name:?}",
        )));
    }
    Ok(symbol)
}

fn append_index(representation: AtomView<'_>, index: &Atom) -> PyResult<Atom> {
    let AtomView::Fun(function) = representation else {
        return Err(PyValueError::new_err(
            "expected a stripped representation pattern",
        ));
    };
    if !function.get_symbol().has_tag(&SPENSO_TAG.representation) || function.get_nargs() != 1 {
        return Err(PyValueError::new_err(
            "expected a stripped representation pattern with one dimension argument",
        ));
    }
    Ok(FunctionBuilder::new(function.get_symbol())
        .add_args(function.iter())
        .add_arg(index)
        .finish())
}

fn index_representation(representation: Atom, index: Atom) -> PyResult<Atom> {
    let AtomView::Fun(function) = representation.as_view() else {
        return Err(PyValueError::new_err(
            "expected a stripped representation pattern",
        ));
    };
    if function.get_symbol() == AIND_SYMBOLS.dind {
        let arguments = function.iter().collect::<Vec<_>>();
        let [inner] = arguments.as_slice() else {
            return Err(PyValueError::new_err(
                "malformed dual representation pattern",
            ));
        };
        return append_index(*inner, &index).map(|inner| AIND_SYMBOLS.dual(inner));
    }
    append_index(representation.as_view(), &index)
}

fn built_in_pattern(
    structure: &Canonicalized<ExplicitKey<AbstractIndex>>,
    logical_ports: Vec<Atom>,
) -> Atom {
    let ports = structure.layout().logical_to_canonical(&logical_ports);
    fixed_tensor(
        structure
            .canonical()
            .global_name
            .expect("built-in tensor structures have names"),
        Vec::new(),
        ports,
    )
}

fn minkowski_port(dimension: &Atom, index: &Atom) -> Atom {
    Minkowski {}.to_symbolic([dimension, index])
}

fn bispinor_port(dimension: &Atom, index: &Atom) -> Atom {
    Bispinor {}.to_symbolic([dimension, index])
}

fn adjoint_port(dimension: &Atom, index: &Atom) -> Atom {
    ColorAdjoint {}.to_symbolic([dimension, index])
}

fn fundamental_port(dimension: &Atom, index: &Atom) -> Atom {
    ColorFundamental {}.to_symbolic([dimension, index])
}

fn antifundamental_port(dimension: &Atom, index: &Atom) -> Atom {
    ColorAntiFundamental {}.to_symbolic([dimension, index])
}

/// A Symbolica expression representing one tensor port in a rewrite pattern.
///
/// Port patterns may name an exact representation or constrain a wildcard
/// representation head by its Spenso duality tags. Omitting `index` produces a
/// stripped representation, suitable for compact vector and trace patterns.
///
/// Examples
/// --------
/// >>> import symbolica as sp
/// >>> from symbolica.community.spenso import PortPattern, Representation
/// >>> D_, i_ = sp.S("D_", "i_")
/// >>> exact = PortPattern.exact(Representation.mink(4), i_)
/// >>> generic = PortPattern.self_dual("R_", D_, i_)
/// >>> stripped = PortPattern.dualizable("C_", D_)
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    frozen,
    extends = PythonExpression,
    from_py_object,
    name = "PortPattern",
    module = "symbolica.community.spenso"
)]
#[derive(Clone, Copy)]
pub struct PortPattern;

impl PortPattern {
    fn from_atom(py: Python<'_>, atom: Atom) -> PyResult<Py<Self>> {
        Py::new(py, (Self, PythonExpression { expr: atom }))
    }

    fn tagged<E: std::fmt::Display>(
        py: Python<'_>,
        name: String,
        dimension: ConvertibleToExpression,
        index: Option<ConvertibleToExpression>,
        symbol: impl FnOnce(&str) -> Result<Symbol, E>,
        dual: bool,
    ) -> PyResult<Py<Self>> {
        let head = wildcard_head(&name, symbol)?;
        let mut arguments = vec![dimension.to_expression().expr];
        arguments.extend(index.map(|value| value.to_expression().expr));
        let atom = head.atom_with_args(arguments);
        Self::from_atom(py, if dual { AIND_SYMBOLS.dual(atom) } else { atom })
    }
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "python_stubgen"), remove_gen_stub)]
#[pymethods]
impl PortPattern {
    /// Match one exact representation, optionally carrying `index`.
    #[staticmethod]
    #[pyo3(signature = (rep, index=None))]
    fn exact(
        py: Python<'_>,
        rep: &SpensoRepresentation,
        index: Option<ConvertibleToExpression>,
    ) -> PyResult<Py<Self>> {
        let index = index.map(|value| value.to_expression().expr);
        Self::from_atom(py, rep.representation.to_symbolic(index))
    }

    /// Match any Spenso representation head.
    ///
    /// `name` is the reusable Symbolica wildcard name and must end in one
    /// underscore, for example `"R_"`.
    #[staticmethod]
    #[pyo3(signature = (name, dimension, index=None))]
    fn any(
        py: Python<'_>,
        name: String,
        dimension: ConvertibleToExpression,
        index: Option<ConvertibleToExpression>,
    ) -> PyResult<Py<Self>> {
        Self::tagged(
            py,
            name,
            dimension,
            index,
            |name| symbolica::try_symbol!(name, tag = &SPENSO_TAG.representation),
            false,
        )
    }

    /// Match a self-dual representation head.
    #[staticmethod]
    #[pyo3(signature = (name, dimension, index=None))]
    fn self_dual(
        py: Python<'_>,
        name: String,
        dimension: ConvertibleToExpression,
        index: Option<ConvertibleToExpression>,
    ) -> PyResult<Py<Self>> {
        Self::tagged(
            py,
            name,
            dimension,
            index,
            |name| {
                symbolica::try_symbol!(
                    name,
                    tags = [&SPENSO_TAG.representation, &SPENSO_TAG.self_dual]
                )
            },
            false,
        )
    }

    /// Match a dualizable representation in its base or dual orientation.
    #[staticmethod]
    #[pyo3(signature = (name, dimension, index=None, *, dual=false))]
    fn dualizable(
        py: Python<'_>,
        name: String,
        dimension: ConvertibleToExpression,
        index: Option<ConvertibleToExpression>,
        dual: bool,
    ) -> PyResult<Py<Self>> {
        Self::tagged(
            py,
            name,
            dimension,
            index,
            |name| {
                symbolica::try_symbol!(
                    name,
                    tags = [&SPENSO_TAG.representation, &SPENSO_TAG.dualizable]
                )
            },
            dual,
        )
    }
}

impl ModuleInit for PortPattern {}

/// A tagged Symbolica expression intended for tensor rewrite rules.
///
/// `args` contains scalar tensor arguments and `ports` contains structural
/// syntax. The two named sections are concatenated in that order. Since this is
/// a pattern expression rather than a concrete tensor, either section may
/// contain ordinary or sequence-wildcard Symbolica expressions.
///
/// Examples
/// --------
/// >>> import symbolica as sp
/// >>> from symbolica.community.spenso import (
/// ...     PortPattern, TensorExpression, TensorName, TensorPattern,
/// ... )
/// >>> k_, D_, mu_, i_, j_, rest___ = sp.S(
/// ...     "k_", "D_", "mu_", "i_", "j_", "rest___",
/// ... )
/// >>> fixed = TensorPattern(
/// ...     TensorName("A"),
/// ...     args=[k_],
/// ...     ports=[PortPattern.any("R_", D_, i_), rest___],
/// ... )
/// >>> generic = TensorPattern.any("T_", ports=[rest___])
/// >>> target = TensorExpression.gamma(4)("mu", "i", "j").to_expression()
/// >>> rule = TensorPattern.gamma(D_, mu_, i_, j_)
/// >>> target.replace(rule, 0)
/// 0
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    frozen,
    extends = PythonExpression,
    from_py_object,
    name = "TensorPattern",
    module = "symbolica.community.spenso"
)]
#[derive(Clone, Copy)]
pub struct TensorPattern;

impl TensorPattern {
    fn from_atom(py: Python<'_>, atom: Atom) -> PyResult<Py<Self>> {
        Py::new(py, (Self, PythonExpression { expr: atom }))
    }

    fn wildcard<E: std::fmt::Display>(
        py: Python<'_>,
        name: String,
        args: Vec<PatternArgument>,
        ports: Vec<PatternPort>,
        symbol: impl FnOnce(&str) -> Result<Symbol, E>,
    ) -> PyResult<Py<Self>> {
        let head = wildcard_head(&name, symbol)?;
        let arguments = expression_atoms(args)
            .into_iter()
            .chain(port_atoms(ports))
            .collect::<Vec<_>>();
        Self::from_atom(py, head.atom_with_args(arguments))
    }

    fn built_in(
        py: Python<'_>,
        structure: &Canonicalized<ExplicitKey<AbstractIndex>>,
        logical_ports: Vec<Atom>,
    ) -> PyResult<Py<Self>> {
        Self::from_atom(py, built_in_pattern(structure, logical_ports))
    }
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "python_stubgen"), remove_gen_stub)]
#[pymethods]
impl TensorPattern {
    /// Create a pattern for the fixed tensor `head`.
    ///
    /// Scalar `args` are always emitted before structural `ports`.
    #[new]
    #[gen_stub(override_return_type(type_repr = "TensorPattern"))]
    #[pyo3(signature = (head, *, args=Vec::new(), ports=Vec::new()))]
    fn new(
        head: SpensoName,
        args: Vec<PatternArgument>,
        ports: Vec<PatternPort>,
    ) -> (TensorPattern, PythonExpression) {
        (
            TensorPattern,
            PythonExpression {
                expr: fixed_tensor(head.name, expression_atoms(args), port_atoms(ports)),
            },
        )
    }

    /// Match any tensor-tagged head.
    ///
    /// `name` must end in one underscore, for example `"T_"`.
    #[staticmethod]
    #[pyo3(signature = (name, *, args=Vec::new(), ports=Vec::new()))]
    fn any(
        py: Python<'_>,
        name: String,
        args: Vec<PatternArgument>,
        ports: Vec<PatternPort>,
    ) -> PyResult<Py<Self>> {
        Self::wildcard(py, name, args, ports, |name| {
            symbolica::try_symbol!(name, tag = &SPENSO_TAG.tensor)
        })
    }

    /// Match any rank-one tensor head.
    ///
    /// Pattern construction intentionally does not impose a concrete arity, so
    /// a sequence wildcard may describe the scalar arguments or structural port.
    #[staticmethod]
    #[pyo3(signature = (name, *, args=Vec::new(), ports=Vec::new()))]
    fn vector(
        py: Python<'_>,
        name: String,
        args: Vec<PatternArgument>,
        ports: Vec<PatternPort>,
    ) -> PyResult<Py<Self>> {
        Self::wildcard(py, name, args, ports, |name| {
            symbolica::try_symbol!(name, tags = [&SPENSO_TAG.tensor, &SPENSO_TAG.rank1])
        })
    }

    /// Match a metric tensor in logical index order.
    #[staticmethod]
    fn g(
        py: Python<'_>,
        rep_pattern: RepresentationPattern,
        i: ConvertibleToExpression,
        j: ConvertibleToExpression,
    ) -> PyResult<Py<Self>> {
        let i = index_representation(rep_pattern.0.clone(), i.to_expression().expr)?;
        let j = index_representation(rep_pattern.0, j.to_expression().expr)?;
        Self::from_atom(py, fixed_tensor(ETS.metric, Vec::new(), vec![i, j]))
    }

    /// Match a musical-isomorphism tensor in logical index order.
    #[staticmethod]
    fn flat(
        py: Python<'_>,
        rep_pattern: RepresentationPattern,
        i: ConvertibleToExpression,
        j: ConvertibleToExpression,
    ) -> PyResult<Py<Self>> {
        let i = index_representation(rep_pattern.0.clone(), i.to_expression().expr)?;
        let j = index_representation(rep_pattern.0, j.to_expression().expr)?;
        Self::from_atom(py, fixed_tensor(ETS.flat, Vec::new(), vec![i, j]))
    }

    /// Match a gamma matrix. Arguments use logical `(mu, i, j)` order.
    #[staticmethod]
    fn gamma(
        py: Python<'_>,
        minkowski_dimension: ConvertibleToExpression,
        mu: ConvertibleToExpression,
        i: ConvertibleToExpression,
        j: ConvertibleToExpression,
    ) -> PyResult<Py<Self>> {
        let dimension = minkowski_dimension.to_expression().expr;
        let logical_ports = vec![
            minkowski_port(&dimension, &mu.to_expression().expr),
            bispinor_port(&Atom::num(4), &i.to_expression().expr),
            bispinor_port(&Atom::num(4), &j.to_expression().expr),
        ];
        Self::built_in(py, &AGS.gamma_strct::<AbstractIndex>(4), logical_ports)
    }

    /// Match a gamma-five matrix in logical `(i, j)` order.
    #[staticmethod]
    fn gamma5(
        py: Python<'_>,
        spinor_dimension: ConvertibleToExpression,
        i: ConvertibleToExpression,
        j: ConvertibleToExpression,
    ) -> PyResult<Py<Self>> {
        let dimension = spinor_dimension.to_expression().expr;
        let logical_ports = vec![
            bispinor_port(&dimension, &i.to_expression().expr),
            bispinor_port(&dimension, &j.to_expression().expr),
        ];
        Self::built_in(py, &AGS.gamma5_strct::<AbstractIndex>(4), logical_ports)
    }

    /// Match a left-chiral projector in logical `(i, j)` order.
    #[staticmethod]
    fn projm(
        py: Python<'_>,
        spinor_dimension: ConvertibleToExpression,
        i: ConvertibleToExpression,
        j: ConvertibleToExpression,
    ) -> PyResult<Py<Self>> {
        let dimension = spinor_dimension.to_expression().expr;
        let logical_ports = vec![
            bispinor_port(&dimension, &i.to_expression().expr),
            bispinor_port(&dimension, &j.to_expression().expr),
        ];
        Self::built_in(py, &AGS.projm_strct::<AbstractIndex>(4), logical_ports)
    }

    /// Match a right-chiral projector in logical `(i, j)` order.
    #[staticmethod]
    fn projp(
        py: Python<'_>,
        spinor_dimension: ConvertibleToExpression,
        i: ConvertibleToExpression,
        j: ConvertibleToExpression,
    ) -> PyResult<Py<Self>> {
        let dimension = spinor_dimension.to_expression().expr;
        let logical_ports = vec![
            bispinor_port(&dimension, &i.to_expression().expr),
            bispinor_port(&dimension, &j.to_expression().expr),
        ];
        Self::built_in(py, &AGS.projp_strct::<AbstractIndex>(4), logical_ports)
    }

    /// Match a sigma matrix. Arguments use logical `(mu, nu, i, j)` order.
    #[staticmethod]
    fn sigma(
        py: Python<'_>,
        minkowski_dimension: ConvertibleToExpression,
        mu: ConvertibleToExpression,
        nu: ConvertibleToExpression,
        i: ConvertibleToExpression,
        j: ConvertibleToExpression,
    ) -> PyResult<Py<Self>> {
        let dimension = minkowski_dimension.to_expression().expr;
        let logical_ports = vec![
            minkowski_port(&dimension, &mu.to_expression().expr),
            minkowski_port(&dimension, &nu.to_expression().expr),
            bispinor_port(&Atom::num(4), &i.to_expression().expr),
            bispinor_port(&Atom::num(4), &j.to_expression().expr),
        ];
        Self::built_in(py, &AGS.sigma_strct::<AbstractIndex>(4), logical_ports)
    }

    /// Match a color structure constant in logical `(a, b, c)` order.
    #[staticmethod]
    fn f(
        py: Python<'_>,
        adjoint_dimension: ConvertibleToExpression,
        a: ConvertibleToExpression,
        b: ConvertibleToExpression,
        c: ConvertibleToExpression,
    ) -> PyResult<Py<Self>> {
        let dimension = adjoint_dimension.to_expression().expr;
        let logical_ports = vec![
            adjoint_port(&dimension, &a.to_expression().expr),
            adjoint_port(&dimension, &b.to_expression().expr),
            adjoint_port(&dimension, &c.to_expression().expr),
        ];
        Self::built_in(py, &CS.f_strct::<AbstractIndex>(1), logical_ports)
    }

    /// Match a color generator in logical `(a, i, j)` order.
    #[staticmethod]
    fn t(
        py: Python<'_>,
        adjoint_dimension: ConvertibleToExpression,
        fundamental_dimension: ConvertibleToExpression,
        a: ConvertibleToExpression,
        i: ConvertibleToExpression,
        j: ConvertibleToExpression,
    ) -> PyResult<Py<Self>> {
        let adjoint_dimension = adjoint_dimension.to_expression().expr;
        let fundamental_dimension = fundamental_dimension.to_expression().expr;
        let logical_ports = vec![
            adjoint_port(&adjoint_dimension, &a.to_expression().expr),
            fundamental_port(&fundamental_dimension, &i.to_expression().expr),
            antifundamental_port(&fundamental_dimension, &j.to_expression().expr),
        ];
        Self::built_in(py, &CS.t_strct::<AbstractIndex>(1, 1), logical_ports)
    }
}

impl ModuleInit for TensorPattern {}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    PortPattern::init(module)?;
    TensorPattern::init(module)
}

#[cfg(test)]
mod tests {
    use super::*;
    use pyo3::types::{PyDict, PyList};
    use spenso::structure::dimension::Dimension;
    use symbolica::symbol;

    use crate::expression::TensorExpression;

    fn wildcard(name: &str) -> PythonExpression {
        PythonExpression::from(Atom::var(symbol!(name)))
    }

    fn assert_replaces(
        target: &Bound<'_, PyAny>,
        pattern: &Bound<'_, PyAny>,
        replacement: i64,
    ) -> PyResult<()> {
        let replaced = target
            .call_method1("replace", (pattern, replacement))?
            .extract::<PythonExpression>()?;
        assert_eq!(replaced.expr, Atom::num(replacement));
        Ok(())
    }

    #[test]
    fn built_in_patterns_translate_logical_ports_to_canonical_atom_order() {
        let dimension = Atom::var(symbol!("D_"));
        let mu = Atom::var(symbol!("mu_"));
        let i = Atom::var(symbol!("i_"));
        let j = Atom::var(symbol!("j_"));
        let logical = vec![
            minkowski_port(&dimension, &mu),
            bispinor_port(&Atom::num(4), &i),
            bispinor_port(&Atom::num(4), &j),
        ];
        let pattern = built_in_pattern(&AGS.gamma_strct::<AbstractIndex>(4), logical);
        let AtomView::Fun(function) = pattern.as_view() else {
            panic!("expected a gamma function pattern")
        };
        let ports = function.iter().collect::<Vec<_>>();

        assert_eq!(function.get_symbol(), AGS.gamma);
        assert!(
            matches!(ports[0], AtomView::Fun(port) if port.get_symbol().has_tag(&SPENSO_TAG.self_dual))
        );
        assert!(
            matches!(ports[2], AtomView::Fun(port) if port.get_symbol() == spenso::structure::representation::LibraryRep::from(Minkowski {}).symbol())
        );
    }

    #[test]
    fn dual_representation_is_indexed_inside_its_wrapper() {
        let representation = ColorAntiFundamental {}.to_symbolic([Atom::var(symbol!("D_"))]);
        let indexed = index_representation(representation, Atom::var(symbol!("i_"))).unwrap();
        let AtomView::Fun(dual) = indexed.as_view() else {
            panic!("expected a dual wrapper")
        };
        let inner = dual.iter().next().unwrap();
        let AtomView::Fun(inner) = inner else {
            panic!("expected an indexed representation")
        };

        assert_eq!(dual.get_symbol(), AIND_SYMBOLS.dind);
        assert_eq!(inner.get_nargs(), 2);
    }

    #[test]
    fn wildcard_heads_require_single_wildcards() {
        assert!(
            wildcard_head("T_", |name| symbolica::try_symbol!(
                name,
                tag = &SPENSO_TAG.tensor
            ))
            .is_ok()
        );
        assert!(
            wildcard_head("T", |name| symbolica::try_symbol!(
                name,
                tag = &SPENSO_TAG.tensor
            ))
            .is_err()
        );
        assert!(
            wildcard_head("T___", |name| symbolica::try_symbol!(
                name,
                tag = &SPENSO_TAG.tensor
            ))
            .is_err()
        );
        assert!(
            wildcard_head("bad head_", |name| symbolica::try_symbol!(
                name,
                tag = &SPENSO_TAG.tensor
            ))
            .is_err()
        );
    }

    #[test]
    fn python_pattern_builders_work_in_variadic_replacements() {
        idenso::representations::initialize();
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let dimension = PythonExpression::from(Atom::var(symbol!("pattern_dimension_")));
            let index = PythonExpression::from(Atom::var(symbol!("pattern_index_")));
            let tail = PythonExpression::from(Atom::var(symbol!("pattern_ports___")));

            let representation = Py::new(
                py,
                SpensoRepresentation {
                    representation: Minkowski {}.new_rep(Dimension::Concrete(4)).cast(),
                },
            )?;
            let port_type = py.get_type::<PortPattern>();
            let exact = port_type.call_method1("exact", (representation, index.clone()))?;
            let any = port_type.call_method1(
                "any",
                ("pattern_rep_", dimension.clone(), index.clone()),
            )?;
            let self_dual = port_type.call_method1(
                "self_dual",
                (
                    "pattern_self_dual_",
                    dimension.clone(),
                    index.clone(),
                ),
            )?;
            let stripped = port_type.call_method1(
                "dualizable",
                ("pattern_dualizable_", dimension.clone()),
            )?;
            let dual_kwargs = PyDict::new(py);
            dual_kwargs.set_item("dual", true)?;
            let dual = port_type.call_method(
                "dualizable",
                (
                    "pattern_dualizable_dual_",
                    dimension.clone(),
                    index.clone(),
                ),
                Some(&dual_kwargs),
            )?;
            for pattern in [&exact, &any, &self_dual, &stripped, &dual] {
                assert!(pattern.is_instance_of::<PortPattern>());
                assert!(pattern.is_instance_of::<PythonExpression>());
            }
            let dual = dual.extract::<PyRef<'_, PortPattern>>()?;
            assert!(
                matches!(dual.as_super().expr.as_view(), AtomView::Fun(function) if function.get_symbol() == AIND_SYMBOLS.dind)
            );

            let tensor_type = py.get_type::<TensorPattern>();
            let fixed_head = Py::new(
                py,
                SpensoName {
                    name: SPENSO_TAG.tensor_symbol("pattern_fixed_tensor"),
                },
            )?;
            let fixed_args = PyList::new(py, [PythonExpression::from(Atom::num(17))])?;
            let fixed_ports = PyList::empty(py);
            fixed_ports.append(&exact)?;
            let fixed_kwargs = PyDict::new(py);
            fixed_kwargs.set_item("args", fixed_args)?;
            fixed_kwargs.set_item("ports", fixed_ports)?;
            let fixed = tensor_type.call((fixed_head,), Some(&fixed_kwargs))?;
            let fixed = fixed.extract::<PyRef<'_, TensorPattern>>()?;
            let AtomView::Fun(fixed_function) = fixed.as_super().expr.as_view() else {
                panic!("expected a fixed-head function pattern")
            };
            let fixed_arguments = fixed_function.iter().collect::<Vec<_>>();
            assert!(matches!(fixed_arguments[0], AtomView::Num(_)));
            assert!(matches!(fixed_arguments[1], AtomView::Fun(_)));

            let ports = PyList::empty(py);
            ports.append(&self_dual)?;
            ports.append(tail)?;
            let wildcard_kwargs = PyDict::new(py);
            wildcard_kwargs.set_item("ports", ports)?;
            let any_tensor = tensor_type.call_method(
                "any",
                ("pattern_tensor_",),
                Some(&wildcard_kwargs),
            )?;
            let vector = tensor_type.call_method(
                "vector",
                ("pattern_vector_",),
                Some(&wildcard_kwargs),
            )?;
            for pattern in [&any_tensor, &vector] {
                assert!(pattern.is_instance_of::<TensorPattern>());
                assert!(pattern.is_instance_of::<PythonExpression>());
            }
            let vector = vector.extract::<PyRef<'_, TensorPattern>>()?;
            let AtomView::Fun(vector_function) = vector.as_super().expr.as_view() else {
                panic!("expected a rank-one function pattern")
            };
            assert!(vector_function.get_symbol().has_tag(&SPENSO_TAG.rank1));

            let concrete_dimension = Atom::num(4);
            let target = PythonExpression::from(fixed_tensor(
                SPENSO_TAG.tensor_symbol("pattern_replacement_target"),
                Vec::new(),
                vec![
                    minkowski_port(&concrete_dimension, &Atom::var(symbol!("mu"))),
                    bispinor_port(&concrete_dimension, &Atom::var(symbol!("i"))),
                ],
            ));
            let target = target.into_pyobject(py)?;
            let replaced = target
                .call_method1("replace", (&any_tensor, 23))?
                .extract::<PythonExpression>()?;
            assert_eq!(replaced.expr, Atom::num(23));

            Ok(())
        })
        .unwrap();
    }

    #[test]
    fn fixed_head_and_port_patterns_work_in_replacements() {
        idenso::representations::initialize();
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let port_type = py.get_type::<PortPattern>();
            let tensor_type = py.get_type::<TensorPattern>();
            let dimension = Atom::num(4);
            let fixed_symbol = SPENSO_TAG.tensor_symbol("pattern_fixed_replacement_target");
            let fixed_head = Py::new(py, SpensoName { name: fixed_symbol })?;
            let target = PythonExpression::from(fixed_tensor(
                fixed_symbol,
                vec![Atom::num(17)],
                vec![minkowski_port(
                    &dimension,
                    &Atom::var(symbol!("pattern_fixed_index")),
                )],
            ))
            .into_pyobject(py)?;
            let representation = Py::new(
                py,
                SpensoRepresentation {
                    representation: Minkowski {}.new_rep(Dimension::Concrete(4)).cast(),
                },
            )?;

            let exact = port_type
                .call_method1("exact", (representation, wildcard("pattern_exact_index_")))?;
            let args = PyList::new(py, [wildcard("pattern_exact_arg_")])?;
            let ports = PyList::new(py, [&exact])?;
            let kwargs = PyDict::new(py);
            kwargs.set_item("args", args)?;
            kwargs.set_item("ports", ports)?;
            let pattern = tensor_type.call((&fixed_head,), Some(&kwargs))?;
            assert_replaces(&target, &pattern, 41)?;

            let any = port_type.call_method1(
                "any",
                (
                    "pattern_any_rep_",
                    wildcard("pattern_any_dimension_"),
                    wildcard("pattern_any_index_"),
                ),
            )?;
            let args = PyList::new(py, [wildcard("pattern_any_arg_")])?;
            let ports = PyList::new(py, [&any])?;
            let kwargs = PyDict::new(py);
            kwargs.set_item("args", args)?;
            kwargs.set_item("ports", ports)?;
            let pattern = tensor_type.call((&fixed_head,), Some(&kwargs))?;
            assert_replaces(&target, &pattern, 43)?;

            let color_symbol = SPENSO_TAG.tensor_symbol("pattern_color_port_target");
            let color_head = Py::new(py, SpensoName { name: color_symbol })?;
            let base_target = PythonExpression::from(fixed_tensor(
                color_symbol,
                Vec::new(),
                vec![fundamental_port(
                    &Atom::num(3),
                    &Atom::var(symbol!("pattern_base_index")),
                )],
            ))
            .into_pyobject(py)?;
            let dual_target = PythonExpression::from(fixed_tensor(
                color_symbol,
                Vec::new(),
                vec![antifundamental_port(
                    &Atom::num(3),
                    &Atom::var(symbol!("pattern_dual_index")),
                )],
            ))
            .into_pyobject(py)?;

            let base = port_type.call_method1(
                "dualizable",
                (
                    "pattern_base_rep_",
                    wildcard("pattern_base_dimension_"),
                    wildcard("pattern_base_index_"),
                ),
            )?;
            let ports = PyList::new(py, [&base])?;
            let kwargs = PyDict::new(py);
            kwargs.set_item("ports", ports)?;
            let pattern = tensor_type.call((&color_head,), Some(&kwargs))?;
            assert_replaces(&base_target, &pattern, 47)?;

            let dual_kwargs = PyDict::new(py);
            dual_kwargs.set_item("dual", true)?;
            let dual = port_type.call_method(
                "dualizable",
                (
                    "pattern_dual_rep_",
                    wildcard("pattern_dual_dimension_"),
                    wildcard("pattern_dual_index_"),
                ),
                Some(&dual_kwargs),
            )?;
            let ports = PyList::new(py, [&dual])?;
            let kwargs = PyDict::new(py);
            kwargs.set_item("ports", ports)?;
            let pattern = tensor_type.call((&color_head,), Some(&kwargs))?;
            assert_replaces(&dual_target, &pattern, 53)?;

            let vector_symbol = SPENSO_TAG.rank_one_tensor_symbol("pattern_stripped_vector_target");
            let vector_target = PythonExpression::from(fixed_tensor(
                vector_symbol,
                vec![Atom::num(5)],
                vec![ColorFundamental {}.to_symbolic([Atom::num(3)])],
            ))
            .into_pyobject(py)?;
            let stripped = port_type.call_method1(
                "dualizable",
                (
                    "pattern_stripped_rep_",
                    wildcard("pattern_stripped_dimension_"),
                ),
            )?;
            let args = PyList::new(py, [wildcard("pattern_vector_arg_")])?;
            let ports = PyList::new(py, [&stripped])?;
            let kwargs = PyDict::new(py);
            kwargs.set_item("args", args)?;
            kwargs.set_item("ports", ports)?;
            let pattern =
                tensor_type.call_method("vector", ("pattern_rank_one_head_",), Some(&kwargs))?;
            assert_replaces(&vector_target, &pattern, 59)?;
            Ok(())
        })
        .unwrap();
    }

    #[test]
    fn every_builtin_pattern_shortcut_matches_its_typed_factory() {
        idenso::representations::initialize();
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let expression_type = py.get_type::<TensorExpression>();
            let pattern_type = py.get_type::<TensorPattern>();
            let representation = Py::new(
                py,
                SpensoRepresentation {
                    representation: Minkowski {}.new_rep(Dimension::Concrete(4)).cast(),
                },
            )?;

            macro_rules! assert_builtin {
                ($name:literal, $factory_args:expr, $indices:expr, $pattern_args:expr, $value:expr) => {{
                    let target = expression_type
                        .call_method1($name, $factory_args)?
                        .call1($indices)?
                        .call_method0("to_expression")?;
                    let pattern = pattern_type.call_method1($name, $pattern_args)?;
                    assert_replaces(&target, &pattern, $value)?;
                }};
            }

            assert_builtin!(
                "g",
                (&representation,),
                ("pattern_g_i", "pattern_g_j"),
                (
                    &representation,
                    wildcard("pattern_g_i_"),
                    wildcard("pattern_g_j_")
                ),
                61
            );
            assert_builtin!(
                "flat",
                (&representation,),
                ("pattern_flat_i", "pattern_flat_j"),
                (
                    &representation,
                    wildcard("pattern_flat_i_"),
                    wildcard("pattern_flat_j_")
                ),
                67
            );
            assert_builtin!(
                "gamma",
                (7,),
                ("pattern_gamma_mu", "pattern_gamma_i", "pattern_gamma_j"),
                (
                    wildcard("pattern_gamma_dimension_"),
                    wildcard("pattern_gamma_mu_"),
                    wildcard("pattern_gamma_i_"),
                    wildcard("pattern_gamma_j_")
                ),
                71
            );
            assert_builtin!(
                "gamma5",
                (6,),
                ("pattern_gamma5_i", "pattern_gamma5_j"),
                (
                    wildcard("pattern_gamma5_dimension_"),
                    wildcard("pattern_gamma5_i_"),
                    wildcard("pattern_gamma5_j_")
                ),
                73
            );
            assert_builtin!(
                "projm",
                (6,),
                ("pattern_projm_i", "pattern_projm_j"),
                (
                    wildcard("pattern_projm_dimension_"),
                    wildcard("pattern_projm_i_"),
                    wildcard("pattern_projm_j_")
                ),
                79
            );
            assert_builtin!(
                "projp",
                (6,),
                ("pattern_projp_i", "pattern_projp_j"),
                (
                    wildcard("pattern_projp_dimension_"),
                    wildcard("pattern_projp_i_"),
                    wildcard("pattern_projp_j_")
                ),
                83
            );
            assert_builtin!(
                "sigma",
                (7,),
                (
                    "pattern_sigma_mu",
                    "pattern_sigma_nu",
                    "pattern_sigma_i",
                    "pattern_sigma_j"
                ),
                (
                    wildcard("pattern_sigma_dimension_"),
                    wildcard("pattern_sigma_mu_"),
                    wildcard("pattern_sigma_nu_"),
                    wildcard("pattern_sigma_i_"),
                    wildcard("pattern_sigma_j_")
                ),
                89
            );
            assert_builtin!(
                "f",
                (8,),
                ("pattern_f_a", "pattern_f_b", "pattern_f_c"),
                (
                    wildcard("pattern_f_dimension_"),
                    wildcard("pattern_f_a_"),
                    wildcard("pattern_f_b_"),
                    wildcard("pattern_f_c_")
                ),
                97
            );
            assert_builtin!(
                "t",
                (8, 3),
                ("pattern_t_a", "pattern_t_i", "pattern_t_j"),
                (
                    wildcard("pattern_t_adjoint_dimension_"),
                    wildcard("pattern_t_fundamental_dimension_"),
                    wildcard("pattern_t_a_"),
                    wildcard("pattern_t_i_"),
                    wildcard("pattern_t_j_")
                ),
                101
            );
            Ok(())
        })
        .unwrap();
    }

    #[cfg(feature = "python_stubgen")]
    #[test]
    fn typed_tensor_stubs_preserve_public_signatures() {
        use pyo3_stub_gen::{
            StubInfo,
            generate::MethodType,
            type_info::{ParameterDefault, ParameterKind, PyMethodsInfo},
        };

        fn parameters<T: 'static>(method_name: &str) -> Vec<&'static str> {
            pyo3_stub_gen::inventory::iter::<PyMethodsInfo>
                .into_iter()
                .filter(|info| (info.struct_id)() == std::any::TypeId::of::<T>())
                .flat_map(|info| info.methods)
                .find(|method| method.name == method_name)
                .unwrap_or_else(|| panic!("{method_name} is missing from generated stubs"))
                .parameters
                .iter()
                .map(|parameter| parameter.name)
                .collect()
        }

        let methods = pyo3_stub_gen::inventory::iter::<PyMethodsInfo>
            .into_iter()
            .filter(|info| (info.struct_id)() == std::any::TypeId::of::<TensorPattern>())
            .flat_map(|info| info.methods)
            .collect::<Vec<_>>();
        let constructor = methods
            .iter()
            .find(|method| method.r#type == MethodType::New)
            .expect("TensorPattern must expose a constructor in generated stubs");
        assert_eq!((constructor.r#return)().name, "TensorPattern");

        for method_name in ["__new__", "any", "vector"] {
            let method = methods
                .iter()
                .find(|method| method.name == method_name)
                .unwrap_or_else(|| panic!("TensorPattern.{method_name} is missing from stubs"));
            for parameter_name in ["args", "ports"] {
                let parameter = method
                    .parameters
                    .iter()
                    .find(|parameter| parameter.name == parameter_name)
                    .unwrap_or_else(|| {
                        panic!("TensorPattern.{method_name} is missing {parameter_name}")
                    });
                assert_eq!(parameter.kind, ParameterKind::KeywordOnly);
                let ParameterDefault::Expr(format) = &parameter.default else {
                    panic!("TensorPattern.{method_name}.{parameter_name} has no list default")
                };
                assert_eq!(format(), "[]");
            }
        }

        for (method, expected) in [
            ("g", &["rep"][..]),
            ("flat", &["rep"]),
            ("gamma", &["minkowski_dimension"]),
            ("gamma5", &["spinor_dimension"]),
            ("projm", &["spinor_dimension"]),
            ("projp", &["spinor_dimension"]),
            ("sigma", &["minkowski_dimension"]),
            ("f", &["adjoint_dimension"]),
            ("t", &["adjoint_dimension", "fundamental_dimension"]),
        ] {
            assert_eq!(parameters::<TensorExpression>(method), expected);
        }
        assert_eq!(parameters::<PortPattern>("exact"), ["rep", "index"]);
        for (method, expected) in [
            ("g", &["rep_pattern", "i", "j"][..]),
            ("flat", &["rep_pattern", "i", "j"]),
            ("gamma", &["minkowski_dimension", "mu", "i", "j"]),
            ("gamma5", &["spinor_dimension", "i", "j"]),
            ("projm", &["spinor_dimension", "i", "j"]),
            ("projp", &["spinor_dimension", "i", "j"]),
            ("sigma", &["minkowski_dimension", "mu", "nu", "i", "j"]),
            ("f", &["adjoint_dimension", "a", "b", "c"]),
            (
                "t",
                &["adjoint_dimension", "fundamental_dimension", "a", "i", "j"],
            ),
        ] {
            assert_eq!(parameters::<TensorPattern>(method), expected);
        }

        let stubs = StubInfo::from_project_root(
            "symbolica.core".to_owned(),
            std::path::PathBuf::from("unused-stub-output"),
        )
        .expect("the downstream community stub generator must accept this inventory");
        let rendered = stubs
            .modules
            .get("symbolica.community.spenso")
            .expect("the Spenso community stub module must be generated")
            .to_string();
        for class in ["PortPattern", "TensorExpression", "TensorPattern"] {
            assert!(
                rendered.contains(&format!("class {class}")),
                "generated Spenso stub is missing {class}"
            );
        }
        for signature in [
            "def gamma(minkowski_dimension:",
            "def gamma5(spinor_dimension:",
            "def t(adjoint_dimension:",
            "def exact(rep:",
            "def g(rep_pattern:",
        ] {
            assert!(
                rendered.contains(signature),
                "generated Spenso stub is missing {signature}"
            );
        }
    }

    #[test]
    fn pattern_public_surface_has_runtime_docstrings() {
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            for (pattern_type, methods) in [
                (
                    py.get_type::<PortPattern>(),
                    &["exact", "any", "self_dual", "dualizable"][..],
                ),
                (
                    py.get_type::<TensorPattern>(),
                    &[
                        "any", "vector", "g", "flat", "gamma", "gamma5", "projm", "projp", "sigma",
                        "f", "t",
                    ][..],
                ),
            ] {
                let class_documentation = pattern_type
                    .getattr("__doc__")?
                    .extract::<Option<String>>()?;
                assert!(class_documentation.is_some_and(|doc| !doc.trim().is_empty()));
                for method in methods {
                    let documentation = pattern_type
                        .getattr(method)?
                        .getattr("__doc__")?
                        .extract::<Option<String>>()?;
                    assert!(
                        documentation.is_some_and(|doc| !doc.trim().is_empty()),
                        "{}.{} is missing its Python docstring",
                        pattern_type.name()?,
                        method,
                    );
                }
            }
            Ok(())
        })
        .unwrap();
    }
}
