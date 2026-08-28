//! Typed Python values for configuring Linnest's Typst renderer.
//!
//! This module deliberately models a closed subset of Typst's native value
//! language. Strings are always data, never source, and executable Typst
//! values can only enter through an explicit module export reference.

use std::collections::{BTreeMap, BTreeSet};
use std::path::{Path, PathBuf};

use pyo3::class::gc::{PyTraverseError, PyVisit};
use pyo3::exceptions::{PyOverflowError, PyTypeError, PyValueError};
use pyo3::prelude::*;
use pyo3::sync::PyOnceLock;
use pyo3::types::{PyAny, PyBool, PyDict, PyFloat, PyInt, PyList, PyModule, PyString, PyTuple};
use pyo3_stub_gen::derive::{gen_stub_pyclass, gen_stub_pyclass_enum, gen_stub_pymethods};
use serde::{Deserialize, Serialize};

use crate::drawing::{DrawingKind, PyEdgeDrawing, PyHalfEdgeDrawing, PyNodeDrawing};

const MAX_NATIVE_DEPTH: usize = 64;
static AUTO_SINGLETON: PyOnceLock<Py<PyAuto>> = PyOnceLock::new();
static INHERIT_SINGLETON: PyOnceLock<Py<PyInherit>> = PyOnceLock::new();

fn finite(value: f64, name: &str) -> PyResult<f64> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(PyValueError::new_err(format!(
            "{name} must be finite, got {value}"
        )))
    }
}

fn non_negative(value: f64, name: &str) -> PyResult<f64> {
    let value = finite(value, name)?;
    if value >= 0.0 {
        Ok(value)
    } else {
        Err(PyValueError::new_err(format!(
            "{name} must be non-negative, got {value}"
        )))
    }
}

fn identifier(value: &str, what: &str) -> PyResult<String> {
    let mut chars = value.chars();
    let Some(first) = chars.next() else {
        return Err(PyValueError::new_err(format!("{what} cannot be empty")));
    };
    if !(first == '_' || first.is_alphabetic())
        || !chars.all(|c| c == '_' || c == '-' || c.is_alphanumeric())
    {
        return Err(PyValueError::new_err(format!(
            "invalid Typst {what} {value:?}"
        )));
    }
    Ok(value.to_owned())
}

fn option_name(value: &str) -> PyResult<String> {
    identifier(value, "named argument")
}

/// Renderer-owned module dependency collected from typed native values.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize)]
#[serde(tag = "kind", content = "value", rename_all = "kebab-case")]
pub(crate) enum TypstModuleSource {
    File(PathBuf),
    Package(String),
}

/// One deterministic import required by a serialized Typst expression.
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct TypstImport {
    pub(crate) alias: String,
    pub(crate) source: TypstModuleSource,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case")]
enum SymbolKind {
    Value,
    Content,
    Function,
}

impl SymbolKind {
    fn name(self) -> &'static str {
        match self {
            Self::Value => "value",
            Self::Content => "content",
            Self::Function => "function",
        }
    }
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct ModuleSymbol {
    module: TypstModuleSource,
    path: Vec<String>,
    kind: SymbolKind,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct CallArguments {
    positional: Vec<NativeValue>,
    named: BTreeMap<String, NativeValue>,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(tag = "kind", content = "value", rename_all = "kebab-case")]
enum TypstExpression {
    Symbol(ModuleSymbol),
    Call {
        callee: Box<TypstExpression>,
        arguments: CallArguments,
    },
    Bind {
        function: Box<TypstExpression>,
        arguments: CallArguments,
    },
}

impl TypstExpression {
    fn collect_modules(&self, modules: &mut BTreeSet<TypstModuleSource>) {
        match self {
            Self::Symbol(symbol) => {
                modules.insert(symbol.module.clone());
            }
            Self::Call { callee, arguments }
            | Self::Bind {
                function: callee,
                arguments,
            } => {
                callee.collect_modules(modules);
                arguments.collect_modules(modules);
            }
        }
    }
}

impl CallArguments {
    fn collect_modules(&self, modules: &mut BTreeSet<TypstModuleSource>) {
        for value in &self.positional {
            value.collect_modules(modules);
        }
        for value in self.named.values() {
            value.collect_modules(modules);
        }
    }

    fn validate(&self, depth: usize) -> PyResult<()> {
        for value in &self.positional {
            value.validate(depth + 1)?;
            if matches!(value, NativeValue::Inherit) {
                return Err(PyValueError::new_err(
                    "INHERIT is not valid as a positional Typst argument",
                ));
            }
        }
        for (name, value) in &self.named {
            identifier(name, "named argument")?;
            value.validate(depth + 1)?;
            if matches!(value, NativeValue::Inherit) {
                return Err(PyValueError::new_err(
                    "INHERIT is not valid in serialized Typst arguments",
                ));
            }
        }
        Ok(())
    }
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case")]
enum LengthUnit {
    Pt,
    Mm,
    Cm,
    In,
    Em,
}

impl LengthUnit {
    fn parse(value: &str) -> PyResult<Self> {
        match value {
            "pt" => Ok(Self::Pt),
            "mm" => Ok(Self::Mm),
            "cm" => Ok(Self::Cm),
            "in" => Ok(Self::In),
            "em" => Ok(Self::Em),
            _ => Err(PyValueError::new_err(format!(
                "unknown length unit {value:?}; expected pt, mm, cm, in, or em"
            ))),
        }
    }

    fn name(&self) -> &'static str {
        match self {
            Self::Pt => "pt",
            Self::Mm => "mm",
            Self::Cm => "cm",
            Self::In => "in",
            Self::Em => "em",
        }
    }
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(tag = "kind", content = "value", rename_all = "kebab-case")]
enum ColorValue {
    Named(String),
    Hex(String),
    Rgba(u8, u8, u8, Option<u8>),
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(tag = "kind", content = "value", rename_all = "kebab-case")]
enum DashValue {
    Named(String),
    Pattern {
        array: Vec<NativeValue>,
        phase: Option<Box<NativeValue>>,
    },
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct TextValue {
    text: String,
    options: BTreeMap<String, NativeValue>,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
enum MathScript {
    Symbol(String),
    Integer(i64),
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct MathValue {
    symbol: String,
    subscript: Option<MathScript>,
    superscript: Option<MathScript>,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum EnumKind {
    LayoutAlgorithm,
    LayoutDirection,
    LabelLayout,
    RankAlignment,
    LayoutNodes,
    Placement,
    Compass,
    Routing,
    RoutePoints,
    Anchor,
    Pattern,
    EdgeLengthResolution,
    DebugLevel,
    StrokeCap,
    StrokeJoin,
    TextStyle,
    DashPattern,
    MarkSymbol,
    MarkPosition,
    MarkOrientation,
    MarkDirection,
}

impl EnumKind {
    const fn name(self) -> &'static str {
        match self {
            Self::LayoutAlgorithm => "LayoutAlgorithm",
            Self::LayoutDirection => "LayoutDirection",
            Self::LabelLayout => "LabelLayout",
            Self::RankAlignment => "RankAlignment",
            Self::LayoutNodes => "LayoutNodes",
            Self::Placement => "Placement",
            Self::Compass => "Compass",
            Self::Routing => "Routing",
            Self::RoutePoints => "RoutePoints",
            Self::Anchor => "Anchor",
            Self::Pattern => "Pattern",
            Self::EdgeLengthResolution => "EdgeLengthResolution",
            Self::DebugLevel => "DebugLevel",
            Self::StrokeCap => "StrokeCap",
            Self::StrokeJoin => "StrokeJoin",
            Self::TextStyle => "TextStyle",
            Self::DashPattern => "DashPattern",
            Self::MarkSymbol => "MarkSymbol",
            Self::MarkPosition => "MarkPosition",
            Self::MarkOrientation => "MarkOrientation",
            Self::MarkDirection => "MarkDirection",
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(tag = "enum", content = "value", rename_all = "kebab-case")]
enum NativeEnum {
    LayoutAlgorithm(PyLayoutAlgorithm),
    LayoutDirection(PyLayoutDirection),
    LabelLayout(PyLabelLayout),
    RankAlignment(PyRankAlignment),
    LayoutNodes(PyLayoutNodes),
    Placement(PyPlacement),
    Compass(PyCompass),
    Routing(PyRouting),
    RoutePoints(PyRoutePoints),
    Anchor(PyAnchor),
    Pattern(PyPattern),
    EdgeLengthResolution(PyEdgeLengthResolution),
    DebugLevel(PyDebugLevel),
    StrokeCap(PyStrokeCap),
    StrokeJoin(PyStrokeJoin),
    TextStyle(PyTextStyle),
    DashPattern(PyDashPattern),
    MarkSymbol(PyMarkSymbol),
    MarkPosition(PyMarkPosition),
    MarkOrientation(PyMarkOrientation),
    MarkDirection(PyMarkDirection),
}

impl NativeEnum {
    const fn kind(self) -> EnumKind {
        match self {
            Self::LayoutAlgorithm(_) => EnumKind::LayoutAlgorithm,
            Self::LayoutDirection(_) => EnumKind::LayoutDirection,
            Self::LabelLayout(_) => EnumKind::LabelLayout,
            Self::RankAlignment(_) => EnumKind::RankAlignment,
            Self::LayoutNodes(_) => EnumKind::LayoutNodes,
            Self::Placement(_) => EnumKind::Placement,
            Self::Compass(_) => EnumKind::Compass,
            Self::Routing(_) => EnumKind::Routing,
            Self::RoutePoints(_) => EnumKind::RoutePoints,
            Self::Anchor(_) => EnumKind::Anchor,
            Self::Pattern(_) => EnumKind::Pattern,
            Self::EdgeLengthResolution(_) => EnumKind::EdgeLengthResolution,
            Self::DebugLevel(_) => EnumKind::DebugLevel,
            Self::StrokeCap(_) => EnumKind::StrokeCap,
            Self::StrokeJoin(_) => EnumKind::StrokeJoin,
            Self::TextStyle(_) => EnumKind::TextStyle,
            Self::DashPattern(_) => EnumKind::DashPattern,
            Self::MarkSymbol(_) => EnumKind::MarkSymbol,
            Self::MarkPosition(_) => EnumKind::MarkPosition,
            Self::MarkOrientation(_) => EnumKind::MarkOrientation,
            Self::MarkDirection(_) => EnumKind::MarkDirection,
        }
    }

    fn to_py(self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        let value = match self {
            Self::LayoutAlgorithm(value) => Py::new(py, value)?.into_any(),
            Self::LayoutDirection(value) => Py::new(py, value)?.into_any(),
            Self::LabelLayout(value) => Py::new(py, value)?.into_any(),
            Self::RankAlignment(value) => Py::new(py, value)?.into_any(),
            Self::LayoutNodes(value) => Py::new(py, value)?.into_any(),
            Self::Placement(value) => Py::new(py, value)?.into_any(),
            Self::Compass(value) => Py::new(py, value)?.into_any(),
            Self::Routing(value) => Py::new(py, value)?.into_any(),
            Self::RoutePoints(value) => Py::new(py, value)?.into_any(),
            Self::Anchor(value) => Py::new(py, value)?.into_any(),
            Self::Pattern(value) => Py::new(py, value)?.into_any(),
            Self::EdgeLengthResolution(value) => Py::new(py, value)?.into_any(),
            Self::DebugLevel(value) => Py::new(py, value)?.into_any(),
            Self::StrokeCap(value) => Py::new(py, value)?.into_any(),
            Self::StrokeJoin(value) => Py::new(py, value)?.into_any(),
            Self::TextStyle(value) => Py::new(py, value)?.into_any(),
            Self::DashPattern(value) => Py::new(py, value)?.into_any(),
            Self::MarkSymbol(value) => Py::new(py, value)?.into_any(),
            Self::MarkPosition(value) => Py::new(py, value)?.into_any(),
            Self::MarkOrientation(value) => Py::new(py, value)?.into_any(),
            Self::MarkDirection(value) => Py::new(py, value)?.into_any(),
        };
        Ok(value)
    }
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(tag = "type", content = "value", rename_all = "kebab-case")]
enum NativeValue {
    Inherit,
    None,
    Auto,
    Bool(bool),
    Int(i64),
    Float(f64),
    String(String),
    Enum(NativeEnum),
    Array(Vec<Self>),
    Dict(BTreeMap<String, Self>),
    Length(f64, LengthUnit),
    Ratio(f64),
    RelativeLength {
        ratio: Option<f64>,
        length: Option<(f64, LengthUnit)>,
    },
    Angle(f64, String),
    Fraction(f64),
    Color(ColorValue),
    Stroke(BTreeMap<String, Self>),
    Dash(DashValue),
    Insets(BTreeMap<String, Self>),
    Mark(BTreeMap<String, Self>),
    Text(TextValue),
    Math(MathValue),
    Expression(TypstExpression),
}

impl NativeValue {
    fn collect_modules(&self, modules: &mut BTreeSet<TypstModuleSource>) {
        match self {
            Self::Array(values) => {
                for value in values {
                    value.collect_modules(modules);
                }
            }
            Self::Dict(values)
            | Self::Stroke(values)
            | Self::Insets(values)
            | Self::Mark(values) => {
                for value in values.values() {
                    value.collect_modules(modules);
                }
            }
            Self::Dash(DashValue::Pattern { array, phase }) => {
                for value in array {
                    value.collect_modules(modules);
                }
                if let Some(phase) = phase {
                    phase.collect_modules(modules);
                }
            }
            Self::Text(value) => {
                for option in value.options.values() {
                    option.collect_modules(modules);
                }
            }
            Self::Expression(expression) => expression.collect_modules(modules),
            _ => {}
        }
    }

    fn to_clinnet(
        &self,
        aliases: &BTreeMap<TypstModuleSource, String>,
    ) -> PyResult<clinnet::TypstValue> {
        use clinnet::TypstValue as Value;

        Ok(match self {
            Self::Inherit => {
                return Err(PyValueError::new_err(
                    "INHERIT cannot be serialized outside an option dictionary",
                ));
            }
            Self::None => Value::None,
            Self::Auto => Value::Auto,
            Self::Bool(value) => Value::Bool(*value),
            Self::Int(value) => Value::Integer(*value),
            Self::Float(value) => Value::Float(*value),
            Self::String(value) => Value::String(value.clone()),
            Self::Enum(value) => value.to_clinnet(),
            Self::Array(values) => Value::Array(
                values
                    .iter()
                    .map(|value| value.to_clinnet(aliases))
                    .collect::<PyResult<_>>()?,
            ),
            Self::Dict(values) | Self::Stroke(values) | Self::Mark(values) => {
                Value::Dictionary(clinnet_dict(values, aliases)?)
            }
            Self::Insets(values) => {
                if let Some(value) = values.get("all") {
                    if !matches!(value, Self::Inherit) {
                        return value.to_clinnet(aliases);
                    }
                }
                Value::Dictionary(clinnet_dict(values, aliases)?)
            }
            Self::Length(value, unit) => Value::Length(*value, unit.to_clinnet()),
            Self::Ratio(value) => Value::Ratio(*value),
            Self::RelativeLength { ratio, length } => Value::RelativeLength {
                ratio: *ratio,
                length: length
                    .as_ref()
                    .map(|(value, unit)| (*value, unit.to_clinnet())),
            },
            Self::Angle(value, unit) => Value::Angle(
                *value,
                match unit.as_str() {
                    "deg" => clinnet::TypstAngleUnit::Deg,
                    "rad" => clinnet::TypstAngleUnit::Rad,
                    _ => {
                        return Err(PyValueError::new_err(format!(
                            "unsupported Typst angle unit {unit:?}"
                        )));
                    }
                },
            ),
            Self::Fraction(value) => Value::Fraction(*value),
            Self::Color(ColorValue::Named(value)) => Value::NamedColor(value.clone()),
            Self::Color(ColorValue::Hex(value)) => Value::HexColor(value.clone()),
            Self::Color(ColorValue::Rgba(red, green, blue, alpha)) => {
                Value::Rgba(*red, *green, *blue, *alpha)
            }
            Self::Dash(DashValue::Named(value)) => Value::String(value.clone()),
            Self::Dash(DashValue::Pattern { array, phase }) => {
                let array = Value::Array(
                    array
                        .iter()
                        .map(|value| value.to_clinnet(aliases))
                        .collect::<PyResult<_>>()?,
                );
                if let Some(phase) = phase {
                    Value::Dictionary(BTreeMap::from([
                        ("array".to_owned(), array),
                        ("phase".to_owned(), phase.to_clinnet(aliases)?),
                    ]))
                } else {
                    array
                }
            }
            Self::Text(value) => Value::Text {
                text: value.text.clone(),
                options: clinnet_dict(&value.options, aliases)?,
            },
            Self::Math(value) => Value::Math {
                symbol: value.symbol.clone(),
                subscript: value.subscript.as_ref().map(MathScript::to_clinnet),
                superscript: value.superscript.as_ref().map(MathScript::to_clinnet),
            },
            Self::Expression(expression) => expression.to_clinnet(aliases)?,
        })
    }

    fn clinnet_config(&self) -> PyResult<(clinnet::TypstConfig, Vec<TypstImport>)> {
        let mut modules = BTreeSet::new();
        self.collect_modules(&mut modules);
        let imports = modules
            .into_iter()
            .enumerate()
            .map(|(index, source)| TypstImport {
                alias: format!("_linnet_m{index}"),
                source,
            })
            .collect::<Vec<_>>();
        let aliases = imports
            .iter()
            .map(|import| (import.source.clone(), import.alias.clone()))
            .collect();
        let clinnet::TypstValue::Dictionary(fields) = self.to_clinnet(&aliases)? else {
            return Err(PyTypeError::new_err(
                "render configuration must be a dictionary",
            ));
        };
        let config = clinnet::TypstConfig::new(fields)
            .map_err(|error| PyValueError::new_err(error.to_string()))?;
        Ok((config, imports))
    }

    fn kind(&self) -> &'static str {
        match self {
            Self::Inherit => "INHERIT",
            Self::None => "none",
            Self::Auto => "auto",
            Self::Bool(_) => "bool",
            Self::Int(_) => "int",
            Self::Float(_) => "float",
            Self::String(_) => "string",
            Self::Enum(value) => value.kind().name(),
            Self::Array(_) => "array",
            Self::Dict(_) => "dictionary",
            Self::Length(..) => "length",
            Self::Ratio(_) => "ratio",
            Self::RelativeLength { .. } => "relative length",
            Self::Angle(..) => "angle",
            Self::Fraction(_) => "fraction",
            Self::Color(_) => "color",
            Self::Stroke(_) => "stroke",
            Self::Dash(_) => "dash",
            Self::Insets(_) => "insets",
            Self::Mark(_) => "mark",
            Self::Text(_) => "text content",
            Self::Math(_) => "math content",
            Self::Expression(_) => "module expression",
        }
    }
}

impl NativeEnum {
    fn to_clinnet(self) -> clinnet::TypstValue {
        let value = match self {
            Self::DebugLevel(value) => {
                return clinnet::TypstValue::Integer(i64::from(value.level()));
            }
            Self::LayoutAlgorithm(value) => value.typst_name(),
            Self::LayoutDirection(value) => value.typst_name(),
            Self::LabelLayout(value) => value.typst_name(),
            Self::RankAlignment(value) => value.typst_name(),
            Self::LayoutNodes(value) => value.typst_name(),
            Self::Placement(value) => value.typst_name(),
            Self::Compass(value) => value.typst_name(),
            Self::Routing(value) => value.typst_name(),
            Self::RoutePoints(value) => value.typst_name(),
            Self::Anchor(value) => value.typst_name(),
            Self::Pattern(value) => value.typst_name(),
            Self::EdgeLengthResolution(value) => value.typst_name(),
            Self::StrokeCap(value) => value.typst_name(),
            Self::StrokeJoin(value) => value.typst_name(),
            Self::TextStyle(value) => value.typst_name(),
            Self::DashPattern(value) => value.typst_name(),
            Self::MarkSymbol(value) => value.typst_name(),
            Self::MarkPosition(value) => value.typst_name(),
            Self::MarkOrientation(value) => value.typst_name(),
            Self::MarkDirection(value) => value.typst_name(),
        };
        clinnet::TypstValue::String(value.to_owned())
    }
}

impl LengthUnit {
    fn to_clinnet(&self) -> clinnet::TypstLengthUnit {
        match self {
            Self::Pt => clinnet::TypstLengthUnit::Pt,
            Self::Mm => clinnet::TypstLengthUnit::Mm,
            Self::Cm => clinnet::TypstLengthUnit::Cm,
            Self::In => clinnet::TypstLengthUnit::In,
            Self::Em => clinnet::TypstLengthUnit::Em,
        }
    }
}

impl MathScript {
    fn to_clinnet(&self) -> clinnet::TypstMathScript {
        match self {
            Self::Symbol(value) => clinnet::TypstMathScript::Symbol(value.clone()),
            Self::Integer(value) => clinnet::TypstMathScript::Integer(*value),
        }
    }
}

impl TypstExpression {
    fn to_clinnet(
        &self,
        aliases: &BTreeMap<TypstModuleSource, String>,
    ) -> PyResult<clinnet::TypstValue> {
        Ok(match self {
            Self::Symbol(symbol) => clinnet::TypstValue::ModuleSymbol {
                alias: aliases.get(&symbol.module).cloned().ok_or_else(|| {
                    PyValueError::new_err("internal error: missing Typst module alias")
                })?,
                path: symbol.path.clone(),
            },
            Self::Call { callee, arguments } => clinnet::TypstValue::Call {
                callee: Box::new(callee.to_clinnet(aliases)?),
                positional: clinnet_values(&arguments.positional, aliases)?,
                named: clinnet_dict(&arguments.named, aliases)?,
            },
            Self::Bind {
                function,
                arguments,
            } => clinnet::TypstValue::Bind {
                function: Box::new(function.to_clinnet(aliases)?),
                positional: clinnet_values(&arguments.positional, aliases)?,
                named: clinnet_dict(&arguments.named, aliases)?,
            },
        })
    }
}

fn clinnet_values(
    values: &[NativeValue],
    aliases: &BTreeMap<TypstModuleSource, String>,
) -> PyResult<Vec<clinnet::TypstValue>> {
    values
        .iter()
        .map(|value| value.to_clinnet(aliases))
        .collect()
}

fn clinnet_dict(
    values: &BTreeMap<String, NativeValue>,
    aliases: &BTreeMap<TypstModuleSource, String>,
) -> PyResult<BTreeMap<String, clinnet::TypstValue>> {
    values
        .iter()
        .filter(|(_, value)| !matches!(value, NativeValue::Inherit))
        .map(|(name, value)| Ok((name.clone(), value.to_clinnet(aliases)?)))
        .collect()
}

/// Type of the `AUTO` sentinel, which requests automatic selection.
#[gen_stub_pyclass]
#[pyclass(skip_from_py_object, frozen, name = "Auto")]
#[derive(Clone, Copy, Debug)]
struct PyAuto;

pyo3_stub_gen::module_variable!("linnet_py", "AUTO", PyAuto);

#[gen_stub_pymethods]
#[pymethods]
impl PyAuto {
    fn __repr__(&self) -> &'static str {
        "AUTO"
    }

    fn __copy__(&self, py: Python<'_>) -> Py<PyAuto> {
        auto_singleton(py)
    }

    fn __deepcopy__(&self, py: Python<'_>, _memo: &Bound<'_, PyAny>) -> Py<PyAuto> {
        auto_singleton(py)
    }
}

/// Type of the `INHERIT` sentinel, which preserves a lower-precedence setting.
#[gen_stub_pyclass]
#[pyclass(skip_from_py_object, frozen, name = "Inherit")]
#[derive(Clone, Copy, Debug)]
struct PyInherit;

pyo3_stub_gen::module_variable!("linnet_py", "INHERIT", PyInherit);

#[gen_stub_pymethods]
#[pymethods]
impl PyInherit {
    fn __repr__(&self) -> &'static str {
        "INHERIT"
    }

    fn __copy__(&self, py: Python<'_>) -> Py<PyInherit> {
        inherit_singleton(py)
    }

    fn __deepcopy__(&self, py: Python<'_>, _memo: &Bound<'_, PyAny>) -> Py<PyInherit> {
        inherit_singleton(py)
    }
}

macro_rules! typst_string_enum {
    (
        $(#[$meta:meta])*
        $rust:ident, $python:literal, $native:ident {
            $($variant:ident => $value:literal),+ $(,)?
        }
    ) => {
        $(#[$meta])*
        #[gen_stub_pyclass_enum]
        #[pyclass(from_py_object, eq, eq_int, name = $python)]
        #[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
        #[serde(rename_all = "kebab-case")]
        enum $rust {
            $($variant),+
        }

        impl $rust {
            const fn typst_name(self) -> &'static str {
                match self {
                    $(Self::$variant => $value),+
                }
            }
        }

        impl From<$rust> for NativeEnum {
            fn from(value: $rust) -> Self {
                Self::$native(value)
            }
        }
    };
}

typst_string_enum! {
    /// Algorithm used by a Linnest layout pass.
    PyLayoutAlgorithm, "LayoutAlgorithm", LayoutAlgorithm {
        Force => "force",
        Anneal => "anneal",
        Tree => "tree",
        Dot => "dot",
        StableLayered => "stable-layered",
    }
}

typst_string_enum! {
    /// Rank direction used by tree and layered layouts.
    PyLayoutDirection, "LayoutDirection", LayoutDirection {
        Down => "down",
        Right => "right",
    }
}

typst_string_enum! {
    /// Strategy used to position edge labels.
    PyLabelLayout, "LabelLayout", LabelLayout {
        Normal => "normal",
        DanglingTangent => "dangling-tangent",
        FixedLength => "fixed-length",
    }
}

typst_string_enum! {
    /// Alignment within a layered-layout rank.
    PyRankAlignment, "RankAlignment", RankAlignment {
        Center => "center",
        Start => "start",
        Left => "left",
        End => "end",
        Right => "right",
    }
}

typst_string_enum! {
    /// Whether a pass lays out or preserves its selected nodes.
    PyLayoutNodes, "LayoutNodes", LayoutNodes {
        Layout => "layout",
        Fixed => "fixed",
    }
}

typst_string_enum! {
    /// Coordinate-constraint mode for an explicit placement.
    PyPlacement, "Placement", Placement {
        Pin => "pin",
        Start => "start",
    }
}

typst_string_enum! {
    /// A finite DOT compass point.
    PyCompass, "Compass", Compass {
        N => "n",
        NE => "ne",
        E => "e",
        SE => "se",
        S => "s",
        SW => "sw",
        W => "w",
        NW => "nw",
        Center => "c",
        Any => "_",
    }
}

typst_string_enum! {
    /// Route construction used for an anchored edge.
    PyRouting, "Routing", Routing {
        EdgePosition => "edge-pos",
        Direct => "direct",
        HobbyThrough => "hobby-through",
        StraightThrough => "straight-through",
    }
}

typst_string_enum! {
    /// Whether layout-provided route points participate in edge routing.
    PyRoutePoints, "RoutePoints", RoutePoints {
        Ignore => "ignore",
        Through => "through",
    }
}

typst_string_enum! {
    /// A finite CeTZ box or content anchor.
    PyAnchor, "Anchor", Anchor {
        Center => "center",
        North => "north",
        NorthEast => "north-east",
        East => "east",
        SouthEast => "south-east",
        South => "south",
        SouthWest => "south-west",
        West => "west",
        NorthWest => "north-west",
    }
}

typst_string_enum! {
    /// Kurvst path decoration pattern.
    PyPattern, "Pattern", Pattern {
        Wave => "wave",
        Zigzag => "zigzag",
        Coil => "coil",
        Normal => "normal",
        Curve => "curve",
    }
}

typst_string_enum! {
    /// How fixed and relative edge-length limits are combined.
    PyEdgeLengthResolution, "EdgeLengthResolution", EdgeLengthResolution {
        Min => "min",
        Shorter => "shorter",
        Max => "max",
        Longer => "longer",
        Length => "length",
        Fixed => "fixed",
        Ratio => "ratio",
        Relative => "relative",
        Disabled => "none",
        Full => "full",
    }
}

/// Finite debug detail accepted by Linnest's drawing layer.
#[gen_stub_pyclass_enum]
#[pyclass(from_py_object, eq, eq_int, name = "DebugLevel")]
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case")]
enum PyDebugLevel {
    Off,
    Canvas,
    EdgePositions,
}

impl PyDebugLevel {
    const fn level(self) -> u8 {
        match self {
            Self::Off => 0,
            Self::Canvas => 1,
            Self::EdgePositions => 2,
        }
    }
}

impl From<PyDebugLevel> for NativeEnum {
    fn from(value: PyDebugLevel) -> Self {
        Self::DebugLevel(value)
    }
}

typst_string_enum! {
    /// Endpoint shape for a Typst stroke.
    PyStrokeCap, "StrokeCap", StrokeCap {
        Butt => "butt",
        Round => "round",
        Square => "square",
    }
}

typst_string_enum! {
    /// Join shape for a Typst stroke.
    PyStrokeJoin, "StrokeJoin", StrokeJoin {
        Miter => "miter",
        Round => "round",
        Bevel => "bevel",
    }
}

typst_string_enum! {
    /// Finite Typst text slant.
    PyTextStyle, "TextStyle", TextStyle {
        Normal => "normal",
        Italic => "italic",
        Oblique => "oblique",
    }
}

typst_string_enum! {
    /// Built-in Typst stroke dash pattern.
    PyDashPattern, "DashPattern", DashPattern {
        Solid => "solid",
        Dotted => "dotted",
        DenselyDotted => "densely-dotted",
        LooselyDotted => "loosely-dotted",
        Dashed => "dashed",
        DenselyDashed => "densely-dashed",
        LooselyDashed => "loosely-dashed",
        DashDotted => "dash-dotted",
        DenselyDashDotted => "densely-dash-dotted",
        LooselyDashDotted => "loosely-dash-dotted",
    }
}

typst_string_enum! {
    /// Built-in arrowhead symbol used by CeTZ edge marks.
    PyMarkSymbol, "MarkSymbol", MarkSymbol {
        Barbed => "barbed",
        Straight => "straight",
    }
}

typst_string_enum! {
    /// Location of a mark on a rendered edge path.
    PyMarkPosition, "MarkPosition", MarkPosition {
        End => "end",
        Center => "center",
        CenterIfDangling => "center-if-dangling",
    }
}

typst_string_enum! {
    /// Frame used to orient an edge mark.
    PyMarkOrientation, "MarkOrientation", MarkOrientation {
        Path => "path",
        Edge => "edge",
    }
}

typst_string_enum! {
    /// Direction of a mark along an edge path.
    PyMarkDirection, "MarkDirection", MarkDirection {
        Forward => "forward",
        Backward => "backward",
    }
}

/// A Typst length such as `2pt` or `1.2em`.
#[gen_stub_pyclass]
#[pyclass(from_py_object, frozen, name = "Length")]
#[derive(Clone, Debug, PartialEq)]
struct PyLength {
    value: f64,
    unit: LengthUnit,
}

#[gen_stub_pymethods]
#[pymethods]
impl PyLength {
    #[new]
    fn new(
        value: f64,
        #[gen_stub(override_type(type_repr="typing.Literal['pt', 'mm', 'cm', 'in', 'em']", imports=("typing")))]
        unit: &str,
    ) -> PyResult<Self> {
        Ok(Self {
            value: finite(value, "length")?,
            unit: LengthUnit::parse(unit)?,
        })
    }

    #[staticmethod]
    fn pt(value: f64) -> PyResult<Self> {
        Self::new(value, "pt")
    }

    #[staticmethod]
    fn mm(value: f64) -> PyResult<Self> {
        Self::new(value, "mm")
    }

    #[staticmethod]
    fn cm(value: f64) -> PyResult<Self> {
        Self::new(value, "cm")
    }

    #[staticmethod]
    fn inches(value: f64) -> PyResult<Self> {
        Self::new(value, "in")
    }

    #[staticmethod]
    fn em(value: f64) -> PyResult<Self> {
        Self::new(value, "em")
    }

    #[getter]
    fn value(&self) -> f64 {
        self.value
    }

    #[getter]
    fn unit(&self) -> &'static str {
        self.unit.name()
    }

    fn __repr__(&self) -> String {
        format!("Length({}, {:?})", self.value, self.unit.name())
    }
}

/// A Typst ratio expressed in percent.
#[gen_stub_pyclass]
#[pyclass(from_py_object, frozen, name = "Ratio")]
#[derive(Clone, Copy, Debug, PartialEq)]
struct PyRatio {
    percent: f64,
}

#[gen_stub_pymethods]
#[pymethods]
impl PyRatio {
    #[new]
    fn new(percent: f64) -> PyResult<Self> {
        Ok(Self {
            percent: finite(percent, "ratio")?,
        })
    }

    #[staticmethod]
    fn from_fraction(value: f64) -> PyResult<Self> {
        Self::new(finite(value, "ratio fraction")? * 100.0)
    }

    #[getter]
    fn percent(&self) -> f64 {
        self.percent
    }

    fn __repr__(&self) -> String {
        format!("Ratio({})", self.percent)
    }
}

/// A sum of a Typst ratio and length.
#[gen_stub_pyclass]
#[pyclass(from_py_object, frozen, name = "RelativeLength")]
#[derive(Clone, Debug, PartialEq)]
struct PyRelativeLength {
    ratio: Option<f64>,
    length: Option<(f64, LengthUnit)>,
}

#[gen_stub_pymethods]
#[pymethods]
impl PyRelativeLength {
    #[new]
    #[pyo3(signature = (ratio=None, length=None))]
    fn new(ratio: Option<&PyRatio>, length: Option<&PyLength>) -> PyResult<Self> {
        if ratio.is_none() && length.is_none() {
            return Err(PyValueError::new_err(
                "RelativeLength needs a ratio, a length, or both",
            ));
        }
        Ok(Self {
            ratio: ratio.map(|value| value.percent),
            length: length.map(|value| (value.value, value.unit.clone())),
        })
    }

    fn __repr__(&self) -> String {
        format!(
            "RelativeLength(ratio={:?}, length={:?})",
            self.ratio,
            self.length
                .as_ref()
                .map(|(value, unit)| format!("{}{}", value, unit.name()))
        )
    }
}

/// A Typst angle in degrees or radians.
#[gen_stub_pyclass]
#[pyclass(from_py_object, frozen, name = "Angle")]
#[derive(Clone, Debug, PartialEq)]
struct PyAngle {
    value: f64,
    unit: String,
}

#[gen_stub_pymethods]
#[pymethods]
impl PyAngle {
    #[new]
    fn new(
        value: f64,
        #[gen_stub(override_type(type_repr="typing.Literal['deg', 'rad']", imports=("typing")))]
        unit: &str,
    ) -> PyResult<Self> {
        if !matches!(unit, "deg" | "rad") {
            return Err(PyValueError::new_err("angle unit must be 'deg' or 'rad'"));
        }
        Ok(Self {
            value: finite(value, "angle")?,
            unit: unit.to_owned(),
        })
    }

    #[staticmethod]
    fn degrees(value: f64) -> PyResult<Self> {
        Self::new(value, "deg")
    }

    #[staticmethod]
    fn radians(value: f64) -> PyResult<Self> {
        Self::new(value, "rad")
    }

    fn __repr__(&self) -> String {
        format!("Angle({}, {:?})", self.value, self.unit)
    }
}

/// A Typst fractional track size such as `1fr`.
#[gen_stub_pyclass]
#[pyclass(from_py_object, frozen, name = "Fraction")]
#[derive(Clone, Copy, Debug, PartialEq)]
struct PyFraction {
    value: f64,
}

#[gen_stub_pymethods]
#[pymethods]
impl PyFraction {
    #[new]
    fn new(value: f64) -> PyResult<Self> {
        Ok(Self {
            value: non_negative(value, "fraction")?,
        })
    }

    #[getter]
    fn value(&self) -> f64 {
        self.value
    }

    fn __repr__(&self) -> String {
        format!("Fraction({})", self.value)
    }
}

/// A safe Typst color value.
#[gen_stub_pyclass]
#[pyclass(from_py_object, frozen, name = "Color")]
#[derive(Clone, Debug, PartialEq)]
struct PyColor {
    value: ColorValue,
}

#[gen_stub_pymethods]
#[pymethods]
impl PyColor {
    #[new]
    fn new(value: &str) -> PyResult<Self> {
        const NAMED: &[&str] = &[
            "black", "gray", "silver", "white", "navy", "blue", "aqua", "teal", "eastern",
            "purple", "fuchsia", "maroon", "red", "orange", "yellow", "olive", "green",
        ];
        let value = if let Some(digits) = value.strip_prefix('#') {
            if !matches!(digits.len(), 3 | 4 | 6 | 8)
                || !digits.chars().all(|c| c.is_ascii_hexdigit())
            {
                return Err(PyValueError::new_err(format!(
                    "invalid hexadecimal color {value:?}"
                )));
            }
            ColorValue::Hex(value.to_ascii_lowercase())
        } else if NAMED.contains(&value) {
            ColorValue::Named(value.to_owned())
        } else {
            return Err(PyValueError::new_err(format!(
                "unknown named color {value:?}; use a hexadecimal color or Color.rgb"
            )));
        };
        Ok(Self { value })
    }

    #[staticmethod]
    #[pyo3(signature = (red, green, blue, alpha=None))]
    fn rgb(red: u8, green: u8, blue: u8, alpha: Option<u8>) -> Self {
        Self {
            value: ColorValue::Rgba(red, green, blue, alpha),
        }
    }

    fn __repr__(&self) -> String {
        match &self.value {
            ColorValue::Named(value) | ColorValue::Hex(value) => {
                format!("Color({value:?})")
            }
            ColorValue::Rgba(red, green, blue, None) => {
                format!("Color.rgb({red}, {green}, {blue})")
            }
            ColorValue::Rgba(red, green, blue, Some(alpha)) => {
                format!("Color.rgb({red}, {green}, {blue}, {alpha})")
            }
        }
    }
}

#[derive(Clone, Copy)]
enum ValueRule {
    Any,
    Bool,
    String,
    NonNegativeInt,
    Number,
    NonNegativeNumber,
    Point,
    Angle,
    Enum(EnumKind),
    EnumOrFunction(EnumKind),
    NodeIndices,
    NodeGroups,
    HedgeSelection,
    DrawSubgraphs,
    Dictionary,
    Unit,
    Length,
    Paint,
    Stroke,
    Dash,
    Content,
    ContentOrFunction,
    Function,
    Style,
    StyleLayers,
    Radius,
    Padding,
    PatternOrStyleLayers,
}

impl ValueRule {
    fn accepts(self, value: &NativeValue) -> bool {
        let expression_kind_of = |value: &NativeValue, expected| match value {
            NativeValue::Expression(TypstExpression::Symbol(symbol)) => symbol.kind == expected,
            NativeValue::Expression(TypstExpression::Bind { .. }) => {
                expected == SymbolKind::Function
            }
            NativeValue::Expression(TypstExpression::Call { .. }) => {
                expected != SymbolKind::Function
            }
            _ => false,
        };
        let expression_kind = |expected| expression_kind_of(value, expected);
        match self {
            Self::Any => true,
            Self::Bool => {
                matches!(value, NativeValue::Bool(_)) || expression_kind(SymbolKind::Value)
            }
            Self::String => {
                matches!(value, NativeValue::String(_)) || expression_kind(SymbolKind::Value)
            }
            Self::NonNegativeInt => {
                matches!(value, NativeValue::Int(number) if *number >= 0)
                    || expression_kind(SymbolKind::Value)
            }
            Self::Number => {
                matches!(value, NativeValue::Int(_) | NativeValue::Float(_))
                    || expression_kind(SymbolKind::Value)
            }
            Self::NonNegativeNumber => {
                matches!(value, NativeValue::Int(number) if *number >= 0)
                    || matches!(value, NativeValue::Float(number) if *number >= 0.0)
                    || expression_kind(SymbolKind::Value)
            }
            Self::Point => {
                matches!(value, NativeValue::Array(values) if values.len() == 2 && values.iter().all(|value| Self::Number.accepts(value)))
                    || matches!(value, NativeValue::Dict(values) if {
                        values.len() == 2
                            && values.get("x").is_some_and(|value| Self::Number.accepts(value))
                            && values.get("y").is_some_and(|value| Self::Number.accepts(value))
                    })
                    || expression_kind(SymbolKind::Value)
            }
            Self::Angle => {
                matches!(
                    value,
                    NativeValue::Int(_) | NativeValue::Float(_) | NativeValue::Angle(..)
                ) || expression_kind(SymbolKind::Value)
            }
            Self::Enum(expected) => {
                matches!(value, NativeValue::Enum(value) if value.kind() == expected)
                    || expression_kind(SymbolKind::Value)
            }
            Self::EnumOrFunction(expected) => {
                matches!(value, NativeValue::Enum(value) if value.kind() == expected)
                    || expression_kind(SymbolKind::Value)
                    || expression_kind(SymbolKind::Function)
            }
            Self::NodeIndices => {
                matches!(value, NativeValue::Array(values) if values.iter().all(|value| matches!(value, NativeValue::Int(index) if *index >= 0)))
                    || expression_kind(SymbolKind::Value)
            }
            Self::NodeGroups => {
                matches!(value, NativeValue::Array(groups) if groups.iter().all(|group| {
                    matches!(group, NativeValue::Array(values) if values.iter().all(|value| matches!(value, NativeValue::Int(index) if *index >= 0)))
                        || expression_kind_of(group, SymbolKind::Value)
                        || expression_kind_of(group, SymbolKind::Function)
                })) || expression_kind(SymbolKind::Value)
            }
            Self::HedgeSelection => {
                matches!(value, NativeValue::Array(values) if values.iter().all(|value| matches!(value, NativeValue::Bool(_))))
                    || expression_kind(SymbolKind::Value)
                    || expression_kind(SymbolKind::Function)
            }
            Self::DrawSubgraphs => {
                matches!(value, NativeValue::Array(values) if {
                    values.iter().all(|value| matches!(value, NativeValue::Bool(_)))
                        || values.iter().all(|value| {
                            matches!(value, NativeValue::Array(bits) if bits.iter().all(|bit| matches!(bit, NativeValue::Bool(_))))
                                || expression_kind_of(value, SymbolKind::Value)
                                || expression_kind_of(value, SymbolKind::Function)
                                || matches!(value, NativeValue::Dict(record) if {
                                    record.keys().all(|key| matches!(key.as_str(), "subgraph" | "edge-style"))
                                        && record.get("subgraph").is_some_and(|selection| {
                                            matches!(selection, NativeValue::Array(bits) if bits.iter().all(|bit| matches!(bit, NativeValue::Bool(_))))
                                                || expression_kind_of(selection, SymbolKind::Value)
                                                || expression_kind_of(selection, SymbolKind::Function)
                                        })
                                        && record.get("edge-style").is_none_or(|style| {
                                            matches!(style, NativeValue::None)
                                                || Self::Style.accepts(style)
                                        })
                                })
                        })
                }) || expression_kind(SymbolKind::Value)
                    || expression_kind(SymbolKind::Function)
            }
            Self::Dictionary => {
                matches!(value, NativeValue::Dict(_)) || expression_kind(SymbolKind::Value)
            }
            Self::Unit => {
                matches!(
                    value,
                    NativeValue::Int(_)
                        | NativeValue::Float(_)
                        | NativeValue::Length(..)
                        | NativeValue::Ratio(_)
                        | NativeValue::RelativeLength { .. }
                ) || expression_kind(SymbolKind::Value)
            }
            Self::Length => {
                matches!(
                    value,
                    NativeValue::Int(_)
                        | NativeValue::Float(_)
                        | NativeValue::Length(..)
                        | NativeValue::Ratio(_)
                        | NativeValue::RelativeLength { .. }
                ) || expression_kind(SymbolKind::Value)
            }
            Self::Paint => {
                matches!(value, NativeValue::Color(_)) || expression_kind(SymbolKind::Value)
            }
            Self::Stroke => {
                matches!(
                    value,
                    NativeValue::Stroke(_)
                        | NativeValue::Color(_)
                        | NativeValue::Length(..)
                        | NativeValue::Int(_)
                        | NativeValue::Float(_)
                ) || expression_kind(SymbolKind::Value)
            }
            Self::Dash => {
                matches!(value, NativeValue::Dash(_)) || expression_kind(SymbolKind::Value)
            }
            Self::Content => {
                matches!(
                    value,
                    NativeValue::String(_) | NativeValue::Text(_) | NativeValue::Math(_)
                ) || expression_kind(SymbolKind::Content)
            }
            Self::ContentOrFunction => {
                Self::Content.accepts(value) || expression_kind(SymbolKind::Function)
            }
            Self::Function => expression_kind(SymbolKind::Function),
            Self::Style => {
                matches!(value, NativeValue::Dict(_))
                    || expression_kind(SymbolKind::Value)
                    || expression_kind(SymbolKind::Function)
            }
            Self::StyleLayers => {
                matches!(value, NativeValue::Dict(_))
                    || matches!(value, NativeValue::Array(values) if values.iter().all(|value| matches!(value, NativeValue::Dict(_))))
                    || expression_kind(SymbolKind::Value)
                    || expression_kind(SymbolKind::Function)
            }
            Self::Radius => {
                matches!(value, NativeValue::Int(_) | NativeValue::Float(_))
                    || matches!(value, NativeValue::Array(values) if values.iter().all(|value| matches!(value, NativeValue::Int(_) | NativeValue::Float(_))))
                    || expression_kind(SymbolKind::Value)
            }
            Self::Padding => {
                matches!(
                    value,
                    NativeValue::Int(_)
                        | NativeValue::Float(_)
                        | NativeValue::Array(_)
                        | NativeValue::Dict(_)
                        | NativeValue::Insets(_)
                ) || expression_kind(SymbolKind::Value)
            }
            Self::PatternOrStyleLayers => {
                matches!(value, NativeValue::Enum(value) if value.kind() == EnumKind::Pattern)
                    || Self::StyleLayers.accepts(value)
            }
        }
    }

    fn name(self) -> &'static str {
        match self {
            Self::Any => "a supported native Typst value",
            Self::Bool => "a bool",
            Self::String => "a string or module value",
            Self::NonNegativeInt => "a non-negative int",
            Self::Number => "a number",
            Self::NonNegativeNumber => "a non-negative number",
            Self::Point => "a two-number point dictionary or array, or module value",
            Self::Angle => "an Angle, number, or module value",
            Self::Enum(kind) => kind.name(),
            Self::EnumOrFunction(kind) => match kind {
                EnumKind::EdgeLengthResolution => "EdgeLengthResolution or a module function",
                _ => "an enum value or module function",
            },
            Self::NodeIndices => "an array of non-negative node indices",
            Self::NodeGroups => "an array of non-negative node-index arrays or module functions",
            Self::HedgeSelection => "a boolean half-edge selection array or module function",
            Self::DrawSubgraphs => {
                "a boolean half-edge selection, module function, or array of typed selections"
            }
            Self::Dictionary => "a dictionary",
            Self::Unit => "a number, length, ratio, or relative length",
            Self::Length => "a number, length, ratio, or relative length",
            Self::Paint => "a Color or module value",
            Self::Stroke => "a Stroke, color, thickness, or module value",
            Self::Dash => "a Dash or module value",
            Self::Content => "text, math, or module content",
            Self::ContentOrFunction => "text, math, module content, or a label callback",
            Self::Function => "a module function or bound module function",
            Self::Style => "a style dictionary, module value, or module function",
            Self::StyleLayers => {
                "a style dictionary, array of style dictionaries, or module function"
            }
            Self::Radius => "a number, numeric array, or module value",
            Self::Padding => "a number, array, dictionary, Insets, or module value",
            Self::PatternOrStyleLayers => {
                "a Pattern, style dictionary, style layers, or module function"
            }
        }
    }
}

#[derive(Clone, Copy)]
struct FieldSpec {
    python: &'static str,
    typst: &'static str,
    rule: ValueRule,
    allow_none: bool,
    allow_auto: bool,
}

impl FieldSpec {
    const fn new(python: &'static str, typst: &'static str, rule: ValueRule) -> Self {
        Self {
            python,
            typst,
            rule,
            allow_none: false,
            allow_auto: false,
        }
    }

    const fn none(mut self) -> Self {
        self.allow_none = true;
        self
    }

    const fn auto(mut self) -> Self {
        self.allow_auto = true;
        self
    }

    fn validate(self, value: &NativeValue) -> PyResult<()> {
        if matches!(value, NativeValue::Inherit)
            || (self.allow_none && matches!(value, NativeValue::None))
            || (self.allow_auto && matches!(value, NativeValue::Auto))
            || self.rule.accepts(value)
        {
            Ok(())
        } else {
            Err(PyTypeError::new_err(format!(
                "{} must be {}, not {}",
                self.python,
                self.rule.name(),
                value.kind()
            )))
        }
    }
}

fn parse_options(
    kwargs: Option<&Bound<'_, PyDict>>,
    fields: &[FieldSpec],
) -> PyResult<BTreeMap<String, NativeValue>> {
    let mut out = fields
        .iter()
        .map(|field| (field.typst.to_owned(), NativeValue::Inherit))
        .collect::<BTreeMap<_, _>>();
    let Some(kwargs) = kwargs else {
        return Ok(out);
    };
    for (key, value) in kwargs.iter() {
        let key = key
            .extract::<String>()
            .map_err(|_| PyTypeError::new_err("Typst option names must be Python strings"))?;
        let field = fields
            .iter()
            .find(|field| field.python == key)
            .ok_or_else(|| PyTypeError::new_err(format!("unknown Typst option {key:?}")))?;
        let value = native_from_py(&value, 0)?;
        field.validate(&value)?;
        if !matches!(value, NativeValue::Inherit) {
            out.insert(field.typst.to_owned(), value);
        }
    }
    Ok(out)
}

fn options_repr(name: &str, values: &[&BTreeMap<String, NativeValue>]) -> String {
    let count = values
        .iter()
        .flat_map(|value| value.values())
        .filter(|value| !matches!(value, NativeValue::Inherit))
        .count();
    format!("{name}(fields={count})")
}

const STROKE_FIELDS: &[FieldSpec] = &[
    FieldSpec::new("paint", "paint", ValueRule::Paint).none(),
    FieldSpec::new("thickness", "thickness", ValueRule::Length).none(),
    FieldSpec::new("cap", "cap", ValueRule::Enum(EnumKind::StrokeCap)).none(),
    FieldSpec::new("join", "join", ValueRule::Enum(EnumKind::StrokeJoin)).none(),
    FieldSpec::new("dash", "dash", ValueRule::Dash).none(),
    FieldSpec::new("miter_limit", "miter-limit", ValueRule::NonNegativeNumber),
];

/// A typed Typst stroke dictionary.
#[gen_stub_pyclass]
#[pyclass(from_py_object, frozen, name = "Stroke")]
#[derive(Clone, Debug, PartialEq)]
struct PyStroke {
    values: BTreeMap<String, NativeValue>,
}

#[gen_stub_pymethods]
#[pymethods]
impl PyStroke {
    #[new]
    #[pyo3(signature = (**kwargs))]
    #[gen_stub(skip)]
    fn new(kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<Self> {
        Ok(Self {
            values: parse_options(kwargs, STROKE_FIELDS)?,
        })
    }

    fn __repr__(&self) -> String {
        options_repr("Stroke", &[&self.values])
    }
}

/// A named or explicit Typst dash pattern.
#[gen_stub_pyclass]
#[pyclass(from_py_object, frozen, name = "Dash")]
#[derive(Clone, Debug, PartialEq)]
struct PyDash {
    value: DashValue,
}

#[gen_stub_pymethods]
#[pymethods]
impl PyDash {
    #[new]
    fn new(pattern: PyDashPattern) -> Self {
        Self {
            value: DashValue::Named(pattern.typst_name().to_owned()),
        }
    }

    #[staticmethod]
    #[pyo3(signature = (values, *, phase=None))]
    fn pattern(
        #[gen_stub(override_type(type_repr = "_LengthArray"))] values: &Bound<'_, PyAny>,
        #[gen_stub(override_type(type_repr = "_LengthValue | None"))] phase: Option<
            &Bound<'_, PyAny>,
        >,
    ) -> PyResult<Self> {
        let NativeValue::Array(array) = native_from_py(values, 0)? else {
            return Err(PyTypeError::new_err(
                "Dash.pattern values must be a list or tuple",
            ));
        };
        if array.is_empty() {
            return Err(PyValueError::new_err(
                "Dash.pattern needs at least one dash or gap",
            ));
        }
        for value in &array {
            if !ValueRule::Length.accepts(value) {
                return Err(PyTypeError::new_err(
                    "Dash.pattern entries must be numbers or lengths",
                ));
            }
        }
        let phase = phase
            .map(|value| native_from_py(value, 0).map(Box::new))
            .transpose()?;
        if let Some(phase) = &phase {
            if !ValueRule::Length.accepts(phase) {
                return Err(PyTypeError::new_err(
                    "Dash.pattern phase must be a number or length",
                ));
            }
        }
        Ok(Self {
            value: DashValue::Pattern { array, phase },
        })
    }

    fn __repr__(&self) -> String {
        match &self.value {
            DashValue::Named(name) => format!("Dash({name:?})"),
            DashValue::Pattern { array, .. } => format!("Dash.pattern(len={})", array.len()),
        }
    }
}

const INSET_FIELDS: &[FieldSpec] = &[
    FieldSpec::new("all", "all", ValueRule::Length),
    FieldSpec::new("x", "x", ValueRule::Length),
    FieldSpec::new("y", "y", ValueRule::Length),
    FieldSpec::new("left", "left", ValueRule::Length),
    FieldSpec::new("right", "right", ValueRule::Length),
    FieldSpec::new("top", "top", ValueRule::Length),
    FieldSpec::new("bottom", "bottom", ValueRule::Length),
];

/// Typed CeTZ/Typst inset values.
#[gen_stub_pyclass]
#[pyclass(from_py_object, frozen, name = "Insets")]
#[derive(Clone, Debug, PartialEq)]
struct PyInsets {
    values: BTreeMap<String, NativeValue>,
}

#[gen_stub_pymethods]
#[pymethods]
impl PyInsets {
    #[new]
    #[pyo3(signature = (**kwargs))]
    #[gen_stub(skip)]
    fn new(kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<Self> {
        let values = parse_options(kwargs, INSET_FIELDS)?;
        let explicit = values
            .iter()
            .filter(|(_, value)| !matches!(value, NativeValue::Inherit))
            .collect::<Vec<_>>();
        if explicit.is_empty() {
            return Err(PyValueError::new_err(
                "Insets needs at least one inset value",
            ));
        }
        if explicit.iter().any(|(name, _)| name.as_str() == "all") && explicit.len() > 1 {
            return Err(PyValueError::new_err(
                "Insets(all=...) cannot be combined with directional values",
            ));
        }
        Ok(Self { values })
    }

    fn __repr__(&self) -> String {
        options_repr("Insets", &[&self.values])
    }
}

const MARK_FIELDS: &[FieldSpec] = &[
    FieldSpec::new("start", "start", ValueRule::Any).none(),
    FieldSpec::new("end", "end", ValueRule::Any).none(),
    FieldSpec::new("fill", "fill", ValueRule::Paint).none(),
    FieldSpec::new("stroke", "stroke", ValueRule::Stroke).none(),
    FieldSpec::new("scale", "scale", ValueRule::NonNegativeNumber),
    FieldSpec::new("anchor", "anchor", ValueRule::Enum(EnumKind::Anchor)),
    FieldSpec::new("shorten_to", "shorten-to", ValueRule::Length)
        .none()
        .auto(),
];

/// Typed CeTZ mark configuration.
#[gen_stub_pyclass]
#[pyclass(from_py_object, frozen, name = "Mark")]
#[derive(Clone, Debug, PartialEq)]
struct PyMark {
    values: BTreeMap<String, NativeValue>,
}

#[gen_stub_pymethods]
#[pymethods]
impl PyMark {
    #[new]
    #[pyo3(signature = (**kwargs))]
    #[gen_stub(skip)]
    fn new(kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<Self> {
        let values = parse_options(kwargs, MARK_FIELDS)?;
        if values
            .values()
            .all(|value| matches!(value, NativeValue::Inherit))
        {
            return Err(PyValueError::new_err("Mark needs at least one option"));
        }
        Ok(Self { values })
    }

    #[staticmethod]
    fn barbed() -> Self {
        Self {
            values: BTreeMap::from([(
                "end".to_owned(),
                NativeValue::Enum(PyMarkSymbol::Barbed.into()),
            )]),
        }
    }

    #[staticmethod]
    fn straight() -> Self {
        Self {
            values: BTreeMap::from([(
                "end".to_owned(),
                NativeValue::Enum(PyMarkSymbol::Straight.into()),
            )]),
        }
    }

    fn __repr__(&self) -> String {
        options_repr("Mark", &[&self.values])
    }
}

const TEXT_FIELDS: &[FieldSpec] = &[
    FieldSpec::new("fill", "fill", ValueRule::Paint),
    FieldSpec::new("size", "size", ValueRule::Length),
    FieldSpec::new("weight", "weight", ValueRule::Any),
    FieldSpec::new("style", "style", ValueRule::Enum(EnumKind::TextStyle)),
];

/// Literal text content with optional typed text styling.
#[gen_stub_pyclass]
#[pyclass(from_py_object, frozen, name = "TextLabel")]
#[derive(Clone, Debug, PartialEq)]
struct PyTextLabel {
    value: TextValue,
}

#[gen_stub_pymethods]
#[pymethods]
impl PyTextLabel {
    #[new]
    #[pyo3(signature = (text, **kwargs))]
    #[gen_stub(skip)]
    fn new(text: String, kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<Self> {
        Ok(Self {
            value: TextValue {
                text,
                options: parse_options(kwargs, TEXT_FIELDS)?,
            },
        })
    }

    #[getter]
    fn text(&self) -> &str {
        &self.value.text
    }

    fn __repr__(&self) -> String {
        format!("TextLabel({:?})", self.value.text)
    }
}

fn math_token(value: &str, what: &str) -> PyResult<String> {
    if value.is_empty()
        || !value
            .chars()
            .all(|c| c == '_' || c == '-' || c.is_alphanumeric())
    {
        Err(PyValueError::new_err(format!(
            "{what} must contain only letters, digits, '_' or '-'"
        )))
    } else {
        Ok(value.to_owned())
    }
}

fn math_script(value: &Bound<'_, PyAny>, what: &str) -> PyResult<MathScript> {
    if !value.is_instance_of::<PyBool>() {
        if let Ok(value) = value.extract::<i64>() {
            return Ok(MathScript::Integer(value));
        }
    }
    let value = value
        .extract::<String>()
        .map_err(|_| PyTypeError::new_err(format!("{what} must be a symbol string or integer")))?;
    Ok(MathScript::Symbol(math_token(&value, what)?))
}

/// Safe mathematical identifier content, optionally with scripts.
#[gen_stub_pyclass]
#[pyclass(from_py_object, frozen, name = "MathSymbol")]
#[derive(Clone, Debug, PartialEq)]
struct PyMathSymbol {
    value: MathValue,
}

#[gen_stub_pymethods]
#[pymethods]
impl PyMathSymbol {
    #[new]
    #[pyo3(signature = (symbol, *, subscript=None, superscript=None))]
    fn new(
        symbol: &str,
        #[gen_stub(override_type(type_repr="builtins.str | builtins.int | None", imports=("builtins")))]
        subscript: Option<&Bound<'_, PyAny>>,
        #[gen_stub(override_type(type_repr="builtins.str | builtins.int | None", imports=("builtins")))]
        superscript: Option<&Bound<'_, PyAny>>,
    ) -> PyResult<Self> {
        Ok(Self {
            value: MathValue {
                symbol: math_token(symbol, "math symbol")?,
                subscript: subscript
                    .map(|value| math_script(value, "math subscript"))
                    .transpose()?,
                superscript: superscript
                    .map(|value| math_script(value, "math superscript"))
                    .transpose()?,
            },
        })
    }

    fn __repr__(&self) -> String {
        format!(
            "MathSymbol({:?}, subscript={:?}, superscript={:?})",
            self.value.symbol, self.value.subscript, self.value.superscript
        )
    }
}

/// A local or package Typst module whose exports can be referenced safely.
#[gen_stub_pyclass]
#[pyclass(from_py_object, frozen, name = "TypstModule")]
#[derive(Clone, Debug, PartialEq, Eq)]
struct PyTypstModule {
    source: TypstModuleSource,
}

#[gen_stub_pymethods]
#[pymethods]
impl PyTypstModule {
    #[new]
    fn new(
        #[gen_stub(override_type(type_repr="builtins.str | os.PathLike[builtins.str]", imports=("builtins", "os")))]
        path: PathBuf,
    ) -> Self {
        Self {
            source: TypstModuleSource::File(path),
        }
    }

    #[staticmethod]
    fn file(
        #[gen_stub(override_type(type_repr="builtins.str | os.PathLike[builtins.str]", imports=("builtins", "os")))]
        path: PathBuf,
    ) -> Self {
        Self::new(path)
    }

    #[staticmethod]
    fn package(specification: String) -> PyResult<Self> {
        clinnet::TypstModule::package("_linnet_validate", specification.clone()).map_err(|_| {
            PyValueError::new_err("Typst package must look like '@namespace/name:version'")
        })?;
        Ok(Self {
            source: TypstModuleSource::Package(specification),
        })
    }

    #[gen_stub(override_return_type(type_repr = "_TypstValueRef"))]
    fn value(&self, export: &str) -> PyResult<PyTypstRef> {
        self.reference(export, SymbolKind::Value)
    }

    #[gen_stub(override_return_type(type_repr = "_TypstContentRef"))]
    fn content(&self, export: &str) -> PyResult<PyTypstRef> {
        self.reference(export, SymbolKind::Content)
    }

    #[gen_stub(override_return_type(type_repr = "_TypstFunctionRef"))]
    fn function(&self, export: &str) -> PyResult<PyTypstRef> {
        self.reference(export, SymbolKind::Function)
    }

    fn __repr__(&self) -> String {
        match &self.source {
            TypstModuleSource::File(path) => format!("TypstModule({path:?})"),
            TypstModuleSource::Package(specification) => {
                format!("TypstModule.package({specification:?})")
            }
        }
    }
}

impl PyTypstModule {
    fn reference(&self, export: &str, kind: SymbolKind) -> PyResult<PyTypstRef> {
        let path = export
            .split('.')
            .map(|part| identifier(part, "module export"))
            .collect::<PyResult<Vec<_>>>()?;
        Ok(PyTypstRef {
            expression: TypstExpression::Symbol(ModuleSymbol {
                module: self.source.clone(),
                path,
                kind,
            }),
        })
    }
}

fn call_arguments(
    args: &Bound<'_, PyTuple>,
    kwargs: Option<&Bound<'_, PyDict>>,
) -> PyResult<CallArguments> {
    let positional = args
        .iter()
        .map(|value| native_from_py(&value, 0))
        .collect::<PyResult<Vec<_>>>()?;
    if positional
        .iter()
        .any(|value| matches!(value, NativeValue::Inherit))
    {
        return Err(PyTypeError::new_err(
            "INHERIT is only valid as an option value",
        ));
    }
    let mut named = BTreeMap::new();
    if let Some(kwargs) = kwargs {
        for (name, value) in kwargs.iter() {
            let name = option_name(&name.extract::<String>()?)?;
            let value = native_from_py(&value, 0)?;
            if !matches!(value, NativeValue::Inherit) {
                named.insert(name, value);
            }
        }
    }
    Ok(CallArguments { positional, named })
}

/// A typed export from a Typst module.
#[gen_stub_pyclass]
#[pyclass(from_py_object, frozen, name = "TypstRef")]
#[derive(Clone, Debug, PartialEq)]
struct PyTypstRef {
    expression: TypstExpression,
}

#[gen_stub_pymethods]
#[pymethods]
impl PyTypstRef {
    #[pyo3(signature = (*args, **kwargs))]
    fn call(
        &self,
        #[gen_stub(override_type(type_repr = "_NativeValue"))] args: &Bound<'_, PyTuple>,
        #[gen_stub(override_type(type_repr = "_NativeValue"))] kwargs: Option<&Bound<'_, PyDict>>,
    ) -> PyResult<PyTypstCall> {
        let TypstExpression::Symbol(symbol) = &self.expression else {
            unreachable!("TypstRef always contains a symbol")
        };
        if symbol.kind != SymbolKind::Function {
            return Err(PyTypeError::new_err(format!(
                "cannot call a Typst {} reference",
                symbol.kind.name()
            )));
        }
        Ok(PyTypstCall {
            expression: TypstExpression::Call {
                callee: Box::new(self.expression.clone()),
                arguments: call_arguments(args, kwargs)?,
            },
        })
    }

    #[pyo3(signature = (*args, **kwargs))]
    fn bind(
        &self,
        #[gen_stub(override_type(type_repr = "_NativeValue"))] args: &Bound<'_, PyTuple>,
        #[gen_stub(override_type(type_repr = "_NativeValue"))] kwargs: Option<&Bound<'_, PyDict>>,
    ) -> PyResult<PyTypstBind> {
        let TypstExpression::Symbol(symbol) = &self.expression else {
            unreachable!("TypstRef always contains a symbol")
        };
        if symbol.kind != SymbolKind::Function {
            return Err(PyTypeError::new_err(format!(
                "cannot bind a Typst {} reference",
                symbol.kind.name()
            )));
        }
        Ok(PyTypstBind {
            expression: TypstExpression::Bind {
                function: Box::new(self.expression.clone()),
                arguments: call_arguments(args, kwargs)?,
            },
        })
    }

    fn __repr__(&self) -> String {
        let TypstExpression::Symbol(symbol) = &self.expression else {
            unreachable!("TypstRef always contains a symbol")
        };
        format!(
            "TypstRef({:?}, {})",
            symbol.path.join("."),
            symbol.kind.name()
        )
    }
}

/// A call to an explicitly imported Typst function.
#[gen_stub_pyclass]
#[pyclass(from_py_object, frozen, name = "TypstCall")]
#[derive(Clone, Debug, PartialEq)]
struct PyTypstCall {
    expression: TypstExpression,
}

#[gen_stub_pymethods]
#[pymethods]
impl PyTypstCall {
    fn __repr__(&self) -> &'static str {
        "TypstCall(...)"
    }
}

/// A Typst function partially applied through its native `.with` method.
#[gen_stub_pyclass]
#[pyclass(from_py_object, frozen, name = "TypstBind")]
#[derive(Clone, Debug, PartialEq)]
struct PyTypstBind {
    expression: TypstExpression,
}

#[gen_stub_pymethods]
#[pymethods]
impl PyTypstBind {
    #[pyo3(signature = (*args, **kwargs))]
    fn call(
        &self,
        #[gen_stub(override_type(type_repr = "_NativeValue"))] args: &Bound<'_, PyTuple>,
        #[gen_stub(override_type(type_repr = "_NativeValue"))] kwargs: Option<&Bound<'_, PyDict>>,
    ) -> PyResult<PyTypstCall> {
        Ok(PyTypstCall {
            expression: TypstExpression::Call {
                callee: Box::new(self.expression.clone()),
                arguments: call_arguments(args, kwargs)?,
            },
        })
    }

    fn __repr__(&self) -> &'static str {
        "TypstBind(...)"
    }
}

fn native_from_py(value: &Bound<'_, PyAny>, depth: usize) -> PyResult<NativeValue> {
    if depth > MAX_NATIVE_DEPTH {
        return Err(PyValueError::new_err(format!(
            "native Typst value exceeds maximum nesting depth {MAX_NATIVE_DEPTH}"
        )));
    }
    if value.is_none() {
        return Ok(NativeValue::None);
    }
    if value.extract::<PyRef<'_, PyInherit>>().is_ok() {
        return Ok(NativeValue::Inherit);
    }
    if value.extract::<PyRef<'_, PyAuto>>().is_ok() {
        return Ok(NativeValue::Auto);
    }
    macro_rules! extract_enum {
        ($($enum:ty),+ $(,)?) => {
            $(
                if let Ok(value) = value.extract::<$enum>() {
                    return Ok(NativeValue::Enum(value.into()));
                }
            )+
        };
    }
    extract_enum!(
        PyLayoutAlgorithm,
        PyLayoutDirection,
        PyLabelLayout,
        PyRankAlignment,
        PyLayoutNodes,
        PyPlacement,
        PyCompass,
        PyRouting,
        PyRoutePoints,
        PyAnchor,
        PyPattern,
        PyEdgeLengthResolution,
        PyDebugLevel,
        PyStrokeCap,
        PyStrokeJoin,
        PyTextStyle,
        PyDashPattern,
        PyMarkSymbol,
        PyMarkPosition,
        PyMarkOrientation,
        PyMarkDirection,
    );
    if let Ok(value) = value.extract::<bool>() {
        return Ok(NativeValue::Bool(value));
    }
    if value.cast::<PyInt>().is_ok() {
        return value.extract::<i64>().map(NativeValue::Int).map_err(|_| {
            PyOverflowError::new_err("Typst integers must fit in a signed 64-bit integer")
        });
    }
    if value.cast::<PyFloat>().is_ok() {
        let value = value.extract::<f64>()?;
        return Ok(NativeValue::Float(finite(value, "Typst float")?));
    }
    if let Ok(value) = value.extract::<String>() {
        return Ok(NativeValue::String(value));
    }
    if let Ok(value) = value.extract::<PyRef<'_, PyLength>>() {
        return Ok(NativeValue::Length(value.value, value.unit.clone()));
    }
    if let Ok(value) = value.extract::<PyRef<'_, PyRatio>>() {
        return Ok(NativeValue::Ratio(value.percent));
    }
    if let Ok(value) = value.extract::<PyRef<'_, PyRelativeLength>>() {
        return Ok(NativeValue::RelativeLength {
            ratio: value.ratio,
            length: value.length.clone(),
        });
    }
    if let Ok(value) = value.extract::<PyRef<'_, PyAngle>>() {
        return Ok(NativeValue::Angle(value.value, value.unit.clone()));
    }
    if let Ok(value) = value.extract::<PyRef<'_, PyFraction>>() {
        return Ok(NativeValue::Fraction(value.value));
    }
    if let Ok(value) = value.extract::<PyRef<'_, PyColor>>() {
        return Ok(NativeValue::Color(value.value.clone()));
    }
    if let Ok(value) = value.extract::<PyRef<'_, PyStroke>>() {
        return Ok(NativeValue::Stroke(value.values.clone()));
    }
    if let Ok(value) = value.extract::<PyRef<'_, PyDash>>() {
        return Ok(NativeValue::Dash(value.value.clone()));
    }
    if let Ok(value) = value.extract::<PyRef<'_, PyInsets>>() {
        return Ok(NativeValue::Insets(value.values.clone()));
    }
    if let Ok(value) = value.extract::<PyRef<'_, PyMark>>() {
        return Ok(NativeValue::Mark(value.values.clone()));
    }
    if let Ok(value) = value.extract::<PyRef<'_, PyTextLabel>>() {
        return Ok(NativeValue::Text(value.value.clone()));
    }
    if let Ok(value) = value.extract::<PyRef<'_, PyMathSymbol>>() {
        return Ok(NativeValue::Math(value.value.clone()));
    }
    if let Ok(value) = value.extract::<PyRef<'_, PyTypstRef>>() {
        return Ok(NativeValue::Expression(value.expression.clone()));
    }
    if let Ok(value) = value.extract::<PyRef<'_, PyTypstCall>>() {
        return Ok(NativeValue::Expression(value.expression.clone()));
    }
    if let Ok(value) = value.extract::<PyRef<'_, PyTypstBind>>() {
        return Ok(NativeValue::Expression(value.expression.clone()));
    }
    if let Ok(value) = value.cast::<PyList>() {
        return Ok(NativeValue::Array(
            value
                .iter()
                .map(|value| native_from_py(&value, depth + 1))
                .collect::<PyResult<Vec<_>>>()?,
        ));
    }
    if let Ok(value) = value.cast::<PyTuple>() {
        return Ok(NativeValue::Array(
            value
                .iter()
                .map(|value| native_from_py(&value, depth + 1))
                .collect::<PyResult<Vec<_>>>()?,
        ));
    }
    if let Ok(value) = value.cast::<PyDict>() {
        let mut values = BTreeMap::new();
        for (key, value) in value.iter() {
            let key = key.extract::<String>().map_err(|_| {
                PyTypeError::new_err("native Typst dictionary keys must be strings")
            })?;
            values.insert(key, native_from_py(&value, depth + 1)?);
        }
        return Ok(NativeValue::Dict(values));
    }
    Err(PyTypeError::new_err(format!(
        "unsupported native Typst value {}; Python callables and arbitrary objects cannot cross the Typst process",
        value.get_type().name()?
    )))
}

const GRAPH_STYLE_FIELDS: &[FieldSpec] = &[
    FieldSpec::new("scope", "scope", ValueRule::Dictionary),
    FieldSpec::new("unit", "unit", ValueRule::Unit),
    FieldSpec::new("node_label", "node-label", ValueRule::ContentOrFunction)
        .none()
        .auto(),
    FieldSpec::new("node_label_style", "node-label-style", ValueRule::Style),
    FieldSpec::new("node_style", "node-style", ValueRule::Style).none(),
    FieldSpec::new("edge_label", "edge-label", ValueRule::ContentOrFunction).none(),
    FieldSpec::new("edge_label_style", "edge-label-style", ValueRule::Style),
];

#[derive(Debug, Default)]
enum SelectorSetting {
    #[default]
    Inherit,
    None,
    Value(Py<PyAny>),
}

impl Clone for SelectorSetting {
    fn clone(&self) -> Self {
        match self {
            Self::Inherit => Self::Inherit,
            Self::None => Self::None,
            Self::Value(value) => Python::attach(|py| Self::Value(value.clone_ref(py))),
        }
    }
}

impl SelectorSetting {
    fn from_python(value: &Bound<'_, PyAny>, name: &str) -> PyResult<Self> {
        if value.extract::<PyRef<'_, PyInherit>>().is_ok() {
            Ok(Self::Inherit)
        } else if value.is_none() {
            Ok(Self::None)
        } else if value.is_callable() {
            Ok(Self::Value(value.clone().unbind()))
        } else {
            Err(PyTypeError::new_err(format!(
                "{name} must be a Python callable, None, or INHERIT"
            )))
        }
    }

    fn overlay(&self, overlay: &Self) -> Self {
        match overlay {
            Self::Inherit => self.clone(),
            _ => overlay.clone(),
        }
    }

    fn resolved(&self, py: Python<'_>) -> Option<Py<PyAny>> {
        match self {
            Self::Value(value) => Some(value.clone_ref(py)),
            Self::Inherit | Self::None => None,
        }
    }
}

#[derive(Clone, Debug, Default)]
struct SelectorSettings {
    node: SelectorSetting,
    edge: SelectorSetting,
    source: SelectorSetting,
    sink: SelectorSetting,
}

impl SelectorSettings {
    fn is_inherit(&self) -> bool {
        [&self.node, &self.edge, &self.source, &self.sink]
            .into_iter()
            .all(|setting| matches!(setting, SelectorSetting::Inherit))
    }

    fn is_none(&self) -> bool {
        [&self.node, &self.edge, &self.source, &self.sink]
            .into_iter()
            .all(|setting| matches!(setting, SelectorSetting::None))
    }

    fn overlay(&self, overlay: &Self) -> Self {
        Self {
            node: self.node.overlay(&overlay.node),
            edge: self.edge.overlay(&overlay.edge),
            source: self.source.overlay(&overlay.source),
            sink: self.sink.overlay(&overlay.sink),
        }
    }

    fn resolved(&self, py: Python<'_>) -> SelectorCallbacks {
        SelectorCallbacks {
            node: self.node.resolved(py),
            edge: self.edge.resolved(py),
            source: self.source.resolved(py),
            sink: self.sink.resolved(py),
        }
    }

    fn traverse(&self, visit: PyVisit<'_>) -> Result<(), PyTraverseError> {
        for selector in [&self.node, &self.edge, &self.source, &self.sink] {
            if let SelectorSetting::Value(selector) = selector {
                visit.call(selector)?;
            }
        }
        Ok(())
    }

    fn clear(&mut self) {
        self.node = SelectorSetting::None;
        self.edge = SelectorSetting::None;
        self.source = SelectorSetting::None;
        self.sink = SelectorSetting::None;
    }
}

/// Python callbacks evaluated against live graph element views before Typst serialization.
#[derive(Debug)]
pub(crate) struct SelectorCallbacks {
    pub(crate) node: Option<Py<PyAny>>,
    pub(crate) edge: Option<Py<PyAny>>,
    pub(crate) source: Option<Py<PyAny>>,
    pub(crate) sink: Option<Py<PyAny>>,
}

/// Per-render Python callbacks returning typed drawing patches.
#[gen_stub_pyclass]
#[pyclass(skip_from_py_object, name = "DrawingSelectors")]
#[derive(Clone, Debug, Default)]
struct PyDrawingSelectors {
    settings: SelectorSettings,
}

#[gen_stub_pymethods]
#[pymethods]
impl PyDrawingSelectors {
    #[new]
    #[pyo3(signature = (**kwargs))]
    #[gen_stub(skip)]
    fn new(kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<Self> {
        let mut settings = SelectorSettings::default();
        if let Some(kwargs) = kwargs {
            for (key, value) in kwargs.iter() {
                let key = key.extract::<String>()?;
                let setting = SelectorSetting::from_python(&value, &key)?;
                match key.as_str() {
                    "node" => settings.node = setting,
                    "edge" => settings.edge = setting,
                    "source" => settings.source = setting,
                    "sink" => settings.sink = setting,
                    _ => {
                        return Err(PyTypeError::new_err(format!(
                            "unknown DrawingSelectors option {key:?}"
                        )))
                    }
                }
            }
        }
        Ok(Self { settings })
    }

    fn __repr__(&self) -> &'static str {
        "DrawingSelectors(...)"
    }

    #[gen_stub(skip)]
    fn __traverse__(&self, visit: PyVisit<'_>) -> Result<(), PyTraverseError> {
        self.settings.traverse(visit)
    }

    #[gen_stub(skip)]
    fn __clear__(&mut self) {
        self.settings.clear();
    }
}

/// Options applied by `linnest.graph.style` before layout measurement.
#[gen_stub_pyclass]
#[pyclass(from_py_object, frozen, name = "GraphStyleOptions")]
#[derive(Clone, Debug, Default, PartialEq)]
struct PyGraphStyleOptions {
    values: BTreeMap<String, NativeValue>,
}

#[gen_stub_pymethods]
#[pymethods]
impl PyGraphStyleOptions {
    #[new]
    #[pyo3(signature = (**kwargs))]
    #[gen_stub(skip)]
    fn new(kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<Self> {
        Ok(Self {
            values: parse_options(kwargs, GRAPH_STYLE_FIELDS)?,
        })
    }

    fn __repr__(&self) -> String {
        options_repr("GraphStyleOptions", &[&self.values])
    }
}

const LAYOUT_FIELDS: &[FieldSpec] = &[
    FieldSpec::new("subgraph", "subgraph", ValueRule::HedgeSelection).none(),
    FieldSpec::new("viewport_width", "viewport-w", ValueRule::NonNegativeNumber),
    FieldSpec::new(
        "viewport_height",
        "viewport-h",
        ValueRule::NonNegativeNumber,
    ),
    FieldSpec::new("tree_dx", "tree-dx", ValueRule::Number),
    FieldSpec::new("tree_dy", "tree-dy", ValueRule::Number),
    FieldSpec::new("steps", "steps", ValueRule::NonNegativeInt),
    FieldSpec::new("seed", "seed", ValueRule::NonNegativeInt),
    FieldSpec::new("step", "step", ValueRule::Number),
    FieldSpec::new("step_shrink", "step-shrink", ValueRule::Number),
    FieldSpec::new("cool", "cool", ValueRule::Number),
    FieldSpec::new("accept_floor", "accept-floor", ValueRule::Number),
    FieldSpec::new("early_tolerance", "early-tol", ValueRule::NonNegativeNumber),
    FieldSpec::new("temperature", "temp", ValueRule::NonNegativeNumber),
    FieldSpec::new("delta", "delta", ValueRule::NonNegativeNumber),
    FieldSpec::new("beta", "beta", ValueRule::Number),
    FieldSpec::new("spring_strength", "k-spring", ValueRule::Number),
    FieldSpec::new("centering_strength", "g-center", ValueRule::Number),
    FieldSpec::new("epochs", "epochs", ValueRule::NonNegativeInt),
    FieldSpec::new("crossing_penalty", "crossing-penalty", ValueRule::Number),
    FieldSpec::new("dangling_repulsion", "gamma-dangling", ValueRule::Number),
    FieldSpec::new("edge_edge_repulsion", "gamma-ee", ValueRule::Number),
    FieldSpec::new("directional_force", "directional-force", ValueRule::Number),
    FieldSpec::new(
        "label_length_scale",
        "label-length-scale",
        ValueRule::Number,
    ),
    FieldSpec::new("label_spring", "label-spring", ValueRule::Number),
    FieldSpec::new("label_charge", "label-charge", ValueRule::Number),
    FieldSpec::new("label_steps", "label-steps", ValueRule::NonNegativeInt),
    FieldSpec::new(
        "label_layout",
        "label-layout",
        ValueRule::Enum(EnumKind::LabelLayout),
    ),
    FieldSpec::new("label_step", "label-step", ValueRule::Number),
    FieldSpec::new(
        "label_early_tolerance",
        "label-early-tol",
        ValueRule::NonNegativeNumber,
    ),
    FieldSpec::new(
        "label_max_delta_scale",
        "label-max-delta-scale",
        ValueRule::NonNegativeNumber,
    ),
    FieldSpec::new("edge_vertex_repulsion", "gamma-ev", ValueRule::Number),
    FieldSpec::new("epsilon", "eps", ValueRule::NonNegativeNumber),
    FieldSpec::new("incremental_energy", "incremental-energy", ValueRule::Bool),
    FieldSpec::new(
        "algorithm",
        "layout-algo",
        ValueRule::Enum(EnumKind::LayoutAlgorithm),
    ),
    FieldSpec::new(
        "nodes",
        "layout-nodes",
        ValueRule::Enum(EnumKind::LayoutNodes),
    ),
    FieldSpec::new(
        "direction",
        "layout-direction",
        ValueRule::Enum(EnumKind::LayoutDirection),
    ),
    FieldSpec::new(
        "rank_align",
        "rank-align",
        ValueRule::Enum(EnumKind::RankAlignment),
    ),
    FieldSpec::new("roots", "layout-roots", ValueRule::NodeIndices),
    FieldSpec::new("rank_same", "rank-same", ValueRule::NodeGroups),
    FieldSpec::new("route_edge_weight", "route-edge-weight", ValueRule::Number),
    FieldSpec::new("route_exit_weight", "route-exit-weight", ValueRule::Number),
    FieldSpec::new(
        "route_label_width_scale",
        "route-label-width-scale",
        ValueRule::Number,
    ),
    FieldSpec::new(
        "route_label_width_cap",
        "route-label-width-cap",
        ValueRule::Number,
    ),
    FieldSpec::new("z_spring", "z-spring", ValueRule::Number),
    FieldSpec::new("z_spring_growth", "z-spring-growth", ValueRule::Number),
    FieldSpec::new("length_scale", "length-scale", ValueRule::Number),
];

/// One or more ordered Linnest layout passes.
#[gen_stub_pyclass]
#[pyclass(from_py_object, frozen, name = "LayoutOptions")]
#[derive(Clone, Debug, PartialEq)]
struct PyLayoutOptions {
    passes: Vec<BTreeMap<String, NativeValue>>,
}

#[gen_stub_pymethods]
#[pymethods]
impl PyLayoutOptions {
    #[new]
    #[pyo3(signature = (**kwargs))]
    #[gen_stub(skip)]
    fn new(kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<Self> {
        Ok(Self {
            passes: vec![parse_options(kwargs, LAYOUT_FIELDS)?],
        })
    }

    #[staticmethod]
    fn sequence(passes: Vec<PyLayoutOptions>) -> PyResult<Self> {
        if passes.is_empty() {
            return Err(PyValueError::new_err(
                "LayoutOptions.sequence needs at least one pass",
            ));
        }
        Ok(Self {
            passes: passes
                .into_iter()
                .flat_map(|options| options.passes)
                .collect(),
        })
    }

    #[pyo3(signature = (**kwargs))]
    #[gen_stub(skip)]
    fn then(&self, kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<Self> {
        let mut passes = self.passes.clone();
        passes.push(parse_options(kwargs, LAYOUT_FIELDS)?);
        Ok(Self { passes })
    }

    #[getter]
    fn pass_count(&self) -> usize {
        self.passes.len()
    }

    fn __repr__(&self) -> String {
        format!("LayoutOptions(passes={})", self.passes.len())
    }
}

impl PyLayoutOptions {
    fn native(&self) -> NativeValue {
        NativeValue::Array(self.passes.iter().cloned().map(NativeValue::Dict).collect())
    }
}

const DRAW_FIELDS: &[FieldSpec] = &[
    FieldSpec::new("scope", "scope", ValueRule::Dictionary),
    FieldSpec::new("unit", "unit", ValueRule::Unit).auto(),
    FieldSpec::new("title", "title", ValueRule::Content)
        .none()
        .auto(),
    FieldSpec::new("subgraph", "subgraph", ValueRule::DrawSubgraphs).none(),
    FieldSpec::new("debug", "debug", ValueRule::Enum(EnumKind::DebugLevel)),
    FieldSpec::new("node_radius", "node-radius", ValueRule::Radius).auto(),
    FieldSpec::new(
        "node_min_radius",
        "node-min-radius",
        ValueRule::NonNegativeNumber,
    ),
    FieldSpec::new(
        "node_label_padding",
        "node-label-padding",
        ValueRule::NonNegativeNumber,
    ),
    FieldSpec::new("node_fill", "node-fill", ValueRule::Paint),
    FieldSpec::new("node_stroke", "node-stroke", ValueRule::Stroke),
    FieldSpec::new("node_outset", "node-outset", ValueRule::Number).auto(),
    FieldSpec::new("node_label_style", "node-label-style", ValueRule::Style),
    FieldSpec::new("node_style", "node-style", ValueRule::Style).none(),
    FieldSpec::new("node_label", "node-label", ValueRule::ContentOrFunction)
        .none()
        .auto(),
    FieldSpec::new("draw_node", "draw-node", ValueRule::Function).auto(),
    FieldSpec::new("edge_stroke", "edge-stroke", ValueRule::Stroke),
    FieldSpec::new("edge_offset", "edge-offset", ValueRule::Number),
    FieldSpec::new("edge_length", "edge-length", ValueRule::Number).none(),
    FieldSpec::new("edge_ratio", "edge-ratio", ValueRule::Number).none(),
    FieldSpec::new(
        "edge_resolve_length",
        "edge-resolve-length",
        ValueRule::EnumOrFunction(EnumKind::EdgeLengthResolution),
    ),
    FieldSpec::new(
        "edge_accuracy",
        "edge-accuracy",
        ValueRule::NonNegativeNumber,
    ),
    FieldSpec::new("edge_optimize", "edge-optimize", ValueRule::Bool),
    FieldSpec::new("source_style", "source-style", ValueRule::StyleLayers).none(),
    FieldSpec::new("sink_style", "sink-style", ValueRule::StyleLayers).none(),
    FieldSpec::new("edge_label", "edge-label", ValueRule::ContentOrFunction).none(),
    FieldSpec::new("edge_label_style", "edge-label-style", ValueRule::Style).none(),
    FieldSpec::new("edge_omega", "edge-omega", ValueRule::Number),
    FieldSpec::new(
        "edge_trim_accuracy",
        "edge-trim-accuracy",
        ValueRule::NonNegativeNumber,
    ),
    FieldSpec::new("padding", "padding", ValueRule::Padding).none(),
    FieldSpec::new(
        "debug_edge_radius",
        "debug-edge-radius",
        ValueRule::NonNegativeNumber,
    ),
    FieldSpec::new("debug_edge_fill", "debug-edge-fill", ValueRule::Paint),
    FieldSpec::new("debug_edge_stroke", "debug-edge-stroke", ValueRule::Stroke),
    FieldSpec::new(
        "debug_edge_label_fill",
        "debug-edge-label-fill",
        ValueRule::Paint,
    ),
    FieldSpec::new(
        "subgraph_edge_style",
        "subgraph-edge-style",
        ValueRule::Style,
    ),
    FieldSpec::new(
        "subgraph_edge_underlay",
        "subgraph-edge-underlay",
        ValueRule::Bool,
    ),
];

/// Full typed option surface for `linnest.draw`.
#[gen_stub_pyclass]
#[pyclass(from_py_object, frozen, name = "DrawOptions")]
#[derive(Clone, Debug, Default, PartialEq)]
struct PyDrawOptions {
    values: BTreeMap<String, NativeValue>,
}

#[gen_stub_pymethods]
#[pymethods]
impl PyDrawOptions {
    #[new]
    #[pyo3(signature = (**kwargs))]
    #[gen_stub(skip)]
    fn new(kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<Self> {
        Ok(Self {
            values: parse_options(kwargs, DRAW_FIELDS)?,
        })
    }

    fn __repr__(&self) -> String {
        options_repr("DrawOptions", &[&self.values])
    }
}

#[derive(Clone, Debug, Default)]
enum PathSetting {
    #[default]
    Inherit,
    None,
    Value(PathBuf),
}

impl PathSetting {
    fn overlay(&self, overlay: &Self) -> Self {
        match overlay {
            Self::Inherit => self.clone(),
            _ => overlay.clone(),
        }
    }

    fn resolved(&self) -> Option<PathBuf> {
        match self {
            Self::Value(path) => Some(path.clone()),
            Self::Inherit | Self::None => None,
        }
    }
}

fn deep_overlay(
    base: &BTreeMap<String, NativeValue>,
    overlay: &BTreeMap<String, NativeValue>,
) -> BTreeMap<String, NativeValue> {
    let mut merged = base.clone();
    for (key, value) in overlay {
        if matches!(value, NativeValue::Inherit) {
            continue;
        }
        match (merged.get(key), value) {
            (Some(NativeValue::Dict(base)), NativeValue::Dict(overlay)) => {
                merged.insert(key.clone(), NativeValue::Dict(deep_overlay(base, overlay)));
            }
            _ => {
                merged.insert(key.clone(), value.clone());
            }
        }
    }
    merged
}

/// Complete typed rendering configuration.
#[gen_stub_pyclass]
#[pyclass(skip_from_py_object, name = "RenderConfig")]
#[derive(Clone, Debug, Default)]
pub(crate) struct PyRenderConfig {
    template: PathSetting,
    values: BTreeMap<String, NativeValue>,
    selectors: SelectorSettings,
}

#[gen_stub_pymethods]
#[pymethods]
impl PyRenderConfig {
    #[new]
    #[pyo3(signature = (**kwargs))]
    #[gen_stub(skip)]
    fn new(kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<Self> {
        let mut config = Self::default();
        let Some(kwargs) = kwargs else {
            return Ok(config);
        };
        for (key, value) in kwargs.iter() {
            let key = key.extract::<String>()?;
            config.assign(&key, &value)?;
        }
        Ok(config)
    }

    fn overlay(&self, overlay: &PyRenderConfig) -> Self {
        self.merged(overlay)
    }

    #[getter]
    #[gen_stub(override_return_type(type_repr = "_TemplatePath"))]
    fn template(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        self.path_setting(py, &self.template)
    }

    #[setter]
    fn set_template(
        &mut self,
        #[gen_stub(override_type(type_repr = "_TemplatePath"))] value: &Bound<'_, PyAny>,
    ) -> PyResult<()> {
        self.assign("template", value)
    }

    #[getter]
    #[gen_stub(override_return_type(type_repr = "_AutoOptionalStaticContent"))]
    fn title(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        self.native_setting(py, "title")
    }

    #[setter]
    fn set_title(
        &mut self,
        #[gen_stub(override_type(type_repr = "_AutoOptionalStaticContent"))] value: &Bound<
            '_,
            PyAny,
        >,
    ) -> PyResult<()> {
        self.assign("title", value)
    }

    #[getter]
    #[gen_stub(override_return_type(type_repr = "_RenderStyle"))]
    fn style(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        match self.values.get("style") {
            None => Ok(inherit(py)),
            Some(NativeValue::None) => Ok(py.None()),
            Some(NativeValue::Dict(values)) => Ok(Py::new(
                py,
                PyGraphStyleOptions {
                    values: values.clone(),
                },
            )?
            .into_any()),
            Some(_) => Err(PyValueError::new_err(
                "internal error: RenderConfig style has an invalid value",
            )),
        }
    }

    #[setter]
    fn set_style(
        &mut self,
        #[gen_stub(override_type(type_repr = "_RenderStyle"))] value: &Bound<'_, PyAny>,
    ) -> PyResult<()> {
        self.assign("style", value)
    }

    #[getter]
    #[gen_stub(override_return_type(type_repr = "_RenderLayouts"))]
    fn layouts(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        match self.values.get("layouts") {
            None => Ok(inherit(py)),
            Some(NativeValue::None) => Ok(py.None()),
            Some(NativeValue::Array(passes)) => {
                let passes = passes
                    .iter()
                    .map(|pass| match pass {
                        NativeValue::Dict(values) => Ok(values.clone()),
                        _ => Err(PyValueError::new_err(
                            "internal error: RenderConfig layouts contains an invalid pass",
                        )),
                    })
                    .collect::<PyResult<Vec<_>>>()?;
                Ok(Py::new(py, PyLayoutOptions { passes })?.into_any())
            }
            Some(_) => Err(PyValueError::new_err(
                "internal error: RenderConfig layouts has an invalid value",
            )),
        }
    }

    #[setter]
    fn set_layouts(
        &mut self,
        #[gen_stub(override_type(type_repr = "_RenderLayouts"))] value: &Bound<'_, PyAny>,
    ) -> PyResult<()> {
        self.assign("layouts", value)
    }

    #[getter]
    #[gen_stub(override_return_type(type_repr = "_RenderDrawing"))]
    fn drawing(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        match self.values.get("draw") {
            None => Ok(inherit(py)),
            Some(NativeValue::None) => Ok(py.None()),
            Some(NativeValue::Dict(values)) => Ok(Py::new(
                py,
                PyDrawOptions {
                    values: values.clone(),
                },
            )?
            .into_any()),
            Some(_) => Err(PyValueError::new_err(
                "internal error: RenderConfig drawing has an invalid value",
            )),
        }
    }

    #[setter]
    fn set_drawing(
        &mut self,
        #[gen_stub(override_type(type_repr = "_RenderDrawing"))] value: &Bound<'_, PyAny>,
    ) -> PyResult<()> {
        self.assign("drawing", value)
    }

    #[getter]
    #[gen_stub(override_return_type(type_repr = "_RenderSelectors"))]
    fn selectors(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        if self.selectors.is_inherit() {
            Ok(inherit(py))
        } else if self.selectors.is_none() {
            Ok(py.None())
        } else {
            Ok(Py::new(
                py,
                PyDrawingSelectors {
                    settings: self.selectors.clone(),
                },
            )?
            .into_any())
        }
    }

    #[setter]
    fn set_selectors(
        &mut self,
        #[gen_stub(override_type(type_repr = "_RenderSelectors"))] value: &Bound<'_, PyAny>,
    ) -> PyResult<()> {
        self.assign("selectors", value)
    }

    #[getter]
    #[gen_stub(override_return_type(type_repr = "_TemplateOptions"))]
    fn template_options(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        self.native_setting(py, "options")
    }

    #[setter]
    fn set_template_options(
        &mut self,
        #[gen_stub(override_type(type_repr = "_TemplateOptions"))] value: &Bound<'_, PyAny>,
    ) -> PyResult<()> {
        self.assign("template_options", value)
    }

    fn __repr__(&self) -> String {
        format!(
            "RenderConfig(template={:?}, fields={})",
            self.template.resolved(),
            self.values.len()
        )
    }

    #[gen_stub(skip)]
    fn __traverse__(&self, visit: PyVisit<'_>) -> Result<(), PyTraverseError> {
        self.selectors.traverse(visit)
    }

    #[gen_stub(skip)]
    fn __clear__(&mut self) {
        self.selectors.clear();
    }
}

pyo3_stub_gen::inventory::submit! {
    pyo3_stub_gen::derive::gen_methods_from_python! { r#"
        class PyStroke:
            def __new__(cls, *, paint: _OptionalPaint = ..., thickness: _OptionalLengthValue = ..., cap: _OptionalStrokeCap = ..., join: _OptionalStrokeJoin = ..., dash: _OptionalDashValue = ..., miter_limit: _Number = ...) -> Stroke: ...
    "# }
}

pyo3_stub_gen::inventory::submit! {
    pyo3_stub_gen::derive::gen_methods_from_python! { r#"
        class PyInsets:
            def __new__(cls, *, all: _LengthValue = ..., x: _LengthValue = ..., y: _LengthValue = ..., left: _LengthValue = ..., right: _LengthValue = ..., top: _LengthValue = ..., bottom: _LengthValue = ...) -> Insets: ...
    "# }
}

pyo3_stub_gen::inventory::submit! {
    pyo3_stub_gen::derive::gen_methods_from_python! { r#"
        class PyMark:
            def __new__(cls, *, start: _NativeValue = ..., end: _NativeValue = ..., fill: _OptionalPaint = ..., stroke: _OptionalStrokeValue = ..., scale: _Number = ..., anchor: _MarkAnchor = ..., shorten_to: _MarkShorten = ...) -> Mark: ...
    "# }
}

pyo3_stub_gen::inventory::submit! {
    pyo3_stub_gen::derive::gen_methods_from_python! { r#"
        class PyTextLabel:
            def __new__(cls, text: str, *, fill: _Paint = ..., size: _LengthValue = ..., weight: _NativeValue = ..., style: _TextStyleValue = ...) -> TextLabel: ...
    "# }
}

pyo3_stub_gen::inventory::submit! {
    pyo3_stub_gen::derive::gen_methods_from_python! { r#"
        import typing

        class PyGraphStyleOptions:
            def __new__(cls, *, scope: _Dictionary = ..., unit: _LengthValue = ..., node_label: _AutoOptionalContent = ..., node_label_style: _Style = ..., node_style: _OptionalStyle = ..., edge_label: _OptionalContent = ..., edge_label_style: _Style = ...) -> GraphStyleOptions: ...
    "# }
}

pyo3_stub_gen::inventory::submit! {
    pyo3_stub_gen::derive::gen_methods_from_python! { r#"
        class PyDrawingSelectors:
            def __new__(cls, *, node: _NodeDrawingSelector = ..., edge: _EdgeDrawingSelector = ..., source: _HalfEdgeDrawingSelector = ..., sink: _HalfEdgeDrawingSelector = ...) -> DrawingSelectors: ...
    "# }
}

pyo3_stub_gen::inventory::submit! {
    pyo3_stub_gen::derive::gen_methods_from_python! { r#"
        import typing

        class PyLayoutOptions:
            def __new__(cls, *, subgraph: _OptionalHedgeSelection = ..., viewport_width: _Number = ..., viewport_height: _Number = ..., tree_dx: _Number = ..., tree_dy: _Number = ..., steps: _Integer = ..., seed: _Integer = ..., step: _Number = ..., step_shrink: _Number = ..., cool: _Number = ..., accept_floor: _Number = ..., early_tolerance: _Number = ..., temperature: _Number = ..., delta: _Number = ..., beta: _Number = ..., spring_strength: _Number = ..., centering_strength: _Number = ..., epochs: _Integer = ..., crossing_penalty: _Number = ..., dangling_repulsion: _Number = ..., edge_edge_repulsion: _Number = ..., directional_force: _Number = ..., label_length_scale: _Number = ..., label_spring: _Number = ..., label_charge: _Number = ..., label_steps: _Integer = ..., label_layout: _LabelLayoutValue = ..., label_step: _Number = ..., label_early_tolerance: _Number = ..., label_max_delta_scale: _Number = ..., edge_vertex_repulsion: _Number = ..., epsilon: _Number = ..., incremental_energy: _Boolean = ..., algorithm: _LayoutAlgorithmValue = ..., nodes: _LayoutNodesValue = ..., direction: _LayoutDirectionValue = ..., rank_align: _RankAlignmentValue = ..., roots: _NodeIndices = ..., rank_same: _NodeGroups = ..., route_edge_weight: _Number = ..., route_exit_weight: _Number = ..., route_label_width_scale: _Number = ..., route_label_width_cap: _Number = ..., z_spring: _Number = ..., z_spring_growth: _Number = ..., length_scale: _Number = ...) -> LayoutOptions: ...
            def then(self, *, subgraph: _OptionalHedgeSelection = ..., viewport_width: _Number = ..., viewport_height: _Number = ..., tree_dx: _Number = ..., tree_dy: _Number = ..., steps: _Integer = ..., seed: _Integer = ..., step: _Number = ..., step_shrink: _Number = ..., cool: _Number = ..., accept_floor: _Number = ..., early_tolerance: _Number = ..., temperature: _Number = ..., delta: _Number = ..., beta: _Number = ..., spring_strength: _Number = ..., centering_strength: _Number = ..., epochs: _Integer = ..., crossing_penalty: _Number = ..., dangling_repulsion: _Number = ..., edge_edge_repulsion: _Number = ..., directional_force: _Number = ..., label_length_scale: _Number = ..., label_spring: _Number = ..., label_charge: _Number = ..., label_steps: _Integer = ..., label_layout: _LabelLayoutValue = ..., label_step: _Number = ..., label_early_tolerance: _Number = ..., label_max_delta_scale: _Number = ..., edge_vertex_repulsion: _Number = ..., epsilon: _Number = ..., incremental_energy: _Boolean = ..., algorithm: _LayoutAlgorithmValue = ..., nodes: _LayoutNodesValue = ..., direction: _LayoutDirectionValue = ..., rank_align: _RankAlignmentValue = ..., roots: _NodeIndices = ..., rank_same: _NodeGroups = ..., route_edge_weight: _Number = ..., route_exit_weight: _Number = ..., route_label_width_scale: _Number = ..., route_label_width_cap: _Number = ..., z_spring: _Number = ..., z_spring_growth: _Number = ..., length_scale: _Number = ...) -> LayoutOptions: ...
    "# }
}

pyo3_stub_gen::inventory::submit! {
    pyo3_stub_gen::derive::gen_methods_from_python! { r#"
        import typing

        class PyDrawOptions:
            def __new__(cls, *, scope: _Dictionary = ..., unit: _AutoLengthValue = ..., title: _AutoOptionalStaticContent = ..., subgraph: _DrawSubgraphs = ..., debug: _DebugValue = ..., node_radius: _AutoRadius = ..., node_min_radius: _Number = ..., node_label_padding: _Number = ..., node_fill: _Paint = ..., node_stroke: _StrokeValue = ..., node_outset: _AutoNumber = ..., node_label_style: _Style = ..., node_style: _OptionalStyle = ..., node_label: _AutoOptionalContent = ..., draw_node: _AutoFunction = ..., edge_stroke: _StrokeValue = ..., edge_offset: _Number = ..., edge_length: _OptionalNumber = ..., edge_ratio: _OptionalNumber = ..., edge_resolve_length: _EdgeLengthResolver = ..., edge_accuracy: _Number = ..., edge_optimize: _Boolean = ..., source_style: _OptionalStyleLayers = ..., sink_style: _OptionalStyleLayers = ..., edge_label: _OptionalContent = ..., edge_label_style: _OptionalStyle = ..., edge_omega: _Number = ..., edge_trim_accuracy: _Number = ..., padding: _OptionalPadding = ..., debug_edge_radius: _Number = ..., debug_edge_fill: _Paint = ..., debug_edge_stroke: _StrokeValue = ..., debug_edge_label_fill: _Paint = ..., subgraph_edge_style: _Style = ..., subgraph_edge_underlay: _Boolean = ...) -> DrawOptions: ...
    "# }
}

pyo3_stub_gen::inventory::submit! {
    pyo3_stub_gen::derive::gen_methods_from_python! { r#"
        class PyRenderConfig:
            def __new__(cls, *, template: _TemplatePath = ..., title: _AutoOptionalStaticContent = ..., style: _RenderStyle = ..., layouts: _RenderLayouts = ..., drawing: _RenderDrawing = ..., selectors: _RenderSelectors = ..., template_options: _TemplateOptions = ...) -> RenderConfig: ...
    "# }
}

impl PyRenderConfig {
    fn assign(&mut self, key: &str, value: &Bound<'_, PyAny>) -> PyResult<()> {
        match key {
            "template" => {
                self.template = if value.extract::<PyRef<'_, PyInherit>>().is_ok() {
                    PathSetting::Inherit
                } else if value.is_none() {
                    PathSetting::None
                } else {
                    PathSetting::Value(value.extract::<PathBuf>().map_err(|_| {
                        PyTypeError::new_err("template must be a path, None, or INHERIT")
                    })?)
                };
            }
            "title" => {
                let title = native_from_py(value, 0)?;
                title.validate(0)?;
                if matches!(title, NativeValue::Inherit) {
                    self.values.remove("title");
                } else if matches!(title, NativeValue::None | NativeValue::Auto)
                    || ValueRule::Content.accepts(&title)
                {
                    self.values.insert("title".to_owned(), title);
                } else {
                    return Err(PyTypeError::new_err(
                        "title must be static text/content, None, AUTO, or INHERIT",
                    ));
                }
            }
            "style" => {
                if value.extract::<PyRef<'_, PyInherit>>().is_ok() {
                    self.values.remove("style");
                } else if value.is_none() {
                    self.values.insert("style".to_owned(), NativeValue::None);
                } else {
                    let options =
                        value
                            .extract::<PyRef<'_, PyGraphStyleOptions>>()
                            .map_err(|_| {
                                PyTypeError::new_err(
                                    "style must be GraphStyleOptions, None, or INHERIT",
                                )
                            })?;
                    self.values.insert(
                        "style".to_owned(),
                        NativeValue::Dict(options.values.clone()),
                    );
                }
            }
            "layouts" => {
                if value.extract::<PyRef<'_, PyInherit>>().is_ok() {
                    self.values.remove("layouts");
                } else if value.is_none() {
                    self.values.insert("layouts".to_owned(), NativeValue::None);
                } else {
                    let options = value.extract::<PyRef<'_, PyLayoutOptions>>().map_err(|_| {
                        PyTypeError::new_err("layouts must be LayoutOptions, None, or INHERIT")
                    })?;
                    self.values.insert("layouts".to_owned(), options.native());
                }
            }
            "drawing" => {
                if value.extract::<PyRef<'_, PyInherit>>().is_ok() {
                    self.values.remove("draw");
                } else if value.is_none() {
                    self.values.insert("draw".to_owned(), NativeValue::None);
                } else {
                    let options = value.extract::<PyRef<'_, PyDrawOptions>>().map_err(|_| {
                        PyTypeError::new_err("drawing must be DrawOptions, None, or INHERIT")
                    })?;
                    self.values
                        .insert("draw".to_owned(), NativeValue::Dict(options.values.clone()));
                }
            }
            "selectors" => {
                if value.extract::<PyRef<'_, PyInherit>>().is_ok() {
                    self.selectors = SelectorSettings::default();
                } else if value.is_none() {
                    self.selectors = SelectorSettings {
                        node: SelectorSetting::None,
                        edge: SelectorSetting::None,
                        source: SelectorSetting::None,
                        sink: SelectorSetting::None,
                    };
                } else {
                    let options =
                        value
                            .extract::<PyRef<'_, PyDrawingSelectors>>()
                            .map_err(|_| {
                                PyTypeError::new_err(
                                    "selectors must be DrawingSelectors, None, or INHERIT",
                                )
                            })?;
                    self.selectors = options.settings.clone();
                }
            }
            "template_options" => {
                if value.extract::<PyRef<'_, PyInherit>>().is_ok() {
                    self.values.remove("options");
                } else if value.is_none() {
                    self.values.insert("options".to_owned(), NativeValue::None);
                } else {
                    let options = native_from_py(value, 0)?;
                    if !matches!(options, NativeValue::Dict(_)) {
                        return Err(PyTypeError::new_err(
                            "template_options must be a dictionary, None, or INHERIT",
                        ));
                    }
                    options.validate(0)?;
                    self.values.insert("options".to_owned(), options);
                }
            }
            _ => {
                return Err(PyTypeError::new_err(format!(
                    "unknown RenderConfig option {key:?}"
                )));
            }
        }
        Ok(())
    }

    fn path_value(&self, py: Python<'_>, path: &Path) -> PyResult<Py<PyAny>> {
        Ok(path.to_path_buf().into_pyobject(py)?.into_any().unbind())
    }

    fn path_setting(&self, py: Python<'_>, setting: &PathSetting) -> PyResult<Py<PyAny>> {
        match setting {
            PathSetting::Inherit => Ok(inherit(py)),
            PathSetting::None => Ok(py.None()),
            PathSetting::Value(path) => self.path_value(py, path),
        }
    }

    fn native_setting(&self, py: Python<'_>, key: &str) -> PyResult<Py<PyAny>> {
        self.values
            .get(key)
            .map_or_else(|| Ok(inherit(py)), |value| native_to_py(py, value))
    }

    fn merged(&self, overlay: &Self) -> Self {
        Self {
            template: self.template.overlay(&overlay.template),
            values: deep_overlay(&self.values, &overlay.values),
            selectors: self.selectors.overlay(&overlay.selectors),
        }
    }

    fn native(&self, elements: NativeValue) -> NativeValue {
        let mut values = self.values.clone();
        values.insert(
            "schema".to_owned(),
            NativeValue::String("linnet-render-config".to_owned()),
        );
        values.insert("version".to_owned(), NativeValue::Int(1));
        values.insert("elements".to_owned(), elements);
        NativeValue::Dict(values)
    }
}

/// Resolved renderer transport plus its closed native Typst configuration.
#[derive(Debug)]
pub(crate) struct RenderConfigTransport {
    pub(crate) config: clinnet::TypstConfig,
    pub(crate) imports: Vec<TypstImport>,
    pub(crate) template: Option<PathBuf>,
    pub(crate) selectors: SelectorCallbacks,
}

impl NativeValue {
    fn validate(&self, depth: usize) -> PyResult<()> {
        if depth > MAX_NATIVE_DEPTH {
            return Err(PyValueError::new_err(format!(
                "native Typst value exceeds maximum nesting depth {MAX_NATIVE_DEPTH}"
            )));
        }
        match self {
            Self::Float(value) | Self::Ratio(value) | Self::Length(value, _) => {
                finite(*value, self.kind())?;
            }
            Self::Angle(value, unit) => {
                finite(*value, self.kind())?;
                if !matches!(unit.as_str(), "deg" | "rad") {
                    return Err(PyValueError::new_err(format!(
                        "angle unit must be 'deg' or 'rad', not {unit:?}"
                    )));
                }
            }
            Self::Fraction(value) => {
                non_negative(*value, "fraction")?;
            }
            Self::RelativeLength { ratio, length } => {
                if ratio.is_none() && length.is_none() {
                    return Err(PyValueError::new_err(
                        "decoded RelativeLength has no ratio or length",
                    ));
                }
                if let Some(ratio) = ratio {
                    finite(*ratio, "relative ratio")?;
                }
                if let Some((length, _)) = length {
                    finite(*length, "relative length")?;
                }
            }
            Self::Array(values) => {
                for value in values {
                    value.validate(depth + 1)?;
                    if matches!(value, Self::Inherit) {
                        return Err(PyValueError::new_err(
                            "INHERIT is not valid inside a Typst array",
                        ));
                    }
                }
            }
            Self::Dict(values) => {
                for value in values.values() {
                    value.validate(depth + 1)?;
                }
            }
            Self::Stroke(values) => {
                validate_option_values(values, STROKE_FIELDS, "Stroke", depth)?;
            }
            Self::Insets(values) => {
                validate_option_values(values, INSET_FIELDS, "Insets", depth)?;
                let explicit = values
                    .iter()
                    .filter(|(_, value)| !matches!(value, Self::Inherit))
                    .collect::<Vec<_>>();
                if explicit.is_empty() {
                    return Err(PyValueError::new_err(
                        "Insets needs at least one inset value",
                    ));
                }
                if explicit.iter().any(|(name, _)| name.as_str() == "all") && explicit.len() > 1 {
                    return Err(PyValueError::new_err(
                        "Insets(all=...) cannot be combined with directional values",
                    ));
                }
            }
            Self::Mark(values) => {
                validate_option_values(values, MARK_FIELDS, "Mark", depth)?;
                if values.values().all(|value| matches!(value, Self::Inherit)) {
                    return Err(PyValueError::new_err("Mark needs at least one option"));
                }
            }
            Self::Color(ColorValue::Named(value) | ColorValue::Hex(value)) => {
                PyColor::new(value)?;
            }
            Self::Dash(DashValue::Pattern { array, phase }) => {
                if array.is_empty() {
                    return Err(PyValueError::new_err(
                        "decoded Dash pattern cannot be empty",
                    ));
                }
                for value in array {
                    value.validate(depth + 1)?;
                    if !ValueRule::Length.accepts(value) {
                        return Err(PyValueError::new_err(
                            "decoded Dash contains a non-length entry",
                        ));
                    }
                }
                if let Some(phase) = phase {
                    phase.validate(depth + 1)?;
                    if !ValueRule::Length.accepts(phase) {
                        return Err(PyValueError::new_err("decoded Dash phase is not a length"));
                    }
                }
            }
            Self::Dash(DashValue::Named(value)) => {
                if !matches!(
                    value.as_str(),
                    "solid"
                        | "dotted"
                        | "densely-dotted"
                        | "loosely-dotted"
                        | "dashed"
                        | "densely-dashed"
                        | "loosely-dashed"
                        | "dash-dotted"
                        | "densely-dash-dotted"
                        | "loosely-dash-dotted"
                ) {
                    return Err(PyValueError::new_err(format!(
                        "unknown Typst dash pattern {value:?}"
                    )));
                }
            }
            Self::Text(value) => {
                validate_option_values(&value.options, TEXT_FIELDS, "TextLabel", depth)?;
            }
            Self::Math(value) => {
                math_token(&value.symbol, "math symbol")?;
                if let Some(MathScript::Symbol(subscript)) = &value.subscript {
                    math_token(subscript, "math subscript")?;
                }
                if let Some(MathScript::Symbol(superscript)) = &value.superscript {
                    math_token(superscript, "math superscript")?;
                }
            }
            Self::Expression(expression) => expression.validate(depth + 1)?,
            Self::Inherit
            | Self::None
            | Self::Auto
            | Self::Bool(_)
            | Self::Int(_)
            | Self::String(_)
            | Self::Enum(_)
            | Self::Color(ColorValue::Rgba(..)) => {}
        }
        Ok(())
    }

    fn reject_dot_module_references(&self) -> PyResult<()> {
        let mut modules = BTreeSet::new();
        self.collect_modules(&mut modules);
        if modules.is_empty() {
            Ok(())
        } else {
            Err(PyValueError::new_err(
                "DotCodec.linnest() does not serialize Typst module references; use a custom DotCodec to reconstruct trusted executable values explicitly",
            ))
        }
    }
}

fn validate_option_values(
    values: &BTreeMap<String, NativeValue>,
    fields: &[FieldSpec],
    kind: &str,
    depth: usize,
) -> PyResult<()> {
    for (name, value) in values {
        let field = fields
            .iter()
            .find(|field| field.typst == name)
            .ok_or_else(|| PyValueError::new_err(format!("unknown {kind} option {name:?}")))?;
        value.validate(depth + 1)?;
        field.validate(value)?;
    }
    Ok(())
}

impl TypstExpression {
    fn validate(&self, depth: usize) -> PyResult<()> {
        if depth > MAX_NATIVE_DEPTH {
            return Err(PyValueError::new_err(format!(
                "native Typst value exceeds maximum nesting depth {MAX_NATIVE_DEPTH}"
            )));
        }
        match self {
            Self::Symbol(symbol) => {
                if symbol.path.is_empty() {
                    return Err(PyValueError::new_err(
                        "decoded Typst module reference has no export",
                    ));
                }
                for part in &symbol.path {
                    identifier(part, "module export")?;
                }
                if let TypstModuleSource::Package(package) = &symbol.module {
                    clinnet::TypstModule::package("_linnet_validate", package.clone()).map_err(
                        |_| PyValueError::new_err("decoded Typst package reference is invalid"),
                    )?;
                }
            }
            Self::Call { callee, arguments } => {
                if !matches!(
                    callee.as_ref(),
                    Self::Symbol(ModuleSymbol {
                        kind: SymbolKind::Function,
                        ..
                    }) | Self::Bind { .. }
                ) {
                    return Err(PyValueError::new_err(
                        "a Typst call must target a module function or binding",
                    ));
                }
                callee.validate(depth + 1)?;
                arguments.validate(depth + 1)?;
            }
            Self::Bind {
                function,
                arguments,
            } => {
                if !matches!(
                    function.as_ref(),
                    Self::Symbol(ModuleSymbol {
                        kind: SymbolKind::Function,
                        ..
                    })
                ) {
                    return Err(PyValueError::new_err(
                        "a Typst binding must target an imported module function",
                    ));
                }
                function.validate(depth + 1)?;
                arguments.validate(depth + 1)?;
            }
        }
        Ok(())
    }
}

fn native_to_py(py: Python<'_>, value: &NativeValue) -> PyResult<Py<PyAny>> {
    let value = match value {
        NativeValue::Inherit => inherit_singleton(py).into_any(),
        NativeValue::None => py.None(),
        NativeValue::Auto => auto_singleton(py).into_any(),
        NativeValue::Bool(value) => PyBool::new(py, *value).to_owned().into_any().unbind(),
        NativeValue::Int(value) => PyInt::new(py, *value).into_any().unbind(),
        NativeValue::Float(value) => PyFloat::new(py, *value).into_any().unbind(),
        NativeValue::String(value) => PyString::new(py, value).into_any().unbind(),
        NativeValue::Enum(value) => value.to_py(py)?,
        NativeValue::Array(values) => {
            let values = values
                .iter()
                .map(|value| native_to_py(py, value))
                .collect::<PyResult<Vec<_>>>()?;
            PyTuple::new(py, values)?.into_any().unbind()
        }
        NativeValue::Dict(values) => {
            let dictionary = PyDict::new(py);
            for (key, value) in values {
                dictionary.set_item(key, native_to_py(py, value)?)?;
            }
            dictionary.into_any().unbind()
        }
        NativeValue::Length(value, unit) => Py::new(
            py,
            PyLength {
                value: *value,
                unit: unit.clone(),
            },
        )?
        .into_any(),
        NativeValue::Ratio(percent) => Py::new(py, PyRatio { percent: *percent })?.into_any(),
        NativeValue::RelativeLength { ratio, length } => Py::new(
            py,
            PyRelativeLength {
                ratio: *ratio,
                length: length.clone(),
            },
        )?
        .into_any(),
        NativeValue::Angle(value, unit) => Py::new(
            py,
            PyAngle {
                value: *value,
                unit: unit.clone(),
            },
        )?
        .into_any(),
        NativeValue::Fraction(value) => Py::new(py, PyFraction { value: *value })?.into_any(),
        NativeValue::Color(value) => Py::new(
            py,
            PyColor {
                value: value.clone(),
            },
        )?
        .into_any(),
        NativeValue::Stroke(values) => Py::new(
            py,
            PyStroke {
                values: values.clone(),
            },
        )?
        .into_any(),
        NativeValue::Dash(value) => Py::new(
            py,
            PyDash {
                value: value.clone(),
            },
        )?
        .into_any(),
        NativeValue::Insets(values) => Py::new(
            py,
            PyInsets {
                values: values.clone(),
            },
        )?
        .into_any(),
        NativeValue::Mark(values) => Py::new(
            py,
            PyMark {
                values: values.clone(),
            },
        )?
        .into_any(),
        NativeValue::Text(value) => Py::new(
            py,
            PyTextLabel {
                value: value.clone(),
            },
        )?
        .into_any(),
        NativeValue::Math(value) => Py::new(
            py,
            PyMathSymbol {
                value: value.clone(),
            },
        )?
        .into_any(),
        NativeValue::Expression(expression @ TypstExpression::Symbol(_)) => Py::new(
            py,
            PyTypstRef {
                expression: expression.clone(),
            },
        )?
        .into_any(),
        NativeValue::Expression(expression @ TypstExpression::Call { .. }) => Py::new(
            py,
            PyTypstCall {
                expression: expression.clone(),
            },
        )?
        .into_any(),
        NativeValue::Expression(expression @ TypstExpression::Bind { .. }) => Py::new(
            py,
            PyTypstBind {
                expression: expression.clone(),
            },
        )?
        .into_any(),
    };
    Ok(value)
}

const DOT_NATIVE_SCHEMA: &str = "linnet-typst-native";
const DOT_NATIVE_VERSION: u8 = 1;

#[derive(Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct DotNativeEnvelope {
    schema: String,
    version: u8,
    value: NativeValue,
}

fn auto_singleton(py: Python<'_>) -> Py<PyAuto> {
    AUTO_SINGLETON
        .get_or_init(py, || {
            Py::new(py, PyAuto).expect("allocating the zero-sized AUTO sentinel cannot fail")
        })
        .clone_ref(py)
}

fn inherit_singleton(py: Python<'_>) -> Py<PyInherit> {
    INHERIT_SINGLETON
        .get_or_init(py, || {
            Py::new(py, PyInherit).expect("allocating the zero-sized INHERIT sentinel cannot fail")
        })
        .clone_ref(py)
}

/// Return the shared INHERIT sentinel for graph-owned option dictionaries.
pub(crate) fn inherit(py: Python<'_>) -> Py<PyAny> {
    inherit_singleton(py).into_any()
}

const NODE_DRAWING_FIELDS: &[FieldSpec] = &[
    FieldSpec::new("label", "label", ValueRule::Content).none(),
    FieldSpec::new("shift", "shift", ValueRule::Point).none(),
    FieldSpec::new("rank", "rank", ValueRule::NonNegativeInt).none(),
    FieldSpec::new("minimum_size", "minimum_size", ValueRule::NonNegativeNumber).none(),
    FieldSpec::new("maximum_size", "maximum_size", ValueRule::NonNegativeNumber).none(),
    FieldSpec::new("style", "style", ValueRule::Style).none(),
    FieldSpec::new("label_style", "label_style", ValueRule::Style).none(),
];

const EDGE_DRAWING_FIELDS: &[FieldSpec] = &[
    FieldSpec::new("label", "label", ValueRule::Content).none(),
    FieldSpec::new("label_position", "label_position", ValueRule::Point).none(),
    FieldSpec::new("label_offset", "label_offset", ValueRule::Number).none(),
    FieldSpec::new("label_angle", "label_angle", ValueRule::Angle).none(),
    FieldSpec::new("bend", "bend", ValueRule::Angle).none(),
    FieldSpec::new("routing", "routing", ValueRule::Enum(EnumKind::Routing)).none(),
    FieldSpec::new(
        "minimum_length",
        "minimum_length",
        ValueRule::NonNegativeInt,
    )
    .none(),
    FieldSpec::new("same_rank", "same_rank", ValueRule::Bool).none(),
    FieldSpec::new("style", "style", ValueRule::StyleLayers).none(),
    FieldSpec::new("label_style", "label_style", ValueRule::Style).none(),
    FieldSpec::new("decoration", "decoration", ValueRule::PatternOrStyleLayers).none(),
];

const HALF_EDGE_DRAWING_FIELDS: &[FieldSpec] = &[
    FieldSpec::new("label", "label", ValueRule::Content).none(),
    FieldSpec::new("statement", "statement", ValueRule::String).none(),
    FieldSpec::new("port_label", "port_label", ValueRule::String).none(),
    FieldSpec::new("compass", "compass", ValueRule::Enum(EnumKind::Compass)).none(),
    FieldSpec::new("anchor", "anchor", ValueRule::Enum(EnumKind::Anchor))
        .none()
        .auto(),
    FieldSpec::new("routing", "routing", ValueRule::Enum(EnumKind::Routing)).none(),
    FieldSpec::new("style", "style", ValueRule::StyleLayers).none(),
];

fn drawing_field_spec(kind: &str, field: &str) -> Option<FieldSpec> {
    let fields = match kind {
        "node" | "NodeDrawing" => NODE_DRAWING_FIELDS,
        "edge" | "EdgeDrawing" => EDGE_DRAWING_FIELDS,
        "half-edge" | "HalfEdgeDrawing" => HALF_EDGE_DRAWING_FIELDS,
        _ => return None,
    };
    fields
        .iter()
        .copied()
        .find(|candidate| candidate.python == field)
}

/// Validate that a Python value belongs to the closed native Typst model.
pub(crate) fn validate_native(value: &Bound<'_, PyAny>) -> PyResult<()> {
    native_from_py(value, 0)?.validate(0)
}

fn drawing_enum_rule(field: &str) -> Option<(EnumKind, bool, bool)> {
    match field {
        "compass" => Some((EnumKind::Compass, true, false)),
        "routing" | "route" => Some((EnumKind::Routing, true, false)),
        "anchor" | "source-anchor" | "sink-anchor" => Some((EnumKind::Anchor, true, true)),
        "pattern" => Some((EnumKind::Pattern, true, false)),
        "debug" => Some((EnumKind::DebugLevel, false, false)),
        "label_layout" | "label-layout" => Some((EnumKind::LabelLayout, false, false)),
        "route_points" | "route-points" => Some((EnumKind::RoutePoints, true, false)),
        "edge_resolve_length" | "edge-resolve-length" | "resolve_length" | "resolve-length" => {
            Some((EnumKind::EdgeLengthResolution, true, false))
        }
        "mark_position" | "mark-position" => Some((EnumKind::MarkPosition, false, false)),
        "mark_orientation" | "mark-orientation" => Some((EnumKind::MarkOrientation, false, false)),
        "mark_direction" | "mark-direction" => Some((EnumKind::MarkDirection, false, false)),
        "cap" => Some((EnumKind::StrokeCap, true, false)),
        "join" => Some((EnumKind::StrokeJoin, true, false)),
        _ => None,
    }
}

fn validate_placement(value: &NativeValue) -> PyResult<()> {
    match value {
        NativeValue::Inherit | NativeValue::None | NativeValue::Expression(_) => Ok(()),
        NativeValue::Enum(NativeEnum::Placement(PyPlacement::Start)) => Ok(()),
        NativeValue::Enum(NativeEnum::Placement(PyPlacement::Pin)) => Err(PyValueError::new_err(
            "pin placement must constrain x, y, ref, dx, or dy",
        )),
        NativeValue::Dict(values) => {
            const FIELDS: &[&str] = &["x", "y", "ref", "dx", "dy", "mode", "x-mode", "y-mode"];
            if let Some(field) = values.keys().find(|field| !FIELDS.contains(&field.as_str())) {
                return Err(PyTypeError::new_err(format!(
                    "unknown placement field {field:?}"
                )));
            }
            if let Some(mode) = values.get("mode") {
                validate_drawing_enum("mode", mode, (EnumKind::Placement, false, false))?;
            }
            for field in ["x-mode", "y-mode"] {
                if let Some(mode) = values.get(field) {
                    validate_drawing_enum(field, mode, (EnumKind::Placement, true, false))?;
                }
            }
            for field in ["x", "y"] {
                if let Some(coordinate) = values.get(field) {
                    let valid_group = match coordinate {
                        NativeValue::Dict(group) => {
                            group.keys().all(|key| matches!(key.as_str(), "kind" | "name" | "side"))
                                && matches!(group.get("kind"), Some(NativeValue::String(kind)) if kind == "group")
                                && matches!(group.get("name"), Some(NativeValue::String(_)))
                                && group.get("side").is_none_or(|side| {
                                    matches!(side, NativeValue::None)
                                        || matches!(side, NativeValue::String(side) if matches!(side.as_str(), "+" | "-" | "positive" | "negative"))
                                })
                        }
                        _ => false,
                    };
                    if !matches!(coordinate, NativeValue::None | NativeValue::Int(_) | NativeValue::Float(_) | NativeValue::Expression(_))
                        && !valid_group
                    {
                        return Err(PyTypeError::new_err(format!(
                            "placement {field} must be a number, group dictionary, module value, or None"
                        )));
                    }
                }
            }
            if let Some(reference) = values.get("ref") {
                let valid = match reference {
                    NativeValue::None | NativeValue::Expression(_) => true,
                    NativeValue::Int(value) => *value >= 0,
                    _ => false,
                };
                if !valid {
                    return Err(PyTypeError::new_err(
                        "placement ref must be a non-negative node index, module value, or None",
                    ));
                }
            }
            for field in ["dx", "dy"] {
                if let Some(offset) = values.get(field) {
                    if !matches!(offset, NativeValue::None | NativeValue::Int(_) | NativeValue::Float(_) | NativeValue::Expression(_)) {
                        return Err(PyTypeError::new_err(format!(
                            "placement {field} must be a number, module value, or None"
                        )));
                    }
                }
            }
            let pin = matches!(
                values.get("mode"),
                None | Some(NativeValue::Enum(NativeEnum::Placement(PyPlacement::Pin)))
            );
            if pin
                && !["x", "y", "ref", "dx", "dy"]
                    .iter()
                    .any(|field| values.get(*field).is_some_and(|value| !matches!(value, NativeValue::None)))
            {
                return Err(PyValueError::new_err(
                    "pin placement must constrain x, y, ref, dx, or dy",
                ));
            }
            Ok(())
        }
        _ => Err(PyTypeError::new_err(format!(
            "drawing field \"placement\" must be Placement, a typed placement dictionary, a module value, None, or INHERIT, not {}",
            value.kind()
        ))),
    }
}

fn validate_drawing_enum(
    field: &str,
    value: &NativeValue,
    rule: (EnumKind, bool, bool),
) -> PyResult<()> {
    let (expected, allow_none, allow_auto) = rule;
    if matches!(value, NativeValue::Inherit)
        || (allow_none && matches!(value, NativeValue::None))
        || (allow_auto && matches!(value, NativeValue::Auto))
        || matches!(value, NativeValue::Enum(value) if value.kind() == expected)
    {
        Ok(())
    } else {
        Err(PyTypeError::new_err(format!(
            "drawing field {field:?} must be {}, not {}",
            expected.name(),
            value.kind()
        )))
    }
}

fn validate_nested_drawing_enums(value: &NativeValue) -> PyResult<()> {
    match value {
        NativeValue::Array(values) => {
            for value in values {
                validate_nested_drawing_enums(value)?;
            }
        }
        NativeValue::Dict(values)
        | NativeValue::Stroke(values)
        | NativeValue::Insets(values)
        | NativeValue::Mark(values) => {
            for (field, value) in values {
                if let Some(rule) = drawing_enum_rule(field) {
                    validate_drawing_enum(field, value, rule)?;
                }
                validate_nested_drawing_enums(value)?;
            }
        }
        _ => {}
    }
    Ok(())
}

/// Validate finite drawing fields without treating escaped strings as enum values.
pub(crate) fn validate_drawing_field(
    kind: &str,
    field: &str,
    value: &Bound<'_, PyAny>,
) -> PyResult<()> {
    if !matches!(
        kind,
        "node" | "NodeDrawing" | "edge" | "EdgeDrawing" | "half-edge" | "HalfEdgeDrawing"
    ) {
        return Err(PyValueError::new_err(format!(
            "unknown drawing element kind {kind:?}"
        )));
    }
    let native = native_from_py(value, 0)?;
    native.validate(0)?;
    if field == "placement" {
        validate_placement(&native)?;
    } else {
        drawing_field_spec(kind, field)
            .ok_or_else(|| {
                PyValueError::new_err(format!("unknown built-in {kind} drawing field {field:?}"))
            })?
            .validate(&native)?;
    }
    validate_nested_drawing_enums(&native)
}

/// Convert a placement mode shorthand to the graph's native placement dictionary.
pub(crate) fn normalize_placement(py: Python<'_>, value: &Bound<'_, PyAny>) -> PyResult<Py<PyAny>> {
    let value = native_from_py(value, 0)?;
    let value = match value {
        NativeValue::Enum(mode) if mode.kind() == EnumKind::Placement => NativeValue::Dict(
            BTreeMap::from([("mode".to_owned(), NativeValue::Enum(mode))]),
        ),
        value => value,
    };
    native_to_py(py, &value)
}

/// Recursively snapshot a supported native value without retaining mutable aliases.
pub(crate) fn copy_native(py: Python<'_>, value: &Py<PyAny>) -> PyResult<Py<PyAny>> {
    let value = native_from_py(value.bind(py), 0)?;
    value.validate(0)?;
    native_to_py(py, &value)
}

/// Convert a typed Compass value to its DOT keyword.
pub(crate) fn compass_keyword(value: &Bound<'_, PyAny>) -> PyResult<Option<String>> {
    if value.is_none() || value.extract::<PyRef<'_, PyInherit>>().is_ok() {
        return Ok(None);
    }
    value
        .extract::<PyCompass>()
        .map(|value| Some(value.typst_name().to_owned()))
        .map_err(|_| PyTypeError::new_err("compass must be Compass, None, or INHERIT"))
}

/// Convert one validated DOT compass keyword to its typed Compass value.
pub(crate) fn compass_value(py: Python<'_>, value: &str) -> PyResult<Py<PyAny>> {
    let value = match value {
        "n" => PyCompass::N,
        "ne" => PyCompass::NE,
        "e" => PyCompass::E,
        "se" => PyCompass::SE,
        "s" => PyCompass::S,
        "sw" => PyCompass::SW,
        "w" => PyCompass::W,
        "nw" => PyCompass::NW,
        "c" => PyCompass::Center,
        "_" => PyCompass::Any,
        _ => {
            return Err(PyValueError::new_err(format!(
                "unknown DOT compass keyword {value:?}"
            )))
        }
    };
    Ok(Py::new(py, value)?.into_any())
}

/// Encode one closed native value for a reserved DOT drawing attribute.
pub(crate) fn encode_dot_native(_py: Python<'_>, value: &Bound<'_, PyAny>) -> PyResult<String> {
    let value = native_from_py(value, 0)?;
    value.validate(0)?;
    value.reject_dot_module_references()?;
    serde_json::to_string(&DotNativeEnvelope {
        schema: DOT_NATIVE_SCHEMA.to_owned(),
        version: DOT_NATIVE_VERSION,
        value,
    })
    .map_err(|error| PyValueError::new_err(format!("failed to encode Typst value: {error}")))
}

/// Decode a reserved DOT drawing attribute into a closed native Python value.
pub(crate) fn decode_dot_native(py: Python<'_>, text: &str) -> PyResult<Py<PyAny>> {
    let envelope: DotNativeEnvelope = serde_json::from_str(text)
        .map_err(|error| PyValueError::new_err(format!("invalid encoded Typst value: {error}")))?;
    if envelope.schema != DOT_NATIVE_SCHEMA || envelope.version != DOT_NATIVE_VERSION {
        return Err(PyValueError::new_err(format!(
            "unsupported encoded Typst value schema {:?} version {}",
            envelope.schema, envelope.version
        )));
    }
    envelope.value.validate(0)?;
    envelope.value.reject_dot_module_references()?;
    native_to_py(py, &envelope.value)
}

/// Create the sparse default configuration owned by every Python graph.
pub(crate) fn default_render_config(py: Python<'_>) -> PyResult<Py<PyAny>> {
    Ok(Py::new(py, PyRenderConfig::default())?.into_any())
}

/// Reject non-RenderConfig graph defaults at the assignment boundary.
pub(crate) fn validate_render_config(py: Python<'_>, config: &Py<PyAny>) -> PyResult<()> {
    config
        .bind(py)
        .extract::<PyRef<'_, PyRenderConfig>>()
        .map(|_| ())
        .map_err(|_| PyTypeError::new_err("expected RenderConfig"))
}

/// Snapshot a RenderConfig and reject all other objects.
pub(crate) fn render_config_copy(py: Python<'_>, config: &Py<PyAny>) -> PyResult<Py<PyAny>> {
    let config = config
        .bind(py)
        .extract::<PyRef<'_, PyRenderConfig>>()
        .map_err(|_| PyTypeError::new_err("expected RenderConfig"))?
        .clone();
    Ok(Py::new(py, config)?.into_any())
}

fn effective_render_config(
    py: Python<'_>,
    base: &Py<PyAny>,
    overlay: Option<&Bound<'_, PyAny>>,
) -> PyResult<PyRenderConfig> {
    let base = base
        .bind(py)
        .extract::<PyRef<'_, PyRenderConfig>>()
        .map_err(|_| PyTypeError::new_err("base render configuration must be RenderConfig"))?
        .clone();
    let Some(overlay) = overlay else {
        return Ok(base);
    };
    let overlay = overlay
        .extract::<PyRef<'_, PyRenderConfig>>()
        .map_err(|_| PyTypeError::new_err("config override must be RenderConfig"))?;
    Ok(base.merged(&overlay))
}

/// Resolve renderer settings and serialize a base/overlay configuration.
pub(crate) fn render_config_transport(
    py: Python<'_>,
    base: &Py<PyAny>,
    overlay: Option<&Bound<'_, PyAny>>,
    elements: &Bound<'_, PyDict>,
) -> PyResult<RenderConfigTransport> {
    let config = effective_render_config(py, base, overlay)?;
    let elements = native_from_py(elements.as_any(), 0)?;
    elements.validate(0)?;
    let (native_config, imports) = config.native(elements).clinnet_config()?;
    Ok(RenderConfigTransport {
        config: native_config,
        imports,
        template: config.template.resolved(),
        selectors: config.selectors.resolved(py),
    })
}

/// Evaluate one explicit selector and snapshot its closed native result.
pub(crate) fn evaluate_selector(
    py: Python<'_>,
    selector: &Py<PyAny>,
    element: &Bound<'_, PyAny>,
    kind: DrawingKind,
) -> PyResult<Option<Py<PyDict>>> {
    let result = selector.bind(py).call1((element,))?;
    if result.is_none() {
        return Ok(None);
    }
    let expected = match kind {
        DrawingKind::Node => "NodeDrawing",
        DrawingKind::Edge => "EdgeDrawing",
        DrawingKind::HalfEdge => "HalfEdgeDrawing",
    };
    let wrong_type = || {
        PyTypeError::new_err(format!(
            "{} drawing selector must return {expected} or None",
            kind.name(),
        ))
    };
    let values = match kind {
        DrawingKind::Node => result
            .extract::<PyRef<'_, PyNodeDrawing>>()
            .map_err(|_| wrong_type())?
            .clone_values(py)?,
        DrawingKind::Edge => result
            .extract::<PyRef<'_, PyEdgeDrawing>>()
            .map_err(|_| wrong_type())?
            .clone_values(py)?,
        DrawingKind::HalfEdge => result
            .extract::<PyRef<'_, PyHalfEdgeDrawing>>()
            .map_err(|_| wrong_type())?
            .clone_values(py)?,
    };
    Ok(Some(values))
}

/// Register the typed Typst surface and its two singleton sentinels.
pub(crate) fn register_typst_api(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PyAuto>()?;
    module.add_class::<PyInherit>()?;
    module.add_class::<PyLayoutAlgorithm>()?;
    module.add_class::<PyLayoutDirection>()?;
    module.add_class::<PyLabelLayout>()?;
    module.add_class::<PyRankAlignment>()?;
    module.add_class::<PyLayoutNodes>()?;
    module.add_class::<PyPlacement>()?;
    module.add_class::<PyCompass>()?;
    module.add_class::<PyRouting>()?;
    module.add_class::<PyRoutePoints>()?;
    module.add_class::<PyAnchor>()?;
    module.add_class::<PyPattern>()?;
    module.add_class::<PyEdgeLengthResolution>()?;
    module.add_class::<PyDebugLevel>()?;
    module.add_class::<PyStrokeCap>()?;
    module.add_class::<PyStrokeJoin>()?;
    module.add_class::<PyTextStyle>()?;
    module.add_class::<PyDashPattern>()?;
    module.add_class::<PyMarkSymbol>()?;
    module.add_class::<PyMarkPosition>()?;
    module.add_class::<PyMarkOrientation>()?;
    module.add_class::<PyMarkDirection>()?;
    module.add_class::<PyLength>()?;
    module.add_class::<PyRatio>()?;
    module.add_class::<PyRelativeLength>()?;
    module.add_class::<PyAngle>()?;
    module.add_class::<PyFraction>()?;
    module.add_class::<PyColor>()?;
    module.add_class::<PyStroke>()?;
    module.add_class::<PyDash>()?;
    module.add_class::<PyInsets>()?;
    module.add_class::<PyMark>()?;
    module.add_class::<PyTextLabel>()?;
    module.add_class::<PyMathSymbol>()?;
    module.add_class::<PyTypstModule>()?;
    module.add_class::<PyTypstRef>()?;
    module.add_class::<PyTypstCall>()?;
    module.add_class::<PyTypstBind>()?;
    module.add_class::<PyGraphStyleOptions>()?;
    module.add_class::<PyDrawingSelectors>()?;
    module.add_class::<PyLayoutOptions>()?;
    module.add_class::<PyDrawOptions>()?;
    module.add_class::<PyRenderConfig>()?;
    module.add("AUTO", auto_singleton(module.py()))?;
    module.add("INHERIT", inherit_singleton(module.py()))?;
    Ok(())
}
