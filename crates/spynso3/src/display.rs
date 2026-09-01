use pyo3::{
    exceptions::{PyImportError, PyRuntimeError, PyValueError},
    prelude::*,
    types::{PyBytes, PyDict},
};
use spenso::{
    algebra::complex::RealOrComplexRef,
    iterators::IteratableTensor,
    shadowing::symbolica_utils::SpensoPrintSettings,
    structure::{
        HasStructure, TensorStructure,
        abstract_index::AbstractIndex,
        permuted::PermuteTensor,
        representation::{IndexRow, LibraryRep, RepName, RepresentationClass},
        slot::{IsAbstractSlot, Slot},
    },
    tensors::{data::SparseOrDense, parametric::AtomViewOrConcrete},
};
use symbolica::{
    api::python::{PythonExpression, PythonFormattedOutput},
    atom::Atom,
    domains::SelfRing,
    printer::{PrintOptions, PrintState},
};
use tymbolica_atom_payload::{AttachmentSet, encode_atom_render_tree};

use crate::{Spensor, composition::StructuredAtom};

#[cfg(feature = "python_stubgen")]
use pyo3_stub_gen::derive::{gen_stub_pyclass, gen_stub_pyfunction};

const RENDER_TYP: &str = include_str!("../typst/render.typ");
const NOTATION_TYP: &str = include_str!("../typst/notation.typ");

/// Presentation settings shared by Typst source, HTML, and SVG rendering.
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    frozen,
    from_py_object,
    name = "DisplaySettings",
    module = "symbolica.community.spenso"
)]
#[derive(Clone, Debug, PartialEq)]
pub struct DisplaySettings {
    #[pyo3(get)]
    tensor_layout: String,
    #[pyo3(get)]
    show_dimensions: bool,
    #[pyo3(get)]
    parentheses: bool,
    #[pyo3(get)]
    commas: Option<bool>,
    #[pyo3(get)]
    symbol_scripts: bool,
    #[pyo3(get)]
    index_gap: String,
    #[pyo3(get)]
    factor_gap: String,
}

impl Default for DisplaySettings {
    fn default() -> Self {
        Self {
            tensor_layout: "ports".to_owned(),
            show_dimensions: false,
            parentheses: true,
            commas: None,
            symbol_scripts: true,
            index_gap: "0.08em".to_owned(),
            factor_gap: "0.12em".to_owned(),
        }
    }
}

fn validate_typst_length(value: &str, field: &str) -> PyResult<()> {
    let number = ["pt", "mm", "cm", "in", "em", "%"]
        .into_iter()
        .find_map(|unit| value.strip_suffix(unit));
    let valid_number = number
        .filter(|number| !number.is_empty())
        .and_then(|number| number.parse::<f64>().ok())
        .is_some_and(f64::is_finite);
    if !valid_number {
        return Err(PyValueError::new_err(format!(
            "{field} must be a simple Typst length such as 0.08em"
        )));
    }
    Ok(())
}

#[pymethods]
impl DisplaySettings {
    #[new]
    #[pyo3(signature = (
        tensor_layout = "ports",
        show_dimensions = false,
        parentheses = true,
        commas = None,
        symbol_scripts = true,
        index_gap = "0.08em",
        factor_gap = "0.12em",
    ))]
    fn new(
        tensor_layout: &str,
        show_dimensions: bool,
        parentheses: bool,
        commas: Option<bool>,
        symbol_scripts: bool,
        index_gap: &str,
        factor_gap: &str,
    ) -> PyResult<Self> {
        if !matches!(tensor_layout, "ports" | "schoonschip" | "call") {
            return Err(PyValueError::new_err(
                "tensor_layout must be 'ports', 'schoonschip', or 'call'",
            ));
        }
        validate_typst_length(index_gap, "index_gap")?;
        validate_typst_length(factor_gap, "factor_gap")?;
        Ok(Self {
            tensor_layout: tensor_layout.to_owned(),
            show_dimensions,
            parentheses,
            commas,
            symbol_scripts,
            index_gap: index_gap.to_owned(),
            factor_gap: factor_gap.to_owned(),
        })
    }

    #[staticmethod]
    fn ports() -> Self {
        Self::default()
    }

    #[staticmethod]
    fn schoonschip() -> Self {
        Self {
            tensor_layout: "schoonschip".to_owned(),
            ..Self::default()
        }
    }

    #[staticmethod]
    fn call() -> Self {
        Self {
            tensor_layout: "call".to_owned(),
            commas: Some(true),
            ..Self::default()
        }
    }

    fn __repr__(&self) -> String {
        format!(
            "DisplaySettings(tensor_layout={:?}, show_dimensions={})",
            self.tensor_layout, self.show_dimensions
        )
    }
}

pub(crate) fn resolved_settings(
    show_dimensions: Option<bool>,
    settings: Option<&DisplaySettings>,
) -> DisplaySettings {
    let mut resolved = settings.cloned().unwrap_or_default();
    if let Some(show_dimensions) = show_dimensions {
        resolved.show_dimensions = show_dimensions;
    }
    resolved
}

fn has_renderer_only_settings(settings: &DisplaySettings) -> bool {
    let defaults = DisplaySettings::default();
    settings.tensor_layout != defaults.tensor_layout
        || settings.index_gap != defaults.index_gap
        || settings.factor_gap != defaults.factor_gap
}

fn reject_renderer_only_settings(settings: &DisplaySettings, method: &str) -> PyResult<()> {
    if has_renderer_only_settings(settings) {
        return Err(PyValueError::new_err(format!(
            "{method} emits ports-style source; Schoonschip, call, and custom spacing layouts require to_html or to_svg"
        )));
    }
    Ok(())
}

pub(crate) fn validate_plain_source_settings(settings: &DisplaySettings) -> PyResult<()> {
    reject_renderer_only_settings(settings, "format_tensor")
}

pub(crate) fn validate_typst_source_settings(settings: &DisplaySettings) -> PyResult<()> {
    reject_renderer_only_settings(settings, "to_typst")
}

#[derive(Clone, Copy)]
enum TensorDisplayMode {
    Plain,
    Latex,
    Typst,
}

fn display_options(mode: TensorDisplayMode, show_dimensions: bool) -> PrintOptions {
    let mut presentation = match mode {
        TensorDisplayMode::Plain => SpensoPrintSettings::compact(),
        TensorDisplayMode::Latex | TensorDisplayMode::Typst => SpensoPrintSettings::typst(),
    };
    presentation.with_dim = show_dimensions;

    let mut options = match mode {
        TensorDisplayMode::Plain => presentation.nice_symbolica(),
        TensorDisplayMode::Latex => PrintOptions {
            custom_print_mode: presentation.into(),
            ..PrintOptions::latex()
        },
        TensorDisplayMode::Typst => PrintOptions {
            custom_print_mode: presentation.into(),
            ..PrintOptions::typst()
        },
    };
    options.color_builtin_symbols = false;
    options.color_top_level_sum = false;
    options.terms_on_new_line = false;
    options.max_line_length = None;
    options
}

fn display_options_with_settings(
    mode: TensorDisplayMode,
    settings: &DisplaySettings,
) -> PrintOptions {
    let mut presentation = match mode {
        TensorDisplayMode::Plain => SpensoPrintSettings::compact(),
        TensorDisplayMode::Latex | TensorDisplayMode::Typst => SpensoPrintSettings::typst(),
    };
    presentation.with_dim = settings.show_dimensions;
    presentation.parens = settings.parentheses;
    presentation.commas = settings.commas.unwrap_or(settings.tensor_layout == "call");
    presentation.symbol_scripts = settings.symbol_scripts;

    let mut options = match mode {
        TensorDisplayMode::Plain => presentation.nice_symbolica(),
        TensorDisplayMode::Latex => PrintOptions {
            custom_print_mode: presentation.into(),
            ..PrintOptions::latex()
        },
        TensorDisplayMode::Typst => PrintOptions {
            custom_print_mode: presentation.into(),
            ..PrintOptions::typst()
        },
    };
    options.color_builtin_symbols = false;
    options.color_top_level_sum = false;
    options.terms_on_new_line = false;
    options.max_line_length = None;
    options
}

fn format_atom_with_settings(
    atom: &Atom,
    mode: TensorDisplayMode,
    settings: &DisplaySettings,
) -> String {
    atom.format_string(
        &display_options_with_settings(mode, settings),
        PrintState::new(),
    )
}

fn format_atom_with_mode(atom: &Atom, mode: TensorDisplayMode, show_dimensions: bool) -> String {
    atom.format_string(&display_options(mode, show_dimensions), PrintState::new())
}

fn format_structured_with_mode(
    value: &StructuredAtom,
    mode: TensorDisplayMode,
    show_dimensions: bool,
) -> String {
    format_atom_with_mode(&value.presentation_atom(), mode, show_dimensions)
}

fn format_structured_with_settings(
    value: &StructuredAtom,
    mode: TensorDisplayMode,
    settings: &DisplaySettings,
) -> String {
    format_atom_with_settings(&value.presentation_atom(), mode, settings)
}

pub(crate) fn structured_to_typst_with_settings(
    value: &StructuredAtom,
    settings: &DisplaySettings,
) -> String {
    format_structured_with_settings(value, TensorDisplayMode::Typst, settings)
}

pub(crate) fn format_structured_settings(
    value: &StructuredAtom,
    settings: &DisplaySettings,
) -> String {
    format_structured_with_settings(value, TensorDisplayMode::Plain, settings)
}

#[cfg(test)]
fn format_atom(atom: &Atom, show_dimensions: bool) -> String {
    format_atom_with_mode(atom, TensorDisplayMode::Plain, show_dimensions)
}

pub(crate) fn format_structured(value: &StructuredAtom, show_dimensions: bool) -> String {
    format_structured_with_mode(value, TensorDisplayMode::Plain, show_dimensions)
}

pub(crate) fn structured_to_latex(value: &StructuredAtom, show_dimensions: bool) -> String {
    let body = format_structured_with_mode(value, TensorDisplayMode::Latex, show_dimensions);
    format!("$${body}$$")
}

// Typst supplies semantic HTML when the optional renderer is installed;
// Symbolica's LaTeX remains the notebook fallback otherwise.
pub(crate) fn format_atom_output_rich(
    py: Python<'_>,
    atom: &Atom,
    settings: &DisplaySettings,
    notation_source: Option<&str>,
) -> PythonFormattedOutput {
    PythonFormattedOutput {
        text: format_atom_with_settings(atom, TensorDisplayMode::Plain, settings),
        html: atom_to_html(py, atom, settings, notation_source).ok(),
        latex: Some({
            let body = format_atom_with_settings(atom, TensorDisplayMode::Latex, settings);
            format!("$${body}$$")
        }),
    }
}

pub(crate) fn format_structured_output_rich(
    py: Python<'_>,
    value: &StructuredAtom,
    settings: &DisplaySettings,
    notation_source: Option<&str>,
) -> PythonFormattedOutput {
    format_atom_output_rich(py, &value.presentation_atom(), settings, notation_source)
}

fn portable_attachments(atom: &Atom) -> Result<AttachmentSet, String> {
    spenso::portable_payload::attachments_for_atom(atom).map_err(|error| error.to_string())
}

fn typst_settings_source(settings: &DisplaySettings) -> String {
    let commas = match settings.commas {
        Some(value) => value.to_string(),
        None => "none".to_owned(),
    };
    format!(
        concat!(
            "(\n",
            "  tensor-layout: {:?},\n",
            "  with-dim: {},\n",
            "  parens: {},\n",
            "  commas: {},\n",
            "  symbol-scripts: {},\n",
            "  index-gap: {},\n",
            "  factor-gap: {},\n",
            ")"
        ),
        settings.tensor_layout,
        settings.show_dimensions,
        settings.parentheses,
        commas,
        settings.symbol_scripts,
        settings.index_gap,
        settings.factor_gap,
    )
}

fn typst_main_source(settings: &DisplaySettings) -> String {
    format!(
        concat!(
            "#import \"notation.typ\" as tensor-notation\n",
            "#set page(width: auto, height: auto, margin: 4pt)\n",
            "#let tree = cbor(read(\"tree.cbor\", encoding: none))\n",
            "#let settings = {}\n",
            "#let visual = tensor-notation.render(\n",
            "  tree,\n",
            "  notation: tensor-notation.default-notation(settings: settings),\n",
            ")\n",
            "$ #visual $\n"
        ),
        typst_settings_source(settings),
    )
}

fn render_atom_tree(atom: &Atom) -> PyResult<Vec<u8>> {
    let attachments = portable_attachments(atom).map_err(PyRuntimeError::new_err)?;
    encode_atom_render_tree(atom, &attachments)
        .map_err(|error| PyRuntimeError::new_err(error.to_string()))
}

fn compile_typst(
    py: Python<'_>,
    main_source: &str,
    format: &str,
    notation_source: Option<&str>,
    tree: Option<&[u8]>,
) -> PyResult<Vec<u8>> {
    let typst = PyModule::import(py, "typst").map_err(|_| {
        PyImportError::new_err(
            "Typst rendering requires the optional dependency; install gammaloop[typst-display]",
        )
    })?;
    let files = PyDict::new(py);
    files.set_item("main.typ", PyBytes::new(py, main_source.as_bytes()))?;
    files.set_item("render.typ", PyBytes::new(py, RENDER_TYP.as_bytes()))?;
    files.set_item(
        "notation.typ",
        PyBytes::new(py, notation_source.unwrap_or(NOTATION_TYP).as_bytes()),
    )?;
    if let Some(tree) = tree {
        files.set_item("tree.cbor", PyBytes::new(py, tree))?;
    }
    let kwargs = PyDict::new(py);
    kwargs.set_item("format", format)?;
    kwargs.set_item("pretty", true)?;
    typst
        .getattr("compile")?
        .call((files,), Some(&kwargs))?
        .extract::<Vec<u8>>()
}

fn extract_html_fragment(document: &str) -> Result<String, &'static str> {
    let body_start = document
        .find("<body")
        .and_then(|start| document[start..].find('>').map(|offset| start + offset + 1))
        .ok_or("Typst HTML output has no body")?;
    let body_end = document[body_start..]
        .find("</body>")
        .map(|offset| body_start + offset)
        .ok_or("Typst HTML output has no closing body")?;

    let mut fragment = String::new();
    let mut rest = document;
    while let Some(style_start) = rest.find("<style") {
        let style = &rest[style_start..];
        let Some(style_end) = style.find("</style>") else {
            break;
        };
        fragment.push_str(&style[..style_end + "</style>".len()]);
        rest = &style[style_end + "</style>".len()..];
    }
    fragment.push_str(document[body_start..body_end].trim());
    Ok(fragment)
}

pub(crate) fn atom_to_html(
    py: Python<'_>,
    atom: &Atom,
    settings: &DisplaySettings,
    notation_source: Option<&str>,
) -> PyResult<String> {
    let tree = render_atom_tree(atom)?;
    let source = typst_main_source(settings);
    let html = compile_typst(py, &source, "html", notation_source, Some(&tree))?;
    let html = String::from_utf8(html).map_err(|error| {
        PyRuntimeError::new_err(format!("Typst returned invalid UTF-8: {error}"))
    })?;
    extract_html_fragment(&html).map_err(PyRuntimeError::new_err)
}

pub(crate) fn atom_to_svg(
    py: Python<'_>,
    atom: &Atom,
    settings: &DisplaySettings,
    notation_source: Option<&str>,
) -> PyResult<String> {
    let tree = render_atom_tree(atom)?;
    let source = typst_main_source(settings);
    let svg = compile_typst(py, &source, "svg", notation_source, Some(&tree))?;
    String::from_utf8(svg)
        .map_err(|error| PyRuntimeError::new_err(format!("Typst returned invalid UTF-8: {error}")))
}

pub(crate) fn structured_to_html(
    py: Python<'_>,
    value: &StructuredAtom,
    settings: &DisplaySettings,
    notation_source: Option<&str>,
) -> PyResult<String> {
    atom_to_html(py, &value.presentation_atom(), settings, notation_source)
}

pub(crate) fn structured_to_svg(
    py: Python<'_>,
    value: &StructuredAtom,
    settings: &DisplaySettings,
    notation_source: Option<&str>,
) -> PyResult<String> {
    atom_to_svg(py, &value.presentation_atom(), settings, notation_source)
}

fn format_tensor_value(
    value: AtomViewOrConcrete<'_, RealOrComplexRef<'_, f64>>,
    mode: TensorDisplayMode,
    settings: &DisplaySettings,
) -> String {
    match value {
        AtomViewOrConcrete::Atom(atom) => {
            format_atom_with_settings(&atom.to_owned(), mode, settings)
        }
        AtomViewOrConcrete::Concrete(RealOrComplexRef::Real(value)) => value.to_string(),
        AtomViewOrConcrete::Concrete(RealOrComplexRef::Complex(value)) => value.to_string(),
    }
}

fn format_tensor_interface(
    tensor: &Spensor,
    mode: TensorDisplayMode,
    settings: &DisplaySettings,
) -> String {
    logical_tensor_slots(tensor)
        .into_iter()
        .map(|slot| format_atom_with_settings(&slot.to_atom(), mode, settings))
        .collect::<Vec<_>>()
        .join(" ")
}

fn format_tensor_interface_rows(
    tensor: &Spensor,
    settings: &DisplaySettings,
) -> (Vec<String>, Vec<String>) {
    let mut top = Vec::new();
    let mut bottom = Vec::new();
    for slot in logical_tensor_slots(tensor) {
        let representation = slot.rep_name();
        let row = representation
            .metadata()
            .map(|metadata| {
                if metadata.class == RepresentationClass::Dualizable && representation.is_dual() {
                    metadata.index_row.opposite()
                } else {
                    metadata.index_row
                }
            })
            .unwrap_or(IndexRow::Top);
        let source = format_atom_with_settings(&slot.to_atom(), TensorDisplayMode::Typst, settings);
        match row {
            IndexRow::Top => top.push(source),
            IndexRow::Bottom => bottom.push(source),
        }
    }
    (top, bottom)
}

fn logical_tensor_slots(tensor: &Spensor) -> Vec<Slot<LibraryRep, AbstractIndex>> {
    let canonical = tensor
        .tensor
        .structure
        .structure()
        .external_structure_iter()
        .collect::<Vec<_>>();
    let representation_sorted = tensor.tensor.index_permutation.apply_slice_inv(&canonical);
    tensor
        .tensor
        .rep_permutation
        .apply_slice_inv(&representation_sorted)
}

fn concrete_tensor_body(
    tensor: &Spensor,
    mode: TensorDisplayMode,
    settings: &DisplaySettings,
) -> String {
    let dense = tensor
        .tensor
        .structure
        .clone()
        .to_dense()
        .permute_inds(&tensor.tensor.index_permutation.inverse())
        .permute_inds(&tensor.tensor.rep_permutation.inverse());
    let values = dense
        .iter_flat()
        .map(|(_, value)| format_tensor_value(value, mode, settings))
        .collect::<Vec<_>>();
    let shape = logical_tensor_slots(tensor)
        .into_iter()
        .map(|slot| usize::try_from(slot.dim()))
        .collect::<Result<Vec<_>, _>>()
        .ok();

    match (mode, shape.as_deref()) {
        (_, Some([])) => values.first().cloned().unwrap_or_else(|| "0".to_string()),
        (TensorDisplayMode::Plain, Some([_])) => format!("[{}]", values.join(", ")),
        (TensorDisplayMode::Plain, Some([_, columns])) if *columns > 0 => values
            .chunks(*columns)
            .map(|row| format!("[{}]", row.join(", ")))
            .collect::<Vec<_>>()
            .join("\n"),
        (TensorDisplayMode::Plain, _) => format!("[{}]", values.join(", ")),
        (TensorDisplayMode::Typst, Some([_])) => format!("vec({})", values.join(",")),
        (TensorDisplayMode::Typst, Some([_, columns])) if *columns > 0 => {
            let rows = values
                .chunks(*columns)
                .map(|row| row.join(","))
                .collect::<Vec<_>>();
            format!("mat({})", rows.join(";"))
        }
        (TensorDisplayMode::Typst, _) => {
            format!(r#"op("Tensor")({})"#, values.join(","))
        }
        (TensorDisplayMode::Latex, Some([_])) => {
            format!(
                r#"\begin{{pmatrix}}{}\end{{pmatrix}}"#,
                values.join(r#" \\ "#)
            )
        }
        (TensorDisplayMode::Latex, Some([_, columns])) if *columns > 0 => {
            let rows = values
                .chunks(*columns)
                .map(|row| row.join(" & "))
                .collect::<Vec<_>>();
            format!(
                r#"\begin{{pmatrix}}{}\end{{pmatrix}}"#,
                rows.join(r#" \\ "#)
            )
        }
        (TensorDisplayMode::Latex, _) => {
            format!(
                r#"\operatorname{{Tensor}}\left({}\right)"#,
                values.join(",")
            )
        }
    }
}

fn concrete_render_project(
    tensor: &Spensor,
    settings: &DisplaySettings,
) -> PyResult<(String, Vec<u8>)> {
    let descriptor = tensor.descriptor.presentation_atom();
    let tree = render_atom_tree(&descriptor)?;
    let body = concrete_tensor_body(tensor, TensorDisplayMode::Typst, settings);
    let source = format!(
        concat!(
            "#import \"notation.typ\" as tensor-notation\n",
            "#set page(width: auto, height: auto, margin: 4pt)\n",
            "#let tree = cbor(read(\"tree.cbor\", encoding: none))\n",
            "#let settings = {}\n",
            "#let descriptor = tensor-notation.render(\n",
            "  tree,\n",
            "  notation: tensor-notation.default-notation(settings: settings),\n",
            ")\n",
            "$ #descriptor = {} $\n"
        ),
        typst_settings_source(settings),
        body,
    );
    Ok((source, tree))
}

fn format_concrete_with_settings(
    tensor: &Spensor,
    mode: TensorDisplayMode,
    settings: &DisplaySettings,
) -> String {
    let body = concrete_tensor_body(tensor, mode, settings);
    if matches!(mode, TensorDisplayMode::Typst) {
        let (top, bottom) = format_tensor_interface_rows(tensor, settings);
        let mut rows = Vec::with_capacity(2);
        if !top.is_empty() {
            rows.push(format!("t:({})", top.join(" ")));
        }
        if !bottom.is_empty() {
            rows.push(format!("b:({})", bottom.join(" ")));
        }
        return if rows.is_empty() {
            body
        } else {
            format!("attach({body},{})", rows.join(","))
        };
    }

    let interface = format_tensor_interface(tensor, mode, settings);
    if interface.is_empty() {
        return body;
    }

    match mode {
        TensorDisplayMode::Plain => format!("Tensor_({interface})\n{body}"),
        TensorDisplayMode::Typst => unreachable!("Typst rendering returned above"),
        TensorDisplayMode::Latex => format!("$${body}_{{{interface}}}$$"),
    }
}

pub(crate) fn format_concrete_tensor(tensor: &Spensor, show_dimensions: bool) -> String {
    let settings = resolved_settings(Some(show_dimensions), None);
    format_concrete_with_settings(tensor, TensorDisplayMode::Plain, &settings)
}

pub(crate) fn format_concrete_tensor_with_settings(
    tensor: &Spensor,
    settings: &DisplaySettings,
) -> String {
    format_concrete_with_settings(tensor, TensorDisplayMode::Plain, settings)
}

pub(crate) fn concrete_tensor_to_typst(tensor: &Spensor, settings: &DisplaySettings) -> String {
    format_concrete_with_settings(tensor, TensorDisplayMode::Typst, settings)
}

pub(crate) fn concrete_tensor_to_html(
    py: Python<'_>,
    tensor: &Spensor,
    settings: &DisplaySettings,
    notation_source: Option<&str>,
) -> PyResult<String> {
    let (source, tree) = concrete_render_project(tensor, settings)?;
    let html = compile_typst(py, &source, "html", notation_source, Some(&tree))?;
    let html = String::from_utf8(html).map_err(|error| {
        PyRuntimeError::new_err(format!("Typst returned invalid UTF-8: {error}"))
    })?;
    extract_html_fragment(&html).map_err(PyRuntimeError::new_err)
}

pub(crate) fn concrete_tensor_to_svg(
    py: Python<'_>,
    tensor: &Spensor,
    settings: &DisplaySettings,
    notation_source: Option<&str>,
) -> PyResult<String> {
    let (source, tree) = concrete_render_project(tensor, settings)?;
    let svg = compile_typst(py, &source, "svg", notation_source, Some(&tree))?;
    String::from_utf8(svg)
        .map_err(|error| PyRuntimeError::new_err(format!("Typst returned invalid UTF-8: {error}")))
}

pub(crate) fn format_concrete_tensor_output(
    tensor: &Spensor,
    show_dimensions: bool,
) -> PythonFormattedOutput {
    let settings = resolved_settings(Some(show_dimensions), None);
    PythonFormattedOutput {
        text: format_concrete_with_settings(tensor, TensorDisplayMode::Plain, &settings),
        html: None,
        latex: Some(format_concrete_with_settings(
            tensor,
            TensorDisplayMode::Latex,
            &settings,
        )),
    }
}

pub(crate) fn format_concrete_tensor_output_rich(
    py: Python<'_>,
    tensor: &Spensor,
    settings: &DisplaySettings,
    notation_source: Option<&str>,
) -> PythonFormattedOutput {
    PythonFormattedOutput {
        text: format_concrete_with_settings(tensor, TensorDisplayMode::Plain, settings),
        html: concrete_tensor_to_html(py, tensor, settings, notation_source).ok(),
        latex: Some(format_concrete_with_settings(
            tensor,
            TensorDisplayMode::Latex,
            settings,
        )),
    }
}

/// Format a tensor expression using compact Spenso notation.
#[cfg_attr(
    feature = "python_stubgen",
    gen_stub_pyfunction(module = "symbolica.community.spenso")
)]
#[pyfunction]
#[pyo3(signature = (expression, show_dimensions = None, *, settings = None))]
fn format_tensor(
    expression: &PythonExpression,
    show_dimensions: Option<bool>,
    settings: Option<PyRef<'_, DisplaySettings>>,
) -> PyResult<String> {
    let settings = resolved_settings(show_dimensions, settings.as_deref());
    validate_plain_source_settings(&settings)?;
    Ok(format_atom_with_settings(
        &expression.expr,
        TensorDisplayMode::Plain,
        &settings,
    ))
}

/// Format a tensor expression as Typst math source.
#[cfg_attr(
    feature = "python_stubgen",
    gen_stub_pyfunction(module = "symbolica.community.spenso")
)]
#[pyfunction]
#[pyo3(signature = (expression, show_dimensions = None, *, settings = None))]
fn to_typst(
    expression: &PythonExpression,
    show_dimensions: Option<bool>,
    settings: Option<PyRef<'_, DisplaySettings>>,
) -> PyResult<String> {
    let settings = resolved_settings(show_dimensions, settings.as_deref());
    validate_typst_source_settings(&settings)?;
    Ok(format_atom_with_settings(
        &expression.expr,
        TensorDisplayMode::Typst,
        &settings,
    ))
}

/// Render a tensor expression to semantic HTML through the optional Typst runtime.
#[cfg_attr(
    feature = "python_stubgen",
    gen_stub_pyfunction(module = "symbolica.community.spenso")
)]
#[pyfunction]
#[pyo3(signature = (expression, show_dimensions = None, *, settings = None, notation_source = None))]
fn to_html(
    py: Python<'_>,
    expression: &PythonExpression,
    show_dimensions: Option<bool>,
    settings: Option<PyRef<'_, DisplaySettings>>,
    notation_source: Option<String>,
) -> PyResult<String> {
    let settings = resolved_settings(show_dimensions, settings.as_deref());
    atom_to_html(py, &expression.expr, &settings, notation_source.as_deref())
}

/// Render a tensor expression to an SVG string through the optional Typst runtime.
#[cfg_attr(
    feature = "python_stubgen",
    gen_stub_pyfunction(module = "symbolica.community.spenso")
)]
#[pyfunction]
#[pyo3(signature = (expression, show_dimensions = None, *, settings = None, notation_source = None))]
fn to_svg(
    py: Python<'_>,
    expression: &PythonExpression,
    show_dimensions: Option<bool>,
    settings: Option<PyRef<'_, DisplaySettings>>,
    notation_source: Option<String>,
) -> PyResult<String> {
    let settings = resolved_settings(show_dimensions, settings.as_deref());
    atom_to_svg(py, &expression.expr, &settings, notation_source.as_deref())
}

/// Build Symbolica's rich display wrapper for a tensor expression.
#[cfg_attr(
    feature = "python_stubgen",
    gen_stub_pyfunction(module = "symbolica.community.spenso")
)]
#[pyfunction]
#[pyo3(signature = (expression, show_dimensions = None, *, settings = None, notation_source = None))]
fn formatted(
    py: Python<'_>,
    expression: &PythonExpression,
    show_dimensions: Option<bool>,
    settings: Option<PyRef<'_, DisplaySettings>>,
    notation_source: Option<String>,
) -> PythonFormattedOutput {
    let settings = resolved_settings(show_dimensions, settings.as_deref());
    format_atom_output_rich(py, &expression.expr, &settings, notation_source.as_deref())
}

pub(crate) fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<DisplaySettings>()?;
    m.add_function(wrap_pyfunction!(format_tensor, m)?)?;
    m.add_function(wrap_pyfunction!(to_typst, m)?)?;
    m.add_function(wrap_pyfunction!(to_html, m)?)?;
    m.add_function(wrap_pyfunction!(to_svg, m)?)?;
    m.add_function(wrap_pyfunction!(formatted, m)?)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use spenso::{
        network::tags::SPENSO_TAG,
        structure::{
            abstract_index::AbstractIndex,
            dimension::Dimension,
            representation::{LibraryRep, Representation},
        },
        vector_symbol,
    };
    use symbolica::{atom::AtomCore, function, symbol};

    fn vector() -> Atom {
        let representation = Representation {
            rep: LibraryRep::SelfDual(1),
            dim: Dimension::Concrete(4),
        };
        function!(
            vector_symbol!("display_test_vector"),
            Atom::num(1),
            representation
                .slot::<AbstractIndex, _>(AbstractIndex::from(symbol!("mu")))
                .to_atom()
        )
    }

    #[test]
    fn dimensions_are_opt_in_and_formatting_does_not_change_the_atom() {
        let atom = vector();
        let canonical = atom.to_canonical_string();

        let compact = format_atom(&atom, false);
        let dimensioned = format_atom(&atom, true);

        assert_ne!(compact, dimensioned);
        assert!(!compact.contains('4'));
        assert!(dimensioned.contains('4'));
        assert_eq!(atom.to_canonical_string(), canonical);
    }

    #[test]
    fn typst_automatically_uses_spenso_printers() {
        let atom = vector();
        let automatic = atom.printer(PrintOptions::typst()).to_string();

        assert!(!automatic.contains("spenso::"));
        assert!(symbol!("display_test_vector").has_tag(&SPENSO_TAG.rank1));
    }

    #[test]
    fn document_settings_cover_each_tensor_layout() {
        for layout in ["ports", "schoonschip", "call"] {
            let settings = DisplaySettings::new(
                layout,
                true,
                false,
                Some(layout == "call"),
                false,
                "0.04em",
                "0.1em",
            )
            .unwrap();
            let source = typst_settings_source(&settings);
            assert!(source.contains(&format!("tensor-layout: {layout:?}")));
            assert!(source.contains("with-dim: true"));
            assert!(source.contains("parens: false"));
        }
    }

    #[test]
    fn spacing_requires_a_finite_number_and_supported_typst_unit() {
        for valid in ["0pt", "0.08em", "-1.5mm", ".25%"] {
            assert!(validate_typst_length(valid, "gap").is_ok(), "{valid}");
        }
        for invalid in ["", "10", "1garbage", "NaNem", "infpt"] {
            assert!(validate_typst_length(invalid, "gap").is_err(), "{invalid}");
        }
    }

    #[test]
    fn standalone_typst_source_rejects_renderer_only_settings() {
        assert!(validate_typst_source_settings(&DisplaySettings::default()).is_ok());
        assert!(validate_plain_source_settings(&DisplaySettings::default()).is_ok());
        assert!(validate_typst_source_settings(&DisplaySettings::schoonschip()).is_err());
        assert!(validate_plain_source_settings(&DisplaySettings::call()).is_err());
        let mut spaced = DisplaySettings::default();
        spaced.index_gap = "0.2em".to_owned();
        assert!(validate_typst_source_settings(&spaced).is_err());
    }

    #[test]
    fn generated_project_reads_the_portable_render_tree_as_binary() {
        let source = typst_main_source(&DisplaySettings::default());
        assert!(source.contains("cbor(read(\"tree.cbor\", encoding: none))"));
        assert!(source.contains("tensor-notation.render"));
    }

    #[test]
    fn html_fragment_keeps_styles_and_body_mathml() {
        let html = "<html><head><style>.x{color:red}</style></head><body><math><mi>x</mi></math></body></html>";
        assert_eq!(
            extract_html_fragment(html).unwrap(),
            "<style>.x{color:red}</style><math><mi>x</mi></math>"
        );
    }
}
