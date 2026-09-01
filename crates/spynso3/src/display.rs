use std::{cmp::Reverse, collections::BinaryHeap, fmt::Write};

use pyo3::{
    exceptions::{PyImportError, PyRuntimeError, PyValueError},
    prelude::*,
    types::PyDict,
};
use spenso::{
    algebra::complex::RealOrComplexRef,
    iterators::IteratableTensor,
    shadowing::symbolica_utils::SpensoPrintSettings,
    structure::{
        TensorDataLayout,
        partial::PartialStructureExt,
        representation::{IndexRow, RepName, RepresentationClass},
        slot::IsAbstractSlot,
    },
    tensors::{
        complex::RealOrComplexTensor,
        data::{DataTensor, GetTensorData},
        parametric::{AtomViewOrConcrete, ParamOrConcrete},
    },
};
use std::path::Path;
use symbolica::{
    api::python::{PythonExpression, PythonFormattedOutput},
    atom::Atom,
    domains::SelfRing,
    printer::{AnsiHtmlFormatter, PrintOptions, PrintState},
};
use tabled::{
    builder::Builder,
    settings::{Alignment, Style},
};
use tymbolica_atom_payload::{AttachmentSet, encode_atom_render_tree};

use crate::{
    Spensor,
    composition::{self, StructuredAtom},
};

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

const MAX_DISPLAY_ELEMENTS: usize = 100;
const DISPLAY_EDGE: usize = 3;
const DISPLAY_SLICE_EDGE: usize = 1;

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum DisplayIndex {
    Index(usize),
    Ellipsis,
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
    let compile = typst.getattr("compile")?;
    // typst-py treats every value in its virtual-file mapping as UTF-8 source,
    // so binary render trees must be routed through an actual temporary file.
    let temporary_directory = PyModule::import(py, "tempfile")?
        .getattr("TemporaryDirectory")?
        .call0()?;
    let root_string: String = temporary_directory.getattr("name")?.extract()?;
    let root = Path::new(&root_string);
    let write_file = |name: &str, contents: &[u8]| {
        std::fs::write(root.join(name), contents).map_err(|error| {
            PyRuntimeError::new_err(format!("could not prepare Typst input {name}: {error}"))
        })
    };
    write_file("main.typ", main_source.as_bytes())?;
    write_file("render.typ", RENDER_TYP.as_bytes())?;
    write_file(
        "notation.typ",
        notation_source.unwrap_or(NOTATION_TYP).as_bytes(),
    )?;
    if let Some(tree) = tree {
        write_file("tree.cbor", tree)?;
    }
    let kwargs = PyDict::new(py);
    kwargs.set_item("format", format)?;
    kwargs.set_item("pretty", true)?;
    kwargs.set_item("root", &root_string)?;
    let result = compile
        .call((root.join("main.typ"),), Some(&kwargs))
        .and_then(|output| output.extract::<Vec<u8>>());
    let cleanup = temporary_directory.call_method0("cleanup");
    match result {
        Ok(output) => {
            cleanup?;
            Ok(output)
        }
        Err(error) => {
            let _ = cleanup;
            Err(error)
        }
    }
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
    tensor
        .descriptor
        .interface
        .logical_slots()
        .into_iter()
        .map(composition::port_atom)
        .map(|port| format_atom_with_settings(&port, mode, settings))
        .collect::<Vec<_>>()
        .join(" ")
}

fn tensor_value_at_storage(
    tensor: &Spensor,
    index: usize,
) -> Option<AtomViewOrConcrete<'_, RealOrComplexRef<'_, f64>>> {
    let index = index.into();
    match &tensor.tensor {
        ParamOrConcrete::Concrete(RealOrComplexTensor::Real(DataTensor::Dense(data))) => data
            .get_ref_linear(index)
            .map(|value| AtomViewOrConcrete::Concrete(RealOrComplexRef::Real(value))),
        ParamOrConcrete::Concrete(RealOrComplexTensor::Real(DataTensor::Sparse(data))) => {
            Some(AtomViewOrConcrete::Concrete(RealOrComplexRef::Real(
                data.get_ref_linear(index).unwrap_or(&data.zero),
            )))
        }
        ParamOrConcrete::Concrete(RealOrComplexTensor::Complex(DataTensor::Dense(data))) => data
            .get_ref_linear(index)
            .map(|value| AtomViewOrConcrete::Concrete(RealOrComplexRef::Complex(value))),
        ParamOrConcrete::Concrete(RealOrComplexTensor::Complex(DataTensor::Sparse(data))) => {
            Some(AtomViewOrConcrete::Concrete(RealOrComplexRef::Complex(
                data.get_ref_linear(index).unwrap_or(&data.zero),
            )))
        }
        ParamOrConcrete::Param(data) => match &data.tensor {
            DataTensor::Dense(data) => data
                .get_ref_linear(index)
                .map(|value| AtomViewOrConcrete::Atom(value.as_view())),
            DataTensor::Sparse(data) => Some(AtomViewOrConcrete::Atom(
                data.get_ref_linear(index).unwrap_or(&data.zero).as_view(),
            )),
        },
    }
}

fn tensor_sparse_default(
    tensor: &Spensor,
) -> Option<AtomViewOrConcrete<'_, RealOrComplexRef<'_, f64>>> {
    match &tensor.tensor {
        ParamOrConcrete::Concrete(RealOrComplexTensor::Real(DataTensor::Sparse(data))) => Some(
            AtomViewOrConcrete::Concrete(RealOrComplexRef::Real(&data.zero)),
        ),
        ParamOrConcrete::Concrete(RealOrComplexTensor::Complex(DataTensor::Sparse(data))) => Some(
            AtomViewOrConcrete::Concrete(RealOrComplexRef::Complex(&data.zero)),
        ),
        ParamOrConcrete::Param(data) => match &data.tensor {
            DataTensor::Sparse(data) => Some(AtomViewOrConcrete::Atom(data.zero.as_view())),
            DataTensor::Dense(_) => None,
        },
        ParamOrConcrete::Concrete(
            RealOrComplexTensor::Real(DataTensor::Dense(_))
            | RealOrComplexTensor::Complex(DataTensor::Dense(_)),
        ) => None,
    }
}

fn tensor_is_sparse(tensor: &Spensor) -> bool {
    tensor_sparse_default(tensor).is_some()
}

struct ConcreteTensorView<'a> {
    tensor: &'a Spensor,
    layout: TensorDataLayout,
    mode: TensorDisplayMode,
    settings: &'a DisplaySettings,
}

impl<'a> ConcreteTensorView<'a> {
    fn new(
        tensor: &'a Spensor,
        mode: TensorDisplayMode,
        settings: &'a DisplaySettings,
    ) -> Option<Self> {
        Some(Self {
            tensor,
            layout: crate::tensor_data_layout(&tensor.descriptor.interface).ok()?,
            mode,
            settings,
        })
    }

    fn value(&self, logical: &[usize]) -> String {
        self.layout
            .logical_expanded_to_storage_flat(logical)
            .ok()
            .and_then(|storage| tensor_value_at_storage(self.tensor, storage))
            .map(|value| format_tensor_value(value, self.mode, self.settings))
            .unwrap_or_else(|| "?".to_string())
    }

    fn truncated(&self) -> bool {
        self.layout.size() > MAX_DISPLAY_ELEMENTS
    }
}

struct SparseEntry {
    logical: Vec<usize>,
    value: String,
}

struct SparsePreview {
    entries: Vec<SparseEntry>,
    stored: usize,
    truncated: bool,
}

fn sparse_preview(view: &ConcreteTensorView<'_>) -> SparsePreview {
    let edge = MAX_DISPLAY_ELEMENTS / 2;
    let mut first = BinaryHeap::with_capacity(edge + 1);
    let mut last = BinaryHeap::with_capacity(edge + 1);
    let mut stored = 0;

    for (storage, _) in view.tensor.tensor.iter_flat() {
        let storage = usize::from(storage);
        if let Ok(logical) = view.layout.storage_flat_to_logical_flat(storage) {
            let entry = (logical, storage);
            stored += 1;
            first.push(entry);
            if first.len() > edge {
                first.pop();
            }
            last.push(Reverse(entry));
            if last.len() > edge {
                last.pop();
            }
        }
    }

    let mut selected = first
        .into_iter()
        .chain(last.into_iter().map(|Reverse(entry)| entry))
        .collect::<Vec<_>>();
    selected.sort_unstable();
    selected.dedup();
    let entries = selected
        .into_iter()
        .filter_map(|(logical_flat, storage)| {
            Some(SparseEntry {
                logical: expanded_from_flat(logical_flat, view.layout.logical_shape()),
                value: tensor_value_at_storage(view.tensor, storage)
                    .map(|value| format_tensor_value(value, view.mode, view.settings))?,
            })
        })
        .collect::<Vec<_>>();

    SparsePreview {
        entries,
        stored,
        truncated: stored > MAX_DISPLAY_ELEMENTS,
    }
}

fn sparse_displayed_indices(preview: &SparsePreview) -> Vec<DisplayIndex> {
    let mut visible = (0..preview.entries.len())
        .map(DisplayIndex::Index)
        .collect::<Vec<_>>();
    if preview.truncated {
        visible.insert(
            (MAX_DISPLAY_ELEMENTS / 2).min(visible.len()),
            DisplayIndex::Ellipsis,
        );
    }
    visible
}

fn displayed_indices(length: usize, edge: usize, truncate: bool) -> Vec<DisplayIndex> {
    if !truncate || length <= edge.saturating_mul(2) {
        return (0..length).map(DisplayIndex::Index).collect();
    }

    (0..edge)
        .map(DisplayIndex::Index)
        .chain(std::iter::once(DisplayIndex::Ellipsis))
        .chain((length - edge..length).map(DisplayIndex::Index))
        .collect()
}

fn expanded_from_flat(mut flat: usize, shape: &[usize]) -> Vec<usize> {
    let mut expanded = vec![0; shape.len()];
    for (axis, &dimension) in shape.iter().enumerate().rev() {
        if dimension > 0 {
            expanded[axis] = flat % dimension;
            flat /= dimension;
        }
    }
    expanded
}

fn logical_shape_plain(shape: &[usize]) -> String {
    match shape {
        [] => "()".to_string(),
        [dimension] => format!("({dimension},)"),
        _ => format!(
            "({})",
            shape
                .iter()
                .map(usize::to_string)
                .collect::<Vec<_>>()
                .join(", ")
        ),
    }
}

fn logical_shape_typst(shape: &[usize]) -> String {
    shape
        .iter()
        .map(usize::to_string)
        .collect::<Vec<_>>()
        .join(" times ")
}

fn logical_shape_latex(shape: &[usize]) -> String {
    shape
        .iter()
        .map(usize::to_string)
        .collect::<Vec<_>>()
        .join(r"\times")
}

fn slice_label(prefix: &[usize]) -> String {
    let mut positions = prefix.iter().map(usize::to_string).collect::<Vec<_>>();
    positions.extend([":".to_string(), ":".to_string()]);
    format!("[{}]", positions.join(", "))
}

fn plain_vector(view: &ConcreteTensorView<'_>) -> String {
    let [dimension] = view.layout.logical_shape() else {
        return "[?]".to_string();
    };
    let values = displayed_indices(*dimension, DISPLAY_EDGE, view.truncated())
        .into_iter()
        .map(|index| match index {
            DisplayIndex::Index(index) => view.value(&[index]),
            DisplayIndex::Ellipsis => "…".to_string(),
        })
        .collect::<Vec<_>>();
    format!("[{}]", values.join(", "))
}

fn matrix_indices(view: &ConcreteTensorView<'_>) -> (Vec<DisplayIndex>, Vec<DisplayIndex>) {
    let shape = view.layout.logical_shape();
    let rows = shape[shape.len() - 2];
    let columns = shape[shape.len() - 1];
    (
        displayed_indices(rows, DISPLAY_EDGE, view.truncated()),
        displayed_indices(columns, DISPLAY_EDGE, view.truncated()),
    )
}

fn plain_matrix(view: &ConcreteTensorView<'_>, prefix: &[usize]) -> String {
    let (rows, columns) = matrix_indices(view);
    if columns.is_empty() {
        return rows
            .into_iter()
            .map(|row| match row {
                DisplayIndex::Index(_) => "[]",
                DisplayIndex::Ellipsis => "⋮",
            })
            .collect::<Vec<_>>()
            .join("\n");
    }

    let mut table = Builder::new();
    for row in rows {
        let values = match row {
            DisplayIndex::Ellipsis => vec!["⋮".to_string(); columns.len()],
            DisplayIndex::Index(row) => columns
                .iter()
                .map(|column| match column {
                    DisplayIndex::Index(column) => {
                        let mut logical = prefix.to_vec();
                        logical.extend([row, *column]);
                        view.value(&logical)
                    }
                    DisplayIndex::Ellipsis => "…".to_string(),
                })
                .collect(),
        };
        table.push_record(values);
    }

    let mut table = table.build();
    table.with(Style::blank()).with(Alignment::right());
    table
        .to_string()
        .lines()
        .map(|row| format!("[{row}]"))
        .collect::<Vec<_>>()
        .join("\n")
}

fn typst_vector(view: &ConcreteTensorView<'_>) -> String {
    let [dimension] = view.layout.logical_shape() else {
        return "vec(?)".to_string();
    };
    let values = displayed_indices(*dimension, DISPLAY_EDGE, view.truncated())
        .into_iter()
        .map(|index| match index {
            DisplayIndex::Index(index) => view.value(&[index]),
            DisplayIndex::Ellipsis => "dots.v".to_string(),
        })
        .collect::<Vec<_>>();
    format!("vec({})", values.join(","))
}

fn typst_matrix(view: &ConcreteTensorView<'_>, prefix: &[usize]) -> String {
    let (rows, columns) = matrix_indices(view);
    if rows.is_empty() || columns.is_empty() {
        return "mat()".to_string();
    }
    let rows = rows
        .into_iter()
        .map(|row| match row {
            DisplayIndex::Ellipsis => columns
                .iter()
                .map(|_| "dots.v")
                .collect::<Vec<_>>()
                .join(","),
            DisplayIndex::Index(row) => columns
                .iter()
                .map(|column| match column {
                    DisplayIndex::Index(column) => {
                        let mut logical = prefix.to_vec();
                        logical.extend([row, *column]);
                        view.value(&logical)
                    }
                    DisplayIndex::Ellipsis => "dots.h".to_string(),
                })
                .collect::<Vec<_>>()
                .join(","),
        })
        .collect::<Vec<_>>();
    format!("mat({})", rows.join(";"))
}

fn latex_vector(view: &ConcreteTensorView<'_>) -> String {
    let [dimension] = view.layout.logical_shape() else {
        return "?".to_string();
    };
    let values = displayed_indices(*dimension, DISPLAY_EDGE, view.truncated())
        .into_iter()
        .map(|index| match index {
            DisplayIndex::Index(index) => view.value(&[index]),
            DisplayIndex::Ellipsis => r"\vdots".to_string(),
        })
        .collect::<Vec<_>>();
    format!(r"\begin{{pmatrix}}{}\end{{pmatrix}}", values.join(r" \\ "))
}

fn latex_matrix(view: &ConcreteTensorView<'_>, prefix: &[usize]) -> String {
    let (rows, columns) = matrix_indices(view);
    if rows.is_empty() || columns.is_empty() {
        return r"\begin{pmatrix}\end{pmatrix}".to_string();
    }
    let rows = rows
        .into_iter()
        .map(|row| match row {
            DisplayIndex::Ellipsis => columns
                .iter()
                .map(|_| r"\vdots")
                .collect::<Vec<_>>()
                .join(" & "),
            DisplayIndex::Index(row) => columns
                .iter()
                .map(|column| match column {
                    DisplayIndex::Index(column) => {
                        let mut logical = prefix.to_vec();
                        logical.extend([row, *column]);
                        view.value(&logical)
                    }
                    DisplayIndex::Ellipsis => r"\cdots".to_string(),
                })
                .collect::<Vec<_>>()
                .join(" & "),
        })
        .collect::<Vec<_>>();
    format!(r"\begin{{pmatrix}}{}\end{{pmatrix}}", rows.join(r" \\ "))
}

fn high_rank_body(view: &ConcreteTensorView<'_>) -> String {
    let shape = view.layout.logical_shape();
    let prefix_shape = &shape[..shape.len() - 2];
    let slice_count = prefix_shape.iter().product();
    let slices = displayed_indices(slice_count, DISPLAY_SLICE_EDGE, view.truncated());

    match view.mode {
        TensorDisplayMode::Plain => slices
            .into_iter()
            .map(|slice| match slice {
                DisplayIndex::Ellipsis => "…".to_string(),
                DisplayIndex::Index(slice) => {
                    let prefix = expanded_from_flat(slice, prefix_shape);
                    format!("{}\n{}", slice_label(&prefix), plain_matrix(view, &prefix))
                }
            })
            .collect::<Vec<_>>()
            .join("\n\n"),
        TensorDisplayMode::Typst => {
            let slices = slices
                .into_iter()
                .map(|slice| match slice {
                    DisplayIndex::Ellipsis => "dots.v".to_string(),
                    DisplayIndex::Index(slice) => {
                        let prefix = expanded_from_flat(slice, prefix_shape);
                        format!(
                            r#"attach({},t:op("{}"))"#,
                            typst_matrix(view, &prefix),
                            slice_label(&prefix)
                        )
                    }
                })
                .collect::<Vec<_>>();
            format!(
                r#"op("Tensor")_({})({})"#,
                logical_shape_typst(shape),
                slices.join(",")
            )
        }
        TensorDisplayMode::Latex => {
            let slices = slices
                .into_iter()
                .map(|slice| match slice {
                    DisplayIndex::Ellipsis => r"\vdots".to_string(),
                    DisplayIndex::Index(slice) => {
                        let prefix = expanded_from_flat(slice, prefix_shape);
                        format!(
                            r"\text{{slice }}{} \\ {}",
                            slice_label(&prefix),
                            latex_matrix(view, &prefix)
                        )
                    }
                })
                .collect::<Vec<_>>();
            format!(
                r"\operatorname{{Tensor}}_{{{}}}\left(\begin{{array}}{{c}}{}\end{{array}}\right)",
                logical_shape_latex(shape),
                slices.join(r" \\[0.6em] ")
            )
        }
    }
}

fn sparse_body(view: &ConcreteTensorView<'_>) -> String {
    let preview = sparse_preview(view);
    let entries = &preview.entries;
    let default = tensor_sparse_default(view.tensor)
        .map(|value| format_tensor_value(value, view.mode, view.settings))
        .unwrap_or_else(|| "0".to_string());
    let visible = sparse_displayed_indices(&preview);
    let shape = view.layout.logical_shape();

    match view.mode {
        TensorDisplayMode::Plain => {
            let stored = preview.stored;
            let formatted = visible
                .into_iter()
                .map(|entry| match entry {
                    DisplayIndex::Ellipsis => "  …".to_string(),
                    DisplayIndex::Index(entry) => format!(
                        "  [{}]: {}",
                        entries[entry]
                            .logical
                            .iter()
                            .map(usize::to_string)
                            .collect::<Vec<_>>()
                            .join(", "),
                        entries[entry].value
                    ),
                })
                .collect::<Vec<_>>();
            format!(
                "Sparse(shape={}, stored={}, default={default}) {{\n{}\n}}",
                logical_shape_plain(shape),
                stored,
                formatted.join("\n")
            )
        }
        TensorDisplayMode::Typst => {
            let mut formatted = vec![format!(r#"op("default")={default}"#)];
            formatted.extend(visible.into_iter().map(|entry| match entry {
                DisplayIndex::Ellipsis => "dots.v".to_string(),
                DisplayIndex::Index(entry) => format!(
                    r#"op("[{}]")={}"#,
                    entries[entry]
                        .logical
                        .iter()
                        .map(usize::to_string)
                        .collect::<Vec<_>>()
                        .join(","),
                    entries[entry].value
                ),
            }));
            format!(
                r#"op("SparseTensor")_({})({})"#,
                logical_shape_typst(shape),
                formatted.join(",")
            )
        }
        TensorDisplayMode::Latex => {
            let formatted = visible
                .into_iter()
                .map(|entry| match entry {
                    DisplayIndex::Ellipsis => r"\vdots".to_string(),
                    DisplayIndex::Index(entry) => format!(
                        r"[{}]&\mapsto&{}",
                        entries[entry]
                            .logical
                            .iter()
                            .map(usize::to_string)
                            .collect::<Vec<_>>()
                            .join(","),
                        entries[entry].value
                    ),
                })
                .collect::<Vec<_>>();
            format!(
                r"\operatorname{{SparseTensor}}_{{{}}}\left\{{\begin{{array}}{{rcl}}{}\end{{array}};\ {default}\ \text{{otherwise}}\right.",
                logical_shape_latex(shape),
                formatted.join(r" \\ ")
            )
        }
    }
}

fn format_tensor_interface_rows(
    tensor: &Spensor,
    settings: &DisplaySettings,
) -> (Vec<String>, Vec<String>) {
    let mut top = Vec::new();
    let mut bottom = Vec::new();
    for slot in tensor.descriptor.interface.logical_slots() {
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
        let port = composition::port_atom(slot);
        let source = format_atom_with_settings(&port, TensorDisplayMode::Typst, settings);
        match row {
            IndexRow::Top => top.push(source),
            IndexRow::Bottom => bottom.push(source),
        }
    }
    (top, bottom)
}

fn concrete_tensor_body(
    tensor: &Spensor,
    mode: TensorDisplayMode,
    settings: &DisplaySettings,
) -> String {
    let Some(view) = ConcreteTensorView::new(tensor, mode, settings) else {
        return "Tensor(<invalid logical layout>)".to_string();
    };

    if tensor_is_sparse(tensor) && view.truncated() {
        return sparse_body(&view);
    }

    match (mode, view.layout.logical_shape()) {
        (_, []) => view.value(&[]),
        (TensorDisplayMode::Plain, [_]) => plain_vector(&view),
        (TensorDisplayMode::Plain, [_, _]) => plain_matrix(&view, &[]),
        (TensorDisplayMode::Typst, [_]) => typst_vector(&view),
        (TensorDisplayMode::Typst, [_, _]) => typst_matrix(&view, &[]),
        (TensorDisplayMode::Latex, [_]) => latex_vector(&view),
        (TensorDisplayMode::Latex, [_, _]) => latex_matrix(&view, &[]),
        (_, _) => high_rank_body(&view),
    }
}

fn escape_html(value: &str) -> String {
    AnsiHtmlFormatter::escape_html(value)
}

fn html_table(view: &ConcreteTensorView<'_>, prefix: &[usize], caption: Option<&str>) -> String {
    let (rows, columns) = matrix_indices(view);
    let mut html = String::from(
        r#"<table style="border-collapse:collapse;font-family:ui-monospace,SFMono-Regular,Menlo,Consolas,monospace;font-size:.92em">"#,
    );
    if let Some(caption) = caption {
        let _ = write!(
            html,
            r#"<caption style="caption-side:top;text-align:left;padding:0 0 .25em;color:inherit"><code>{}</code></caption>"#,
            escape_html(caption)
        );
    }
    html.push_str("<tbody>");
    if rows.is_empty() || columns.is_empty() {
        html.push_str(
            r#"<tr><td style="padding:.25em .55em;border:1px solid rgba(127,127,127,.35)"><code>∅</code></td></tr>"#,
        );
    } else {
        for row in rows {
            match row {
                DisplayIndex::Ellipsis => {
                    let _ = write!(
                        html,
                        r#"<tr><td colspan="{}" style="padding:.1em .55em;text-align:center;border:1px solid rgba(127,127,127,.35)"><code>⋮</code></td></tr>"#,
                        columns.len()
                    );
                }
                DisplayIndex::Index(row) => {
                    html.push_str("<tr>");
                    for column in &columns {
                        let value = match column {
                            DisplayIndex::Index(column) => {
                                let mut logical = prefix.to_vec();
                                logical.extend([row, *column]);
                                escape_html(&view.value(&logical))
                            }
                            DisplayIndex::Ellipsis => "&hellip;".to_string(),
                        };
                        let _ = write!(
                            html,
                            r#"<td style="padding:.25em .55em;text-align:right;border:1px solid rgba(127,127,127,.35)"><code>{value}</code></td>"#
                        );
                    }
                    html.push_str("</tr>");
                }
            }
        }
    }
    html.push_str("</tbody></table>");
    html
}

fn html_vector(view: &ConcreteTensorView<'_>) -> String {
    let [dimension] = view.layout.logical_shape() else {
        return "<code>?</code>".to_string();
    };
    let mut html = String::from(
        r#"<table style="border-collapse:collapse;font-family:ui-monospace,SFMono-Regular,Menlo,Consolas,monospace;font-size:.92em"><tbody><tr>"#,
    );
    for index in displayed_indices(*dimension, DISPLAY_EDGE, view.truncated()) {
        let value = match index {
            DisplayIndex::Index(index) => escape_html(&view.value(&[index])),
            DisplayIndex::Ellipsis => "&hellip;".to_string(),
        };
        let _ = write!(
            html,
            r#"<td style="padding:.25em .55em;text-align:right;border:1px solid rgba(127,127,127,.35)"><code>{value}</code></td>"#
        );
    }
    html.push_str("</tr></tbody></table>");
    html
}

fn html_high_rank(view: &ConcreteTensorView<'_>) -> String {
    let shape = view.layout.logical_shape();
    let prefix_shape = &shape[..shape.len() - 2];
    let slice_count = prefix_shape.iter().product();
    let mut html = String::from(
        r#"<div style="display:flex;flex-wrap:wrap;align-items:flex-start;gap:.8em 1.2em">"#,
    );
    for slice in displayed_indices(slice_count, DISPLAY_SLICE_EDGE, view.truncated()) {
        match slice {
            DisplayIndex::Ellipsis => html.push_str(
                r#"<div style="align-self:center;padding:.5em"><code>&hellip;</code></div>"#,
            ),
            DisplayIndex::Index(slice) => {
                let prefix = expanded_from_flat(slice, prefix_shape);
                html.push_str(&html_table(view, &prefix, Some(&slice_label(&prefix))));
            }
        }
    }
    html.push_str("</div>");
    html
}

fn html_sparse(view: &ConcreteTensorView<'_>) -> String {
    let preview = sparse_preview(view);
    let entries = &preview.entries;
    let default = tensor_sparse_default(view.tensor)
        .map(|value| {
            escape_html(&format_tensor_value(
                value,
                TensorDisplayMode::Plain,
                view.settings,
            ))
        })
        .unwrap_or_else(|| "0".to_string());
    let visible = sparse_displayed_indices(&preview);
    let mut html = format!(
        r#"<div style="margin-bottom:.35em;opacity:.75">stored={} &middot; default=<code>{default}</code></div><table style="border-collapse:collapse;font-family:ui-monospace,SFMono-Regular,Menlo,Consolas,monospace;font-size:.92em"><thead><tr><th scope="col" style="padding:.2em .55em;text-align:left;border-bottom:1px solid rgba(127,127,127,.5)">logical index</th><th scope="col" style="padding:.2em .55em;text-align:right;border-bottom:1px solid rgba(127,127,127,.5)">value</th></tr></thead><tbody>"#,
        preview.stored
    );
    for entry in visible {
        match entry {
            DisplayIndex::Ellipsis => html.push_str(
                r#"<tr><td colspan="2" style="padding:.1em .55em;text-align:center"><code>&hellip;</code></td></tr>"#,
            ),
            DisplayIndex::Index(entry) => {
                let coordinate = format!(
                    "[{}]",
                    entries[entry]
                        .logical
                        .iter()
                        .map(usize::to_string)
                        .collect::<Vec<_>>()
                        .join(", ")
                );
                let _ = write!(
                    html,
                    r#"<tr><td style="padding:.2em .55em;text-align:left;border-bottom:1px solid rgba(127,127,127,.2)"><code>{}</code></td><td style="padding:.2em .55em;text-align:right;border-bottom:1px solid rgba(127,127,127,.2)"><code>{}</code></td></tr>"#,
                    escape_html(&coordinate),
                    escape_html(&entries[entry].value)
                );
            }
        }
    }
    html.push_str("</tbody></table>");
    html
}

#[allow(dead_code)]
fn concrete_tensor_to_table_html(tensor: &Spensor, settings: &DisplaySettings) -> String {
    let Some(view) = ConcreteTensorView::new(tensor, TensorDisplayMode::Plain, settings) else {
        return r#"<div data-spenso-tensor><code>Tensor(&lt;invalid logical layout&gt;)</code></div>"#
            .to_string();
    };
    let interface = escape_html(&format_tensor_interface(
        tensor,
        TensorDisplayMode::Plain,
        settings,
    ));
    let shape = escape_html(&logical_shape_plain(view.layout.logical_shape()));
    let mut html = String::from(
        r#"<div data-spenso-tensor style="display:inline-block;max-width:100%;color:inherit">"#,
    );
    html.push_str(r#"<div style="margin-bottom:.45em"><strong>Tensor</strong>"#);
    if !interface.is_empty() {
        let _ = write!(html, "<sub>({interface})</sub>");
    }
    let _ = write!(
        html,
        r#" <span style="opacity:.65">shape={shape}</span></div>"#
    );

    if tensor_is_sparse(tensor) && view.truncated() {
        html.push_str(&html_sparse(&view));
    } else {
        match view.layout.logical_shape() {
            [] => {
                let _ = write!(
                    html,
                    "<div><code>{}</code></div>",
                    escape_html(&view.value(&[]))
                );
            }
            [_] => html.push_str(&html_vector(&view)),
            [_, _] => html.push_str(&html_table(&view, &[], None)),
            _ => html.push_str(&html_high_rank(&view)),
        }
    }
    html.push_str("</div>");
    html
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
        return match mode {
            TensorDisplayMode::Latex => format!("$${body}$$"),
            TensorDisplayMode::Plain | TensorDisplayMode::Typst => body,
        };
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
        network::{parsing::ShadowedStructure, tags::SPENSO_TAG},
        structure::{
            OrderedStructure,
            abstract_index::AbstractIndex,
            concrete_index::FlatIndex,
            dimension::Dimension,
            partial::{PartialIndex, PartialStructure, PartialStructureExt},
            representation::{ExtendibleReps, LibraryRep, RepName, Representation},
            slot::IsAbstractSlot,
        },
        tensors::data::{DenseTensor, SparseTensor},
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
    fn display_selection_keeps_both_edges() {
        assert_eq!(
            displayed_indices(10, 2, true),
            vec![
                DisplayIndex::Index(0),
                DisplayIndex::Index(1),
                DisplayIndex::Ellipsis,
                DisplayIndex::Index(8),
                DisplayIndex::Index(9),
            ]
        );
        assert_eq!(
            displayed_indices(4, 2, true),
            (0..4).map(DisplayIndex::Index).collect::<Vec<_>>()
        );
    }

    #[test]
    fn leading_slice_coordinates_are_row_major() {
        assert_eq!(expanded_from_flat(0, &[2, 3]), vec![0, 0]);
        assert_eq!(expanded_from_flat(4, &[2, 3]), vec![1, 1]);
        assert_eq!(slice_label(&[1, 1]), "[1, 1, :, :]");
    }

    #[test]
    fn tensor_html_escapes_cell_content() {
        assert_eq!(escape_html("x < y & z"), "x &lt; y &amp; z");
    }

    #[test]
    fn scalar_latex_is_a_complete_notebook_math_fragment() {
        let interface = PartialStructure::from_logical_slots([]);
        let tensor = Spensor::scalar_with_descriptor(
            2.,
            StructuredAtom::new(Atom::Zero, interface),
            None,
            Vec::new(),
        );

        let settings = DisplaySettings::default();
        assert_eq!(
            format_concrete_with_settings(&tensor, TensorDisplayMode::Latex, &settings),
            "$$2$$"
        );
    }

    #[test]
    fn sparse_preview_keeps_bounded_logical_edges() {
        let euc = ExtendibleReps::EUCLIDEAN.new_rep(Dimension::Concrete(1_000));
        let interface = PartialStructure::from_logical_slots([
            euc.slot(PartialIndex::Explicit(AbstractIndex::Normal(0)))
        ]);
        let structure = OrderedStructure::new(
            interface
                .logical_slots()
                .into_iter()
                .map(|slot| slot.rep().slot(AbstractIndex::Normal(0)))
                .collect(),
        )
        .map_canonical(|structure| ShadowedStructure {
            structure,
            global_name: None,
            additional_args: None,
        })
        .into_canonical();
        let tensor = Spensor::from_storage_with_descriptor(
            SparseTensor {
                elements: (0usize..=200)
                    .map(|index| (FlatIndex::from(index), index as f64))
                    .collect(),
                zero: -1.,
                structure,
            }
            .into(),
            StructuredAtom::new(Atom::Zero, interface),
            None,
            Vec::new(),
        );
        let settings = DisplaySettings::default();
        let view = ConcreteTensorView::new(&tensor, TensorDisplayMode::Plain, &settings).unwrap();
        let preview = sparse_preview(&view);

        assert_eq!(preview.stored, 201);
        assert!(preview.truncated);
        assert_eq!(preview.entries.len(), MAX_DISPLAY_ELEMENTS);
        assert_eq!(
            preview
                .entries
                .iter()
                .map(|entry| entry.logical[0])
                .collect::<Vec<_>>(),
            (0..50).chain(151..=200).collect::<Vec<_>>()
        );
    }

    #[test]
    fn mixed_representation_tensor_displays_in_logical_order() {
        let mink = ExtendibleReps::MINKOWSKI.new_rep(Dimension::Concrete(2));
        let euc = ExtendibleReps::EUCLIDEAN.new_rep(Dimension::Concrete(3));
        let interface = PartialStructure::from_logical_slots([
            mink.slot(PartialIndex::Explicit(AbstractIndex::Normal(0))),
            euc.slot(PartialIndex::Explicit(AbstractIndex::Normal(1))),
        ]);
        let tensor = OrderedStructure::new(
            interface
                .logical_slots()
                .into_iter()
                .enumerate()
                .map(|(index, slot)| slot.rep().slot(AbstractIndex::Normal(index)))
                .collect(),
        )
        .map_canonical(|structure| ShadowedStructure {
            structure,
            global_name: None,
            additional_args: None,
        })
        .map_canonical(|structure| {
            DenseTensor::from_storage_data(vec![1., 4000., 20., 50., 300., 6.], structure)
                .unwrap()
                .into()
        })
        .into_canonical();
        let tensor = Spensor::from_storage_with_descriptor(
            tensor,
            StructuredAtom::new(Atom::Zero, interface),
            None,
            Vec::new(),
        );

        let settings = DisplaySettings::default();
        let plain = format_concrete_tensor(&tensor, false);
        assert!(
            plain.contains("[    1   20   300 ]\n[ 4000   50     6 ]"),
            "unexpected plain tensor display:\n{plain}"
        );
        assert!(concrete_tensor_to_typst(&tensor, &settings).contains("mat(1,20,300;4000,50,6)"));
        let html = concrete_tensor_to_table_html(&tensor, &settings);
        let positions = [1, 20, 300, 4000, 50, 6]
            .into_iter()
            .map(|value| html.find(&format!("<code>{value}</code>")).unwrap())
            .collect::<Vec<_>>();
        assert!(positions.windows(2).all(|pair| pair[0] < pair[1]));
    }

    #[test]
    fn concrete_tensor_interface_uses_descriptor_indices() {
        let euc = ExtendibleReps::EUCLIDEAN.new_rep(Dimension::Concrete(2));
        let interface = PartialStructure::from_logical_slots([
            euc.slot(PartialIndex::Explicit(AbstractIndex::from(symbol!("row")))),
            euc.slot(PartialIndex::Explicit(AbstractIndex::from(symbol!(
                "column"
            )))),
        ]);
        let storage = OrderedStructure::new(vec![
            euc.slot(AbstractIndex::Open {
                owner: 101,
                axis: 0,
            }),
            euc.slot(AbstractIndex::Open {
                owner: 101,
                axis: 1,
            }),
        ])
        .map_canonical(|structure| ShadowedStructure {
            structure,
            global_name: None,
            additional_args: None,
        })
        .map_canonical(|structure| {
            DenseTensor::from_storage_data(vec![1., 2., 3., 4.], structure)
                .unwrap()
                .into()
        })
        .into_canonical();
        let tensor = Spensor::from_storage_with_descriptor(
            storage,
            StructuredAtom::new(Atom::Zero, interface),
            None,
            Vec::new(),
        );

        let settings = DisplaySettings::default();
        for rendered in [
            format_concrete_tensor(&tensor, false),
            concrete_tensor_to_typst(&tensor, &settings),
            concrete_tensor_to_table_html(&tensor, &settings),
            format_concrete_with_settings(&tensor, TensorDisplayMode::Latex, &settings),
        ] {
            let row = rendered.find("row").unwrap();
            let column = rendered.find("column").unwrap();
            assert!(row < column, "unexpected descriptor order: {rendered}");
            assert!(
                !rendered.contains("open") && !rendered.contains("101"),
                "storage index leaked: {rendered}"
            );
        }
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
        let spaced = DisplaySettings {
            index_gap: "0.2em".to_owned(),
            ..DisplaySettings::default()
        };
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
