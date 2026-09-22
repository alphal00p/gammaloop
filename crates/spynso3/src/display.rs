use std::{
    cmp::Reverse,
    collections::{BTreeMap, BinaryHeap, HashMap, HashSet},
    fmt::Write,
};

use pyo3::{
    exceptions::{PyImportError, PyRuntimeError, PyValueError},
    prelude::*,
    types::PyDict,
};
use spenso::{
    algebra::complex::RealOrComplexRef,
    iterators::IteratableTensor,
    network::tags::SPENSO_TAG,
    portable_payload::register_math_display_symbol,
    shadowing::symbolica_utils::SpensoPrintSettings,
    structure::{
        TensorDataLayout,
        partial::PartialStructureExt,
        representation::{
            IndexDisplay, IndexPalette, IndexRow, RepName, RepresentationClass,
            RepresentationMetadata,
        },
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
    atom::{Atom, AtomCore, AtomView, FunctionBuilder, Symbol},
    domains::SelfRing,
    printer::{AnsiHtmlFormatter, PrintOptions, PrintState},
};
use symbolica_typst_atom_payload::{AttachmentSet, encode_atom_render_tree};
use tabled::{
    builder::Builder,
    settings::{Alignment, Style},
};

use crate::{
    Spensor,
    composition::{self, StructuredAtom},
};

#[cfg(feature = "python_stubgen")]
use pyo3_stub_gen::derive::{gen_stub_pyclass, gen_stub_pyfunction, gen_stub_pymethods};

const RENDER_TYP: &str = include_str!("../typst/render.typ");
const NOTATION_TYP: &str = include_str!("../typst/notation.typ");
const NOTEBOOK_STYLE: &str = include_str!("../typst/notebook.css");

/// Presentation settings shared by Typst source, HTML, and SVG rendering.
///
/// ``index_style="alphabet"`` assigns representation-specific letters to graph
/// indices within each displayed expression. ``"graph"`` keeps graph identities
/// in their subscripts; ``"raw"`` preserves the original index notation. All
/// styles leave the underlying expression and tensor interface unchanged.
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
    index_style: String,
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
            index_style: "alphabet".to_owned(),
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

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
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
        index_style = "alphabet",
    ))]
    #[allow(clippy::too_many_arguments)] // Each Python display setting is independently optional.
    fn new(
        tensor_layout: &str,
        show_dimensions: bool,
        parentheses: bool,
        commas: Option<bool>,
        symbol_scripts: bool,
        index_gap: &str,
        factor_gap: &str,
        index_style: &str,
    ) -> PyResult<Self> {
        if !matches!(tensor_layout, "ports" | "schoonschip" | "call") {
            return Err(PyValueError::new_err(
                "tensor_layout must be 'ports', 'schoonschip', or 'call'",
            ));
        }
        if !matches!(index_style, "alphabet" | "graph" | "raw") {
            return Err(PyValueError::new_err(
                "index_style must be 'alphabet', 'graph', or 'raw'",
            ));
        }
        validate_typst_length(index_gap, "index_gap")?;
        validate_typst_length(factor_gap, "factor_gap")?;
        Ok(Self {
            tensor_layout: tensor_layout.to_owned(),
            index_style: index_style.to_owned(),
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
            "DisplaySettings(tensor_layout={:?}, show_dimensions={}, index_style={:?})",
            self.tensor_layout, self.show_dimensions, self.index_style
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

/// Aliases belong to one display, never to the algebraic index namespace.
/// Rich renderers receive these as notation records over the original Atom tree.
#[derive(Default)]
struct IndexAliases {
    entries: HashMap<(Symbol, Atom), (usize, IndexDisplay)>,
}

impl IndexAliases {
    fn for_descriptor(descriptor: &StructuredAtom, style: &str) -> Self {
        let mut bundle = FunctionBuilder::new(spenso::structure::abstract_index::AIND_SYMBOLS.aind)
            .add_arg(&descriptor.atom);
        for slot in descriptor.interface.logical_slots() {
            bundle = bundle.add_arg(composition::port_atom(slot));
        }
        Self::for_atom(&bundle.finish(), style)
    }

    fn for_atom(atom: &Atom, style: &str) -> Self {
        if style == "raw" {
            return Self::default();
        }
        let mut slots = Vec::new();
        let _ = atom.replace_map(|value, _, output| {
            let AtomView::Fun(function) = value else {
                return;
            };
            if function.get_symbol().is_scalar() {
                **output = value.to_owned();
            } else if function.get_symbol().has_tag(&SPENSO_TAG.representation)
                && function.get_nargs() == 2
            {
                slots.push((
                    function.get_symbol(),
                    function.iter().nth(1).unwrap().to_owned(),
                ));
                **output = value.to_owned();
            }
        });
        let mut occupied: HashMap<Symbol, HashSet<String>> = HashMap::new();
        let mut pending = BTreeMap::new();
        for (representation, index) in slots {
            let Some(metadata) = RepresentationMetadata::from_symbol(representation) else {
                continue;
            };
            let IndexPalette::Cyclic { .. } = metadata.index_palette else {
                continue;
            };
            if Self::compound(index.as_view()).is_some() {
                pending.insert(
                    (
                        representation.get_name().to_owned(),
                        index.to_canonical_string(),
                    ),
                    (representation, index, metadata.index_palette),
                );
            } else {
                let label = usize::try_from(index.as_view())
                    .ok()
                    .and_then(|position| metadata.index_palette.resolve(position))
                    .or_else(|| match index.as_view() {
                        AtomView::Var(variable) => IndexDisplay::from_symbol(variable.get_symbol())
                            .or_else(|| {
                                IndexDisplay::symbol(variable.get_symbol().get_stripped_name()).ok()
                            }),
                        _ => None,
                    });
                if let Some(label) = label {
                    occupied
                        .entry(representation)
                        .or_default()
                        .insert(Self::label_key(&label));
                }
            }
        }
        let mut aliases = Self::default();
        for (_, (representation, index, palette)) in pending {
            let IndexPalette::Cyclic { start, .. } = &palette else {
                unreachable!()
            };
            let mut position = *start;
            let display = if style == "graph" {
                let (head, arguments) = Self::compound(index.as_view()).unwrap();
                let name = match head.get_name() {
                    "gammalooprs::hedge" => "h",
                    "gammalooprs::edge" => "e",
                    "gammalooprs::vertex" => "v",
                    name => name,
                };
                let suffix = match (name, arguments.as_slice()) {
                    ("h" | "e" | "v", [owner, 1]) => format!("{name}{owner}"),
                    ("h" | "e" | "v", [owner, local]) => format!("{name}{owner}.{local}"),
                    (_, arguments) => format!(
                        "{name}({})",
                        arguments
                            .iter()
                            .map(usize::to_string)
                            .collect::<Vec<_>>()
                            .join(",")
                    ),
                };
                let Some(base) = palette.resolve(*start) else {
                    continue;
                };
                let Ok(suffix) = IndexDisplay::text(suffix) else {
                    continue;
                };
                Some(base.with_bottom(suffix))
            } else {
                let used = occupied.entry(representation).or_default();
                loop {
                    let Some(display) = palette.resolve(position) else {
                        break None;
                    };
                    if used.insert(Self::label_key(&display)) {
                        break Some(display);
                    }
                    let Some(next) = position
                        .checked_add(1)
                        .filter(|next| *next <= i64::MAX as usize)
                    else {
                        break None;
                    };
                    position = next;
                }
            };
            // Exhaustion is possible for a custom palette starting at i64::MAX.
            let Some(display) = display else { continue };
            aliases
                .entries
                .insert((representation, index), (position, display));
        }
        aliases
    }

    fn compound(index: AtomView<'_>) -> Option<(Symbol, Vec<usize>)> {
        let AtomView::Fun(function) = index else {
            return None;
        };
        if !function.get_symbol().has_tag(&SPENSO_TAG.index)
            || !(1..=2).contains(&function.get_nargs())
        {
            return None;
        }
        let arguments = function
            .iter()
            .map(usize::try_from)
            .collect::<Result<Vec<_>, _>>()
            .ok()?;
        Some((function.get_symbol(), arguments))
    }

    fn label_key(display: &IndexDisplay) -> String {
        // Unicode spellings and Typst's named Greek symbols have the same glyph.
        let mut key = display.to_typst_source();
        for (glyph, name) in [
            ('α', "alpha"),
            ('β', "beta"),
            ('γ', "gamma"),
            ('δ', "delta"),
            ('ε', "epsilon"),
            ('ζ', "zeta"),
            ('η', "eta"),
            ('θ', "theta"),
            ('ι', "iota"),
            ('κ', "kappa"),
            ('λ', "lambda"),
            ('μ', "mu"),
            ('ν', "nu"),
            ('ξ', "xi"),
            ('ο', "omicron"),
            ('π', "pi"),
            ('ρ', "rho"),
            ('σ', "sigma"),
            ('τ', "tau"),
            ('υ', "upsilon"),
            ('φ', "phi"),
            ('χ', "chi"),
            ('ψ', "psi"),
            ('ω', "omega"),
        ] {
            key = key.replace(glyph, name);
        }
        key
    }

    fn presentation_atom(&self, atom: &Atom, style: &str) -> Atom {
        if self.entries.is_empty() {
            return atom.clone();
        }
        atom.replace_map(|value, _, output| {
            let AtomView::Fun(function) = value else {
                return;
            };
            if function.get_symbol().is_scalar() {
                **output = value.to_owned();
                return;
            }
            if !function.get_symbol().has_tag(&SPENSO_TAG.representation)
                || function.get_nargs() != 2
            {
                return;
            }
            let mut arguments = function.iter();
            let dimension = arguments.next().unwrap();
            let index = arguments.next().unwrap();
            let Some((position, display)) =
                self.entries.get(&(function.get_symbol(), index.to_owned()))
            else {
                return;
            };
            let replacement = if style == "graph" {
                let Ok(symbol) = register_math_display_symbol(display, "spenso::display") else {
                    return;
                };
                Atom::var(symbol)
            } else {
                Atom::num(*position as i64)
            };
            **output = FunctionBuilder::new(function.get_symbol())
                .add_arg(dimension)
                .add_arg(replacement)
                .finish();
        })
    }

    fn typst_source(&self, style: &str) -> String {
        let mut records = self
            .entries
            .iter()
            .map(|((representation, index), (position, display))| {
                let (head, arguments) = Self::compound(index.as_view()).unwrap();
                let arguments = arguments
                    .iter()
                    .map(|argument| format!("{argument},"))
                    .collect::<String>();
                let label = if style == "graph" {
                    format!("${}$", display.to_typst_source())
                } else {
                    position.to_string()
                };
                format!(
                    "({:?},{:?},({arguments}),{label}),",
                    representation.get_name(),
                    head.get_name()
                )
            })
            .collect::<Vec<_>>();
        records.sort();
        format!("({})", records.join(""))
    }
}

fn format_atom_with_settings(
    atom: &Atom,
    mode: TensorDisplayMode,
    settings: &DisplaySettings,
) -> String {
    let aliases = IndexAliases::for_atom(atom, &settings.index_style);
    let presentation = aliases.presentation_atom(atom, &settings.index_style);
    presentation.format_string(
        &display_options_with_settings(mode, settings),
        PrintState::new(),
    )
}

fn format_atom_with_mode(atom: &Atom, mode: TensorDisplayMode, show_dimensions: bool) -> String {
    let aliases = IndexAliases::for_atom(atom, "alphabet");
    aliases
        .presentation_atom(atom, "alphabet")
        .format_string(&display_options(mode, show_dimensions), PrintState::new())
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

fn typst_settings_source(settings: &DisplaySettings, atom: &Atom) -> String {
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
            "  index-aliases: {},\n",
            ")"
        ),
        settings.tensor_layout,
        settings.show_dimensions,
        settings.parentheses,
        commas,
        settings.symbol_scripts,
        settings.index_gap,
        settings.factor_gap,
        IndexAliases::for_atom(atom, &settings.index_style).typst_source(&settings.index_style),
    )
}

fn typst_main_source(settings: &DisplaySettings, atom: &Atom) -> String {
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
        typst_settings_source(settings, atom),
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
    let source = typst_main_source(settings, atom);
    let html = compile_typst(py, &source, "html", notation_source, Some(&tree))?;
    let html = String::from_utf8(html).map_err(|error| {
        PyRuntimeError::new_err(format!("Typst returned invalid UTF-8: {error}"))
    })?;
    let fragment = extract_html_fragment(&html).map_err(PyRuntimeError::new_err)?;
    // Use the page's math fonts so each standalone fragment stays small.
    Ok(format!(
        "<style>{NOTEBOOK_STYLE}</style><div data-spenso-math>{fragment}</div>"
    ))
}

pub(crate) fn atom_to_svg(
    py: Python<'_>,
    atom: &Atom,
    settings: &DisplaySettings,
    notation_source: Option<&str>,
) -> PyResult<String> {
    let tree = render_atom_tree(atom)?;
    let source = typst_main_source(settings, atom);
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
    let aliases = IndexAliases::for_descriptor(&tensor.descriptor, &settings.index_style);
    tensor
        .descriptor
        .interface
        .logical_slots()
        .into_iter()
        .map(composition::port_atom)
        .map(|port| aliases.presentation_atom(&port, &settings.index_style))
        .map(|port| {
            port.format_string(
                &display_options_with_settings(mode, settings),
                PrintState::new(),
            )
        })
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
    let aliases = IndexAliases::for_descriptor(&tensor.descriptor, &settings.index_style);
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
        let port = aliases.presentation_atom(&port, &settings.index_style);
        let source = port.format_string(
            &display_options_with_settings(TensorDisplayMode::Typst, settings),
            PrintState::new(),
        );
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
        typst_settings_source(settings, &descriptor),
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

    fn graph_index(name: &str, arguments: &[usize]) -> Atom {
        use symbolica::atom::{NamespacedSymbol, SymbolBuilder};
        let name = NamespacedSymbol::parse(name);
        let head = Symbol::get_symbol(name.clone()).unwrap_or_else(|| {
            SymbolBuilder::new(name)
                .with_tags([SPENSO_TAG.index.clone()])
                .build()
                .unwrap()
        });
        let mut index = FunctionBuilder::new(head);
        for argument in arguments {
            index = index.add_arg(Atom::num(*argument as i64));
        }
        index.finish()
    }

    fn indexed_test_tensor(indices: impl IntoIterator<Item = Atom>) -> Atom {
        let mut tensor = FunctionBuilder::new(spenso::tensor_symbol!("display_index_tests::T"));
        let mink = LibraryRep::from(spenso::structure::representation::Minkowski {});
        for index in indices {
            tensor = tensor.add_arg(mink.to_symbolic([Atom::num(4), index]));
        }
        tensor.finish()
    }

    #[test]
    fn index_styles_are_validated_and_reported_in_settings() {
        for style in ["alphabet", "graph", "raw"] {
            let settings =
                DisplaySettings::new("ports", false, true, None, true, "0.08em", "0.12em", style)
                    .unwrap();
            assert_eq!(settings.index_style, style);
            assert!(settings.__repr__().contains(style));
        }
        assert!(
            DisplaySettings::new(
                "ports", false, true, None, true, "0.08em", "0.12em", "invalid"
            )
            .is_err()
        );
        assert_eq!(DisplaySettings::default().index_style, "alphabet");
    }

    #[test]
    fn compound_index_aliases_share_contractions_without_changing_exact_atoms() {
        let first = graph_index("display_index_tests::hedge", &[0, 1]);
        let second = graph_index("display_index_tests::hedge", &[1, 1]);
        let tensor = indexed_test_tensor([first.clone(), second.clone()]);
        let atom = tensor.clone()
            * indexed_test_tensor([first.clone()])
            * indexed_test_tensor([second.clone()]);
        let original = atom.to_canonical_string();
        let exact_tree = render_atom_tree(&atom).unwrap();
        let aliases = IndexAliases::for_atom(&atom, "alphabet");
        assert_eq!(aliases.entries.len(), 2);
        let presentation = aliases.presentation_atom(&atom, "alphabet");
        let expected = indexed_test_tensor([Atom::num(1), Atom::num(2)])
            * indexed_test_tensor([Atom::num(1)])
            * indexed_test_tensor([Atom::num(2)]);
        assert_eq!(presentation, expected);
        for mode in [
            TensorDisplayMode::Plain,
            TensorDisplayMode::Latex,
            TensorDisplayMode::Typst,
        ] {
            assert!(!format_atom_with_mode(&atom, mode, false).contains("hedge"));
        }
        let source = typst_main_source(&DisplaySettings::default(), &atom);
        assert!(source.contains("index-aliases:"));
        assert!(source.contains("display_index_tests::hedge"));
        assert_eq!(atom.to_canonical_string(), original);
        assert_eq!(render_atom_tree(&atom).unwrap(), exact_tree);
        assert_ne!(render_atom_tree(&presentation).unwrap(), exact_tree);
    }

    #[test]
    fn compound_index_alphabet_reserves_numeric_manual_and_wrapped_labels() {
        let manual = register_math_display_symbol(
            &IndexDisplay::symbol("rho").unwrap(),
            "display_index_tests",
        )
        .unwrap();
        let mut indices = vec![
            Atom::num(1),
            Atom::var(symbol!("nu")),
            Atom::var(manual),
            Atom::num(5),
        ];
        indices.extend((0..5).map(|i| graph_index("display_index_tests::hedge", &[i, 1])));
        let atom = indexed_test_tensor(indices);
        let aliases = IndexAliases::for_atom(&atom, "alphabet");
        let mut positions = aliases
            .entries
            .values()
            .map(|(position, _)| *position)
            .collect::<Vec<_>>();
        positions.sort();
        assert_eq!(positions, [4, 6, 7, 8, 9]);
        let raw = DisplaySettings {
            index_style: "raw".to_owned(),
            ..DisplaySettings::default()
        };
        assert!(format_atom_with_settings(&atom, TensorDisplayMode::Typst, &raw).contains("hedge"));
    }

    #[test]
    fn compound_index_aliases_preserve_heads_arity_and_whole_sum_identity() {
        let indices = [
            graph_index("display_index_tests::edge", &[0, 1]),
            graph_index("display_index_tests::hedge", &[0, 1]),
            graph_index("display_index_tests::vertex", &[0, 1]),
            graph_index("display_index_tests::hedge", &[0]),
            graph_index("display_index_tests::hedge", &[0, 0]),
        ];
        let atom = indices
            .into_iter()
            .map(|index| indexed_test_tensor([index]))
            .fold(Atom::Zero, |left, right| left + right);
        let aliases = IndexAliases::for_atom(&atom, "alphabet");
        assert_eq!(aliases.entries.len(), 5);
        let presentation = aliases.presentation_atom(&atom, "alphabet");
        let AtomView::Add(sum) = presentation.as_view() else {
            panic!("distinct terms combined while formatting")
        };
        assert_eq!(sum.get_nargs(), 5);
    }

    #[test]
    fn compound_index_graph_labels_are_unambiguous_and_use_each_representations_palette() {
        let atom = indexed_test_tensor([
            graph_index("gammalooprs::hedge", &[4, 1]),
            graph_index("gammalooprs::edge", &[4, 1]),
            graph_index("gammalooprs::vertex", &[4, 1]),
            graph_index("gammalooprs::hedge", &[4]),
            graph_index("gammalooprs::hedge", &[4, 0]),
            graph_index("gammalooprs::hedge", &[4, 2]),
        ]);
        let settings = DisplaySettings {
            index_style: "graph".to_owned(),
            ..DisplaySettings::default()
        };
        let aliases = IndexAliases::for_atom(&atom, "graph");
        let labels = aliases
            .entries
            .values()
            .map(|(_, display)| display.to_native_string())
            .collect::<HashSet<_>>();
        assert_eq!(labels.len(), 6);
        for expected in [
            "mu_(h4)",
            "mu_(e4)",
            "mu_(v4)",
            "mu_(h(4))",
            "mu_(h4.0)",
            "mu_(h4.2)",
        ] {
            assert!(labels.contains(expected), "missing {expected}: {labels:?}");
        }
        for mode in [
            TensorDisplayMode::Plain,
            TensorDisplayMode::Latex,
            TensorDisplayMode::Typst,
        ] {
            let rendered = format_atom_with_settings(&atom, mode, &settings);
            assert!(
                !rendered.contains("spenso_index_"),
                "display symbol escaped into output: {rendered}"
            );
            assert!(
                rendered.contains("h4") && !rendered.contains("h4.1"),
                "graph identity absent: {rendered}"
            );
        }
        assert!(typst_settings_source(&settings, &atom).contains("$attach(mu,b:upright(\"h4\"))$"));
    }

    #[test]
    fn compound_index_aliases_share_dual_palettes_and_leave_scalar_metadata_opaque() {
        let rep = LibraryRep::new_dual_with_index_palette(
            "display_index_tests::custom",
            IndexPalette::cyclic(1, [IndexDisplay::symbol("zeta").unwrap()]).unwrap(),
        )
        .unwrap();
        let index = graph_index("display_index_tests::hedge", &[0, 1]);
        let slot = rep.to_symbolic([Atom::num(4), index.clone()]);
        let dual = rep.dual().to_symbolic([Atom::num(4), index.clone()]);
        let hidden = indexed_test_tensor([graph_index("display_index_tests::hedge", &[9, 1])]);
        let scalar = function!(symbol!("display_index_tests::scalar"; Scalar), hidden);
        let tensor = function!(
            spenso::tensor_symbol!("display_index_tests::Dual"),
            slot,
            dual
        ) * scalar.clone();
        let aliases = IndexAliases::for_atom(&tensor, "alphabet");
        assert_eq!(aliases.entries.len(), 1);
        assert_eq!(
            aliases.entries.values().next().unwrap().1,
            IndexDisplay::symbol("zeta").unwrap()
        );
        let presentation = aliases.presentation_atom(&tensor, "alphabet");
        assert!(
            presentation
                .to_canonical_string()
                .contains(&scalar.to_canonical_string())
        );
        let formatted = format_atom_with_mode(&tensor, TensorDisplayMode::Typst, false);
        assert!(formatted.matches("zeta").count() >= 2);
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
                "alphabet",
            )
            .unwrap();
            let source = typst_settings_source(&settings, &Atom::Zero);
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
        let source = typst_main_source(&DisplaySettings::default(), &Atom::Zero);
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
