//! Local tensor notation for ordinary Symbolica expressions.
//!
//! The interactive callbacks replace bounded semantic leaves, leaving arithmetic
//! and paging to Symbolica. Static exports walk the complete expression through
//! the same notation owner and decline when their explicit limits are exceeded.
//!
//! Use [`render`] to offer tensor notation to an interactive printer, [`expression`]
//! to export an entire expression, and [`fragment`] to embed a printer's output.
//! Each entry point has an example demonstrating its distinct boundary.
use super::{
    COMPONENT_LABELS_SUBSCRIPT_TAG, ChainOrientation, CompactPolarity, RepresentationValue,
    SPENSO_TAG, TENSOR_PRINT_CALLBACK_TAG, TENSOR_PRINT_HEAD_PREFIX,
};
use crate::network::library::symbolic::ETS;
use crate::network::tags;
use crate::shadowing::symbolica_utils::SpensoPrintSettings;
use crate::shadowing::{CYCLIC, SYM};
use crate::structure::abstract_index::AIND_SYMBOLS;
use crate::structure::representation::{
    IndexDisplay, IndexRow, RepresentationClass, RepresentationMetadata,
};
use crate::structure::slot::{SlotMatch, SlotMatcher};
use std::cell::Cell;
use std::num::ParseFloatError;
use symbolica::atom::{Atom, AtomCore, AtomView, Symbol, UserData};
use symbolica::coefficient::CoefficientView;
use symbolica::printer::{AnsiHtmlFormatter, PrintMode, PrintOptions, PrintState, PrintUserData};
use symbolica::state::State;
use thiserror::Error;

const MAX_BYTES: usize = 4096;
const MAX_EXPORT_NODES: usize = 1_000_000;
const MAX_EXPORT_DEPTH: usize = 256;
const MAX_EXPORT_BYTES: usize = 16 * 1024 * 1024;

fn fits(atom: AtomView<'_>, remaining: &mut usize, depth: usize) -> bool {
    fn visit(atom: AtomView<'_>, remaining: &mut usize, bytes: &mut usize, depth: usize) -> bool {
        if *remaining == 0 || depth > 32 || atom.get_byte_size() > MAX_BYTES {
            return false;
        }
        *remaining -= 1;
        let size = atom.get_symbol().map_or_else(
            || {
                if matches!(atom, AtomView::Num(_)) {
                    atom.get_byte_size()
                } else {
                    2
                }
            },
            |symbol| symbol.get_name().len() + 2,
        );
        let Some(left) = bytes.checked_sub(size) else {
            return false;
        };
        *bytes = left;
        atom.children()
            .all(|child| visit(child, remaining, bytes, depth + 1))
    }
    let mut bytes = MAX_BYTES;
    visit(atom, remaining, &mut bytes, depth)
}

fn source(atom: AtomView<'_>, options: &PrintOptions) -> Option<String> {
    if !fits(atom, &mut 256, 0) {
        return None;
    }
    let atom = source_label(atom);
    let native = atom
        .printer(PrintOptions {
            mode: PrintMode::MathML,
            max_terms: Some(500),
            ..options.clone()
        })
        .to_string();
    fragment(native)
}

fn source_label(atom: AtomView<'_>) -> AtomView<'_> {
    // Reversible index symbols retain their original local label as Atom data.
    // Display simple payloads without renaming the index or inventing aliases.
    // Nested payloads remain opaque, avoiding recursive metadata expansion.
    if let AtomView::Var(variable) = atom
        && variable.get_symbol().has_tag(&SPENSO_TAG.index)
        && let UserData::Atom(payload) = variable.get_symbol().get_data()
        && let AtomView::Fun(label) = payload.as_view()
        && label.get_symbol().has_tag(&SPENSO_TAG.index)
        && fits(payload.as_view(), &mut 64, 0)
        && label.iter().all(|argument| match argument {
            AtomView::Num(_) => true,
            AtomView::Var(variable) => matches!(variable.get_symbol().get_data(), UserData::None),
            _ => false,
        })
    {
        payload.as_view()
    } else {
        atom
    }
}

/// Prepare a complete, path-free MathML fragment for embedding in a semantic leaf.
///
/// Native and custom printer output must fit within 4096 bytes and contain no
/// `math` root or omission placeholder. Remove actual `data-symbolica-*`
/// attributes while preserving token text and other attributes, including local
/// visual slash and open-axis annotations. Malformed opening attributes decline.
/// `None` means this fragment cannot be embedded under that contract; callers
/// offering optional notation may defer to their ordinary structural printer.
/// This removes viewer metadata, not arbitrary untrusted HTML.
///
/// # Examples
///
/// ```
/// use spenso::network::tags::mathml;
///
/// let native = r#"<mi data-symbolica-path="0">x</mi>"#.to_owned();
/// assert_eq!(mathml::fragment(native).as_deref(), Some("<mi>x</mi>"));
/// ```
pub fn fragment(value: String) -> Option<String> {
    bounded_fragment(value, MAX_BYTES)
}

fn bounded_fragment(mut value: String, max_bytes: usize) -> Option<String> {
    if value.len() > max_bytes || value.contains("<mtext>⋯</mtext>") || value.contains("<math") {
        return None;
    }
    let bytes = value.as_bytes();
    let mut cursor = 0;
    let mut attributes = Vec::new();
    while let Some(start) = value[cursor..].find('<').map(|offset| cursor + offset) {
        let tail = &value[start..];
        if tail.starts_with("<!--") {
            cursor = start + tail.find("-->")? + 3;
            continue;
        }
        if tail.starts_with("<![CDATA[") {
            cursor = start + tail.find("]]>")? + 3;
            continue;
        }
        if tail.starts_with("<?") {
            cursor = start + tail.find("?>")? + 2;
            continue;
        }
        cursor = start + 1;
        if matches!(bytes.get(cursor), Some(b'/' | b'!')) {
            cursor = start + tail.find('>')? + 1;
            continue;
        }
        // Token text is opaque. Only names at opening-tag attribute boundaries
        // can own viewer metadata; quoted values may themselves contain `>`.
        while let Some(byte) = bytes.get(cursor)
            && !byte.is_ascii_whitespace()
            && !matches!(byte, b'/' | b'>')
        {
            cursor += 1;
        }
        loop {
            let whitespace = cursor;
            while bytes.get(cursor).is_some_and(u8::is_ascii_whitespace) {
                cursor += 1;
            }
            if bytes.get(cursor) == Some(&b'>') {
                cursor += 1;
                break;
            }
            if bytes.get(cursor..cursor + 2) == Some(b"/>".as_slice()) {
                cursor += 2;
                break;
            }
            let name = cursor;
            while let Some(byte) = bytes.get(cursor)
                && !byte.is_ascii_whitespace()
                && !matches!(byte, b'=' | b'/' | b'>')
            {
                cursor += 1;
            }
            if cursor == name {
                return None;
            }
            let symbolic = value[name..cursor].starts_with("data-symbolica-");
            while bytes.get(cursor).is_some_and(u8::is_ascii_whitespace) {
                cursor += 1;
            }
            if bytes.get(cursor) != Some(&b'=') {
                return None;
            }
            cursor += 1;
            while bytes.get(cursor).is_some_and(u8::is_ascii_whitespace) {
                cursor += 1;
            }
            let quote = *bytes.get(cursor)?;
            if !matches!(quote, b'\'' | b'"') {
                return None;
            }
            cursor += 1;
            cursor += value[cursor..].find(char::from(quote))? + 1;
            if symbolic {
                attributes.push(whitespace..cursor);
            }
        }
    }
    for attribute in attributes.into_iter().rev() {
        value.replace_range(attribute, "");
    }
    value.trim_start().starts_with('<').then_some(value)
}

fn head(symbol: Symbol, options: Option<&PrintOptions>) -> Option<String> {
    if let Some(source) = symbol.get_tags().iter().find_map(|tag| {
        tag.strip_prefix(TENSOR_PRINT_HEAD_PREFIX)?
            .strip_prefix("mathml:")
    }) {
        return (source.len() + "<mrow></mrow>".len() <= MAX_BYTES)
            .then(|| format!("<mrow>{source}</mrow>"));
    }
    if let Some(display) = IndexDisplay::from_symbol(symbol) {
        return display.to_mathml();
    }
    let label = symbol
        .get_tags()
        .iter()
        .find_map(|tag| tag.strip_prefix("spenso::tensor-label:"))
        .unwrap_or_else(|| {
            if options.is_some_and(|options| {
                !options.hide_all_namespaces
                    && options.hide_namespace.as_deref() != Some(symbol.get_namespace())
                    && !State::BUILTIN_SYMBOLS.contains(&symbol)
            }) {
                symbol.get_name()
            } else {
                symbol.get_stripped_name()
            }
        });
    (label.len() <= 1024).then(|| format!("<mi>{}</mi>", AnsiHtmlFormatter::escape_html(label)))
}

fn group(body: String) -> String {
    format!("<mrow><mo>(</mo>{body}<mo>)</mo></mrow>")
}

fn arguments(values: &[AtomView<'_>], options: &PrintOptions) -> Option<String> {
    Some(
        values
            .iter()
            .map(|value| source(*value, options))
            .collect::<Option<Vec<_>>>()?
            .join("<mo>,</mo>"),
    )
}

fn labelled_head(
    symbol: Symbol,
    labels: &[AtomView<'_>],
    options: &PrintOptions,
) -> Option<String> {
    let head = head(symbol, None)?;
    if labels.is_empty() {
        return Some(head);
    }
    let labels = arguments(labels, options)?;
    Some(
        if symbol.has_tag(&SPENSO_TAG.rank1) || symbol.has_tag(COMPONENT_LABELS_SUBSCRIPT_TAG) {
            format!("<msub>{head}<mrow>{labels}</mrow></msub>")
        } else {
            format!("<mrow>{head}{}</mrow>", group(labels))
        },
    )
}

pub(crate) fn tensor(
    atom: AtomView<'_>,
    options: &PrintOptions,
    state: &PrintState,
) -> Option<String> {
    if !fits(atom, &mut 64, 0) {
        return None;
    }
    let symbol = atom.get_symbol()?;
    if let AtomView::Var(_) = atom {
        return symbol
            .has_tag(&SPENSO_TAG.tensor)
            .then(|| head(symbol, None))?;
    }
    let AtomView::Fun(function) = atom else {
        return None;
    };
    let values = function.iter().collect::<Vec<_>>();
    if symbol == ETS.metric {
        let [first, second] = values.as_slice() else {
            return None;
        };
        let first = tags::tensor_slot(*first)?;
        let second = tags::tensor_slot(*second)?;
        if first.representation != second.representation || first.dimension != second.dimension {
            return None;
        }
    }
    if symbol.has_tag(&SPENSO_TAG.index) {
        let label = symbol
            .get_tags()
            .iter()
            .find_map(|tag| tag.strip_prefix("spenso::index-label:"))?;
        if values.is_empty() || label.len() > 128 {
            return None;
        }
        let values = if matches!(label, "s" | "t")
            && values.len() == 2
            && tags::natural_index(values[1]) == Some(1)
        {
            &values[..1]
        } else {
            &values
        };
        let body = format!(
            "<msub><mi>{}</mi><mrow>{}</mrow></msub>",
            AnsiHtmlFormatter::escape_html(label),
            arguments(values, options)?
        );
        return (body.len() <= MAX_BYTES).then_some(body);
    }
    if !symbol.has_tag(&SPENSO_TAG.tensor)
        || values.iter().any(|value| {
            tags::is_chain_marker(*value, SPENSO_TAG.chain_in)
                || tags::is_chain_marker(*value, SPENSO_TAG.chain_out)
        })
    {
        return None;
    }
    let mut labels = Vec::new();
    let mut columns = Vec::new();
    let open_count = values
        .iter()
        .filter(|value| {
            tags::representation_value(**value).is_some()
                || tags::tensor_slot(**value).is_some_and(|slot| {
                    matches!(slot.index, AtomView::Fun(index)
                        if index.get_symbol() == AIND_SYMBOLS.openind && index.get_nargs() == 2)
                })
        })
        .count();
    let mut axis = 0;
    for value in &values {
        let slot = tags::tensor_slot(*value);
        let representation = tags::representation_value(*value);
        if slot.is_some() || representation.is_some() {
            let (index, row) = if let Some(slot) = slot {
                let open = matches!(slot.index, AtomView::Fun(index)
                    if index.get_symbol() == AIND_SYMBOLS.openind && index.get_nargs() == 2);
                ((!open).then_some(slot.index), slot.row)
            } else {
                (None, representation?.row)
            };
            columns.push((
                if let Some(index) = index {
                    source(index, options)?
                } else if open_count > 1 {
                    format!("<msub><mo>□</mo><mn>{axis}</mn></msub>")
                } else {
                    "<mo>□</mo>".to_owned()
                },
                row,
            ));
            axis += 1;
        } else if let AtomView::Fun(component) = value
            && component.get_symbol() == AIND_SYMBOLS.cind
        {
            if values.last().copied() != Some(*value)
                || !columns.is_empty()
                || component
                    .iter()
                    .any(|index| tags::natural_index(index).is_none())
            {
                return None;
            }
            columns.push((
                arguments(&component.iter().collect::<Vec<_>>(), options)?,
                IndexRow::Top,
            ));
        } else {
            labels.push(*value);
        }
    }
    let mut body = labelled_head(symbol, &labels, options)?;
    let upper = columns.iter().any(|(_, row)| *row == IndexRow::Top);
    if !columns.is_empty() {
        let mut scripts = String::new();
        for (index, row) in columns {
            let pair = match row {
                IndexRow::Top => format!("<none/><mrow>{index}</mrow>"),
                IndexRow::Bottom => format!("<mrow>{index}</mrow><none/>"),
            };
            scripts.push_str(&pair);
        }
        body = format!("<mmultiscripts>{body}{scripts}</mmultiscripts>");
    }
    if state.in_exp_base && upper {
        body = group(body);
    }
    (body.len() <= MAX_BYTES).then_some(body)
}

/// Standalone spaces and slots retain their family and dimension. Inside a
/// tensor those qualifiers belong to the tensor's type, so only its index is
/// shown by `tensor` instead.
pub(crate) fn representation(
    atom: AtomView<'_>,
    options: &PrintOptions,
    state: &PrintState,
) -> Option<String> {
    if !fits(atom, &mut 64, 0) {
        return None;
    }
    let (symbol, dimension, index, class, row) = if let Some(slot) = tags::tensor_slot(atom) {
        (
            slot.representation,
            slot.dimension,
            Some(slot.index),
            RepresentationMetadata::from_symbol(slot.representation)?.class,
            slot.row,
        )
    } else {
        let value = tags::representation_value(atom)?;
        (value.symbol, value.dimension, None, value.class, value.row)
    };
    // A bare symbol defers to Symbolica, retaining its namespace/style options
    // without calling this representation printer recursively on the function.
    let mut body = format!(
        "<msub>{}<mrow>{}</mrow></msub>",
        source(Atom::var(symbol).as_view(), options)?,
        source(dimension, options)?,
    );
    let dual = matches!(tags::unwrap_tensor_display(atom), AtomView::Fun(wrapper)
        if wrapper.get_symbol() == AIND_SYMBOLS.dind);
    if dual && class == RepresentationClass::Dualizable {
        body = format!("<msup>{body}<mo>*</mo></msup>");
    }
    if let Some(index) = index {
        let index = source(index, options)?;
        let pair = match row {
            IndexRow::Top => format!("<none/><mrow>{index}</mrow>"),
            IndexRow::Bottom => format!("<mrow>{index}</mrow><none/>"),
        };
        body = format!("<mmultiscripts>{body}{pair}</mmultiscripts>");
    }
    if state.in_exp_base
        && ((dual && class == RepresentationClass::Dualizable)
            || (index.is_some() && row == IndexRow::Top))
    {
        body = group(body);
    }
    (body.len() <= MAX_BYTES).then_some(body)
}

fn vector<'a>(
    atom: AtomView<'a>,
    options: &PrintOptions,
) -> Option<(String, RepresentationValue<'a>)> {
    let AtomView::Fun(function) = tags::unwrap_tensor_display(atom) else {
        return None;
    };
    let symbol = function.get_symbol();
    if !symbol.has_tag(&SPENSO_TAG.tensor)
        || !symbol.has_tag(&SPENSO_TAG.rank1)
        || symbol.has_tag(TENSOR_PRINT_CALLBACK_TAG)
    {
        return None;
    }
    let mut representation = None;
    let mut labels = Vec::new();
    for value in function.iter() {
        if tags::tensor_slot(value).is_some() {
            return None;
        }
        if let Some(value) = tags::representation_value(value) {
            if representation.replace(value).is_some() {
                return None;
            }
        } else {
            labels.push(value);
        }
    }
    Some((labelled_head(symbol, &labels, options)?, representation?))
}

pub(super) fn dot(
    atom: AtomView<'_>,
    options: &PrintOptions,
    state: &PrintState,
) -> Option<String> {
    if !fits(atom, &mut 64, 0) {
        return None;
    }
    let AtomView::Fun(function) = atom else {
        return None;
    };
    if function.get_nargs() != 2 {
        return None;
    }
    let mut values = function.iter();
    let (a, a_rep) = vector(values.next()?, options)?;
    let (b, b_rep) = vector(values.next()?, options)?;
    if a_rep.symbol != b_rep.symbol
        || a_rep.dimension != b_rep.dimension
        || !matches!(
            a_rep.class,
            RepresentationClass::SelfDual | RepresentationClass::InlineMetric
        )
    {
        return None;
    }
    let mut body = format!("<mrow>{a}<mo>·</mo>{b}</mrow>");
    if state.in_product || state.in_exp_base {
        body = group(body);
    }
    (body.len() <= MAX_BYTES).then_some(body)
}

/// Render one complete, bounded semantic tensor leaf without a `<math>` root.
///
/// The index resolver receives the representation and the original index Atom;
/// returning `None` keeps its representation-owned palette or ordinary label.
/// Arithmetic roots remain owned by Symbolica's structural expression viewer.
/// `None` declines this notation offer for an unsupported leaf, invalid settings,
/// or a leaf outside the interactive bounds. Use [`expression`] when a complete
/// export or a diagnostic is required.
///
/// # Examples
///
/// ```
/// use spenso::network::tags::mathml;
/// use spenso::tensor_symbol;
/// use symbolica::atom::{Atom, AtomCore};
/// use symbolica::printer::{PrintOptions, PrintState};
///
/// let tensor = Atom::var(tensor_symbol!("mathml_example::T"));
/// let html = mathml::render(
///     tensor.as_view(), &PrintOptions::mathml(), &PrintState::new(), &|_, _| None,
/// );
/// assert_eq!(html.as_deref(), Some("<mi>T</mi>"));
/// ```
pub fn render(
    atom: AtomView<'_>,
    options: &PrintOptions,
    state: &PrintState,
    index: &dyn for<'a> Fn(Symbol, AtomView<'a>) -> Option<String>,
) -> Option<String> {
    if !fits(atom, &mut 256, 0) {
        return None;
    }
    let context = Render::new(options, index).ok()?;
    let atom = tags::unwrap_tensor_display(atom);
    fragment(context.leaf(atom, state).or_else(|| {
        let symbol = atom.get_symbol()?;
        if context.owns(atom) {
            return None;
        }
        (symbol.get_print_function()?)(atom, options, state)
    })?)
}

/// Render a complete static MathML fragment using the interactive tensor notation.
///
/// No arithmetic or function argument is paged or replaced by an omission marker.
/// The source is limited to 1,000,000 nodes and depth 256, and generated fragments
/// share a 16 MiB rendering budget.
/// Native numeric leaves retain their existing 4 KiB/precision bounds.
/// The index resolver has the same contract as [`render`].
///
/// # Errors
///
/// Returns [`ExportError::InvalidSettings`] for an unsupported layout or invalid
/// spacing, [`ExportError::SourceLimit`] for source bounds, and
/// [`ExportError::IncompleteOutput`] when a complete supported fragment cannot be
/// produced within the rendering bounds. No partial expression is returned.
///
/// # Examples
///
/// ```
/// use spenso::network::tags::mathml;
/// use spenso::network::tags::mathml::ExportError;
/// use symbolica::atom::{Atom, AtomCore};
/// use symbolica::printer::{PrintOptions, PrintState};
///
/// let value = Atom::num(7);
/// let html = mathml::expression(
///     value.as_view(), &PrintOptions::mathml(), &PrintState::new(), &|_, _| None,
/// )?;
/// assert!(html.contains("<mn>7</mn>"));
/// # Ok::<(), ExportError>(())
/// ```
pub fn expression(
    atom: AtomView<'_>,
    options: &PrintOptions,
    state: &PrintState,
    index: &dyn for<'a> Fn(Symbol, AtomView<'a>) -> Option<String>,
) -> Result<String, ExportError> {
    let mut remaining = MAX_EXPORT_NODES;
    if !Render::expression_fits(atom, &mut remaining, 0) {
        return Err(ExportError::SourceLimit);
    }
    let options = PrintOptions {
        mode: PrintMode::MathML,
        ..options.clone()
    };
    let mut context = Render::new(&options, index)?;
    context.export = Some(Export {
        remaining: Cell::new(MAX_EXPORT_BYTES),
        state: Cell::new(*state),
        failed: Cell::new(false),
    });
    let body = context
        .visual(atom, state)
        .ok_or(ExportError::IncompleteOutput)?;
    bounded_fragment(body, MAX_EXPORT_BYTES).ok_or(ExportError::IncompleteOutput)
}

/// Why a complete MathML export could not be produced.
///
/// See [`expression`] for successful export and error propagation. Semantic
/// callback declines are instead represented by [`render`]'s optional result.
#[derive(Clone, Debug, Error, PartialEq, Eq)]
pub enum ExportError {
    /// An unsupported layout or malformed/non-finite spacing setting.
    #[error("invalid MathML rendering setting `{key}`")]
    InvalidSettings {
        /// The key in [`PrintOptions::custom_print_mode`] that could not be used.
        key: &'static str,
        /// The numeric parsing failure, when spacing failed to parse.
        #[source]
        source: Option<ParseFloatError>,
    },
    /// The input exceeds the source node, depth, byte, or symbol-name limits.
    #[error("MathML export source exceeds the supported size or depth limits")]
    SourceLimit,
    /// A complete fragment could not be rendered within its output constraints.
    ///
    /// These include cumulative allocation and native numeric bounds, malformed
    /// fragments, and unsupported index display metadata. See [`expression`].
    #[error("MathML export could not produce a complete fragment within the rendering limits")]
    IncompleteOutput,
}

struct Export {
    remaining: Cell<usize>,
    state: Cell<PrintState>,
    failed: Cell<bool>,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Layout {
    Ports,
    Schoonschip,
    Call,
}

struct Render<'a> {
    options: &'a PrintOptions,
    index: &'a dyn for<'b> Fn(Symbol, AtomView<'b>) -> Option<String>,
    settings: SpensoPrintSettings,
    layout: Layout,
    index_gap: String,
    factor_gap: String,
    export: Option<Export>,
}

impl<'a> Render<'a> {
    fn expression_fits(atom: AtomView<'_>, remaining: &mut usize, depth: usize) -> bool {
        if *remaining == 0
            || depth > MAX_EXPORT_DEPTH
            || atom.get_byte_size() > MAX_EXPORT_BYTES
            || atom
                .get_symbol()
                .is_some_and(|symbol| symbol.get_name().len() > MAX_BYTES)
        {
            return false;
        }
        *remaining -= 1;
        atom.children()
            .all(|child| Self::expression_fits(child, remaining, depth + 1))
    }

    fn fragment(&self, value: String) -> Option<String> {
        let Some(export) = &self.export else {
            return fragment(value);
        };
        let Some(remaining) = export.remaining.get().checked_sub(value.len()) else {
            export.failed.set(true);
            return None;
        };
        if export.failed.get() {
            return None;
        }
        export.remaining.set(remaining);
        bounded_fragment(value, MAX_EXPORT_BYTES)
    }

    fn child_state(&self) -> PrintState {
        self.export.as_ref().map_or_else(PrintState::new, |export| {
            let state = export.state.get();
            PrintState {
                in_sum: false,
                in_product: false,
                in_exp: false,
                in_exp_base: false,
                level: state.level.saturating_add(1),
                ..state
            }
        })
    }

    fn owns(&self, atom: AtomView<'_>) -> bool {
        atom.get_symbol().is_some_and(|symbol| {
            symbol.has_tag(&SPENSO_TAG.tensor)
                || symbol.has_tag(&SPENSO_TAG.index)
                || symbol.has_tag(&SPENSO_TAG.representation)
                || [
                    SPENSO_TAG.dot,
                    SPENSO_TAG.chain,
                    SPENSO_TAG.trace,
                    SPENSO_TAG.bracket,
                    AIND_SYMBOLS.dind,
                ]
                .contains(&symbol)
        })
    }

    fn new(
        options: &'a PrintOptions,
        index: &'a dyn for<'b> Fn(Symbol, AtomView<'b>) -> Option<String>,
    ) -> Result<Self, ExportError> {
        let text = |key: &str| match options.custom_print_mode.get(key) {
            Some(PrintUserData::String(value)) => Some(value.as_str()),
            _ => None,
        };
        let layout = match text("spenso.tensor_layout").unwrap_or("ports") {
            "ports" => Layout::Ports,
            "schoonschip" => Layout::Schoonschip,
            "call" => Layout::Call,
            _ => {
                return Err(ExportError::InvalidSettings {
                    key: "spenso.tensor_layout",
                    source: None,
                });
            }
        };
        let settings = SpensoPrintSettings::resolve(options)
            .map(|resolved| resolved.presentation)
            .unwrap_or_else(|| {
                let mut settings = SpensoPrintSettings::typst();
                settings.commas = layout == Layout::Call;
                settings
            });
        let gap = |key, default| {
            let value = text(key).unwrap_or(default);
            let (number, unit) = ["pt", "mm", "cm", "in", "em", "%", "ex", "px"]
                .into_iter()
                .find_map(|unit| value.strip_suffix(unit).map(|number| (number, unit)))
                .ok_or(ExportError::InvalidSettings { key, source: None })?;
            let number = number
                .parse::<f64>()
                .map_err(|source| ExportError::InvalidSettings {
                    key,
                    source: Some(source),
                })?;
            if !number.is_finite() {
                return Err(ExportError::InvalidSettings { key, source: None });
            }
            Ok(format!("{number}{unit}"))
        };
        Ok(Self {
            options,
            index,
            settings,
            layout,
            index_gap: gap("spenso.index_gap", "0.08em")?,
            factor_gap: gap("spenso.factor_gap", "0.12em")?,
            export: None,
        })
    }

    fn space(gap: &str) -> String {
        format!("<mspace width=\"{gap}\"/>")
    }

    fn fence(body: String, left: &str, right: &str) -> String {
        format!("<mrow><mo>{left}</mo>{body}<mo>{right}</mo></mrow>")
    }

    fn head(&self, symbol: Symbol) -> Option<String> {
        if symbol
            .get_tags()
            .iter()
            .any(|tag| tag.starts_with("spenso::print-head:mathml:"))
            || IndexDisplay::from_symbol(symbol).is_some()
        {
            return self.fragment(head(symbol, None)?);
        }
        let body = match symbol.get_name() {
            "spenso::gamma" => "<mi>γ</mi>".to_owned(),
            "spenso::gamma0" => "<msub><mi>γ</mi><mn>0</mn></msub>".to_owned(),
            "spenso::gamma5" => "<msub><mi>γ</mi><mn>5</mn></msub>".to_owned(),
            "spenso::projp" => "<msup><mi>ℙ</mi><mo>+</mo></msup>".to_owned(),
            "spenso::projm" => "<msup><mi>ℙ</mi><mo>−</mo></msup>".to_owned(),
            _ => head(symbol, None)?,
        };
        if self.export.is_some() {
            self.fragment(body)
        } else {
            Some(body)
        }
    }

    fn leaf(&self, atom: AtomView<'_>, state: &PrintState) -> Option<String> {
        let atom = tags::unwrap_tensor_display(atom);
        let symbol = atom.get_symbol()?;
        // Static scalar variables go straight to Symbol::format(), which owns
        // both their callback and native name. A declined callback runs once.
        if symbol.has_tag(TENSOR_PRINT_CALLBACK_TAG)
            && !(self.export.is_some() && matches!(atom, AtomView::Var(_)) && !self.owns(atom))
            && let Some(custom) = tags::tensor_custom_print(atom, self.options, state)
        {
            if let Some(custom) = self.fragment(custom) {
                return Some(custom);
            }
            if !symbol.has_tag(&SPENSO_TAG.tensor) {
                return None;
            }
            // An unusable override does not hide a bounded tensor's arguments.
            // Continue through its existing semantic notation without recalling
            // the callable through the global printer.
        }
        if symbol == SPENSO_TAG.dot {
            return self.dot(atom, state);
        }
        if symbol == SPENSO_TAG.chain {
            return self.chain(atom, state);
        }
        if symbol == SPENSO_TAG.trace {
            return self.trace(atom, state);
        }
        if symbol == SPENSO_TAG.bracket {
            let AtomView::Fun(function) = atom else {
                return None;
            };
            if function.get_nargs() == 0 {
                return None;
            }
            let body = self.factors(&function.iter().collect::<Vec<_>>())?;
            return Some(if state.in_exp_base || state.in_exp {
                group(body)
            } else {
                body
            });
        }
        if tags::tensor_slot(atom).is_some() || tags::representation_value(atom).is_some() {
            return self.representation(atom, state);
        }
        if !symbol.has_tag(&SPENSO_TAG.tensor) && !symbol.has_tag(&SPENSO_TAG.index) {
            return None;
        }
        if let AtomView::Var(_) = atom {
            if symbol.has_tag(&SPENSO_TAG.index) {
                return if let Some(display) = IndexDisplay::from_symbol(symbol) {
                    self.fragment(display.to_mathml()?)
                } else if self.export.is_some() {
                    let label = source_label(atom);
                    if label != atom {
                        self.visual(label, state)
                    } else {
                        let mut body = String::new();
                        symbol.format(self.options, *state, &mut body).ok()?;
                        self.fragment(body)
                    }
                } else {
                    source(atom, self.options)
                };
            }
            let body = self.head(symbol)?;
            return Some(if state.in_exp_base && body.contains("<msup") {
                group(body)
            } else {
                body
            });
        }
        let AtomView::Fun(function) = atom else {
            return None;
        };
        let arguments = function.iter().collect::<Vec<_>>();
        if symbol.has_tag(&SPENSO_TAG.index) {
            let label = symbol
                .get_tags()
                .iter()
                .find_map(|tag| tag.strip_prefix("spenso::index-label:"))?;
            if arguments.is_empty() || label.len() > 128 {
                return None;
            }
            let arguments = if matches!(label, "s" | "t")
                && arguments.len() == 2
                && tags::natural_index(arguments[1]) == Some(1)
            {
                &arguments[..1]
            } else {
                &arguments
            };
            return Some(format!(
                "<msub><mi>{}</mi><mrow>{}</mrow></msub>",
                AnsiHtmlFormatter::escape_html(label),
                self.arguments(arguments, "<mo>.</mo>")?
            ));
        }
        if symbol == ETS.metric {
            let [first, second] = arguments.as_slice() else {
                return None;
            };
            let mut matcher = SlotMatcher::default();
            let mut representation = |argument| {
                let argument = tags::unwrap_tensor_display(argument);
                if tags::tensor_slot(argument).is_none()
                    && tags::representation_value(argument).is_none()
                {
                    return None;
                }
                match matcher.classify(argument) {
                    SlotMatch::Explicit(slot) => Some(slot.representation()),
                    _ => matcher.compact_representation(argument),
                }
            };
            if !representation(*first)?.matches(representation(*second)?) {
                return None;
            }
        }
        if let Some(AtomView::Fun(component)) = arguments.last()
            && component.get_symbol() == AIND_SYMBOLS.cind
        {
            let labels = &arguments[..arguments.len() - 1];
            let coordinates = component.iter().collect::<Vec<_>>();
            if coordinates
                .iter()
                .any(|index| tags::natural_index(*index).is_none())
            {
                return None;
            }
            let mut base = self.head(symbol)?;
            if !labels.is_empty() {
                let labels = self.arguments(labels, "<mo>,</mo>")?;
                base = if symbol.has_tag(COMPONENT_LABELS_SUBSCRIPT_TAG) {
                    format!("<msub>{base}<mrow>{labels}</mrow></msub>")
                } else {
                    format!("<mrow>{base}{}</mrow>", group(labels))
                };
            }
            if coordinates.is_empty() {
                return Some(base);
            }
            let coordinates = self.arguments(&coordinates, "<mo>,</mo>")?;
            if self.settings.array_components {
                return Some(format!(
                    "<mrow>{base}{}</mrow>",
                    Self::fence(coordinates, "[", "]")
                ));
            }
            let body = format!("<msup>{base}<mrow>{coordinates}</mrow></msup>");
            return Some(if state.in_exp_base { group(body) } else { body });
        }
        let mut base = self.head(symbol)?;
        let selected = if symbol.get_name() == "spenso::gamma" && arguments.len() == 3 {
            if let Some((label, _)) = self.vector(arguments[2]) {
                // The Lorentz port is consumed by the momentum, never by an AUTO square.
                base = format!(
                    "<mrow data-spenso-cancel=\"updiagonal\" style=\"background-image:linear-gradient(to bottom right,transparent calc(50% - 0.025em),currentColor calc(50% - 0.025em),currentColor calc(50% + 0.025em),transparent calc(50% + 0.025em))\">{label}</mrow>"
                );
                &arguments[..2]
            } else {
                &arguments
            }
        } else {
            &arguments
        };
        self.tensor_body(base, selected, state)
    }

    fn visual(&self, atom: AtomView<'_>, state: &PrintState) -> Option<String> {
        let Some(export) = &self.export else {
            return self.visual_body(atom, state);
        };
        let previous = export.state.replace(*state);
        let body = self.visual_body(atom, state);
        export.state.set(previous);
        if export.failed.get() {
            return None;
        }
        body
    }

    fn visual_body(&self, atom: AtomView<'_>, state: &PrintState) -> Option<String> {
        let atom = tags::unwrap_tensor_display(atom);
        let body = match atom {
            AtomView::Add(sum) => {
                let mut term_state = self.child_state();
                term_state.in_sum = self.export.is_some();
                let body =
                    self.sequence(&sum.iter().collect::<Vec<_>>(), "<mo>+</mo>", &term_state)?;
                if state.in_product || state.in_exp_base {
                    group(body)
                } else {
                    format!("<mrow>{body}</mrow>")
                }
            }
            AtomView::Mul(product) => {
                let body = self.factors(&product.iter().collect::<Vec<_>>())?;
                if state.in_exp_base { group(body) } else { body }
            }
            AtomView::Pow(power) => {
                let (base, exponent) = power.get_base_exp();
                // Match Symbolica's native numeric power-base rule: a bare
                // nonnegative integer is unambiguous, other numbers need fences.
                let numeric_fence = matches!(base, AtomView::Num(number)
                    if !matches!(number.get_coeff_view(), CoefficientView::Natural(r, 1, 0, 1) if r >= 0));
                let mut base_state = self.child_state();
                base_state.in_exp_base = true;
                let mut base = self.visual(base, &base_state)?;
                if numeric_fence {
                    base = group(base);
                }
                let mut exponent_state = self.child_state();
                exponent_state.in_exp = true;
                let exponent = self.visual(exponent, &exponent_state)?;
                let body = format!("<msup><mrow>{base}</mrow><mrow>{exponent}</mrow></msup>");
                if state.in_exp_base { group(body) } else { body }
            }
            _ => {
                if let Some(body) = self.leaf(atom, state) {
                    body
                } else if self.export.is_some() {
                    self.structural(atom, state)?
                } else if self.owns(atom) {
                    return None;
                } else if let Some(custom) = atom
                    .get_symbol()
                    .and_then(|symbol| (symbol.get_print_function()?)(atom, self.options, state))
                {
                    fragment(custom)?
                } else {
                    source(atom, self.options)?
                }
            }
        };
        self.fragment(body)
    }

    fn structural(&self, atom: AtomView<'_>, state: &PrintState) -> Option<String> {
        match atom {
            AtomView::Num(number) => {
                if !fits(atom, &mut 256, 0)
                    || self
                        .options
                        .precision
                        .is_some_and(|precision| precision > MAX_BYTES)
                {
                    return None;
                }
                let polynomial = matches!(
                    number.get_coeff_view(),
                    CoefficientView::RationalPolynomial(_)
                );
                let body = if polynomial {
                    // A Num root has no pageable children. Disable the separate
                    // plain coefficient term limit too, so exotic numeric
                    // domains either remain exact or fail the atomic byte cap.
                    fragment(
                        atom.printer(PrintOptions {
                            max_terms: None,
                            ..self.options.clone()
                        })
                        .to_string(),
                    )?
                } else {
                    source(atom, self.options)?
                };
                Some(
                    if state.in_product
                        && (polynomial
                            || !number.get_coeff_view().is_real()
                            || number.get_coeff_view().to_owned().is_negative())
                    {
                        group(body)
                    } else {
                        body
                    },
                )
            }
            AtomView::Var(variable) if !self.owns(atom) => {
                let mut body = String::new();
                variable
                    .get_symbol()
                    .format(self.options, *state, &mut body)
                    .ok()?;
                self.fragment(body)
            }
            AtomView::Var(variable) => self.head(variable.get_symbol()),
            AtomView::Fun(function) => {
                let symbol = function.get_symbol();
                // Semantic callbacks were already tried by leaf(). Ordinary
                // callables receive the original Atom and surrounding state.
                if !self.owns(atom)
                    && !symbol.has_tag(TENSOR_PRINT_CALLBACK_TAG)
                    && let Some(custom) = symbol.get_print_function()
                    && let Some(custom) = custom(atom, self.options, state)
                    && let Some(custom) = self.fragment(custom)
                {
                    return Some(custom);
                }
                let arguments = function.iter().collect::<Vec<_>>();
                let count = arguments.len();
                if symbol == Symbol::DERIVATIVE
                    && count >= 3
                    && count % 2 == 1
                    && arguments[count / 2].get_symbol().is_some()
                {
                    // Native der(d1,...,dn,f,x1,...,xn) notation, without its
                    // paging limits. Recurse through the original head too:
                    // a function-valued head must retain all its children.
                    let count = count / 2;
                    let mut head_state = self.child_state();
                    head_state.in_exp_base = true;
                    let head = self.visual(arguments[count], &head_state)?;
                    let orders = self.arguments(&arguments[..count], "<mo>,</mo>")?;
                    let arguments = self.arguments(&arguments[count + 1..], "<mo>,</mo>")?;
                    let body = format!(
                        "<mrow><msup><mrow>{head}</mrow>{}</msup><mo>&#x2061;</mo>{}</mrow>",
                        group(orders),
                        group(arguments)
                    );
                    return Some(if state.in_exp_base { group(body) } else { body });
                }
                let head = if !self.owns(atom) && symbol.get_print_function().is_none() {
                    let mut head = String::new();
                    symbol.format(self.options, *state, &mut head).ok()?;
                    self.fragment(head)?
                } else if !self.owns(atom) {
                    // The function callback already declined. Keep its source
                    // namespace through the bounded semantic head owner rather
                    // than retrying the callback on a bare variable.
                    self.fragment(head(symbol, Some(self.options))?)?
                } else {
                    self.head(symbol)?
                };
                let arguments = self.arguments(&arguments, "<mo>,</mo>")?;
                Some(format!("<mrow>{head}{}</mrow>", group(arguments)))
            }
            _ => None,
        }
    }

    fn arguments(&self, arguments: &[AtomView<'_>], separator: &str) -> Option<String> {
        self.sequence(arguments, separator, &self.child_state())
    }

    fn sequence(
        &self,
        arguments: &[AtomView<'_>],
        separator: &str,
        state: &PrintState,
    ) -> Option<String> {
        let mut output = String::new();
        for (position, argument) in arguments.iter().enumerate() {
            if position != 0 {
                output.push_str(separator);
            }
            output.push_str(&self.visual(*argument, state)?);
            if output.len() > self.export.as_ref().map_or(MAX_BYTES, |_| MAX_EXPORT_BYTES) {
                return None;
            }
        }
        Some(output)
    }

    fn factors(&self, factors: &[AtomView<'_>]) -> Option<String> {
        let mut state = self.child_state();
        state.in_product = true;
        let output = self.sequence(factors, &Self::space(&self.factor_gap), &state)?;
        Some(format!("<mrow>{output}</mrow>"))
    }

    fn index(&self, symbol: Symbol, index: AtomView<'_>) -> Option<String> {
        if let Some(alias) = (self.index)(symbol, index) {
            return self.fragment(alias);
        }
        if let Some(index) = tags::natural_index(index)
            && let Some(display) = RepresentationMetadata::from_symbol(symbol)?
                .index_palette
                .resolve(index)
        {
            return self.fragment(display.to_mathml()?);
        }
        if let AtomView::Var(variable) = index
            && variable.get_symbol().has_tag(&SPENSO_TAG.index)
            && let Some(display) = IndexDisplay::from_symbol(variable.get_symbol())
        {
            return self.fragment(display.to_mathml()?);
        }
        self.visual(index, &self.child_state())
    }

    fn qualified(&self, source: String, symbol: Symbol, dimension: AtomView<'_>) -> Option<String> {
        if !self.settings.with_dim {
            return Some(source);
        }
        let label = RepresentationMetadata::from_symbol(symbol)?
            .label
            .to_mathml()?;
        let dimension = self.visual(dimension, &self.child_state())?;
        Some(format!(
            "<msup><mrow>{source}</mrow><msub><mrow>{label}</mrow><mrow>{dimension}</mrow></msub></msup>"
        ))
    }

    fn port(&self, rep: RepresentationValue<'_>) -> Option<String> {
        let shape = match (rep.class, rep.polarity) {
            (RepresentationClass::InlineMetric, _) => "■",
            (RepresentationClass::SelfDual, _) => "●",
            (_, CompactPolarity::Bra) => "▶",
            (_, CompactPolarity::Ket) => "◀",
        };
        self.qualified(format!("<mo>{shape}</mo>"), rep.symbol, rep.dimension)
    }

    fn open_columns(&self, arguments: &[AtomView<'_>]) -> Option<Vec<Option<(String, IndexRow)>>> {
        let open = |argument| {
            tags::representation_value(argument).is_some() || tags::tensor_slot(argument).is_some_and(|slot| {
            matches!(slot.index, AtomView::Fun(index) if index.get_symbol() == AIND_SYMBOLS.openind && index.get_nargs() == 2)
        })
        };
        let numbered = arguments.iter().filter(|argument| open(**argument)).count() > 1;
        let mut axis = 0;
        arguments
            .iter()
            .map(|argument| {
                let slot = tags::tensor_slot(*argument);
                let rep = tags::representation_value(*argument);
                let position = axis;
                if slot.is_some() || rep.is_some() {
                    axis += 1;
                }
                if !open(*argument) {
                    return Some(None);
                }
                let (symbol, dimension, row) = if let Some(slot) = slot {
                    (slot.representation, slot.dimension, slot.row)
                } else {
                    let rep = rep?;
                    (rep.symbol, rep.dimension, rep.row)
                };
                let source = if numbered {
                    format!("<msub><mo>□</mo><mn>{position}</mn></msub>")
                } else {
                    "<mo>□</mo>".to_owned()
                };
                let source = self.qualified(source, symbol, dimension)?;
                let metadata = RepresentationMetadata::from_symbol(symbol)?;
                let dual = metadata.class == RepresentationClass::Dualizable
                    && row != metadata.index_row;
                let title = format!(
                    "Unresolved axis {position} · {}({dimension}){}",
                    symbol.get_name(),
                    if dual { " · dual" } else { "" },
                );
                let title = AnsiHtmlFormatter::escape_html(&title);
                Some(Some((
                    format!("<mrow title=\"{title}\" data-spenso-open-axis=\"{position}\">{source}</mrow>"),
                    row,
                )))
            })
            .collect()
    }

    fn attachment(&self, base: String, columns: &[(String, IndexRow)]) -> String {
        if columns.is_empty() {
            return base;
        }
        let mut scripts = String::new();
        for (position, (source, row)) in columns.iter().enumerate() {
            if position != 0 {
                let gap = Self::space(&self.index_gap);
                scripts.push_str(&format!("<mrow>{gap}</mrow><mrow>{gap}</mrow>"));
            }
            scripts.push_str(&match row {
                IndexRow::Top => format!("<none/><mrow>{source}</mrow>"),
                IndexRow::Bottom => format!("<mrow>{source}</mrow><none/>"),
            });
        }
        format!("<mmultiscripts><mrow>{base}</mrow>{scripts}</mmultiscripts>")
    }

    fn call(&self, base: String, columns: &[(String, IndexRow)]) -> String {
        if columns.is_empty() {
            return base;
        }
        let gap = Self::space(&self.index_gap);
        let comma_separator = self.settings.commas.then(|| format!("<mo>,</mo>{gap}"));
        let separator = comma_separator.as_deref().unwrap_or(&gap);
        let row = |row| {
            columns
                .iter()
                .filter(|(_, side)| *side == row)
                .map(|(source, _)| source.as_str())
                .collect::<Vec<_>>()
                .join(separator)
        };
        let (top, bottom) = (row(IndexRow::Top), row(IndexRow::Bottom));
        let arguments = if bottom.is_empty() {
            top
        } else if top.is_empty() {
            format!("<mo>;</mo>{gap}{bottom}")
        } else {
            format!("{top}{gap}<mo>;</mo>{gap}{bottom}")
        };
        format!("<mrow>{base}{}</mrow>", group(arguments))
    }

    fn vector<'b>(&self, atom: AtomView<'b>) -> Option<(String, RepresentationValue<'b>)> {
        let AtomView::Fun(function) = tags::unwrap_tensor_display(atom) else {
            return None;
        };
        let symbol = function.get_symbol();
        if !symbol.has_tag(&SPENSO_TAG.tensor)
            || !symbol.has_tag(&SPENSO_TAG.rank1)
            || symbol.has_tag(TENSOR_PRINT_CALLBACK_TAG)
        {
            return None;
        }
        let mut rep = None;
        let mut labels = Vec::new();
        let mut matcher = SlotMatcher::default();
        for argument in function.iter() {
            if tags::tensor_slot(argument).is_some()
                || matches!(matcher.classify(argument), SlotMatch::Explicit(_))
            {
                return None;
            }
            if let Some(candidate) = tags::representation_value(argument) {
                if rep.replace(candidate).is_some() {
                    return None;
                }
            } else {
                labels.push(argument);
            }
        }
        let mut label = self.head(symbol)?;
        if !labels.is_empty() {
            label = if self.settings.symbol_scripts {
                let columns = labels
                    .iter()
                    .map(|label| {
                        Some((self.visual(*label, &self.child_state())?, IndexRow::Bottom))
                    })
                    .collect::<Option<Vec<_>>>()?;
                self.attachment(label, &columns)
            } else {
                let separator = if self.settings.commas {
                    "<mo>,</mo>".to_owned()
                } else {
                    Self::space("0.15em")
                };
                format!(
                    "<mrow>{label}{}</mrow>",
                    group(self.arguments(&labels, &separator)?)
                )
            };
        }
        Some((label, rep?))
    }

    fn tensor_body(
        &self,
        mut base: String,
        arguments: &[AtomView<'_>],
        state: &PrintState,
    ) -> Option<String> {
        let head_upper = base.contains("<msup");
        let mut columns = Vec::new();
        let mut ordinary = Vec::new();
        let mut bras = Vec::new();
        let mut kets = Vec::new();
        let markers = tags::chain_markers(arguments);
        // A malformed marker pair must not become a visible index or label.
        if markers.is_none()
            && arguments.iter().any(|argument| {
                tags::is_chain_marker(*argument, SPENSO_TAG.chain_in)
                    || tags::is_chain_marker(*argument, SPENSO_TAG.chain_out)
            })
        {
            return None;
        }
        let open_columns = self.open_columns(arguments)?;
        for ((position, argument), open_column) in
            arguments.iter().copied().enumerate().zip(open_columns)
        {
            if markers
                .is_some_and(|markers| position == markers.input || position == markers.output)
            {
                continue;
            }
            if let Some(column) = open_column {
                columns.push(column);
            } else if let Some(slot) = tags::tensor_slot(argument) {
                let source = self.index(slot.representation, slot.index)?;
                columns.push((
                    self.qualified(source, slot.representation, slot.dimension)?,
                    slot.row,
                ));
            } else if let Some((label, rep)) = self.vector(argument) {
                let source = if self.layout == Layout::Ports {
                    self.port(rep)?
                } else {
                    self.qualified(
                        format!("<mstyle mathvariant=\"bold\">{label}</mstyle>"),
                        rep.symbol,
                        rep.dimension,
                    )?
                };
                columns.push((source, rep.row));
                if self.layout == Layout::Ports {
                    match rep.polarity {
                        CompactPolarity::Bra => bras.push(label),
                        CompactPolarity::Ket => kets.push(label),
                    }
                }
            } else {
                let source = self.visual(argument, &self.child_state())?;
                if self.settings.symbol_scripts {
                    columns.push((source, IndexRow::Bottom));
                } else {
                    ordinary.push(source);
                }
            }
        }
        if !ordinary.is_empty() {
            let separator = if self.settings.commas {
                "<mo>,</mo>".to_owned()
            } else {
                Self::space("0.15em")
            };
            base = format!("<mrow>{base}{}</mrow>", group(ordinary.join(&separator)));
        }
        let upper = columns.iter().any(|(_, row)| *row == IndexRow::Top);
        base = if self.layout == Layout::Call {
            self.call(base, &columns)
        } else {
            self.attachment(base, &columns)
        };
        let transposed =
            markers.is_some_and(|markers| markers.orientation == ChainOrientation::Transposed);
        if transposed {
            base = format!("<msup><mrow>{base}</mrow><mi mathvariant=\"normal\">T</mi></msup>");
        }
        if self.layout == Layout::Ports {
            base = Self::ports(base, &bras, &kets);
        }
        if state.in_exp_base
            && (upper || transposed || head_upper || !bras.is_empty() || !kets.is_empty())
        {
            base = group(base);
        }
        self.fragment(base)
    }

    fn ports(base: String, bras: &[String], kets: &[String]) -> String {
        let mut body = String::new();
        if !bras.is_empty() {
            body.push_str(&Self::fence(bras.join("<mo>,</mo>"), "⟨", "|"));
        }
        body.push_str(&base);
        if !kets.is_empty() {
            body.push_str(&Self::fence(kets.join("<mo>,</mo>"), "|", "⟩"));
        }
        if bras.is_empty() && kets.is_empty() {
            body
        } else {
            format!("<mrow>{body}</mrow>")
        }
    }

    fn representation(&self, atom: AtomView<'_>, state: &PrintState) -> Option<String> {
        let (symbol, dimension, index, class, row) = if let Some(slot) = tags::tensor_slot(atom) {
            (
                slot.representation,
                slot.dimension,
                Some(slot.index),
                RepresentationMetadata::from_symbol(slot.representation)?.class,
                slot.row,
            )
        } else {
            let rep = tags::representation_value(atom)?;
            (rep.symbol, rep.dimension, None, rep.class, rep.row)
        };
        let label = RepresentationMetadata::from_symbol(symbol)?
            .label
            .to_mathml()?;
        let dimension = self.visual(dimension, &self.child_state())?;
        let mut body = format!("<msub><mrow>{label}</mrow><mrow>{dimension}</mrow></msub>");
        let dual = matches!(tags::unwrap_tensor_display(atom), AtomView::Fun(wrapper) if wrapper.get_symbol() == AIND_SYMBOLS.dind)
            && class == RepresentationClass::Dualizable;
        if dual {
            body = format!("<msup>{body}<mo>*</mo></msup>");
        }
        if let Some(index) = index {
            body = self.attachment(body, &[(self.index(symbol, index)?, row)]);
        }
        Some(
            if state.in_exp_base && (dual || (index.is_some() && row == IndexRow::Top)) {
                group(body)
            } else {
                body
            },
        )
    }

    fn dot(&self, atom: AtomView<'_>, _state: &PrintState) -> Option<String> {
        let AtomView::Fun(function) = atom else {
            return None;
        };
        if function.get_nargs() != 2 {
            return None;
        }
        let mut arguments = function.iter();
        let (a, left) = self.vector(arguments.next()?)?;
        let (b, right) = self.vector(arguments.next()?)?;
        if left.symbol != right.symbol
            || left.dimension != right.dimension
            || !matches!(
                left.class,
                RepresentationClass::SelfDual | RepresentationClass::InlineMetric
            )
        {
            return None;
        }
        Some(group(format!(
            "<mrow>{a}<mo lspace=\"0\" rspace=\"0\">·</mo>{b}</mrow>"
        )))
    }

    fn chain(&self, atom: AtomView<'_>, state: &PrintState) -> Option<String> {
        let AtomView::Fun(function) = atom else {
            return None;
        };
        let arguments = function.iter().collect::<Vec<_>>();
        let [start, end, factors @ ..] = arguments.as_slice() else {
            return None;
        };
        let mut body = self.factors(factors)?;
        if self.settings.parens {
            body = Self::fence(body, "[", "]");
        }
        let mut columns = Vec::new();
        let mut labels = [None, None];
        let open = self.open_columns(&[*start, *end])?;
        for (position, (endpoint, open_column)) in [*start, *end].into_iter().zip(open).enumerate()
        {
            if tags::is_chain_marker(endpoint, SPENSO_TAG.chain_in)
                || tags::is_chain_marker(endpoint, SPENSO_TAG.chain_out)
            {
                continue;
            }
            if let Some((label, rep)) = self.vector(endpoint) {
                labels[position] = Some(label);
                columns.push((self.port(rep)?, rep.row));
            } else if let Some(column) = open_column {
                columns.push(column);
            } else if let Some(slot) = tags::tensor_slot(endpoint) {
                let source = self.index(slot.representation, slot.index)?;
                columns.push((
                    self.qualified(source, slot.representation, slot.dimension)?,
                    slot.row,
                ));
            } else {
                labels[position] = Some(self.visual(endpoint, &self.child_state())?);
            }
        }
        body = self.attachment(body, &columns);
        body = Self::ports(body, labels[0].as_slice(), labels[1].as_slice());
        Some(if state.in_exp_base { group(body) } else { body })
    }

    fn trace(&self, atom: AtomView<'_>, _state: &PrintState) -> Option<String> {
        let AtomView::Fun(function) = atom else {
            return None;
        };
        let arguments = function.iter().collect::<Vec<_>>();
        let [representation, payload @ ..] = arguments.as_slice() else {
            return None;
        };
        let factors = match payload {
            [] => Vec::new(),
            [AtomView::Fun(cyclic)] if cyclic.get_symbol() == *CYCLIC => {
                cyclic.iter().collect::<Vec<_>>()
            }
            [projector @ AtomView::Fun(function)] if function.get_symbol() == *SYM => {
                vec![*projector]
            }
            _ => return None,
        };
        let mut head = "<mi mathvariant=\"normal\">Tr</mi>".to_owned();
        if self.settings.with_dim {
            head = format!(
                "<msub>{head}<mrow>{}</mrow></msub>",
                self.visual(*representation, &self.child_state())?
            );
        }
        if factors.is_empty() {
            return Some(head);
        }
        Some(format!(
            "<mrow>{head}{}</mrow>",
            group(self.factors(&factors)?)
        ))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::network::tags::tensor_print;
    use crate::structure::dimension::Dimension;
    use crate::structure::representation::{IndexPalette, LibraryRep, Minkowski};
    use crate::{dind, euc, lor, mink, tensor_symbol, vector_symbol};
    use ahash::HashMap;
    use std::sync::atomic::{AtomicUsize, Ordering};
    use std::sync::{Arc, Mutex};
    use symbolica::atom::{FunctionBuilder, NamespacedSymbol, SymbolAttribute, SymbolBuilder};
    use symbolica::{function, symbol, wrap_symbol};

    fn options() -> PrintOptions {
        PrintOptions {
            mode: PrintMode::MathML,
            hide_all_namespaces: true,
            ..PrintOptions::new()
        }
    }

    fn scoped(atom: &Atom, options: &PrintOptions) -> String {
        render(atom.as_view(), options, &PrintState::new(), &|_, _| None).unwrap()
    }

    #[test]
    fn static_expression_keeps_all_factors_past_native_arithmetic_and_page_limits() {
        let head = tensor_symbol!("static_mathml::Wide");
        let index = symbol!("static_mathml::wide_index");
        let mut value = Atom::num(1);
        for label in 0..320 {
            value *= function!(head, label, mink!(4, Atom::var(index)));
        }
        let original = value.as_view().get_data().to_vec();
        let mut options = options();
        options.max_terms = Some(1);
        let html = expression(value.as_view(), &options, &PrintState::new(), &|_, _| {
            Some("<mi mathcolor=\"red\">α</mi>".to_owned())
        })
        .unwrap();
        assert_eq!(html.matches("<mi>Wide</mi>").count(), 320);
        assert_eq!(html.matches("<mi mathcolor=\"red\">α</mi>").count(), 320);
        for label in 0..320 {
            assert!(
                html.contains(&format!("<mn>{label}</mn>")),
                "missing {label}"
            );
        }
        assert!(html.len() > MAX_BYTES);
        assert!(!html.contains('⋯') && !html.contains("data-symbolica-"));
        assert_eq!(value.as_view().get_data(), original);
        assert!(render(value.as_view(), &options, &PrintState::new(), &|_, _| None).is_none());
    }

    #[test]
    fn static_expression_recurses_through_unsupported_calls_without_losing_tensors() {
        let head = tensor_symbol!("static_mathml::Nested");
        let plain = symbol!("static_mathml::outer");
        let tensor = function!(head, mink!(4, symbol!("static_mathml::mu")));
        let mut value = FunctionBuilder::new(plain)
            .add_arg(&tensor)
            .add_args((0..80).map(Atom::num))
            .finish();
        for _ in 0..40 {
            value = function!(plain, value);
        }
        let original = value.as_view().get_data().to_vec();
        let html = expression(value.as_view(), &options(), &PrintState::new(), &|_, _| {
            None
        })
        .unwrap();
        assert_eq!(html.matches("<mi>outer</mi>").count(), 41);
        assert!(html.contains("<mmultiscripts>") && html.contains("<mi>Nested</mi>"));
        assert!(html.contains("<mi>mu</mi>"));
        for label in 0..80 {
            assert!(
                html.contains(&format!("<mn>{label}</mn>")),
                "missing {label}"
            );
        }
        assert!(!html.contains('⋯') && !html.contains("data-symbolica-"));
        assert_eq!(value.as_view().get_data(), original);

        // Declined semantic calls use the same recursive structural fallback.
        let malformed = function!(SPENSO_TAG.trace, mink!(4), &tensor, &tensor);
        let html = expression(
            malformed.as_view(),
            &options(),
            &PrintState::new(),
            &|_, _| None,
        )
        .unwrap();
        assert_eq!(html.matches("<mi>Nested</mi>").count(), 2);
        assert!(html.contains("<mi>mink</mi>") && !html.contains('⋯'));
    }

    #[test]
    fn static_expression_uses_native_single_and_multivariate_derivative_notation() {
        let f = symbol!("static_mathml::DerivativeFunction");
        let x = symbol!("static_mathml::derivative_x");
        let y = symbol!("static_mathml::derivative_y");
        let z = symbol!("static_mathml::derivative_z");
        let single = function!(f, Atom::var(x)).derivative(x);
        let multiple = function!(
            Symbol::DERIVATIVE,
            2,
            0,
            3,
            Atom::var(f),
            Atom::var(x),
            Atom::var(y),
            Atom::var(z)
        );
        for (value, orders, arguments) in [
            (single, vec![1], vec![x]),
            (multiple, vec![2, 0, 3], vec![x, y, z]),
        ] {
            assert_eq!(value.as_view().get_symbol(), Some(Symbol::DERIVATIVE));
            let original = value.as_view().get_data().to_vec();
            let html = expression(value.as_view(), &options(), &PrintState::new(), &|_, _| {
                None
            })
            .unwrap();
            assert!(
                html.starts_with("<mrow><msup><mrow><mi>DerivativeFunction</mi></mrow>"),
                "{html}"
            );
            let (script, arguments_html) = html.split_once("</msup><mo>&#x2061;</mo>").unwrap();
            let order_positions = orders
                .iter()
                .map(|order| script.find(&format!("<mn>{order}</mn>")).unwrap())
                .collect::<Vec<_>>();
            let argument_positions = arguments
                .iter()
                .map(|argument| {
                    arguments_html
                        .find(&format!("<mi>{}</mi>", argument.get_stripped_name()))
                        .unwrap()
                })
                .collect::<Vec<_>>();
            assert!(
                order_positions.windows(2).all(|pair| pair[0] < pair[1]),
                "{html}"
            );
            assert!(
                argument_positions.windows(2).all(|pair| pair[0] < pair[1]),
                "{html}"
            );
            assert_eq!(script.matches("<mo>,</mo>").count(), orders.len() - 1);
            assert_eq!(
                arguments_html.matches("<mo>,</mo>").count(),
                arguments.len() - 1
            );
            assert!(!html.contains("<mi>der</mi>") && !html.contains('⋯'));
            assert_eq!(value.as_view().get_data(), original);
            let mut state = PrintState::new();
            state.in_exp_base = true;
            assert_eq!(
                expression(value.as_view(), &options(), &state, &|_, _| None),
                Ok(group(html))
            );
        }
    }

    #[test]
    fn static_expression_derivatives_keep_rich_orders_heads_arguments_and_aliases() {
        let head = SymbolBuilder::new(wrap_symbol!("static_mathml::DerivativeHead"))
            .with_tags([
                SPENSO_TAG.tensor.as_str(),
                "spenso::print-head:mathml:<mi mathcolor=\"purple\">𝒜</mi>",
            ])
            .with_print_function(tensor_print)
            .build()
            .unwrap();
        let tensor = tensor_symbol!("static_mathml::DerivativeArgument");
        let index = symbol!("static_mathml::derivative_index");
        let rich = function!(tensor, mink!(4, Atom::var(index)));
        // The native derivative shape admits any symbolic head. Its source
        // children and non-scalar orders must also survive static rendering.
        let value = function!(Symbol::DERIVATIVE, &rich, function!(head, &rich), &rich);
        let original = value.as_view().get_data().to_vec();
        let html = expression(
            value.as_view(),
            &options(),
            &PrintState::new(),
            &|_, original_index| {
                assert_eq!(original_index, Atom::var(index).as_view());
                Some("<mi mathcolor=\"orange\">α</mi>".to_owned())
            },
        )
        .unwrap();
        assert!(html.contains("<mi mathcolor=\"purple\">𝒜</mi>"), "{html}");
        assert_eq!(html.matches("<mi>DerivativeArgument</mi>").count(), 3);
        assert_eq!(html.matches("<mi mathcolor=\"orange\">α</mi>").count(), 3);
        assert!(html.contains("<mmultiscripts>") && html.contains("<mo>&#x2061;</mo>"));
        assert!(!html.contains('⋯') && !html.contains("data-symbolica-"));
        assert_eq!(value.as_view().get_data(), original);
    }

    #[test]
    fn static_expression_invalid_derivative_shapes_keep_structural_arguments() {
        let f = symbol!("static_mathml::InvalidDerivativeFunction");
        let same_name = symbol!("static_mathml::der");
        let tensor = tensor_symbol!("static_mathml::InvalidDerivativeArgument");
        let rich = function!(
            tensor,
            mink!(4, symbol!("static_mathml::invalid_derivative_index"))
        );
        for (value, occurrences) in [
            (function!(Symbol::DERIVATIVE), 0),
            (function!(Symbol::DERIVATIVE, &rich), 1),
            (function!(Symbol::DERIVATIVE, 1, &rich), 1),
            (function!(Symbol::DERIVATIVE, 1, 2, &rich, &rich), 2),
            (function!(Symbol::DERIVATIVE, 1, 7, &rich), 1),
            (function!(same_name, 1, Atom::var(f), &rich), 1),
        ] {
            let original = value.as_view().get_data().to_vec();
            let html = expression(value.as_view(), &options(), &PrintState::new(), &|_, _| {
                None
            })
            .unwrap();
            assert!(html.contains("<mi>der</mi>"), "{html}");
            assert!(
                !html.contains("<mo>&#x2061;</mo>") && !html.contains('⋯'),
                "{html}"
            );
            assert_eq!(
                html.matches("<mi>InvalidDerivativeArgument</mi>").count(),
                occurrences
            );
            assert_eq!(value.as_view().get_data(), original);
        }
    }

    #[test]
    fn static_expression_keeps_numeric_signs_and_fences_products_and_power_bases() {
        let options = options();
        let exponent = Atom::var(symbol!("static_mathml::exponent"));
        let factor = Atom::var(symbol!("static_mathml::factor"));
        let large = Atom::num(2).pow(80);
        for scalar in [
            Atom::num(-2),
            Atom::num((-2, 3)),
            Atom::num(-0.5),
            Atom::num((-2, 3)) + Atom::i() * Atom::num((5, 7)),
            Atom::num((2, 3)) - Atom::i() * Atom::num((5, 7)),
            Atom::num((-2, 3)) - Atom::i() * Atom::num((5, 7)),
            -Atom::i(),
            large - Atom::i() * Atom::num((11, 13)),
        ] {
            assert!(matches!(scalar.as_view(), AtomView::Num(_)));
            let native = source(scalar.as_view(), &options).unwrap();
            let html =
                expression(scalar.as_view(), &options, &PrintState::new(), &|_, _| None).unwrap();
            assert_eq!(html, native);
            let product = &factor * &scalar;
            let original = product.as_view().get_data().to_vec();
            let html = expression(product.as_view(), &options, &PrintState::new(), &|_, _| {
                None
            })
            .unwrap();
            assert!(html.contains(&group(native.clone())), "{html}");
            assert_eq!(product.as_view().get_data(), original);
            let power = scalar.pow(&exponent);
            let html =
                expression(power.as_view(), &options, &PrintState::new(), &|_, _| None).unwrap();
            assert!(
                html.contains(&format!("<msup><mrow>{}</mrow>", group(native))),
                "{html}"
            );
        }
        let positive = Atom::num((1, 2)).pow(&exponent);
        let html = expression(positive.as_view(), &options, &PrintState::new(), &|_, _| {
            None
        })
        .unwrap();
        assert!(html.contains(&group(
            source(Atom::num((1, 2)).as_view(), &options).unwrap()
        )));
    }

    #[test]
    fn static_expression_preserves_callback_options_state_and_shared_leaf_notation() {
        let states = Arc::new(Mutex::new(Vec::new()));
        let callback_states = Arc::clone(&states);
        let callback = SymbolBuilder::new(wrap_symbol!("static_mathml::callback"))
            .with_print_function(move |_, options, state| {
                assert!(options.mode.is_mathml());
                assert!(matches!(
                    options.custom_print_mode.get("static-test"),
                    Some(PrintUserData::Integer(17))
                ));
                callback_states.lock().unwrap().push(*state);
                Some("<mi data-symbolica-path=\"hidden\">custom</mi>".to_owned())
            })
            .build()
            .unwrap();
        let mut options = options();
        options
            .custom_print_mode
            .insert("static-test".to_owned(), PrintUserData::Integer(17));
        let mut state = PrintState::new();
        state.level = 17;
        state.in_exp_base = true;
        state.in_product = true;
        let value = Atom::var(callback);
        assert_eq!(
            expression(value.as_view(), &options, &state, &|_, _| None).as_deref(),
            Ok("<mi>custom</mi>")
        );
        let seen = states.lock().unwrap();
        assert_eq!(seen.len(), 1);
        assert_eq!(seen[0].level, 17);
        assert!(seen[0].in_exp_base && seen[0].in_product);
        drop(seen);
        states.lock().unwrap().clear();
        let nested = function!(SPENSO_TAG.bracket, value.pow(2));
        expression(nested.as_view(), &options, &state, &|_, _| None).unwrap();
        let seen = states.lock().unwrap();
        assert_eq!(seen.len(), 1);
        assert_eq!(seen[0].level, 19);
        assert!(seen[0].in_exp_base && !seen[0].in_product);
        drop(seen);
        states.lock().unwrap().clear();
        let sum = &value + Atom::num(1);
        expression(sum.as_view(), &options, &state, &|_, _| None).unwrap();
        let seen = states.lock().unwrap();
        assert_eq!(seen.len(), 1);
        assert!(seen[0].in_sum && seen[0].level == 18);
        drop(seen);

        let head = tensor_symbol!("static_mathml::Shared");
        let value = function!(head, 3, mink!(4, symbol!("static_mathml::shared_index")));
        for layout in ["ports", "schoonschip", "call"] {
            options.custom_print_mode.insert(
                "spenso.tensor_layout".to_owned(),
                PrintUserData::String(layout.to_owned()),
            );
            assert_eq!(
                expression(value.as_view(), &options, &state, &|_, _| None).unwrap(),
                render(value.as_view(), &options, &state, &|_, _| None).unwrap()
            );
        }
    }

    #[test]
    fn static_expression_declined_scalar_callback_runs_once_per_render() {
        let calls = Arc::new(AtomicUsize::new(0));
        let callback_calls = Arc::clone(&calls);
        let symbol = SymbolBuilder::new(wrap_symbol!("static_mathml::declined_scalar"))
            .with_attributes(vec![SymbolAttribute::Scalar])
            .with_tags([TENSOR_PRINT_CALLBACK_TAG])
            .with_print_function(move |atom, options, state| {
                assert!(matches!(atom, AtomView::Var(_)));
                assert!(options.mode.is_mathml());
                assert!(state.in_exp_base);
                assert_eq!(state.level, 23);
                (callback_calls.fetch_add(1, Ordering::SeqCst) != 0)
                    .then(|| "<mi>later-render</mi>".to_owned())
            })
            .build()
            .unwrap();
        let value = Atom::var(symbol);
        let original = value.as_view().get_data().to_vec();
        let mut state = PrintState::new();
        state.level = 23;
        state.in_exp_base = true;
        assert_eq!(
            expression(value.as_view(), &options(), &state, &|_, _| None).as_deref(),
            Ok("<mi>declined_scalar</mi>")
        );
        assert_eq!(calls.load(Ordering::SeqCst), 1);
        assert_eq!(
            expression(value.as_view(), &options(), &state, &|_, _| None).as_deref(),
            Ok("<mi>later-render</mi>")
        );
        assert_eq!(calls.load(Ordering::SeqCst), 2);
        assert_eq!(value.as_view().get_data(), original);
    }

    #[test]
    fn static_expression_declined_function_callbacks_keep_source_namespaces_once() {
        let calls = Arc::new(AtomicUsize::new(0));
        let symbols = [
            ("static_mathml_declined_left::γ", true),
            ("static_mathml_declined_right::γ", false),
        ]
        .map(|(name, tagged)| {
            let calls = Arc::clone(&calls);
            SymbolBuilder::new(NamespacedSymbol::parse(name))
                .with_tags(if tagged {
                    vec![TENSOR_PRINT_CALLBACK_TAG]
                } else {
                    Vec::new()
                })
                .with_print_function(move |atom, options, state| {
                    assert!(matches!(atom, AtomView::Fun(function) if function.get_nargs() == 1));
                    assert!(options.mode.is_mathml());
                    assert_eq!(state.level, 19);
                    calls.fetch_add(1, Ordering::SeqCst);
                    None
                })
                .build()
                .unwrap()
        });
        let argument = Atom::var(symbol!("static_mathml::declined_function_argument"));
        let mut state = PrintState::new();
        state.level = 19;
        for (hide_all_namespaces, hide_namespace) in [
            (false, None),
            (true, None),
            (false, Some("static_mathml_declined_left")),
            (false, Some("static_mathml_declined_right")),
        ] {
            let options = PrintOptions {
                hide_all_namespaces,
                hide_namespace: hide_namespace.map(Into::into),
                ..options()
            };
            let rendered = symbols.map(|symbol| {
                calls.store(0, Ordering::SeqCst);
                let value = function!(symbol, &argument);
                let original = value.as_view().get_data().to_vec();
                let html = expression(value.as_view(), &options, &state, &|_, _| None).unwrap();
                assert_eq!(calls.load(Ordering::SeqCst), 1);
                let qualified =
                    !hide_all_namespaces && hide_namespace != Some(symbol.get_namespace());
                let label = if qualified {
                    symbol.get_name()
                } else {
                    symbol.get_stripped_name()
                };
                assert!(html.contains(&format!("<mi>{label}</mi>")), "{html}");
                assert!(
                    html.contains("<mi>declined_function_argument</mi>"),
                    "{html}"
                );
                assert_eq!(html.contains(symbol.get_namespace()), qualified, "{html}");
                assert_eq!(value.as_view().get_data(), original);
                html
            });
            if hide_all_namespaces {
                assert_eq!(rendered[0], rendered[1]);
            } else {
                assert_ne!(rendered[0], rendered[1]);
            }
        }
    }

    #[test]
    fn static_expression_keeps_native_greek_names_namespaces_and_numbers() {
        for hide_all_namespaces in [false, true] {
            let options = PrintOptions {
                hide_all_namespaces,
                ..options()
            };
            for value in [
                Atom::var(symbol!("static_mathml::α")),
                Atom::var(symbol!("static_mathml::γ")),
                Atom::var(symbol!("static_mathml::native_scalar")),
                Atom::var(Symbol::PI),
                Atom::num(0),
                Atom::num(-2),
                Atom::num((1, 3)),
                Atom::num(2.5),
                Atom::num((2, 3)) - Atom::i() * Atom::num((5, 7)),
            ] {
                let native = source(value.as_view(), &options).unwrap();
                assert_eq!(
                    expression(value.as_view(), &options, &PrintState::new(), &|_, _| None),
                    Ok(native)
                );
            }
        }
    }

    #[test]
    fn static_expression_declines_limits_without_returning_partial_output() {
        let plain = symbol!("static_mathml::depth");
        let mut value = Atom::var(plain);
        for _ in 0..=MAX_EXPORT_DEPTH {
            value = function!(plain, value);
        }
        assert_eq!(
            expression(value.as_view(), &options(), &PrintState::new(), &|_, _| {
                None
            }),
            Err(ExportError::SourceLimit)
        );
        assert!(!Render::expression_fits(
            function!(plain, 1, 2).as_view(),
            &mut 2,
            0
        ));

        let callback = SymbolBuilder::new(wrap_symbol!("static_mathml::oversized"))
            .with_tags([SPENSO_TAG.tensor.as_str(), TENSOR_PRINT_CALLBACK_TAG])
            .with_print_function(|_, _, _| {
                Some(format!("<mtext>{}</mtext>", "x".repeat(MAX_EXPORT_BYTES)))
            })
            .build()
            .unwrap();
        let value = function!(callback, 1);
        let original = value.as_view().get_data().to_vec();
        assert_eq!(
            expression(value.as_view(), &options(), &PrintState::new(), &|_, _| {
                None
            }),
            Err(ExportError::IncompleteOutput)
        );
        assert_eq!(value.as_view().get_data(), original);
        assert!(fragment(format!("<mtext>{}</mtext>", "x".repeat(MAX_BYTES))).is_none());
    }

    #[test]
    fn static_expression_reports_invalid_settings_and_retains_parse_errors() {
        let value = Atom::var(tensor_symbol!("static_mathml::InvalidSettings"));
        let state = PrintState::new();
        for (key, setting, parse_error) in [
            ("spenso.tensor_layout", "unknown", false),
            ("spenso.index_gap", "1unknown", false),
            ("spenso.factor_gap", "bad-em", true),
            ("spenso.index_gap", "infem", false),
        ] {
            let mut options = options();
            options
                .custom_print_mode
                .insert(key.to_owned(), PrintUserData::String(setting.to_owned()));
            let error = expression(value.as_view(), &options, &state, &|_, _| None).unwrap_err();
            let ExportError::InvalidSettings {
                key: actual,
                source,
            } = error
            else {
                panic!("expected invalid settings for {key}={setting}");
            };
            assert_eq!(actual, key);
            assert_eq!(source.is_some(), parse_error, "parse source for {setting}");
            assert!(render(value.as_view(), &options, &state, &|_, _| None).is_none());
        }
    }

    fn gamma() -> Symbol {
        tags::register_tensor_symbol(
            NamespacedSymbol::parse("spenso::gamma"),
            vec![SymbolAttribute::Linear],
            false,
        )
        .unwrap()
    }

    #[test]
    fn fragment_removes_only_opening_tag_metadata_and_preserves_quoted_text() {
        let value = r#"<mrow title="quoted > data-symbolica-path='keep'" data-symbolica-path = "0>1" class='kept' data-symbolica-hidden='yes'><mtext>literal data-symbolica-path="keep" text</mtext><annotation encoding="text/plain">literal data-symbolica-hidden='keep'</annotation><mo title='data-symbolica-path="keep" >'>&gt;</mo><annotation><![CDATA[<mrow data-symbolica-path="keep">]]></annotation></mrow>"#;
        let expected = r#"<mrow title="quoted > data-symbolica-path='keep'" class='kept'><mtext>literal data-symbolica-path="keep" text</mtext><annotation encoding="text/plain">literal data-symbolica-hidden='keep'</annotation><mo title='data-symbolica-path="keep" >'>&gt;</mo><annotation><![CDATA[<mrow data-symbolica-path="keep">]]></annotation></mrow>"#;
        assert_eq!(fragment(value.to_owned()).as_deref(), Some(expected));
        assert!(fragment(format!("<mtext>{}</mtext>", "x".repeat(MAX_BYTES))).is_none());
    }

    #[test]
    fn scoped_bracket_fences_numeric_power_bases_using_native_mathml_rule() {
        let options = options();
        let exponent = Atom::var(symbol!("scoped_mathml::numeric_exponent"));
        for (base, fenced) in [
            (Atom::num(-2), true),
            (Atom::num((1, 2)), true),
            (Atom::num((-1, 2)), true),
            (Atom::num(2), false),
        ] {
            let value = function!(SPENSO_TAG.bracket, base.pow(&exponent));
            let original = value.as_view().get_data().to_vec();
            let mut base = source(base.as_view(), &options).unwrap();
            if fenced {
                base = group(base);
            }
            let html = scoped(&value, &options);
            assert!(
                html.contains(&format!("<msup><mrow>{base}</mrow>")),
                "{html}"
            );
            assert_eq!(value.as_view().get_data(), original);
        }
    }

    #[test]
    fn scoped_metrics_accept_compatible_open_and_explicit_ports_only() {
        let rep = LibraryRep::new_dual("scoped_mathml::metric_open_space").unwrap();
        let open = |axis| function!(AIND_SYMBOLS.openind, 17, axis);
        let valid = [
            (mink!(4), mink!(4), 2),
            (euc!(3), euc!(3), 2),
            (
                function!(rep.symbol(), 5),
                dind!(function!(rep.symbol(), 5)),
                2,
            ),
            (mink!(4, open(0)), mink!(4, open(1)), 2),
            (mink!(4, symbol!("scoped_mathml::metric_i")), mink!(4), 1),
            (
                function!(rep.symbol(), 5, symbol!("scoped_mathml::metric_j")),
                dind!(function!(rep.symbol(), 5)),
                1,
            ),
        ];
        for (first, second, ports) in valid {
            let value = ETS.metric(first, second);
            let original = value.as_view().get_data().to_vec();
            let html = scoped(&value, &options());
            assert_eq!(
                html.matches("data-spenso-open-axis=").count(),
                ports,
                "{html}"
            );
            assert!(
                html.contains("data-spenso-open-axis=\"0\"")
                    || html.contains("data-spenso-open-axis=\"1\""),
                "{html}"
            );
            assert!(html.contains("title=\"Unresolved axis"), "{html}");
            if ports == 2 {
                assert!(
                    html.contains("data-spenso-open-axis=\"0\"")
                        && html.contains("data-spenso-open-axis=\"1\""),
                    "{html}"
                );
            }
            assert_eq!(value.as_view().get_data(), original);
        }
        for (first, second) in [
            (mink!(4), euc!(4)),
            (mink!(4), mink!(5)),
            (function!(rep.symbol(), 5), function!(rep.symbol(), 5)),
            (
                dind!(function!(rep.symbol(), 5)),
                dind!(function!(rep.symbol(), 5)),
            ),
            (
                function!(rep.symbol(), 5),
                dind!(function!(rep.symbol(), 6)),
            ),
        ] {
            let value = ETS.metric(first, second);
            assert!(
                render(value.as_view(), &options(), &PrintState::new(), &|_, _| {
                    None
                })
                .is_none()
            );
        }
        let open_metric = ETS.metric(mink!(4), mink!(4));
        assert!(tensor(open_metric.as_view(), &options(), &PrintState::new()).is_none());
    }

    #[test]
    fn scoped_gamma_consumes_only_lorentz_port_and_styles_its_slash() {
        let spinor = LibraryRep::new_self_dual_with_index_palette_and_row(
            "scoped_mathml::spinor",
            IndexPalette::Numeric,
            IndexRow::Bottom,
        )
        .unwrap();
        let p = vector_symbol!("scoped_mathml::p");
        let momentum = function!(p, 7, mink!(4));
        for (left, right, ports) in [
            (
                function!(spinor.symbol(), 4),
                function!(spinor.symbol(), 4),
                2,
            ),
            (
                function!(spinor.symbol(), 4, symbol!("scoped_mathml::a")),
                function!(spinor.symbol(), 4),
                1,
            ),
            (
                Atom::var(SPENSO_TAG.chain_in),
                Atom::var(SPENSO_TAG.chain_out),
                0,
            ),
        ] {
            let value = function!(gamma(), left, right, &momentum);
            let original = value.as_view().get_data().to_vec();
            let html = scoped(&value, &options());
            assert_eq!(html.matches('□').count(), ports, "{html}");
            assert!(html.contains("data-spenso-cancel=\"updiagonal\""), "{html}");
            assert!(html.contains("linear-gradient"), "{html}");
            assert!(html.contains("currentColor"), "{html}");
            assert!(html.contains("<mn>7</mn>"), "{html}");
            assert!(!html.contains("data-symbolica-"), "{html}");
            assert!(
                !html.contains("<mi>in</mi>") && !html.contains("<mi>out</mi>"),
                "{html}"
            );
            assert_eq!(value.as_view().get_data(), original);
        }
        let indexed = function!(
            gamma(),
            function!(spinor.symbol(), 4, symbol!("scoped_mathml::a")),
            function!(spinor.symbol(), 4, symbol!("scoped_mathml::b")),
            mink!(4, symbol!("scoped_mathml::mu"))
        );
        let html = scoped(&indexed, &options());
        assert!(html.contains("<mi>γ</mi>"), "{html}");
        for name in ["a", "b", "mu"] {
            assert!(html.contains(&format!("<mi>{name}</mi>")), "{html}");
        }
        assert!(!html.contains("data-spenso-cancel"), "{html}");
    }

    #[test]
    fn scoped_tensor_layouts_keep_bra_ket_polarity_and_representation_rows() {
        let rep = LibraryRep::new_dual("scoped_mathml::ports_space").unwrap();
        let u = vector_symbol!("scoped_mathml::u");
        let v = vector_symbol!("scoped_mathml::v");
        let value = function!(
            tensor_symbol!("scoped_mathml::A"),
            function!(u, 1, dind!(function!(rep.symbol(), 3))),
            function!(v, 2, function!(rep.symbol(), 3))
        );
        let mut options = options();
        let html = scoped(&value, &options);
        assert!(
            html.contains("<mo>⟨</mo>") && html.contains("<mo>⟩</mo>"),
            "{html}"
        );
        assert!(
            html.contains("<mo>▶</mo>") && html.contains("<mo>◀</mo>"),
            "{html}"
        );
        assert!(html.contains("<mrow><mo>▶</mo></mrow><none/>"), "{html}");
        assert!(html.contains("<none/><mrow><mo>◀</mo></mrow>"), "{html}");
        for layout in ["schoonschip", "call"] {
            options.custom_print_mode.insert(
                "spenso.tensor_layout".to_owned(),
                PrintUserData::String(layout.to_owned()),
            );
            let html = scoped(&value, &options);
            assert!(html.contains("mathvariant=\"bold\""), "{html}");
            assert!(
                !html.contains('⟨')
                    && !html.contains('⟩')
                    && !html.contains('▶')
                    && !html.contains('◀'),
                "{html}"
            );
            if layout == "call" {
                assert!(html.contains("<mo>;</mo>"), "{html}");
            }
        }
    }

    #[test]
    fn scoped_chain_trace_and_transposition_hide_internal_wiring() {
        let mu = symbol!("scoped_mathml::mu");
        let nu = symbol!("scoped_mathml::nu");
        let factor = function!(
            gamma(),
            Atom::var(SPENSO_TAG.chain_in),
            Atom::var(SPENSO_TAG.chain_out),
            mink!(4, Atom::var(mu))
        );
        let reversed = function!(
            gamma(),
            Atom::var(SPENSO_TAG.chain_out),
            Atom::var(SPENSO_TAG.chain_in),
            mink!(4, Atom::var(nu))
        );
        let value = function!(
            SPENSO_TAG.chain,
            mink!(4, symbol!("scoped_mathml::i")),
            mink!(4, symbol!("scoped_mathml::j")),
            &factor,
            &reversed
        );
        let mut options = options();
        options.custom_print_mode.insert(
            "spenso.factor_gap".to_owned(),
            PrintUserData::String("20pt".to_owned()),
        );
        let html = scoped(&value, &options);
        assert!(
            html.contains("<mo>[</mo>") && html.contains("<mo>]</mo>"),
            "{html}"
        );
        assert!(html.contains("mathvariant=\"normal\">T</mi>"), "{html}");
        assert!(html.contains("width=\"20pt\""), "{html}");
        for name in ["i", "j", "mu", "nu"] {
            assert!(html.contains(&format!("<mi>{name}</mi>")), "{html}");
        }
        assert!(
            !html.contains("<mi>in</mi>") && !html.contains("<mi>out</mi>"),
            "{html}"
        );
        let trace = function!(
            SPENSO_TAG.trace,
            mink!(4),
            function!(*CYCLIC, &factor, &reversed)
        );
        let html = scoped(&trace, &options);
        assert!(html.contains(">Tr</mi>"), "{html}");
        assert!(
            !html.contains("cyclic") && !html.contains("<mi>in</mi>"),
            "{html}"
        );
        assert!(html.contains("width=\"20pt\""), "{html}");
        let malformed = function!(SPENSO_TAG.trace, mink!(4), &factor, &reversed);
        assert!(
            render(
                malformed.as_view(),
                &options,
                &PrintState::new(),
                &|_, _| None
            )
            .is_none()
        );
    }

    #[test]
    fn scoped_aliases_receive_original_indices_and_symbolic_dimensions_survive() {
        let head = tensor_symbol!("scoped_mathml::Aliases");
        let original = symbol!("scoped_mathml::original_index");
        let dimension = symbol!("scoped_mathml::D");
        let value = function!(
            head,
            mink!(Dimension::Symbolic(dimension.into()), Atom::var(original))
        );
        let mut options = options();
        let mut settings = SpensoPrintSettings::typst();
        settings.with_dim = true;
        options.custom_print_mode.extend(HashMap::from(settings));
        let calls = Cell::new(0);
        let html = render(
            value.as_view(),
            &options,
            &PrintState::new(),
            &|rep, index| {
                assert_eq!(rep, LibraryRep::from(Minkowski {}).symbol());
                assert_eq!(index, Atom::var(original).as_view());
                calls.set(calls.get() + 1);
                Some("<msub><mi>α</mi><mn>9</mn></msub>".to_owned())
            },
        )
        .unwrap();
        assert_eq!(calls.get(), 1);
        assert!(
            html.contains("<mi>α</mi>") && html.contains("<mi>D</mi>"),
            "{html}"
        );
        assert!(
            !html.contains("original_index") && !html.contains("data-symbolica-"),
            "{html}"
        );
    }

    #[test]
    fn scoped_components_scripts_grouping_and_all_spacing_units() {
        let head = SymbolBuilder::new(wrap_symbol!("scoped_mathml::component"))
            .with_tags([SPENSO_TAG.tensor.as_str(), COMPONENT_LABELS_SUBSCRIPT_TAG])
            .with_print_function(tensor_print)
            .build()
            .unwrap();
        let component = function!(head, 7, function!(AIND_SYMBOLS.cind, 2, 3));
        let mut options = options();
        let mut state = PrintState::new();
        let html = scoped(&component, &options);
        assert!(
            html.contains("<msub><mi>component</mi>") && html.contains("<msup>"),
            "{html}"
        );
        state.in_exp_base = true;
        assert_eq!(
            render(component.as_view(), &options, &state, &|_, _| None),
            Some(group(html))
        );
        let mut settings = SpensoPrintSettings::typst();
        settings.array_components = true;
        options.custom_print_mode.extend(HashMap::from(settings));
        let html = scoped(&component, &options);
        assert!(
            html.contains("<mo>[</mo>") && !html.contains("<msup>"),
            "{html}"
        );
        let indexed = function!(
            head,
            mink!(4, symbol!("scoped_mathml::i")),
            mink!(4, symbol!("scoped_mathml::j"))
        );
        for gap in ["-0.02em", "20pt", "1mm", "2cm", "0.1in", "25%"] {
            options.custom_print_mode.insert(
                "spenso.index_gap".to_owned(),
                PrintUserData::String(gap.to_owned()),
            );
            assert!(
                scoped(&indexed, &options).contains(&format!("width=\"{gap}\"")),
                "{gap}"
            );
        }
        settings.symbol_scripts = false;
        settings.commas = true;
        settings.parens = false;
        options.custom_print_mode.extend(HashMap::from(settings));
        let labels = function!(head, 5, 6, mink!(4, symbol!("scoped_mathml::i")));
        let html = scoped(&labels, &options);
        let positions = [
            "<mo>(</mo>",
            "<mn>5</mn>",
            "<mo>,</mo>",
            "<mn>6</mn>",
            "<mo>)</mo>",
        ]
        .map(|token| html.find(token).expect("ordered scalar argument token"));
        assert!(positions.windows(2).all(|pair| pair[0] < pair[1]), "{html}");
        assert!(!html.contains("<msub>"), "{html}");
        let chain = function!(
            SPENSO_TAG.chain,
            Atom::var(SPENSO_TAG.chain_in),
            Atom::var(SPENSO_TAG.chain_out),
            &labels
        );
        assert!(!scoped(&chain, &options).contains("<mo>[</mo>"));
    }

    #[test]
    fn scoped_semantic_leaves_keep_native_callbacks_and_reject_unbounded_work() {
        let callback = symbol!(
            "scoped_mathml::custom_scalar",
            print = |_, options, state| {
                assert!(options.mode.is_mathml());
                Some(
                    if state.in_exp_base {
                        "<mi>power</mi>"
                    } else {
                        "<mi>ordinary</mi>"
                    }
                    .to_owned(),
                )
            }
        );
        let mut state = PrintState::new();
        state.in_exp_base = true;
        assert_eq!(
            render(
                Atom::var(callback).as_view(),
                &options(),
                &state,
                &|_, _| None
            ),
            Some("<mi>power</mi>".to_owned())
        );
        let nested = function!(SPENSO_TAG.bracket, Atom::var(callback).pow(2));
        let nested = scoped(&nested, &options());
        assert!(nested.contains("<mi>power</mi>"), "{nested}");
        let head = tensor_symbol!("scoped_mathml::Bounded");
        let value = function!(head, mink!(4, symbol!("scoped_mathml::i")));
        let huge_alias =
            |_, _: AtomView<'_>| Some(format!("<mtext>{}</mtext>", "x".repeat(MAX_BYTES)));
        assert!(render(value.as_view(), &options(), &state, &huge_alias).is_none());
        let many = FunctionBuilder::new(head)
            .add_args((0..256).map(Atom::num))
            .finish();
        assert!(render(many.as_view(), &options(), &state, &|_, _| None).is_none());
        let coordinates = FunctionBuilder::new(AIND_SYMBOLS.cind)
            .add_args((0..70).map(Atom::num))
            .finish();
        let component = function!(head, coordinates);
        assert!(render(component.as_view(), &options(), &state, &|_, _| None).is_some());
        let plain = symbol!("scoped_mathml::nested");
        let mut deep = Atom::var(plain);
        for _ in 0..34 {
            deep = function!(plain, deep);
        }
        let deep = function!(head, deep);
        assert!(render(deep.as_view(), &options(), &state, &|_, _| None).is_none());
        let arithmetic = value + Atom::num(1);
        assert!(render(arithmetic.as_view(), &options(), &state, &|_, _| None).is_none());
        assert_eq!(
            fragment("<mrow data-symbolica-path=\"0:1\"><mi>x</mi></mrow>".to_owned()),
            Some("<mrow><mi>x</mi></mrow>".to_owned())
        );
        assert!(fragment("<math><mi>x</mi></math>".to_owned()).is_none());
    }

    #[test]
    fn scoped_oversized_tensor_callable_falls_back_once_and_keeps_arguments() {
        let calls = Arc::new(AtomicUsize::new(0));
        let callback_calls = Arc::clone(&calls);
        let head = SymbolBuilder::new(wrap_symbol!("scoped_mathml::oversized_callable"))
            .with_tags([SPENSO_TAG.tensor.as_str(), TENSOR_PRINT_CALLBACK_TAG])
            .with_print_function(move |_, options, _| {
                assert!(options.mode.is_mathml());
                callback_calls.fetch_add(1, Ordering::SeqCst);
                Some(format!("<mtext>{}</mtext>", "x".repeat(5000)))
            })
            .build()
            .unwrap();
        let value = function!(head, 1, 2);
        let original = value.as_view().get_data().to_vec();
        let html = scoped(&value, &options());
        assert_eq!(calls.load(Ordering::SeqCst), 1);
        for token in ["<mi>oversized_callable</mi>", "<mn>1</mn>", "<mn>2</mn>"] {
            assert!(html.contains(token), "{html}");
        }
        assert!(!html.contains("<mtext>") && !html.contains('⋯'), "{html}");
        assert!(html.len() <= MAX_BYTES);
        assert_eq!(value.as_view().get_data(), original);
    }

    #[test]
    fn scoped_special_heads_and_dot_power_bases_keep_their_grouping() {
        let mut state = PrintState::new();
        state.in_exp_base = true;
        for (name, label) in [
            ("gamma0", "<mn>0</mn>"),
            ("gamma5", "<mn>5</mn>"),
            ("projp", "<mo>+</mo>"),
            ("projm", "<mo>−</mo>"),
        ] {
            let symbol = SymbolBuilder::new(NamespacedSymbol::parse(&format!("spenso::{name}")))
                .with_tags([SPENSO_TAG.tensor.clone()])
                .build()
                .unwrap();
            let value = function!(
                symbol,
                Atom::var(SPENSO_TAG.chain_in),
                Atom::var(SPENSO_TAG.chain_out)
            );
            let html = scoped(&value, &options());
            assert!(html.contains(label), "{html}");
            let power = render(value.as_view(), &options(), &state, &|_, _| None).unwrap();
            if name.starts_with("proj") {
                assert_eq!(power, group(html));
            }
        }
        let p = vector_symbol!("scoped_mathml::dot_p");
        let q = vector_symbol!("scoped_mathml::dot_q");
        let value = function!(
            SPENSO_TAG.dot,
            function!(p, 1, mink!(4)),
            function!(q, 2, mink!(4))
        );
        let html = scoped(&value, &options());
        assert!(html.starts_with("<mrow><mo>(</mo>"), "{html}");
        assert!(!html.contains('□'), "{html}");
        assert_eq!(
            render(value.as_view(), &options(), &state, &|_, _| None),
            Some(html)
        );
    }

    #[test]
    fn tensor_scripts_keep_rows_arguments_and_power_grouping() {
        let head = tensor_symbol!("native_mathml_print::T");
        let value = function!(
            head,
            symbol!("native_mathml_print::parameter"),
            mink!(4, symbol!("native_mathml_print::mu")),
            dind!(lor!(4, symbol!("native_mathml_print::nu")))
        );
        let options = options();
        let mut state = PrintState::new();
        let html = tensor(value.as_view(), &options, &state).unwrap();
        assert!(html.contains("<mi>parameter</mi>"), "{html}");
        assert!(html.contains("<none/><mrow><mi>mu</mi></mrow>"), "{html}");
        assert!(html.contains("<mrow><mi>nu</mi></mrow><none/>"), "{html}");
        state.in_exp_base = true;
        assert_eq!(
            tensor(value.as_view(), &options, &state).unwrap(),
            group(html)
        );

        let oversized = FunctionBuilder::new(head)
            .add_args((0..65).map(Atom::num))
            .finish();
        assert!(tensor(oversized.as_view(), &options, &state).is_none());
    }

    #[test]
    fn dots_only_hide_compatible_contracted_representations() {
        let p = vector_symbol!("native_mathml_print::p");
        let q = vector_symbol!("native_mathml_print::q");
        let options = options();
        let mut state = PrintState::new();
        let valid = function!(
            SPENSO_TAG.dot,
            function!(p, 2, mink!(4)),
            function!(q, 3, mink!(4))
        );
        let html = dot(valid.as_view(), &options, &state).unwrap();
        assert!(html.contains("<msub><mi>p</mi>"), "{html}");
        assert!(html.contains("<mn>2</mn>"), "{html}");
        assert!(html.contains("<mo>·</mo>"), "{html}");
        assert!(!html.contains("□"), "{html}");
        for (in_product, in_exp_base) in [(true, false), (false, true), (true, true)] {
            state.in_product = in_product;
            state.in_exp_base = in_exp_base;
            assert_eq!(
                dot(valid.as_view(), &options, &state).unwrap(),
                group(html.clone())
            );
        }
        let incompatible = function!(
            SPENSO_TAG.dot,
            function!(p, 2, mink!(4)),
            function!(q, 3, mink!(5))
        );
        assert!(dot(incompatible.as_view(), &options, &state).is_none());
    }

    #[test]
    fn custom_head_fragments_remain_one_script_base() {
        let symbol = SymbolBuilder::new(wrap_symbol!("native_mathml_print::multi_root_head"))
            .with_tags([
                SPENSO_TAG.tensor.as_str(),
                SPENSO_TAG.rank1.as_str(),
                "spenso::print-head:mathml:<mi>A</mi><mo>′</mo>",
            ])
            .with_user_data(
                IndexDisplay::math("overline", vec![IndexDisplay::symbol("B").unwrap()])
                    .unwrap()
                    .symbol_user_data(),
            )
            .with_print_function(tensor_print)
            .build()
            .unwrap();
        let value = function!(symbol, 7, mink!(4, symbol!("native_mathml_print::mu")));
        let original = value.as_view().get_data().to_vec();
        let head = "<mrow><mi>A</mi><mo>′</mo></mrow>";
        // The native printer subscripts the vector label; scoped and static
        // notation give labels and representation slots their own columns.
        for (html, prefix) in [
            (
                tensor(value.as_view(), &options(), &PrintState::new()).unwrap(),
                "<mmultiscripts><msub>",
            ),
            (scoped(&value, &options()), "<mmultiscripts><mrow>"),
            (
                expression(value.as_view(), &options(), &PrintState::new(), &|_, _| {
                    None
                })
                .unwrap(),
                "<mmultiscripts><mrow>",
            ),
        ] {
            assert!(html.starts_with(&format!("{prefix}{head}")), "{html}");
            assert_eq!(html.matches(head).count(), 1, "{html}");
            assert_eq!(html.matches("<mi>A</mi>").count(), 1, "{html}");
            assert_eq!(html.matches("<mo>′</mo>").count(), 1, "{html}");
            assert_eq!(html.matches("<mn>7</mn>").count(), 1, "{html}");
            assert_eq!(html.matches("<mi>mu</mi>").count(), 1, "{html}");
            assert!(!html.contains("<mi>B</mi>"), "{html}");
        }
        assert_eq!(value.as_view().get_data(), original);
    }

    #[test]
    fn structured_tensor_heads_use_existing_mathml_bounds() {
        let display = IndexDisplay::Sequence(vec![IndexDisplay::Symbol("<".repeat(128)); 10]);
        let symbol =
            SymbolBuilder::new(wrap_symbol!("native_mathml_print::bounded_structured_head"))
                .with_tags([SPENSO_TAG.tensor.as_str()])
                .with_user_data(display.symbol_user_data())
                .build()
                .unwrap();
        assert_eq!(IndexDisplay::from_symbol(symbol), Some(display));
        assert!(head(symbol, None).is_none());
        assert!(
            render(
                Atom::var(symbol).as_view(),
                &options(),
                &PrintState::new(),
                &|_, _| None,
            )
            .is_none()
        );
    }

    #[test]
    fn dot_compaction_defers_to_callable_owned_vector_notation() {
        let symbol = SymbolBuilder::new(wrap_symbol!("native_mathml_print::callable_vector"))
            .with_tags([
                SPENSO_TAG.tensor.as_str(),
                SPENSO_TAG.rank1.as_str(),
                TENSOR_PRINT_CALLBACK_TAG,
            ])
            .with_print_function(|_, options, _| {
                options
                    .mode
                    .is_mathml()
                    .then(|| "<mi>custom-vector</mi>".to_owned())
            })
            .build()
            .unwrap();
        let p = vector_symbol!("native_mathml_print::ordinary_vector");
        let value = function!(
            SPENSO_TAG.dot,
            function!(symbol, 2, mink!(4)),
            function!(p, 3, mink!(4))
        );
        let options = options();
        assert!(dot(value.as_view(), &options, &PrintState::new()).is_none());
        let html = value.printer(options).to_string();
        assert!(html.contains("<mi>custom-vector</mi>"), "{html}");
    }

    #[test]
    fn native_metrics_keep_original_indices_variance_and_source_bytes() {
        let rep = LibraryRep::new_dual("native_mathml_print::metric_space").unwrap();
        let mu = symbol!("native_mathml_print::metric_mu");
        let nu = symbol!("native_mathml_print::metric_nu");
        let first = function!(rep.symbol(), 3, mu);
        let second = function!(rep.symbol(), 3, nu);
        let options = options();
        for (left, right, first_lower, second_lower) in [
            (first.clone(), second.clone(), false, false),
            (dind!(first.clone()), second.clone(), true, false),
            (dind!(first.clone()), dind!(second.clone()), true, true),
        ] {
            let value = ETS.metric(left, right);
            let bytes = value.as_view().get_data().to_vec();
            let html = value.printer(options.clone()).to_string();
            assert!(html.contains("<mi>g</mi>"), "{html}");
            assert!(html.contains("<mi>metric_mu</mi>"), "{html}");
            assert!(html.contains("<mi>metric_nu</mi>"), "{html}");
            let [first, second] =
                [("metric_mu", first_lower), ("metric_nu", second_lower)].map(|(label, lower)| {
                    let index = format!("<mrow><mi>{label}</mi></mrow>");
                    if lower {
                        format!("{index}<none/>")
                    } else {
                        format!("<none/>{index}")
                    }
                });
            // The metric is symmetric: either slot order is canonical.
            assert!(
                html.contains(&format!("<mi>g</mi>{first}{second}</mmultiscripts>"))
                    || html.contains(&format!("<mi>g</mi>{second}{first}</mmultiscripts>")),
                "{html}"
            );
            assert!(!html.contains("metric_space"), "{html}");
            assert!(!html.contains("data-symbolica-path"), "{html}");
            assert_eq!(value.as_view().get_data(), bytes);
        }

        let value = ETS.metric(mink!(4, Atom::var(mu)), mink!(4, Atom::var(nu)));
        let mut state = PrintState::new();
        let html = tensor(value.as_view(), &options, &state).unwrap();
        state.in_exp_base = true;
        assert_eq!(tensor(value.as_view(), &options, &state), Some(group(html)));
    }

    #[test]
    fn native_metric_dispatch_defers_unsupported_calls() {
        let p = vector_symbol!("native_mathml_print::metric_p");
        let q = vector_symbol!("native_mathml_print::metric_q");
        let options = options();

        let mu = symbol!("native_mathml_print::unsupported_mu");
        let nu = symbol!("native_mathml_print::unsupported_nu");
        for value in [
            function!(ETS.metric, mink!(4, Atom::var(mu))),
            ETS.metric(mink!(4, Atom::var(mu)), mink!(5, Atom::var(nu))),
            ETS.metric(mink!(4, Atom::var(mu)), lor!(4, Atom::var(nu))),
            ETS.metric(function!(p, 2, mink!(4)), function!(q, 3, mink!(4))),
            ETS.metric(function!(p, 2, mink!(4)), function!(q, 3, mink!(5))),
        ] {
            assert!(tensor(value.as_view(), &options, &PrintState::new()).is_none());
            let html = value.printer(options.clone()).to_string();
            assert!(html.contains("<mi>g</mi>"), "{html}");
        }
        let ordinary_g = symbol!("native_mathml_print::g");
        let value = function!(ordinary_g, mink!(4, Atom::var(mu)), mink!(4, Atom::var(nu)));
        assert!(tensor(value.as_view(), &options, &PrintState::new()).is_none());
    }

    #[test]
    fn standalone_representations_keep_family_dimension_and_duality() {
        let rep = LibraryRep::new_dual("native_mathml_print::standalone_space").unwrap();
        let options = options();
        let dimension = symbol!("native_mathml_print::D");
        let index = symbol!("native_mathml_print::standalone_i");
        for value in [
            function!(rep.symbol(), dimension),
            function!(rep.symbol(), 3, index),
            dind!(function!(rep.symbol(), 3, index)),
        ] {
            let bytes = value.as_view().get_data().to_vec();
            let html = value.printer(options.clone()).to_string();
            assert!(html.contains("<mi>standalone_space</mi>"), "{html}");
            assert!(html.contains("<msub>"), "{html}");
            assert!(!html.contains("data-symbolica-path"), "{html}");
            if value == function!(rep.symbol(), dimension) {
                assert!(html.contains("<mi>D</mi>"), "{html}");
                assert!(!html.contains("<mmultiscripts>"), "{html}");
            } else {
                assert!(html.contains("<mn>3</mn>"), "{html}");
                assert!(html.contains("<mi>standalone_i</mi>"), "{html}");
            }
            assert_eq!(value.as_view().get_data(), bytes);
        }
        let dual = dind!(function!(rep.symbol(), 3, index));
        let html = representation(dual.as_view(), &options, &PrintState::new()).unwrap();
        assert!(html.contains("<mo>*</mo>"), "{html}");
        assert!(html.contains("</mrow><none/>"), "{html}");
        let mut power = PrintState::new();
        power.in_exp_base = true;
        assert_eq!(
            representation(dual.as_view(), &options, &power),
            Some(group(html))
        );

        let coad = LibraryRep::new_self_dual("coad").unwrap();
        for (value, family, dim) in [
            (
                function!(LibraryRep::from(Minkowski {}).symbol(), dimension),
                "mink",
                "<mi>D</mi>",
            ),
            (function!(coad.symbol(), 8), "coad", "<mn>8</mn>"),
        ] {
            let html = value.printer(options.clone()).to_string();
            assert!(html.contains(&format!("<mi>{family}</mi>")), "{html}");
            assert!(html.contains(dim), "{html}");
            assert!(!html.contains('*'), "{html}");
        }

        let value = function!(rep.symbol(), 3);
        assert!(representation(value.as_view(), &options, &PrintState::new()).is_some());
        for value in [
            Atom::var(rep.symbol()),
            function!(rep.symbol()),
            function!(rep.symbol(), 3, index, index),
            function!(AIND_SYMBOLS.dind, symbol!("native_mathml_print::opaque")),
        ] {
            assert!(representation(value.as_view(), &options, &PrintState::new()).is_none());
        }
    }
}
