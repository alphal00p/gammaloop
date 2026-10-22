//! Local tensor notation for ordinary Symbolica expressions.
//!
//! These callbacks only replace bounded semantic leaves. Symbolica retains
//! ownership of arithmetic, responsive layouts, paging, and unsupported calls.
use super::*;
use crate::structure::slot::{SlotMatch, SlotMatcher};
use symbolica::printer::{AnsiHtmlFormatter, PrintMode, PrintUserData};

const MAX_BYTES: usize = 4096;

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
    // Reversible index symbols retain their original local label as Atom data.
    // Display simple payloads without renaming the index or inventing aliases.
    // Nested payloads remain opaque, avoiding recursive metadata expansion.
    let atom = if let AtomView::Var(variable) = atom
        && variable.get_symbol().has_tag(&SPENSO_TAG.index)
        && let symbolica::atom::UserData::Atom(payload) = variable.get_symbol().get_data()
        && let AtomView::Fun(label) = payload.as_view()
        && label.get_symbol().has_tag(&SPENSO_TAG.index)
        && fits(payload.as_view(), &mut 64, 0)
        && label.iter().all(|argument| match argument {
            AtomView::Num(_) => true,
            AtomView::Var(variable) => matches!(
                variable.get_symbol().get_data(),
                symbolica::atom::UserData::None
            ),
            _ => false,
        }) {
        payload.as_view()
    } else {
        atom
    };
    let native = atom
        .printer(PrintOptions {
            mode: PrintMode::MathML,
            max_terms: Some(500),
            ..options.clone()
        })
        .to_string();
    fragment(native)
}

/// Prepare a complete, path-free MathML fragment for embedding in a semantic leaf.
///
/// Native and custom printer output must fit within 4096 bytes and contain no
/// `math` root or omission placeholder. Remove actual `data-symbolica-*`
/// attributes while preserving token text and other attributes, including local
/// visual slash and open-axis annotations. Malformed opening attributes decline.
pub fn fragment(mut value: String) -> Option<String> {
    if value.len() > MAX_BYTES || value.contains("<mtext>⋯</mtext>") || value.contains("<math") {
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

fn head(symbol: Symbol) -> Option<String> {
    if let Some(source) = symbol.get_tags().iter().find_map(|tag| {
        tag.strip_prefix(TENSOR_PRINT_HEAD_PREFIX)?
            .strip_prefix("mathml:")
    }) {
        return (source.len() + "<mrow></mrow>".len() <= MAX_BYTES)
            .then(|| format!("<mrow>{source}</mrow>"));
    }
    let label = symbol
        .get_tags()
        .iter()
        .find_map(|tag| tag.strip_prefix("spenso::tensor-label:"))
        .unwrap_or_else(|| symbol.get_stripped_name());
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
    let head = head(symbol)?;
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
        return symbol.has_tag(&SPENSO_TAG.tensor).then(|| head(symbol))?;
    }
    let AtomView::Fun(function) = atom else {
        return None;
    };
    let values = function.iter().collect::<Vec<_>>();
    if symbol == crate::network::library::symbolic::ETS.metric {
        let [first, second] = values.as_slice() else {
            return None;
        };
        let first = tensor_slot(*first)?;
        let second = tensor_slot(*second)?;
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
            && natural_index(values[1]) == Some(1)
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
            is_chain_marker(*value, SPENSO_TAG.chain_in)
                || is_chain_marker(*value, SPENSO_TAG.chain_out)
        })
    {
        return None;
    }
    let mut labels = Vec::new();
    let mut columns = Vec::new();
    let open_count = values
        .iter()
        .filter(|value| {
            representation_value(**value).is_some()
                || tensor_slot(**value).is_some_and(|slot| {
                    matches!(slot.index, AtomView::Fun(index)
                        if index.get_symbol() == AIND_SYMBOLS.openind && index.get_nargs() == 2)
                })
        })
        .count();
    let mut axis = 0;
    for value in &values {
        let slot = tensor_slot(*value);
        let representation = representation_value(*value);
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
                || component.iter().any(|index| natural_index(index).is_none())
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
    let (symbol, dimension, index, class, row) = if let Some(slot) = tensor_slot(atom) {
        (
            slot.representation,
            slot.dimension,
            Some(slot.index),
            RepresentationMetadata::from_symbol(slot.representation)?.class,
            slot.row,
        )
    } else {
        let value = representation_value(atom)?;
        (value.symbol, value.dimension, None, value.class, value.row)
    };
    // A bare symbol defers to Symbolica, retaining its namespace/style options
    // without calling this representation printer recursively on the function.
    let mut body = format!(
        "<msub>{}<mrow>{}</mrow></msub>",
        source(Atom::var(symbol).as_view(), options)?,
        source(dimension, options)?,
    );
    let dual = matches!(unwrap_tensor_display(atom), AtomView::Fun(wrapper)
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
    let AtomView::Fun(function) = unwrap_tensor_display(atom) else {
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
        if tensor_slot(value).is_some() {
            return None;
        }
        if let Some(value) = representation_value(value) {
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
pub fn render(
    atom: AtomView<'_>,
    options: &PrintOptions,
    state: &PrintState,
    index: &dyn for<'a> Fn(Symbol, AtomView<'a>) -> Option<String>,
) -> Option<String> {
    if !fits(atom, &mut 256, 0) {
        return None;
    }
    let context = Render::new(options, index)?;
    let atom = unwrap_tensor_display(atom);
    fragment(context.leaf(atom, state).or_else(|| {
        let symbol = atom.get_symbol()?;
        if context.owns(atom) {
            return None;
        }
        (symbol.get_print_function()?)(atom, options, state)
    })?)
}

#[derive(Clone, Copy, PartialEq, Eq)]
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
}

impl<'a> Render<'a> {
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
    ) -> Option<Self> {
        let text = |key: &str| match options.custom_print_mode.get(key) {
            Some(PrintUserData::String(value)) => Some(value.as_str()),
            _ => None,
        };
        let layout = match text("spenso.tensor_layout").unwrap_or("ports") {
            "ports" => Layout::Ports,
            "schoonschip" => Layout::Schoonschip,
            "call" => Layout::Call,
            _ => return None,
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
                .find_map(|unit| value.strip_suffix(unit).map(|number| (number, unit)))?;
            let number = number.parse::<f64>().ok()?;
            number.is_finite().then(|| format!("{number}{unit}"))
        };
        Some(Self {
            options,
            index,
            settings,
            layout,
            index_gap: gap("spenso.index_gap", "0.08em")?,
            factor_gap: gap("spenso.factor_gap", "0.12em")?,
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
        {
            return fragment(head(symbol)?);
        }
        Some(match symbol.get_name() {
            "spenso::gamma" => "<mi>γ</mi>".to_owned(),
            "spenso::gamma0" => "<msub><mi>γ</mi><mn>0</mn></msub>".to_owned(),
            "spenso::gamma5" => "<msub><mi>γ</mi><mn>5</mn></msub>".to_owned(),
            "spenso::projp" => "<msup><mi>ℙ</mi><mo>+</mo></msup>".to_owned(),
            "spenso::projm" => "<msup><mi>ℙ</mi><mo>−</mo></msup>".to_owned(),
            _ => head(symbol)?,
        })
    }

    fn leaf(&self, atom: AtomView<'_>, state: &PrintState) -> Option<String> {
        let atom = unwrap_tensor_display(atom);
        let symbol = atom.get_symbol()?;
        if symbol.has_tag(TENSOR_PRINT_CALLBACK_TAG)
            && let Some(custom) = tensor_custom_print(atom, self.options, state)
        {
            if let Some(custom) = fragment(custom) {
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
        if tensor_slot(atom).is_some() || representation_value(atom).is_some() {
            return self.representation(atom, state);
        }
        if !symbol.has_tag(&SPENSO_TAG.tensor) && !symbol.has_tag(&SPENSO_TAG.index) {
            return None;
        }
        if let AtomView::Var(_) = atom {
            if symbol.has_tag(&SPENSO_TAG.index) {
                return if let Some(display) = IndexDisplay::from_symbol(symbol) {
                    fragment(display.to_mathml()?)
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
                && natural_index(arguments[1]) == Some(1)
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
        if symbol == crate::network::library::symbolic::ETS.metric {
            let [first, second] = arguments.as_slice() else {
                return None;
            };
            let mut matcher = SlotMatcher::default();
            let mut representation = |argument| {
                let argument = unwrap_tensor_display(argument);
                if tensor_slot(argument).is_none() && representation_value(argument).is_none() {
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
                .any(|index| natural_index(*index).is_none())
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
        let atom = unwrap_tensor_display(atom);
        let body = match atom {
            AtomView::Add(sum) => {
                let body = self.arguments(&sum.iter().collect::<Vec<_>>(), "<mo>+</mo>")?;
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
                let mut base_state = PrintState::new();
                base_state.in_exp_base = true;
                let mut base = self.visual(base, &base_state)?;
                if numeric_fence {
                    base = group(base);
                }
                let mut exponent_state = PrintState::new();
                exponent_state.in_exp = true;
                let exponent = self.visual(exponent, &exponent_state)?;
                let body = format!("<msup><mrow>{base}</mrow><mrow>{exponent}</mrow></msup>");
                if state.in_exp_base { group(body) } else { body }
            }
            _ => {
                if let Some(body) = self.leaf(atom, state) {
                    body
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
        fragment(body)
    }

    fn arguments(&self, arguments: &[AtomView<'_>], separator: &str) -> Option<String> {
        let mut output = String::new();
        for (position, argument) in arguments.iter().enumerate() {
            if position != 0 {
                output.push_str(separator);
            }
            output.push_str(&self.visual(*argument, &PrintState::new())?);
            if output.len() > MAX_BYTES {
                return None;
            }
        }
        Some(output)
    }

    fn factors(&self, factors: &[AtomView<'_>]) -> Option<String> {
        let mut state = PrintState::new();
        state.in_product = true;
        let mut output = String::new();
        for (position, factor) in factors.iter().enumerate() {
            if position != 0 {
                output.push_str(&Self::space(&self.factor_gap));
            }
            output.push_str(&self.visual(*factor, &state)?);
            if output.len() > MAX_BYTES {
                return None;
            }
        }
        Some(format!("<mrow>{output}</mrow>"))
    }

    fn index(&self, symbol: Symbol, index: AtomView<'_>) -> Option<String> {
        if let Some(alias) = (self.index)(symbol, index) {
            return fragment(alias);
        }
        if let Some(index) = natural_index(index)
            && let Some(display) = RepresentationMetadata::from_symbol(symbol)?
                .index_palette
                .resolve(index)
        {
            return fragment(display.to_mathml()?);
        }
        if let AtomView::Var(variable) = index
            && variable.get_symbol().has_tag(&SPENSO_TAG.index)
            && let Some(display) = IndexDisplay::from_symbol(variable.get_symbol())
        {
            return fragment(display.to_mathml()?);
        }
        self.visual(index, &PrintState::new())
    }

    fn qualified(&self, source: String, symbol: Symbol, dimension: AtomView<'_>) -> Option<String> {
        if !self.settings.with_dim {
            return Some(source);
        }
        let label = RepresentationMetadata::from_symbol(symbol)?
            .label
            .to_mathml()?;
        let dimension = self.visual(dimension, &PrintState::new())?;
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
            representation_value(argument).is_some() || tensor_slot(argument).is_some_and(|slot| {
            matches!(slot.index, AtomView::Fun(index) if index.get_symbol() == AIND_SYMBOLS.openind && index.get_nargs() == 2)
        })
        };
        let numbered = arguments.iter().filter(|argument| open(**argument)).count() > 1;
        let mut axis = 0;
        arguments
            .iter()
            .map(|argument| {
                let slot = tensor_slot(*argument);
                let rep = representation_value(*argument);
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
        let separator = if self.settings.commas {
            format!("<mo>,</mo>{gap}")
        } else {
            gap.clone()
        };
        let row = |row| {
            columns
                .iter()
                .filter(|(_, side)| *side == row)
                .map(|(source, _)| source.as_str())
                .collect::<Vec<_>>()
                .join(&separator)
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
        let AtomView::Fun(function) = unwrap_tensor_display(atom) else {
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
            if tensor_slot(argument).is_some()
                || matches!(matcher.classify(argument), SlotMatch::Explicit(_))
            {
                return None;
            }
            if let Some(candidate) = representation_value(argument) {
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
                    .map(|label| Some((self.visual(*label, &PrintState::new())?, IndexRow::Bottom)))
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
        let markers = chain_markers(arguments);
        // A malformed marker pair must not become a visible index or label.
        if markers.is_none()
            && arguments.iter().any(|argument| {
                is_chain_marker(*argument, SPENSO_TAG.chain_in)
                    || is_chain_marker(*argument, SPENSO_TAG.chain_out)
            })
        {
            return None;
        }
        let open_columns = self.open_columns(arguments)?;
        for (position, argument) in arguments.iter().copied().enumerate() {
            if markers
                .is_some_and(|markers| position == markers.input || position == markers.output)
            {
                continue;
            }
            if let Some(column) = &open_columns[position] {
                columns.push(column.clone());
            } else if let Some(slot) = tensor_slot(argument) {
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
                let source = self.visual(argument, &PrintState::new())?;
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
        fragment(base)
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
        let (symbol, dimension, index, class, row) = if let Some(slot) = tensor_slot(atom) {
            (
                slot.representation,
                slot.dimension,
                Some(slot.index),
                RepresentationMetadata::from_symbol(slot.representation)?.class,
                slot.row,
            )
        } else {
            let rep = representation_value(atom)?;
            (rep.symbol, rep.dimension, None, rep.class, rep.row)
        };
        let label = RepresentationMetadata::from_symbol(symbol)?
            .label
            .to_mathml()?;
        let dimension = self.visual(dimension, &PrintState::new())?;
        let mut body = format!("<msub><mrow>{label}</mrow><mrow>{dimension}</mrow></msub>");
        let dual = matches!(unwrap_tensor_display(atom), AtomView::Fun(wrapper) if wrapper.get_symbol() == AIND_SYMBOLS.dind)
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
        for (position, endpoint) in [*start, *end].into_iter().enumerate() {
            if is_chain_marker(endpoint, SPENSO_TAG.chain_in)
                || is_chain_marker(endpoint, SPENSO_TAG.chain_out)
            {
                continue;
            }
            if let Some((label, rep)) = self.vector(endpoint) {
                labels[position] = Some(label);
                columns.push((self.port(rep)?, rep.row));
            } else if let Some(column) = &open[position] {
                columns.push(column.clone());
            } else if let Some(slot) = tensor_slot(endpoint) {
                let source = self.index(slot.representation, slot.index)?;
                columns.push((
                    self.qualified(source, slot.representation, slot.dimension)?,
                    slot.row,
                ));
            } else {
                labels[position] = Some(self.visual(endpoint, &PrintState::new())?);
            }
        }
        body = self.attachment(body, &columns);
        body = Self::ports(
            body,
            &labels[0].iter().cloned().collect::<Vec<_>>(),
            &labels[1].iter().cloned().collect::<Vec<_>>(),
        );
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
                self.visual(*representation, &PrintState::new())?
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
    use crate::{
        dind, euc, lor, mink,
        structure::representation::{LibraryRep, Minkowski},
    };
    use symbolica::{function, printer::PrintMode};

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

    fn gamma() -> Symbol {
        register_tensor_symbol(
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
        use crate::network::library::symbolic::ETS;

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
            crate::structure::representation::IndexPalette::Numeric,
            IndexRow::Bottom,
        )
        .unwrap();
        let p = crate::vector_symbol!("scoped_mathml::p");
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
        let u = crate::vector_symbol!("scoped_mathml::u");
        let v = crate::vector_symbol!("scoped_mathml::v");
        let value = function!(
            crate::tensor_symbol!("scoped_mathml::A"),
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
        let head = crate::tensor_symbol!("scoped_mathml::Aliases");
        let original = symbol!("scoped_mathml::original_index");
        let dimension = symbol!("scoped_mathml::D");
        let value = function!(
            head,
            mink!(
                crate::structure::dimension::Dimension::Symbolic(dimension.into()),
                Atom::var(original)
            )
        );
        let mut options = options();
        let mut settings = SpensoPrintSettings::typst();
        settings.with_dim = true;
        options
            .custom_print_mode
            .extend(ahash::HashMap::from(settings));
        let calls = std::cell::Cell::new(0);
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
        let head = SymbolBuilder::new(symbolica::wrap_symbol!("scoped_mathml::component"))
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
        options
            .custom_print_mode
            .extend(ahash::HashMap::from(settings));
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
        options
            .custom_print_mode
            .extend(ahash::HashMap::from(settings));
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
        let head = crate::tensor_symbol!("scoped_mathml::Bounded");
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
        use std::sync::{
            Arc,
            atomic::{AtomicUsize, Ordering},
        };

        let calls = Arc::new(AtomicUsize::new(0));
        let callback_calls = Arc::clone(&calls);
        let head = SymbolBuilder::new(symbolica::wrap_symbol!("scoped_mathml::oversized_callable"))
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
        let p = crate::vector_symbol!("scoped_mathml::dot_p");
        let q = crate::vector_symbol!("scoped_mathml::dot_q");
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
        let head = crate::tensor_symbol!("native_mathml_print::T");
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
        let p = crate::vector_symbol!("native_mathml_print::p");
        let q = crate::vector_symbol!("native_mathml_print::q");
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
        let symbol = SymbolBuilder::new(symbolica::wrap_symbol!(
            "native_mathml_print::multi_root_head"
        ))
        .with_tags([
            SPENSO_TAG.tensor.as_str(),
            SPENSO_TAG.rank1.as_str(),
            "spenso::print-head:mathml:<mi>A</mi><mo>′</mo>",
        ])
        .with_print_function(tensor_print)
        .build()
        .unwrap();
        let value = function!(symbol, 7, mink!(4, symbol!("native_mathml_print::mu")));
        let html = tensor(value.as_view(), &options(), &PrintState::new()).unwrap();
        assert!(
            html.starts_with("<mmultiscripts><msub><mrow><mi>A</mi><mo>′</mo></mrow>"),
            "{html}"
        );
    }

    #[test]
    fn dot_compaction_defers_to_callable_owned_vector_notation() {
        let symbol = SymbolBuilder::new(symbolica::wrap_symbol!(
            "native_mathml_print::callable_vector"
        ))
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
        let p = crate::vector_symbol!("native_mathml_print::ordinary_vector");
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
        use crate::network::library::symbolic::ETS;

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
        use crate::network::library::symbolic::ETS;

        let p = crate::vector_symbol!("native_mathml_print::metric_p");
        let q = crate::vector_symbol!("native_mathml_print::metric_q");
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
