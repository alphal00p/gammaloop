//! Local tensor notation for ordinary Symbolica expressions.
//!
//! These callbacks only replace bounded semantic leaves. Symbolica retains
//! ownership of arithmetic, responsive layouts, paging, and unsupported calls.
use super::*;
use symbolica::printer::{AnsiHtmlFormatter, mathml::MathMlCursorCache};

const MAX_BYTES: usize = 4096;

fn fits(atom: AtomView<'_>, remaining: &mut usize, depth: usize) -> bool {
    if *remaining == 0 || depth > 12 || atom.get_byte_size() > 2048 {
        return false;
    }
    *remaining -= 1;
    atom.children()
        .all(|child| fits(child, remaining, depth + 1))
}

fn source(atom: AtomView<'_>, options: &PrintOptions) -> Option<String> {
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
    let page = atom.render_mathml_layout_cached(
        options,
        0,
        500,
        0,
        &[],
        &MathMlCursorCache::new(atom),
        None,
    );
    if page.is_partial() || page.html.len() > MAX_BYTES {
        return None;
    }
    let mut fragment = page.fragment().into_owned();
    // This complete label is atomic in its caller. Its standalone paths would
    // otherwise collide with the surrounding expression's expansion targets.
    const PATH_PREFIX: &str = "<mrow data-symbolica-path=\"";
    while let Some(start) = fragment.find(PATH_PREFIX) {
        let value = start + PATH_PREFIX.len();
        let end = value + fragment[value..].find('"')? + 1;
        fragment.replace_range(start + "<mrow".len()..end, "");
    }
    Some(fragment)
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

pub(super) fn tensor(
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{dind, lor, mink};
    use symbolica::{function, printer::PrintMode};

    fn options() -> PrintOptions {
        PrintOptions {
            mode: PrintMode::MathML,
            hide_all_namespaces: true,
            ..PrintOptions::new()
        }
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
}
