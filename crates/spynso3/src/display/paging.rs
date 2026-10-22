//! Tensor notation leaves in Symbolica's responsive MathML renderer.
use super::*;
use std::collections::VecDeque;
use symbolica::api::python::MathMLViewer;
use symbolica::printer::{PrintMode, PrintUserData};

pub(crate) const DEFAULT_PAGE_SIZE: usize = 25;
const MAX_LEAF_NODES: usize = 256;
const MAX_LEAF_BYTES: usize = 4096;
const MAX_ALIASES: usize = 4096;

fn fits(atom: AtomView<'_>, budget: &mut usize, depth: usize) -> bool {
    if *budget == 0
        || depth > 32
        || atom.get_byte_size() > MAX_LEAF_BYTES
        || atom
            .get_symbol()
            .is_some_and(|symbol| symbol.get_name().len() > MAX_LEAF_BYTES)
    {
        return false;
    }
    *budget -= 1;
    atom.children().all(|child| fits(child, budget, depth + 1))
}

pub(crate) fn is_large(value: &SymbolicTensor<PartialStructure>) -> bool {
    !fits(value.expression().as_view(), &mut 1000, 0)
}

type LeafCacheKey = (Atom, [u16; 10]);

impl DisplaySettings {
    /// Tensor presentation options shared by static exports and scoped paging.
    pub(crate) fn mathml_options(&self) -> PrintOptions {
        let mut options = display_options_with_settings(TensorDisplayMode::Typst, self);
        options.mode = PrintMode::MathML;
        for (key, value) in [
            ("spenso.tensor_layout", &self.tensor_layout),
            ("spenso.index_gap", &self.index_gap),
            ("spenso.factor_gap", &self.factor_gap),
        ] {
            options
                .custom_print_mode
                .insert(key.into(), PrintUserData::String(value.clone()));
        }
        options
    }
}

/// Complete export uses the same tensor notation and source-label owner as
/// paging. Its explicit export cap reports an error rather than dropping terms.
pub(crate) fn atom_to_html(
    py: Python<'_>,
    atom: &Atom,
    settings: &DisplaySettings,
) -> PyResult<String> {
    let aliases = IndexAliases::for_atom(atom, &settings.index_style);
    html(py, atom, settings, &aliases)
}

pub(crate) fn structured_to_html(
    py: Python<'_>,
    value: &SymbolicTensor<PartialStructure>,
    settings: &DisplaySettings,
) -> PyResult<String> {
    let aliases = IndexAliases::for_descriptor(value, &settings.index_style);
    html(py, &value.presentation_atom(), settings, &aliases)
}

fn html(
    py: Python<'_>,
    atom: &Atom,
    settings: &DisplaySettings,
    aliases: &IndexAliases,
) -> PyResult<String> {
    let options = settings.mathml_options();
    let invalid_label = std::cell::Cell::new(false);
    let body = spenso::network::tags::mathml::expression(
        atom.as_view(),
        &options,
        &PrintState::new(),
        &|representation, index| {
            let key = (representation, index.to_owned());
            if let Some((_, display)) = aliases.entries.get(&key) {
                let rendered = display.to_mathml();
                invalid_label.set(invalid_label.get() || rendered.is_none());
                return rendered;
            }
            aliases.source_indices.get(&key).and_then(|source| {
                let rendered = spenso::network::tags::mathml::expression(
                    source.as_view(),
                    &options,
                    &PrintState::new(),
                    &|_, _| None,
                );
                invalid_label.set(invalid_label.get() || rendered.is_none());
                rendered
            })
        },
    );
    if let Some(error) = PyErr::take(py) {
        return Err(error);
    }
    let body = body.filter(|_| !invalid_label.get()).ok_or_else(|| {
        PyValueError::new_err(
            "complete MathML export exceeds its resource limit or contains an unsupported label; use paged() to inspect it",
        )
    })?;
    Ok(notebook_html_fragment(&format!(
        "<math xmlns=\"http://www.w3.org/1998/Math/MathML\"><mrow>{body}</mrow></math>"
    )))
}

struct TensorLeaves {
    settings: DisplaySettings,
    aliases: IndexAliases,
    // Responsive alternatives share direct fragments, with bounded retention.
    cache: VecDeque<(LeafCacheKey, Option<String>)>,
}

impl TensorLeaves {
    fn render(
        &mut self,
        py: Python<'_>,
        value: &SymbolicTensor<PartialStructure>,
        view: AtomView<'_>,
        state: &PrintState,
    ) -> PyResult<Option<String>> {
        let Some(symbol) = view.get_symbol() else {
            return Ok(None);
        };
        if !fits(view, &mut { MAX_LEAF_NODES }, 0) {
            return Ok(None);
        }
        // Arithmetic and ordinary calls belong to Symbolica. Only semantic
        // tensor leaves enter the scoped semantic printer.
        let semantic = symbol.has_tag(&SPENSO_TAG.tensor)
            || symbol.has_tag("tensor")
            || symbol.has_tag(&SPENSO_TAG.index)
            || symbol.has_tag(&SPENSO_TAG.representation)
            || symbol.has_tag(TENSOR_PRINT_CALLBACK_TAG)
            || symbol
                .get_tags()
                .iter()
                .any(|tag| tag.starts_with(TENSOR_PRINT_HEAD_PREFIX))
            || matches!(
                symbol.get_name(),
                "spenso::bracket"
                    | "spenso::cas"
                    | "spenso::idx"
                    | "spenso::gram"
                    | "spenso::gamma"
                    | "spenso::gamma0"
                    | "spenso::gamma5"
                    | "spenso::projp"
                    | "spenso::projm"
                    | "spenso::dot"
                    | "spenso::chain"
                    | "spenso::trace"
                    | "spenso::Nc"
            );
        if !semantic {
            return Ok(None);
        }
        let root = value.expression().as_view().get_data();
        // Custom printers receive the complete state. Reusing a fragment must
        // preserve every observable flag and level, including exponent context.
        let key = (
            view.to_owned(),
            [
                u16::from(state.in_sum),
                u16::from(state.in_product),
                u16::from(state.suppress_one),
                u16::from(state.in_exp),
                u16::from(state.in_exp_base),
                u16::from(state.top_level_add_child),
                u16::from(state.superscript),
                state.level,
                state.bracket_level,
                state.indentation_level,
            ],
        );
        if let Some(index) = self
            .cache
            .iter()
            .position(|(candidate, _)| *candidate == key)
        {
            let entry = self.cache.remove(index).unwrap();
            let result = entry.1.clone();
            self.cache.push_back(entry);
            return Ok(result);
        }
        // The inherited Python expression may own a byte-for-byte Atom clone.
        // Equality is bounded by fits above and does not scan a large root.
        let atom = if view.get_data() == root {
            value.presentation_atom()
        } else {
            value.presentation_fragment(view)
        };
        let mut aliases = IndexAliases::for_atom(&atom, &self.settings.index_style);
        // Reserve the public free ports without walking the original expression.
        let mut bundle = FunctionBuilder::new(spenso::structure::abstract_index::AIND_SYMBOLS.aind)
            .add_arg(&atom);
        for slot in value.structure().logical_slots() {
            bundle = bundle.add_arg(composition::port_atom(slot));
        }
        let discovered = IndexAliases::for_atom(&bundle.finish(), &self.settings.index_style);
        for (key, (mut position, mut display)) in discovered.entries {
            if let Some(existing) = self.aliases.entries.get(&key) {
                if discovered
                    .reserved_labels
                    .get(&key.0)
                    .is_some_and(|labels| labels.contains(&IndexAliases::label_key(&existing.1)))
                {
                    // An explicit internal label may first appear on a later
                    // page. Keep it distinct by retaining this index's source
                    // spelling here, without renaming already displayed parts.
                    aliases.entries.remove(&key);
                    continue;
                }
                if let Some(alias) = aliases.entries.get_mut(&key) {
                    *alias = existing.clone();
                }
                continue;
            }
            let mut retain = self.aliases.entries.len() < MAX_ALIASES;
            if self.settings.index_style == "alphabet" {
                let metadata = RepresentationMetadata::from_symbol(key.0);
                while retain {
                    let label = IndexAliases::label_key(&display);
                    if !self
                        .aliases
                        .reserved_labels
                        .get(&key.0)
                        .is_some_and(|used| used.contains(&label))
                        && !discovered
                            .reserved_labels
                            .get(&key.0)
                            .is_some_and(|used| used.contains(&label))
                    {
                        break;
                    }
                    let Some(next_position) = position.and_then(|p| p.checked_add(1)) else {
                        retain = false;
                        break;
                    };
                    let Some(next) = metadata
                        .as_ref()
                        .and_then(|metadata| metadata.index_palette.resolve(next_position))
                    else {
                        retain = false;
                        break;
                    };
                    position = Some(next_position);
                    display = next;
                }
            }
            if retain {
                self.aliases
                    .reserved_labels
                    .entry(key.0)
                    .or_default()
                    .insert(IndexAliases::label_key(&display));
                let alias = (position, display);
                if let Some(local) = aliases.entries.get_mut(&key) {
                    *local = alias.clone();
                }
                self.aliases.entries.insert(key, alias);
            } else {
                // Keep established names stable. New keys retain their source
                // spelling once the bounded viewer alias inventory is full.
                aliases.entries.remove(&key);
            }
        }
        let options = self.settings.mathml_options();
        // Aliases remain borrowed viewer state, rather than serialized options
        // or symbols inserted into the expression being paged.
        let result = spenso::network::tags::mathml::render(
            atom.as_view(),
            &options,
            state,
            &|representation, index| {
                let key = (representation, index.to_owned());
                aliases
                    .entries
                    .get(&key)
                    .and_then(|(_, display)| display.to_mathml())
                    .or_else(|| {
                        let source = aliases.source_indices.get(&key)?;
                        if !fits(source.as_view(), &mut { MAX_LEAF_NODES }, 0) {
                            return None;
                        }
                        Some(source.printer(options.clone()).to_string())
                    })
            },
        );
        if let Some(error) = PyErr::take(py) {
            return Err(error);
        }
        // Oversized calls defer to structural rendering so their arguments can
        // still be paged and inspected, rather than becoming an atomic loop.
        let result = result.filter(|html| html.len() <= MAX_LEAF_BYTES);
        self.cache.push_back((key, result.clone()));
        if self.cache.len() > 128 {
            self.cache.pop_front();
        }
        Ok(result)
    }
}

pub(crate) fn viewer(
    py: Python<'_>,
    expression: Py<PythonExpression>,
    value: Arc<SymbolicTensor<PartialStructure>>,
    settings: DisplaySettings,
    notation: Option<String>,
    size: usize,
) -> PyResult<Py<PyAny>> {
    if notation.is_some() {
        return Err(PyValueError::new_err(
            "The direct MathML pager does not accept Typst notation_source; use to_typst_project() or to_svg() for custom Typst exports",
        ));
    }
    let mut leaves = TensorLeaves {
        settings,
        aliases: IndexAliases::default(),
        cache: VecDeque::new(),
    };
    MathMLViewer::expression(
        expression,
        PrintOptions::new()
            .multiplication_operator('·')
            .num_exp_as_superscript(true),
        Some(Box::new(move |py, atom, _options, state| {
            leaves.render(py, &value, atom, state)
        })),
    )
    .page_size(size)
    .into_python(py)
}

pub(crate) fn static_html(
    py: Python<'_>,
    expression: Py<PythonExpression>,
    value: Arc<SymbolicTensor<PartialStructure>>,
    settings: DisplaySettings,
    notation: Option<String>,
) -> PyResult<String> {
    viewer(py, expression, value, settings, notation, DEFAULT_PAGE_SIZE)?
        .call_method0(py, "_repr_html_")?
        .extract(py)
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn semantic_leaf_work_is_bounded_before_notation_callbacks() {
        let small = symbolica::parse!("f(x,y)");
        assert!(fits(small.as_view(), &mut { MAX_LEAF_NODES }, 0));
        assert!(!fits(small.as_view(), &mut 1, 0));
        assert!(!fits(small.as_view(), &mut { MAX_LEAF_NODES }, 33));
    }
}
