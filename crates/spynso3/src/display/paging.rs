//! Tensor notation leaves in Symbolica's responsive MathML renderer.
use super::*;
use std::{collections::VecDeque, ffi::CString, time::Instant};
use symbolica::{
    api::python::{MathMLMarker, MathMLPage, MathMLPageRequest, MathMLPageSource, MathMLViewer},
    printer::mathml::{LAYOUTS, MathMlCursorCache},
};

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

// Retain bounded native cursor bookmarks over the original immutable owner.
self_cell::self_cell!(
    struct CachedTensor {
        owner: Arc<SymbolicTensor<PartialStructure>>,
        #[covariant]
        dependent: MathMlCursorCache,
    }
);

struct PageSource {
    tensor: CachedTensor,
    leaves: TensorLeaves,
}
struct TensorLeaves {
    settings: DisplaySettings,
    notation: Option<String>,
    aliases: IndexAliases,
    // Responsive alternatives share compiled leaves, with bounded retention.
    cache: VecDeque<((Atom, bool), Option<String>)>,
}

impl TensorLeaves {
    fn render(
        &mut self,
        py: Python<'_>,
        value: &SymbolicTensor<PartialStructure>,
        view: AtomView<'_>,
        state: &PrintState,
        deadline: Instant,
    ) -> PyResult<Option<String>> {
        let Some(symbol) = view.get_symbol() else {
            return Ok(None);
        };
        if !fits(view, &mut { MAX_LEAF_NODES }, 0) {
            return Ok(None);
        }
        // Arithmetic and ordinary calls belong to Symbolica. Only semantic
        // tensor leaves, registered labels and custom notation enter Typst.
        let semantic = self.notation.is_some()
            || symbol.has_tag(&SPENSO_TAG.tensor)
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
        let labelled =
            matches!(view, AtomView::Var(_)) && !labels::typst_heads(&view.to_owned()).is_empty();
        if !(semantic || labelled) {
            return Ok(None);
        }
        let root = value.expression().as_view().get_data();
        let key = (view.to_owned(), state.in_exp_base);
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
        let atom =
            if view.get_data().as_ptr() == root.as_ptr() && view.get_byte_size() == root.len() {
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
        let mut settings = typst_settings_source(&self.settings, &Atom::zero());
        settings.push_str(&format!(
            "\n#let settings = (..settings, index-aliases: {})\n",
            aliases.typst_source(&self.settings.index_style)
        ));
        settings.push_str(&typst_custom_print_source(&self.settings, &atom)?);
        let attachments = portable_attachments(&atom).map_err(PyRuntimeError::new_err)?;
        let tree = symbolica_typst_plugin::payload::atom_render_tree_value(&atom, &attachments)
            .map_err(|error| PyRuntimeError::new_err(error.to_string()))?;
        let mut bytes = Vec::new();
        ciborium::into_writer(&tree, &mut bytes)
            .map_err(|error| PyRuntimeError::new_err(error.to_string()))?;
        let payload = PyDict::new(py);
        payload.set_item("tree", PyBytes::new(py, &bytes))?;
        payload.set_item("settings", settings)?;
        payload.set_item("notation", self.notation.as_deref().unwrap_or(NOTATION_TYP))?;
        payload.set_item("power_base", state.in_exp_base)?;
        let remaining = deadline
            .saturating_duration_since(Instant::now())
            .as_secs_f64();
        if remaining == 0.0 {
            return Err(PyRuntimeError::new_err(
                "Tensor notation exceeded the 10-second deadline",
            ));
        }
        payload.set_item("timeout", remaining)?;
        let result: Option<String> = module(py)?
            .getattr("compile_leaf")?
            .call1((payload,))?
            .extract()?;
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

impl MathMLPageSource for PageSource {
    fn page(&mut self, py: Python<'_>, request: MathMLPageRequest<'_>) -> PyResult<MathMLPage> {
        if !(1..=500).contains(&request.size)
            || request.layout >= LAYOUTS.len()
            || request.path.len() > 64
            || request
                .expansions
                .iter()
                .any(|(path, size)| path.len() > 64 || !(1..=500).contains(size))
        {
            return Err(PyValueError::new_err("invalid tensor page request"));
        }
        let cursors = self.tensor.borrow_dependent();
        let atom = cursors
            .locate(request.path)
            .ok_or_else(|| PyValueError::new_err("unknown subexpression"))?;
        let deadline = Instant::now() + std::time::Duration::from_secs(10);
        let mut error = None;
        let mut leaf = |atom: AtomView<'_>, _options: &PrintOptions, state: &PrintState| {
            if error.is_some() {
                return None;
            }
            match self
                .leaves
                .render(py, self.tensor.borrow_owner(), atom, state, deadline)
            {
                Ok(fragment) => fragment,
                Err(failure) => {
                    error = Some(failure);
                    None
                }
            }
        };
        let rendered = atom.render_mathml_layout_cached(
            &PrintOptions::new()
                .multiplication_operator('·')
                .num_exp_as_superscript(true),
            request.start,
            request.size,
            request.layout,
            request.expansions,
            cursors,
            Some(&mut leaf),
        );
        let mut page = MathMLPage::new(
            format!(
                "<style>{NOTEBOOK_STYLE}</style>{}",
                rendered
                    .html
                    .replacen("<math ", "<math data-spenso-math=\"\" ", 1)
            ),
            rendered.start,
            rendered.end,
            rendered.total,
        );
        if let Some(error) = error {
            page.html = format!(
                "<p role=\"alert\">Page rendering failed: {}</p>",
                AnsiHtmlFormatter::escape_html(
                    &error.to_string().chars().take(2048).collect::<String>()
                )
            );
            return Ok(page);
        }
        let (kind, unit) = match atom {
            AtomView::Fun(_) if rendered.custom => ("Function", "expressions"),
            AtomView::Add(_) => ("Sum", "terms"),
            AtomView::Mul(_) => ("Product", "factors"),
            AtomView::Fun(_) => ("Function", "arguments"),
            AtomView::Pow(_) => ("Power", "expressions"),
            _ => ("Value", "expressions"),
        };
        page.kind = kind.into();
        page.unit = unit.into();
        page.markers = rendered
            .markers
            .into_iter()
            .map(|marker| {
                let mut path = request.path.to_vec();
                path.extend(marker.path);
                MathMLMarker::new(path, marker.label, marker.start)
            })
            .collect();
        Ok(page)
    }
}

fn module(py: Python<'_>) -> PyResult<Bound<'_, PyModule>> {
    let modules = py.import("sys")?.getattr("modules")?;
    let name = "_spenso_paging";
    if let Ok(module) = modules.get_item(name) {
        return Ok(module.cast_into::<PyModule>()?);
    }
    let code = CString::new(include_str!("../../typst/paging.py")).unwrap();
    let module = PyModule::from_code(py, &code, c"spenso_paging.py", c"_spenso_paging")?;
    module.add_function(wrap_pyfunction!(compile_typst, &module)?)?;
    modules.set_item(name, &module)?;
    Ok(module)
}

pub(crate) fn viewer(
    py: Python<'_>,
    value: Arc<SymbolicTensor<PartialStructure>>,
    settings: DisplaySettings,
    notation: Option<String>,
    size: usize,
) -> PyResult<Py<PyAny>> {
    let source = PageSource {
        tensor: CachedTensor::new(value, |value| {
            MathMlCursorCache::new(value.expression().as_view())
        }),
        leaves: TensorLeaves {
            settings,
            notation,
            aliases: IndexAliases::default(),
            cache: VecDeque::new(),
        },
    };
    MathMLViewer::new(source)
        .page_size(size)
        .layout_count(LAYOUTS.len())
        .supports_expansion(true)
        .into_python(py)
}

pub(crate) fn static_html(
    py: Python<'_>,
    value: Arc<SymbolicTensor<PartialStructure>>,
    settings: DisplaySettings,
    notation: Option<String>,
) -> PyResult<String> {
    viewer(py, value, settings, notation, DEFAULT_PAGE_SIZE)?
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
