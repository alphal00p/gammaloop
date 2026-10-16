//! Bounded presentation trees. Ellipses are render-tree nodes, never algebraic atoms.
use super::*;
use ciborium::Value;
use std::{collections::BTreeMap, ffi::CString};
use symbolica::atom::representation::ListIterator;

const MAX_NODES: usize = 2_000;
const MAX_DEPTH: usize = 32;
const MAX_BYTES: usize = 64 * 1024;
pub(crate) const DEFAULT_PAGE_SIZE: usize = 25;

fn map(fields: impl IntoIterator<Item = (&'static str, Value)>) -> Value {
    Value::Map(
        fields
            .into_iter()
            .map(|(k, v)| (Value::Text(k.into()), v))
            .collect(),
    )
}
fn text(s: &str) -> Value {
    Value::Text(s.into())
}
fn children(atom: AtomView<'_>) -> Option<ListIterator<'_>> {
    match atom {
        AtomView::Add(a) => Some(a.iter()),
        AtomView::Mul(a) => Some(a.iter()),
        AtomView::Pow(a) => Some(a.iter()),
        AtomView::Fun(a) => Some(a.iter()),
        _ => None,
    }
}
fn fits(atom: AtomView<'_>, nodes: &mut usize, depth: usize) -> bool {
    if depth > MAX_DEPTH
        || *nodes == 0
        || atom.get_byte_size() > MAX_BYTES
        || atom
            .get_symbol()
            .is_some_and(|s| s.get_name().len() > MAX_BYTES)
    {
        return false;
    }
    *nodes -= 1;
    children(atom).is_none_or(|mut c| c.all(|a| fits(a, nodes, depth + 1)))
}

fn large_sum(a: AtomView<'_>, limit: usize, budget: &mut usize) -> bool {
    large_sum_at(a, limit, budget, 0)
}
fn large_sum_at(a: AtomView<'_>, limit: usize, budget: &mut usize, depth: usize) -> bool {
    if *budget == 0 || depth > MAX_DEPTH {
        return true;
    }
    *budget -= 1;
    if matches!(a, AtomView::Add(s) if s.get_nargs()>limit) {
        return true;
    }
    children(a).is_some_and(|mut c| c.any(|n| large_sum_at(n, limit, budget, depth + 1)))
}

pub(crate) fn is_large(value: &SymbolicTensor<PartialStructure>) -> bool {
    let a = value.expression().as_view();
    let mut nodes = MAX_NODES;
    !fits(a, &mut nodes, 0) || large_sum(a, 100, &mut { MAX_NODES })
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct Span {
    start: usize,
    len: usize,
}
#[derive(Clone, Copy)]
struct Cursor {
    focus: Span,
    sum: Span,
    next: Span,
    index: usize,
}

#[pyclass]
struct PageSource {
    value: Arc<SymbolicTensor<PartialStructure>>,
    settings: DisplaySettings,
    notation: Option<String>,
    spans: BTreeMap<usize, Span>,
    cursors: Vec<Cursor>,
    aliases: IndexAliases,
}

struct Builder<'a> {
    root: &'a [u8],
    value: &'a SymbolicTensor<PartialStructure>,
    settings: &'a DisplaySettings,
    aliases: &'a mut IndexAliases,
    spans: &'a mut BTreeMap<usize, Span>,
    cursor: Option<Cursor>,
    focus: Span,
    next: Option<Cursor>,
    active: Option<Span>,
    selected_page: Option<Value>,
    selection: Option<String>,
    start: usize,
    end: usize,
    total: usize,
    unit: &'static str,
    limit: usize,
    summaries: bool,
    nodes: usize,
    bytes: usize,
    attachments: Vec<Value>,
    sources: Vec<String>,
    page_aliases: IndexAliases,
    holes: Vec<(usize, String)>,
}

impl Builder<'_> {
    fn span(&self, a: AtomView<'_>) -> Span {
        Span {
            start: a.get_data().as_ptr() as usize - self.root.as_ptr() as usize,
            len: a.get_byte_size(),
        }
    }
    fn hole(&mut self, a: AtomView<'_>, label: &str) -> Value {
        let span = self.span(a);
        self.spans.insert(span.start, span);
        let index = self.holes.iter().position(|(id, _)| *id == span.start);
        let index = index.or_else(|| {
            if self.holes.len() >= 32 || span == self.focus {
                return None;
            }
            let index = self.holes.len();
            let label = match a {
                AtomView::Add(sum) => {
                    format!("Open sum [{}] · {} terms", index + 1, sum.get_nargs())
                }
                AtomView::Mul(product) => format!(
                    "Open product [{}] · {} factors",
                    index + 1,
                    product.get_nargs()
                ),
                _ => format!("{} [{}]", label, index + 1),
            };
            self.holes.push((span.start, label));
            Some(index)
        });
        let label = index.map_or_else(|| "⋯".to_owned(), |i| format!("⋯[{}]", i + 1));
        let marker = map([("kind", text("omission")), ("label", text(&label))]);
        if matches!(a, AtomView::Add(_)) {
            map([("kind", text("sum")), ("terms", Value::Array(vec![marker]))])
        } else {
            marker
        }
    }
    fn leaf(&mut self, a: AtomView<'_>, remaining: usize) -> PyResult<Value> {
        self.nodes = remaining;
        self.bytes -= a.get_byte_size();
        let atom = self.value.presentation_fragment(a);
        // Alias ownership lasts for the viewer, rather than restarting on every page.
        let discovered = IndexAliases::for_atom(&atom, &self.settings.index_style);
        for (key, (position, display)) in discovered.entries {
            if let Some(existing) = self.aliases.entries.get(&key) {
                self.page_aliases.entries.insert(key, existing.clone());
                continue;
            }
            let mut position = position;
            let mut display = display;
            if self.settings.index_style == "alphabet"
                && let (Some(mut p), Some(metadata)) =
                    (position, RepresentationMetadata::from_symbol(key.0))
            {
                while self
                    .aliases
                    .entries
                    .iter()
                    .any(|((rep, _), (_, d))| *rep == key.0 && *d == display)
                {
                    p = p.saturating_add(1);
                    let Some(next) = metadata.index_palette.resolve(p) else {
                        break;
                    };
                    display = next;
                }
                position = Some(p);
            }
            self.page_aliases
                .entries
                .insert(key.clone(), (position, display.clone()));
            self.aliases.entries.insert(key, (position, display));
        }
        self.sources
            .push(typst_custom_print_source(self.settings, &atom)?);
        let attachments = portable_attachments(&atom).map_err(PyRuntimeError::new_err)?;
        let Value::Map(mut envelope) =
            symbolica_typst_plugin::payload::atom_render_tree_value(&atom, &attachments)
                .map_err(|e| PyRuntimeError::new_err(e.to_string()))?
        else {
            unreachable!()
        };
        let mut result = None;
        for (k, v) in envelope.drain(..) {
            if k == text("root") {
                result = Some(v);
            } else if k == text("attachments")
                && let Value::Array(a) = v
            {
                self.attachments.extend(a);
            }
        }
        Ok(result.expect("render tree root"))
    }
    fn sequence(
        &mut self,
        span: Span,
        mut iter: ListIterator<'_>,
        depth: usize,
        kind: &'static str,
        field: &'static str,
        unit: &'static str,
    ) -> PyResult<Value> {
        if self.active.is_some() || self.cursor.is_some_and(|c| c.sum != span) {
            return Ok(self.hole(
                AtomView::from(&self.root[span.start..span.start + span.len]),
                "Open subexpression",
            ));
        }
        self.active = Some(span);
        self.unit = unit;
        self.total = iter.len();
        // Keep a compact expression outline separate from the selected inner
        // sum. Wrapping a sum inside its original fences makes those fences
        // stretch across the whole page and hides the surrounding factors.
        let marker = if span != self.focus {
            let original = AtomView::from(&self.root[span.start..span.start + span.len]);
            let marker = self.hole(original, "Open selected sum");
            self.selection = Some(
                self.holes
                    .iter()
                    .position(|(id, _)| *id == span.start)
                    .map_or_else(|| "Selected sum".to_owned(), |i| format!("Sum [{}]", i + 1)),
            );
            Some(marker)
        } else {
            None
        };
        let mut offset = 0;
        // Resume from a previously validated byte boundary, without rescanning terms.
        let mut rest = None;
        if let Some(cursor) = self.cursor {
            offset = cursor.index;
            rest = Some(&self.root[cursor.next.start..span.start + span.len]);
        }
        self.start = offset;
        self.end = offset;
        let mut terms = Vec::new();
        if offset > 0 {
            terms.push(map([("kind", text("omission"))]));
        }
        for _ in 0..self.limit {
            if self.end >= self.total {
                break;
            }
            let term = if let Some(data) = rest {
                // The next term's length is retained in the cursor; following
                // terms are read using the original list iterator's bounded skip.
                let original = AtomView::from(&self.root[span.start..span.start + span.len]);
                let original_iter = children(original).expect("sequence");
                // Construct a list cursor from a validated suffix (no data copying).
                // SAFETY: cursors are created only from this immutable
                // sum's iterator; every update advances by that atom's length.
                let mut suffix = unsafe { original_iter.resume_at(data, self.total - self.end) };
                let t = suffix.next().unwrap();
                rest = Some(&data[t.get_byte_size()..]);
                t
            } else {
                iter.next().unwrap()
            };
            if self.nodes < 16 || self.bytes < term.get_byte_size().min(256) {
                let next = self.span(term);
                self.next = Some(Cursor {
                    focus: self.focus,
                    sum: span,
                    next,
                    index: self.end,
                });
                break;
            }
            terms.push(self.build(term, depth + 1)?);
            self.end += 1;
        }
        if self.end < self.total {
            if self.next.is_none() {
                let next = if let Some(data) = rest {
                    Span {
                        start: data.as_ptr() as usize - self.root.as_ptr() as usize,
                        len: 0,
                    }
                } else {
                    self.span(iter.next().unwrap())
                };
                self.next = Some(Cursor {
                    focus: self.focus,
                    sum: span,
                    next,
                    index: self.end,
                });
            }
            terms.push(map([("kind", text("omission"))]));
        }
        let page = map([("kind", text(kind)), (field, Value::Array(terms))]);
        if let Some(marker) = marker {
            self.selected_page = Some(page);
            Ok(marker)
        } else {
            Ok(page)
        }
    }
    fn build(&mut self, a: AtomView<'_>, depth: usize) -> PyResult<Value> {
        if depth > MAX_DEPTH || self.nodes == 0 {
            return Ok(self.hole(a, "Open subexpression"));
        }
        let span = self.span(a);
        if self.summaries && span != self.focus && children(a).is_some() {
            return Ok(self.hole(a, "Open oversized term"));
        }
        let mut remaining = self.nodes;
        let too_many = large_sum(a, self.limit.min(100), &mut { MAX_NODES });
        let target = self.cursor.is_some_and(|c| c.sum == span);
        if !(self.summaries && children(a).is_some())
            && !too_many
            && !target
            && a.get_byte_size() <= self.bytes
            && fits(a, &mut remaining, depth)
        {
            return self.leaf(a, remaining);
        }
        self.nodes -= 1;
        match a {
            AtomView::Add(sum) => self.sequence(span, sum.iter(), depth, "sum", "terms", "Terms"),

            AtomView::Mul(m) => {
                // Avoid a single product with arbitrarily many factor placeholders.
                if m.get_nargs() > 32 {
                    if span == self.focus {
                        return self.sequence(
                            span,
                            m.iter(),
                            depth,
                            "product",
                            "factors",
                            "Factors",
                        );
                    }
                    return Ok(self.hole(a, "Open product factors"));
                }
                let mut factors = Vec::new();
                for c in m.iter() {
                    factors.push(self.build(c, depth + 1)?);
                }
                Ok(map([
                    ("kind", text("product")),
                    ("factors", Value::Array(factors)),
                ]))
            }
            AtomView::Pow(p) => {
                let (b, e) = p.get_base_exp();
                Ok(map([
                    ("kind", text("power")),
                    ("base", self.build(b, depth + 1)?),
                    ("exponent", self.build(e, depth + 1)?),
                ]))
            }
            AtomView::Fun(f) if span == self.focus => {
                self.sequence(span, f.iter(), depth, "inspection", "items", "Arguments")
            }
            _ => Ok(self.hole(a, "Open arguments")),
        }
    }
}

#[pymethods]
impl PageSource {
    #[pyo3(signature = (focus, cursor, limit, summaries=false))]
    fn page(
        &mut self,
        py: Python<'_>,
        focus: usize,
        cursor: Option<usize>,
        limit: usize,
        summaries: bool,
    ) -> PyResult<Py<PyDict>> {
        if !(1..=500).contains(&limit) {
            return Err(PyValueError::new_err("page size must be 1..500"));
        }
        let focus = *self
            .spans
            .get(&focus)
            .ok_or_else(|| PyValueError::new_err("unknown subexpression"))?;
        let cursor = cursor
            .map(|id| {
                self.cursors
                    .get(id)
                    .copied()
                    .ok_or_else(|| PyValueError::new_err("unknown cursor"))
            })
            .transpose()?;
        if cursor.is_some_and(|c| c.focus != focus) {
            return Err(PyValueError::new_err(
                "cursor belongs to another subexpression",
            ));
        }
        let root = self.value.expression().as_view().get_data();
        let atom = AtomView::from(&root[focus.start..focus.start + focus.len]);
        let mut builder = Builder {
            root,
            value: &self.value,
            settings: &self.settings,
            aliases: &mut self.aliases,
            spans: &mut self.spans,
            cursor,
            focus,
            next: None,
            active: None,
            selected_page: None,
            selection: None,
            start: 0,
            end: 0,
            total: 0,
            unit: "Terms",
            limit,
            summaries,
            nodes: MAX_NODES,
            bytes: MAX_BYTES,
            attachments: Vec::new(),
            sources: Vec::new(),
            page_aliases: IndexAliases::default(),
            holes: Vec::new(),
        };
        let node = builder.build(atom, 0)?;
        let result = PyDict::new(py);
        result.set_item("start", builder.start)?;
        result.set_item("end", builder.end)?;
        result.set_item("total", builder.total)?;
        result.set_item("unit", builder.unit)?;
        result.set_item("selection", &builder.selection)?;
        result.set_item("holes", &builder.holes)?;
        let next = builder.next.map(|c| {
            let id = self.cursors.len();
            self.cursors.push(c);
            id
        });
        result.set_item("next", next)?;
        let tree = map([
            ("kind", text("atom-render-tree")),
            ("root", node),
            ("page", builder.selected_page.unwrap_or(Value::Null)),
            ("attachments", Value::Array(builder.attachments)),
        ]);
        let mut bytes = Vec::new();
        ciborium::into_writer(&tree, &mut bytes)
            .map_err(|e| PyRuntimeError::new_err(e.to_string()))?;
        result.set_item("tree", PyBytes::new(py, &bytes))?;
        let mut source = typst_settings_source(&self.settings, &Atom::zero());
        source.push_str(&format!(
            "\n#let settings = (..settings, index-aliases: {})\n",
            builder
                .page_aliases
                .typst_source(&self.settings.index_style)
        ));
        // Merge callbacks from the bounded fragments, rather than losing earlier overrides.
        for s in builder.sources {
            source.push_str("\n#let saved-heads = settings.at(\"print-heads\", default: (:))\n#let saved-calls = settings.at(\"print-calls\", default: ())\n");
            source.push_str(&s);
            source.push_str("#let settings = (..settings, print-heads: (: ..saved-heads, ..settings.print-heads), print-calls: saved-calls + settings.print-calls)\n");
        }
        result.set_item("settings", source)?;
        result.set_item("notation", self.notation.as_deref().unwrap_or(NOTATION_TYP))?;
        result.set_item("render", RENDER_TYP)?;
        Ok(result.unbind())
    }
}

fn module(py: Python<'_>) -> PyResult<Bound<'_, PyModule>> {
    // Cache the module in sys.modules; all viewers share code, not expression state.
    let modules = py.import("sys")?.getattr("modules")?;
    let name = "_spenso_paging";
    if let Ok(m) = modules.get_item(name) {
        return Ok(m.cast_into::<PyModule>()?);
    }
    let code = CString::new(include_str!("../../typst/paging.py")).unwrap();
    let m = PyModule::from_code(py, &code, c"spenso_paging.py", c"_spenso_paging")?;
    m.add_function(wrap_pyfunction!(compile_typst, &m)?)?;
    m.setattr("WIDGET_ESM", include_str!("../../typst/paging.js"))?;
    m.setattr("NOTEBOOK_STYLE", NOTEBOOK_STYLE)?;
    modules.set_item(name, &m)?;
    Ok(m)
}
pub(crate) fn viewer(
    py: Python<'_>,
    value: Arc<SymbolicTensor<PartialStructure>>,
    settings: DisplaySettings,
    notation: Option<String>,
    size: usize,
) -> PyResult<Py<PyAny>> {
    if ![25, 100, 250, 500].contains(&size) {
        return Err(PyValueError::new_err(
            "page_size must be 25, 100, 250, or 500",
        ));
    }
    let len = value.expression().as_view().get_byte_size();
    let source = PageSource {
        value,
        settings,
        notation,
        spans: BTreeMap::from([(0, Span { start: 0, len })]),
        cursors: Vec::new(),
        aliases: IndexAliases::default(),
    };
    Ok(module(py)?
        .getattr("Pager")?
        .call1((Py::new(py, source)?, size))?
        .unbind())
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
    fn budgets_reject_deep_wide_and_large_atoms() {
        let small = symbolica::parse!("x+y");
        assert!(fits(small.as_view(), &mut { MAX_NODES }, 0));
        assert!(!fits(small.as_view(), &mut 1, 0));
        assert!(!fits(small.as_view(), &mut { MAX_NODES }, MAX_DEPTH + 1));
        let huge = Atom::add_many(
            (0..3000)
                .map(|n| symbolica::parse!("x").pow(n))
                .collect::<Vec<_>>(),
        );
        assert!(!fits(huge.as_view(), &mut { MAX_NODES }, 0));
    }
    #[test]
    fn continuation_resumes_the_exact_immutable_sum_suffix() {
        let atom = symbolica::parse!("a+b+c+d+e");
        let AtomView::Add(sum) = atom.as_view() else {
            panic!("sum")
        };
        let all = sum.iter().collect::<Vec<_>>();
        let root = atom.as_view().get_data();
        let start = all[2].get_data().as_ptr() as usize - root.as_ptr() as usize;
        // SAFETY: the offset comes from the original immutable iterator.
        let resumed = unsafe { sum.iter().resume_at(&root[start..], all.len() - 2) };
        assert_eq!(resumed.collect::<Vec<_>>(), all[2..]);
    }
}
