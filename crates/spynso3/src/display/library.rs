//! Library discovery uses the same symbolic printer and component viewer as tensors.
use std::sync::atomic::{AtomicU64, Ordering};

use pyo3::types::{PyString, PyTuple};
use spenso::structure::dimension::Dimension;

use super::*;
use crate::{expression::TensorExpression, library::SpensorLibrary};

const STYLE: &str = include_str!("../../typst/library.css");
const COPY: &str = include_str!("../../typst/library-copy.js");
const MAX_ENTRIES: usize = 32;
const MAX_HTML_BYTES: usize = 8_000_000;
static NEXT_CATALOGUE: AtomicU64 = AtomicU64::new(0);

struct Recipe {
    access: String,
    expression: String,
    indices: Vec<usize>,
}

impl Recipe {
    fn quote(py: Python<'_>, text: &str) -> PyResult<String> {
        PyString::new(py, text).repr()?.extract()
    }

    fn argument(py: Python<'_>, atom: &Atom) -> PyResult<String> {
        if let Ok(number) = i64::try_from(atom.as_view()) {
            Ok(number.to_string())
        } else if let AtomView::Var(variable) = atom.as_view() {
            Ok(format!(
                "S({})",
                Self::quote(py, variable.get_symbol().get_name())?
            ))
        } else {
            Ok(format!(
                "E({})",
                Self::quote(py, &atom.to_canonical_string())?
            ))
        }
    }

    fn new(py: Python<'_>, tensor: &Spensor) -> PyResult<Self> {
        let name = tensor.descriptor_name.expect("library tensors are named");
        let mut arguments = tensor
            .descriptor_args
            .iter()
            .map(|arg| Self::argument(py, arg))
            .collect::<PyResult<Vec<_>>>()?;
        let mut counters = HashMap::new();
        let mut indices = Vec::new();
        for slot in tensor.descriptor.interface.logical_slots() {
            let rep = slot.rep();
            let base = rep.rep.base();
            let dimension = match rep.dim {
                Dimension::Concrete(dimension) => dimension.to_string(),
                _ => Self::argument(py, &rep.dim.to_symbolic())?,
            };
            let constructor = match base.symbol().get_name() {
                "spenso::euc" => "euc",
                "spenso::mink" => "mink",
                "spenso::bis" => "bis",
                "spenso::cof" => "cof",
                "spenso::coad" => "coad",
                "spenso::cos" => "cos",
                _ => "",
            };
            let mut source = if constructor.is_empty() {
                format!(
                    "Representation({}, {dimension}, is_self_dual={})",
                    Self::quote(py, base.symbol().get_name())?,
                    if base.is_self_dual() { "True" } else { "False" }
                )
            } else {
                format!("Representation.{constructor}({dimension})")
            };
            if !rep.rep.is_base() {
                source.push_str(".dual()");
            }
            arguments.push(source);
            let next = counters.entry(base).or_insert_with(|| {
                match base.metadata().map(|metadata| metadata.index_palette) {
                    Some(IndexPalette::Cyclic { start, .. }) => start,
                    _ => 0,
                }
            });
            indices.push(*next);
            *next += 1;
        }
        let call = if arguments.is_empty() {
            String::new()
        } else {
            format!("\n    {},\n", arguments.join(",\n    "))
        };
        let access = format!(
            "from symbolica import E, S\nfrom symbolica.community.spenso import Representation, TensorName\n\nkey = TensorName({})({call})\ntensor = library[key]",
            Self::quote(py, name.get_name())?,
        );
        let labels = indices
            .iter()
            .map(usize::to_string)
            .collect::<Vec<_>>()
            .join(", ");
        let expression = format!("# After the Access example\nexpr = key({labels})\nexpr");
        Ok(Self {
            access,
            expression,
            indices,
        })
    }

    fn html(&self, id: &str) -> String {
        let print = "# After the Expression example\nexpr                         # Notebook display\nprint(expr)                  # Text\nprint(expr.to_typst())        # Typst source\nprint(expr.to_latex())        # LaTeX source";
        let mut tabs = String::new();
        let mut panels = String::new();
        for (i, (label, code)) in [
            ("Access", self.access.as_str()),
            ("Expression", self.expression.as_str()),
            ("Print", print),
        ]
        .into_iter()
        .enumerate()
        {
            tabs.push_str(&format!("<label><input type=\"radio\" name=\"{id}-recipe\" value=\"{i}\" {}><span>{label}</span></label>", if i == 0 { "checked" } else { "" }));
            panels.push_str(&format!("<div class=\"sl-code\" data-recipe=\"{i}\"><button type=\"button\" class=\"sl-copy\" onclick=\"{}\" aria-label=\"Copy {label} example\">Copy</button><pre><code>{}</code></pre></div>", escape_html(COPY), escape_html(code)));
        }
        format!(
            "<details class=\"sl-python\"><summary>Python</summary><div class=\"sl-recipes\"><div class=\"sl-tabs\" role=\"group\" aria-label=\"Python examples\">{tabs}</div>{panels}<p>Use your library variable in place of <code>library</code>.</p></div></details>"
        )
    }
}

pub(crate) fn to_html(
    py: Python<'_>,
    library: &SpensorLibrary,
    settings: &DisplaySettings,
) -> PyResult<String> {
    let id = NEXT_CATALOGUE.fetch_add(1, Ordering::Relaxed);
    let root = format!("spenso-library-{id}");
    let selector = format!("[data-spenso-library=\"{id}\"]");
    let keys = library.keys(py)?;
    let count = keys.len();
    let mut chips = String::new();
    let mut panels = String::new();
    let mut rules = String::new();
    let mut shown = 0;
    for (index, key) in keys.into_iter().take(MAX_ENTRIES).enumerate() {
        if panels.len() >= MAX_HTML_BYTES {
            break;
        }
        let tensor = library.__getitem__(py, key.bind(py).extract()?)?;
        let mut tensor = tensor.borrow(py).clone();
        let recipe = Recipe::new(py, &tensor)?;
        // Index a snapshot for presentation only; the stored key stays unresolved.
        let indexed = key
            .bind(py)
            .call1(PyTuple::new(py, &recipe.indices)?)?
            .extract::<Py<TensorExpression>>()?;
        tensor.descriptor = TensorExpression::structured(&indexed.borrow(py));
        let name = tensor.descriptor_name.expect("library tensors are named");
        let title = name.get_name();
        let svg = structured_to_svg(py, &tensor.descriptor, settings, None)?;
        let shape = tensor
            .descriptor
            .interface
            .logical_slots()
            .into_iter()
            .map(|slot| slot.rep().dim.to_string())
            .collect::<Vec<_>>()
            .join(" × ");
        let shape = if shape.is_empty() { "Scalar" } else { &shape };
        let control = format!("{root}-{index}");
        chips.push_str(&format!("<label class=\"sl-chip\" title=\"{}\"><input type=\"radio\" name=\"{root}\" value=\"{index}\" aria-label=\"{} · {}\" aria-controls=\"{control}\"><span class=\"sl-formula\">{svg}</span><span class=\"sl-shape\">{}</span></label>", escape_html(title), escape_html(title), escape_html(shape), escape_html(shape)));
        rules.push_str(&format!("{selector}:has(.sl-chip input[value=\"{index}\"]:checked) > .sl-panel[data-entry=\"{index}\"]{{display:block}}"));
        let viewer = concrete_tensor_to_html(py, &tensor, settings, None)?;
        panels.push_str(&format!("<section class=\"sl-panel\" id=\"{control}\" data-entry=\"{index}\"><div class=\"sl-panel-heading\"><code>{}</code><label class=\"sl-close\"><input type=\"radio\" name=\"{root}\" value=\"closed-{index}\" aria-label=\"Close tensor\"><span aria-hidden=\"true\">×</span></label></div>{}{viewer}</section>", escape_html(title), recipe.html(&control)));
        shown += 1;
    }
    let empty = if count == 0 {
        "<p class=\"sl-empty\">No stored tensors. Register a named Tensor to add it to this library.</p>"
    } else {
        ""
    };
    let limited = if shown < count {
        format!(
            "<p class=\"sl-note\">Showing {shown} of {count} tensors. Use library.keys() or library.items() to access every entry.</p>"
        )
    } else {
        String::new()
    };
    let factories = library.library.generic_len();
    Ok(format!(
        "<style>{STYLE}{rules}</style><div data-spenso-library=\"{id}\"><div class=\"sl-heading\"><code>TensorLibrary</code><span>{count} stored tensors</span></div><div class=\"sl-catalogue\">{chips}</div>{empty}{panels}{limited}<details class=\"sl-factories\"><summary>{factories} dimension-dependent factories</summary><p>Components are generated for a concrete signature on lookup. Factories are excluded from stored entries.</p><code>library[TensorExpression.g(Representation.mink(4))]</code></details></div>"
    ))
}
