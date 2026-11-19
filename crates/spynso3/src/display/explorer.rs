//! Bounded, read-only component snapshots for self-contained notebook exploration.
use crate::Spensor;
use crate::display;
use crate::display::{
    ConcreteTensorView, DisplayIndex, DisplaySettings, IndexAliases, NOTEBOOK_STYLE,
    TensorDisplayMode,
};
use pyo3::exceptions::PyValueError;
use pyo3::types::{PyAnyMethods, PyDict, PyDictMethods, PyModule};
use pyo3::{PyResult, Python};
use spenso::algebra::complex::RealOrComplexRef;
use spenso::structure::partial::PartialStructureExt;
use spenso::tensors::parametric::AtomViewOrConcrete;
use std::collections::HashMap;
use std::collections::hash_map::Entry as HashMapEntry;
use std::mem;
use symbolica::atom::Atom;

const STYLE: &str = include_str!("../../typst/explorer.css");
const SCRIPT: &str = include_str!("../../typst/explorer.js");
const MAX_COMPONENTS: usize = 512;
const MAX_FORMULA_BYTES: usize = 4_096;
const MAX_HTML_BYTES: usize = 2_000_000;

type Component = (usize, String, String);
type Entry = (Vec<usize>, Component);

fn payload_bytes(value: &AtomViewOrConcrete<'_, RealOrComplexRef<'_, f64>>) -> usize {
    match value {
        AtomViewOrConcrete::Atom(atom) => atom.get_byte_size(),
        AtomViewOrConcrete::Concrete(RealOrComplexRef::Real(value)) => mem::size_of_val(*value),
        AtomViewOrConcrete::Concrete(RealOrComplexRef::Complex(value)) => mem::size_of_val(*value),
    }
}

struct Renderer<'a, 'py> {
    py: Python<'py>,
    settings: &'a DisplaySettings,
    formulas: HashMap<Atom, String>,
}

impl Renderer<'_, '_> {
    fn component(
        &mut self,
        value: AtomViewOrConcrete<'_, RealOrComplexRef<'_, f64>>,
    ) -> PyResult<Component> {
        let bytes = payload_bytes(&value);
        if bytes > MAX_FORMULA_BYTES {
            let text = format!("Expression ({bytes} bytes; full formula omitted from preview)");
            let html = format!("<code>{text}</code>");
            return Ok((bytes, text, html));
        }
        let formula = match &value {
            AtomViewOrConcrete::Atom(atom) => {
                let html = match self.formulas.entry((*atom).to_owned()) {
                    HashMapEntry::Occupied(entry) => entry.into_mut(),
                    HashMapEntry::Vacant(entry) => {
                        let html =
                            display::atom_to_html(self.py, entry.key(), self.settings, None)?;
                        entry.insert(html)
                    }
                };
                // The cached formula and each component record own their markup.
                Some(html.clone())
            }
            AtomViewOrConcrete::Concrete(_) => None,
        };
        let html = match formula {
            Some(formula) => formula,
            None => display::format_tensor_value_html(self.py, &value, self.settings)?,
        };
        let plain = display::format_tensor_value(value, TensorDisplayMode::Plain, self.settings);
        Ok((bytes, plain, html))
    }
}

impl ConcreteTensorView<'_> {
    fn axis_labels(&self) -> Vec<String> {
        let aliases =
            IndexAliases::for_descriptor(&self.tensor.descriptor, &self.settings.index_style);
        self.tensor
            .descriptor
            .structure()
            .logical_slots()
            .into_iter()
            .map(|slot| aliases.port_label(slot, self.settings))
            .collect()
    }

    fn explorer_coordinates(&self) -> (Vec<Vec<usize>>, usize) {
        if display::tensor_is_sparse(self.tensor) {
            let preview = display::sparse_preview(self);
            (
                preview
                    .entries
                    .into_iter()
                    .map(|entry| entry.logical)
                    .collect(),
                preview.stored,
            )
        } else {
            let size = self.layout.size();
            let coordinates =
                display::displayed_indices(size, MAX_COMPONENTS / 2, size > MAX_COMPONENTS)
                    .into_iter()
                    .filter_map(|index| match index {
                        DisplayIndex::Index(index) => Some(display::expanded_from_flat(
                            index,
                            self.layout.logical_shape(),
                        )),
                        DisplayIndex::Ellipsis => None,
                    })
                    .collect();
            (coordinates, size)
        }
    }
}

pub(super) fn to_html(
    py: Python<'_>,
    tensor: &Spensor,
    settings: &DisplaySettings,
) -> PyResult<String> {
    let view = ConcreteTensorView::new(tensor, TensorDisplayMode::Plain, settings)
        .ok_or_else(|| PyValueError::new_err("invalid tensor logical layout"))?;
    let (coordinates, stored) = view.explorer_coordinates();
    let mut renderer = Renderer {
        py,
        settings,
        formulas: HashMap::new(),
    };
    let default = display::tensor_sparse_default(tensor)
        .map(|value| renderer.component(value))
        .transpose()?;
    let mut entries: Vec<Entry> = Vec::with_capacity(coordinates.len());
    let mut html_bytes = 0;
    for coordinate in coordinates {
        let storage = view
            .layout
            .logical_expanded_to_storage_flat(&coordinate)
            .map_err(|error| PyValueError::new_err(error.to_string()))?;
        let value = display::tensor_value_at_storage(tensor, storage)
            .ok_or_else(|| PyValueError::new_err("missing stored tensor component"))?;
        let component = renderer.component(value)?;
        html_bytes += component.1.len() + component.2.len();
        if html_bytes > MAX_HTML_BYTES {
            break;
        }
        entries.push((coordinate, component));
    }
    let payload = PyDict::new(py);
    payload.set_item("shape", view.layout.logical_shape())?;
    payload.set_item("stored", stored)?;
    payload.set_item("complete", entries.len() == stored)?;
    payload.set_item("sparse", default.is_some())?;
    payload.set_item("default", default)?;
    payload.set_item("entries", entries)?;
    payload.set_item("axes", view.axis_labels())?;
    let descriptor = tensor.descriptor.presentation_atom();
    let heading = if descriptor.as_view().get_byte_size() <= MAX_FORMULA_BYTES {
        display::structured_to_html(py, &tensor.descriptor, settings, None)?
    } else {
        "<span>Tensor</span>".to_owned()
    };
    let json: String = PyModule::import(py, "json")?
        .call_method1("dumps", (payload,))?
        .extract()?;
    // The data is JSON, never executable source. Escape '<' even inside srcdoc:
    // tensor names and custom printer output must not close the data script.
    let json = json.replace('<', "\\u003c");
    let document = format!(
        concat!(
            "<!doctype html><html><head><meta charset=\"utf-8\">",
            "<meta name=\"viewport\" content=\"width=device-width,initial-scale=1\">",
            "<meta http-equiv=\"Content-Security-Policy\" content=\"default-src 'none'; ",
            "script-src 'unsafe-inline'; style-src 'unsafe-inline'; img-src data:;\">",
            "<style>{NOTEBOOK_STYLE}\n{STYLE}</style></head><body>",
            "<main id=\"spenso-explorer\"><div class=\"tensor-metadata\">",
            "<code class=\"python-type\">Tensor</code><span id=\"shape\"></span></div>",
            "<section aria-labelledby=\"expression-label\">",
            "<div id=\"expression-label\">Expression</div>",
            "<div class=\"descriptor\">{heading}</div></section>",
            "<div id=\"controls\"></div>",
            "<div id=\"detail\" aria-live=\"polite\"></div>",
            "<div id=\"plot\"></div><div id=\"legend\"></div>",
            "<div id=\"status\"></div></main>",
            "<script id=\"tensor-data\" type=\"application/json\">{json}</script>",
            "<script>{SCRIPT}</script></body></html>"
        ),
        STYLE = STYLE,
        NOTEBOOK_STYLE = NOTEBOOK_STYLE,
        SCRIPT = SCRIPT,
        heading = heading,
        json = json,
    );
    // Notebook HTML formatters do not execute script tags. A sandboxed srcdoc
    // keeps this portable across marimo and Jupyter, without a server or widget dependency.
    // Explicitly forward the page's theme as well: WebKit does not reliably
    // propagate color-scheme into sandboxed frames. The sandbox stays opaque.
    let height = if view.layout.logical_shape().len() == 3 {
        760
    } else {
        520
    };
    Ok(format!(
        r#"<iframe data-spenso-explorer title="Tensor component explorer" sandbox="allow-scripts" style="display:block;width:100%;max-width:100%;height:{height}px;min-height:240px;resize:vertical;overflow:auto;border:0;color-scheme:inherit" onload="{}" srcdoc="{}"></iframe>"#,
        display::escape_html(include_str!("../../typst/explorer-host.js")),
        display::escape_html(&document)
    ))
}

#[cfg(test)]
mod tests {
    use super::*;
    use spenso::algebra::complex::Complex;
    use symbolica::atom::AtomCore;
    use symbolica::parse;

    #[test]
    fn memory_shading_uses_encoded_bytes_or_native_numeric_storage() {
        let small = Atom::num(1);
        let large = parse!("(x+y+z)^4").expand();
        let real = 1e300;
        let complex = Complex::new(1., 2.);
        assert_eq!(
            payload_bytes(&AtomViewOrConcrete::Atom(small.as_view())),
            small.as_view().get_byte_size()
        );
        assert!(
            payload_bytes(&AtomViewOrConcrete::Atom(large.as_view()))
                > payload_bytes(&AtomViewOrConcrete::Atom(small.as_view()))
        );
        assert_eq!(
            payload_bytes(&AtomViewOrConcrete::Concrete(RealOrComplexRef::Real(&real))),
            8
        );
        assert_eq!(
            payload_bytes(&AtomViewOrConcrete::Concrete(RealOrComplexRef::Complex(
                &complex
            ))),
            16
        );
    }
}
