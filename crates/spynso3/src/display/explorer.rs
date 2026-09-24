//! Bounded, read-only component snapshots for self-contained notebook exploration.
use super::*;

const STYLE: &str = include_str!("../../typst/explorer.css");
const SCRIPT: &str = include_str!("../../typst/explorer.js");
const MAX_COMPONENTS: usize = 512;
const MAX_HTML_BYTES: usize = 2_000_000;
const MAX_FORMULA_BYTES: usize = 4_096;

type Component = (usize, String, String);
type Entry = (Vec<usize>, Component);

fn payload_bytes(value: &AtomViewOrConcrete<'_, RealOrComplexRef<'_, f64>>) -> usize {
    match value {
        AtomViewOrConcrete::Atom(atom) => atom.get_byte_size(),
        AtomViewOrConcrete::Concrete(RealOrComplexRef::Real(value)) => size_of_val(*value),
        AtomViewOrConcrete::Concrete(RealOrComplexRef::Complex(value)) => size_of_val(*value),
    }
}

struct Renderer<'a, 'py> {
    py: Python<'py>,
    settings: &'a DisplaySettings,
    notation_source: Option<&'a str>,
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
            return Ok((bytes, text.clone(), format!("<code>{text}</code>")));
        }
        let formula = match &value {
            AtomViewOrConcrete::Atom(atom) => {
                let atom = (*atom).to_owned();
                if !self.formulas.contains_key(&atom) {
                    let svg = atom_to_svg(self.py, &atom, self.settings, self.notation_source)?;
                    self.formulas.insert(atom.clone(), svg);
                }
                Some(self.formulas[&atom].clone())
            }
            AtomViewOrConcrete::Concrete(_) => None,
        };
        let plain = format_tensor_value(value, TensorDisplayMode::Plain, self.settings);
        let html = formula.unwrap_or_else(|| format!("<code>{}</code>", escape_html(&plain)));
        Ok((bytes, plain, html))
    }
}

impl ConcreteTensorView<'_> {
    fn axis_labels(&self) -> Vec<String> {
        use spenso::structure::slot::{SlotMatch, SlotMatcher};

        let mut matcher = SlotMatcher::default();

        let aliases =
            IndexAliases::for_descriptor(&self.tensor.descriptor, &self.settings.index_style);
        self.tensor
            .descriptor
            .interface
            .logical_slots()
            .into_iter()
            .map(|slot| {
                let port = composition::port_atom(slot);
                let port = aliases.presentation_atom(&port, &self.settings.index_style);
                let SlotMatch::Explicit(view) = matcher.classify(port.as_view()) else {
                    return port.format_string(
                        &display_options_with_settings(TensorDisplayMode::Plain, self.settings),
                        PrintState::new(),
                    );
                };
                let index = view.index();
                let display = usize::try_from(index)
                    .ok()
                    .and_then(|position| {
                        slot.rep_name().metadata()?.index_palette.resolve(position)
                    })
                    .or_else(|| match index {
                        AtomView::Var(variable) => IndexDisplay::from_symbol(variable.get_symbol()),
                        _ => None,
                    });
                display
                    .map(|label| label.to_native_string())
                    .unwrap_or_else(|| {
                        index.to_owned().format_string(
                            &display_options_with_settings(TensorDisplayMode::Plain, self.settings),
                            PrintState::new(),
                        )
                    })
            })
            .collect()
    }

    fn explorer_coordinates(&self) -> (Vec<Vec<usize>>, usize) {
        if tensor_is_sparse(self.tensor) {
            let preview = sparse_preview(self);
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
            let coordinates = displayed_indices(size, MAX_COMPONENTS / 2, size > MAX_COMPONENTS)
                .into_iter()
                .filter_map(|index| match index {
                    DisplayIndex::Index(index) => {
                        Some(expanded_from_flat(index, self.layout.logical_shape()))
                    }
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
    notation_source: Option<&str>,
) -> PyResult<String> {
    let view = ConcreteTensorView::new(tensor, TensorDisplayMode::Plain, settings)
        .ok_or_else(|| PyValueError::new_err("invalid tensor logical layout"))?;
    let (coordinates, stored) = view.explorer_coordinates();
    let mut renderer = Renderer {
        py,
        settings,
        notation_source,
        formulas: HashMap::new(),
    };
    let default = tensor_sparse_default(tensor)
        .map(|value| renderer.component(value))
        .transpose()?;
    let mut entries: Vec<Entry> = Vec::with_capacity(coordinates.len());
    let mut html_bytes = 0;
    for coordinate in coordinates {
        let storage = view
            .layout
            .logical_expanded_to_storage_flat(&coordinate)
            .map_err(|error| PyValueError::new_err(error.to_string()))?;
        let value = tensor_value_at_storage(tensor, storage)
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
        atom_to_svg(py, &descriptor, settings, notation_source)?
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
            "<style>{STYLE}</style></head><body>",
            "<main id=\"spenso-explorer\"><div class=\"descriptor\">{heading}</div>",
            "<div id=\"shape\"></div><div id=\"controls\"></div>",
            "<div id=\"detail\" aria-live=\"polite\"></div>",
            "<div id=\"plot\"></div><div id=\"legend\"></div>",
            "<div id=\"status\"></div></main>",
            "<script id=\"tensor-data\" type=\"application/json\">{json}</script>",
            "<script>{SCRIPT}</script></body></html>"
        ),
        STYLE = STYLE,
        SCRIPT = SCRIPT,
        heading = heading,
        json = json,
    );
    // Notebook HTML formatters do not execute script tags. A sandboxed srcdoc
    // keeps this portable across marimo and Jupyter, without a server or widget dependency.
    let height = if view.layout.logical_shape().len() == 3 {
        760
    } else {
        520
    };
    Ok(format!(
        r#"<iframe data-spenso-explorer title="Tensor component explorer" sandbox="allow-scripts" style="width:100%;height:{height}px;min-height:240px;resize:vertical;overflow:auto;border:0;color-scheme:light dark" srcdoc="{}"></iframe>"#,
        escape_html(&document)
    ))
}

#[cfg(test)]
mod tests {
    use super::*;
    use spenso::algebra::complex::Complex;

    #[test]
    fn memory_shading_uses_encoded_bytes_or_native_numeric_storage() {
        let small = Atom::num(1);
        let large = symbolica::parse!("(x+y+z)^4").expand();
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
