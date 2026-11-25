use linnest::svg::HighlightStyle;
use pyo3::exceptions::PyValueError;
use pyo3::types::{
    PyAny, PyAnyMethods, PySequence, PySequenceMethods, PyString, PyTuple, PyTupleMethods,
};
use pyo3::{Bound, PyResult};

/// Resolve a Python selection or an ordered sequence of selections and
/// `(selection, style)` tuples. The caller validates each selection's owner.
/// See [`linnest::svg::Scene::apply_highlights`] for the rendering contract.
pub fn extract_highlights<'py, T>(
    highlight: &Bound<'py, PyAny>,
    mut selection: impl FnMut(&Bound<'py, PyAny>) -> PyResult<T>,
) -> PyResult<Vec<(T, HighlightStyle)>> {
    if let Some(sequence) = highlight
        .cast::<PySequence>()
        .ok()
        .filter(|_| !highlight.is_instance_of::<PyString>())
    {
        (0..sequence.len()?)
            .map(|index| {
                let entry = sequence.get_item(index)?;
                if let Some(pair) = entry.cast::<PyTuple>().ok().filter(|pair| pair.len() == 2) {
                    let style = pair.get_item(1)?.extract::<String>()?;
                    let style = match style.as_str() {
                        "amber" => HighlightStyle::Amber,
                        "acted-on" => HighlightStyle::ActedOn,
                        "middle" => HighlightStyle::Middle,
                        "treated" => HighlightStyle::Treated,
                        _ => return Err(PyValueError::new_err(
                            "highlight style must be 'amber', 'acted-on', 'middle', or 'treated'",
                        )),
                    };
                    Ok((selection(&pair.get_item(0)?)?, style))
                } else {
                    Ok((selection(&entry)?, HighlightStyle::Amber))
                }
            })
            .collect()
    } else {
        Ok(vec![(selection(highlight)?, HighlightStyle::Amber)])
    }
}

/// CSS for the layered palette, scoped to the supplied internal SVG class.
/// See [`linnest::svg::Scene::apply_highlights`] for the named stroke treatments.
pub fn highlight_css(svg: &str, class: &str, variable_prefix: &str) -> String {
    let mut css = String::new();
    let mut light = String::new();
    let mut dark = String::new();
    for (role, fallback, dark_color) in HighlightStyle::THEME {
        let attribute = format!("stroke=\"{fallback}\"");
        if svg.contains(&attribute) {
            css.push_str(&format!(
                ".{class} [{attribute}]{{stroke:var(--{variable_prefix}-{role},{fallback});}}"
            ));
            light.push_str(&format!("--{variable_prefix}-{role}:{fallback};"));
            dark.push_str(&format!("--{variable_prefix}-{role}:{dark_color};"));
        }
    }
    if light.is_empty() {
        return css;
    }
    css.push_str(&format!(
        r#".{class}{{{light}}}
@media(prefers-color-scheme:dark){{.{class}{{{dark}}}}}
.{class}[data-theme="dark"],:is(.dark,[data-theme="dark"],[data-jp-theme-light="false"]) .{class}{{{dark}}}
.{class}[data-theme="light"],:is(.light,[data-theme="light"],[data-jp-theme-light="true"]) .{class}{{{light}}}"#
    ));
    css
}
