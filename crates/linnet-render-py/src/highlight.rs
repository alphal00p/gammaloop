use std::fmt::Write;

use linnest::svg::{
    CountertermScheme, HighlightDetails, HighlightOverlap, HighlightStyle, SlotPaint,
    UnknownCountertermScheme, UnknownHighlightOverlap, UnknownHighlightStyle,
};
use pyo3::exceptions::PyValueError;
use pyo3::types::{
    PyAny, PyAnyMethods, PyDict, PyDictMethods, PySequence, PySequenceMethods, PyString, PyTuple,
    PyTupleMethods,
};
use pyo3::{Bound, PyResult};
use HighlightStyle::{ActedOn, Amber, Middle, Treated};

/// A Python selection resolved by its owning graph. See [`extract_highlights`].
#[derive(Debug, PartialEq)]
pub struct SelectedHighlight<T> {
    /// Selection validated by the caller's callback.
    pub selection: T,
    /// Explicit or inferred style.
    pub style: HighlightStyle,
    /// Hover text and counterterm treatment.
    pub details: HighlightDetails,
}

/// Parse a Python overlap keyword. See [`extract_highlights`].
///
/// # Errors
/// Raises ValueError for an unknown mode.
pub fn highlight_overlap(mode: &str) -> PyResult<HighlightOverlap> {
    mode.parse()
        .map_err(|error: UnknownHighlightOverlap| PyValueError::new_err(error.to_string()))
}

/// Resolve a selection, a styled pair, or an ordered sequence of either.
/// The callback validates the owner of each selection. One, two, or three
/// unstyled entries receive amber, blue/red, or blue/violet/red respectively.
/// See [`linnest::svg::Scene::apply_highlights`] for the rendering contract.
///
/// # Errors
/// Returns a Python exception for an invalid selection, name, or sequence.
pub fn extract_highlights<'py, T>(
    highlight: &Bound<'py, PyAny>,
    overlap: HighlightOverlap,
    mut selection: impl FnMut(&Bound<'py, PyAny>) -> PyResult<T>,
) -> PyResult<Vec<SelectedHighlight<T>>> {
    let pair = |entry: &Bound<'py, PyAny>| -> PyResult<_> {
        let Ok(pair) = entry.cast::<PyTuple>() else {
            return Ok(None);
        };
        if pair.len() != 2 {
            return Ok(None);
        }
        let options = pair.get_item(1)?;
        let mut details = HighlightDetails::default();
        let role = if options.is_instance_of::<PyString>() {
            Some(options.extract::<String>()?)
        } else if let Ok(fields) = options.cast::<PyDict>() {
            for (key, value) in fields.iter() {
                match key.extract::<String>()?.as_str() {
                    "style" => {}
                    "label" => details.label = value.extract()?,
                    "description" => details.description = value.extract()?,
                    "contracted" => details.contracted = value.extract()?,
                    "scheme" => {
                        if let Some(name) = value.extract::<Option<String>>()? {
                            details.scheme =
                                name.parse().map_err(|error: UnknownCountertermScheme| {
                                    PyValueError::new_err(error.to_string())
                                })?;
                        }
                    }
                    key => {
                        return Err(PyValueError::new_err(format!(
                            "unknown highlight field {key:?}"
                        )))
                    }
                }
            }
            fields
                .get_item("style")?
                .map(|value| value.extract::<String>())
                .transpose()?
        } else {
            return Ok(None);
        };
        let style = role
            .map(|role| {
                role.parse().map_err(|error: UnknownHighlightStyle| {
                    PyValueError::new_err(error.to_string())
                })
            })
            .transpose()?;
        if !details.contracted && details.scheme != CountertermScheme::Unspecified {
            return Err(PyValueError::new_err(
                "scheme requires a contracted highlight layer",
            ));
        }
        let style = if details.contracted {
            if style.is_some_and(|style| style != Treated) {
                return Err(PyValueError::new_err(
                    "a contracted highlight layer requires style 'treated'",
                ));
            }
            Some(Treated)
        } else {
            style
        };
        Ok(Some((pair.get_item(0)?, style, details)))
    };
    let entries = if let Some(entry) = pair(highlight)? {
        vec![entry]
    } else if let Some(sequence) = highlight
        .cast::<PySequence>()
        .ok()
        .filter(|_| !highlight.is_instance_of::<PyString>())
    {
        (0..sequence.len()?)
            .map(|index| {
                let entry = sequence.get_item(index)?;
                Ok(match pair(&entry)? {
                    Some(entry) => entry,
                    None => (entry, None, HighlightDetails::default()),
                })
            })
            .collect::<PyResult<Vec<_>>>()?
    } else {
        vec![(highlight.clone(), None, HighlightDetails::default())]
    };
    let defaults: &[HighlightStyle] = match (
        overlap.band_spacing(),
        entries
            .iter()
            .filter(|(_, style, _)| style.is_none())
            .count(),
    ) {
        _ if overlap == HighlightOverlap::Lanes && entries.len() > 1 => &[],
        (Some(_), _) | (None, 0) => &[],
        (None, 1) => &[Amber],
        (None, 2) => &[ActedOn, Treated],
        (None, 3) => &[ActedOn, Middle, Treated],
        _ if overlap != HighlightOverlap::Casings => &[],
        _ => {
            return Err(PyValueError::new_err(
                "provide explicit highlight styles for all but at most three layers",
            ))
        }
    };
    let mut defaults = defaults.iter().copied();
    entries
        .into_iter()
        .enumerate()
        .map(|(index, (entry, style, details))| {
            Ok(SelectedHighlight {
                selection: selection(&entry)?,
                details,
                style: style
                    .unwrap_or_else(|| defaults.next().unwrap_or(HighlightStyle::Slot(index + 1))),
            })
        })
        .collect()
}

/// One CSS palette role. See [`theme_svg`].
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ThemePaint<'a> {
    role: &'a str,
    light: &'a str,
    dark: &'a str,
}

impl<'a> ThemePaint<'a> {
    /// Associate a CSS role with its two paints. See [`theme_svg`].
    pub const fn new(role: &'a str, light: &'a str, dark: &'a str) -> Self {
        Self { role, light, dark }
    }
}

/// Byte formatting of public SVG exports. See [`theme_svg`].
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ThemeFormat {
    /// Legacy Linnet subgraph style block.
    Linnet,
    /// Compact FeynKit style block.
    FeynKit,
}

/// Apply shared light, dark, and notebook theme selectors to an SVG.
/// Palette rows and class names are trusted renderer inputs.
///
/// ```
/// use linnet_render_py::{ThemeFormat, ThemePaint, theme_svg};
/// let svg = theme_svg("<svg xmlns=\"http://www.w3.org/2000/svg\"></svg>".into(),
///     &[ThemePaint::new("ink", "#000000", "#ffffff")], &[], ThemeFormat::Linnet);
/// assert!(svg.contains("data-jp-theme-light"));
/// ```
pub fn theme_svg(
    svg: String,
    palette: &[ThemePaint<'_>],
    highlights: &[HighlightStyle],
    format: ThemeFormat,
) -> String {
    let linnet = matches!(format, ThemeFormat::Linnet);
    let (class, prefix) = if linnet {
        ("linnet-subgraph", "linnet")
    } else {
        ("feynkit-diagram-svg", "feynkit")
    };
    let gap = if linnet { " " } else { "" };
    let separator = if linnet { "\n" } else { "" };
    let mut styles = String::new();
    let mut light = String::new();
    let mut dark = String::new();
    let mut rules = String::new();
    let named = HighlightStyle::LAYERED.into_iter().filter_map(|style| {
        highlights.contains(&style).then_some(())?;
        Some(ThemePaint::new(
            style.role(),
            style.light_paint()?,
            style.dark_paint()?,
        ))
    });
    let slots = highlights
        .iter()
        .filter_map(|style| match *style {
            HighlightStyle::Slot(slot) => Some((style.to_string(), SlotPaint::new(slot))),
            _ => None,
        })
        .collect::<Vec<_>>();
    let slot_rows = slots
        .iter()
        .map(|(role, colors)| ThemePaint::new(role, &colors.light, &colors.dark));
    for paint in palette.iter().copied().chain(named).chain(slot_rows) {
        let variable = format!("--{prefix}-{}", paint.role);
        write!(light, "{variable}:{gap}{};{gap}", paint.light).expect("writing a string");
        write!(dark, "{variable}:{gap}{};{gap}", paint.dark).expect("writing a string");
        for attribute in ["fill", "stroke"] {
            let value = if linnet {
                format!("var({variable})")
            } else {
                format!("var({variable},{})", paint.light)
            };
            write!(rules, ".{class} [{attribute}=\"{}\"]{gap}{{{gap}{attribute}:{gap}{value};{gap}}}{separator}", paint.light).expect("writing a string");
        }
    }
    if svg.contains("data-linnet-counterterm=") || svg.contains("data-linnet-paper=") {
        let paper = format!("--{prefix}-counterterm-paper");
        write!(light, "{paper}:#ffffff;{gap}").expect("writing a string");
        write!(dark, "{paper}:#1b1d23;{gap}").expect("writing a string");
        write!(
            rules,
            ".{class} [data-linnet-counterterm]:not([data-scheme=on-shell]) > circle{{fill:var({paper},#ffffff);}}.{class} [data-linnet-paper]{{stroke:var({paper},#ffffff);}}{separator}"
        )
        .expect("writing a string");
    }
    if !linnet {
        styles.push_str(&rules);
    }
    writeln!(styles, ".{class}{gap}{{{gap}{light}}}").expect("writing a string");
    if linnet {
        styles.push_str(&rules);
    }
    if linnet {
        writeln!(
            styles,
            "@media (prefers-color-scheme: dark) {{\n  .{class} {{ {dark}}}\n}}"
        )
        .expect("writing a string");
    } else {
        writeln!(
            styles,
            "@media(prefers-color-scheme:dark){{.{class}{{{dark}}}}}"
        )
        .expect("writing a string");
    }
    for (theme, jupyter, variables) in [("dark", "false", &dark), ("light", "true", &light)] {
        let join = if linnet { "\n" } else { "" };
        let declarations = if linnet {
            format!("\n  {}\n", variables.trim_end())
        } else {
            variables.to_owned()
        };
        writeln!(styles, ".{class}[data-theme=\"{theme}\"],{join}:is(.{theme},{gap}[data-theme=\"{theme}\"],{gap}[data-jp-theme-light=\"{jupyter}\"]) .{class}{gap}{{{declarations}}}").expect("writing a string");
    }
    let svg_gap = if linnet { " " } else { "" };
    svg.replacen("<svg ", &format!("<svg class=\"{class}\" "), 1)
        .replacen('>', &format!(">{svg_gap}<style>{styles}</style>"), 1)
}
