use std::fmt::Write;

use linnest::svg::{HighlightStyle, UnknownHighlightStyle};
use pyo3::exceptions::PyValueError;
use pyo3::types::{
    PyAny, PyAnyMethods, PySequence, PySequenceMethods, PyString, PyTuple, PyTupleMethods,
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
    mut selection: impl FnMut(&Bound<'py, PyAny>) -> PyResult<T>,
) -> PyResult<Vec<SelectedHighlight<T>>> {
    let pair = |entry: &Bound<'py, PyAny>| -> PyResult<_> {
        let Ok(pair) = entry.cast::<PyTuple>() else {
            return Ok(None);
        };
        if pair.len() == 2 && pair.get_item(1)?.is_instance_of::<PyString>() {
            let role: String = pair.get_item(1)?.extract()?;
            let style = role
                .parse()
                .map_err(|error: UnknownHighlightStyle| PyValueError::new_err(error.to_string()))?;
            Ok(Some((pair.get_item(0)?, style)))
        } else {
            Ok(None)
        }
    };
    let entries = if let Some((entry, style)) = pair(highlight)? {
        vec![(entry, Some(style))]
    } else if let Some(sequence) = highlight
        .cast::<PySequence>()
        .ok()
        .filter(|_| !highlight.is_instance_of::<PyString>())
    {
        (0..sequence.len()?)
            .map(|index| {
                let entry = sequence.get_item(index)?;
                Ok(match pair(&entry)? {
                    Some((entry, style)) => (entry, Some(style)),
                    None => (entry, None),
                })
            })
            .collect::<PyResult<Vec<_>>>()?
    } else {
        vec![(highlight.clone(), None)]
    };
    let defaults: &[HighlightStyle] =
        match entries.iter().filter(|(_, style)| style.is_none()).count() {
            0 => &[],
            1 => &[Amber],
            2 => &[ActedOn, Treated],
            3 => &[ActedOn, Middle, Treated],
            _ => {
                return Err(PyValueError::new_err(
                    "provide explicit highlight styles for all but at most three layers",
                ))
            }
        };
    let mut defaults = defaults.iter().copied();
    entries
        .into_iter()
        .map(|(entry, style)| {
            Ok(SelectedHighlight {
                selection: selection(&entry)?,
                style: style
                    .or_else(|| defaults.next())
                    .expect("one default per unstyled entry"),
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
    for paint in palette.iter().copied().chain(named) {
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
