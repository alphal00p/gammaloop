use itertools::Itertools;
use linnet::half_edge::involution::{EdgeIndex, Orientation};
use spenso::algebra::algebraic_traits::IsZero;
use spenso::algebra::complex::Complex;
use symbolica::numerical_integration::Sample;

use crate::{
    integrate::discrete_axis_label,
    utils::{self, F},
};

#[derive(Clone, Copy, Debug, Eq, PartialEq, Default)]
pub(crate) struct TextStyle {
    pub(crate) color: Option<TextColor>,
    pub(crate) bold: bool,
    pub(crate) dimmed: bool,
}

impl TextStyle {
    pub(crate) const PLAIN: Self = Self {
        color: None,
        bold: false,
        dimmed: false,
    };

    pub(crate) const fn with_color(mut self, color: TextColor) -> Self {
        self.color = Some(color);
        self
    }

    pub(crate) const fn bold(mut self) -> Self {
        self.bold = true;
        self
    }

    pub(crate) const fn dimmed(mut self) -> Self {
        self.dimmed = true;
        self
    }

    pub(crate) const fn green() -> Self {
        Self::PLAIN.with_color(TextColor::Green)
    }

    pub(crate) const fn blue() -> Self {
        Self::PLAIN.with_color(TextColor::Blue)
    }

    pub(crate) const fn red() -> Self {
        Self::PLAIN.with_color(TextColor::Red)
    }

    pub(crate) const fn yellow() -> Self {
        Self::PLAIN.with_color(TextColor::Yellow)
    }

    pub(crate) const fn pink() -> Self {
        Self::PLAIN.with_color(TextColor::Pink)
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) enum TextColor {
    Green,
    Blue,
    Red,
    Yellow,
    Pink,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub(crate) struct TextSpan {
    pub(crate) text: String,
    pub(crate) style: TextStyle,
}

impl TextSpan {
    pub(crate) fn new(text: impl Into<String>, style: TextStyle) -> Self {
        Self {
            text: text.into(),
            style,
        }
    }
}

#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub(crate) struct StyledText {
    pub(crate) spans: Vec<TextSpan>,
}

impl StyledText {
    pub(crate) fn new() -> Self {
        Self::default()
    }

    pub(crate) fn plain(text: impl Into<String>) -> Self {
        Self::styled(text, TextStyle::PLAIN)
    }

    pub(crate) fn styled(text: impl Into<String>, style: TextStyle) -> Self {
        Self {
            spans: vec![TextSpan::new(text, style)],
        }
    }

    pub(crate) fn push_span(&mut self, span: TextSpan) {
        self.spans.push(span);
    }

    pub(crate) fn push_text(&mut self, text: impl Into<String>, style: TextStyle) {
        self.push_span(TextSpan::new(text, style));
    }

    pub(crate) fn append(&mut self, other: Self) {
        self.spans.extend(other.spans);
    }

    pub(crate) fn to_plain_string(&self) -> String {
        self.spans.iter().map(|span| span.text.as_str()).collect()
    }

    pub(crate) fn is_empty(&self) -> bool {
        self.spans.iter().all(|span| span.text.is_empty())
    }
}

impl From<String> for StyledText {
    fn from(value: String) -> Self {
        Self::plain(value)
    }
}

impl From<&str> for StyledText {
    fn from(value: &str) -> Self {
        Self::plain(value)
    }
}

#[derive(Clone, Debug)]
pub(crate) struct DisplayField<T> {
    pub(crate) raw: T,
    pub(crate) display: StyledText,
}

impl<T> DisplayField<T> {
    pub(crate) fn new(raw: T, display: StyledText) -> Self {
        Self { raw, display }
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) enum UncertaintyNotation {
    Dynamic,
    Scientific,
}

fn format_max_eval_coordinate(value: F<f64>) -> String {
    let formatted = format!("{:.16e}", value.0);
    let Some((mantissa, exponent)) = formatted.rsplit_once('e') else {
        return formatted;
    };
    let exponent = exponent.parse::<i32>().unwrap_or_default();
    format!("{mantissa}e{exponent:+03}")
}

fn format_max_eval_coordinates(xs: &[F<f64>]) -> String {
    if xs.len() <= 3 {
        return format!(
            "[ {} ]",
            xs.iter()
                .map(|value| format_max_eval_coordinate(*value))
                .join(" ")
        );
    }

    let rows = xs
        .chunks(3)
        .map(|chunk| {
            chunk
                .iter()
                .map(|value| format_max_eval_coordinate(*value))
                .join(" ")
        })
        .join("\n");
    format!("[\n{rows} ]")
}

fn append_max_eval_sample_parts(
    sample: &Sample<F<f64>>,
    axis_labels: &[String],
    depth: usize,
    parts: &mut Vec<String>,
) {
    match sample {
        Sample::Continuous(_, xs) => {
            parts.push(format!("xs: {}", format_max_eval_coordinates(xs)));
        }
        Sample::Discrete(_, index, Some(nested_sample)) => {
            parts.push(format!(
                "{}: {}",
                discrete_axis_label(axis_labels, depth),
                index
            ));
            append_max_eval_sample_parts(nested_sample, axis_labels, depth + 1, parts);
        }
        Sample::Discrete(_, index, None) => {
            parts.push(format!(
                "{}: {}",
                discrete_axis_label(axis_labels, depth),
                index
            ));
        }
        Sample::Uniform(_, indices, xs) => {
            for (offset, index) in indices.iter().enumerate() {
                parts.push(format!(
                    "{}: {}",
                    discrete_axis_label(axis_labels, depth + offset),
                    index
                ));
            }
            parts.push(format!("xs: {}", format_max_eval_coordinates(xs)));
        }
    }
}

pub(crate) fn format_max_eval_sample(
    sample: &Sample<F<f64>>,
    axis_labels: &[String],
    prefix_path: &[usize],
) -> String {
    let mut parts = prefix_path
        .iter()
        .enumerate()
        .map(|(depth, index)| format!("{}: {}", discrete_axis_label(axis_labels, depth), index))
        .collect_vec();
    append_max_eval_sample_parts(sample, axis_labels, prefix_path.len(), &mut parts);
    if parts.is_empty() {
        String::from("N/A")
    } else {
        parts.join(", ")
    }
}

pub(crate) fn format_significant_percentage(
    value: f64,
    significant_digits: usize,
    scientific_threshold: Option<(f64, f64)>,
) -> String {
    if !value.is_finite() {
        return "None".to_string();
    }

    if value == 0.0 {
        return format!("{:.*}%", significant_digits.saturating_sub(1), 0.0);
    }

    let abs_value = value.abs();
    if scientific_threshold.is_some_and(|(upper, lower)| abs_value >= upper || abs_value < lower) {
        return format!("{:.*e}%", significant_digits.saturating_sub(1), value);
    }

    let exponent = abs_value.log10().floor() as i32;
    let decimals = (significant_digits as i32 - exponent - 1).max(0) as usize;
    format!("{value:.decimals$}%")
}

fn should_use_scientific_uncertainty_notation(avg: F<f64>, err: F<f64>) -> bool {
    let avg_abs = avg.abs().0;
    let err_abs = err.abs().0;
    avg_abs >= 1e6
        || (avg_abs != 0.0 && avg_abs < 1e-5)
        || (avg.is_zero() && !(1e-4..1e5).contains(&err_abs))
}

pub(crate) fn format_uncertainty_with_notation(
    avg: F<f64>,
    err: F<f64>,
    notation: UncertaintyNotation,
) -> String {
    if !matches!(notation, UncertaintyNotation::Scientific)
        && !should_use_scientific_uncertainty_notation(avg, err)
    {
        return utils::format_uncertainty(avg, err);
    }

    if !avg.0.is_finite() || !err.0.is_finite() {
        return utils::format_uncertainty(avg, err);
    }

    let exponent = if avg.is_non_zero() {
        avg.abs().0.log10().floor() as i32
    } else if err.is_non_zero() {
        err.abs().0.log10().floor() as i32
    } else {
        0
    };
    let scale = 10_f64.powi(exponent);
    let mantissa = utils::format_uncertainty(F(avg.0 / scale), F(err.0 / scale));
    format!("{mantissa}e{exponent}")
}

pub(crate) fn format_signed_uncertainty(
    avg: F<f64>,
    err: F<f64>,
    notation: UncertaintyNotation,
) -> String {
    let formatted = format_uncertainty_with_notation(avg, err, notation);
    if avg.0.is_sign_negative() {
        formatted
    } else {
        format!("+{formatted}")
    }
}

pub(crate) fn format_abbreviated_count(value: usize) -> String {
    if value < 1_000 {
        return value.to_string();
    }

    let value = value as f64;
    if value < 1_000_000.0 {
        return format!("{:.2}K", value / 1_000.0);
    }
    if value < 1_000_000_000.0 {
        return format!("{:.2}M", value / 1_000_000.0);
    }
    if value < 1_000_000_000_000.0 {
        return format!("{:.2}B", value / 1_000_000_000.0);
    }

    format!("{:.3}T", value / 1_000_000_000_000.0)
}

pub(crate) fn format_iteration_points(points: usize) -> String {
    format_abbreviated_count(points)
}

pub(crate) fn format_total_points(points: usize) -> String {
    format_abbreviated_count(points)
}

pub(crate) fn graph_group_description<I>(graph_names: I) -> String
where
    I: IntoIterator<Item = String>,
{
    let graph_names = graph_names.into_iter().collect_vec();
    if graph_names.len() <= 1 {
        return graph_names.into_iter().next().unwrap_or_default();
    }

    format!("[{}]", graph_names.join(","))
}

pub(crate) fn orientation_description(orientation: &[Orientation]) -> String {
    orientation
        .iter()
        .map(|orientation| match *orientation {
            Orientation::Default => '+',
            Orientation::Reversed => '-',
            Orientation::Undirected => '0',
        })
        .collect()
}

pub(crate) fn lmb_channel_description(loop_edges: &[EdgeIndex]) -> String {
    format!(
        "({})",
        loop_edges
            .iter()
            .map(|edge_id| edge_id.0.to_string())
            .join(",")
    )
}

pub(crate) fn styled_orientation_description(description: &str) -> StyledText {
    let mut styled = StyledText::new();
    for sign in description.chars() {
        match sign {
            '+' => styled.push_text("+", TextStyle::green().bold()),
            '-' => styled.push_text("-", TextStyle::red().bold()),
            '0' => styled.push_text("0", TextStyle::PLAIN.dimmed()),
            other => styled.push_text(other.to_string(), TextStyle::PLAIN),
        }
    }
    styled
}

pub(crate) fn styled_graph_description(description: &str) -> StyledText {
    if let Some(inner) = description
        .strip_prefix('[')
        .and_then(|trimmed| trimmed.strip_suffix(']'))
    {
        let names = inner
            .split(',')
            .filter(|name| !name.is_empty())
            .collect_vec();
        if names.is_empty() {
            return StyledText::styled(description, TextStyle::green().bold());
        }

        let mut styled = StyledText::plain("[");
        for (index, name) in names.into_iter().enumerate() {
            if index > 0 {
                styled.push_text(",", TextStyle::PLAIN);
            }
            let style = if index == 0 {
                TextStyle::green().bold()
            } else {
                TextStyle::blue().bold()
            };
            styled.push_text(name, style);
        }
        styled.push_text("]", TextStyle::PLAIN);
        return styled;
    }

    StyledText::styled(description, TextStyle::green().bold())
}

pub(crate) fn styled_bin_description(axis_label: &str, description: &str) -> StyledText {
    match axis_label {
        "orientation" => styled_orientation_description(description),
        "graph" => styled_graph_description(description),
        _ => StyledText::styled(description, TextStyle::green().bold()),
    }
}

pub(crate) fn styled_signed_complex(value: Complex<F<f64>>, error: Complex<F<f64>>) -> StyledText {
    let rendered = format!(
        "( {}, {} )",
        format_signed_uncertainty(value.re, error.re, UncertaintyNotation::Scientific),
        format_signed_uncertainty(value.im, error.im, UncertaintyNotation::Scientific)
    );
    StyledText::styled(rendered, TextStyle::blue().bold())
}
