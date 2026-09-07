use std::fs;
use std::path::PathBuf;
use std::process::Command;

use clinnet::TypstRenderer;

fn svg_paths(svg: &str) -> Vec<(usize, &str)> {
    let mut paths = Vec::new();
    let mut cursor = 0;
    while let Some(relative) = svg[cursor..].find("<path ") {
        let start = cursor + relative;
        let end = start + svg[start..].find("/>").unwrap() + 2;
        paths.push((start, &svg[start..end]));
        cursor = end;
    }
    paths
}

fn svg_attr<'a>(tag: &'a str, name: &str) -> Option<&'a str> {
    let needle = format!(r#"{name}=""#);
    let start = tag.find(&needle)? + needle.len();
    let rest = &tag[start..];
    Some(&rest[..rest.find('"')?])
}

fn paths_with_attr<'a>(svg: &'a str, name: &str, value: &str) -> Vec<(usize, &'a str)> {
    svg_paths(svg)
        .into_iter()
        .filter(|(_, tag)| svg_attr(tag, name) == Some(value))
        .collect()
}

fn translation(value: &str) -> (f64, f64) {
    let mut coordinates = value
        .strip_prefix("translate(")
        .unwrap()
        .split([' ', ',', ')'])
        .filter(|part| !part.is_empty())
        .map(|part| part.parse().unwrap());
    (
        coordinates.next().unwrap(),
        coordinates.next().unwrap_or(0.0),
    )
}

fn own_translation(tag: &str) -> (f64, f64) {
    svg_attr(tag, "transform").map_or((0.0, 0.0), translation)
}

fn preceding_group_translation(svg: &str, before: usize) -> (f64, f64) {
    let key = r#"<g transform=""#;
    let start = svg[..before].rfind(key).unwrap() + key.len();
    let rest = &svg[start..];
    translation(&rest[..rest.find('"').unwrap()])
}

fn svg_numbers(value: &str) -> Vec<f64> {
    value
        .split(|character: char| {
            !(character.is_ascii_digit() || matches!(character, '.' | '-' | '+' | 'e' | 'E'))
        })
        .filter_map(|part| part.parse().ok())
        .collect()
}

// The fixture deliberately uses horizontal cubic paths, so inspecting their
// final x delta is enough and avoids coupling this test to a general SVG parser.
fn horizontal_span(tag: &str) -> (f64, f64) {
    let values = svg_numbers(svg_attr(tag, "d").unwrap());
    assert!(values[0].abs() < 1e-9 && values[1].abs() < 1e-9);
    let start = own_translation(tag).0;
    let end = start + values[values.len() - 2];
    (start.min(end), start.max(end))
}

fn merged_spans(mut spans: Vec<(f64, f64)>) -> Vec<(f64, f64)> {
    spans.sort_by(|left, right| left.0.total_cmp(&right.0));
    let mut merged: Vec<(f64, f64)> = Vec::new();
    for span in spans {
        if let Some(last) = merged.last_mut().filter(|last| span.0 <= last.1 + 1e-3) {
            last.1 = last.1.max(span.1);
        } else {
            merged.push(span);
        }
    }
    merged
}

fn stroke_spans(svg: &str, color: &str) -> Vec<(f64, f64)> {
    merged_spans(
        paths_with_attr(svg, "stroke", color)
            .into_iter()
            .map(|(_, tag)| horizontal_span(tag))
            .collect(),
    )
}

fn assert_close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "expected {expected} ± {tolerance}, got {actual}"
    );
}

#[test]
fn public_linnest_layout_and_drawing_behavior_is_observable() {
    let configured_typst = std::env::var_os("TYPST_TEST_EXECUTABLE").map(PathBuf::from);
    let typst = configured_typst
        .clone()
        .unwrap_or_else(|| PathBuf::from("typst"));
    let version = match Command::new(&typst).arg("--version").output() {
        Ok(version) => version,
        Err(error) if configured_typst.is_some() => {
            panic!("configured Typst executable failed: {error}")
        }
        Err(_) => return,
    };
    if !version.status.success() {
        if configured_typst.is_some() {
            panic!(
                "configured Typst executable returned {status}",
                status = version.status
            );
        }
        return;
    }
    let version = String::from_utf8_lossy(&version.stdout);
    let probe_svg = version.starts_with("typst 0.15.");
    if configured_typst.is_some() && !probe_svg {
        panic!("SVG behavior fixture requires Typst 0.15.x, found {version}");
    }

    let base = tempfile::tempdir().unwrap();
    let renderer = TypstRenderer::new(base.path()).typst_executable(typst);
    renderer.check_version().unwrap();
    renderer.stage_default_assets().unwrap();
    fs::write(
        base.path().join(".clinnet/templates/map-style.typ"),
        include_str!("../../linnest/typst/examples/map-style.typ"),
    )
    .unwrap();

    let fixture = base
        .path()
        .join(".clinnet/templates/linnest-public-behavior.typ");
    fs::write(
        &fixture,
        include_str!("resources/linnest-public-behavior.typ"),
    )
    .unwrap();
    let output = base.path().join("linnest-public-behavior.svg");
    renderer.compile_template(&fixture, &output, &[]).unwrap();
    let svg = fs::read_to_string(output).unwrap();
    // The Typst assertions above remain useful on newer versions. The SVG
    // serializer probes below are intentionally pinned to Typst 0.15.
    if !probe_svg {
        return;
    }

    let shaft = stroke_spans(&svg, "#16a34a");
    let shaft_reference = stroke_spans(&svg, "#0ea5e9");
    let arrow_head = paths_with_attr(&svg, "fill", "#d119e6");
    assert_eq!(shaft.len(), 1);
    assert_eq!(shaft_reference.len(), 1);
    assert_eq!(arrow_head.len(), 1);
    let shaft_span = shaft[0];
    let shaft_reference_span = shaft_reference[0];
    assert_close(shaft_span.0 - shaft_reference_span.0, 10.0, 1e-3);
    assert_close(shaft_span.1, own_translation(arrow_head[0].1).0, 1e-3);

    let crossing = stroke_spans(&svg, "#dc2626");
    assert_eq!(crossing.len(), 2);
    let crossing_gap = crossing[1].0 - crossing[0].1;
    let crossing_span = crossing[1].1 - crossing[0].0;
    assert_close(crossing_gap / crossing_span, 1.0 / 6.0, 5e-4);
    let gap_center = (crossing[0].1 + crossing[1].0) / 2.0;
    assert_close(
        (gap_center - crossing[0].0) / crossing_span,
        2.0 / 6.0,
        5e-4,
    );
    assert!(!stroke_spans(&svg, "#2563eb").is_empty());
    let crossing_mark = paths_with_attr(&svg, "fill", "#7e22ce");
    let crossing_reference_mark = paths_with_attr(&svg, "fill", "#9333ea");
    assert_eq!(crossing_mark.len(), 1);
    assert_eq!(crossing_reference_mark.len(), 1);
    assert_close(
        own_translation(crossing_mark[0].1).0,
        own_translation(crossing_reference_mark[0].1).0,
        1e-3,
    );

    let left_label = paths_with_attr(&svg, "fill", "#ea580c");
    let right_label = paths_with_attr(&svg, "fill", "#0891b2");
    assert_eq!(left_label.len(), 1);
    assert_eq!(right_label.len(), 1);
    let label_reference = stroke_spans(&svg, "#475569");
    let left_carrier_paths = paths_with_attr(&svg, "stroke", "#c2410c");
    let right_carrier_paths = paths_with_attr(&svg, "stroke", "#0e7490");
    let left_carrier = stroke_spans(&svg, "#c2410c");
    let right_carrier = stroke_spans(&svg, "#0e7490");
    assert_eq!(label_reference.len(), 1);
    assert_eq!(left_carrier.len(), 1);
    assert_eq!(right_carrier.len(), 1);
    let label_span = label_reference[0].1 - label_reference[0].0;
    let left_center = (left_carrier[0].0 + left_carrier[0].1) / 2.0;
    let right_center = (right_carrier[0].0 + right_carrier[0].1) / 2.0;
    assert_close((right_center - left_center) / label_span, 2.0 / 6.0, 5e-4);

    let left_label_position = preceding_group_translation(&svg, left_label[0].0);
    let right_label_position = preceding_group_translation(&svg, right_label[0].0);
    assert_close(left_label_position.0 + 1.5, left_center, 0.02);
    assert_close(right_label_position.0 + 1.5, right_center, 0.02);
    assert_close(
        own_translation(left_carrier_paths[0].1).1 - (left_label_position.1 + 1.5),
        2.5,
        0.1,
    );
    assert_close(
        own_translation(right_carrier_paths[0].1).1 - (right_label_position.1 + 1.5),
        2.5,
        0.1,
    );

    let momentum_reference = paths_with_attr(&svg, "stroke", "#78716c");
    let reference_y = own_translation(momentum_reference[0].1).1;
    let reference_spans = stroke_spans(&svg, "#78716c");
    assert_eq!(reference_spans.len(), 1);
    let reference_center = (reference_spans[0].0 + reference_spans[0].1) / 2.0;
    for (shaft_color, label_color, side, label_shift) in [
        ("#86198f", "#a21caf", 1.0, 5.0),
        ("#075985", "#0369a1", -1.0, -10.0),
    ] {
        let shafts = paths_with_attr(&svg, "stroke", shaft_color);
        let labels = paths_with_attr(&svg, "fill", label_color);
        assert_eq!(stroke_spans(&svg, shaft_color).len(), 1);
        assert_eq!(labels.len(), 1);
        let shaft_y = own_translation(shafts[0].1).1;
        let label_position = preceding_group_translation(&svg, labels[0].0);
        assert_close(reference_y - shaft_y, side * 8.0, 0.02);
        assert_close(shaft_y - (label_position.1 + 1.5), side * 6.0, 0.02);
        assert_close(label_position.0 + 1.5, reference_center + label_shift, 0.02);
    }

    for color in ["#f97316", "#65a30d", "#0f766e"] {
        let spans = stroke_spans(&svg, color);
        assert_eq!(spans.len(), 1, "expected one visible {color} edge");
        assert_close(spans[0].1 - spans[0].0, 40.0, 0.02);
    }
    for color in ["#84cc16", "#ca8a04", "#4d7c0f", "#be123c"] {
        assert!(
            paths_with_attr(&svg, "stroke", color).is_empty(),
            "unexpected {color} fallback or overridden style"
        );
    }
}
