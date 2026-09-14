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

    fs::write(
        base.path()
            .join(".clinnet/templates/curved-arrow-behavior.typ"),
        include_str!("resources/curved-arrow-behavior.typ"),
    )
    .unwrap();

    fs::write(
        base.path()
            .join(".clinnet/templates/weighted-cut-behavior.typ"),
        include_str!("resources/weighted-cut-behavior.typ"),
    )
    .unwrap();

    fs::write(
        base.path()
            .join(".clinnet/templates/named-map-behavior.typ"),
        include_str!("resources/named-map-behavior.typ"),
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

    for (fill_color, stroke_color, edge_color, node_color, unit) in [
        ("#f472b6", "#9f1239", "#b45309", "#1e3a8a", 10.0),
        ("#22d3ee", "#6d28d9", "#92400e", "#172554", 20.0),
    ] {
        let overlays = paths_with_attr(&svg, "fill", fill_color);
        let strokes = paths_with_attr(&svg, "stroke", stroke_color);
        let edges = paths_with_attr(&svg, "stroke", edge_color);
        let nodes = paths_with_attr(&svg, "fill", node_color);
        assert_eq!(overlays.len(), 1);
        assert_eq!(strokes.len(), 1);
        assert_eq!(nodes.len(), 2);
        assert!(!edges.is_empty());
        let (overlay_offset, overlay) = overlays[0];
        assert!(edges.iter().all(|(offset, _)| *offset < overlay_offset));
        assert!(nodes.iter().all(|(offset, _)| *offset < overlay_offset));
        assert!(overlay_offset < strokes[0].0);

        let edge_spans = merged_spans(edges.iter().map(|(_, tag)| horizontal_span(tag)).collect());
        assert_eq!(edge_spans.len(), 1);
        let span = horizontal_span(strokes[0].1);
        assert_close(span.1 - span.0, 4.0 * unit, 1e-3);
        assert_close(span.0, edge_spans[0].0, 1e-3);
        assert_close(span.1, edge_spans[0].1, 1e-3);
        let overlay_position = own_translation(overlay);
        let stroke_y = own_translation(strokes[0].1).1;
        assert_close(span.0 - overlay_position.0, unit, 1e-3);
        assert_close(stroke_y - overlay_position.1, unit, 1e-3);
        for (_, edge) in edges {
            assert_close(own_translation(edge).1, stroke_y, 1e-3);
        }
        for ((_, node), center_x) in nodes.iter().zip([span.0, span.1]) {
            let position = own_translation(node);
            assert_close(position.0 + 0.5 * unit, center_x, 1e-3);
            assert_close(position.1 + 0.5 * unit, stroke_y, 1e-3);
        }
    }

    let fixture = base
        .path()
        .join(".clinnet/templates/style-defaults-behavior.typ");
    let mut style_outputs = std::collections::BTreeMap::new();
    for mode in [
        "defaults",
        "draw-edge",
        "record-edge",
        "draw-endpoint",
        "record-endpoint",
        "half-endpoint",
        "delegate",
        "hidden",
        "decoration",
        "carrier-clean",
        "carrier-hostile",
        "carrier-removed",
        "ordinary",
        "carrier-source",
        "carrier-sink",
        "carrier-auto",
        "carrier-left",
        "carrier-auto-right",
        "carrier-right",
    ] {
        let output = base.path().join(format!("style-{mode}.svg"));
        fs::write(
            &fixture,
            include_str!("resources/style-defaults-behavior.typ")
                .replace("default: \"defaults\"", &format!("default: \"{mode}\"")),
        )
        .unwrap();
        renderer.compile_template(&fixture, &output, &[]).unwrap();
        style_outputs.insert(mode, fs::read_to_string(output).unwrap());
    }
    // Observe precedence in the final painted paths, through public graph.style/draw.
    for (mode, colors) in [
        ("defaults", vec!["#2563eb"]),
        ("draw-edge", vec!["#16a34a"]),
        ("record-edge", vec!["#ea580c"]),
        ("draw-endpoint", vec!["#0891b2", "#ea580c"]),
        ("record-endpoint", vec!["#9333ea", "#ea580c"]),
        ("half-endpoint", vec!["#dc2626", "#ea580c"]),
        ("hidden", vec![]),
        ("decoration", vec!["#16a34a", "#0891b2", "#9333ea"]),
    ] {
        let svg = &style_outputs[mode];
        for color in [
            "#2563eb", "#16a34a", "#ea580c", "#0891b2", "#9333ea", "#dc2626",
        ] {
            assert_eq!(
                !paths_with_attr(svg, "stroke", color).is_empty(),
                colors.contains(&color),
                "unexpected {color} in {mode}",
            );
        }
    }
    assert_eq!(style_outputs["defaults"], style_outputs["delegate"]);
    // Painting overrides and subgraph underlays cannot expose a label-only carrier.
    assert_eq!(
        style_outputs["carrier-clean"],
        style_outputs["carrier-hostile"]
    );
    // Unlaid-out curved graphs choose the same side for carrier and label.
    assert_eq!(style_outputs["carrier-auto"], style_outputs["carrier-left"]);
    assert_eq!(
        style_outputs["carrier-auto-right"],
        style_outputs["carrier-right"]
    );
    // Removing carriers restores the ordinary measured label for all three edge kinds.
    assert_eq!(style_outputs["carrier-removed"], style_outputs["ordinary"]);
    for mode in [
        "carrier-clean",
        "carrier-removed",
        "carrier-source",
        "carrier-sink",
    ] {
        let svg = &style_outputs[mode];
        assert_eq!(paths_with_attr(svg, "fill", "#ea580c").len(), 3, "{mode}");
        if mode == "carrier-source" || mode == "carrier-sink" {
            assert_eq!(paths_with_attr(svg, "stroke", "#2563eb").len(), 2, "{mode}");
        }
    }

    let templates = base.path().join(".clinnet/templates");
    fs::create_dir_all(templates.join("impl")).unwrap();
    fs::write(
        templates.join("gamma-physics-edge-style.typ"),
        include_str!("../../../assets/embedded/drawing/templates/physics-edge-style.typ"),
    )
    .unwrap();
    fs::write(
        templates.join("impl/physics-edge-style.typ"),
        include_str!("../../../assets/embedded/drawing/templates/impl/physics-edge-style.typ"),
    )
    .unwrap();
    fs::write(
        templates.join("gamma-layout-core.typ"),
        include_str!("../../../assets/embedded/drawing/templates/layout-core.typ"),
    )
    .unwrap();
    let fixture = templates.join("gamma-momentum-behavior.typ");
    let mut momentum_outputs = std::collections::BTreeMap::new();
    let mut label_positions = std::collections::BTreeMap::new();
    for mode in [
        "short",
        "long",
        "arrow-shift",
        "label-shift",
        "label-default",
        "anchor-center",
        "right",
        "no-mark",
        "reverse",
        "undirected",
        "incoming",
        "outgoing",
        "fallback",
        "ordinary",
        "curved-short",
        "curved-long",
        "default-mark",
    ] {
        let output = base.path().join(format!("momentum-{mode}.svg"));
        fs::write(
            &fixture,
            include_str!("resources/gamma-momentum-behavior.typ")
                .replace("default: \"short\"", &format!("default: \"{mode}\"")),
        )
        .unwrap();
        renderer.compile_template(&fixture, &output, &[]).unwrap();
        let svg = fs::read_to_string(output).unwrap();
        let labels = paths_with_attr(&svg, "fill", "#ea580c");
        // The native label replaces generated particle/q content and appears once.
        assert_eq!(labels.len(), 1, "{mode}");
        label_positions.insert(mode, preceding_group_translation(&svg, labels[0].0));
        assert_eq!(
            paths_with_attr(&svg, "fill", "#16a34a").len(),
            usize::from(mode != "undirected"),
            "{mode}"
        );
        assert_eq!(
            paths_with_attr(&svg, "fill", "#dc2626").len(),
            usize::from(!["no-mark", "fallback", "ordinary", "default-mark"].contains(&mode)),
            "{mode}"
        );
        momentum_outputs.insert(mode, svg);
    }
    // Test the GammaLoop callbacks themselves, independently of the Xbox example.
    for mode in [
        "long",
        "arrow-shift",
        "label-default",
        "no-mark",
        "reverse",
        "undirected",
    ] {
        assert_eq!(label_positions["short"], label_positions[mode], "{mode}");
    }
    assert_eq!(
        label_positions["curved-short"],
        label_positions["curved-long"]
    );
    assert_eq!(momentum_outputs["short"], momentum_outputs["label-default"]);
    assert_eq!(momentum_outputs["fallback"], momentum_outputs["ordinary"]);
    assert_close(
        label_positions["label-shift"].0 - label_positions["short"].0,
        -15.0,
        0.02,
    );
    for (mode, length, shift, side, gap) in [
        ("short", 6.0, 5.0, -1.0, 6.0),
        ("long", 24.0, 5.0, -1.0, 6.0),
        ("arrow-shift", 6.0, 10.0, -1.0, 6.0),
        ("label-shift", 6.0, 5.0, -1.0, 6.0),
        ("anchor-center", 6.0, 5.0, -1.0, 4.5),
        ("right", 6.0, 5.0, 1.0, 6.0),
        ("incoming", 6.0, 5.0, -1.0, 6.0),
        ("outgoing", 6.0, 5.0, -1.0, 6.0),
    ] {
        let svg = &momentum_outputs[mode];
        let reference = stroke_spans(svg, "#2563eb");
        let shafts = stroke_spans(svg, "#dc2626");
        assert_eq!(reference.len(), 1, "{mode}");
        assert_eq!(shafts.len(), 1, "{mode}");
        assert_close(shafts[0].1 - shafts[0].0, length, 0.02);
        assert_close(
            (shafts[0].0 + shafts[0].1 - reference[0].0 - reference[0].1) / 2.0,
            shift,
            0.02,
        );
        let shaft_y = own_translation(paths_with_attr(svg, "stroke", "#dc2626")[0].1).1;
        let reference_y = own_translation(paths_with_attr(svg, "stroke", "#2563eb")[0].1).1;
        assert_close(shaft_y - reference_y, side * 6.2, 0.02);
        assert_close(label_positions[mode].1 + 1.5 - shaft_y, side * gap, 0.02);
    }
    let forward = &momentum_outputs["short"];
    let reversed = &momentum_outputs["reverse"];
    let forward_fermion = paths_with_attr(forward, "fill", "#16a34a")[0].1;
    let reversed_fermion = paths_with_attr(reversed, "fill", "#16a34a")[0].1;
    // The triangle's first relative move locates its tip: right for forward,
    // left for reversed. Momentum remains source-to-sink in every orientation.
    assert!(svg_numbers(svg_attr(forward_fermion, "d").unwrap())[2] > 2.0);
    assert_close(
        svg_numbers(svg_attr(reversed_fermion, "d").unwrap())[2],
        0.0,
        1e-6,
    );
    for mode in ["reverse", "undirected", "incoming", "outgoing"] {
        assert_eq!(
            svg_attr(
                paths_with_attr(&momentum_outputs[mode], "fill", "#dc2626")[0].1,
                "d"
            ),
            svg_attr(paths_with_attr(forward, "fill", "#dc2626")[0].1, "d")
        );
    }
    let default_marks = paths_with_attr(&momentum_outputs["default-mark"], "stroke", "#dc2626");
    assert_eq!(default_marks.len(), 2);
    assert!(
        default_marks
            .iter()
            .all(|(_, tag)| svg_attr(tag, "stroke-width") == Some("0.4"))
    );

    let mut outside_positions = std::collections::BTreeMap::new();
    for mode in [
        "outside-incoming",
        "outside-incoming-ordinary",
        "outside-outgoing",
        "outside-outgoing-ordinary",
        "outside-unmatched",
        "outside-unmatched-ordinary",
        "outside-incoming-shift",
        "outside-outgoing-label-shift",
        "outside-outgoing-anchor",
    ] {
        let output = base.path().join(format!("momentum-{mode}.svg"));
        fs::write(
            &fixture,
            include_str!("resources/gamma-momentum-behavior.typ")
                .replace("default: \"short\"", &format!("default: \"{mode}\"")),
        )
        .unwrap();
        renderer.compile_template(&fixture, &output, &[]).unwrap();
        let svg = fs::read_to_string(output).unwrap();
        let labels = paths_with_attr(&svg, "fill", "#ea580c");
        assert_eq!(labels.len(), 1, "{mode}");
        let label_position = preceding_group_translation(&svg, labels[0].0);
        let reference = stroke_spans(&svg, "#2563eb")[0];
        let reference_y = own_translation(paths_with_attr(&svg, "stroke", "#2563eb")[0].1).1;
        let position = (
            label_position.0 - reference.0,
            label_position.1 - reference_y,
        );
        outside_positions.insert(mode, position);
        assert_eq!(
            paths_with_attr(&svg, "fill", "#dc2626").len(),
            usize::from(!mode.ends_with("-ordinary")),
            "{mode}"
        );
        if mode == "outside-incoming" {
            // The fixture's label is a rectangular box around actual particle/q content.
            let width = svg_numbers(svg_attr(labels[0].1, "d").unwrap())[3];
            assert_close(position.0 + width, -10.0, 0.02);
        } else if mode == "outside-outgoing" || mode == "outside-unmatched" {
            assert_close(label_position.0 - reference.1, 10.0, 0.02);
        }
    }
    // Ordered external labels retain the ordinary outside position when arrows
    // are enabled, including cut tags with no matching counterpart.
    for (mode, ordinary) in [
        ("outside-incoming", "outside-incoming-ordinary"),
        ("outside-outgoing", "outside-outgoing-ordinary"),
        ("outside-unmatched", "outside-unmatched-ordinary"),
    ] {
        assert_close(
            outside_positions[mode].0,
            outside_positions[ordinary].0,
            1e-6,
        );
        assert_close(
            outside_positions[mode].1,
            outside_positions[ordinary].1,
            1e-6,
        );
    }
    // Explicit momentum placement retains the path carrier. An arrow shift
    // still supplies the omitted label shift, even for prepared external legs.
    for mode in [
        "outside-incoming-shift",
        "outside-outgoing-label-shift",
        "outside-outgoing-anchor",
    ] {
        assert!(outside_positions[mode].1 < outside_positions["outside-outgoing"].1 - 5.0);
    }

    fs::write(
        templates.join("epemttbar.dot"),
        include_str!("../../../tests/resources/graphs/epemttbar.dot"),
    )
    .unwrap();
    let fixture = templates.join("gamma-cross-section-behavior.typ");
    fs::write(
        &fixture,
        include_str!("resources/gamma-cross-section-behavior.typ"),
    )
    .unwrap();
    renderer
        .compile_template(&fixture, base.path().join("gamma-cross-section.pdf"), &[])
        .unwrap();
    fs::write(
        templates.join("gamma-debug.dot"),
        r#"digraph {
            a [pos="0,0!"]; b [pos="4,0!"]; incoming [style=invis];
            incoming -> a [particle=fermion, pos="-2,0!"];
            a -> b [particle=fermion];
        }"#,
    )
    .unwrap();
    let fixture = templates.join("gamma-debug-behavior.typ");
    fs::write(&fixture, include_str!("resources/gamma-debug-behavior.typ")).unwrap();
    renderer
        .compile_template(&fixture, base.path().join("gamma-debug.pdf"), &[])
        .unwrap();
}

#[test]
fn public_weighted_cut_rejects_invalid_selections_and_stale_topology() {
    let configured_typst = std::env::var_os("TYPST_TEST_EXECUTABLE").map(PathBuf::from);
    let typst = configured_typst
        .clone()
        .unwrap_or_else(|| PathBuf::from("typst"));
    match Command::new(&typst).arg("--version").output() {
        Ok(version) if version.status.success() => {}
        result if configured_typst.is_some() => {
            panic!("configured Typst executable failed: {result:?}")
        }
        _ => return,
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
    let fixture = base.path().join(".clinnet/templates/invalid-cut.typ");
    let output = base.path().join("invalid-cut.svg");
    let prelude = r#"
#set page(width: auto, height: auto)
#import "crates/linnest/typst/src/lib.typ": graph, subgraph, layout, draw
#import graph: node, edge, source, sink
#import "map-style.typ": momentum
#let g = graph.build({
  node(<a>); node(<b>); node(<c>)
  edge(<e>, source(<a>), sink(<b>))
  edge(<f>, source(<b>), sink(<c>))
  edge(<out>, source(<c>))
})
#let left = subgraph.select(g, source: (<e>,))
#let right = subgraph.select(g, sink: (<e>,))
"#;
    // Prove the shared setup and installed Wasm work before accepting any failure.
    fs::write(
        &fixture,
        format!(
            "{prelude}\n#let opened = graph.cut(g, left: subgraph.with-data(g, left, hedge: (winding: 2)), right: right)\n#assert(graph.boundaries(opened).len() == 4)"
        ),
    )
    .unwrap();
    renderer.compile_template(&fixture, &output, &[]).unwrap();

    let mut cases = vec![
        (
            "overlap".to_owned(),
            "#let _ = graph.cut(g, left: left, right: left)".to_owned(),
            "disjoint, paired half-edges",
        ),
        (
            "not-inverse".to_owned(),
            "#let _ = graph.cut(g, left: left, right: subgraph.select(g, sink: (<f>,)))".to_owned(),
            "right must be the involution of left",
        ),
        (
            "dangling".to_owned(),
            "#let _ = graph.cut(g, left: subgraph.select(g, source: (<out>,)), right: right)".to_owned(),
            "select exactly one half of each paired edge on each side",
        ),
        (
            "mismatched-winding".to_owned(),
            "#let _ = graph.cut(g, left: subgraph.with-data(g, left, hedge: (winding: 2)), right: subgraph.with-data(g, right, hedge: (winding: 3)))".to_owned(),
            "winding must be a positive integer agreeing on both sides",
        ),
        (
            "unknown-name".to_owned(),
            "#let _ = subgraph.select(g, edges: (<missing>,))".to_owned(),
            "subgraph.select: unknown edge",
        ),
        (
            "unknown-hedge".to_owned(),
            "#let _ = subgraph.select(g, hedges: (99,))".to_owned(),
            "subgraph.select: unknown half-edge ID",
        ),
        (
            "unselected-annotation".to_owned(),
            "#let _ = subgraph.hedge-data(left, subgraph.hedges(right).first())".to_owned(),
            "subgraph.hedge-data: half-edge ID is not selected",
        ),
        (
            "raw-bytes".to_owned(),
            "#let _ = graph.edges(g, subgraph: left.bytes)".to_owned(),
            "subgraph: expected a Linnest subgraph object",
        ),
    ];
    for (options, expected) in [
        ("1", "momentum: expected named options"),
        ("gap: 0.2", "momentum: unknown arrow option gap"),
        ("label: none", "momentum: label must be a dictionary"),
        ("label: 1", "momentum: label must be a dictionary"),
        ("label: ()", "momentum: label must be a dictionary"),
        (
            "label: (length: 1)",
            "momentum: unknown label option length",
        ),
    ] {
        cases.push((
            format!("momentum-{options}"),
            format!("#let _ = momentum({options})"),
            expected,
        ));
    }
    // Name maps validate every supplied key, including explicit no-op entries.
    for (kind, key, missing) in [("node", "a", "A"), ("edge", "e", "E")] {
        for name in [missing, "0", "*", "default"] {
            cases.push((
                format!("map-{kind}-unknown-{name}"),
                format!("#let _ = graph.map(g, {kind}: (\"{name}\": none))"),
                "no record named",
            ));
        }
        for invalid in ["0", "false", "auto", "\"patch\"", "[patch]", "()"] {
            cases.push((
                format!("map-{kind}-invalid-mapper-{invalid}"),
                format!("#let _ = graph.map(graph.build(), {kind}: {invalid})"),
                "expected none, a callback, or a name-keyed dictionary",
            ));
            cases.push((
                format!("map-{kind}-invalid-entry-{invalid}"),
                format!("#let _ = graph.map(g, {kind}: ({key}: {invalid}))"),
                "must be none, a dictionary, or a callback",
            ));
            for mapper in [
                format!("_ => {invalid}"),
                format!("({key}: _ => {invalid})"),
            ] {
                cases.push((
                    format!("map-{kind}-invalid-return-{mapper}"),
                    format!("#let _ = graph.map(g, {kind}: {mapper})"),
                    "callback must return none or a dictionary",
                ));
            }
        }
        cases.push((
            format!("map-{kind}-unknown-on-empty-graph"),
            format!("#let _ = graph.map(graph.build(), {kind}: ({key}: none))"),
            "no record named",
        ));
    }
    for kind in ["graph", "source", "sink"] {
        cases.push((
            format!("map-{kind}-rejects-dictionary"),
            format!("#let _ = graph.map(g, {kind}: (:))"),
            "expected none or a callback",
        ));
    }
    for name in ["D1", "D1.00", "D1.3", "d1.0"] {
        cases.push((
            format!("map-unknown-cut-fragment-{name}"),
            format!(
                r#"#let master = graph.build({{
  node(<a>); node(<b>)
  edge(<D1>, source(<a>), sink(<b>))
}})
#let opened = graph.cut(master,
  left: subgraph.with-data(master, subgraph.select(master, source: (<D1>,)), hedge: (winding: 2)),
  right: subgraph.select(master, sink: (<D1>,)),
)
#let _ = graph.map(opened, edge: ("{name}": none))"#
            ),
            "graph.map edge: no record named",
        ));
    }
    for winding in ["0", "-1", "1.5", "\"2\""] {
        cases.push((
            format!("invalid-winding-{winding}"),
            format!(
                "#let _ = graph.cut(g, left: subgraph.with-data(g, left, hedge: (winding: {winding})), right: right)"
            ),
            "winding must be a positive integer agreeing on both sides",
        ));
    }
    for (left_winding, right_winding) in [("2", "2.0"), ("2.0", "2")] {
        cases.push((
            format!("mixed-winding-types-{left_winding}-{right_winding}"),
            format!(
                "#let _ = graph.cut(g, left: subgraph.with-data(g, left, hedge: (winding: {left_winding})), right: subgraph.with-data(g, right, hedge: (winding: {right_winding})))"
            ),
            "winding must be a positive integer agreeing on both sides",
        ));
    }
    // All-false masks have no out-of-range selected IDs to expose a size mismatch.
    for size in [4, 6] {
        cases.push((
            format!("wrong-mask-size-{size}-for-5-hedges"),
            format!(
                r#"#let other = graph.build({{
  node(<other>)
  for _ in range({size}) {{ edge(source(<other>)) }}
}})
#let malformed = left + (bytes: subgraph.select(other).bytes)
#assert(malformed.topology == left.topology and subgraph.hedges(malformed) == ())
#let _ = graph.edges(g, subgraph: malformed)"#
            ),
            "subgraph: mask size does not match graph",
        ));
    }
    // Identical sizes, names, and hedge IDs are insufficient: e's source moved.
    let changed_source = r#"
#let g = graph.build({
  node(<a>); node(<b>); node(<c>)
  edge(<e>, source(<c>), sink(<b>))
  edge(<f>, source(<b>), sink(<c>))
  edge(<out>, source(<c>))
})
"#;
    for consumer in [
        "#let _ = graph.nodes(g, subgraph: left)",
        "#let _ = graph.edges(g, subgraph: left)",
        "#let _ = layout(g, subgraph: left)",
        "#draw(g, subgraph: left)",
        "#let _ = subgraph.with-data(g, left)",
        "#let _ = subgraph.complement(g, left)",
        "#let _ = graph.cut(g, left: left, right: right)",
    ] {
        cases.push((
            format!("stale-source: {consumer}"),
            format!("{changed_source}\n{consumer}"),
            "subgraph: topology does not match graph",
        ));
    }
    for (name, body, expected) in cases {
        fs::write(&fixture, format!("{prelude}\n{body}")).unwrap();
        let error = renderer
            .compile_template(&fixture, &output, &[])
            .expect_err(&format!("{name} unexpectedly compiled"));
        assert!(
            error.to_string().contains(expected),
            "{name}: expected diagnostic {expected:?}, got {error}"
        );
    }
}

#[test]
fn public_gammaloop_momentum_geometry_ignores_label_visibility() {
    let configured_typst = std::env::var_os("TYPST_TEST_EXECUTABLE").map(PathBuf::from);
    let typst = configured_typst
        .clone()
        .unwrap_or_else(|| PathBuf::from("typst"));
    match Command::new(&typst).arg("--version").output() {
        Ok(version) if version.status.success() => {
            let version = String::from_utf8_lossy(&version.stdout);
            if !version.starts_with("typst 0.15.") {
                assert!(
                    configured_typst.is_none(),
                    "SVG fixture requires Typst 0.15.x, found {version}"
                );
                return;
            }
        }
        result if configured_typst.is_some() => {
            panic!("configured Typst executable failed: {result:?}")
        }
        _ => return,
    }
    let base = tempfile::tempdir().unwrap();
    let renderer = TypstRenderer::new(base.path()).typst_executable(typst);
    renderer.check_version().unwrap();
    renderer.stage_default_assets().unwrap();
    let templates = base.path().join(".clinnet/templates");
    fs::create_dir_all(templates.join("impl")).unwrap();
    fs::write(
        templates.join("gamma-physics-edge-style.typ"),
        include_str!("../../../assets/embedded/drawing/templates/physics-edge-style.typ"),
    )
    .unwrap();
    fs::write(
        templates.join("impl/physics-edge-style.typ"),
        include_str!("../../../assets/embedded/drawing/templates/impl/physics-edge-style.typ"),
    )
    .unwrap();
    let fixture = templates.join("gamma-momentum-side-behavior.typ");
    let mut baseline: Option<Vec<(String, (f64, f64))>> = None;
    for mode in [
        "full-above",
        "full-below",
        "particle-above",
        "particle-below",
        "none-above",
        "none-below",
    ] {
        fs::write(
            &fixture,
            include_str!("resources/gamma-momentum-side-behavior.typ")
                .replace("default: \"full-above\"", &format!("default: \"{mode}\"")),
        )
        .unwrap();
        let output = base.path().join(format!("momentum-side-{mode}.svg"));
        renderer.compile_template(&fixture, &output, &[]).unwrap();
        let svg = fs::read_to_string(output).unwrap();
        let particle = paths_with_attr(&svg, "stroke", "#2563eb");
        assert!(!particle.is_empty(), "{mode}");
        assert_eq!(paths_with_attr(&svg, "fill", "#dc2626").len(), 9, "{mode}");
        let origin = own_translation(particle[0].1);
        // Text changes the canvas bounds. Compare actual shafts, arrowheads,
        // and particle paths relative to the first underlying particle path.
        let geometry = svg_paths(&svg)
            .into_iter()
            .filter(|(_, tag)| {
                matches!(svg_attr(tag, "stroke"), Some("#2563eb" | "#dc2626"))
                    || svg_attr(tag, "fill") == Some("#dc2626")
            })
            .map(|(_, tag)| {
                let position = own_translation(tag);
                (
                    svg_attr(tag, "d").unwrap().to_owned(),
                    (position.0 - origin.0, position.1 - origin.1),
                )
            })
            .collect::<Vec<_>>();
        if let Some(expected) = &baseline {
            assert_eq!(geometry.len(), expected.len(), "{mode}");
            for (actual, expected) in geometry.iter().zip(expected) {
                assert_eq!(actual.0, expected.0, "{mode}");
                assert_close(actual.1.0, expected.1.0, 1e-4);
                assert_close(actual.1.1, expected.1.1, 1e-4);
            }
        } else {
            baseline = Some(geometry);
        }
    }
}

#[test]
fn public_gammaloop_external_order_tracks_physical_indices() {
    let configured_typst = std::env::var_os("TYPST_TEST_EXECUTABLE").map(PathBuf::from);
    let typst = configured_typst
        .clone()
        .unwrap_or_else(|| PathBuf::from("typst"));
    match Command::new(&typst).arg("--version").output() {
        Ok(version) if version.status.success() => {}
        result if configured_typst.is_some() => {
            panic!("configured Typst executable failed: {result:?}")
        }
        _ => return,
    }
    let base = tempfile::tempdir().unwrap();
    let renderer = TypstRenderer::new(base.path()).typst_executable(typst);
    renderer.check_version().unwrap();
    renderer.stage_default_assets().unwrap();
    let templates = base.path().join(".clinnet/templates");
    fs::create_dir_all(templates.join("impl")).unwrap();
    for (name, source) in [
        (
            "gamma-physics-edge-style.typ",
            include_str!("../../../assets/embedded/drawing/templates/physics-edge-style.typ"),
        ),
        (
            "impl/physics-edge-style.typ",
            include_str!("../../../assets/embedded/drawing/templates/impl/physics-edge-style.typ"),
        ),
        (
            "gamma-layout-core.typ",
            include_str!("../../../assets/embedded/drawing/templates/layout-core.typ"),
        ),
        (
            "gamma-external-order-behavior.typ",
            include_str!("resources/gamma-external-order-behavior.typ"),
        ),
    ] {
        fs::write(templates.join(name), source).unwrap();
    }
    renderer
        .compile_template(
            templates.join("gamma-external-order-behavior.typ"),
            base.path().join("gamma-external-order.pdf"),
            &[],
        )
        .unwrap();
}
