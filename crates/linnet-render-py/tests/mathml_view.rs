//! Native receipt regression tests and an opt-in, finite browser producer.
//!
//! Root coordinates compilation/execution. The ignored producer retains every
//! prepared view until its own browser receipt has completed native rendering.
use std::ffi::CString;

use linnest::svg::{
    Dash, Details, EdgeDrawing, LabelMetrics, LabelSource, NodeDrawing, Scene, Stroke,
};
use linnet_render_py::PyDiagramView;
use pyo3::exceptions::PyValueError;
use pyo3::types::{PyAnyMethods, PyDict, PyDictMethods, PyModule, PyModuleMethods};
use pyo3::{pyfunction, wrap_pyfunction, Py, PyResult, Python};
use roxmltree::Document;
use serde_json::{json, Value};

fn fixture_scene(marker: &str) -> Scene {
    // Scale the certified EC embedding uniformly instead of constraining its
    // nodes independently of routed carriers. A four-unit median route is
    // 54pt, leaving label room with zero solver steps and no label feedback.
    let node = |index| json!({"index": index, "name": format!("vertex-{index}")});
    let stroke = Stroke {
        paint: "#246813".into(),
        width: 0.7,
        dash: Dash::Solid,
        round_cap: true,
    };
    let labels = [
        ("<mmultiscripts><mi>T</mi><mi>ν</mi><mi>μ</mi><mprescripts/><none/><mi>ρ</mi></mmultiscripts>", "$attach(T, tl: rho, br: nu, tr: mu)$"),
        ("<mrow><msubsup><mi>α</mi><mi>i</mi><mn>2</mn></msubsup><mo>+</mo><mi>β</mi></mrow>", "$alpha_i^2 + beta$"),
        ("<mrow><mo>⟨</mo><msub><mi>V</mi><mi>μ</mi></msub><mo>⊗</mo><msup><mi>W</mi><mi>ν</mi></msup><mo>⟩</mo></mrow>", "$⟨ V_mu ⊗ W^nu ⟩$"),
        ("<mfrac><msup><mi>p</mi><mn>2</mn></msup><mrow><mi>m</mi><mo>+</mo><mi>ε</mi></mrow></mfrac>", "$frac(p^2, m + epsilon)$"),
        ("<msub><mi>ℓ</mi><mn>0</mn></msub>", "$ell_0$"),
    ];
    Scene {
        graph: serde_json::from_value(json!({
            "name": format!("mathml-{marker}"),
            "nodes": [node(0), node(1), node(2)],
            "edges": [
                {"id": 0, "source": {"node": 0, "id": 0}, "sink": {"node": 0, "id": 1}},
                {"id": 1, "source": {"node": 0, "id": 2}, "sink": {"node": 1, "id": 3}},
                {"id": 2, "source": {"node": 0, "id": 4}, "sink": {"node": 1, "id": 5}},
                {"id": 3, "sink": {"node": 0, "id": 6}},
                {"id": 4, "source": {"node": 1, "id": 7}}
            ]
        }))
        .unwrap(),
        nodes: (0..3)
            .map(|index| NodeDrawing {
                radius: 0.2,
                label_padding: Some(0.3),
                label: Some(index),
                rectangular: true,
                fill: "#edf3e7".into(),
                stroke: stroke.clone(),
                details: [
                    ("fixture", marker.to_owned()),
                    ("vertex", index.to_string()),
                ]
                .into_iter()
                .collect(),
            })
            .collect(),
        edges: (0..5)
            .map(|index| EdgeDrawing {
                stroke: Stroke {
                    dash: if index == 2 {
                        Dash::Dashed(3.0, 2.0)
                    } else {
                        Dash::Solid
                    },
                    ..stroke.clone()
                },
                pattern: None,
                flow: Some(index != 2),
                momentum: false,
                momentum_offset: EdgeDrawing::DEFAULT_MOMENTUM_OFFSET,
                label: Some(if index == 0 { 4 } else { 3 }),
                details: [
                    ("fixture", marker.to_owned()),
                    ("particle", format!("p{index}")),
                ]
                .into_iter()
                .collect::<Details>(),
            })
            .collect(),
        preamble: "#set text(size: 9pt)".into(),
        title: Some(LabelSource::text(&format!(
            "Native {marker}: α < β & tensor slots"
        ))),
        pages: labels
            .into_iter()
            .map(|(mathml, typst)| LabelSource {
                mathml: mathml.into(),
                typst: Some(typst.into()),
            })
            .collect(),
        layout: serde_json::from_value(
            json!({"layout-algo": "impred", "impred-steps": 0, "impred-spacing": 4.0}),
        )
        .unwrap(),
        layout_passes: Vec::new(),
        layout_edges: None,
        label_feedback: false,
    }
}

#[pyfunction]
fn fresh_view(py: Python<'_>, marker: &str) -> PyResult<PyDiagramView> {
    PyDiagramView::from_scene(py, fixture_scene(marker), None)
}

/// A second native oracle uses the browser's actual measurements, never guessed
/// dimensions. This checks that the retained view reaches Scene::render intact.
#[pyfunction]
fn native_oracle(marker: &str, receipt: &str) -> PyResult<String> {
    let receipt: Value =
        serde_json::from_str(receipt).map_err(|error| PyValueError::new_err(error.to_string()))?;
    let metrics = receipt["metrics"]
        .as_array()
        .ok_or_else(|| PyValueError::new_err("metrics"))?
        .iter()
        .map(|entry| {
            serde_json::from_value::<LabelMetrics>(json!({
                "width": entry["width"], "height": entry["height"], "baseline": entry["baseline"]
            }))
            .map_err(|error| PyValueError::new_err(error.to_string()))
        })
        .collect::<PyResult<Vec<_>>>()?;
    let scene = fixture_scene(marker);
    let typeset = scene
        .typeset_mathml(&metrics)
        .map_err(|error| PyValueError::new_err(error.to_string()))?;
    scene.render(&typeset).map_err(PyValueError::new_err)
}

fn run_python(entry: &str) -> PyResult<()> {
    Python::initialize();
    Python::attach(|py| {
        let module = PyModule::new(py, "mathml_native_fixture")?;
        linnet_render_py::register(&module)?;
        module.add_function(wrap_pyfunction!(fresh_view, &module)?)?;
        module.add_function(wrap_pyfunction!(native_oracle, &module)?)?;
        let locals = PyDict::new(py);
        locals.set_item("native", module)?;
        locals.set_item("javascript", include_str!("../src/view.js"))?;
        locals.set_item("adapter_source", include_str!("../src/view.py"))?;
        locals.set_item(
            "test_directory",
            concat!(env!("CARGO_MANIFEST_DIR"), "/tests"),
        )?;
        let source = CString::new(include_str!("mathml_view_support.py")).unwrap();
        py.run(&source, Some(&locals), Some(&locals))?;
        locals.get_item(entry)?.unwrap().call0()?;
        Ok(())
    })
}

#[test]
fn mathml_view_receipt_protocol_and_native_capture() {
    run_python("verify_native").unwrap();
}

#[test]
fn mathml_view_fixture_seed_scale_reaches_native_scene_geometry() {
    let configured_scale = fixture_scene("scale").layout["impred-spacing"]
        .as_f64()
        .unwrap();
    assert_eq!(configured_scale, 4.0);
    // A label-free structural probe checks actual native nodes and external
    // endpoints, without invented font metrics or a browser/Python adapter.
    let positions = |scale| {
        let mut scene = fixture_scene("scale");
        scene.layout.insert("impred-spacing".into(), json!(scale));
        scene.title = None;
        scene.pages.clear();
        for node in &mut scene.nodes {
            node.label = None;
        }
        for edge in &mut scene.edges {
            edge.label = None;
        }
        let typeset = scene.typeset_mathml(&[]).unwrap();
        let svg = scene.render(&typeset).unwrap();
        let document = Document::parse(&svg).unwrap();
        let mut centers: Vec<[f64; 2]> = (0..3)
            .map(|id| {
                let id = id.to_string();
                let shape = document
                    .descendants()
                    .find(|node| node.attribute("data-linnet-node") == Some(id.as_str()))
                    .unwrap();
                let value = |name| shape.attribute(name).unwrap().parse::<f64>().unwrap();
                [
                    value("x") + value("width") / 2.0,
                    value("y") + value("height") / 2.0,
                ]
            })
            .collect();
        for (hedge, free_start) in [(6, true), (7, false)] {
            let prefix = format!("#linnet-halfedge-{hedge}?");
            let carrier = document
                .descendants()
                .find(|node| {
                    node.attribute("data-linnet-carrier") == Some("true")
                        && node
                            .attribute("href")
                            .is_some_and(|href| href.starts_with(&prefix))
                })
                .unwrap();
            let data = carrier
                .first_element_child()
                .unwrap()
                .attribute("d")
                .unwrap();
            let coordinates: Vec<f64> = data
                .split(|character: char| {
                    character.is_ascii_alphabetic() || character.is_whitespace() || character == ','
                })
                .filter(|part| !part.is_empty())
                .map(|part| part.parse().unwrap())
                .collect();
            let actual = if free_start {
                &coordinates[..2]
            } else {
                &coordinates[coordinates.len() - 2..]
            };
            centers.push([actual[0], actual[1]]);
        }
        centers
    };
    let baseline = positions(1.0);
    let scaled = positions(configured_scale);
    assert!(baseline[0] != baseline[1]);
    for (before, after) in baseline.iter().zip(&scaled) {
        for axis in 0..2 {
            let expected = (before[axis] - baseline[0][axis]) * configured_scale;
            let actual = after[axis] - scaled[0][axis];
            // Each serialized coordinate rounds to 1e-4pt. Scaling a
            // difference propagates two such rounding intervals per drawing.
            assert!(
                (actual - expected).abs() <= (1.0 + configured_scale) * 1e-4,
                "EC scale did not reach native geometry: {actual} versus {expected}"
            );
        }
    }

    // Actual 9px SVG-context Chromium 151 receipt from attempt-4. This locks
    // the motivating allocated-box failure before another browser matrix.
    let metrics: Vec<LabelMetrics> = [
        (169.1256606471119, 9.000000108719297, 6.999999925613113),
        (13.623100090234727, 9.468749794005543, 7.46875032615789),
        (24.18759448279095, 10.343749874114499, 7.406250320435822),
        (44.42009539099762, 12.999999759673134, 8.812500091553092),
        (16.92194110133243, 17.765624301907675, 11.515625160217912),
        (8.250032226688386, 9.000000108719297, 6.999999925613113),
    ]
    .into_iter()
    .map(|(width, height, baseline)| LabelMetrics {
        width,
        height,
        baseline,
    })
    .collect();
    let scene = fixture_scene("chromium-350-a");
    let svg = scene
        .render(&scene.typeset_mathml(&metrics).unwrap())
        .unwrap();
    let document = Document::parse(&svg).unwrap();
    let boxes: Vec<[f64; 4]> = document
        .descendants()
        .filter(|node| node.tag_name().name() == "foreignObject")
        .map(|node| {
            let transform = node.parent().unwrap().attribute("transform").unwrap();
            let mut coordinates = transform
                .strip_prefix("translate(")
                .unwrap()
                .strip_suffix(')')
                .unwrap()
                .split_whitespace()
                .map(|value| value.parse::<f64>().unwrap());
            let (x, y) = (coordinates.next().unwrap(), coordinates.next().unwrap());
            [
                x,
                y,
                x + node.attribute("width").unwrap().parse::<f64>().unwrap(),
                y + node.attribute("height").unwrap().parse::<f64>().unwrap(),
            ]
        })
        .collect();
    assert_eq!(boxes.len(), 9);
    for (index, a) in boxes.iter().enumerate() {
        for (other, b) in boxes.iter().enumerate().skip(index + 1) {
            assert!(
                a[2].min(b[2]) <= a[0].max(b[0]) || a[3].min(b[3]) <= a[1].max(b[1]),
                "archived actual browser label rectangles overlap: {index}/{other}: {a:?} / {b:?}"
            );
        }
    }
}

#[test]
fn mathml_view_scene_and_collection_ownership() {
    Python::initialize();
    Python::attach(|py| -> PyResult<()> {
        let scenes = PyDiagramView::from_scenes(py, vec![fixture_scene("channel-a"), fixture_scene("channel-b")], None)?;
        let first = PyDiagramView::from_scene(py, fixture_scene("child-a"), None)?;
        let second = PyDiagramView::from_scene(py, fixture_scene("child-b"), None)?;
        let collection = PyDiagramView::from_views(vec![first.clone(), second.clone()], Some(
            "<section data-collection='native'><!--linnet-svg-1--><!--linnet-svg-0--></section>".into(),
        ))?;
        let locals = PyDict::new(py);
        locals.set_item("scenes", Py::new(py, scenes)?)?;
        locals.set_item("first", Py::new(py, first)?)?;
        locals.set_item("second", Py::new(py, second)?)?;
        locals.set_item("collection", Py::new(py, collection)?)?;
        let source = CString::new(include_str!("mathml_view_collections.py")).unwrap();
        py.run(&source, Some(&locals), Some(&locals))
    }).unwrap();
}

#[test]
#[ignore = "requires root GO and the explicit browser client; finite local producer"]
fn mathml_view_browser_producer() {
    run_python("serve_browser").unwrap();
}

#[test]
#[ignore = "explicit optional Python Typst vector-export verification"]
fn mathml_view_outlined_exports() {
    run_python("verify_outlined").unwrap();
}
