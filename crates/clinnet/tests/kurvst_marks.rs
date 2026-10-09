use std::{fs, path::PathBuf, process::Command};

use clinnet::TypstRenderer;
use kurvst::marks::MarkGeometrySpec;
use serde_json::json;

#[test]
fn mark_geometry_matches_native_through_the_typst_wasm_boundary() {
    let typst = std::env::var_os("TYPST_TEST_EXECUTABLE")
        .map_or_else(|| PathBuf::from("typst"), PathBuf::from);
    let version = Command::new(&typst).arg("--version").output().unwrap();
    assert!(
        version.status.success()
            && String::from_utf8_lossy(&version.stdout).starts_with("typst 0.15."),
        "mark parity tests require Typst 0.15"
    );
    let base = tempfile::tempdir().unwrap();
    let renderer = TypstRenderer::new(base.path()).typst_executable(typst);
    renderer.stage_default_assets().unwrap();
    let root = base.path().join(".clinnet/templates");
    let shapes = [
        "triangle", "straight", "stealth", "round", "tikz", "barb", "hooks", "bar", "bracket",
        "circle", "square", "diamond", "rays",
    ];
    let mut templates: Vec<_> = shapes
        .iter()
        .map(|shape| {
            json!({
                "mark": {"shape": shape},
                "context": {"units-per-pt": 1.0, "line-thickness": 1.0}
            })
        })
        .collect();
    templates.push(json!({
        "mark": {
            "shape": "combine",
            "parts": [
                {"shape": "bar"},
                {"gap": {"points": 2.0, "ratio": 0.5}},
                {"shape": "circle"}
            ],
            "fit": "bend"
        },
        "context": {
            "units-per-pt": 0.5, "line-thickness": 0.7,
            "shaft-stroke": {"cap": "round", "join": "bevel", "miter-limit": 3.0}
        }
    }));
    let carriers = json!([
        {"path": {"elements": [
            {"kind": "move", "start": [0.0, 0.0]},
            {"kind": "line", "end": [100.0, 0.0]}
        ]}},
        {"path": {"elements": [
            {"kind": "move", "start": [0.0, 0.0]},
            {"kind": "cubic", "control-start": [5.0, 20.0],
             "control-end": [30.0, -20.0], "end": [40.0, 20.0]}
        ]}},
        {"path": {"elements": [
            {"kind": "move", "start": [0.0, 0.0]},
            {"kind": "line", "end": [100.0, 0.0]}
        ]}},
        {"path": {"elements": [
            {"kind": "move", "start": [0.0, 0.0]},
            {"kind": "line", "end": [2.0, 0.0]}
        ]}}
    ]);
    let candidate_placements: Vec<_> = (0..templates.len())
        .map(|template| {
            json!({
                "template": template, "carrier": usize::from(template == 13),
                "station": {"kind": "end"}
            })
        })
        .collect();
    // Exercise the shared painted-center rule across the complete catalogue
    // and composite, both numeric station forms, both directions and shifts.
    // The composite's distinct shaft context needs its own carrier.
    let mut centered_placements = Vec::new();
    for template in 0..templates.len() {
        for direction in ["forward", "backward"] {
            for station in [
                json!({"kind": "ratio", "value": 0.5}),
                json!({"kind": "distance", "value": 50.0}),
            ] {
                for shift in [-3.0, 3.0] {
                    centered_placements.push(json!({
                        "template": template, "carrier": if template == 13 { 2 } else { 0 },
                        "station": station, "direction": direction, "shift": shift
                    }));
                }
            }
        }
    }
    for direction in ["forward", "backward"] {
        centered_placements.push(json!({
            "template": 0, "carrier": 3, "direction": direction,
            "station": {"kind": "ratio", "value": 0.5}
        }));
    }
    for (mode, placements) in [
        ("candidates", json!(candidate_placements)),
        (
            "selected",
            json!([
                {"template": 0, "carrier": 0, "station": {"kind": "start"},
                 "direction": "backward"},
                {"template": 0, "carrier": 0, "station": {"kind": "end"}},
                {"template": 3, "carrier": 0,
                 "station": {"kind": "ratio", "value": 0.5}}
            ]),
        ),
        ("selected", json!(centered_placements)),
    ] {
        let request = json!({
            "mode": mode,
            "templates": templates,
            "carriers": carriers,
            "placements": placements
        });
        let native: MarkGeometrySpec = serde_json::from_value(request.clone()).unwrap();
        let output = native.geometry().unwrap();
        fs::write(
            root.join("request.json"),
            serde_json::to_vec(&request).unwrap(),
        )
        .unwrap();
        let mut expected = Vec::new();
        ciborium::ser::into_writer(&output, &mut expected).unwrap();
        fs::write(root.join("expected.cbor"), expected).unwrap();
        let input = root.join("parity.typ");
        fs::write(
            &input,
            r#"
#let engine = plugin("crates/kurvst/typst/kurvst.wasm")
#let actual = cbor(engine.mark_geometry(cbor.encode(json("request.json"))))
#let expected = cbor(read("expected.cbor", encoding: none))
#let first-difference(a, b, path: "root") = {
  if a == b { return none }
  let coordinate = (path.ends-with(".end") or path.contains(".elements[")
    or path.contains(".tip[") or path.contains(".back[")
    or path.contains(".shaft-contact[") or path.contains(".outline-bounds[")
    or path.contains(".footprint-bounds["))
  if coordinate and type(a) in (int, float) and type(b) in (int, float) {
    if calc.abs(a - b) <= 1e-12 * calc.max(1, calc.abs(a), calc.abs(b)) {
      return none
    }
    return (path: path, actual: repr(a), expected: repr(b))
  }
  if type(a) == array and type(b) == array and a.len() == b.len() {
    for i in range(a.len()) {
      let difference = first-difference(a.at(i), b.at(i), path: path + "[" + str(i) + "]")
      if difference != none { return difference }
    }
    return none
  } else if type(a) == dictionary and type(b) == dictionary and a.keys() == b.keys() {
    for key in a.keys() {
      let difference = first-difference(a.at(key), b.at(key), path: path + "." + key)
      if difference != none { return difference }
    }
    return none
  }
  (path: path, actual: repr(a), expected: repr(b))
}
#assert.eq(
  first-difference(actual, expected),
  none,
  message: "native/Wasm geometry within 1e-12; structure and styles exact",
)
Mark geometry parity passed.
"#,
        )
        .unwrap();
        renderer
            .compile_template(&input, base.path().join(format!("{mode}.pdf")), &[])
            .unwrap();
    }
}

#[test]
fn public_mark_constructors_and_batch_geometry_behave() {
    let typst = std::env::var_os("TYPST_TEST_EXECUTABLE")
        .map_or_else(|| PathBuf::from("typst"), PathBuf::from);
    let base = tempfile::tempdir().unwrap();
    let renderer = TypstRenderer::new(base.path()).typst_executable(typst);
    renderer.check_version().unwrap();
    renderer.stage_default_assets().unwrap();
    let root = base.path().join(".clinnet/templates");
    for (name, source) in [
        (
            "mark-contract.typ",
            include_str!("../../kurvst/typst/tests/mark-contract.typ"),
        ),
        (
            "mark-geometry.typ",
            include_str!("../../kurvst/typst/tests/mark-geometry.typ"),
        ),
    ] {
        let source = source
            .replace("../src/mark.typ", "crates/kurvst/typst/src/mark.typ")
            .replace("../src/lib.typ", "crates/kurvst/typst/src/lib.typ");
        let input = root.join(name);
        fs::write(&input, &source).unwrap();
        let output = base.path().join(name).with_extension("pdf");
        renderer.compile_template(&input, &output, &[]).unwrap();
        if name == "mark-contract.typ" {
            for mode in [
                "positional",
                "unsupported",
                "tikz-arc",
                "symbol",
                "scale",
                "anchor",
                "shorten-to",
                "cetz",
                "n",
                "n-fraction",
                "align",
                "fit",
                "nonfinite",
                "bad-context",
                "thickness",
                "child-fit",
                "missing-parts",
                "named-parts",
                "fill-content",
                "fill-function",
                "stroke-dictionary",
                "stroke-thickness",
                "stroke-cap",
                "stroke-join",
                "stroke-dash",
                "stroke-miter",
            ] {
                fs::write(
                    &input,
                    source.replace("default: \"valid\"", &format!("default: \"{mode}\"")),
                )
                .unwrap();
                let error = renderer.compile_template(&input, &output, &[]).unwrap_err();
                assert!(error.to_string().contains("kurvst.mark"), "{mode}: {error}");
            }
        }
    }
}
