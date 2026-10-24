//! Shared graph presentation bindings, independent of Python graph construction.
/// Bundled offline math font shared by graph and tensor presentation.
///
/// ```no_run
/// use linnet_render_py::data_assets::{STIX_TWO_MATH, STIX_TWO_MATH_LICENSE};
/// use std::error::Error;
/// let font = STIX_TWO_MATH.bytes()?;
/// assert!(!font.is_empty());
/// assert!(STIX_TWO_MATH_LICENSE.contains("SIL OPEN FONT LICENSE"));
/// # Ok::<(), Box<dyn Error>>(())
/// ```
pub mod data_assets {
    /// License text accompanying the embedded font; see the [`crate::data_assets`] example.
    pub const STIX_TWO_MATH_LICENSE: &str = include_str!("../fonts/STIX-Two-OFL.txt");
    include!(concat!(env!("OUT_DIR"), "/data_assets.rs"));
}
#[cfg(feature = "python_stubgen")]
use pyo3::types::{PyAnyMethods, PyStringMethods, PyTypeMethods};
use pyo3::types::{PyModule, PyModuleMethods};
#[cfg(feature = "python_stubgen")]
use pyo3::Python;
use pyo3::{Bound, PyResult};
pub mod drawing;
mod scene;
mod snapshot;
mod svg;
pub mod typst;
pub mod typst_export;
mod view;
pub use scene::register_scene_selectors;
pub use snapshot::PyDiagramRender;
pub use typst::{PyAuto as Auto, PyRenderSettings, PyRenderSettings as RenderSettings};
pub use view::PyDiagramView;

pub fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    drawing::register(module)?;
    typst::register_typst_api(module)?;
    module.add_class::<PyDiagramRender>()?;
    module.add_class::<PyDiagramView>()?;
    Ok(())
}

pub struct GraphModule;
impl symbolica::api::python::SymbolicaCommunityModule for GraphModule {
    fn get_name() -> String {
        "graph".into()
    }
    fn register_module(module: &Bound<'_, PyModule>) -> PyResult<()> {
        register(module)
    }
}

#[cfg(feature = "python_stubgen")]
fn stub_info() -> pyo3_stub_gen::Result<pyo3_stub_gen::StubInfo> {
    pyo3_stub_gen::StubInfo::from_project_root(
        "symbolica".into(),
        std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")),
    )
}
#[cfg(feature = "python_stubgen")]
const STUB_TYPE_ALIASES: &str = r#"class _TypstValueRef(typing.Protocol): ...
class _TypstContentRef(typing.Protocol): ...
class _TypstFunctionRef(typing.Protocol):
    def call(self, *args: _NativeValue, **kwargs: _NativeValue) -> TypstCall: ...
    def bind(self, *args: _NativeValue, **kwargs: _NativeValue) -> TypstBind: ...

_FunctionExpression: typing.TypeAlias = _TypstFunctionRef | TypstBind
_TypstExpression: typing.TypeAlias = _TypstValueRef | _TypstContentRef | _FunctionExpression | TypstCall
_NativeArray: typing.TypeAlias = builtins.list["_NativeValue"] | builtins.tuple["_NativeValue", ...]
_NativeDict: typing.TypeAlias = builtins.dict[builtins.str, "_NativeValue"]
_NativeValue: typing.TypeAlias = (
    None
    | builtins.bool
    | builtins.int
    | builtins.float
    | builtins.str
    | Auto
    | Inherit
    | LayoutAlgorithm
    | LayoutDirection
    | LabelLayout
    | RankAlignment
    | LayoutNodes
    | Placement
    | Compass
    | Routing
    | RoutePoints
    | DanglingTangent
    | Anchor
    | Pattern
    | EdgeLengthResolution
    | DebugLevel
    | StrokeCap
    | StrokeJoin
    | TextStyle
    | DashPattern
    | MarkSymbol
    | MarkPosition
    | MarkOrientation
    | MarkDirection
    | Length
    | Ratio
    | RelativeLength
    | Angle
    | Fraction
    | Color
    | Stroke
    | Dash
    | Insets
    | Mark
    | TextLabel
    | MathSymbol
    | _TypstExpression
    | _NativeArray
    | _NativeDict
)
_ValueExpression: typing.TypeAlias = _TypstValueRef | TypstCall
_Number: typing.TypeAlias = builtins.int | builtins.float | _ValueExpression | Inherit
_Integer: typing.TypeAlias = builtins.int | _ValueExpression | Inherit
_Boolean: typing.TypeAlias = builtins.bool | _ValueExpression | Inherit
_LengthValue: typing.TypeAlias = builtins.int | builtins.float | Length | Ratio | RelativeLength | _ValueExpression | Inherit
_LengthArray: typing.TypeAlias = builtins.list[_LengthValue] | builtins.tuple[_LengthValue, ...]
_Paint: typing.TypeAlias = Color | _ValueExpression | Inherit
_StrokeValue: typing.TypeAlias = Stroke | Color | Length | builtins.int | builtins.float | _ValueExpression | Inherit
_DashValue: typing.TypeAlias = Dash | _ValueExpression | Inherit
_StaticContent: typing.TypeAlias = builtins.str | TextLabel | MathSymbol | _TypstContentRef | TypstCall | Inherit
_Content: typing.TypeAlias = _StaticContent | _FunctionExpression
_Dictionary: typing.TypeAlias = _NativeDict | _ValueExpression | Inherit
_StyleValue: typing.TypeAlias = _NativeDict | _ValueExpression | _FunctionExpression
_Style: typing.TypeAlias = _StyleValue | Inherit
_StyleLayers: typing.TypeAlias = _Style | builtins.list[_NativeDict] | builtins.tuple[_NativeDict, ...]
_PlacementValue: typing.TypeAlias = Placement | _NativeDict | _TypstExpression | None | Inherit
_CompassValue: typing.TypeAlias = Compass | None | Inherit
_RoutingValue: typing.TypeAlias = Routing | None | Inherit
_DanglingTangentValue: typing.TypeAlias = DanglingTangent | Auto | Inherit
_AnchorValue: typing.TypeAlias = Anchor | Auto | None | Inherit
_Radius: typing.TypeAlias = builtins.int | builtins.float | builtins.list[builtins.int | builtins.float] | builtins.tuple[builtins.int | builtins.float, ...] | _ValueExpression | Inherit
_Padding: typing.TypeAlias = builtins.int | builtins.float | _NativeArray | _NativeDict | Insets | _ValueExpression | Inherit
_OptionalString: typing.TypeAlias = builtins.str | None
_OptionalRenderSettings: typing.TypeAlias = RenderSettings | None
_OptionalPaint: typing.TypeAlias = _Paint | None
_OptionalLengthValue: typing.TypeAlias = _LengthValue | None
_OptionalStrokeCap: typing.TypeAlias = StrokeCap | None | Inherit
_OptionalStrokeJoin: typing.TypeAlias = StrokeJoin | None | Inherit
_OptionalDashValue: typing.TypeAlias = _DashValue | None
_OptionalStrokeValue: typing.TypeAlias = _StrokeValue | None
_MarkAnchor: typing.TypeAlias = Anchor | Inherit
_MarkShorten: typing.TypeAlias = _LengthValue | None | Auto
_TextStyleValue: typing.TypeAlias = TextStyle | Inherit
_AutoOptionalStaticContent: typing.TypeAlias = _StaticContent | None | Auto
_OptionalStaticContent: typing.TypeAlias = _StaticContent | None
_AutoOptionalContent: typing.TypeAlias = _Content | None | Auto
_OptionalContent: typing.TypeAlias = _Content | None
_OptionalStyle: typing.TypeAlias = _Style | None
_OptionalStyleLayers: typing.TypeAlias = _StyleLayers | None
_OptionalBoolean: typing.TypeAlias = _Boolean | None
_OptionalInteger: typing.TypeAlias = _Integer | None
_DrawingCoordinate: typing.TypeAlias = builtins.int | builtins.float | _ValueExpression
_DrawingPointDict = typing.TypedDict(
    "_DrawingPointDict",
    {"x": _DrawingCoordinate, "y": _DrawingCoordinate},
)
_DrawingPoint: typing.TypeAlias = builtins.list[_DrawingCoordinate] | builtins.tuple[_DrawingCoordinate, _DrawingCoordinate] | _DrawingPointDict | _ValueExpression | None | Inherit
_DrawingPoints: typing.TypeAlias = typing.Sequence[builtins.list[_DrawingCoordinate] | builtins.tuple[_DrawingCoordinate, _DrawingCoordinate] | _DrawingPointDict | _ValueExpression] | _ValueExpression | None | Inherit
_DrawingString: typing.TypeAlias = builtins.str | _ValueExpression | None | Inherit
_DrawingAngle: typing.TypeAlias = builtins.int | builtins.float | Angle | _ValueExpression | None | Inherit
_DrawingDecoration: typing.TypeAlias = Pattern | _StyleLayers | None | Inherit
_NodeDrawingSelector: typing.TypeAlias = typing.Callable[[typing.Any], NodeDrawing | None] | None | Inherit
_EdgeDrawingSelector: typing.TypeAlias = typing.Callable[[typing.Any], EdgeDrawing | None] | None | Inherit
_HalfEdgeDrawingSelector: typing.TypeAlias = typing.Callable[[typing.Any], HalfEdgeDrawing | None] | None | Inherit
_HedgeSelection: typing.TypeAlias = builtins.list[builtins.bool] | builtins.tuple[builtins.bool, ...]
_SubgraphExpression: typing.TypeAlias = _ValueExpression | _FunctionExpression
_SubgraphSelection: typing.TypeAlias = _HedgeSelection | _SubgraphExpression
_OptionalHedgeSelection: typing.TypeAlias = _SubgraphSelection | None | Inherit
_NodeIndexArray: typing.TypeAlias = builtins.list[builtins.int] | builtins.tuple[builtins.int, ...]
_NodeIndices: typing.TypeAlias = _NodeIndexArray | _ValueExpression | Inherit
_NodeGroup: typing.TypeAlias = _NodeIndexArray | _SubgraphExpression
_NodeGroups: typing.TypeAlias = builtins.list[_NodeGroup] | builtins.tuple[_NodeGroup, ...] | _ValueExpression | Inherit
_LabelLayoutValue: typing.TypeAlias = LabelLayout | Inherit
_LayoutAlgorithmValue: typing.TypeAlias = LayoutAlgorithm | Inherit
_LayoutNodesValue: typing.TypeAlias = LayoutNodes | Inherit
_LayoutDirectionValue: typing.TypeAlias = LayoutDirection | Inherit
_RankAlignmentValue: typing.TypeAlias = RankAlignment | Inherit
_DrawSubgraphRecord = typing.TypedDict(
    "_DrawSubgraphRecord",
    {"subgraph": _SubgraphSelection},
)
_StyledDrawSubgraphRecord = typing.TypedDict(
    "_StyledDrawSubgraphRecord",
    {"subgraph": _SubgraphSelection, "edge-style": _StyleValue | None},
)
_DrawSubgraphEntry: typing.TypeAlias = _SubgraphSelection | _DrawSubgraphRecord | _StyledDrawSubgraphRecord
_DrawSubgraphs: typing.TypeAlias = _SubgraphSelection | builtins.list[_DrawSubgraphEntry] | builtins.tuple[_DrawSubgraphEntry, ...] | None | Inherit
_AutoLengthValue: typing.TypeAlias = _LengthValue | Auto
_AutoNumber: typing.TypeAlias = _Number | Auto
_DebugValue: typing.TypeAlias = DebugLevel | Inherit
_AutoRadius: typing.TypeAlias = _Radius | Auto
_AutoFunction: typing.TypeAlias = _FunctionExpression | Auto
_OptionalNumber: typing.TypeAlias = _Number | None
_EdgeLengthResolver: typing.TypeAlias = EdgeLengthResolution | _ValueExpression | _FunctionExpression | Inherit
_OptionalPadding: typing.TypeAlias = _Padding | None
_RenderStyle: typing.TypeAlias = GraphStyleOptions | None | Inherit
_RenderLayouts: typing.TypeAlias = LayoutSettings | None | Inherit
_RenderDrawing: typing.TypeAlias = DrawOptions | None | Inherit
_RenderSelectors: typing.TypeAlias = DrawingSelectors | None | Inherit
"#;

/// Render the installed and documented Python surface from one stub inventory.
#[cfg(feature = "python_stubgen")]
pub fn canonical_stub() -> pyo3_stub_gen::Result<String> {
    let info = stub_info()?;
    let module = info
        .modules
        .get("symbolica.community.graph")
        .expect("render StubInfo must contain the graph display module");
    let mut exports = module
        .class
        .values()
        .map(|class| class.name)
        .chain(module.enum_.values().map(|enumeration| enumeration.name))
        .chain(module.function.keys().copied())
        .chain(module.variables.values().map(|variable| variable.name))
        .collect::<Vec<_>>();
    exports.sort_unstable();
    exports.dedup();

    // Rust collection constructors are valid PyO3 defaults but are not Python
    // syntax. Keep the generated stub's conventional unspecified default.
    let generated = module.to_string().replace(" = Vec :: new()", " = ...");
    Python::initialize();
    let import_end_line = Python::attach(|py| -> PyResult<usize> {
        let ast = py.import("ast")?;
        let mut import_end = 0;
        for node in ast
            .call_method1("parse", (&generated,))?
            .getattr("body")?
            .try_iter()?
        {
            let node = node?;
            if matches!(node.get_type().name()?.to_str()?, "Import" | "ImportFrom") {
                import_end = node.getattr("end_lineno")?.extract()?;
            }
        }
        Ok(import_end)
    })?;
    let import_end = generated
        .split_inclusive('\n')
        .take(import_end_line)
        .map(str::len)
        .sum::<usize>();
    let declarations = generated[import_end..].trim_matches('\n');
    let exports = exports
        .into_iter()
        .map(|name| format!("    \"{name}\","))
        .collect::<Vec<_>>()
        .join("\n");
    let mut stub = format!(
        "{}\n__all__ = [\n{exports}\n]\n\n{STUB_TYPE_ALIASES}\n\n{declarations}\n",
        &generated[..import_end]
    );
    // pyo3-stub-gen retains Rust identifiers for classes declared by macros,
    // so normalize those inventory names to the names registered with Python;
    // module variables such as AUTO and INHERIT are already included.
    for (rust, python) in [
        ("PyEdgeLengthResolution", "EdgeLengthResolution"),
        ("PyDirectionBasis", "DirectionBasis"),
        ("PyHalfEdgeDrawing", "HalfEdgeDrawing"),
        ("PyHalfEdgeValue", "HalfEdgeValue"),
        ("PyLayoutAlgorithm", "LayoutAlgorithm"),
        ("PyLayoutDirection", "LayoutDirection"),
        ("PyMarkOrientation", "MarkOrientation"),
        ("PyNodeDrawing", "NodeDrawing"),
        ("PyRankAlignment", "RankAlignment"),
        ("PyEdgeDrawing", "EdgeDrawing"),
        ("PyEdgeValue", "EdgeValue"),
        ("PyHalfEdge", "HalfEdge"),
        ("PyOrientedCut", "OrientedCut"),
        ("PyLabelLayout", "LabelLayout"),
        ("PyMarkDirection", "MarkDirection"),
        ("PyMarkPosition", "MarkPosition"),
        ("PyLayoutNodes", "LayoutNodes"),
        ("PyDashPattern", "DashPattern"),
        ("PyMarkSymbol", "MarkSymbol"),
        ("PyDrawingSelectors", "DrawingSelectors"),
        ("PyRoutePoints", "RoutePoints"),
        ("PyDanglingTangent", "DanglingTangent"),
        ("PyStrokeJoin", "StrokeJoin"),
        ("PyTextStyle", "TextStyle"),
        ("PyNodeValue", "NodeValue"),
        ("PySubgraph", "Subgraph"),
        ("PyTraversalTree", "TraversalTree"),
        ("PyCycle", "Cycle"),
        ("PyCutPartition", "CutPartition"),
        ("PyPlacement", "Placement"),
        ("PyStrokeCap", "StrokeCap"),
        ("PyCompass", "Compass"),
        ("PyPattern", "Pattern"),
        ("PyRouting", "Routing"),
        ("PyAnchor", "Anchor"),
        ("PyEdge", "Edge"),
        ("PyNode", "Node"),
    ] {
        stub = stub.replace(rust, python);
    }
    // Stub generation runs without signature inference, so non-builtin Rust
    // values collapse to `...`. These closed enums have one API-wide default.
    for (type_name, default) in [
        ("NodeStore", "NodeStore.Vec"),
        ("Orientation", "Orientation.Default"),
        ("DirectionBasis", "DirectionBasis.Underlying"),
    ] {
        stub = stub.replace(
            &format!(": {type_name} = ..."),
            &format!(": {type_name} = {default}"),
        );
    }
    stub = stub.replace(
        "name: builtins.str = ..., *, payload:",
        "name: builtins.str = '', *, payload:",
    );
    // Generated docstrings can retain indentation on blank lines; keep the
    // checked-in canonical stub free of trailing whitespace.
    let mut canonical = stub
        .lines()
        .map(str::trim_end)
        .collect::<Vec<_>>()
        .join("\n");
    canonical.push('\n');
    Ok(canonical)
}

#[cfg(test)]
mod tests {
    use super::*;
    use linnest::svg::{Dash, Details, EdgeDrawing, LabelSource, NodeDrawing, Scene, Stroke};
    use pyo3::exceptions::PyRuntimeError;
    use pyo3::types::{
        PyAnyMethods, PyDict, PyDictMethods, PyStringMethods, PyType, PyTypeMethods,
    };
    use pyo3::{Py, PyRef, Python};
    use serde_json::json;
    use std::fs;
    use std::path::PathBuf;
    use symbolica::api::python::SymbolicaCommunityModule;

    #[test]
    fn community_module_registers_display_types_without_graph_construction() {
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let module = PyModule::new(py, "symbolica.community.graph")?;
            GraphModule::register_module(&module)?;
            for name in [
                "DiagramRender",
                "DiagramView",
                "RenderSettings",
                "LayoutSettings",
                "Stroke",
                "DrawOptions",
            ] {
                let class = module.getattr(name)?.cast_into::<PyType>()?;
                assert_eq!(class.module()?.to_str()?, "symbolica.community.graph");
            }
            for name in [
                "Graph",
                "Subgraph",
                "NodeStore",
                "Node",
                "Edge",
                "node",
                "source",
                "sink",
                "edge",
                "build",
            ] {
                assert!(
                    !module.hasattr(name)?,
                    "unexpected generic graph export {name}"
                );
            }
            Ok(())
        })
        .unwrap();
    }

    #[test]
    fn full_document_rendering_is_not_exposed() {
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let module = PyModule::new(py, "symbolica.community.graph")?;
            register(&module)?;
            let locals = PyDict::new(py);
            locals.set_item("graph", module)?;
            py.run(
                c"assert not hasattr(graph.DiagramRender, 'from_sources')",
                Some(&locals),
                Some(&locals),
            )?;
            py.run(
                c"assert not hasattr(graph.DiagramRender, 'typst_source')",
                Some(&locals),
                Some(&locals),
            )?;
            for option in ["template", "source_root", "template_options"] {
                let kwargs = PyDict::new(py);
                kwargs.set_item(option, py.None())?;
                assert!(locals
                    .get_item("graph")?
                    .unwrap()
                    .getattr("RenderSettings")?
                    .call((), Some(&kwargs))
                    .is_err());
            }
            Ok(())
        })
        .unwrap();
    }

    #[test]
    fn native_layout_passes_preserve_labels_paint_and_exports_without_graph_bindings() {
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let module = PyModule::new(py, "symbolica.community.graph")?;
            register(&module)?;
            let locals = PyDict::new(py);
            locals.set_item("graph", module)?;
            py.run(
                cr#"
config = graph.RenderSettings(
    layouts=graph.LayoutSettings(impred_steps=1).then(impred_steps=1),
)
"#,
                Some(&locals),
                Some(&locals),
            )?;
            let settings = locals.get_item("config")?.unwrap();
            let settings = settings.extract::<PyRef<'_, PyRenderSettings>>()?;
            let stroke = Stroke {
                paint: "#246813".into(),
                width: 0.6,
                dash: Dash::Solid,
                round_cap: false,
            };
            let node = NodeDrawing {
                radius: 0.2,
                label_padding: None,
                label: None,
                rectangular: false,
                fill: "#357913".into(),
                stroke: stroke.clone(),
                details: Details::default(),
            };
            let scene = Scene {
                graph: serde_json::from_value(json!({
                    "nodes": [{"index": 0}, {"index": 1}],
                    "edges": [{"id": 0, "source": {"node": 0, "id": 0}, "sink": {"node": 1, "id": 1}}],
                })).unwrap(),
                nodes: vec![
                    NodeDrawing { label: Some(0), ..node.clone() },
                    node,
                ],
                edges: vec![EdgeDrawing {
                    stroke,
                    pattern: None,
                    flow: None,
                    momentum: false,
                    momentum_offset: EdgeDrawing::DEFAULT_MOMENTUM_OFFSET,
                    label: None,
                    details: [("particle", "phi")].into_iter().collect(),
                }],
                preamble: String::new(),
                title: None,
                pages: vec![LabelSource {
                    mathml: "<mi>x</mi>".into(),
                    typst: Some("[$x$]".into()),
                }],
                layout: Default::default(),
                layout_passes: Vec::new(),
                layout_edges: None,
                label_feedback: false,
            };
            let view = PyDiagramView::from_scene(py, scene, Some(&settings))?;
            let request = view.measurement_request(py)?;
            locals.set_item("view", Py::new(py, view.clone())?)?;
            locals.set_item("request", request)?;
            // EXTERNAL METRICS FIXTURE in points; no browser or compiler ran.
            py.run(
                cr#"
assert request["labels"] == [{"id": 0, "mathml": "<mi>x</mi>"}]
receipt = {"id": request["id"], "metrics": [
    {"id": 0, "width": 12.0, "height": 9.0, "baseline": 7.0},
]}
assert "awaiting browser measurement" in repr(view)
"#,
                Some(&locals),
                Some(&locals),
            )?;
            let receipt = locals.get_item("receipt")?.unwrap().cast_into::<PyDict>()?;
            match view.measured(py, None, None) {
                Ok(_) => panic!("unmeasured labels produced a snapshot"),
                Err(error) => assert!(error.to_string().contains("measurement is pending")),
            }
            let result = view.measured(py, Some(&receipt), None)?;
            let svg = result.to_svg()?;
            for value in [
                "<svg", "#357913", "#246813", "particle", "<mi>x</mi>", "foreignObject",
            ] {
                assert!(svg.contains(value), "missing scene data {value}");
            }
            assert_eq!(view.measured(py, None, None)?.to_svg()?, svg);
            let output = tempfile::tempdir()?;
            let result = Py::new(py, result)?;
            let error = result.bind(py)
                .call_method1("save", (output.path().join("scene.pdf"),))
                .unwrap_err();
            assert!(error.is_instance_of::<PyRuntimeError>(py));
            assert!(error.to_string().contains("require outlined labels"));
            Ok(())
        })
        .unwrap();
    }

    #[test]
    fn captured_vector_svg_preserves_pdf_png_and_typst_exports() {
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            // A genuine vector path fixture, independent of label compilation.
            let svg = r##"<svg xmlns="http://www.w3.org/2000/svg" width="24pt" height="18pt" viewBox="0 0 24 18"><path d="M2 2H22V16H2Z" fill="#357913"/></svg>"##;
            let snapshot = PyDiagramRender::new(svg.into(), format!("<figure>{svg}</figure>"));
            let collection = PyDiagramRender::from_diagrams(
                vec![snapshot.clone(), snapshot.clone()],
                format!("<figure>{svg}{svg}</figure>"),
            )?;
            assert_eq!(collection.to_svg_pages()?.len(), 2);
            let result = Py::new(py, snapshot)?;
            let output = tempfile::tempdir()?;
            for suffix in ["pdf", "png", "typ"] {
                let path = output.path().join(format!("vector.{suffix}"));
                result.bind(py).call_method1("save", (&path,))?;
                let bytes = fs::read(path)?;
                match suffix {
                    "pdf" => assert!(bytes.starts_with(b"%PDF")),
                    "png" => assert!(bytes.starts_with(b"\x89PNG\r\n\x1a\n")),
                    _ => assert!(String::from_utf8(bytes).unwrap().contains("#image(bytes(")),
                }
            }
            let collection = Py::new(py, collection)?;
            let pdf = output.path().join("collection.pdf");
            collection.bind(py).call_method1("save", (&pdf,))?;
            assert!(fs::read(pdf)?.starts_with(b"%PDF"));
            let kwargs = PyDict::new(py);
            kwargs.set_item("format", "png")?;
            let pages = collection.bind(py)
                .call_method("save_pages", (output.path().join("pages"),), Some(&kwargs))?
                .extract::<Vec<PathBuf>>()?;
            assert_eq!(pages.len(), 2);
            for page in pages {
                assert!(fs::read(page)?.starts_with(b"\x89PNG\r\n\x1a\n"));
            }
            Ok(())
        })
        .unwrap();
    }
}
