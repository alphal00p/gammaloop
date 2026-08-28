//! Python bindings for Linnet's typed graph, DOT, and Typst drawing APIs.

use pyo3::prelude::*;
use pyo3_stub_gen::define_stub_info_gatherer;

mod dot;
mod drawing;
mod graph;
mod mutations;
mod native_graph;
mod render;
mod topology;
mod typst;

#[pymodule]
fn linnet_py(module: &Bound<'_, PyModule>) -> PyResult<()> {
    graph::register(module)?;
    topology::register(module)?;
    drawing::register(module)?;
    dot::register(module)?;
    typst::register_typst_api(module)?;
    Ok(())
}

define_stub_info_gatherer!(stub_info);

// Rust's `PyAny` boundary is intentionally wider than the Python API: values
// are recursively checked before they can cross into Typst.  Keep the aliases
// here so manual signatures and generated signatures share one closed model.
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
    | Anchor
    | Pattern
    | EdgeLengthResolution
    | DebugLevel
    | StrokeCap
    | StrokeJoin
    | TextStyle
    | RenderMode
    | TypstFields
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
_AnchorValue: typing.TypeAlias = Anchor | Auto | None | Inherit
_Radius: typing.TypeAlias = builtins.int | builtins.float | builtins.list[builtins.int | builtins.float] | builtins.tuple[builtins.int | builtins.float, ...] | _ValueExpression | Inherit
_Padding: typing.TypeAlias = builtins.int | builtins.float | _NativeArray | _NativeDict | Insets | _ValueExpression | Inherit
_FieldNames: typing.TypeAlias = builtins.str | builtins.list[builtins.str] | builtins.tuple[builtins.str, ...] | _ValueExpression | Inherit
_OptionalString: typing.TypeAlias = builtins.str | None
_OptionalHalfEdgeSpec: typing.TypeAlias = HalfEdgeSpec | None
_EndpointTarget: typing.TypeAlias = NodeSpec | Node | builtins.int | builtins.str
_HalfEdgeTarget: typing.TypeAlias = HalfEdge | builtins.int
_GraphItem: typing.TypeAlias = NodeSpec | EdgeSpec
_OptionalGlobalData: typing.TypeAlias = GlobalData | None
_OptionalDotCodec: typing.TypeAlias = DotCodec | None
_OptionalRenderConfig: typing.TypeAlias = RenderConfig | None
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
_DrawingString: typing.TypeAlias = builtins.str | _ValueExpression | None | Inherit
_DrawingScalar: typing.TypeAlias = builtins.str | builtins.int | builtins.float | _ValueExpression | None | Inherit
_DrawingAngle: typing.TypeAlias = builtins.int | builtins.float | Angle | _ValueExpression | None | Inherit
_DrawingDecoration: typing.TypeAlias = Pattern | _StyleLayers | None | Inherit
_NodeSelector: typing.TypeAlias = typing.Callable[[Node], _OptionalStyle] | None | Inherit
_EdgeSelector: typing.TypeAlias = typing.Callable[[Edge], _OptionalStyleLayers] | None | Inherit
_HalfEdgeSelector: typing.TypeAlias = typing.Callable[[HalfEdge], _OptionalStyleLayers] | None | Inherit
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
_TypstFieldsValue: typing.TypeAlias = TypstFields | Inherit
_MomentumMark: typing.TypeAlias = Mark | _NativeDict | _ValueExpression | None | Auto | Inherit
_OptionalFieldNames: typing.TypeAlias = _FieldNames | None
_TemplatePath: typing.TypeAlias = builtins.str | os.PathLike[builtins.str] | None | Inherit
_ExecutablePath: typing.TypeAlias = builtins.str | os.PathLike[builtins.str] | Inherit
_RenderModeValue: typing.TypeAlias = RenderMode | Auto | Inherit
_RenderStyle: typing.TypeAlias = GraphStyleOptions | None | Inherit
_RenderLayouts: typing.TypeAlias = LayoutOptions | None | Inherit
_RenderDrawing: typing.TypeAlias = DrawOptions | None | Inherit
_RenderPhysics: typing.TypeAlias = PhysicsOptions | None | Inherit
"#;

/// Render the installed and documented Python surface from one stub inventory.
pub fn canonical_stub() -> pyo3_stub_gen::Result<String> {
    let info = stub_info()?;
    let module = info
        .modules
        .get("linnet_py")
        .expect("linnet StubInfo must contain the linnet_py module");
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

    let generated = module.to_string();
    let import_end = generated
        .split_inclusive('\n')
        .scan(0, |offset, line| {
            *offset += line.len();
            Some((line, *offset))
        })
        .filter(|(line, _)| {
            let line = line.trim_start();
            line.starts_with("import ") || line.starts_with("from ")
        })
        .map(|(_, offset)| offset)
        .last()
        .expect("linnet stub must import its referenced Python types");
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
        ("PyRenderMode", "RenderMode"),
        ("PyRoutePoints", "RoutePoints"),
        ("PyStrokeJoin", "StrokeJoin"),
        ("PyTextStyle", "TextStyle"),
        ("PyTypstFields", "TypstFields"),
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
    // Rust collection constructors are valid PyO3 defaults but are not Python
    // syntax. Keep the generated stub's conventional unspecified default.
    stub = stub.replace(" = Vec :: new()", " = ...");
    Ok(stub)
}
