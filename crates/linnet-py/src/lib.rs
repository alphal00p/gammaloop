//! Python bindings for Linnet's typed graph, DOT, and Typst drawing APIs.

use pyo3::prelude::*;

mod dot;
mod drawing;
mod graph;
mod mutations;
mod render;
pub use snapshot::PyDiagramRender;
mod scene;
mod snapshot;
pub use graph::PyGraph;
pub use topology::PySubgraph;
pub use typst::PyRenderSettings;
pub use typst::{PyAuto as Auto, PyRenderSettings as RenderSettings};
mod streaming;
mod svg;
mod topology;
mod typst;

pub fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    graph::register(module)?;
    topology::register(module)?;
    drawing::register(module)?;
    dot::register(module)?;
    typst::register_typst_api(module)?;
    module.add_class::<PyDiagramRender>()?;
    streaming::register(module)?;
    Ok(())
}

// Classes take the Community module path whenever `community` is enabled, so a
// standalone module in that build would publish mislabeled classes. Feature
// unification would also export `PyInit_linnet_py` from Symbolica's extension.
#[cfg(all(feature = "standalone", feature = "community"))]
compile_error!("the `standalone` linnet_py module cannot be built with `community`");

/// Graphs, topology, layout, and rendering backed by Linnet.
#[cfg(feature = "standalone")]
#[pymodule]
fn linnet_py(module: &Bound<'_, PyModule>) -> PyResult<()> {
    register(module)
}

/// A Python module that publishes Linnet's classes.
#[derive(Clone, Copy, Debug)]
pub enum PythonModule {
    /// The standalone `linnet-py` distribution.
    Standalone,
    /// The Symbolica Community submodule registered by [`GraphModule`].
    #[cfg(feature = "community")]
    Community,
}

impl PythonModule {
    pub const fn name(self) -> &'static str {
        match self {
            PythonModule::Standalone => "linnet_py",
            #[cfg(feature = "community")]
            PythonModule::Community => "symbolica.community.graph",
        }
    }
}

#[cfg(feature = "python_stubgen")]
pub fn stub_info() -> pyo3_stub_gen::Result<pyo3_stub_gen::StubInfo> {
    // Dependency classes without an explicit module belong to Symbolica;
    // Linnet declarations explicitly name their own module.
    pyo3_stub_gen::StubInfo::from_project_root(
        "symbolica".to_owned(),
        std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")),
    )
}

// Rust's `PyAny` boundary is intentionally wider than the Python API: values
// are recursively checked before they can cross into Typst.  Keep the aliases
// here so manual signatures and generated signatures share one closed model.
#[cfg(feature = "python_stubgen")]
const STUB_TYPE_ALIASES: &str = r#"_LabelKey: typing.TypeAlias = typing.Callable[[typing.Any], _Key] | None

class _TypstValueRef(typing.Protocol): ...
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
_EndpointTarget: typing.TypeAlias = Node | NodeView | builtins.int | builtins.str
_EdgeEndpoint: typing.TypeAlias = HalfEdge | _EndpointTarget | None
_NodeInput: typing.TypeAlias = Node | builtins.str | None
_EdgeInput: typing.TypeAlias = Edge | _EdgeEndpoint
_OptionalOrientation: typing.TypeAlias = Orientation | None
_OptionalNodeDrawing: typing.TypeAlias = NodeDrawing | None
_OptionalEdgeDrawing: typing.TypeAlias = EdgeDrawing | None
_OptionalHalfEdgeDrawing: typing.TypeAlias = HalfEdgeDrawing | None
_HalfEdgeTarget: typing.TypeAlias = HalfEdgeView | builtins.int
_OptionalGlobalData: typing.TypeAlias = GlobalData | None
_OptionalDotCodec: typing.TypeAlias = DotCodec | None
_OptionalRenderSettings: typing.TypeAlias = RenderSettings | None
_OptionalPaint: typing.TypeAlias = _Paint | None
_OptionalLengthValue: typing.TypeAlias = _LengthValue | None
_OptionalStrokeCap: typing.TypeAlias = StrokeCap | None | Inherit
_OptionalStrokeJoin: typing.TypeAlias = StrokeJoin | None | Inherit
_OptionalDashValue: typing.TypeAlias = _DashValue | None
_OptionalStrokeValue: typing.TypeAlias = _StrokeValue | None
_MarkSize: typing.TypeAlias = builtins.str | Length | Ratio | RelativeLength
_MarkPaint: typing.TypeAlias = Color | builtins.str | builtins.bool | None | Auto
_AutoMarkSize: typing.TypeAlias = _MarkSize | Auto
_MarkRatio: typing.TypeAlias = builtins.str | Ratio
_MarkStroke: typing.TypeAlias = _MarkPaint | Stroke
_MarkPhase: typing.TypeAlias = Angle | Auto
_MarkPart: typing.TypeAlias = Mark | _NativeDict | _MarkSize
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
_NodeDrawingSelector: typing.TypeAlias = typing.Callable[[NodeView], NodeDrawing | None] | None | Inherit
_EdgeDrawingSelector: typing.TypeAlias = typing.Callable[[EdgeView], EdgeDrawing | None] | None | Inherit
_HalfEdgeDrawingSelector: typing.TypeAlias = typing.Callable[[HalfEdgeView], HalfEdgeDrawing | None] | None | Inherit
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
_TemplatePath: typing.TypeAlias = builtins.str | os.PathLike[builtins.str] | None | Inherit
_SourceRootPath: typing.TypeAlias = builtins.str | os.PathLike[builtins.str] | None | Inherit
_RenderStyle: typing.TypeAlias = GraphStyleOptions | None | Inherit
_RenderLayouts: typing.TypeAlias = LayoutSettings | None | Inherit
_RenderDrawing: typing.TypeAlias = DrawOptions | None | Inherit
_RenderSelectors: typing.TypeAlias = DrawingSelectors | None | Inherit
_TemplateOptions: typing.TypeAlias = _NativeDict | None | Inherit
"#;

#[cfg(feature = "python_stubgen")]
impl PythonModule {
    /// The module every class is compiled into, and so filed under in the stub
    /// inventory. `#[pyclass]` fixes the qualified type name at compile time,
    /// so each class selects it with
    /// `#[cfg_attr(feature = "community", pyo3(module = ...))]`.
    #[cfg(feature = "community")]
    pub(crate) const COMPILED: Self = Self::Community;
    #[cfg(not(feature = "community"))]
    pub(crate) const COMPILED: Self = Self::Standalone;

    /// Render this module's installed and documented Python surface. Only the
    /// Community build accepts Symbolica expressions as keys of their own kind.
    pub fn stub(self) -> pyo3_stub_gen::Result<String> {
        match self {
            PythonModule::Standalone => render_stub("", "typing.Hashable"),
            #[cfg(feature = "community")]
            PythonModule::Community => render_stub(
                "from symbolica.core import Expression\n",
                "typing.Hashable | Expression",
            ),
        }
    }
}

/// Render the Python surface from the one stub inventory, with the given
/// import and type of graph keys.
#[cfg(feature = "python_stubgen")]
fn render_stub(key_import: &str, key: &str) -> pyo3_stub_gen::Result<String> {
    let info = stub_info()?;
    let module = info
        .modules
        .get(PythonModule::COMPILED.name())
        .expect("linnet StubInfo must contain the linnet module");
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
        "{}{key_import}\n__all__ = [\n{exports}\n]\n\n_Key: typing.TypeAlias = {key}\n{STUB_TYPE_ALIASES}\n\n{declarations}\n",
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
        ("PyHalfEdgeSpec", "HalfEdge"),
        ("PyNodeSpec", "Node"),
        ("PyEdgeSpec", "Edge"),
        ("PyHalfEdge", "HalfEdgeView"),
        ("PyOrientedCut", "OrientedCut"),
        ("PyLabelLayout", "LabelLayout"),
        ("PyMarkDirection", "MarkDirection"),
        ("PyMarkPosition", "MarkPosition"),
        ("PyLayoutNodes", "LayoutNodes"),
        ("PyDashPattern", "DashPattern"),
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
        ("PyEdge", "EdgeView"),
        ("PyNode", "NodeView"),
    ] {
        stub = stub.replace(rust, python);
    }
    // Stub generation runs without signature inference, so non-builtin Rust
    // values collapse to `...`. These closed enums have one API-wide default.
    for (type_name, default) in [
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

/// Generic graphs and graph presentation share a single Community module.
#[cfg(feature = "community")]
pub struct GraphModule;
#[cfg(feature = "community")]
impl symbolica::api::python::SymbolicaCommunityModule for GraphModule {
    fn get_name() -> String {
        "graph".into()
    }
    fn register_module(module: &Bound<'_, PyModule>) -> PyResult<()> {
        register(module)
    }
}

#[cfg(test)]
mod tests {
    use pyo3::{prelude::*, types::PyType};

    #[test]
    fn classes_are_named_after_the_module_of_their_build() {
        let expected = if cfg!(feature = "community") {
            "symbolica.community.graph"
        } else {
            "linnet_py"
        };
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let module = PyModule::new(py, "graph")?;
            super::register(&module)?;
            let mut classes = 0;
            for (name, value) in module.dict() {
                if let Ok(class) = value.cast::<PyType>() {
                    assert_eq!(class.module()?.to_str()?, expected, "{name}");
                    classes += 1;
                }
            }
            assert!(classes > 50, "{classes} classes");
            // CPython's own messages print the compiled type name, not `__module__`.
            let error = module
                .getattr("AUTO")?
                .getattr("missing")
                .unwrap_err()
                .to_string();
            assert!(
                error.contains(&format!("'{expected}.Auto' object has no attribute")),
                "{error}"
            );
            Ok(())
        })
        .unwrap();
    }

    #[cfg(feature = "python_stubgen")]
    #[test]
    fn stub_keys_include_expressions_only_for_the_community_module() {
        use super::PythonModule;

        let stubs = [(PythonModule::Standalone, "typing.Hashable", false)].into_iter();
        #[cfg(feature = "community")]
        let stubs = stubs.chain([(
            PythonModule::Community,
            "typing.Hashable | Expression",
            true,
        )]);
        for (module, key, imports_expression) in stubs {
            let stub = module.stub().unwrap();
            assert_eq!(
                stub.contains("\nfrom symbolica.core import Expression\n"),
                imports_expression
            );
            assert!(stub.contains(&format!("\n_Key: typing.TypeAlias = {key}\n")));
            assert!(stub.contains("typing.Callable[[typing.Any], _Key] | None"));
            assert!(stub.contains("def data(self) -> _Key:"));
            assert_eq!(stub.matches("typing.Hashable").count(), 1);
        }
        // The module path never enters the stub body, so every build renders
        // the checked-in standalone stub.
        assert_eq!(
            PythonModule::Standalone.stub().unwrap(),
            include_str!("../linnet_py.pyi")
        );
    }
}
