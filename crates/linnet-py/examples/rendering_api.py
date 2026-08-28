# ruff: noqa: B018, PLR1711  # Cell outputs and empty returns are Marimo syntax.

import marimo

__generated_with = "0.24.0"
app = marimo.App(width="medium")


@app.cell
def _():
    from dataclasses import dataclass

    import linnet_py as lp
    import marimo as mo

    @dataclass(eq=False)
    class TaskRecord:
        owner: str
        retry_limit: int
        runtime_context: object

    @dataclass(eq=False)
    class ChannelRecord:
        protocol: str
        capacity: int
        runtime_context: object

    @dataclass(eq=False)
    class EndpointRecord:
        slot: str
        runtime_context: object

    return ChannelRecord, EndpointRecord, TaskRecord, lp, mo


@app.cell
def _(mo):
    mo.md(r"""
    # Linnet Python rendering

    A graph's ordinary notebook representation is a model-neutral Linnest
    drawing. `Graph._repr_svg_()` powers the automatic display; `to_svg()` and
    `render(path)` expose the same Typst 0.15 pipeline explicitly.

    Domain-specific conventions are ordinary Python render settings, not modes
    built into Linnet. The example data below is an application workflow: every
    node, edge, and half-edge carries an arbitrary Python dataclass instance
    that is never serialized to Typst.
    """)
    return


@app.cell
def _(lp, mo):
    node_store = mo.ui.dropdown(
        options={
            "Dense Vec": lp.NodeStore.Vec,
            "Identification Forest": lp.NodeStore.Forest,
        },
        value="Dense Vec",
        label="Node store",
    )
    layout_algorithm = mo.ui.dropdown(
        options={
            "Stable layered": lp.LayoutAlgorithm.StableLayered,
            "Force": lp.LayoutAlgorithm.Force,
        },
        value="Force",
        label="Layout",
    )
    custom_theme = mo.ui.checkbox(
        value=False,
        label="Custom Python line theme (optional)",
    )
    controls = mo.hstack(
        [node_store, layout_algorithm, custom_theme],
        justify="start",
        gap=1.5,
    )
    controls
    return custom_theme, layout_algorithm, node_store


@app.cell
def _(
    ChannelRecord,
    EndpointRecord,
    TaskRecord,
    custom_theme,
    layout_algorithm,
    lp,
    node_store,
):
    runtime_context = object()
    receive_data = TaskRecord("platform", 3, runtime_context)
    parse_data = TaskRecord("formats", 2, runtime_context)
    validate_data = TaskRecord("quality", 1, runtime_context)
    normalize_data = TaskRecord("formats", 2, runtime_context)
    quarantine_data = TaskRecord("operations", 5, runtime_context)
    index_data = TaskRecord("search", 3, runtime_context)
    audit_data = TaskRecord("compliance", 4, runtime_context)
    publish_data = TaskRecord("delivery", 2, runtime_context)

    receive = lp.node("receive", data=receive_data)
    parse = lp.node("parse", data=parse_data)
    validate = lp.node("validate", data=validate_data)
    normalize = lp.node("normalize", data=normalize_data)
    quarantine = lp.node("quarantine", data=quarantine_data)
    index = lp.node("index", data=index_data)
    audit = lp.node("audit", data=audit_data)
    publish = lp.node("publish", data=publish_data)

    request_channel = ChannelRecord("https", 256, runtime_context)
    request_port = EndpointRecord("request body", runtime_context)
    graph = lp.build(
        receive,
        parse,
        validate,
        normalize,
        quarantine,
        index,
        audit,
        publish,
        lp.edge(
            lp.sink(receive, data=request_port),
            "requests",
            data=request_channel,
        ),
        lp.edge(
            lp.source(
                receive,
                data=EndpointRecord("accepted documents", runtime_context),
            ),
            "decode queue",
            lp.sink(parse, data=EndpointRecord("raw documents", runtime_context)),
            data=ChannelRecord("memory queue", 128, runtime_context),
        ),
        lp.edge(
            lp.source(parse, data=EndpointRecord("syntax tree", runtime_context)),
            "parsed documents",
            lp.sink(validate, data=EndpointRecord("candidate", runtime_context)),
            data=ChannelRecord("bounded queue", 96, runtime_context),
        ),
        lp.edge(
            lp.source(validate, data=EndpointRecord("valid", runtime_context)),
            "accepted",
            lp.sink(normalize, data=EndpointRecord("source record", runtime_context)),
            data=ChannelRecord("stream", 128, runtime_context),
        ),
        lp.edge(
            lp.source(validate, data=EndpointRecord("rejected", runtime_context)),
            "rejected",
            lp.sink(
                quarantine,
                data=EndpointRecord("failed record", runtime_context),
            ),
            data=ChannelRecord("durable queue", 512, runtime_context),
        ),
        lp.edge(
            lp.source(normalize, data=EndpointRecord("canonical", runtime_context)),
            "index feed",
            lp.sink(index, data=EndpointRecord("document", runtime_context)),
            data=ChannelRecord("batch", 64, runtime_context),
        ),
        lp.edge(
            lp.source(normalize, data=EndpointRecord("change log", runtime_context)),
            "normalization log",
            lp.sink(audit, data=EndpointRecord("normal changes", runtime_context)),
            data=ChannelRecord("append log", 1024, runtime_context),
        ),
        lp.edge(
            lp.source(quarantine, data=EndpointRecord("incident", runtime_context)),
            "incident log",
            lp.sink(audit, data=EndpointRecord("exceptions", runtime_context)),
            data=ChannelRecord("append log", 1024, runtime_context),
        ),
        lp.edge(
            lp.source(index, data=EndpointRecord("index version", runtime_context)),
            "search snapshot",
            lp.sink(publish, data=EndpointRecord("search input", runtime_context)),
            data=ChannelRecord("snapshot", 8, runtime_context),
        ),
        lp.edge(
            lp.source(audit, data=EndpointRecord("approval", runtime_context)),
            "release approval",
            lp.sink(publish, data=EndpointRecord("approval input", runtime_context)),
            data=ChannelRecord("event", 32, runtime_context),
        ),
        lp.edge(
            lp.source(publish, data=EndpointRecord("release", runtime_context)),
            "published documents",
            data=ChannelRecord("https", 256, runtime_context),
        ),
        name="document-intake-workflow",
        node_store=node_store.value,
    )

    if layout_algorithm.value == lp.LayoutAlgorithm.Force:
        _layouts = lp.LayoutOptions(
            algorithm=lp.LayoutAlgorithm.Force,
            direction=lp.LayoutDirection.Right,
            steps=260,
            directional_force=0.55,
            label_steps=60,
        )
    else:
        _layouts = lp.LayoutOptions(
            algorithm=lp.LayoutAlgorithm.StableLayered,
            direction=lp.LayoutDirection.Right,
            label_steps=60,
        )

    graph.render_config = lp.RenderConfig(
        title="Document intake workflow",
        layouts=_layouts,
    )
    if custom_theme.value:
        _channel_colors = {
            "append log": lp.Color("purple"),
            "https": lp.Color("blue"),
            "stream": lp.Color("teal"),
        }

        def _channel_drawing(edge):
            _paint = _channel_colors.get(edge.data.protocol, lp.Color("black"))
            return lp.EdgeDrawing(
                label=lp.TextLabel(edge.name or edge.data.protocol, fill=_paint),
                style={
                    "stroke": lp.Stroke(
                        paint=_paint,
                        thickness=lp.Length.pt(0.7),
                    )
                },
                decoration=(
                    lp.Pattern.Wave if edge.data.protocol == "stream" else lp.INHERIT
                ),
            )

        def _direction_drawing(half_edge):
            _edge = half_edge.edge
            _arrow_half = _edge.sink if _edge.sink is not None else _edge.source
            if _arrow_half is None or half_edge.index != _arrow_half.index:
                return None
            return lp.HalfEdgeDrawing(
                style={
                    "mark": lp.Mark.straight(),
                    "mark-position": lp.MarkPosition.CenterIfDangling,
                    "mark-orientation": lp.MarkOrientation.Edge,
                }
            )

        graph.render_config.selectors = lp.DrawingSelectors(
            edge=_channel_drawing,
            source=_direction_drawing,
            sink=_direction_drawing,
        )
    return graph, receive_data, request_channel, request_port, runtime_context


@app.cell
def _(custom_theme, graph, mo):
    _overlay_note = (
        "The optional theme is enabled. It is assembled entirely in Python "
        "from typed drawing selectors and inspects the workflow payloads."
        if custom_theme.value
        else "No custom theme is active; this is the default generic display."
    )
    mo.md(f"""
    ## Automatic rich display

    The next cell returns the `Graph` itself, so Marimo discovers
    `_repr_svg_()` and performs a fresh render. {_overlay_note}

    `{graph!r}`
    """)
    return


@app.cell
def _(graph):
    graph
    return


@app.cell
def _(graph, mo):
    import time

    _started_at = time.perf_counter()
    svg = graph.to_svg()
    _elapsed = time.perf_counter() - _started_at
    explicit_render = mo.vstack(
        [
            mo.md(
                f"""
                ## Explicit `to_svg()`

                Returned **{len(svg):,} characters** of SVG in
                **{_elapsed:.3f} s**. A sparse per-call `RenderConfig` can
                override this render without mutating `graph.render_config`.
                """
            ),
            mo.Html(svg),
        ]
    )
    explicit_render
    return (svg,)


@app.cell
def _(
    graph,
    mo,
    receive_data,
    request_channel,
    request_port,
    runtime_context,
    svg,
):
    _node_identity = graph.node("receive").data is receive_data
    _edge_identity = graph.edge("requests").data is request_channel
    _half_edge_identity = graph.edge("requests").sink.data is request_port
    _shared_alias = (
        graph.node("receive").data.runtime_context
        is graph.edge("requests").data.runtime_context
        is runtime_context
    )
    mo.md(f"""
    ## Arbitrary data stays in Python

    Only topology and typed drawing configuration cross the renderer boundary.
    The application objects remain attached by identity:

    - nodes / edges / half-edges: **{graph.n_nodes} / {graph.n_edges} / {graph.n_half_edges}**
    - node payload identity preserved: **{_node_identity}**
    - edge payload identity preserved: **{_edge_identity}**
    - half-edge payload identity preserved: **{_half_edge_identity}**
    - shared runtime object remains aliased: **{_shared_alias}**
    - SVG starts with: `{svg[:40]!r}`
    """)
    return


@app.cell
def _():
    return


if __name__ == "__main__":
    app.run()
