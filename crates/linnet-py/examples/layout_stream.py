# /// script
# requires-python = ">=3.10"
# dependencies = [
#     "anywidget==0.9.18",
#     "linnet-py==0.1.0",
#     "marimo==0.24.0",
#     "typst==0.15.0",
# ]
# ///

# ruff: noqa: B018, PLR1711  # Cell outputs and empty returns are Marimo syntax.

import marimo

__generated_with = "0.24.0"
app = marimo.App(width="medium")


@app.cell(hide_code=True)
def _():
    import asyncio
    import sys

    import anywidget
    import linnet_py as lp
    import marimo as mo
    import traitlets

    return anywidget, asyncio, lp, mo, sys, traitlets


@app.cell(hide_code=True)
def _(mo):
    mo.md("""
    # Watch a layout settle

    Edit a graph, open the force controls, and watch the solver work. Each frame
    advances the same Rust simulation; the SVG keeps its nodes and paths and
    updates their coordinates. **Pause** stops the solver between batches.

    This is a geometry preview: blue circles are nodes; small amber points are
    the solver's edge positions. Typst typography and final label placement are
    left to the finished drawing.
    """)
    return


@app.cell(hide_code=True)
def _(mo):
    example = mo.ui.dropdown(
        options={
            "X-box core (K4)": """digraph Xbox {
  a; b; c; d;
  c -> d;
  a -> b;
  c -> a;
  d -> b;
  a -> d;
  b -> c;
}""",
            "Loops and open legs": """digraph Open {
  in [style=invis]; out [style=invis];
  a; b; c;
  in -> a;
  a -> b; a -> b;
  b -> c; c -> a;
  b -> b;
  c -> out;
}""",
        },
        value="X-box core (K4)",
        label="Example",
    )
    restart = mo.ui.button(value=0, on_click=lambda value: value + 1, label="Restart")
    mo.hstack([example, restart], justify="space-between", align="end")
    return example, restart


@app.cell(hide_code=True)
def _(example, mo):
    dot_editor = mo.ui.code_editor(
        value=example.value,
        language="text",
        min_height=220,
        max_height=360,
        debounce=600,
        label="DOT editor · edits restart the layout",
    )
    dot_editor
    return (dot_editor,)


@app.cell(hide_code=True)
def _(mo):
    steps = mo.ui.slider(
        20,
        600,
        step=20,
        value=200,
        show_value=True,
        debounce=True,
        label="Steps per epoch",
    )
    epochs = mo.ui.slider(
        1, 24, value=8, show_value=True, debounce=True, label="Epochs"
    )
    every = mo.ui.slider(
        1, 24, value=4, show_value=True, debounce=True, label="Steps per frame"
    )
    seed = mo.ui.slider(0, 100, value=1, show_value=True, debounce=True, label="Seed")
    step = mo.ui.slider(
        0.005,
        0.2,
        step=0.005,
        value=0.02,
        show_value=True,
        debounce=True,
        label="Step size",
    )
    cool = mo.ui.slider(
        0.5, 1.0, step=0.01, value=0.85, show_value=True, debounce=True, label="Cooling"
    )
    spring_strength = mo.ui.slider(
        0.1,
        20.0,
        step=0.1,
        value=1.0,
        show_value=True,
        debounce=True,
        label="Spring strength",
    )
    repulsion = mo.ui.slider(
        0.0,
        20.0,
        step=0.1,
        value=1.5,
        show_value=True,
        debounce=True,
        label="Repulsion",
    )
    length_scale = mo.ui.slider(
        0.1,
        4.0,
        step=0.1,
        value=1.0,
        show_value=True,
        debounce=True,
        label="Spring length scale",
    )
    depth_scale = mo.ui.slider(
        0.0,
        3.0,
        step=0.1,
        value=1.0,
        show_value=True,
        debounce=True,
        label="Initial depth scale",
    )
    flattening_end = mo.ui.slider(
        0.05,
        1.0,
        step=0.05,
        value=0.5,
        show_value=True,
        debounce=True,
        label="Flattening end",
    )
    mo.accordion(
        {
            "Force controls · change a value to restart": mo.hstack(
                [
                    mo.vstack(
                        [
                            spring_strength,
                            repulsion,
                            length_scale,
                            depth_scale,
                            flattening_end,
                        ]
                    ),
                    mo.vstack([steps, epochs, step, cool, every, seed]),
                ],
                widths="equal",
                gap=2,
            )
        }
    )
    return (
        cool,
        depth_scale,
        epochs,
        every,
        flattening_end,
        length_scale,
        repulsion,
        seed,
        spring_strength,
        step,
        steps,
    )


@app.cell(hide_code=True)
def _(anywidget, traitlets):
    class LayoutView(anywidget.AnyWidget):
        """Retain graph topology in the browser and synchronize position frames."""

        topology = traitlets.Dict().tag(sync=True)
        frame = traitlets.Dict().tag(sync=True)
        paused = traitlets.Bool(False).tag(sync=True)
        _esm = r"""
        export default {
          render({model, el}) {
            const ns = "http://www.w3.org/2000/svg";
            const make = (tag, attrs = {}) => {
              const element = document.createElementNS(ns, tag);
              for (const [key, value] of Object.entries(attrs)) element.setAttribute(key, value);
              return element;
            };
            el.innerHTML = '<div data-linnet-render-ready="stream"></div>';
            const root = el.firstElementChild;
            root.style.cssText = "border:1px solid #cbd5e1;border-radius:12px;overflow:hidden;background:#f8fafc;color:#0f172a";
            const bar = document.createElement("div");
            bar.style.cssText = "display:flex;gap:16px;align-items:center;padding:12px 16px;border-bottom:1px solid #e2e8f0;font:13px system-ui";
            const pause = document.createElement("button");
            pause.textContent = "Pause";
            pause.style.cssText = "padding:5px 14px;border:1px solid #94a3b8;border-radius:6px;background:white;color:#0f172a;cursor:pointer";
            const status = document.createElement("span");
            status.setAttribute("role", "status");
            bar.append(pause, status);
            const svg = make("svg", {viewBox:"0 0 720 420", role:"img", "aria-label":"Streaming force layout"});
            svg.style.cssText = "display:block;width:100%;height:auto;min-height:260px";
            const graph = make("g");
            svg.append(graph);
            root.append(bar, svg);
            const topology = model.get("topology");
            const paths = topology.endpoints.map(() => {
              const path = make("path", {fill:"none", stroke:"#64748b", "stroke-width":1.6});
              graph.append(path);
              return path;
            });
            const edgePoints = topology.endpoints.map(() => {
              const point = make("circle", {r:3, fill:"#d97706"});
              graph.append(point);
              return point;
            });
            const nodes = topology.node_names.map((name) => {
              const group = make("g");
              const circle = make("circle", {r:8, fill:"#2563eb", stroke:"white", "stroke-width":2});
              const label = make("text", {y:25, "text-anchor":"middle", fill:"#0f172a", "font-size":13, "font-family":"system-ui"});
              label.textContent = name;
              group.append(circle, label);
              graph.append(group);
              return group;
            });
            let bounds = null;
            const updatePause = () => {
              pause.textContent = model.get("paused") ? "Resume" : "Pause";
              pause.setAttribute("aria-pressed", String(model.get("paused")));
            };
            const onPause = () => {
              model.set("paused", !model.get("paused"));
              model.save_changes();
            };
            const update = () => {
              const frame = model.get("frame");
              if (!frame.nodes) return;
              root.dataset.iteration = String(frame.iteration);
              root.dataset.done = String(frame.done);
              pause.disabled = frame.done;
              status.textContent = `${frame.done ? "Finished" : "Solving"} · step ${frame.iteration} · Δ ${frame.max_movement.toPrecision(3)} · ${nodes.length} nodes / ${paths.length} edges`;
              const points = [...frame.nodes, ...frame.edges];
              if (!points.length) { status.textContent = "Empty graph · no positions to solve"; return; }
              const xs = points.map(p => p[0]), ys = points.map(p => p[1]);
              const extent = [Math.min(...xs), Math.max(...xs), Math.min(...ys), Math.max(...ys)];
              if (!bounds) bounds = extent;
              else bounds = [Math.min(bounds[0], extent[0]), Math.max(bounds[1], extent[1]), Math.min(bounds[2], extent[2]), Math.max(bounds[3], extent[3])];
              const cx = (bounds[0] + bounds[1]) / 2, cy = (bounds[2] + bounds[3]) / 2;
              const scale = Math.min(580 / Math.max(bounds[1] - bounds[0], 1), 300 / Math.max(bounds[3] - bounds[2], 1));
              const screen = p => [360 + (p[0] - cx) * scale, 200 - (p[1] - cy) * scale];
              const n = frame.nodes.map(screen), e = frame.edges.map(screen);
              nodes.forEach((node, i) => node.setAttribute("transform", `translate(${n[i]})`));
              topology.endpoints.forEach(([source, sink], i) => {
                const p = e[i];
                const a = source === null ? p : n[source];
                const b = sink === null ? p : n[sink];
                let path;
                if (source !== null && source === sink) {
                  const dx = (p[0] - a[0]) * 4 / 3, dy = (p[1] - a[1]) * 4 / 3;
                  path = `M${a} C${a[0]+dx-dy},${a[1]+dy+dx} ${a[0]+dx+dy},${a[1]+dy-dx} ${a}`;
                } else if (source === null || sink === null) {
                  path = `M${a} L${b}`;
                } else {
                  path = `M${a} Q${2*p[0]-(a[0]+b[0])/2},${2*p[1]-(a[1]+b[1])/2} ${b}`;
                }
                paths[i].setAttribute("d", path);
                edgePoints[i].setAttribute("cx", p[0]);
                edgePoints[i].setAttribute("cy", p[1]);
              });
            };
            pause.addEventListener("click", onPause);
            model.on("change:frame", update);
            model.on("change:paused", updatePause);
            update(); updatePause();
            return () => {
              model.off("change:frame", update);
              model.off("change:paused", updatePause);
              pause.removeEventListener("click", onPause);
            };
          }
        };
        """

    return (LayoutView,)


@app.cell(hide_code=True)
def _(
    LayoutView,
    asyncio,
    cool,
    depth_scale,
    dot_editor,
    epochs,
    every,
    flattening_end,
    length_scale,
    lp,
    mo,
    repulsion,
    restart,
    seed,
    spring_strength,
    step,
    steps,
    sys,
):
    _restart = restart.value
    _source = dot_editor.value
    _options = {
        "every": every.value,
        "steps": steps.value,
        "epochs": epochs.value,
        "seed": seed.value,
        "step": step.value,
        "cool": cool.value,
        "spring_strength": spring_strength.value,
        "repulsion": repulsion.value,
        "length_scale": length_scale.value,
        "depth_scale": depth_scale.value,
        "flattening_end": flattening_end.value,
    }

    async def _animate(source=_source, options=_options):
        _thread = mo.current_thread()
        try:
            _stream = lp.LayoutStream.from_dot(source, **options)
            _view = LayoutView(
                topology={
                    "node_names": _stream.node_names,
                    "endpoints": _stream.endpoints,
                }
            )
            mo.output.replace(mo.ui.anywidget(_view))
            while not _thread.should_exit:
                if not _view.paused:
                    _frame = next(_stream, None)
                    if _frame is None:
                        break
                    _view.frame = {
                        "nodes": _frame.nodes,
                        "edges": _frame.edges,
                        "iteration": _frame.iteration,
                        "done": _frame.done,
                        "max_movement": _frame.max_movement,
                    }
                # Give the browser event loop a chance to paint each frame.
                await asyncio.sleep(1 / 30)
        except (TypeError, ValueError, RuntimeError) as _error:
            mo.output.replace(mo.callout(f"DOT/layout error: {_error}", kind="danger"))

    # Marimo manages cancellation when controls change. Pyodide runs an async
    # thread target cooperatively; native Python needs its own event loop.
    _target = (
        _animate if sys.platform == "emscripten" else lambda: asyncio.run(_animate())
    )
    mo.Thread(target=_target, daemon=True).start()
    return


if __name__ == "__main__":
    app.run()
