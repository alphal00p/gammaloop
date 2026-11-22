"""Time one native execution of every enabled cell in an example Marimo notebook."""

import argparse
import ast
import hashlib
import importlib.util
import json
import sys
import threading
from collections import Counter
from contextlib import ExitStack, contextmanager
from pathlib import Path
from time import perf_counter_ns
from unittest.mock import patch

from marimo._ast.app import AppScriptRunner
from marimo._ast.names import SETUP_CELL_NAME, is_internal_cell_name


def load_notebook(path, setup_seconds):
    # app.run() schedules normal cells, but Python executes the setup cell
    # during import. Add a timer inside its with statement, keeping the original
    # filename/line numbers for Marimo's inspect-based cell registration.
    tree = ast.parse(path.read_bytes(), filename=str(path))
    for statement in tree.body:
        if not isinstance(statement, ast.With):
            continue
        for item in statement.items:
            expr = item.context_expr
            target = expr.func if isinstance(expr, ast.Call) else expr
            if ast.unparse(target) == "app.setup":
                timer = ast.parse("with __measure_setup__(): pass").body[0].items[0]
                statement.items.append(timer)
                break

    @contextmanager
    def time_setup():
        start = perf_counter_ns()
        yield
        setup_seconds[SETUP_CELL_NAME] = (perf_counter_ns() - start) / 1e9

    spec = importlib.util.spec_from_file_location("__timed_notebook__", path)
    module = importlib.util.module_from_spec(spec)
    module.__measure_setup__ = time_setup
    sys.modules[spec.name] = module
    code = compile(ast.fix_missing_locations(tree), str(path), "exec")
    exec(code, module.__dict__)  # noqa: S102 - execute the instrumented notebook module
    return module.app


def cell_identities(cells):
    # Marimo's generated four-letter IDs depend on registration order.
    # Named cells and unique definitions survive implementation edits. For
    # displays with the same inputs, compare the same normalized expression.
    definitions = {name for entry in cells for name in entry.cell.defs}
    definitions = {name for name in definitions if not name.startswith("_")}
    keys = {}
    labels = {}
    for entry in cells:
        outputs = sorted(entry.cell.defs & definitions)
        inputs = sorted(entry.cell.refs & definitions)
        if entry.cell_id == SETUP_CELL_NAME:
            key, label = "setup", "Notebook setup"
        elif not is_internal_cell_name(entry.name):
            key, label = "named." + entry.name, entry.name
        elif outputs:
            key, label = "defines." + "+".join(outputs), ", ".join(outputs)
        else:
            key = "uses." + "+".join(inputs) if inputs else "display"
            label = (
                entry.code.strip().splitlines()[0][:100]
                if entry.code.strip()
                else "Empty cell"
            )
        keys[entry.cell_id] = key
        labels[entry.cell_id] = label
    repeats = Counter(keys.values())
    for entry in cells:
        key = keys[entry.cell_id]
        if repeats[key] > 1:
            syntax = ast.dump(ast.parse(entry.code))
            keys[entry.cell_id] = (
                key + "." + hashlib.sha256(syntax.encode()).hexdigest()[:16]
            )
    if len(set(keys.values())) != len(cells):
        raise ValueError(
            "Ambiguous identical cells; give these cells distinct Marimo names"
        )
    return keys, {key: labels[cid] for cid, key in keys.items()}


def measure(path):
    path = path.resolve()
    sys.path.insert(0, str(path.parent))
    seconds = {}
    app = load_notebook(path, seconds)
    app._maybe_initialize()
    enabled = [
        entry
        for entry in app._cell_manager.cell_data()
        if not app._graph.is_disabled(entry.cell_id)
    ]
    if not enabled:
        raise ValueError(f"No enabled cells in {path}")
    keys, labels = cell_identities(enabled)
    scenario, control = {
        "examples/notebooks/feynkit/08_numerical_cff_ltd_marimo.py": (
            "submit-integration-defaults",
            "integration_settings",
        ),
        "crates/linnet-py/examples/physics_render_settings.py": (
            "export-vector-svg",
            "vector_export",
        ),
        "crates/linnet-py/examples/rendering_api.py": (
            "capture-chromium-svg",
            "snapshot_request",
        ),
    }.get("/".join(path.parts[-4:]), ("default", None))
    submitted = {}
    layout_stream = path.parts[-4:] == (
        "crates",
        "linnet-py",
        "examples",
        "layout_stream.py",
    )
    widgets = []
    thread_errors = []

    class Timer:
        name = "notebook-cell-timing"

        def setup(self, cell, glbls):
            if not submitted and control in glbls:
                # This is the next cell's context: Marimo rejects reading
                # a widget in the cell that defines it. Submit before timing.
                widget = glbls[control]
                if control == "integration_settings":
                    widget._update(widget.element._initial_value_frontend)
                    submitted.update(widget.value)
                else:
                    if control == "snapshot_request":
                        # The real widget receives its font metrics before the
                        # user can capture its SVG. Populate the cache via the
                        # view's native headless Chromium implementation.
                        glbls["view"].measured(browser="chromium")
                    widget._update(1)
                    submitted[control] = widget.value
            self.start = perf_counter_ns()

        def teardown(self, cell, glbls, run_result):
            elapsed = perf_counter_ns() - self.start
            if run_result.exception is not None:
                # Marimo normally allows mo.stop(). It is not a completed
                # workload; even a stopped leaf must fail the coverage check.
                return
            cid = cell.cell_id
            if cid in seconds or elapsed <= 0:
                raise RuntimeError(
                    f"Invalid execution count or duration for {keys[cid]}"
                )
            seconds[cid] = elapsed / 1e9
            if layout_stream and "LayoutView" in cell.refs:
                # Keep the actual Marimo cell context alive for its worker.
                # This cell times dispatch; joining/verification is outside
                # the interval and checks the terminal output of the animation.
                for thread in threading.enumerate():
                    if isinstance(thread, glbls["mo"].Thread):
                        thread.join()
                if thread_errors:
                    raise RuntimeError(
                        "Layout animation thread failed"
                    ) from thread_errors[0].exc_value
                frames = [
                    w.frame for w in widgets if isinstance(w, glbls["LayoutView"])
                ]
                if len(frames) != 1 or not frames[0].get("done"):
                    raise RuntimeError(
                        "Layout animation did not produce its final frame"
                    )
                submitted["final_frame"] = True

    class TimedRunner(AppScriptRunner):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            self._evaluator.lifecycles.append(Timer())

    # The pinned evaluator calls its lifecycles around both sync and awaited
    # cell bodies. Leave scheduling, output evaluation and errors to Marimo.
    with ExitStack() as stack:
        stack.enter_context(patch("marimo._ast.app.AppScriptRunner", TimedRunner))
        if layout_stream:
            from ipywidgets import Widget

            stack.enter_context(
                patch.object(threading, "excepthook", thread_errors.append)
            )
            Widget.on_widget_constructed(widgets.append)
            stack.callback(Widget.on_widget_constructed, None)
        app.run()
    missing = keys.keys() - seconds.keys()
    if missing:
        raise RuntimeError(
            f"Enabled cells did not finish in {path}: {[keys[c] for c in sorted(missing)]}"
        )
    if control and not submitted:
        raise RuntimeError(f"Scenario {scenario} never activated {control}")
    if layout_stream and not submitted:
        raise RuntimeError("Missing layout animation cell")
    return {
        "cells": {key: seconds[cid] for cid, key in keys.items()},
        "labels": labels,
        "scenario": "dispatch-layout-animation" if layout_stream else scenario,
        "settings": submitted,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--notebook", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    result = measure(args.notebook)
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")


if __name__ == "__main__":
    main()
