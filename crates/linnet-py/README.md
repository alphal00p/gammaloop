# linnet-py

`linnet-py` exposes Linnet's native half-edge graph, node-store variants,
topology editing, subgraphs, algorithms, arbitrary Python element data, and
typed drawing configuration as a standalone Python package.

Rendering is in-process through `typst` 0.15.0. The wheel embeds Linnest,
Kurvst, CeTZ, and oxifmt, while the source distribution carries their build
inputs. Rendering converts graph topology to a versioned CBOR graph spec and
never serializes arbitrary Python `.data`. Clinnet and a Typst executable are
not runtime dependencies.

The linnet-py, Linnest, and Kurvst sources are MIT-licensed under `LICENSE`;
the distribution also carries the licenses and provenance of its vendored
Typst packages.

```python
import linnet_py as lp

left = lp.node("left", data=object(), label="Left")
right = lp.node("right", data=object(), label="Right")
graph = lp.build(
    left,
    right,
    lp.edge(lp.source(left), "dependency", lp.sink(right)),
)

svg = graph.to_svg()
graph.render("graph.pdf")
```

Generic rendering labels nodes and edges by their structural IDs ($n_i$, $e_i$,
...). Set `lp.DrawOptions(show_half_edge_ids=True)` to add optional endpoint
IDs ($h_i$); explicit drawing labels take precedence.

The package is currently built from the GammaLoop workspace. Its complete API
and development instructions are maintained in the Linnet documentation there.

## Marimo examples

The generic and physics notebooks under `examples/` carry PEP 723 dependencies
for native sandboxing and browser installation. Native editable sessions still
use the workspace package directly:

```console
uvx --from marimo==0.24.0 --with-editable crates/linnet-py \
  marimo edit crates/linnet-py/examples/rendering_api.py
```

Once an Emscripten wheel is available, export both notebooks as editable static
WASM pages without changing their checked-in dependency metadata:

```console
uv run --with marimo==0.24.0 \
  python crates/linnet-py/examples/export_wasm.py \
  --wheel dist/linnet_py-0.1.0-cp310-abi3-pyemscripten_2026_0_wasm32.whl \
  --output dist/linnet-wasm
```

The helper stages the local wheel override temporarily, runs Marimo's strict
`MW` checks, exports in edit mode, and serves the result for an HTTP smoke test.
Omit `--wheel` after publishing the browser wheel. Pass `--browser-smoke` in an
environment with Playwright and Chromium to wait for a real SVG render from
each Pyodide notebook.

The exported pages embed their Python source and open with editable code cells;
browser edits do not modify the checked-in notebooks. Serve the output directory
over HTTP (browsers cannot launch Pyodide from `file://`), for example:

```console
python -m http.server --directory dist/linnet-wasm
```
