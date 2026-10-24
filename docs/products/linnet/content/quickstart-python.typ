#import "../../shared.typ": callout

#let quickstart-python = [
= Using Linnet from Python

Symbolica Community exposes native Linnet graphs as `symbolica.community.graph`.
Build and install the complete Community wheel, or use the repository's notebook host:

// docs-example: syntax
```sh
python -m pip install "maturin>=1.13,<2"
maturin develop --locked --manifest-path examples/notebooks/symbolica-host/Cargo.toml
```

Graph rendering uses native Rust geometry and direct MathML labels. Notebook
display requires `anywidget`; the live browser measures labels before the
drawing is finalized. No Typst compiler is used for this path.

The GammaLoop package and this checkout's isolated notebook host declare
`notebook-display` (`anywidget>=0.9,<1`), `graph-export`
(`playwright>=1.55,<2`) and `typst-export` (`typst>=0.15`) extras. Do not assume
those extras exist in the currently released external Community wheel.
For an existing host environment, install the dependencies you need directly:

// docs-example: syntax
```sh
python -m pip install "anywidget>=0.9,<1"
```

== Parse and inspect a graph

Save this as `linnet_quickstart.py`:

// docs-example: compile linnet-python-quickstart
```python
from symbolica.community import graph as lp

codec = lp.DotCodec.topology()
graph = lp.Graph.from_dot(
    """
    digraph G {
      A;
      B;
      A -> B;
    }
    """,
    codec,
)

whole = graph.full_subgraph()
assert graph.n_nodes == 2
assert graph.n_edges == 1
assert len(graph.nodes_of(whole)) == 2
assert len(graph.edges_of(whole)) == 1

print(graph.to_dot())
```

#context if target() == "html" {
  html.elem("div", attrs: (
    class: "live-notebook",
    "data-notebook": "python_quickstart",
    "aria-label": "Run the Python quickstart",
  ))[
    #html.elem("p", attrs: (class: "live-notebook-fallback"))[
      This editable example loads automatically when it comes into view. It requires JavaScript
      and downloads Python and its dependencies on the first visit. You can also run the
      example above locally.
    ]
  ]
}

Run `python linnet_quickstart.py`. Success means the assertions pass and the graph is printed as
DOT. `Graph` also exposes typed nodes, half-edges, subgraphs, cycles, oriented cuts, and
traversal trees; use the #link("reference/python/")[Python reference] for the exact current
surface.

In a notebook, keep and display the configured view:

```python
from IPython.display import display

view = graph.render()
display(view)
```

After the browser finishes measuring, `snapshot = view.measured()` returns an
immutable drawing. `snapshot.to_svg()`, `snapshot.to_html()` and
`snapshot.save("graph.svg")` preserve its measured boxes and embedded font.
Exporting a pending view raises an error.

Scripts may explicitly request `view.measured(browser="chromium")`. Install
`playwright>=1.55,<2` and run `python -m playwright install chromium` first.
The default SVG contains browser MathML in XHTML `foreignObject` elements;
use `view.outlined()` with optional Python typst-py (distribution `typst>=0.15`)
for vector labels before PDF, PNG or Typst SVG import. Neither an automatic
headless browser nor estimated label bounds are used.

Continue with #link("guides/python-rendering/")[Rendering graphs from Python] to explore a
complete graph with custom drawing selectors. The #link("playground/")[live Python playground]
also lets you edit DOT and tune the layout while watching the solver run.

Use the #link("quickstart/rust/")[Rust guide] when you need Linnet's complete crate API, or
the #link("quickstart/typst/")[Typst guide] when the final result should be a drawing.
]
