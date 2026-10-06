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

Rendering uses the embedded Rust Typst compiler. No separate Python Linnet or Typst
package or external Typst executable is required.

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

Continue with #link("guides/python-rendering/")[Rendering graphs from Python] to explore a
complete graph with custom drawing selectors. The #link("playground/")[live Python playground]
also lets you edit DOT and tune the layout while watching the solver run.

Use the #link("quickstart/rust/")[Rust guide] when you need Linnet's complete crate API, or
the #link("quickstart/typst/")[Typst guide] when the final result should be a drawing.
]
