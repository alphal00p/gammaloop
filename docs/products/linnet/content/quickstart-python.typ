#import "../../shared.typ": callout

#let quickstart-python = [
= Using Linnet from Python

The `linnet_py` extension provides DOT-backed Linnet graphs to Python 3.10 and newer. It is a real
binding with runtime tests, but it is not currently published to PyPI.

#callout("Developer preview: build from source", [
  There is no supported `pip install linnet-py` release yet. The similarly named `linnet` project
  on PyPI is unrelated. Use this page from a GammaLoop checkout until official `linnet-py` wheels
  are published.
])

== Build the extension

Install a stable Rust toolchain and Python 3.10 or newer, then run:

// docs-example: syntax
```sh
git clone https://github.com/alphal00p/gammaloop.git
cd gammaloop
python3.10 -m venv .venv
. .venv/bin/activate
python -m pip install "maturin>=1.7,<2.0"
maturin develop --locked \
  --manifest-path crates/linnet-py/Cargo.toml \
  --features extension-module,abi3-py310
```

== Parse and inspect a graph

Save this as `linnet_quickstart.py`:

// docs-example: compile linnet-python-quickstart
```python
import linnet_py as lp

graph = lp.DotGraph.from_string(
    """
    digraph G {
      A;
      B;
      A -> B;
    }
    """
)

whole = graph.full_filter()
assert graph.n_nodes() == 2
assert graph.n_edges() == 1
assert len(graph.iter_nodes_of(whole)) == 2
assert len(graph.iter_edges_of(whole)) == 1

print(graph.dot())
```

Run `python linnet_quickstart.py`. Success means the assertions pass and the graph is printed as
DOT. `DotGraph` also exposes typed nodes, half-edges, subgraphs, cycles, oriented cuts, and
traversal trees; use the #link("reference/python/")[Python reference] for the exact current
surface.

Use the #link("quickstart/rust/")[Rust guide] when you need Linnet's complete crate API, or
the #link("quickstart/typst/")[Typst guide] when the final result should be a drawing.
]
