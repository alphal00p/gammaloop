# Linnet Python bindings

Native Linnet half-edge graphs for Python: builders, live graph views,
subgraphs, DOT codecs, graph algorithms, layout streams, and rendering with an
embedded Typst compiler. No separate Typst package or executable is required.

The same bindings ship in two distributions:

| Distribution | Import | Graph keys |
| --- | --- | --- |
| `linnet-py` | `import linnet_py` | Python values |
| Symbolica Community wheel | `from symbolica.community import graph` | Python values and Symbolica expressions |

Each build names its classes after its own module, so reprs and pickles never
refer to the other distribution. Only the Community build accepts Symbolica
expressions as keys of their own kind and shares its rendering types with
HepKit diagrams and tensor networks.

```python
import linnet_py as lp

graph = lp.Graph(edges=[lp.Edge(lp.Node("a"), lp.Node("b"), name="ab")])
svg = graph.to_svg()
```

See the [Linnet documentation](https://alphal00p.github.io/gammaloop/linnet/)
for construction, topology, algorithms, DOT, layout streaming, and rendering.

## Licenses

Linnet, Linnest, and Kurvst are MIT-licensed (`LICENSE`). The wheel embeds
Typst and its fonts (`licenses/`) and the CeTZ, oxifmt, and MiTeX Typst
packages, whose licenses and provenance are under `vendor/typst-packages/`.
The wheel's CycloneDX SBOM lists every compiled Rust crate with its license,
including the GPL-2.0-or-later `dot-parser`.

## Development

This crate is a Rust library. Its `standalone` feature adds the `linnet_py`
extension module, which only the `linnet-py` wheel enables. Its `community`
feature instead names the classes after `symbolica.community.graph` for a
Symbolica Community build, which registers `GraphModule`; the two features
cannot be combined.

Build and test the standalone wheel from this directory:

```console
maturin build --release --out dist
python -m pip install dist/linnet_py-*.whl
python -m unittest discover -s tests -p 'test_*.py'
```

The same suite runs against the Community wheel built from
`examples/notebooks/symbolica-host`. It tests whichever module is installed;
set `LINNET_PY_MODULE` to `linnet_py` or `symbolica.community.graph` when both
are. `tests/test_community_rendering.py` needs the Community wheel.

`linnet_py.pyi` is the standalone stub and
`python/symbolica/community/graph/__init__.pyi` the Community stub; they differ
only in the graph key type. `docs/api/python/README.md` describes how to
regenerate and check both.
