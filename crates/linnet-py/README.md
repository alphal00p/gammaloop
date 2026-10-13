# Linnet Python bindings

Native Linnet half-edge graphs for Python: builders, live graph views,
subgraphs, DOT codecs, graph algorithms, layout streams, and rendering with an
embedded Typst compiler; no Typst executable is used. Native drawings, the
default, need no Typst package. The full Typst route (templates, selectors,
layout sequences) imports stock CeTZ 0.5.2 from the local Typst package path
or cache: compile any document that imports `@preview/cetz:0.5.2` once with
the Typst CLI, or set `TYPST_PACKAGE_PATH` to a directory holding
`preview/cetz/0.5.2`.

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

Linnet, Linnest, and Kurvst are MIT-licensed (`LICENSE`); Linnest's layout seed
includes code translated from NetworkX (BSD-3-Clause,
`licenses/networkx-LICENSE`), and Kurvst's mark geometry follows tiptoe (MIT,
`licenses/tiptoe-LICENSE`). The wheel embeds Typst with the permissively
licensed Libertinus Serif, STIX Two Math, and DejaVu Sans Mono fonts
(`licenses/`), and the oxifmt and MiTeX Typst packages, whose licenses and
provenance are under `vendor/typst-packages/`. It does not include CeTZ
(LGPL-3.0): native drawings need no Typst package, and the Typst route reads
stock CeTZ 0.5.2 from the local Typst package path or cache. The wheel's
CycloneDX SBOM lists the Rust crates of its dependency graph with their
licenses.

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

The Typst-route tests need stock CeTZ 0.5.2. The development shell provides it
through `TYPST_PACKAGE_CACHE_PATH`; elsewhere run `scripts/fetch-stock-cetz.sh`
from the repository root and point `TYPST_PACKAGE_PATH` at the absolute path of
its `target/typst-packages`.

The same suite runs against the Community wheel built from
`examples/notebooks/symbolica-host`. It tests whichever module is installed;
set `LINNET_PY_MODULE` to `linnet_py` or `symbolica.community.graph` when both
are. `tests/test_community_rendering.py` needs the Community wheel.

`linnet_py.pyi` is the standalone stub and
`python/symbolica/community/graph/__init__.pyi` the Community stub; they differ
only in the graph key type. `docs/api/python/README.md` describes how to
regenerate and check both.
