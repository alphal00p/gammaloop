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

The package is currently built from the GammaLoop workspace. Its complete API
and development instructions are maintained in the Linnet documentation there.
