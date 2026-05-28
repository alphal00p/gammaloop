# Contributing

## Typst Graph Debugging

When debugging Linnest layout output from Typst, expose the graph records through
Typst metadata and inspect them with `typst query`. The metadata element must be
wrapped in markup, because Typst labels can only be attached in markup mode.

```typst
[#metadata(graph.nodes(g)) <linnest-tree-nodes>]
[#metadata(graph.edges(g)) <linnest-tree-edges>]
```

Then query the values from the same root used for compilation:

```bash
typst query --root crates crates/linnest/typst/examples/tree.typ '<linnest-tree-nodes>' --field value --pretty
typst query --root crates crates/linnest/typst/examples/tree.typ '<linnest-tree-edges>' --field value --pretty
```

This is useful for checking generated `pos`, `label-pos`, `layout-width`,
`layout-height`, `label-width`, and `label-height` values without inferring them
from the rendered PDF.
