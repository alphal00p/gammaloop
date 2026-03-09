<div align="center">
<picture>
  <source media="(prefers-color-scheme: dark)" srcset="https://raw.githubusercontent.com/alphal00p/linnet/6ecc76560c68fa3c988e64f0e46f3009e84cca3d/assets/linnet-logo.webp">
  <img src="https://raw.githubusercontent.com/alphal00p/linnet/6ecc76560c68fa3c988e64f0e46f3009e84cca3d/assets/linnet-logo.webp" width="200" height="200">
</picture>
</div>

[![crates.io](https://img.shields.io/crates/v/linnet.svg)](https://crates.io/crates/linnet)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15494394.svg)](https://doi.org/10.5281/zenodo.15494394)

# Linnet

Linnet is a graph library, specifically designed and developed to represent **tensor networks** and **Feynman diagrams** as used in projects like **[`gammaloop`](https://github.com/alphal00p/gammaloop)** and **[`spenso`](https://github.com/alphal00p/spenso)**.

The choice of a **half-edge data structure** is central to Linnet's design. This structure is exceptionally well-suited for these applications because it intrinsically supports the definition and manipulation of subgraphs that can be cleanly and efficiently "split" along edges (so that node degree is preserved).

All of the graph algorithms, iterators and graph manipulations operate at the level of subgraphs (which can also be the whole graph).

## Key Features

- **Efficient Half-Edge Data Structure:** Utilizes a half-edge representation, enabling efficient graph traversal, easy graph modifications like:
  - node identification
  - subgraph excision (splitting graphs in two)
  - graph joining along sewn along half edges
  - edge contraction
- **Subgraph Manipulation:** Provides capabilities for defining, extracting, and analyzing various types of subgraphs, such as those induced by node properties, connectivity patterns (e.g., cycles, biconnected components), or graph cuts.
- **Graph Drawing & Visualization:** Includes functionalities for generating visual representations of graphs. Provides graph layouting based on simulated annealing of pseudo spring forces.
- **Dot Parsing** Integrates [`dot-parser`](https://codeberg.org/bromind/dot-parser/) curtesy of [Martin Vassor](mailto:martin@vassor.org?subject=[dot-parser]) integrated in a macro.

## Getting Started

To start using Linnet in your Rust project, add it as a dependency in your `Cargo.toml`:

```toml
[dependencies]
linnet = "0.17.0"
```

## AI Agent Context

For AI agents working with this codebase: **See `.agent/README.md` for comprehensive project knowledge and `.agent/notes.md` for session notes.** Key concepts:

- Half-edge data structure is central to everything
- All algorithms operate at the subgraph level
- Designed for tensor networks and Feynman diagrams
- Use `just notes` and `just knowledge` for agent context management

## Acknowledgement

This crate was written by [Lucien Huber](https://github.com/lcnbr/)
