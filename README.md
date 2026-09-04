<div align="center">
<picture>
  <source media="(prefers-color-scheme: dark)" srcset="https://github.com/alphal00p/gammaloop/blob/2ee2ec575fa575c26bdaf89a3e7df41428b879dc/assets/gammalooplogo-dark.svg">
  <img src="https://github.com/alphal00p/gammaloop/blob/2ee2ec575fa575c26bdaf89a3e7df41428b879dc/assets/gammalooplogo-light.svg" width="300">
</picture>

[![Nix-CI](https://nix-ci.com/badge/gh:alphal00p:gammaloop/main?v=2)](https://nix-ci.com/gh:alphal00p:gammaloop/main)
[![GitHub Actions](https://github.com/alphal00p/gammaloop/actions/workflows/continuous-integration.yml/badge.svg?branch=main)](https://github.com/alphal00p/gammaloop/actions/workflows/continuous-integration.yml)
<!--[![crates.io](https://img.shields.io/crates/v/spenso.svg)](https://crates.io/crates/spenso)
[![Build Status](https://github.com/alphal00p/spenso/actions/workflows/ci.yml/badge.svg)](https://github.com/alphal00p/spenso/actions/workflows/ci.yml)
[![codecov](https://codecov.io/github/alphal00p/spenso/graph/badge.svg?token=ST0XA54QSF)](https://codecov.io/github/alphal00p/spenso)
[![Ask DeepWiki](https://deepwiki.com/badge.svg)](https://deepwiki.com/alphal00p/spenso)-->
</div>

# GammaLoop

GammaLoop computes differential collider cross-sections with Local Unitarity. It combines a
stateful command-line application with Rust and Python APIs and the Linnet, Spenso, Idenso, and
Vakint research libraries.

The [live αLoop documentation](https://alphal00p.github.io/gammaloop/) is the user entry point.
Its canonical sources are Typst files under [`docs/products`](docs/products); this README is only
the compatibility index shown by GitHub and package tooling.

## Start here

- [Build, create, inspect, and resume a GammaLoop state](https://alphal00p.github.io/gammaloop/products/gammaloop/latest/tutorial/)
- [Generate processes and manage integration workspaces](https://alphal00p.github.io/gammaloop/products/gammaloop/latest/guides/process-generation/)
- [Inspect events, selectors, and observables from Python or Rust](https://alphal00p.github.io/gammaloop/products/gammaloop/latest/guides/events-and-observables/)
- [Choose the CLI, Rust, or Python interface and enable shell completion](https://alphal00p.github.io/gammaloop/products/gammaloop/latest/reference/interfaces/)
- [Browse the generated API and CLI reference](https://alphal00p.github.io/gammaloop/products/gammaloop/latest/reference/)

For a source checkout, the supported environment and a quick smoke test are:

```bash
nix develop
just build-cli
./gammaloop --help
```

The tutorial explains non-Nix prerequisites, the repository wrapper, state paths, run-card
replay, and the cost of the maintained scientific examples.

## Related products

- [Linnet](https://alphal00p.github.io/gammaloop/products/linnet/latest/): half-edge graphs, Linnest layout, and Clinnet rendering
- [Spenso](https://alphal00p.github.io/gammaloop/products/spenso/latest/): typed tensors and executable tensor networks
- [Idenso](https://alphal00p.github.io/gammaloop/products/idenso/latest/): symbolic tensor identities and rewrites
- [Vakint](https://alphal00p.github.io/gammaloop/products/vakint/latest/): vacuum-integral matching and evaluation

See [www.alphaloop.ch](https://www.alphaloop.ch) for the broader project and literature.

## Contributing

Read [`CONTRIBUTING.typ`](CONTRIBUTING.typ) before editing. Typst is the canonical format for
authored documentation; files named exactly `README.md` and `AGENTS.md` are compatibility
exceptions and should remain concise indexes.
