# Clinnet

`clinnet` is Linnet's command-line companion for recursively rendering DOT files with Typst and
assembling the figures into a grid PDF. The Cargo package is named `clinnet`; its executable is
named `linnet`.

The canonical [Clinnet usage guide](https://alphal00p.github.io/gammaloop/products/linnet/latest/guides/clinnet/) covers
installation, default outputs, template inputs, cache invalidation, and targeted rebuilds. See
the [Clinnet version history](https://alphal00p.github.io/gammaloop/products/linnet/latest/version-history/clinnet/) for
package-specific changes and `linnet --help` for the exact installed command surface.

```bash
cargo install clinnet
cargo install typst-cli --version 0.15.0 --locked
linnet draw examples
```

Contributor policy and documentation checks are in [`CONTRIBUTING.typ`](../../CONTRIBUTING.typ).
