= Standalone Typst packages

These vendored packages support standalone Typst documents and documentation builds.
No Python display (graph, FeynKit or tensor-network rendering) embeds these packages
or their WASM plugins; those displays render SVG natively.

- oxifmt 1.0.0 (`preview/oxifmt/1.0.0`), MIT OR Apache-2.0
  - Source: `https://packages.typst.org/preview/oxifmt-1.0.0.tar.gz`
  - Nix recursive hash: `sha256-RtGKdyiX2kJbUjChPohSGNeYOKVlI2VM0k1uFaEqDC8=`
- MiTeX 0.2.6 (`preview/mitex/0.2.6`), Apache-2.0
  - Source: `https://packages.typst.org/preview/mitex-0.2.6.tar.gz`

oxifmt and MiTeX retain their licenses and manifests in their package directories.

CeTZ (LGPL-3.0-or-later) is not vendored. Linnest and Kurvst draw with stock CeTZ
0.5.2 from Typst Universe, which Typst reads from the local package directories;
Kurvst supplies the arrowhead geometry the former patched CeTZ fork provided.
