= Standalone Typst packages

These vendored packages support standalone Typst documents and documentation builds.
The Python graph renderer does not embed these packages or their WASM plugins.

- Stock CeTZ 0.5.2, LGPL-3.0-or-later, from
  `https://packages.typst.org/preview/cetz-0.5.2.tar.gz`.
  Archive SHA256: `77cf8490114ae04c6e665a11efa691d284a0cadb9719771b5708c1197292f23f`.
  The unmodified archive and license are retained in `archives/`.
  Kurvst owns arrowhead geometry; no local CeTZ patches are applied.
- oxifmt 1.0.0, MIT OR Apache-2.0, from
  `https://packages.typst.org/preview/oxifmt-1.0.0.tar.gz`.
  Nix recursive hash: `sha256-RtGKdyiX2kJbUjChPohSGNeYOKVlI2VM0k1uFaEqDC8=`.
- MiTeX 0.2.6, Apache-2.0, from
  `https://packages.typst.org/preview/mitex-0.2.6.tar.gz`.

The corresponding license files and package manifests remain in each package directory.
